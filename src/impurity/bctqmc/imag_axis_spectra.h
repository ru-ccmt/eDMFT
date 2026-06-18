#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <blitz/array.h>
#include "basics.h"
#include "interpolate.h"
#include "tinyMatmul.h"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"



/*
inline double int_newt_cotes(const Array<double,1>& f, double a, double b, int N) noexcept
{
  // Implements Newton–Cotes integration rules when f is given on a discrete uniform mesh.
  //
  //                 Integrate[ f(x)dx , {x,a,b}]
  // 
  // User give number of points to be used, and interval [a,b], and we use high level 
  // routine to get as precise answer as possible for small price. Should be better than simpson
  double dh = (N>1) ? (b-a)/(N-1) : (b-a);
  if (N>=8){// O(dh^4)
    const double coef[]={17/48.,59/48.,43/48.,49./48.,1.0};
    double dsum=0;
    dsum += coef[0]*(f(0)+f(N-1)) + coef[1]*(f(1)+f(N-2)) + coef[2]*(f(2)+f(N-3)) + coef[3]*(f(3)+f(N-4));
    for (int i=4; i<N-4; ++i) dsum += f(i);
    return dsum*dh;
  }else{
    switch (N){
    case 1:
      return dh*f(0);
    case 2:
      return 0.5*dh*(f(0)+f(1));
    case 3:
      return (f(0)+4*f(1)+f(2))*dh/3.; // O(h^5)
    case 4:
      return (f(0)+3*(f(1)+f(2))+f(3))*3./8*dh; //O(h^5)
    case 5:
      return (7*f(0)+32*f(1)+12*f(2)+32*f(3)+7*f(4))*2./45.*dh; //O(h^7)
    case 6:
      return (19*f(0)+75*f(1)+50*(f(2)+f(3))+75*f(4)+19*f(5))*5./288.*dh; // O(h^7)
    case 7:
      return (41*f(0)+216*f(1)+27*f(2)+272*f(3)+27*f(4)+216*f(5)+41*f(6))/140*dh; //O(h^9)
    }
  }
  return 0;
}
*/

class VolterraSolver{
  // Here we solve DysonLaplace transformation with the double integral equation.
  // As usual, we write in matrix form:
  //                                    G = G0 + G0*Sigma*G
  // and we stay in time, hence all multiplications have to be in time and orbital indices.
  // We also take into account that G(t<0)=G0(t<0)=Sigma(t<0) vanishes, hence only the following
  // integration remains:
  //    G(tau) = G0(tau) + Integrate[G0(tau-t2)*Sigma(t2-t1)*G(t1),{t2,t1,tau},{t1,0,tau}]
  // This is equivalent to
  //    G(tau) = G0(tau) + Integrate[ conv(tau-t)*G(t), {t,0,tau}]
  // where
  //    conv(tau) = Integrate[ G0(tau-t)*Sig(t), {t,0,tau}]
  // We evaluate this double integral below, which is actually just twice single integral and takes O(n^2) time.
  // Notice that conv(0)=0, hence G(tau) does not appear on the right-hand side, but
  // only prior time G(tau-delta). The term with G(tau) has G(tau)*conv(0), which vanishes.
  // Also notice that we take 2^r*N points instead of N points for precision. As it turns out, the error
  // even when r=0 is of the order of 10^-6, which probably means that this is not needed.
  //
  //  This is called:  Volterra equation of the second kind.
  //
public:
  using MatrixRM = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
  struct RuntimeParams{
    // Empirical crossover between the simpler small-block path and the
    // optimized large-block path. This can later be made MPI-size dependent.
    int cmp2_small_cutoff = 4;
    bool print_breakdown = false;
  };
  struct GeoTerm {
    // Instead of computing integration coefficients several times, we precompute them in constructor, and store using this class.
    int l;  // panel of interpolation
    TinyVector<double,4> coeff; // 4 coefficients for interpolation, that are used as f_i*coeff[0]+f_{i+1}*coeff[1]+f''_i*coeff[2]+f''_{i+1}*coeff[3]
  };
  const Array<double,1>& tau;
  std::vector<std::vector<GeoTerm>> wij_geom;  // here we store precomputed part of the Volterra weights
  int nt; // number of time points
  RuntimeParams params;
  VolterraSolver(const Array<double,1>& _tau_, const RuntimeParams& p)
    : tau(_tau_), wij_geom(_tau_.extent(0)*_tau_.extent(0)), nt(_tau_.extent(0)), params(p)
  {
    precompute_geometry(); // precomputing common part of Volterra weights
  }
  VolterraSolver(const Array<double,1>& _tau_) : VolterraSolver(_tau_, RuntimeParams{}) {}
  int geom_index(int i, int j) const { return i*nt + j; }
  static std::pair<int,int> split_range(int irank, int nproc, int leng)
  { // used for inside parallelization when blocks are very large, so that they would block efficient outside parallelization
    int base = leng / nproc;
    int rem  = leng % nproc;
    int iqs  = irank * base + std::min(irank, rem);
    int iqe  = iqs + base + (irank < rem ? 1 : 0);
    return {iqs, iqe};
  }
  TinyVector<double,4> Li(int i, int j, int l, double a0, double b0, double a1, double b1) const
  {// This implements W^0_{ijl} = \integral_a^b C_l(x)*phi_{j,0}(t_i-x), where
    // phi_{j,0}(x) = (t_{j+1}-x)/(t_{j+1}-t_j) comes from the linear approximation of a function in a panel j,
    // and where C_l is cubic spline in panel l, and hence has a form:
    // C_l(x) = C_l * phi_{l,0}(x) + C_{l+1} * phi_{l,1}(x) + C''_l * phi_{l,2}(x) + C''_{l+1} * phi_{l,3}(x).
    // Here phi_{l,0}(x) = (t_{l+1}-x)/(t_{l+1}-t_l)  and
    //      phi_{l,1}(x) = (x-t_l)/(t_{l+1}-t_l)
    //      phi_{l,2}(x) = 1/6*(t_{l+1}-t_l)^2 * ( q^3 - q ) and q=(t_{l+1}-x)/(t_{l+1}-t_l)
    //      phi_{l,3}(x) = 1/6*(t_{l+1}-t_l)^2 * ( p^3 - p ) and p=(x-t_l)/(t_{l+1}-t_l)
    // Using the exact integration, we can write W^0_{ijl} as the following sum:
    // W^0_{ijl} = C_l * L0 + C_{l+1} * L1 + C''_l * L2 + C''_{l+1} * L3
    // and L0, L1, L2, L3 are evaluated below and returned as TinVector<double,4>(L0,L1,L2,L3).
    //  Here the input limits should be fractional, and the original integral from a to b is transformed to
    //    in L_0 and L2:
    //  a0 = (t_{l+1}-b)/(t_{l+1}-t_l)  and b0 = (t_{l+1}-a)/(t_{l+1}-t_l)
    //    and in L1 and L3:
    //  a1 = (a-t_l)/(t_{l+1}-t_l) and b1 = (b-t_l)/(t_{l+1}-t_l)
    //
    // limits:
    // 0<= l <= N-2
    // 0<= j <= N-2
    // 0<= i <= N-1
    // 
    // Notice that parameter gamma = (tau(l+1)-tau(i)+tau(j+1))/(tau(l+1)-tau(l))
    //   is always smaller than dt_j/dt_l = (tau(j+1)-tau(j))/(tau(l+1)-tau(l))
    //
    double dt_l = tau(l+1)-tau(l);
    double dt_j = tau(j+1)-tau(j);
    double gamma = (tau(l+1)-tau(i)+tau(j+1))/dt_l;
    double dt_l_2 = dt_l*dt_l;
    double dt_l_4 = dt_l_2*dt_l_2;
    double b0_2 = b0*b0, a0_2 = a0*a0;
    double b1_2 = b1*b1, a1_2 = a1*a1;
    double L0 = dt_l_2/dt_j*( b0_2 * (1./2.* gamma   -1./3.*b0) - a0_2 * (1./2.* gamma    -1./3.*a0) );
    double L1 = dt_l_2/dt_j*( b1_2 * (1./2.*(gamma-1)+1./3.*b1) - a1_2 * (1./2.*(gamma-1) +1./3.*a1) );
    double L2 = dt_l_4/dt_j*1./6.*( b0_2 * (  gamma   *( 0.25*b0_2 - 0.5) + b0*( 1./3.-1./5.*b0_2)) - a0_2 * (  gamma   *( 0.25*a0_2 - 0.5) + a0*( 1./3.-1./5.*a0_2)) );
    double L3 = dt_l_4/dt_j*1./6.*( b1_2 * ( (gamma-1)*( 0.25*b1_2 - 0.5) + b1*(-1./3.+1./5.*b1_2)) - a1_2 * ( (gamma-1)*( 0.25*a1_2 - 0.5) + a1*(-1./3.+1./5.*a1_2)) );
    return TinyVector<double,4>(L0,L1,L2,L3);
  }
  TinyVector<double,4> Ki(int i, int j, int l, double a0, double b0, double a1, double b1) const
  {// This implements W^1_{ijl} = \integral_a^b C_l(x)*phi_{j-1,1}(t_i-x), where
    // phi_{j-1,1}(x) = (x-t_{j-1})/(t_j-t_{j-1}) comes from the linear approximation of a function in a panel j,
    // and where C_l is cubic spline in panel l, and hence has a form:
    // C_l(x) = C_l * phi_{l,0}(x) + C_{l+1} * phi_{l,1}(x) + C''_l * phi_{l,2}(x) + C''_{l+1} * phi_{l,3}(x).
    // Here phi_{l,0}(x) = (t_{l+1}-x)/(t_{l+1}-t_l)  and
    //      phi_{l,1}(x) = (x-t_l)/(t_{l+1}-t_l)
    //      phi_{l,2}(x) = 1/6*(t_{l+1}-t_l)^2 * ( q^3 - q ) and q=(t_{l+1}-x)/(t_{l+1}-t_l)
    //      phi_{l,3}(x) = 1/6*(t_{l+1}-t_l)^2 * ( p^3 - p ) and p=(x-t_l)/(t_{l+1}-t_l)
    // Using the exact integration, we can write W^1_{ijl} as the following sum:
    // W^1_{ijl} = C_l * K0 + C_{l+1} * K1 + C''_l * K2 + C''_{l+1} * K3
    // and K0, K1, K2, K3 are evaluated below and returned as TinVector<double,4>(K0,K1,K2,K3).
    //  Here the input limits should be fractional, and the original integral from a to b is transformed to
    //    in K_0 and K2:
    //  a0 = (t_{l+1}-b)/(t_{l+1}-t_l)  and b0 = (t_{l+1}-a)/(t_{l+1}-t_l)
    //    and in K1 and K3:
    //  a1 = (a-t_l)/(t_{l+1}-t_l) and b1 = (b-t_l)/(t_{l+1}-t_l)
    //
    // limits:
    // 0<= l <= N-2
    // 1<= j <= N-1
    // 0<= i <= N-1
    //
    // Notice that parameter eps = (tau(i)-tau(j-1)-tau(l))/(tau(l+1)-tau(l))
    //   is always smaller than dt_{j-1}/dt_l = (tau(j)-tau(j-1))/(tau(l+1)-tau(l));
    //
    double dt_l = tau(l+1)-tau(l);
    double dt_j = tau(j)-tau(j-1);
    double eps = (tau(i)-tau(j-1)-tau(l))/dt_l;
    double dt_l_2 = dt_l*dt_l;
    double dt_l_4 = dt_l_2*dt_l_2;
    double b0_2 = b0*b0, a0_2 = a0*a0;
    double b1_2 = b1*b1, a1_2 = a1*a1;
    double K0 = dt_l_2/dt_j*( b0_2 * (1./2.*(eps-1) + 1./3.*b0) - a0_2 * (1./2.* (eps-1) + 1./3.*a0) );
    double K1 = dt_l_2/dt_j*( b1_2 * (1./2.* eps    - 1./3.*b1) - a1_2 * (1./2.* eps     - 1./3.*a1) );
    double K2 = dt_l_4/dt_j*1./6.*( b0_2 * ( (eps-1)*( 0.25*b0_2 - 0.5) + b0*(-1./3.+1./5.*b0_2)) - a0_2 * ( (eps-1)*( 0.25*a0_2 - 0.5) + a0*(-1./3.+1./5.*a0_2)) );
    double K3 = dt_l_4/dt_j*1./6.*( b1_2 * (  eps   *( 0.25*b1_2 - 0.5) + b1*( 1./3.-1./5.*b1_2)) - a1_2 * (  eps   *( 0.25*a1_2 - 0.5) + a1*( 1./3.-1./5.*a1_2)) );
    return TinyVector<double,4>(K0,K1,K2,K3);
  }
  void add_geom_term(int i, int j, int l, const TinyVector<double,4>& coeff)
  {
    wij_geom[geom_index(i,j)].push_back(GeoTerm{l, coeff});
  }
  void precompute_geometry()
  {// Precompute all purely geometrical information entering W(i,j).
    //
    // The expensive part in cmp2 was repeatedly determining:
    //   - which spline panels l overlap with a given (i,j) contribution
    //   - the Li/Ki coefficients for that overlap
    //
    // None of this depends on:
    //   - Sigma
    //   - fconv
    //   - orbital dimension N
    //
    // It depends only on the external mesh tau and on the pair (i,j), so we
    // compute it once in the constructor and reuse it in every cmp2 call.
    //
    // For each (i,j) we store a list of terms:
    //   W(i,j) = sum_term coeff_term[0..3] * spline_panel(l_term)
    //
    // Later cmp2 only needs to loop over these stored terms and call
    // fconv.evaluate_add(l, coeff, Wij(i,j,:)).
    for (int i=1; i<nt; ++i){
      int ai0=i;
      {
        int j=0;
       // For j=0:
        //   W(i,0) = \int_{t_0}^{t_1} dx conv(t_i-x) phi_{0,0}(x)
        //
        // This contributes only the Li-type terms.
         double tij_0 = tau(i)-tau(j+1);
        double tij_1 = tau(i)-tau(j);
        intpar p_l0 = InterpRight(tij_0, ai0, tau);
        intpar p_l1(i-1,1);
        int l0 = p_l0.i;
        int ln = p_l1.i;
        if (l0<ln){
          add_geom_term(i,j,l0, Li(i,j,l0, 0, 1-p_l0.p, p_l0.p, 1));
          for (int l=l0+1; l<ln; ++l) add_geom_term(i,j,l, Li(i,j,l, 0, 1, 0, 1));
          add_geom_term(i,j,ln, Li(i,j,ln, 1-p_l1.p, 1, 0, p_l1.p));
        }else if (l0==ln){
          add_geom_term(i,j,l0, Li(i,j,l0, 1-p_l1.p, 1-p_l0.p, p_l0.p, p_l1.p));
        }
      }
      for (int j=1; j<i; ++j){
       // For 0<j<i:
        //   W(i,j) = \int_{t_j}^{t_{j+1}} dx conv(t_i-x) phi_{j,0}(x)
        //          + \int_{t_{j-1}}^{t_j} dx conv(t_i-x) phi_{j-1,1}(x)
        //
        // The first integral gives Li-type terms.
        // The second integral gives Ki-type terms.
        double tij_0 = tau(i)-tau(j+1);
        double tij_1 = tau(i)-tau(j);
        double tij_2 = tau(i)-tau(j-1);
        intpar p_l0 = InterpRight(tij_0, ai0, tau);
        int ai1 = p_l0.i;
        intpar p_l1 = InterpLeft(tij_1, ai1, tau);
        intpar p_l2 = InterpLeft(tij_2, ai1, tau);
        int l0 = p_l0.i;
        int ln = p_l1.i;
        int lm = p_l2.i;
        if (l0<ln){
          add_geom_term(i,j,l0, Li(i,j,l0, 0, 1-p_l0.p, p_l0.p, 1));
          for (int l=l0+1; l<ln; ++l) add_geom_term(i,j,l, Li(i,j,l, 0, 1, 0, 1));
          add_geom_term(i,j,ln, Li(i,j,ln, 1-p_l1.p, 1, 0, p_l1.p));
        }else if (l0==ln){
          add_geom_term(i,j,l0, Li(i,j,l0, 1-p_l1.p, 1-p_l0.p, p_l0.p, p_l1.p));
        }
        if (ln<lm){
          add_geom_term(i,j,ln, Ki(i,j,ln, 0, 1-p_l1.p, p_l1.p, 1));
          for (int l=ln+1; l<lm; ++l) add_geom_term(i,j,l, Ki(i,j,l, 0, 1, 0, 1));
          add_geom_term(i,j,lm, Ki(i,j,lm, 1-p_l2.p, 1, 0, p_l2.p));
        }else if (ln==lm){
          add_geom_term(i,j,ln, Ki(i,j,ln, 1-p_l2.p, 1-p_l1.p, p_l1.p, p_l2.p));
        }
      }
      {
        int j=i;
        // For j=i:
        //   W(i,i) = \int_{t_{i-1}}^{t_i} dx conv(t_i-x) phi_{i-1,1}(x)
        //
        // This contributes only Ki-type terms.
        double tij_2 = tau(i)-tau(j-1);
        intpar p_l1(0,0);
        int ai1=0;
        intpar p_l2 = InterpLeft(tij_2, ai1, tau);
        int ln = p_l1.i;
        int lm = p_l2.i;
        if (ln<lm){
          add_geom_term(i,j,ln, Ki(i,j,ln, 0, 1-p_l1.p, p_l1.p, 1));
          for (int l=ln+1; l<lm; ++l) add_geom_term(i,j,l, Ki(i,j,l, 0, 1, 0, 1));
          add_geom_term(i,j,lm, Ki(i,j,lm, 1-p_l2.p, 1, 0, p_l2.p));
        }else if (ln==lm){
          add_geom_term(i,j,ln, Ki(i,j,ln, 1-p_l2.p, 1-p_l1.p, p_l1.p, p_l2.p));
        }
      }
    }
  }
  void accumulate_term_matrix_packed(const std::vector<double>& Fpack, const std::vector<double>& F2pack, int N2,
                                     int l, const TinyVector<double,4>& coeff, MatrixRM& W) const
  {
    double* w = W.data();
    const int off0 = l*N2;
    const int off1 = (l+1)*N2;
    for (int k=0; k<N2; ++k){
      w[k] += coeff[0]*Fpack[off0 + k]
           +  coeff[1]*Fpack[off1 + k]
           +  coeff[2]*F2pack[off0 + k]
           +  coeff[3]*F2pack[off1 + k];
    }
  }
  void cmp1(Array<double,1>& Gm, const Array<double,1>& Sm, double Em){
   //---------------------------------------------------------------------------
    // Scalar Σ(t) → G(t)
    //---------------------------------------------------------------------------
    // This is the simplest case in which Sm,Gm are a single function of time Sm(tau), Gm(tau), as well as Em, defined by G0=exp(-Em*tau).
    // Appropriate for 1D blocks
     Spline1D<double> fSig(tau.extent(0));
    for (int i=0; i<tau.extent(0); ++i) fSig[i]=Sm(i);
    fSig.splineIt(tau);
    // We first compute conv(t) = int_0^t G0(t-x) \Sigma(x) dx
    //
    // Because G0(t) = exp(-Em*t) we can derive simple recursive formula:
    //   conv(t_{i+1}) = conv(t_i)*exp(-(t_{i+1}-t_i)*Em) + \int_{t_i}^{t_{i+1}}dx e^{-(t_{i+1}-x)*Em} \Sigma(x)
    // Since \Sigma is given as spline, we can express it in each panel as \Sigma(x)= S_i \phi_{i,0}(x) + S_{i+1} \phi_{i,1}(x) + S''_{i} \phi_{i,2}(x) + S''_{i+1} \phi_{i,3}(x)
    // where \phi_{i,0}(x) = (t_{i+1}-x)/(t_{i+1}-t_i)
    //       \phi_{i,1}(x) = (x-t_i)/(t_{i+1}-t_i)
    //       \phi_{i,2}(x) = 1/6*(t_{i+1}-t_i)^2 * (q^3-q) where q=(t_{i+1}-x)/(t_{i+1}-t_i)
    //       \phi_{i,3}(x) = 1/6*(t_{i+1}-t_i)^2 * (p^3-p) where p=(x-t_i)/(t_{i+1}-t_i)
    // The analytic integration over each panel of the spline gives
    //  \int_{t_i}^{t_{i+1}}dx e^{-(t_{i+1}-x)*Em} \Sigma(x) = Li[0] * S_i + Li[1] * S_{i+1} + Li[2] * S''_{i} + Li[3] * S''_{i+1}
    //
    Array<double,1> convv(tau.extent(0)); convv=0;
    for (int i=0; i<tau.extent(0)-1; ++i){
      double dt = (tau(i+1)-tau(i));
      double a = Em*dt;
      double ema = exp(-a), a2=a*a, a4=a*a*a*a;
      double dt3 = dt*dt*dt;
      TinyVector<double,4> Liw;
      if (fabs(a)>1e-7){
        Liw[0] = dt*(1-ema*(1+a))/a2;
        Liw[1] = dt*(a-1+ema)/a2;
        Liw[2] = -dt3/a4*((1+a+a2/3.)*ema  + a2/6.-1);
        Liw[3] = -dt3/a4*((1-a+a2/3.) + ema*(a2/6.-1));
      }else{ // Taylor expansion for the case of small a, which would otherwise diverge
        Liw[0] = dt*(0.5-a/3.);
        Liw[1] = dt*(0.5-a/6.);
        Liw[2] = dt3*(-1./24.+a/45.);
        Liw[3] = dt3*(-1./24.+a*7./360.);
      }
      convv(i+1) = convv(i)*ema + fSig.evaluate(i,Liw);
    }
    // Next we spline conv(t) = int_0^t G0(t-x) \Sigma(x) dx, so that convolution with G0 and G can be computed next.
    // First we evaluate convolution of conv*G0, i.e.,
    // Rx(t) = \int_0^t conv(t-x) G0(x) = \int_0^t conv(x) * exp(-(t-x)*Em), which can again be handled with the same trick as convolution of Sigma*G0 above.
    //
    // The analytic integration over each panel of the spline gives
    //  dRx_i = \int_{t_i}^{t_{i+1}}dx e^{-(t_{i+1}-x)*Em} conv(x) = Li[0] * conv_i + Li[1] * conv_{i+1} + Li[2] * conv''_{i} + Li[3] * conv''_{i+1}
    //  and Rx_{i+1} = Rx_i *exp(-Em*(t_{i+1}-t_i)) + dRx_i
     Spline1D<double> fconv(tau.extent(0));
    for (int i=0; i<tau.extent(0); ++i) fconv[i]=convv(i);
    fconv.splineIt(tau);
    Array<double,1> Rx(tau.extent(0)); Rx=0;
    for (int i=0; i<tau.extent(0)-1; ++i){
      double dt = (tau(i+1)-tau(i));
      double a = Em*dt;
      double ema = exp(-a), a2=a*a, a4=a*a*a*a;
      double dt3 = dt*dt*dt;
      TinyVector<double,4> Liw;
      if (fabs(a)>1e-7){
        Liw[0] = dt*(1-ema*(1+a))/a2;
        Liw[1] = dt*(a-1+ema)/a2;
        Liw[2] = -dt3/a4*((1+a+a2/3.)*ema  + a2/6.-1);
        Liw[3] = -dt3/a4*((1-a+a2/3.) + ema*(a2/6.-1));
      }else{
        Liw[0] = dt*(0.5-a/3.);
        Liw[1] = dt*(0.5-a/6.);
        Liw[2] = dt3*(-1./24.+a/45.);
        Liw[3] = dt3*(-1./24.+a*7./360.);
      }
      Rx(i+1) = Rx(i)*ema + fconv.evaluate(i,Liw);
    }
   // Finally, we compute weights Wij, which can solve Volterra equation.
    // We write dG(t) = G(t)-G0(t) and write equation for dG(t)
    //   dG(t) = \int_0^\tau dx conv(t-x) dG(x) + Rx(t)
    // We want to compute weights W(i,j), so that we can write:
    //  dG(t_i) = Rx(t_i) + \sum_{j=0}^{j=i} W(i,j) * dG(t_j)
    // which can be solved as
    //  dG(t_i) = (1-W(i,i))^{-1} * ( R_x(t_i) + \sum_{j=0}^{j=i-1} W(i,j) * dG(t_j))
    // and dG(t_i) can be computed from previous dG(t<t_i)
    //
    // We approximate the resulting function dG(t) with linear approximation, i.e., in panel j: dG(x) = dG_j * \phi_{j,0}(x) + dG_{j+1}*\phi_{j,1}(x)
    // This results in
    // dG(t_i) = \sum_{j=0}^{j=i-1} \int_{t_j}^{t_{j+1}}dx conv(t_i-x) dG(x) = \sum_{j=0}^{j=i-1} \int_{t_j}^{t_{j+1}}dx conv(t_i-x) (dG_j * \phi_{j,0}(x) + dG_{j+1} * \phi_{j,1}(x))
    // which can also be written
    // dG(t_i) = \sum_{j=0}^{j=i-1} dG_j * \int_{t_j}^{t_{j+1}}dx conv(t_i-x) \phi_{j,0}(x) + \sum_{j=1}^{j=i} dG_{j} * \int_{t_{j-1}}^{t_{j}}dx conv(t_i-x) \phi_{j-1,1}(x)
    // and therefore (for 0<j<i):
    //              W(i,j) = \int_{t_j}^{t_{j+1}}dx conv(t_i-x) \phi_{j,0}(x) + \int_{t_{j-1}}^{t_{j}}dx conv(t_i-x) \phi_{j-1,1}(x)
    //              W(i,0) = \int_{t_0}^{t_{1}}dx conv(t_i-x) \phi_{0,0}(x)
    //              W(i,i) = \int_{t_{i-1}}^{t_{i}}dx conv(t_i-x) \phi_{i-1,1}(x)
    //
    // Next we take into account that conv(x) is a cubic spline, and hence can also be written in each panel as cubic polinomial
    // We write W(i,j) = \sum_l (Li[i,j,l][0] + Ki[i,j,l][0])* conv_l + (Li[i,j,l][1] + Ki[i,j,l][1]) * conv_{l+1} + (Li[i,j,l][2] + K[i,j,l][2]) * conv''_l + (Li[i,j,l][3] + K[i,j,l][3])* conv''_{l+1}
    //   here Li evaluates the first integral for W(i,j) 
    //    Li[i,j,:] = \int_{t_j}^{t_{j+1}}dx conv(t_i-x) \phi_{j,0}(x)
    //   and Ki evaluates the second integral for W(i,j)
    //    Ki[i,j,:] = \int_{t_{j-1}}^{t_{j}}dx conv(t_i-x) \phi_{j-1,1}(x)
    //  After a long derivation, we obtain Li and Ki as implemented above
    Array<double,2> Wij(tau.extent(0),tau.extent(0)); Wij=0;
    if (true){
      // Reuse the geometry precomputed in the constructor. The scalar cmp1
      // needs exactly the same (i,j,l) panel decomposition as cmp2_large,
      // only the final evaluation is scalar rather than matrix-valued.
      for (int i=1; i<tau.extent(0); ++i){
	for (int j=0; j<=i; ++j){
	  double W_ij = 0.0;
	  const auto& terms = wij_geom[geom_index(i,j)];
	  for (const auto& term : terms) W_ij += fconv.evaluate(term.l, term.coeff);
	  Wij(i,j) = W_ij;
	}
      }
    }else{
      // Original explicit construction kept as fallback/reference.
      for (int i=1; i<tau.extent(0); ++i){
        int ai0=i;
        {
	  int j=0;
	  double tij_0 = tau(i)-tau(j+1); // tau(i)-tau(1)
	  double tij_1 = tau(i)-tau(j);   // tau(i)
	  intpar p_l0 = InterpRight(tij_0, ai0, tau);
	  intpar p_l1(i-1,1); // the tau(i) can be choosen as (i,p=0) or (i-1,p=1);
	  int l0 = p_l0.i;
	  int ln = p_l1.i;
	  double W_ij=0;
	  if (l0<ln){
	    TinyVector<double,4> L1 = Li(i,j,l0, 0, 1-p_l0.p, p_l0.p, 1);
	    W_ij += fconv.evaluate(l0,L1); 
	    for (int l=l0+1; l<ln; ++l){
	      TinyVector<double,4> L = Li(i,j,l, 0, 1, 0, 1);
	      W_ij += fconv.evaluate(l,L); 
	    }
	    TinyVector<double,4> L2 = Li(i,j,ln, 1-p_l1.p, 1, 0, p_l1.p);
	    W_ij += fconv.evaluate(ln,L2); 
	  }else if (l0==ln){
	    TinyVector<double,4> L = Li(i,j,l0, 1-p_l1.p, 1-p_l0.p, p_l0.p, p_l1.p);
	    W_ij += fconv.evaluate(l0,L);
	  }else{
	    cerr<<"ERROR in cmp1 l0>ln"<<endl;// should not happen
	  }
	  Wij(i,j)=W_ij;
        }
        for (int j=1; j<i; ++j){
	  double tij_0 = tau(i)-tau(j+1);
	  double tij_1 = tau(i)-tau(j);
	  double tij_2 = tau(i)-tau(j-1);
	  intpar p_l0 = InterpRight(tij_0, ai0, tau);
	  int ai1 = p_l0.i;
	  intpar p_l1 = InterpLeft(tij_1, ai1, tau);
	  intpar p_l2 = InterpLeft(tij_2, ai1, tau);
	  int l0 = p_l0.i;
	  int ln = p_l1.i;
	  int lm = p_l2.i;
	  double W_ij=0;
	  if (l0<ln){
	    TinyVector<double,4> L1 = Li(i,j,l0, 0, 1-p_l0.p, p_l0.p, 1);
	    W_ij += fconv.evaluate(l0,L1);
	    for (int l=l0+1; l<ln; ++l){
	      TinyVector<double,4> L = Li(i,j,l, 0, 1, 0, 1);
	      W_ij += fconv.evaluate(l,L);
	    }
	    TinyVector<double,4> L2 = Li(i,j,ln, 1-p_l1.p, 1, 0, p_l1.p);
	    W_ij += fconv.evaluate(ln,L2);
	  }else if (l0==ln){
	    TinyVector<double,4> L = Li(i,j,l0, 1-p_l1.p, 1-p_l0.p, p_l0.p, p_l1.p);
	    W_ij += fconv.evaluate(l0,L);
	  }else{
	    cerr<<"ERROR in cmp1 l0>ln"<<endl;// should not happen
	  }
	  if (ln<lm){
	    TinyVector<double,4> K1 = Ki(i,j,ln, 0, 1-p_l1.p, p_l1.p, 1);
	    W_ij += fconv.evaluate(ln,K1);
	    for (int l=ln+1; l<lm; ++l){
	      TinyVector<double,4> K = Ki(i,j,l, 0, 1, 0, 1);
	      W_ij += fconv.evaluate(l,K);
	    }
	    TinyVector<double,4> K2 = Ki(i,j,lm, 1-p_l2.p, 1, 0, p_l2.p);
	    W_ij += fconv.evaluate(lm,K2);
	  }else if (ln==lm){
	    TinyVector<double,4> K = Ki(i,j,ln, 1-p_l2.p, 1-p_l1.p, p_l1.p, p_l2.p);
	    W_ij += fconv.evaluate(ln,K);
	  }else{
	    cerr<<"ERROR in cmp1 ln>lm"<<endl; // should not happen
	  }
	  Wij(i,j) = W_ij;
        }
        {
	  int j=i;
	  double tij_1 = tau(i)-tau(j);   // tij_1=0
	  double tij_2 = tau(i)-tau(j-1); // tij_2=dt_{i-1}
	  intpar p_l1(0,0); // tij_1=0
	  int ai1=0;
	  intpar p_l2 = InterpLeft(tij_2, ai1, tau);
	  int ln = p_l1.i;
	  int lm = p_l2.i;
	  double W_ij=0;
	  if (ln<lm){
	    TinyVector<double,4> K1 = Ki(i,j,ln, 0, 1-p_l1.p, p_l1.p, 1);
	    W_ij += fconv.evaluate(ln,K1);
	    for (int l=ln+1; l<lm; ++l){
	      TinyVector<double,4> K = Ki(i,j,l, 0, 1, 0, 1);
	      W_ij += fconv.evaluate(l,K);
	    }
	    TinyVector<double,4> K2 = Ki(i,j,lm, 1-p_l2.p, 1, 0, p_l2.p);
	    W_ij += fconv.evaluate(lm,K2);
	  }else if (ln==lm){
	    TinyVector<double,4> K = Ki(i,j,ln, 1-p_l2.p, 1-p_l1.p, p_l1.p, p_l2.p);
	    W_ij += fconv.evaluate(ln,K);
	  }else{
	    cerr<<"ERROR in cmp1 ln>lm"<<endl; // should not happen
	  }
	  Wij(i,j) = W_ij;
        }
      }
    }
    Array<double,1> dG(tau.extent(0)); dG=0;
    for (int i=1; i<tau.extent(0); ++i){
      double CG=0;
      for (int j=0; j<i; ++j) CG += Wij(i,j)*dG(j);
      // Note that dG(i) = Rx(i) + \sum_{j=0}^{j=i} Wij(i,j)*dG(j), but now dG(i) appear on the left and right side.
      // When eliminating dG(i) on the right, we get
      //       (1-Wij(i,i))  * dG(i) = Rx(i) + \sum_{j=0}^{j=i-1} Wij(i,j)*dG(j)
      dG(i) = (CG + Rx(i))/(1-Wij(i,i));
    }
    Gm.resize(tau.extent(0));
    for (int i=0; i<tau.extent(0); ++i) Gm(i) = dG(i) + exp(-Em*tau(i));
  }
  void cmp2_small(Array<double,3>& Gm, const Array<double,3>& Sm, const Array<double,1>& Em){
   //---------------------------------------------------------------------------
    // Matrix Σ(t) (N×N) diagonal G0 → G(t) (N×N)      𝑂(N_orb^3 n^2)
    //---------------------------------------------------------------------------
    // Appropriate for blocks of
    //                  Gm(tau,N,N) and Sm(N,N,tau).
    // We work in basis with diagonal G0(N,tau).
    //  temporary objects:
    //        Gr(t,N,N), convr(t,N,N), fr(t,N,N)
    //        in addition to G0(t), Sig(t), 
    // For algorithm, please check cmp1 above, which has a lot of comments.
    const int N = Sm.extent(0);
    Spline1DMany fSig(tau.extent(0),N*N);
    for (int i1=0; i1<N; ++i1)
      for (int i2=0; i2<N; ++i2){
        int icmp=i1*N+i2;
        for (int it=0; it<tau.extent(0); ++it) fSig(it,icmp)=Sm(i1,i2,it);
      }
    fSig.splineIt(tau);
    // We first compute conv(t) = int_0^t G0(t-x) \Sigma(x) dx
    //
    // Because G0(t) = exp(-Em*t) we can derive simple recursive formula:
    //   conv(t_{i+1}) = conv(t_i)*exp(-(t_{i+1}-t_i)*Em) + \int_{t_i}^{t_{i+1}}dx e^{-(t_{i+1}-x)*Em} \Sigma(x)
    // Since \Sigma is given as spline, we can express it in each panel as \Sigma(x)= S_i \phi_{i,0}(x) + S_{i+1} \phi_{i,1}(x) + S''_{i} \phi_{i,2}(x) + S''_{i+1} \phi_{i,3}(x)
    // where \phi_{i,0}(x) = (t_{i+1}-x)/(t_{i+1}-t_i)
    //       \phi_{i,1}(x) = (x-t_i)/(t_{i+1}-t_i)
    //       \phi_{i,2}(x) = 1/6*(t_{i+1}-t_i)^2 * (q^3-q) where q=(t_{i+1}-x)/(t_{i+1}-t_i)
    //       \phi_{i,3}(x) = 1/6*(t_{i+1}-t_i)^2 * (p^3-p) where p=(x-t_i)/(t_{i+1}-t_i)
    // The analytic integration over each panel of the spline gives
    //  \int_{t_i}^{t_{i+1}}dx e^{-(t_{i+1}-x)*Em} \Sigma(x) = Li[0] * S_i + Li[1] * S_{i+1} + Li[2] * S''_{i} + Li[3] * S''_{i+1}
    //
    Spline1DMany fconv(tau.extent(0),N*N);
    for (int it=0; it<tau.extent(0); ++it)
      for (int icmp=0; icmp<N*N; ++icmp) fconv(it,icmp) = 0.0;
    for (int i=0; i<tau.extent(0)-1; ++i){
      double dt = (tau(i+1)-tau(i));
      for (int i1=0; i1<N; ++i1){
        double a = Em(i1)*dt, ema = exp(-a), a2=a*a, a4=a*a*a*a;
        double dt3 = dt*dt*dt;
        TinyVector<double,4> Liw;
        if (fabs(a)>1e-7){
          Liw[0] = dt*(1-ema*(1+a))/a2; Liw[1] = dt*(a-1+ema)/a2;
          Liw[2] = -dt3/a4*((1+a+a2/3.)*ema  + a2/6.-1);
          Liw[3] = -dt3/a4*((1-a+a2/3.) + ema*(a2/6.-1));
        }else{
          Liw[0] = dt*(0.5-a/3.); Liw[1] = dt*(0.5-a/6.);
          Liw[2] = dt3*(-1./24.+a/45.); Liw[3] = dt3*(-1./24.+a*7./360.);
        }
        for (int i2=0; i2<N; ++i2){
          int icmp = i1*N+i2;
          fconv(i+1,icmp) = fconv(i,icmp)*ema + fSig.evaluate(i,icmp,Liw);
        }
      }
    }
    fconv.splineIt(tau);
    std::vector<double> Fpack(tau.extent(0) * N * N), F2pack(tau.extent(0) * N * N);
    for (int l=0; l<tau.extent(0); ++l){
      int off = l * N * N;
      for (int k=0; k<N*N; ++k){ Fpack[off+k]=fconv.F(l,k); F2pack[off+k]=fconv.F2(l,k); }
    }
    // Next we spline conv(t) = int_0^t G0(t-x) \Sigma(x) dx, so that convolution with G0 and G can be computed next.
    // First we evaluate convolution of conv*G0, i.e.,
    // Rx(t) = \int_0^t conv(t-x) G0(x) = \int_0^t conv(x) * exp(-(t-x)*Em), which can again be handled with the same trick as convolution of Sigma*G0 above.
    //
    // The analytic integration over each panel of the spline gives
    //  dRx_i = \int_{t_i}^{t_{i+1}}dx e^{-(t_{i+1}-x)*Em} conv(x) = Li[0] * conv_i + Li[1] * conv_{i+1} + Li[2] * conv''_{i} + Li[3] * conv''_{i+1}
    //  and Rx_{i+1} = Rx_i *exp(-Em*(t_{i+1}-t_i)) + dRx_i
    Array<double,3> Rx(tau.extent(0),N,N); Rx=0;
    for (int i=0; i<tau.extent(0)-1; ++i){
      for (int i2=0; i2<N; ++i2){
        double dt = (tau(i+1)-tau(i));
        double a = Em(i2)*dt, ema = exp(-a), a2=a*a, a4=a*a*a*a;
        double dt3 = dt*dt*dt;
        TinyVector<double,4> Liw;
        if (fabs(a)>1e-7){
          Liw[0] = dt*(1-ema*(1+a))/a2; Liw[1] = dt*(a-1+ema)/a2;
          Liw[2] = -dt3/a4*((1+a+a2/3.)*ema  + a2/6.-1);
          Liw[3] = -dt3/a4*((1-a+a2/3.) + ema*(a2/6.-1));
        }else{
          Liw[0] = dt*(0.5-a/3.); Liw[1] = dt*(0.5-a/6.);
          Liw[2] = dt3*(-1./24.+a/45.); Liw[3] = dt3*(-1./24.+a*7./360.);
        }
        for (int i1=0; i1<N; ++i1){
          int icmp = i1*N+i2;
          Rx(i+1,i1,i2) = Rx(i,i1,i2)*ema + fconv.evaluate(i,icmp,Liw);
        }
      }
    }
   // Finally, we compute weights Wij, which can solve Volterra equation.
    // We write dG(t) = G(t)-G0(t) and write equation for dG(t)
    //   dG(t) = \int_0^\tau dx conv(t-x) dG(x) + Rx(t)
    // We want to compute weights W(i,j), so that we can write:
    //  dG(t_i) = Rx(t_i) + \sum_{j=0}^{j=i} W(i,j) * dG(t_j)
    // which can be solved as
    //  dG(t_i) = (1-W(i,i))^{-1} * ( R_x(t_i) + \sum_{j=0}^{j=i-1} W(i,j) * dG(t_j))
    // and dG(t_i) can be computed from previous dG(t<t_i)
    //
    // We approximate the resulting function dG(t) with linear approximation, i.e., in panel j: dG(x) = dG_j * \phi_{j,0}(x) + dG_{j+1}*\phi_{j,1}(x)
    // This results in
    // dG(t_i) = \sum_{j=0}^{j=i-1} \int_{t_j}^{t_{j+1}}dx conv(t_i-x) dG(x) = \sum_{j=0}^{j=i-1} \int_{t_j}^{t_{j+1}}dx conv(t_i-x) (dG_j * \phi_{j,0}(x) + dG_{j+1} * \phi_{j,1}(x))
    // which can also be written
    // dG(t_i) = \sum_{j=0}^{j=i-1} dG_j * \int_{t_j}^{t_{j+1}}dx conv(t_i-x) \phi_{j,0}(x) + \sum_{j=1}^{j=i} dG_{j} * \int_{t_{j-1}}^{t_{j}}dx conv(t_i-x) \phi_{j-1,1}(x)
    // and therefore (for 0<j<i):
    //              W(i,j) = \int_{t_j}^{t_{j+1}}dx conv(t_i-x) \phi_{j,0}(x) + \int_{t_{j-1}}^{t_{j}}dx conv(t_i-x) \phi_{j-1,1}(x)
    //              W(i,0) = \int_{t_0}^{t_{1}}dx conv(t_i-x) \phi_{0,0}(x)
    //              W(i,i) = \int_{t_{i-1}}^{t_{i}}dx conv(t_i-x) \phi_{i-1,1}(x)
    //
    // Next we take into account that conv(x) is a cubic spline, and hence can also be written in each panel as cubic polinomial
    // We write W(i,j) = \sum_l (Li[i,j,l][0] + Ki[i,j,l][0])* conv_l + (Li[i,j,l][1] + Ki[i,j,l][1]) * conv_{l+1} + (Li[i,j,l][2] + K[i,j,l][2]) * conv''_l + (Li[i,j,l][3] + K[i,j,l][3])* conv''_{l+1}
    //   here Li evaluates the first integral for W(i,j) 
    //    Li[i,j,:] = \int_{t_j}^{t_{j+1}}dx conv(t_i-x) \phi_{j,0}(x)
    //   and Ki evaluates the second integral for W(i,j)
    //    Ki[i,j,:] = \int_{t_{j-1}}^{t_{j}}dx conv(t_i-x) \phi_{j-1,1}(x)
    //  After a long derivation, we obtain Li and Ki as implemented above
    Array<double,3> Wij(tau.extent(0),tau.extent(0),N*N); Wij=0;
    for (int i=1; i<tau.extent(0); ++i)
      for (int j=0; j<=i; ++j){
        MatrixRM Wj = MatrixRM::Zero(N,N);
        const auto& terms = wij_geom[geom_index(i,j)];
        for (const auto& term : terms) accumulate_term_matrix_packed(Fpack, F2pack, N*N, term.l, term.coeff, Wj);
        std::copy(Wj.data(), Wj.data()+N*N, Wij(i,j,Range::all()).dataFirst());
      }
    MatrixRM Kii(N,N), CG(N,N);
    std::vector<MatrixRM> dG(tau.extent(0), MatrixRM::Zero(N,N));
    for (int i=1; i<tau.extent(0); ++i){
      CG.setZero();
      for (int j=0; j<i; ++j){
        auto Wij_view = Wij(i,j,Range::all());
        Eigen::Map<const MatrixRM> Wij_map(Wij_view.dataFirst(), N, N);
        CG.noalias() += Wij_map * dG[j]; // CG[a,b] = \sum_{j<i} Wij[i,j][a,a']*dG[j][a',b]
      }
      {
        auto Rx_view = Rx(i,Range::all(),Range::all());
        Eigen::Map<const MatrixRM> Rx_map(Rx_view.dataFirst(), N, N);
        CG += Rx_map; // CG[a,b] += Rx[i][a,b]
      }
      // Note that dG(i) = Rx(i) + \sum_{j=0}^{j=i} Wij(i,j)*dG(j), but now dG(i) appear on the left and right side.
      // When eliminating dG(i) on the right, we get
      //       (1-Wij(i,i))  * dG(i) = Rx(i) + \sum_{j=0}^{j=i-1} Wij(i,j)*dG(j)
      auto Wii_view = Wij(i,i,Range::all());
      Eigen::Map<const MatrixRM> Wii_map(Wii_view.dataFirst(), N, N);
      Kii.setIdentity();
      Kii -= Wii_map;  // Kii=(1-W(i,i,:,:))
      Eigen::PartialPivLU<MatrixRM> lu(Kii);  // dG[i][a,b] = (Kii^{-1}*CG)_{a,b}
      dG[i] = lu.solve(CG);
    }
    Gm.resize(tau.extent(0),N,N); Gm = 0;
    for (int i=0; i<tau.extent(0); ++i)
      for (int i1=0; i1<N; ++i1){
        for (int i2=0; i2<N; ++i2) Gm(i,i1,i2) = dG[i](i1,i2);
        Gm(i,i1,i1) += exp(-Em(i1)*tau(i));
      }
  }
  void cmp2_large(Array<double,3>& Gm, const Array<double,3>& Sm, const Array<double,1>& Em)
  {
    // Equivalent to cmp2_small, but this is optimized for large blocks NxN with N>params.cmp2_small_cutoff, and also allows inside mpi parallelization.
    //
    const int N = Sm.extent(0);
    Spline1DMany fSig(tau.extent(0),N*N);
    for (int i1=0; i1<N; ++i1)
      for (int i2=0; i2<N; ++i2){
        int icmp=i1*N+i2;
        for (int it=0; it<tau.extent(0); ++it) fSig(it,icmp)=Sm(i1,i2,it);
      }
    fSig.splineIt(tau);
    Spline1DMany fconv(tau.extent(0),N*N);
    for (int it=0; it<tau.extent(0); ++it)
      for (int icmp=0; icmp<N*N; ++icmp) fconv(it,icmp) = 0.0;
    for (int i=0; i<tau.extent(0)-1; ++i){
      double dt = (tau(i+1)-tau(i));
      for (int i1=0; i1<N; ++i1){
        double a = Em(i1)*dt, ema = exp(-a), a2=a*a, a4=a*a*a*a;
        double dt3 = dt*dt*dt;
        TinyVector<double,4> Liw;
        if (fabs(a)>1e-7){
          Liw[0] = dt*(1-ema*(1+a))/a2; Liw[1] = dt*(a-1+ema)/a2;
          Liw[2] = -dt3/a4*((1+a+a2/3.)*ema  + a2/6.-1);
          Liw[3] = -dt3/a4*((1-a+a2/3.) + ema*(a2/6.-1));
        }else{
          Liw[0] = dt*(0.5-a/3.); Liw[1] = dt*(0.5-a/6.);
          Liw[2] = dt3*(-1./24.+a/45.); Liw[3] = dt3*(-1./24.+a*7./360.);
        }
        for (int i2=0; i2<N; ++i2){
          int icmp = i1*N+i2;
          fconv(i+1,icmp) = fconv(i,icmp)*ema + fSig.evaluate(i,icmp,Liw);
        }
      }
    }
    fconv.splineIt(tau);
    std::vector<double> Fpack(tau.extent(0) * N * N), F2pack(tau.extent(0) * N * N);
    for (int l=0; l<tau.extent(0); ++l){
      int off = l * N * N;
      for (int k=0; k<N*N; ++k){ Fpack[off+k]=fconv.F(l,k); F2pack[off+k]=fconv.F2(l,k); }
    }
    Array<double,3> Rx(tau.extent(0),N,N); Rx=0;
    for (int i=0; i<tau.extent(0)-1; ++i){
      for (int i2=0; i2<N; ++i2){
        double dt = (tau(i+1)-tau(i));
        double a = Em(i2)*dt, ema = exp(-a), a2=a*a, a4=a*a*a*a;
        double dt3 = dt*dt*dt;
        TinyVector<double,4> Liw;
        if (fabs(a)>1e-7){
          Liw[0] = dt*(1-ema*(1+a))/a2; Liw[1] = dt*(a-1+ema)/a2;
          Liw[2] = -dt3/a4*((1+a+a2/3.)*ema  + a2/6.-1);
          Liw[3] = -dt3/a4*((1-a+a2/3.) + ema*(a2/6.-1));
        }else{
          Liw[0] = dt*(0.5-a/3.); Liw[1] = dt*(0.5-a/6.);
          Liw[2] = dt3*(-1./24.+a/45.); Liw[3] = dt3*(-1./24.+a*7./360.);
        }
        for (int i1=0; i1<N; ++i1){
          int icmp = i1*N+i2;
          Rx(i+1,i1,i2) = Rx(i,i1,i2)*ema + fconv.evaluate(i,icmp,Liw);
        }
      }
    }
    MatrixRM Kii(N,N), CG(N,N), Wj(N,N), Wii(N,N);
    std::vector<MatrixRM> dG(tau.extent(0), MatrixRM::Zero(N,N));
    for (int i=1; i<tau.extent(0); ++i){
      CG.setZero();
      for (int j=0; j<i; ++j){
        Wj.setZero();
        const auto& terms = wij_geom[geom_index(i,j)];
        for (const auto& term : terms) accumulate_term_matrix_packed(Fpack, F2pack, N*N, term.l, term.coeff, Wj);
        CG.noalias() += Wj * dG[j];
      }
      {
        auto Rx_view = Rx(i,Range::all(),Range::all());
        Eigen::Map<const MatrixRM> Rx_map(Rx_view.dataFirst(), N, N);
        CG += Rx_map;
      }
      Wii.setZero();
      const auto& terms = wij_geom[geom_index(i,i)];
      for (const auto& term : terms) accumulate_term_matrix_packed(Fpack, F2pack, N*N, term.l, term.coeff, Wii);
      Kii.setIdentity();
      Kii -= Wii;
      Eigen::PartialPivLU<MatrixRM> lu(Kii);
      dG[i] = lu.solve(CG);
    }
    Gm.resize(tau.extent(0),N,N); Gm = 0;
    for (int i=0; i<tau.extent(0); ++i)
      for (int i1=0; i1<N; ++i1){
        for (int i2=0; i2<N; ++i2) Gm(i,i1,i2) = dG[i](i1,i2);
        Gm(i,i1,i1) += exp(-Em(i1)*tau(i));
      }
  }
#ifdef _MPI
  void cmp2_large_parallel(Array<double,3>& Gm, const Array<double,3>& Sm, const Array<double,1>& Em,
                           MPI_Comm comm, int comm_rank, int comm_size){
    // Equivalent to cmp2_small, but this is optimized for large blocks NxN with N>params.cmp2_small_cutoff, and also allows inside mpi parallelization.
    //
     const int N = Sm.extent(0);
    Spline1DMany fSig(tau.extent(0),N*N);
    for (int i1=0; i1<N; ++i1)
      for (int i2=0; i2<N; ++i2){
        int icmp=i1*N+i2;
        for (int it=0; it<tau.extent(0); ++it) fSig(it,icmp)=Sm(i1,i2,it);
      }
    fSig.splineIt(tau);
    Spline1DMany fconv(tau.extent(0),N*N);
    for (int it=0; it<tau.extent(0); ++it)
      for (int icmp=0; icmp<N*N; ++icmp) fconv(it,icmp) = 0.0;
    for (int i=0; i<tau.extent(0)-1; ++i){
      double dt = (tau(i+1)-tau(i));
      for (int i1=0; i1<N; ++i1){
        double a = Em(i1)*dt, ema = exp(-a), a2=a*a, a4=a*a*a*a;
        double dt3 = dt*dt*dt;
        TinyVector<double,4> Liw;
        if (fabs(a)>1e-7){
          Liw[0] = dt*(1-ema*(1+a))/a2; Liw[1] = dt*(a-1+ema)/a2;
          Liw[2] = -dt3/a4*((1+a+a2/3.)*ema  + a2/6.-1);
          Liw[3] = -dt3/a4*((1-a+a2/3.) + ema*(a2/6.-1));
        }else{
          Liw[0] = dt*(0.5-a/3.); Liw[1] = dt*(0.5-a/6.);
          Liw[2] = dt3*(-1./24.+a/45.); Liw[3] = dt3*(-1./24.+a*7./360.);
        }
        for (int i2=0; i2<N; ++i2){
          int icmp = i1*N+i2;
          fconv(i+1,icmp) = fconv(i,icmp)*ema + fSig.evaluate(i,icmp,Liw);
        }
      }
    }
    fconv.splineIt(tau);
    std::vector<double> Fpack(tau.extent(0) * N * N), F2pack(tau.extent(0) * N * N);
    for (int l=0; l<tau.extent(0); ++l){
      int off = l * N * N;
      for (int k=0; k<N*N; ++k){ Fpack[off+k]=fconv.F(l,k); F2pack[off+k]=fconv.F2(l,k); }
    }
    Array<double,3> Rx(tau.extent(0),N,N); Rx=0;
    for (int i=0; i<tau.extent(0)-1; ++i){
      for (int i2=0; i2<N; ++i2){
        double dt = (tau(i+1)-tau(i));
        double a = Em(i2)*dt, ema = exp(-a), a2=a*a, a4=a*a*a*a;
        double dt3 = dt*dt*dt;
        TinyVector<double,4> Liw;
        if (fabs(a)>1e-7){
          Liw[0] = dt*(1-ema*(1+a))/a2; Liw[1] = dt*(a-1+ema)/a2;
          Liw[2] = -dt3/a4*((1+a+a2/3.)*ema  + a2/6.-1);
          Liw[3] = -dt3/a4*((1-a+a2/3.) + ema*(a2/6.-1));
        }else{
          Liw[0] = dt*(0.5-a/3.); Liw[1] = dt*(0.5-a/6.);
          Liw[2] = dt3*(-1./24.+a/45.); Liw[3] = dt3*(-1./24.+a*7./360.);
        }
        for (int i1=0; i1<N; ++i1){
          int icmp = i1*N+i2;
          Rx(i+1,i1,i2) = Rx(i,i1,i2)*ema + fconv.evaluate(i,icmp,Liw);
        }
      }
    }
    MatrixRM Kii(N,N), CG(N,N), CG_local(N,N), Wj(N,N), Wii(N,N);
    std::vector<MatrixRM> dG(tau.extent(0), MatrixRM::Zero(N,N));
    for (int i=1; i<tau.extent(0); ++i){
      CG_local.setZero();
      auto [jqs, jqe] = split_range(comm_rank, comm_size, i);
      for (int j=jqs; j<jqe; ++j){
        Wj.setZero();
        const auto& terms = wij_geom[geom_index(i,j)];
        for (const auto& term : terms) accumulate_term_matrix_packed(Fpack, F2pack, N*N, term.l, term.coeff, Wj);
        CG_local.noalias() += Wj * dG[j];
      }
      if (comm_size > 1) MPI_Allreduce(CG_local.data(), CG.data(), N*N, MPI_DOUBLE, MPI_SUM, comm);
      else CG = CG_local;
      {
        auto Rx_view = Rx(i,Range::all(),Range::all());
        Eigen::Map<const MatrixRM> Rx_map(Rx_view.dataFirst(), N, N);
        CG += Rx_map;
      }
      Wii.setZero();
      const auto& terms = wij_geom[geom_index(i,i)];
      for (const auto& term : terms) accumulate_term_matrix_packed(Fpack, F2pack, N*N, term.l, term.coeff, Wii);
      Kii.setIdentity();
      Kii -= Wii;
      Eigen::PartialPivLU<MatrixRM> lu(Kii);
      dG[i] = lu.solve(CG);
    }
    Gm.resize(tau.extent(0),N,N); Gm = 0;
    for (int i=0; i<tau.extent(0); ++i)
      for (int i1=0; i1<N; ++i1){
        for (int i2=0; i2<N; ++i2) Gm(i,i1,i2) = dG[i](i1,i2);
        Gm(i,i1,i1) += exp(-Em(i1)*tau(i));
      }
  }
#endif
  void cmp2(Array<double,3>& Gm, const Array<double,3>& Sm, const Array<double,1>& Em){
    if (Sm.extent(0) <= params.cmp2_small_cutoff) cmp2_small(Gm, Sm, Em);
    else cmp2_large(Gm, Sm, Em);
  }
};

struct SchedulerParams {
  int no_inside_parallel_below = 64;
  int prefer_inside_parallel_at_or_above = 94;
  double block_parallel_rank_fraction = 0.50;
  int min_team_size = 2;
  int max_team_size = 32;
  int practical_team_size_cap = 8;
  int max_dedicated_blocks = 5;
  bool promote_medium_blocks = true;
  bool print_diagnostics = false;
  bool auto_by_mpi_size = true;
  int cmp2_small_cutoff = 4;
};

inline SchedulerParams MakeSchedulerParams(const json& pars, int mpi_size)
{
  SchedulerParams p;
  if (pars.contains("volterra_scheduler_auto") && !pars["volterra_scheduler_auto"].get<bool>())
    p.auto_by_mpi_size = false;
  if (p.auto_by_mpi_size){
    if (mpi_size <= 16){
      p.no_inside_parallel_below = 94;
      p.prefer_inside_parallel_at_or_above = 94;
      p.block_parallel_rank_fraction = 0.50;
      p.practical_team_size_cap = 4;
      p.max_dedicated_blocks = 3;
      p.promote_medium_blocks = false;
      p.cmp2_small_cutoff = 4;
    }else if (mpi_size < 64){
      p.no_inside_parallel_below = 77;
      p.prefer_inside_parallel_at_or_above = 94;
      p.block_parallel_rank_fraction = 0.50;
      p.practical_team_size_cap = 8;
      p.max_dedicated_blocks = 4;
      p.promote_medium_blocks = true;
      p.cmp2_small_cutoff = 4;
    }else{
      p.no_inside_parallel_below = 64;
      p.prefer_inside_parallel_at_or_above = 94;
      p.block_parallel_rank_fraction = 0.50;
      p.practical_team_size_cap = 8;
      p.max_dedicated_blocks = 5;
      p.promote_medium_blocks = true;
      p.cmp2_small_cutoff = 4;
    }
  }
  if (pars.contains("volterra_cmp2_small_cutoff")) p.cmp2_small_cutoff = pars["volterra_cmp2_small_cutoff"].get<int>();
  if (pars.contains("volterra_no_inside_parallel_below")) p.no_inside_parallel_below = pars["volterra_no_inside_parallel_below"].get<int>();
  if (pars.contains("volterra_prefer_inside_parallel_at_or_above")) p.prefer_inside_parallel_at_or_above = pars["volterra_prefer_inside_parallel_at_or_above"].get<int>();
  if (pars.contains("volterra_block_parallel_rank_fraction")) p.block_parallel_rank_fraction = pars["volterra_block_parallel_rank_fraction"].get<double>();
  if (pars.contains("volterra_min_team_size")) p.min_team_size = pars["volterra_min_team_size"].get<int>();
  if (pars.contains("volterra_max_team_size")) p.max_team_size = pars["volterra_max_team_size"].get<int>();
  if (pars.contains("volterra_practical_team_size_cap")) p.practical_team_size_cap = pars["volterra_practical_team_size_cap"].get<int>();
  if (pars.contains("volterra_max_dedicated_blocks")) p.max_dedicated_blocks = pars["volterra_max_dedicated_blocks"].get<int>();
  if (pars.contains("volterra_promote_medium_blocks")) p.promote_medium_blocks = pars["volterra_promote_medium_blocks"].get<bool>();
  if (pars.contains("volterra_print_diagnostics")) p.print_diagnostics = pars["volterra_print_diagnostics"].get<bool>();
  return p;
}

inline std::vector<int> AllocateQuantizedTeamSizes(const std::vector<BlockJob>& jobs,
                                                   int ranks_for_dedicated,
                                                   const SchedulerParams& p)
{
  const int effective_max_team_size = std::min(p.max_team_size, p.practical_team_size_cap);
  std::vector<int> team_sizes(jobs.size(), p.min_team_size);
  int used = p.min_team_size * static_cast<int>(jobs.size());
  if (jobs.empty()) return {};
  if (used > ranks_for_dedicated) return {};
  while (true){
    int best = -1;
    double best_score = -1.0;
    int best_need = 0;
    for (int i=0; i<jobs.size(); ++i){
      int next_size = std::min(team_sizes[i] * 2, effective_max_team_size);
      int need = next_size - team_sizes[i];
      if (need <= 0) continue;
      if (used + need > ranks_for_dedicated) continue;
      double score = jobs[i].cost / team_sizes[i];
      if (score > best_score){ best = i; best_score = score; best_need = need; }
    }
    if (best < 0) break;
    team_sizes[best] += best_need;
    used += best_need;
  }
  return team_sizes;
}

inline void GreedyAssignJobs(std::vector<std::vector<BlockJob>>& assign, const std::vector<BlockJob>& jobs)
{
  std::vector<double> load(assign.size(), 0.0);
  for (const auto& job : jobs){
    int best = 0;
    for (int r=1; r<assign.size(); ++r)
      if (load[r] < load[best]) best = r;
    assign[best].push_back(job);
    load[best] += job.cost;
  }
}

/*
class DysonLaplaceIntegral{
  // Here we solve DysonLaplace transformation with the double integral equation.
  // As usual, we write in matrix form:
  //                                    G = G0 + G0*Sigma*G
  // and we stay in time, hence all multiplications have to be in time and orbital indices.
  // We also take into account that G(t<0)=G0(t<0)=Sigma(t<0) vanishes, hence only the following
  // integration remains:
  //    G(tau) = G0(tau) + Integrate[G0(tau-t2)*Sigma(t2-t1)*G(t1),{t2,t1,tau},{t1,0,tau}]
  // This is equivalent to
  //    G(tau) = G0(tau) + Integrate[ conv(tau-t)*G(t), {t,0,tau}]
  // where
  //    conv(tau) = Integrate[ G0(tau-t)*Sig(t), {t,0,tau}]
  // We evaluate this double integral below, which is actually just twice single integral and takes O(n^2) time.
  // Notice that conv(0)=0, hence G(tau) does not appear on the right-hand side, but
  // only prior time G(tau-delta). The term with G(tau) has G(tau)*conv(0), which vanishes.
  // Also notice that we take 2^r*N points instead of N points for precision. As it turns out, the error
  // even when r=0 is of the order of 10^-6, which probably means that this is not needed.
  //
public:
  int r;
  const Array<double,1>& tau;
  Array<double,1> t;
  Array<double,1> G0, Sig;
  Array<double,1> f, conv, G;
  Array<double,3> fr, convr, Gr;
  Array<intpar,1> ips;
  //
  DysonLaplaceIntegral(const Array<double,1>& _tau_, int _r_) : tau(_tau_), r(_r_)
  {
    if (r==0)
      t.reference(tau);
    else{
      t.resize((1<<r)*(tau.extent(0)-1)+1);
      for (int i=0; i<t.extent(0); ++i)
	t(i) = (tau(tau.extent(0)-1)-tau(0))*i/(t.extent(0)-1.);
    }
    // precompute interpolation coefficients for more dense mesh t when original function is given on mesh tau
    ips.resize(t.extent(0));
    int pos=0;
    for (int i=0; i<t.extent(0); ++i) ips(i) = InterpLeft(t(i),pos,tau);
    Sig.resize(t.extent(0));
    G0.resize(t.extent(0));
    f.resize(t.extent(0));
    conv.resize(t.extent(0));
    G.resize(t.extent(0));
  }
  void cmp1(Array<double,1>& Gm, const Array<double,1>& Sm, double Em){
    //---------------------------------------------------------------------------
    // Scalar Σ(t) → G(t)
    //---------------------------------------------------------------------------
    // This is the simplest case in which Sm,Gm are a single function of time Sm(tau), Gm(tau), as well as Em, defined by G0=exp(-Em*tau).
    // Appropriate for 1D blocks
    // Evaluate G0 on t-mesh
    for (int i=0; i<t.extent(0); ++i) G0(i) = exp(-Em*t(i));
    // Σ(t) on refined t-mesh (either alias or spline‑interpolated)
    if (r==0){
      Sig.reference(Sm);
    }else{
      Spline1D<double> fSig(tau.extent(0));
      for (int i=0; i<tau.extent(0); ++i) fSig[i]=Sm(i);
      fSig.splineIt(tau);
      for (int i=0; i<t.extent(0); ++i) Sig(i) = fSig(ips(i));
    }
    //
    // Here we evaluate the double integral explained above: conv(t) = (G0 * Σ)(t)
    //
    for (int i=0; i<t.extent(0); ++i){
      for (int it=0; it<=i; ++it)
	f(it) = G0(i-it)*Sig(it);
      conv(i) = int_newt_cotes(f, t(0), t(i), i+1);
    }
    // -------- solve Dyson integral for G(t) (scalar) --------
    G=0;
    for (int i=0; i<t.extent(0); ++i){
      for (int it=0; it<=i; ++it)
	f(it) = conv(i-it)*G(it); // notice that the last point is conv(0)*G(i), which is not yet known. However, conv(0)=0, hence G(i) is not needed.
      G(i) = G0(i) + int_newt_cotes(f, t(0), t(i), i+1);
    }
    // down‑sample back to user mesh
    int q=(1<<r);
    Gm.resize(tau.extent(0));
    for (int it=0; it<tau.extent(0); ++it) Gm(it) = G(q*it);
  }
  void cmp2(Array<double,3>& Gm, const Array<double,3>& Sm, const Array<double,1>& Em){
    //---------------------------------------------------------------------------
    // Matrix Σ(t) (N×N) diagonal G0 → G(t) (N×N)      𝑂(N_orb² n²)
    //---------------------------------------------------------------------------
    // Appropriate for blocks of Gm(tau,N,N) and Sm(N,N,tau). We work in basis with diagonal G0.
    //  temporary objects:
    //        Gr(t,N,N), convr(t,N,N), fr(t,N,N)
    //        in addition to G0(t), Sig(t), 
    // 
    // Evaluate G0 on t-mesh
    const int Nt = t.extent(0);
    const int N = Gm.extent(1);
    Gr.resize(Nt,N,N);    Gr=0;   // final result on denser mesh
    convr.resize(Nt,N,N); convr=0;// convolution of G0*Sigma
    // ---------------------------------------------------------------------
    // 1.  Σ interpolation (if needed) + convolution convr = G0 * Σ
    // ---------------------------------------------------------------------
    if (r==0){
      for (int i1=0; i1<N; ++i1){
	for (int i=0; i<t.extent(0); ++i) G0(i) = exp(-Em(i1)*t(i)); // G0(tau_i)(i1,i1)
	for (int i2=0; i2<N; ++i2){
	  // Here we evaluate the double integral explained above
	  for (int i=0; i<Nt; ++i){
	    for (int it=0; it<=i; ++it)
	      f(it) = G0(i-it)*Sm(i1,i2,it); // f = G0 * Sig in matrix form, but notice that G0 is diagonal, hence f(i1,i2)=G0(i1)*Sig(i1,i2)
	    convr(i,i1,i2) = int_newt_cotes(f, t(0), t(i), i+1);
	  }
	}
	Gr(Range::all(),i1,i1) = G0(Range::all()); // we set G = G0
      }
    }else{
      //Sigr.resize(N,N,t.extent(0));
      Spline1D<double> fSig(tau.extent(0));
      for (int i1=0; i1<N; ++i1){
	for (int i=0; i<Nt; ++i) G0(i) = exp(-Em(i1)*t(i)); // G0(tau_i)(i1,i1)
	for (int i2=0; i2<N; ++i2){
	  // first interpolate Sm(i1,i2,tau) on denser mesh t instead of tau, and store in Sig(t)
	  // Σ(i1,i2,t) on refined mesh → Sig
	  for (int k=0; k<tau.extent(0); ++k) fSig[k]=Sm(i1,i2,k);
	  fSig.splineIt(tau);
	  for (int k=0; k<Nt; ++k) Sig(k) = fSig(ips(k));
	  //
	  // Here we evaluate the double integral for convolution convr(tau,i1,i2) explained above
	  for (int i=0; i<Nt; ++i){
	    for (int it=0; it<=i; ++it)
	      f(it) = G0(i-it)*Sig(it); // f = G0 * Sig in matrix form, but notice that G0 is diagonal, hence f(i1,i2)=G0(i1)*Sig(i1,i2)
	    convr(i,i1,i2) = int_newt_cotes(f, t(0), t(i), i+1);
	  }
	}
	Gr(Range::all(),i1,i1) = G0(Range::all()); // we set G = G0 
      }
    }
    // ---------------------------------------------------------------------
    // 2.  Solve Dyson equation: G = G0 + conv * G   (time‑ordered)
    // ---------------------------------------------------------------------
    fr.resize(Nt,N,N);
    for (int i=0; i<t.extent(0); ++i){
      // fr(it ≤ i) = convr(i-it) * Gr(it)
      for (int it=0; it<=i; ++it){ // only need Gr(it<i,:,:) to set up fr(it<=i,:,:)
	auto  fr_view   = fr(it,Range::all(),Range::all());
	//tiny::matmul(fr_view, convr(i-it,Range::all(),Range::all()), Gr(it,Range::all(),Range::all())); // notice that the last point is conv(0)*G(i), which is not yet known. However, conv(0)=0, hence G(i) is not needed.
	tiny::dot(fr_view, convr(i-it,Range::all(),Range::all()), Gr(it,Range::all(),Range::all())); // notice that the last point is conv(0)*G(i), which is not yet known. However, conv(0)=0, hence G(i) is not needed.
      }
      for (int i1=0; i1<N; ++i1)
	for (int i2=0; i2<N; ++i2)
	  Gr(i,i1,i2) += int_newt_cotes(fr(Range::all(),i1,i2), t(0), t(i), i+1); // use fr(it<=i,:,:) to update Gr(i,:,:)
    }
    // ---------------------------------------------------------------------
    // 3.  Down‑sample to user mesh tau
    // ---------------------------------------------------------------------
    int q=(1<<r);
    Gm.resize(tau.size(),N,N);
    for (int it=0; it<tau.extent(0); ++it){
      if (fabs(t(q*it)-tau(it))>1e-6) cout<<"ERROR t["<<(q*it)<<"]="<<t(q*it)<<" tau["<<it<<"]="<<tau(it)<<endl;
      Gm(it,Range::all(),Range::all()) = Gr(q*it,Range::all(),Range::all());
    }
  }
};
 */

class DysonLaplaceIntegral{
  // Here we solve DysonLaplace transformation with the double integral equation.
  // As usual, we write in matrix form:
  //                                    G = G0 + G0*Sigma*G
  // and we stay in time, hence all multiplications have to be in time and orbital indices.
  // We also take into account that G(t<0)=G0(t<0)=Sigma(t<0) vanishes, hence only the following
  // integration remains:
  //    G(tau) = G0(tau) + Integrate[G0(tau-t2)*Sigma(t2-t1)*G(t1),{t2,t1,tau},{t1,0,tau}]
  // This is equivalent to
  //    G(tau) = G0(tau) + Integrate[ conv(tau-t)*G(t), {t,0,tau}]
  // where
  //    conv(tau) = Integrate[ G0(tau-t)*Sig(t), {t,0,tau}]
  // We evaluate this double integral below, which is actually just twice single integral and takes O(n^2) time.
  // Notice that conv(0)=0, hence G(tau) does not appear on the right-hand side, but
  // only prior time G(tau-delta). The term with G(tau) has G(tau)*conv(0), which vanishes.
  // Also notice that we take 2^r*N points instead of N points for precision. As it turns out, the error
  // even when r=0 is of the order of 10^-6, which probably means that this is not needed.
  //
  //  This is called:  Volterra equation of the second kind, 
  //
public:
  int r;
  const Array<double,1>& tau;
  Array<double,1> t;
  Array<double,1> G0, Sig;
  Array<double,1> f, conv, G;
  Array<double,3> fr, convr, Gr;
  Array<intpar,1> ips, ipq;
  //
  DysonLaplaceIntegral(const Array<double,1>& _tau_, int _r_) : tau(_tau_), r(_r_)
  {
    if (r==0){
      t.reference(tau);
      double dh=tau(1)-tau(0);
      for (int i=1; i<tau.size()-1; ++i)
	if (fabs(tau(i+1)-tau(i)-dh)>1e-7) cout<<"ERROR in DysonLaplaceIntegral, we are trying to use non-equidistant mesh on equidistant algorithm"<<endl;
    }else{
      int Nt = tau.extent(0);
      t.resize((1<<r)*(Nt-1)+1);
      for (int i=0; i<t.extent(0); ++i)
	t(i) = (tau(Nt-1)-tau(0))*i/(t.extent(0)-1.);
    }
    // precompute interpolation coefficients for more dense mesh t when original function is given on mesh tau
    ips.resize(t.extent(0));
    int pos=0;
    for (int i=0; i<t.extent(0); ++i) ips(i) = InterpLeft(t(i),pos,tau);
    ipq.resize(tau.extent(0));
    pos=0;
    for (int i=0; i<tau.extent(0); ++i) ipq(i) = InterpLeft(tau(i),pos,t);
    Sig.resize(t.extent(0));
    G0.resize(t.extent(0));
    f.resize(t.extent(0));
    conv.resize(t.extent(0));
    G.resize(t.extent(0));
  }
  void cmp1(Array<double,1>& Gm, const Array<double,1>& Sm, double Em){
    //---------------------------------------------------------------------------
    // Scalar Σ(t) → G(t)
    //---------------------------------------------------------------------------
    // This is the simplest case in which Sm,Gm are a single function of time Sm(tau), Gm(tau), as well as Em, defined by G0=exp(-Em*tau).
    // Appropriate for 1D blocks
    // Evaluate G0 on t-mesh
    static int iwhich=0;
    iwhich += 1;
    
    for (int i=0; i<t.extent(0); ++i) G0(i) = exp(-Em*t(i));
    // Σ(t) on refined t-mesh (either alias or spline‑interpolated)
    //cout<<"t="<<t<<endl;
    //cout<<"G0="<<G0<<endl;
    if (r==0){
      Sig.reference(Sm);
    }else{
      Spline1D<double> fSig(tau.extent(0));
      for (int i=0; i<tau.extent(0); ++i) fSig[i]=Sm(i);
      fSig.splineIt(tau);
      for (int i=0; i<t.extent(0); ++i) Sig(i) = fSig(ips(i));
      //ofstream outd("debug_laplace."+to_string(iwhich)+".dat"); outd.precision(14);
      //for (int it=0; it<t.extent(0); ++it)
      //  outd << t(it) << " " << Sig(it) << endl;
    }
    //cout<<"Sig="<<Sig<<endl;
    //
    // Here we evaluate the double integral explained above: conv(t) = (G0 * Σ)(t)
    //
    for (int i=0; i<t.extent(0); ++i){
      for (int it=0; it<=i; ++it)
	f(it) = G0(i-it)*Sig(it);
      conv(i) = int_newt_cotes(f, t(0), t(i), i+1);
    }
    // -------- solve Dyson integral for G(t) (scalar) --------
    G=0;
    for (int i=0; i<t.extent(0); ++i){
      for (int it=0; it<=i; ++it)
	f(it) = conv(i-it)*G(it); // notice that the last point is conv(0)*G(i), which is not yet known. However, conv(0)=0, hence G(i) is not needed.
      G(i) = G0(i) + int_newt_cotes(f, t(0), t(i), i+1);
    }
    // down‑sample back to user mesh
    //int q=(1<<r);
    //for (int it=0; it<tau.extent(0); ++it) Gm(it) = G(q*it);
    Gm.resize(tau.extent(0));
    {
      Spline1D<double> fGm(t.extent(0));
      for (int i=0; i<t.extent(0); ++i) fGm[i]=G(i);
      fGm.splineIt(t);
      for (int i=0; i<tau.extent(0); ++i) Gm(i) = fGm(ipq(i));
    }
  }
  void cmp2(Array<double,3>& Gm, const Array<double,3>& Sm, const Array<double,1>& Em){
    //---------------------------------------------------------------------------
    // Matrix Σ(t) (N×N) diagonal G0 → G(t) (N×N)      𝑂(N_orb² n²)
    //---------------------------------------------------------------------------
    // Appropriate for blocks of
    //                  Gm(tau,N,N) and Sm(N,N,tau).
    // We work in basis with diagonal G0(N,tau).
    //  temporary objects:
    //        Gr(t,N,N), convr(t,N,N), fr(t,N,N)
    //        in addition to G0(t), Sig(t), 
    //
     // Evaluate G0 on t-mesh
    const int Nt = t.extent(0);
    const int N = Sm.extent(0);
    if (N!=Gm.extent(2) or N!=Sm.extent(0) or N!=Sm.extent(1)) cout<<"ERROR dimensions of Gm="<<Gm.shape()<<" and Sm="<<Sm.shape()<<" are incompatible in DysonLaplaceIntegral.cmp2"<<endl;
    Gr.resize(Nt,N,N);    Gr=0;   // final result on denser mesh
    convr.resize(Nt,N,N); convr=0;// convolution of G0*Sigma
    // ---------------------------------------------------------------------
    // 1.  Σ interpolation (if needed) + convolution convr = G0 * Σ
    // ---------------------------------------------------------------------
    if (r==0){
      for (int i1=0; i1<N; ++i1){
	for (int i=0; i<t.extent(0); ++i) G0(i) = exp(-Em(i1)*t(i)); // G0(tau_i)(i1,i1)
	for (int i2=0; i2<N; ++i2){
	  // Here we evaluate the double integral explained above
	  for (int i=0; i<Nt; ++i){
	    for (int it=0; it<=i; ++it)
	      f(it) = G0(i-it)*Sm(i1,i2,it); // f = G0 * Sig in matrix form, but notice that G0 is diagonal, hence f(i1,i2)=G0(i1)*Sig(i1,i2)
	    convr(i,i1,i2) = int_newt_cotes(f, t(0), t(i), i+1);
	  }
	}
	Gr(Range::all(),i1,i1) = G0(Range::all()); // we set G = G0
      }
    }else{
      Spline1DMany fSig(tau.extent(0), N*N);
      // first interpolate Sm(i1,i2,tau) on uniform mesh t instead of tau, and store in Sig(t)
      // Σ(i1,i2,t) on refined mesh → Sig
      for (int it=0; it<tau.extent(0); ++it)
	for (int i1=0; i1<N; ++i1)
	  for (int i2=0; i2<N; ++i2)
	    fSig(it,i1*N+i2)=Sm(i1,i2,it); // This is self-energy on non-uniform mesh tau(it)
      fSig.splineIt(tau);                  // splining all functions at once.
      //convr(0, Range::all(), Range::all()) = 0.0;
      for (int i1=0; i1<N; ++i1){
	for (int i=0; i<Nt; ++i) G0(i) = exp(-Em(i1)*t(i)); // G0(tau_i)(i1,i1)
	for (int i2=0; i2<N; ++i2){
	  for (int it=0; it<Nt; ++it) Sig(it) = fSig(ips(it),i1*N+i2); // Here we have Sig on uniform mesh t(it).
	  //
	  // Here we evaluate the double integral for convolution convr(tau,i1,i2) explained above
	  for (int i=0; i<Nt; ++i){
	    for (int it=0; it<=i; ++it)
	      f(it) = G0(i-it)*Sig(it); // f = G0 * Sig in matrix form, but notice that G0 is diagonal, hence f(i1,i2)=G0(i1)*Sig(i1,i2)
	    convr(i,i1,i2) = int_newt_cotes(f, t(0), t(i), i+1);
	  }
	}
	Gr(Range::all(),i1,i1) = G0(Range::all()); // we set G = G0 
      }
    }
    // ---------------------------------------------------------------------
    // 2.  Solve Dyson equation: G = G0 + conv * G   (time‑ordered)
    // ---------------------------------------------------------------------
    fr.resize(Nt,N,N);
    for (int i=0; i<t.extent(0); ++i){
      // fr(it ≤ i) = convr(i-it) * Gr(it)
      for (int it=0; it<=i; ++it){ // only need Gr(it<i,:,:) to set up fr(it<=i,:,:)
	auto  fr_view   = fr(it,Range::all(),Range::all());
	//tiny::matmul(fr_view, convr(i-it,Range::all(),Range::all()), Gr(it,Range::all(),Range::all())); // notice that the last point is conv(0)*G(i), which is not yet known. However, conv(0)=0, hence G(i) is not needed.
	tiny::dot(fr_view, convr(i-it,Range::all(),Range::all()), Gr(it,Range::all(),Range::all())); // notice that the last point is conv(0)*G(i), which is not yet known. However, conv(0)=0, hence G(i) is not needed.
      }
      for (int i1=0; i1<N; ++i1)
	for (int i2=0; i2<N; ++i2)
	  Gr(i,i1,i2) += int_newt_cotes(fr(Range::all(),i1,i2), t(0), t(i), i+1); // use fr(it<=i,:,:) to update Gr(i,:,:)
    }
    // ---------------------------------------------------------------------
    // 3.  Down‑sample to user mesh tau
    // ---------------------------------------------------------------------
    Gm.resize(tau.extent(0),N,N);
    {
      Spline1DMany fGm(t.extent(0),N*N);
      for (int i=0; i<t.extent(0); ++i)
	for (int i1=0; i1<N; ++i1)
	  for (int i2=0; i2<N; ++i2)
	    fGm(i,i1*N+i2) = Gr(i,i1,i2); // G on equidistant mesh is interpolated.
      fGm.splineIt(t);                    // with splines on uniform mesh t(it)
      for (int it=0; it<tau.extent(0); ++it) // Now interpolated on non-uniform mesh.
	for (int i1=0; i1<N; ++i1)
	  for (int i2=0; i2<N; ++i2)
	    Gm(it,i1,i2) = fGm(ipq(it),i1*N+i2);
    }
  }
};


void FindLargestEigenvalues(Array<double,2>& gm,
			    const Array<double,2>& Gm, const Array<double,1>& tau, const Cix& cix, int itt, MPIcls& mpi)
{
  gm.resize(cix.Nblock,tau.size());
  gm=0;
  std::vector<BlockJob> jobs; jobs.reserve(cix.Nblock);
  { // this part just figures out which blocks require computation. Some are equivalent and can be skept.
    for (int ib=0; ib<cix.Nblock; ++ib){
      int N = cix.Block[ib].size(); // block size
      jobs.push_back({ib, N, double(N)*N*N}); // {iblock,size,cost}
    }
  }
  std::vector<std::vector<BlockJob>> assign(mpi.size);
  if (mpi.size==1)
    assign[0]=jobs;
  else
    mpi.Assign_jobs(assign, jobs); // Here we actually assign jobs.
  ///////////////////////////////////////////////////////
  
  fsec t_prep_HS=zero_t, t_diag1=zero_t;
  Array<double,1> work(3*cix.Nmaxsize); // for exact diagonalization
  for (int ij=0; ij<assign[mpi.rank].size(); ij++){
    int ib = assign[mpi.rank][ij].id;
    int N = cix.Block[ib].size(); // block size
    if (N==1){ // no need to perform exact diagonalization, because the block is one dimensional
      int _ii_ = cix.Sgm_symmetry[ib](0,0);
      gm(ib,Range::all()) = Gm(_ii_-1,Range::all());
    }else{
      Array<double,2> HS(N,N); 
      Array<double,1> E(N);  // stores eigenvalues
      for (int it=0; it<tau.size(); ++it){
	auto t0 = high_resolution_clock::now();
	HS=0;
	for (int i=0; i<N; ++i){
	  for (int j=0; j<N; ++j){
	    int _ii_ = cix.Sgm_symmetry[ib](i,j);
	    if (_ii_==0) continue;
	    if (_ii_>0)
	      HS(i,j) = Gm(_ii_-1,it);
	    else
	      HS(i,j) = -Gm(-_ii_-1,it);
	  }
	}
	auto t1 = high_resolution_clock::now();
	xsyev(HS.extent(0), HS.data(), HS.extent(1), E.data(), work.data(), work.size(), "N");
	gm(ib,it) = max(E);
	auto t2 = high_resolution_clock::now();
	t_prep_HS += t1-t0;
	t_diag1 += t2-t1;
      }
    }
  }
  mpi.AllReduceSum_simple_array(gm);
  if (mpi.rank==mpi.master){
    fmt::print(*mpi.out,"t_pre={:.4f} t_diag={:.4f}\n",t_prep_HS.count(), t_diag1.count());
    ofstream out("Gpp_blocks."+to_string(itt+1));
    for (int it=0; it<tau.size(); ++it){
      out<<tau(it)<<" ";
      for (int ib=0; ib<cix.Nblock; ++ib) out<<gm(ib,it)<<"  ";
      out<<endl;
    }
  }
}

  struct retVals{
    double lambda0;
    int    followPeak;
    retVals(double l=0.0, int p=-1) : lambda0(l), followPeak(p){};
    operator std::pair<double,int>() const noexcept {return {lambda0, followPeak};}
    retVals(const std::pair<double,int>& pr): lambda0(pr.first), followPeak(pr.second) {}
    retVals& operator=(const std::pair<double,int>& pr){lambda0    = pr.first; followPeak = pr.second; return *this;}
  };

retVals GetLambda0_in_sector(int ib,
			     const Array<double,1>& tau, const Array<double,1>& tau_dh, const Array<double,2>& Sm, const Cix& cix,
			     double peak_posT, Array<double,1>& work, ostream& clog)
{
  retVals r(100,-1);
  int N = cix.Block[ib].size(); // block size
  if (N==1){ // no need to perform exact diagonalization, because the block is one dimensional
    double E = cix.Ham[ib](0,0);
    // this should not need abs, because negative symmetry is allowed only for off-diagonal components.
    int ii = abs(cix.Sgm_symmetry[ib](0,0))-1; // note that cix.Sgm_symmetry starts with 1, hence need to subtract 1.
    double S0 =0;
    for (int it=0; it<tau.extent(0); ++it) S0 += Sm(ii,it)*tau_dh(it);
    double lmbda = E-S0-peak_posT;
    r.lambda0=lmbda;
    r.followPeak=ib+1;
    clog<<"lambda["<<ib+1<<"]="<<lmbda<<" Energy["<<ib+1<<"]="<<E<<" S0="<<S0<<endl;
  }else{
    int counter=0;
    std::map<int,int> indx;
    for (int i=0; i<N; ++i){
      for (int j=0; j<N; ++j){
	int ii = abs(cix.Sgm_symmetry[ib](i,j))-1;
	if (ii>=0){
	  auto [it, inserted] = indx.try_emplace(ii, counter);
	  if (inserted) counter += 1;
	}
      }
    }
    //std::cout << "indx: ";
    //for (const auto& [key, val] : indx) std::cout << "(" << key << "->" << val << ") ";
    //std::cout << "\n";
    // integral of self-energy, which is self-energy at zero frequency
    Array<double,1> S0_pp(counter); S0_pp=0;
    for (const auto& [ii, icount] : indx){
      for (int it=0; it<tau.extent(0); ++it) S0_pp(icount) += Sm(ii,it)*tau_dh(it);
    }
    //cout<<" ingerals are S0_pp"<< S0_pp << " or " << _S0_pp_ << endl;
    Array<double,2> HS(N,N); HS=0;
    for (int i=0; i<N; ++i){
      for (int j=0; j<N; ++j){
	int _ii_ = cix.Sgm_symmetry[ib](i,j);
	int ii = abs(_ii_)-1;
	if (ii>=0){
	  HS(i,j) = S0_pp(indx.at(ii))*sign(_ii_); // We set (i,j) matrix element from precomputed integral.
	  //cout<<"S0["<<ib+1<<",("<<i<<","<<j<<")]="<< HS(i,j) << endl;
	}
      }
    }
    //cout<<"-Sig(om=0)["<<ib<<"]="<<HS<<endl;
    
    Array<double,1> E(N);  // stores eigenvalues
    HS *= -1;
    HS += cix.Ham[ib];
    //cout<<"H+Sig(om=0)["<<ib<<"]="<<HS<<endl;
    //t_prep_HS += t1-t0;
    xsyev(HS.extent(0), HS.data(), HS.extent(1), E.data(), work.data(), work.size(), "N");
    double lmbda = min(E)-peak_posT;
    clog<<"lambda["<<ib+1<<"]="<<lmbda<<" Energy["<<ib+1<<"]-S0="<<E<<endl;
    //cout<<"lambda["<<ib+1<<"]="<<lmbda<<" Energy["<<ib+1<<"]-S0="<<E<<endl;
    r.lambda0=lmbda;
    r.followPeak=ib+1;
  }
  return r;
}


auto DeterminePseudoPropagators(Array<double,2>& Gm,
				const Array<double,1>& tau, const Array<double,1>& tau_dh, const Array<double,2>& Sm,
				const Cix& cix, double beta, double lambda0_fixed, int followPeak_fixed, int itt, MPIcls& mpi,
				const json& pars, double peak_pos=0)
{
  retVals r(100,-1);
  
  bool DEBUG=false;
  bool PrintEigens=false;
  double t_diag1=0, t_diag2=0, t_oneD=0, t_sg_trans=0, t_dl=0, t_g_trans=0, t_set_g=0;
  //map<int,fsec> t_dim;
  //unordered_map<int,int> t_deg;

  ostream& clog = *mpi.out;
  if (tau_dh.extent(0) != tau.extent(0)){ clog<<"ERROR in DeterminePseudoPropagators: tau_dh.size="<<tau_dh.extent(0)<<" while tau.size="<<tau.extent(0); exit(1);}
  SchedulerParams sched = MakeSchedulerParams(pars, mpi.size);
  VolterraSolver::RuntimeParams vparams;
  vparams.cmp2_small_cutoff = sched.cmp2_small_cutoff;
  VolterraSolver volterra(tau, vparams);
  ////////////////////////////////////////////////////////////////////////////
  // distributes blocks across multiple processors according to
  // dimension of each block. The optimized Volterra solver is dominated by
  // matrix-matrix multiplies for larger blocks, hence the cost is closer to N^3.
  ////////////////////////////////////////////////////////////////////////////
  std::vector<BlockJob> jobs; jobs.reserve(cix.Nblock);
  Array<int,1> deg(cix.Nblock); deg=0;
  { // this part just figures out which blocks require computation. Some are equivalent and can be skept.
    Array<int,1> G_done(cix.Npseudo_unique); G_done=0;
    for (int ib=0; ib<cix.Nblock; ++ib){
      int N = cix.Block[ib].size(); // block size
      bool need_to_compute=false;   // check if at least one Gm components is not yet computed.
      for (int i=0; i<N*N; ++i){
	int ii = abs(cix.Sgm_symmetry[ib](i/N,i%N))-1; // use symmetry to not compute the same componente multiple times.
	//if (ii>=0) cout<<"  debug ib="<<ib+1<<" (i,j)="<<(i/N)<<","<<(i%N)<<" ii="<<ii<<" G_done[ii]="<<G_done(ii)<<endl;
	if (ii>=0 and G_done(ii)==0){
	  need_to_compute=true;
	  break;
	}
      }
      if (need_to_compute){
	jobs.push_back({ib, N, double(N)*N*N}); // {iblock,size,cost}
	for (int i=0; i<N; ++i){
	  for (int j=0; j<N; ++j){
	    int ii = abs(cix.Sgm_symmetry[ib](i,j))-1; // use symmetry to not compute the same componente multiple times.
	    if (ii>=0) G_done(ii)=1;
	  }
	}
	deg(ib) += 1;
	//cout<<"ib="<<ib+1<<" deg["<<ib+1<<"]+=1 : "<<deg(ib)<<endl;
      }else{
	int ib_equiv=-1;
	for (int i=0; i<N*N; ++i){
	  int ii = abs(cix.Sgm_symmetry[ib](i/N,i%N))-1;
	  if (ii>=0){
	    ib_equiv = abs(cix.Sgm_indx[ii][0][0]);
	    break;
	  }
	}
	deg(ib_equiv) += 1;
	//cout<<"ib="<<ib+1<<" deg["<<ib_equiv+1<<"]+=1 : "<<deg(ib_equiv)<<endl;
      }
    }
  }

  //clog<<"deg="<<deg<<endl;
  std::vector<std::vector<BlockJob>> assign(mpi.size);
  if (mpi.size==1)
    assign[0]=jobs;
  else
    mpi.Assign_jobs(assign, jobs); // Here we actually assign jobs.
  ///////////////////////////////////////////////////////
  int itlast = tau.size()-1;

  //cout<<"jobs:";
  //for (int l=0; l<jobs.size(); ++l)
  //  cout<<"jobs["<<l<<"]="<< jobs[l] <<endl;
  clog.precision(14);
  Array<double,1> work(3*cix.Nmaxsize); // for exact diagonalization
  if (lambda0_fixed==0){  // we are updating lambda0 at this iteration
    double peak_posT = peak_pos/beta; // if peak is not centered exactly at zero

    if (followPeak_fixed>0){ // lambda0 is taken from one sector, hence all processors compute the single sector where lambda0 is fixed.
      int ib = followPeak_fixed-1;
      r = GetLambda0_in_sector(ib,
			       tau, tau_dh, Sm, cix, peak_posT, work, clog);
    }else{// We need to compute lambda0 in all sectors.
      r.lambda0=100.; r.followPeak=-1;

      for (int ij=0; ij<assign[mpi.rank].size(); ij++){
	int ib = assign[mpi.rank][ij].id;
	auto t0 = get_cpu_time();
	retVals r_ib = GetLambda0_in_sector(ib,
					    tau, tau_dh, Sm, cix, peak_posT, work, clog);
	if (r_ib.lambda0<r.lambda0) r = r_ib;
	auto t1 = get_cpu_time();
	t_diag1 += t1-t0;
      }
      clog<<"  on rank="<<mpi.rank<<" lmbda0 ="<<r.lambda0<<" followPeak="<<r.followPeak<<endl;
      ///////////////////////////////////////////////////////
      if (mpi.size>1) // We have smallest lambda0 in this core, but we need to find the smallest across processors.
	r = mpi.FindSmallest_lambda0(r);
      ///////////////////////////////////////////////////////
    }
    clog<<"   "<<itt<<") lmbda0="<<r.lambda0<<" followPeak="<<r.followPeak<<endl;
  }else{
    clog<<"using previously-fixed lambda0="<<lambda0_fixed<<endl;
    r.lambda0 = lambda0_fixed;
    r.followPeak = 0;
  }

  ///////////////////////////////////////////////////////
  Array<double,1> Prob_eig(cix.Natst); Prob_eig=0;
  Gm.resize(cix.Npseudo_unique,tau.size());
  Gm = 0;
  Array<double,1> E;
  ////// This is helper function solve_block
  auto solve_block = [&](int ib, bool write_output
#ifdef _MPI
			 , bool use_inside_parallel, MPI_Comm team_comm, int team_rank, int team_size
#endif
			 )
  {
    int N = cix.Block[ib].size();
    if (N==1){
      auto t0 = get_cpu_time();
      double Em = cix.Ham[ib](0,0) - r.lambda0;
      int ii = abs(cix.Sgm_symmetry[ib](0,0))-1;
      Array<double,1> Gm_view;
      volterra.cmp1(Gm_view, Sm(ii,Range::all()), Em);
      if (write_output){
	Gm(ii,Range::all()) = Gm_view;
	int i_pseudo = cix.Block[ib](0);
	Prob_eig(i_pseudo-1) = Gm_view(itlast);
      }
      t_oneD += get_cpu_time()-t0;
      return;
    }

    auto t0 = get_cpu_time();
    Array<double,2> H(cix.Ham[ib].copy()); // will contain eigenvectors of Hamiltonian
    E.resize(N);
    xsyev(H.extent(0), H.data(), H.extent(1), E.data(), work.data(), work.size(), "V"); // diagonalizing Hamiltonian
    auto t1 = get_cpu_time();
    t_diag2 += t1-t0;

    Array<double,3> Sgm(N,N,tau.size());
    Array<double,3> Gmm(tau.size(),N,N);
    Array<double,2> Sg1(N,N), Sg2(N,N);
    for (int it=0; it<tau.size(); ++it){
      Sg1=0; // NxN matris of self-energy
      for (int i=0; i<N; ++i){
	for (int j=0; j<N; ++j){
	  int _ii_ = cix.Sgm_symmetry[ib](i,j);
	  int ii = abs(_ii_)-1;
	  if (ii>=0) Sg1(i,j) = Sm(ii,it)*sign(_ii_);
	}
      }
      tiny::dot(Sg2,H,Sg1);             // Sg2 = H*Sg1
      tiny::dot(Sg1,Sg2,H,false,true);  // Sg1 = H*Sg1*H^+
      Sgm(Range::all(),Range::all(),it) = Sg1; // Sgm(:,:,it) <- H*Sg1(it)*H^+
    }
    auto t2 = get_cpu_time();
    t_sg_trans += t2-t1;

    E -= r.lambda0;
#ifdef _MPI
    // If blocks are large, we parallelize inside volterra.cmp2.
    if (use_inside_parallel && team_size>1) volterra.cmp2_large_parallel(Gmm, Sgm, E, team_comm, team_rank, team_size);
    else volterra.cmp2(Gmm, Sgm, E);
#else
    volterra.cmp2(Gmm, Sgm, E);
#endif
    auto t3 = get_cpu_time();
    t_dl += t3-t2;

    if (write_output){
      for (int i=0; i<N; ++i){
	int i_pseudo = cix.Block[ib](i);
	Prob_eig(i_pseudo-1) = Gmm(itlast,i,i);
      }
    }
    // Now transforming G back to non-diagonal basis, i.e., computational basis
    for (int it=0; it<tau.size(); ++it){
      tiny::dot(Sg1,H,Gmm(it,Range::all(),Range::all()),true,false); // Sg1 = H^+*G
      tiny::dot(Sg2,Sg1,H);   // Sg2 = H^+*G*H
      Gmm(it,Range::all(),Range::all()) = Sg2; // Gmm <- H^+*G*H
    }
    auto t4 = get_cpu_time();
    t_g_trans += t4-t3;

    if (write_output){
      for (int i=0; i<N; ++i){
	for (int j=0; j<N; ++j){
	  int _ii_ = cix.Sgm_symmetry[ib](i,j);
	  if (_ii_>0 && bl::all(cix.Sgm_indx[_ii_-1][0] == TinyVector<int,3>(ib,i,j))){
	    int ii = _ii_-1;
	    Gm(ii,Range::all()) = Gmm(Range::all(),i,j);
	  }
	}
      }
    }
    t_set_g += get_cpu_time()-t4;
  };/// end of helper function solve_block
  ////////////////////////////////////////////////
  // parallelization is optimized depending of how many cores are available and how large are blocks.
  std::vector<BlockJob> dedicated_jobs, medium_candidates, small_jobs;
  for (const auto& job : jobs){
    if (job.n >= sched.prefer_inside_parallel_at_or_above) dedicated_jobs.push_back(job);
    else if (sched.promote_medium_blocks && job.n >= sched.no_inside_parallel_below) medium_candidates.push_back(job);
    else small_jobs.push_back(job);
  }
  auto by_cost = [](const BlockJob& a, const BlockJob& b){ return a.cost > b.cost; };
  std::sort(dedicated_jobs.begin(), dedicated_jobs.end(), by_cost);
  std::sort(medium_candidates.begin(), medium_candidates.end(), by_cost);
  if ((int)dedicated_jobs.size() > sched.max_dedicated_blocks) dedicated_jobs.resize(sched.max_dedicated_blocks);
  if (sched.promote_medium_blocks){
    for (const auto& job : medium_candidates){
      if ((int)dedicated_jobs.size() >= sched.max_dedicated_blocks) break;
      dedicated_jobs.push_back(job);
    }
  }
  std::sort(dedicated_jobs.begin(), dedicated_jobs.end(), by_cost);
  std::set<int> dedicated_ids;
  for (const auto& job : dedicated_jobs) dedicated_ids.insert(job.id);
  std::vector<BlockJob> all_small_jobs;
  for (const auto& job : jobs)
    if (!dedicated_ids.count(job.id)) all_small_jobs.push_back(job);

#ifdef _MPI
  int requested_block_workers = all_small_jobs.empty() ? 0 : std::max(1, int(std::ceil(sched.block_parallel_rank_fraction * mpi.size)));
  int ranks_for_dedicated = std::max(0, mpi.size - requested_block_workers);
  std::vector<int> team_sizes;
  if (!dedicated_jobs.empty())
    team_sizes = AllocateQuantizedTeamSizes(dedicated_jobs, ranks_for_dedicated, sched);
  if (team_sizes.empty()){
    dedicated_jobs.clear();
    dedicated_ids.clear();
    all_small_jobs = jobs;
  }
  int total_dedicated_ranks = 0;
  for (int ts : team_sizes) total_dedicated_ranks += ts;
  int n_block_workers = std::max(0, mpi.size - total_dedicated_ranks);
  if (sched.print_diagnostics && mpi.rank==mpi.master){
    clog << "Volterra scheduler size=" << mpi.size
	 << " cmp2_small_cutoff=" << sched.cmp2_small_cutoff
	 << " dedicated_blocks=" << dedicated_jobs.size()
	 << " dedicated_ranks=" << total_dedicated_ranks
	 << " block_workers=" << n_block_workers << endl;
    for (int i=0; i<dedicated_jobs.size(); ++i)
      clog << "  dedicated ib=" << dedicated_jobs[i].id+1 << " N=" << dedicated_jobs[i].n << " team_size=" << team_sizes[i] << endl;
  }

  int rank0 = 0;
  for (int ided=0; ided<dedicated_jobs.size(); ++ided){
    int team_size = team_sizes[ided];
    bool in_team = (mpi.rank >= rank0 && mpi.rank < rank0 + team_size);
    int color = in_team ? (ided+1) : MPI_UNDEFINED;
    MPI_Comm team_comm = MPI_COMM_NULL;
    MPI_Comm_split(MPI_COMM_WORLD, color, mpi.rank, &team_comm);
    if (in_team){
      int team_rank = -1, team_nproc = 0;
      MPI_Comm_rank(team_comm, &team_rank);
      MPI_Comm_size(team_comm, &team_nproc);
      solve_block(dedicated_jobs[ided].id, team_rank==0, true, team_comm, team_rank, team_nproc);
      MPI_Comm_free(&team_comm);
    }
    rank0 += team_size;
  }

  if (n_block_workers>0){
    std::vector<std::vector<BlockJob>> small_assign(n_block_workers);
    GreedyAssignJobs(small_assign, all_small_jobs);
    if (mpi.rank >= total_dedicated_ranks){
      int small_worker_index = mpi.rank - total_dedicated_ranks;
      if (small_worker_index < 0 || small_worker_index >= small_assign.size())
	throw std::runtime_error("Invalid small_worker_index in DeterminePseudoPropagators");
      for (const auto& job : small_assign[small_worker_index])
	solve_block(job.id, true, false, MPI_COMM_NULL, 0, 1);
    }
  }

  if (mpi.size>1){
    mpi.AllReduceSum_simple_array(Gm);
    mpi.AllReduceSum_simple_array(Prob_eig);
  }
#else
  for (const auto& job : jobs) solve_block(job.id, true);
#endif
  
  if (DEBUG and mpi.rank==mpi.master){
    ofstream out2("PPG.dat"); out2.precision(10);
    for (int it=0; it<tau.size(); ++it){
      out2<<tau(it)<<" ";
      for (int ii=0; ii<cix.Npseudo_unique; ++ii) out2<< Gm(ii,it)<<" ";
      out2<<endl;
    }
  }
  fmt::print(clog,"   # t_diag1={:.4f},t_oneD={:.4f} t_diag2={:.4f}s ", t_diag1, t_oneD, t_diag2);
  fmt::print(clog,"t_sg_trans={:.4f}s t_volterra={:.4f}s t_g_trans={:.4f}s t_set_g={:4f}s\n", t_sg_trans, t_dl, t_g_trans, t_set_g);
  //for (auto&& [dim, t] : t_dim | std::views::all) 
  //  fmt::print(clog,"t[{}]={:4f}s ", dim, t.count()/t_deg[dim]);
  //fmt::print(clog,"\n");

  
  if (mpi.rank==mpi.master){
    double small=5e-7;
    // We skipped several states when computing Probability, because they are equivalent to some other state.
    // Now we need to set all these empty values.
    for (int ib=0; ib<cix.Nblock; ++ib){
      int i=0;
      int ii = abs(cix.Sgm_symmetry[ib](i,i))-1;
      if (ii>=0){
	int _ib_ = cix.Sgm_indx[ii][0][0];
	if (_ib_!=ib){
	  for (int i=0; i<cix.Block[ib].size(); ++i){
	    if (Prob_eig(cix.Block[ib](i)-1)==0){
	      Prob_eig(cix.Block[ib](i)-1) = Prob_eig(cix.Block[_ib_](i)-1);
	    }else{
	      clog << "WARN : Found equivalence ib="<<ib+1<<" and ib="<<_ib_+1<<"  P_orig["<<ib+1<<","<<i<<"]="<<Prob_eig(cix.Block[ib](i)-1)<<endl;
	    }
	  }
	}
      }
    }
    double sumQ = sum(Prob_eig);
    Prob_eig *= 1/sumQ;
    /*
    double sm1=0, sm2=0;
    for (int ib=0; ib<cix.Nblock; ++ib){
      int N = cix.Block[ib].size();
      for (int i=0; i<N; ++i){
	int ind = cix.Block[ib](i)-1;
	sm1 += Prob_eig(ind)*deg(ib);
	sm2 += Prob_eig(ind);
	clog<<"sm1["<<ind<<"]="<<sm1<<" sm2["<<ind<<"]="<<sm2<<" deb["<<ib<<"]="<<deg(ib)<<endl;
      }
    }
    */
    /*
    double sumQ = 0;
    for (int ib=0; ib<cix.Nblock; ++ib){
      int N = cix.Block[ib].size();
      for (int i=0; i<N; ++i) sumQ += Prob_eig(cix.Block[ib](i)-1)*deg(ib);
    }
    */
    //double sumQ = sum(Prob_eig);
    clog<<"   "<<itt<<") Q = "<<sumQ<<endl;
    ofstream out("Probability.dat");
    for (int ib=0; ib<cix.Nblock; ++ib){
      int N = cix.Block[ib].size();
      double Ptot=0;
      for (int i=0; i<N; ++i) Ptot += Prob_eig(cix.Block[ib](i)-1);
      if (Ptot>small){
	int nmax=N;
	for (; nmax>0; --nmax){
	  if (abs(Prob_eig(cix.Block[ib](nmax-1)-1)) > 0.2*small) break;
	}
	int Nocc = cix.iN(cix.Block[ib](0)-1);
	fmt::print(clog,"   {:d}) Prob_eig[{:d},N={:d}]={:10.7f} [",itt,ib+1,Nocc,Ptot);
	for (int i=0; i<nmax; ++i)
	  fmt::print(clog,"{:10.7f},",Prob_eig(cix.Block[ib](i)-1));
	clog<<"]"<<endl;
      }
      for (int i=0; i<N; ++i) 
	out<<ib+1<<" "<<i<<" "<<Prob_eig(cix.Block[ib](i)-1)<<endl;
    }
  }
  return r;
}    
auto DeterminProbabilitiesI(const Array<double,1>& tau, const Array<double,2>& Gm, const Cix& cix,
			    double beta, double mu, double lambda0, int followPeak, int itt, MPIcls& mpi, double small=5e-7)
{
  struct retVals{
    double nd, Q, lambda0, Fimp, Epot, TrSigmaG;
    int followPeak;
    Array<double,1> density_matrix, Prob;
  };
  retVals r;
  r.followPeak = followPeak;
  ostream& clog = *mpi.out;
  //
  Array<double,1> Probability(cix.Natst), nalpha(cix.Nbath);
  Probability=0; nalpha = 0;
  //
  int itlast = tau.size()-1;
  double sumQ=0, sum_nd=0, dQ=0;
  r.Epot=0;
  for (int ib=0; ib<cix.Nblock; ++ib){
    int N = cix.Block[ib].size();
    double sumQ1=0;
    for (int i=0; i<N; ++i){
      int i_pseudo = cix.Block[ib](i);
      int ii = cix.Sgm_symmetry[ib](i,i)-1; // note that cix.Sgm_symmetry starts with 1, hence need to subtract 1.
      // diagonal components should not have negative ii. This is only allowed for off-diagonal
      if (ii>=0){
	double Nocc = cix.iN(i_pseudo-1);
	dQ = Gm(ii,itlast);
	sumQ += dQ;
	sumQ1 += dQ;
	Probability(i_pseudo-1) += dQ;
	sum_nd += dQ*Nocc;
	//for (int ibth=0; ibth<cix.Nbath; ++ibth){
	//  int fj = cix.Fi(2*ibth+1,i_pseudo); // This computes psi_{ibth}|i_pseudo> and is this is nonzero, we know that
	//  if (fj!=0)                      // psi^+_{ibth}psi_{ibth}|i_pseudo>=|i_pseudo>
	//    nalpha(ibth) += dQ;
	//}
      }
      for (int ibth=0; ibth<cix.Nbath; ++ibth){
	int j_pseudo = cix.Fi(2*ibth+1,i_pseudo);      // <j_pseudo|F_{ibth}|i_pseudo>
	if (j_pseudo!=0){
	  int k_pseudo = cix.Fi(2*ibth,abs(j_pseudo)); // <k_pseudo|F^+_{ibth}|j_pseudo><j_pseudo|F_{ibth}|i_pseudo>
	  if (k_pseudo!=0){
	    int k = cix.inBlock(abs(k_pseudo));
	    int _ii_ = cix.Sgm_symmetry[ib](k,i);
	    int ii = abs(_ii_)-1;
	    if (ii>=0){
	      double dQm = Gm(ii,itlast);
	      nalpha(ibth) += sign(_ii_*j_pseudo*k_pseudo)*dQm; // dQ[k_pseudo,i_pseudo] * <k_pseudo|F^+_{ibth}|j_pseudo><j_pseudo|F_{ibth}|i_pseudo>
	      //cout<<"ibth="<<ibth<<" i_pseudo="<<i_pseudo<<" j_pseudo="<<j_pseudo<<" k_pseudo="<<k_pseudo<<" ii="<<ii<<" (i,k)=("<<i<<","<<k<<")"<<" dnalpha="<<sign(_ii_*j_pseudo*k_pseudo)*dQm(ii)<<endl;
	    }
	  }
	}
      }
      for (int j=0; j<N; ++j){
	int _ii_ = cix.Sgm_symmetry[ib](i,j);
	int ii = abs(_ii_)-1;
	if (ii>=0) r.Epot += sign(_ii_)*Gm(ii,itlast)*cix.Ham[ib](j,i);
      }
    }
    //clog<<"sumQ["<<ib+1<<"]="<<sumQ1<<endl;
  }
  r.Prob.resize(cix.Npseudo_unique);
  for (int ii=0; ii<cix.Npseudo_unique; ++ii) r.Prob(ii) = Gm(ii,itlast)/sumQ;
  
  Probability *= 1/sumQ;
  r.nd = sum_nd/sumQ;
  nalpha *= 1/sumQ;
  r.Q = sumQ;
  r.lambda0 = lambda0;
  r.Epot *= 1/sumQ;
  r.Epot += r.nd*mu; // this is because r.Epot = Tr(H^{atom}_m P_m), but H^{{atom} was previously modified to H^{atom}-mu*N
  r.Fimp = lambda0-log(r.Q)/beta;

  //cout<<"Ebath="<<cix.Ebath<<" bath_sigind="<<cix.bath_sigind<<" nalpha="<<nalpha<<endl;
  double Epot_singlep=0;
  for (int b=0; b<cix.Nbath; b++){
    int iib = cix.bath_sigind0(b);
    int ib_sg = 1;
    if (iib<0){ iib=-iib; ib_sg=-1;}
    Epot_singlep += (cix.Ebath(iib)+mu)*ib_sg*nalpha(b); // previously we subtracted Ebath-=mu, now we add back, like in Epot
  }
  r.TrSigmaG = r.Epot - Epot_singlep; // This is now Tr(H^{atom}_m P_m) - Tr(E^{imp}_a n_a) and is 1/2Tr(Sigma*G) because we subtract single-particle contribution to potential energy

  r.density_matrix.resize(cix.Nbath_unique_unique); r.density_matrix=0;
  Array<int,1> deg(cix.Nbath_unique_unique); deg=0;
  for (int iorb=0; iorb<cix.Nbath; ++iorb){
    int iib0 = cix.bath_sigind0(iorb);
    double ib_sg = (iib0 < 0 ? -1 : 1);
    int iib   = std::abs(iib0);
    r.density_matrix(iib) += nalpha(iorb)*ib_sg;
    deg(iib) += 1;
  }
  r.density_matrix /= deg;
  
    
  if (mpi.rank==mpi.master){
    clog<<"   "<<itt<<") nd= "<<r.nd<<endl;
    clog<<"   "<<itt<<") Q=  "<<r.Q<<endl;
    for (int ib=0; ib<cix.Nblock; ++ib){
      int N = cix.Block[ib].size();
      double Ptot=0;
      for (int i=0; i<N; ++i) Ptot += Probability(cix.Block[ib](i)-1);
      if (Ptot>small){
	int nmax=N;
	for (; nmax>0; --nmax){
	  if (abs(Probability(cix.Block[ib](nmax-1)-1)) > 0.2*small) break;
	}
	int Nocc = cix.iN(cix.Block[ib](0)-1);
	//clog<<"   "<<itt<<") Prob["<<ib+1<<",N="<<Nocc<<"]=[";
	fmt::print(clog,"   {:d}) Prob[{:d},N={:d}]={:10.7f} [",itt,ib+1,Nocc,Ptot);
	for (int i=0; i<nmax; ++i)
	  fmt::print(clog,"{:10.7f},",Probability(cix.Block[ib](i)-1));
	clog<<"]"<<endl;
      }
    }
    clog<<"   "<<itt<<") Fimp = "<<r.Fimp<<endl;
    clog<<"   "<<itt<<") Epot = "<<r.Epot<<endl;
    clog<<"   "<<itt<<") TrSigmaG = "<<r.TrSigmaG<<endl;
    clog<<"   "<<itt<<") Epot_singlep = "<<Epot_singlep<<endl;
    //clog<<"   occupancies: "<<endl<<"   ";
    clog.setf(ios::fixed);
    clog<<"   "<<itt<<") bath_density_matrix="<< r.density_matrix <<endl;
    clog.unsetf(ios::fixed);
  }
  return r;
}

class DeltaKey{
public:
  int two_N_max_order, Nbath, two_N_max_order_sqr;
  Array<int,1> to_p_3; // This is array for 3^n
  //
  DeltaKey(int _N_max_order_, int _Nbath_) :
    two_N_max_order(2*_N_max_order_),  Nbath(_Nbath_), two_N_max_order_sqr(two_N_max_order*two_N_max_order)
  {}
  int DKey(int ib, int is, int ie) const{
    //return Nbath*(2*N_max_order*ie + is) + ib;
    return is + two_N_max_order*(ie+two_N_max_order*ib);
  }
  auto FromDKey(int ddkey) const{
    struct retVals{
      int ib,is,ie;
      retVals(int _ib_, int _is_, int _ie_): ib(_ib_),is(_is_),ie(_ie_){};
    };
    int is   = ddkey % two_N_max_order;
    int rest = ddkey / two_N_max_order;
    int ie   = rest  % two_N_max_order;
    int ib   = rest  / two_N_max_order;
    return retVals(ib,is,ie);
  }
  int IntKey(int is, int ie){
    return two_N_max_order*ie + is;
  }
  auto FromIntKey(int ddkey){
    struct retVals{
      int is, ie;
      retVals(int _is_, int _ie_) : is(_is_),ie(_ie_){};
    };
    int is = ddkey % two_N_max_order;
    int ie = ddkey / two_N_max_order;
    return retVals(is,ie);
  }
  auto IntwDKey(int ib, int is, int ie) const{
    struct retVals{
      int pkey, dkey;
      retVals(int _pkey_, int _dkey_): pkey(_pkey_), dkey(_dkey_){};
    };
    int pkey = two_N_max_order*ie + is;
    return retVals(pkey, pkey+two_N_max_order_sqr*ib);
  }
};

/////////////////////////////////////////////////////////////////////
// Below is code for smoothing imaginary time data by attenuation  //
/////////////////////////////////////////////////////////////////////
inline int wrap_mod(int i, int N) {
    int r = i % N;
    return (r < 0) ? (r + N) : r;
}
blitz::Array<double,1> gaussian_kernel(int L, double sigma_pts)
{
  // Build a normalized 1D Gaussian kernel of odd length L with width/sigma in *points*
  if (L < 3) L = 3;
  if ((L & 1) == 0) ++L; // make odd
  blitz::Array<double,1> k(L);
  const int mid = L/2;
  double s2 = sigma_pts * sigma_pts;
  double sum = 0.0;
  for (int i = 0; i < L; ++i) {
    double x = double(i - mid);
    k(i) = std::exp(-0.5 * (x*x) / s2);
    sum += k(i);
  }
  k /= sum;
  return k;
}

// Periodic 1D convolution: y = x (*) k    (k length = L odd)
void convolve_periodic(blitz::Array<double,1>& y,
		       const blitz::Array<double,1>& x, const blitz::Array<double,1>& k)
{
  const int N = x.extent(0);
  const int L = k.extent(0);
  const int mid = L/2;
  if (y.extent(0) != N) y.resize(N);
  // y(j) = sum_m k(m) * x(j + (m - mid) mod N)
  for (int j = 0; j < N; ++j) {
    double acc = 0.0;
    for (int m = 0; m < L; ++m) {
      int idx = wrap_mod(j + (m - mid), N);
      acc += k(m) * x(idx);
    }
    y(j) = acc;
  }
}
// Optional: reflected-boundary convolution (if you prefer non-periodic)
int reflect_index(int i, int N) {
  // Reflect into [0, N-1]
  while (i < 0 || i >= N) {
    if (i < 0) i = -1 - i;
    else       i = 2*N - 1 - i;
  }
  return i;
}
void convolve_reflect(blitz::Array<double,1>& y,
		      const blitz::Array<double,1>& x, const blitz::Array<double,1>& k)
{
  const int N = x.extent(0);
  const int L = k.extent(0);
  const int mid = L/2;
  if (y.extent(0) != N) y.resize(N);
  for (int j = 0; j < N; ++j) {
    double acc = 0.0;
    for (int m = 0; m < L; ++m) {
      int idx = reflect_index(j + (m - mid), N);
      acc += k(m) * x(idx);
    }
    y(j) = acc;
  }
}

class SmoothCorrectOrders{
  // Smoothly reduce orders such that each order is smaller in magnitude |G[n]| > T * |G[n-1]| (on average).
  //
  blitz::Array<double,1> ker;
  blitz::Array<double,1> u, us, tmp; // Work buffers (1D along tau)
public:
  //
  // sigma_pts: Gaussian width in grid points for smoothing (5–9 is a good start)
  explicit SmoothCorrectOrders(int Ntau, double sigma_pts = 7.0) : u(Ntau), us(Ntau), tmp(Ntau)
  {
    // Build Gaussian kernel once
    const int L = std::max(3, int(std::ceil(6.0 * sigma_pts)) | 1); // ~6σ, odd
    ker.resize(L);
    ker = gaussian_kernel(L, sigma_pts);
  }
  blitz::Array<double,1> cmp(const blitz::Array<double,1>& G_n, const blitz::Array<double,1>& G_nm, double T = 2.0, double kappa = 1.0, bool periodic = true, bool resmooth_mask = true)
  {
    // Returns attenuation a(τ) in [0,1]
    // G_n(Ntau) :  order n that we are correcting
    // G_nm(Ntau):  order n-1 that we are comparing to
    // T:    threshold factor (2.0 means "more than twice previous order")
    // kappa: how hard to push back once above T (1.0 ~ bring smoothed ratio back to T)
    // periodic: use periodic smoothing in tau (natural for β-interval); set false to use reflect
    //
    const int Ntau   = G_n.extent(0);
    if (Ntau!=u.extent(0)){
      std::cerr<<"ERROR : mismatch with Ntau="<<Ntau<<" previously used Ntau="<<u.extent(0)<<endl;
      assert(Ntau == u.extent(0));
    }
    const double logT = std::log(std::max(1.0, T));
    const double eps  = 1e-14;
    // 1) log-ratio u = log(|G[n]|+eps) - log(|G[n-1]|+eps)
    for (int t = 0; t < Ntau; ++t) {
      double gnm1 = std::abs(G_nm(t)) + eps;
      double gn   = std::abs(G_n(t)) + eps;
      u(t) = std::log(gn) - std::log(gnm1);
    }
    // 2) smooth u along tau
    if (periodic) convolve_periodic(us, u, ker);
    else          convolve_reflect (us, u, ker);
    // 3) excess above threshold and attenuation a = exp(-kappa * max(0, u_s - logT))
    Array<double,1> a(Ntau);
    for (int t = 0; t < Ntau; ++t) {
      double e = us(t) - logT;
      if (e < 0.0) e = 0.0;
      a(t) = std::exp(-kappa * e); // in (0,1]
    }
    // 4) optionally smooth the mask a to be extra gentle
    if (resmooth_mask) {
      if (periodic) convolve_periodic(tmp, a, ker);
      else          convolve_reflect(tmp, a, ker);
      // clamp to [0,1]
      for (int t = 0; t < Ntau; ++t)
	a(t) = std::min(1.0, std::max(0.0, tmp(t)));
    }
    return a;
  }
};
