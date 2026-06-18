#include <iostream>
#include <algorithm>
#include <vector>
#include <deque>
#include <blitz/array.h>
#include "basics.h"
#include "wnormalization.h"
#include "polynomial.h"

using namespace std;
using namespace blitz;


inline double wfactorial(int n){
  double r=1;
  for (int i=2; i<=n; i++) r = r*i;
  return r;
}

inline double ipower(double base, int exp)
{
  switch (exp){
  case 0 : return 1.; break;
  case 1 : return base; break;
  case 2 : return base*base; break;
  case 3 : return base*base*base; break;
  default :
    if (exp<0){
      exp = -exp;
      base = 1./base;
    }
    double result = 1;
    while (exp){
      if (exp & 1) result *= base;
      exp >>= 1;
      base *= base;
    }
    return result;
  }
}

//typedef __float128 dreal;

std::ostream& operator<<(std::ostream& stream, const vector<double>& f)
{
  stream<<"[";
  for (int i=0; i<f.size(); i++) stream<<std::setw(6)<<f[i]<<", ";
  stream<<"]";
  return stream;
}
double IterativeNestedLoop(int depth, const vector<double>& E_r, double E_i)
{
  int max = E_r.size();
  // Initialize the slots to hold the current iteration value for each depth
  //int* slots = (int*)alloca(sizeof(int) * depth);
  Array<int,1> slots(depth);
  slots = 0;

  double ordern = 0;
  
  int index = 0;
  while (true){
    // TODO: Your inner loop code goes here. You can inspect the values in slots
    // slots[0] == most inside loop
    // slots[depth-1] == first loop
    double prod = 1;
    for (int j=0; j<depth; ++j) prod *= (E_r[slots(j)]-E_i);
    // prod at depth=4 : ((E_r[ir]-E_i)*(E_r[jr]-E_i)*(E_r[kr]-E_i)*(E_r[lr]-E_i));
    ordern += 1/prod;
    //cout << slots << " prod=" << prod << " ordern="<< ordern << endl;
    
    // Increment
    slots(0) += 1;
    // Carry
    while (slots(index) == max){
      // Overflow, we're done
      if (index == depth - 1) return ordern;
      //slots(index) = 0; // use this is you want to restart the loop with 0, rather than i_k >= i_{k-1}
      for (int j=0; j<=index; j++) slots(j) = slots(index+1)+1; // This is here because we want i_k >= i_{k-1}
      index += 1;
      slots(index) += 1;
    }
    index = 0;
  }
  return ordern;
}

void GetSortedEnergy(const Array<double,1>& _Et_, vector<double>& Ei)
{ // The order of energies is irrelevant, hence we can sort them for convenience
  int N = _Et_.size();
  //Ei.resize(N);
  for (int i=0; i<N; i++) Ei[i] = _Et_(i);
  std::sort(Ei.begin(),Ei.end());
}
int GetDegIndex(const vector<double>& Ei, vector<vector<int> >& ind_deg)
{// From sorted energies, we can find the degeneracy of each subset of energies.
  // For example, if input energies are [0,1,0,1] the sorted energies Ei=[0,0,1,1] and ind_deg=[[0,1],[2,3]] for degeneries.
  int max_deg=0; // maximum degeneracy for what is needed for small splitting
  int deg=0;
  deque<int> degs;
  for (int i=1; i<Ei.size(); i++){
    //cout<<" E["<<i<<"]="<<Ei[i]<<" E["<<i-1<<"]="<<Ei[i-1]<<endl;
    if ( fabs(Ei[i]-Ei[i-1])<1e-12 ) deg++;
    else{
      degs.push_back(deg+1);
      deg=0;
    }
    max_deg = max(deg,max_deg);
    //cout << " i=" << i << " deg=" << deg << " max_deg=" << max_deg << endl;
  }
  degs.push_back(deg+1);
  //cout << "degs="<< degs << endl;
    
  //vector<vector<int>> ind_deg(degs.size()); 
  ind_deg.resize(degs.size()); // index of degenerate sets of Energies
  int ii=0;
  for (int i=0; i<degs.size(); i++){
    ind_deg[i].resize(degs[i]);
    for (int j=0; j<degs[i]; j++){
      ind_deg[i][j]=ii;
      ii++;
    }
  }
  /*
    cout<<"ind_deg=[";
    for (int i=0; i<ind_deg.size(); ++i){
    cout<<"(";
    for (int j=0; j<ind_deg[i].size(); ++j) cout<<ind_deg[i][j]<<",";
    cout<<"),";
    }
    cout<<"]"<<endl;
  */
  return max_deg;
}

double AnDegInt(const Array<double,1>& Et, double beta)
{ // This routine works for arbitrary set of energies with degeneracies.
  // We are evaluating the integral
  //    \int_0^b dt_n \int_0^{t_n} dt_{n-1} ... \int_0^t_3 dt_2 \int_0^{t_2} dt_1 e^{-E_1 t_1 - E_2(t_2-t_1) - E_3(t_3-t_2) - .. E_n(t_n-t_{n-1}) - E_0(beta-t_n)}
  // which can also be written as:
  //    (\Prod_i \int_0^b dt_i ) theta(0<t_1<t_2<t_3<...t_n<b) \Prod_j e^{-E_j(t_j-t_{j-1})}
  // When there is no degeneracy, the integral is equal to:
  //    \sum_i e^{-b E_i}/(\Prod_{j\ne i} (E_j-E_i) )
  // With degeneracy, it is more complicated, and for a single degeneracy at E_i, we get a contribution
  //    e^{-b E_i}/(\Prod_{E_j\ne E_i} (E_j-E_i) ) * (b - \sum_{Ej!=E_i} 1/(Ej-E_i));
  // For the higher order degeneracy, we worked out the result below.
  vector<double> Ei(Et.extent(0));
  GetSortedEnergy(Et, Ei);
  vector<vector<int> > ind_deg; // index of degenerate sets of Energies
  GetDegIndex(Ei, ind_deg);
  double res=0;
  int N = Ei.size();
  for (int i=0; i<ind_deg.size(); i++){ // over degenerate sets
    double E_i = Ei[ind_deg[i][0]];     // energy of the current degenerate set
    int Ni = ind_deg[i].size();         // size of the degenerate set
    int Nr = N-Ni;                      // number of the rest of the states, which are not in this degenerate set
    vector<double> E_r(Nr);             // Energies of the rest of the states
    int ii=0;
    for (int j=0; j<ind_deg.size(); j++){ // over all the rest degenerate groups
      if (i==j) continue;                 // not this one
      for (int k=0; k<ind_deg[j].size(); k++){
	int ir = ind_deg[j][k];
	E_r[ii] = Ei[ ind_deg[j][k] ];    // save the rest of the states in an array
	ii++;
      }
    }
    //cout << "Now E_i="<< E_i <<" E_r="<< E_r << endl;
    // produces the leading part, which is sufficient when there is no degeneracy
    //   e^{-b E_i}/(\Prod_{E_j\ne E_i} (E_j-E_i) )
    double prod=1;
    for (int ir=0; ir<Nr; ir++) prod *= (E_r[ir]-E_i);
    double first_part = exp(-beta*E_i)/prod;
    //cout << "Prod="<< prod << " first part="<< first_part << endl;
    // second part starts with volume of the degenerate subspace : b^(Ni-1)/(Ni-1)!
    double second_part = ipower(beta,Ni-1)/wfactorial(Ni-1);
    if (E_r.size()>0){
      int m_1_power_n=-1; // (-1)^n
      for (int n=1; n<Ni; n++){
	double ordern = IterativeNestedLoop(n, E_r, E_i);
	second_part += m_1_power_n * ipower(beta,Ni-1-n)/wfactorial(Ni-1-n) * ordern;
	m_1_power_n *= -1;
	// Ni==2 : second_part = b - order1;
	// Ni==3 : second_part = b^2/2! - b*order1 + order2;
	// Ni==4 : second_part = b^3/3! - b^2/2!*order1 + beta*order2 -order3;
	// order1 = \sum_{Ej!=Ei} 1/(Ej-Ei)
	// order2 = \sum_{Ej!=Ei,Ek!=Ei,k>=j} 1/((Ej-Ei)*(Ek-Ei))
	// order3 = \sum_{Ej!=Ei,Ek!=Ei,El!=Ei,k>=j,l>=k} 1/((Ej-Ei)*(Ek-Ei)*(El-Ei))
      }
    }
    //cout << "second part="<< second_part << endl;
    res += first_part * second_part;
  }
  return res;
}
/*
class meassureWeight{
  //The operator() returns the value of a meassuring diagram, which is a function that we know is properly normalized to unity.
  double beta;
  Array<double,1> Et;
public:
  double integral0;
public:
  meassureWeight(){};
  void Init(const Array<double,1>& _Et_, double _beta_){
    beta=_beta_;
    Et.resize(_Et_.size()); Et=_Et_.copy();
    integral0 = AnDegInt(Et,beta);
    //cout<<"Initializing meassureWeight with beta="<<beta<<" n="<<Et.extent(0)<<endl;
  }
  meassureWeight(const Array<double,1>& _Et_, double _beta_) : Et(_Et_), beta(_beta_)
  {
    integral0 = AnDegInt(Et,beta);
  };
  double operator()(const Array<double,1>& times) const
  {
    //cout<<"Evaluating meassureWeight with n="<<Et.size()<<endl;
    int N = Et.size();
    double dexp = Et(N-1)*(beta+times(0)-times(N-1));
    //cout<<" dexp="<<dexp<<" Et[-1]="<<Et(N-1)<<" dt="<<(beta+times(0)-times(N-1))<<endl;
    for (int it=0; it<N-1; it++){
      dexp += Et(it)*(times(it+1)-times(it));
      //cout<<" dexp="<<dexp<<" Et["<<it<<"]="<<Et(it)<<" dt="<<(times(it+1)-times(it))<<endl;
    }
    //cout<<"res="<< (exp(-dexp)/integral0) << endl;
    return exp(-dexp)/integral0;
  }
  
};
*/
class CumulativeInt1D{
  // This class calculates exact integral of picewise function g(x) at an arbitrary point.
  // It implements Int[g(x), {x,-beta,t}] for arbitrary t when g(x) is picewise constant
public:
  double beta, dh;
  double Iall;           // Iall = Int[g(x), {x,0,beta=dh*Nbin}]
  bl::Array<double,1> g; // function g(x) to calculate
  bl::Array<double,1> Ig; // Ig[i] = Int[g(x), {x,0,dh*i}]
  CumulativeInt1D(const bl::Array<double,1>& _g_, double _beta_) : beta(_beta_), g(_g_.extent(0)), Ig(_g_.extent(0))
  {
    int Nbin = _g_.extent(0); // discrete function g(x)
    dh = beta/Nbin;
    for (int i=0; i<Nbin; ++i) g(i) = _g_(i);
    double cum = 0;
    Ig = 0;  // Ig[i] = Int[g(x), {x,0,dh*i}]
    for (int i=0; i<Nbin; ++i){
      Ig(i) = cum;
      cum += g(i)*dh;
    }
    Iall = cum; // Iall = Int[g(x), {x,0,beta=dh*Nbin}]
    //cout << "Ih(N-1)=" << Ih(Nbin-1) << endl;
  }
  double operator()(double t) // Int[g(x), {x,-beta,t}]
  {
    int Nbin = g.extent(0);
    if (t==0){
      return Iall;
    }else if (t < 0){
      int ii = (1.+t/beta)*Nbin;
      if (ii>=Nbin) return Iall; // t/beta is essentially zero, hence ii=Nbin
      //cout << "ii=" << ii << " Nbin=" << Nbin << " Ih=" << Ih(ii) << " h(ii)="<< h(ii)
      return Ig(ii) + g(ii)*(beta+t-ii*dh);
    }else{
      int ii = t/beta*Nbin;
      if (ii >= Nbin) return 0;
      return Iall + Ig(ii) + g(ii)*(t-ii*dh);
    }
  }
};

bl::Array<double,1> NumericInt1DType1(const bl::Array<double,2>& gx, const bl::Array<double,3>& gx_off, double beta, int i_first_time=1, int r=3)
{
  int Nloops = gx.extent(0); // This is N=Nloops-1 in the above equation
  int Nbin = gx.extent(1);   // number of bins
  int Noff = Nloops-1- i_first_time; // Number of off-diagonal functions Noff above
  int Nr = ipower(2,r)+1;                // How many points is used for romberg integration below. It needs to be exact for polynomial of correct degree.
  double dh = beta/Nbin;
  bl::Array<double,1> Ps(Nr), Pt(Nr);
  bl::Array<double,1> t1(Nr);
  bl::Array<double,1> Integ(Nloops);
  Integ=0;
  for (int isk=i_first_time; isk<Nloops; isk++){ // variable j in the above equation
    double dsum = 0.0;
    for (int j=0; j<Nbin; ++j){                        // bins corresponding to variable j
      // equidistant mesh of Nr point betweem [j*dh,....,(j+1)*dh].
      for (int l=0; l<Nr; ++l) t1(l) = j*dh + 1e-6 + (1-2e-6)*dh*l/(Nr-1.); // more precise mesh (for romberg integration) of variable r_j, i.e., the last integral in the above equation
      Pt(bl::Range::all()) = 1.0;                         // Will contain the product of functions \product_{i\ne j and i=1}^{N} .....
      
      for (int ik=i_first_time; ik<Nloops; ++ik){     // Here we compute the product i\ne j, and we use i==ik and j==isk
	if (ik!=isk){
	  CumulativeInt1D I_ci(gx_off(isk,ik,bl::Range::all()), beta); // This class calculates exact integral of picewise function gx_off[j,i](u)*u at an arbitrary point. See discussion on I_ci above.
	  Ps(bl::Range::all()) = 0;                      // Will contain the product of functions \product_{i\ne j and i=1}^{N}
	  for (int i1=0; i1<j+1; i1++){
	    for (int l=0; l<Nr; ++l){             // Here we have t2<t1
	      double t2s = min(    i1*dh, t1(l)); // t2 always smaller than t1 
	      double t2e = min((i1+1)*dh, t1(l)); // t2 always remains smaller than t1
	      double t1t2s = t1(l)-t2s; // t1-t2 at the start, but taking into account t1>t2 hence this is positive
	      double t1t2e = t1(l)-t2e; // t1-t2 at the end,   but taking into account t1>t2 hence this is positive
	      if (t1t2s < 1e-8) t1t2s = 1e-8;
	      if (t1t2e < 1e-8) t1t2e = 1e-8;
	      // Int[gx_off(t2),{t2,-beta,beta-(t1-(i1+1)*dh)}]-Int[gx_off(t2),{t2,-beta,beta-(t1-i1*dh)}]=
	      // t2 is between [i1*dh,i2*dh] but t2<t1
	      // Ps(t1=j*dh+l*dsmall) = gx[i1]*Int[gx_off(t2),{t2,beta-(t1-i1*dh),beta-(t1-(i1+1)*dh)}]
	      // Ps(t1=j*dh+l*dsmall) = gx[i1]*Int[gx_off(t2-t1),{t2,beta+i1*dh,beta+(i+1)*dh}]
	      // always make sure that t1>t2
	      Ps(l) += ( I_ci(beta-t1t2e)-I_ci(beta-t1t2s) ) * gx(ik,i1);
	    }
	  }
	  for (int i1=j; i1<Nbin; i1++){
	    for (int l=0; l<Nr; ++l){                // Here we have t2>t1
	      double t2s = max(      i1*dh, t1(l) ); // t2 always larger than t1
	      double t2e = max( (i1+1.)*dh, t1(l) ); // t2 remains larger than t1
	      if (t2e > beta) t2e = beta;
	      double t2st1 = t2s-t1(l); // t2-t1 at the start, is positive
	      double t2et1 = t2e-t1(l); // t2-t1 at the end, it is still positive
	      if (t2st1 < 1e-8) t2st1 = 1e-8;
	      if (t2et1 < 1e-8) t2et1 = 1e-8;
	      // t2 is between [i1*dh,i2*dh] but t2>t1
	      // Ps(t1=j*dh+l*dsmall) = gx[i1]*Int[gx_off(t2),{t2,t2s-t1,t2e-t1}] = gx[i1]*Int[gx_off(t2-t1),{t2,t2s,t2e}]
	      Ps(l) += ( I_ci(t2et1)-I_ci(t2st1) ) * gx(ik,i1);
	    }
	  }
	  Pt(bl::Range::all()) *= Ps(bl::Range::all()); // Product over all Ps functions. There are Noff products.
	  //cout << "isk=" << isk << " ik=" << ik << " Ps=" << Ps << " Pt=" << Pt << endl;
	}
      }
      double Pt_Int = romberg(Pt, dh);
      dsum += Pt_Int * gx(isk,j);
    }
    Integ(isk) = fabs(dsum);
  }
  return Integ;
}

class TmeassureWeight{
public:
  int Ntime, Nbin;
  double beta;
  double integral0;
  Array<double,1> Et;
  bool self_consistent;
  bl::Array<double,1> gt;
  bl::Array<double,1> t_hist;
  //
  TmeassureWeight(){};
  TmeassureWeight(const Array<double,1>& _Et_, double _beta_, int _Nbin_): Et(_Et_), beta(_beta_), Nbin(_Nbin_), self_consistent(false)
  {
    Ntime = Et.extent(0);
    t_hist.resize(Nbin); t_hist=0;
    integral0 = AnDegInt(Et,beta);
  }
  void Init(const Array<double,1>& _Et_, double _beta_, int _Nbin_){
    self_consistent=false;
    Ntime = _Et_.extent(0);
    beta = _beta_;
    Nbin = _Nbin_;
    Et.resize(Ntime);
    for (int i=0; i<Ntime; ++i) Et(i) = _Et_(i);
    t_hist.resize(Nbin); t_hist=0;
    integral0 = AnDegInt(Et,beta);
  }
  double trst(double x){ return fabs(x); }
  double f1_time(double t) const{
    const double small=1e-16;
    int ip = (t/beta)*Nbin;
    if (ip>=Nbin or ip<0) return 0;
    double res = gt(ip);
    return res > small ? res : small;
  }
  double operator()(const bl::Array<double,1>& times) const
  {
    double PQ_new = 1.0;
    if (! self_consistent){
      double dexp=0;
      for (int it=0; it < Ntime; it++) dexp += Et(it)*(times(it+1)-times(it));
      PQ_new = exp(-dexp)/integral0;
    }else{
      for (int it=0; it < Ntime; it++) PQ_new *= f1_time( times(it+1)-times(it) );
    }
    return PQ_new;
  }
  void Normalize_histogram(double dnrm){
    for (int ip=0; ip<Nbin; ++ip) t_hist(ip) *= dnrm;
  }
  void Recompute(bool SaveData){
    gt.resize(Nbin);
    self_consistent=true;
    gt=0;
    {
      gt(0) = trst( 0.5*(t_hist(0)+t_hist(1)) );
      for (int i=1; i<Nbin-1; i++)
	gt(i) = trst( (t_hist(i-1)+t_hist(i)+t_hist(i+1))/3. );
      gt(Nbin-1) = trst( 0.5*(t_hist(Nbin-2)+t_hist(Nbin-1)) );
    }
    double dh=beta/Nbin;
    // first level of normalization
    double intg1D = sum( gt(bl::Range::all()) ) * dh;
    gt(bl::Range::all()) *= 1./intg1D;
    // second level of correct normalization
    Array<double,2> Kt(Nbin,1);
    for (int i=0; i<Nbin; ++i) Kt(i,0) = gt(i);
    Array<double,2> Rt;
    GiveR(Rt, Kt, gt, dh); // Rx = g*g  -> order=1 and Ntime=2
    int norder=Ntime/2;
    for (int iord=1; iord<norder; ++iord){
      GiveR(Kt, Rt, gt, dh); // Kx = g*g*g
      GiveR(Rt, Kt, gt, dh); // Rx = g*g*g*g -> order 2
    }
    double intgD=0;
    for (int l=0; l<Rt.extent(1); ++l) // polynomial \sum_{m=0}^n Rx(Nbin-1,m)*t^m
      intgD += Rt(Nbin-1,l);
    double nrm = pow(intgD,-1.0/Ntime);
    //cout<<"intgD["<<norder<<"]="<<intgD<<" nrm="<<nrm<<endl;
    gt *= nrm;
    
    if (SaveData){
      double alpha = pow(integral0,1./Ntime);
      //cout<<"integral0="<<integral0<<" alpha="<<alpha<<" Ntime="<<Ntime<<endl;
      //for (int it=0; it<Ntime; it++){
      std::ofstream out( (std::string("tmeasure_weight.")+to_string(Ntime/2)).c_str() );
      for (int i=0; i<Nbin; i++){
	double dt = beta/Nbin*(i+0.5);
	double f0_time = exp(-Et(0)*dt)/alpha;
	out<< dt << " "<< gt(i) << "  "<< f0_time <<endl;
      }
      out.close();
    }
  }
  void Add_to_histogram(double dk_hist, const bl::Array<double,1>& times)
  {
    //cout<<"Adding to histogram with dk_hist="<<dk_hist<<endl;
    double xbeta=times(Ntime)-times(0);
    if (fabs(xbeta-beta)>1e-5) cout<<"ERROR xbeta="<<xbeta<<" Ntime="<<Ntime<<" times="<<times<<" times("<<Ntime+1<<")="<<times(Ntime+1)<<" size="<<times.extent(0)<<endl;
    for (int it=0; it<Ntime; it++){
      double tc = times(it+1)-times(it);
      int ip = tc/beta*Nbin;
      if (ip<0 or ip>=Nbin){
	cout << "ip=" << ip << " it="<<it<<" times="<<times<<" has wrong value as t=" << tc <<" is beyond beta="<< beta << endl;
	exit(1);
      }
      t_hist(ip) += dk_hist;
      //cout<<"t_hist2("<<it<<","<<ip<<")="<<t_hist2(it,it,ip)<<endl;
    }
  }
  double Histogram_Sum()
  {
    double dsum=0;
    for (int ip=0; ip<Nbin; ++ip)
      dsum += t_hist(ip);
    return dsum/Nbin;
  }
  void WriteStatus(ofstream& out)
  {
    if (self_consistent){
      //out << "# gt["<<Nbin<<"]  t_hist["<<Nbin<<"]" << endl;
      for (int i=0; i<Nbin; i++){
	double t = beta/Nbin*(i+0.5);
	out<<t<<" "<<gt(i)<<" "<<t_hist(i)<<endl;
      }
    }else{
      out << "# t_hist["<<Nbin<<"]" << endl;
      for (int i=0; i<Nbin; ++i){
	double t = beta/Nbin*(i+0.5);
	out<<t<<" "<<t_hist(i)<<" "<<t_hist(i)<<endl;
      }
    }
  }
  bool ReadStatus(ifstream& inp){
    std::string line;
    //std::getline (inp,line); // comment gt[Ntime,Nbin]
    gt.resize(Nbin);
    for (int i=0; i<Nbin; i++){
      double x, y, z;
      inp >> x >> y >> z;
      gt(i)=y;
      t_hist(i)=z;
      if (!inp){
	cout << "WARN: reading status file and could not read gt" << endl;
	return false;
      }
      //std::getline (inp,line);
    }
    double dh=beta/Nbin;
    // normalization
    Array<double,2> Kt(Nbin,1);
    for (int i=0; i<Nbin; ++i) Kt(i,0) = gt(i);
    Array<double,2> Rt;
    GiveR(Rt, Kt, gt, dh); // Rx = g*g  -> order=1 and Ntime=2
    int norder=Ntime/2;
    for (int iord=1; iord<norder; ++iord){
      GiveR(Kt, Rt, gt, dh); // Kx = g*g*g
      GiveR(Rt, Kt, gt, dh); // Rx = g*g*g*g -> order 2
    }
    double intgD=0;
    for (int l=0; l<Rt.extent(1); ++l) // polynomial \sum_{m=0}^n Rx(Nbin-1,m)*t^m
      intgD += Rt(Nbin-1,l);
    double nrm = pow(intgD,-1.0/Ntime);
    //cout<<"histogram."<<norder<<" successfully read. Normalization constant["<<norder<<"]="<<intgD<<endl;
    gt *= nrm;
    
    //for (int i=0; i<Nbin; ++i) t_hist(i)=gt(i)*1000;
    
    self_consistent=true;
    return true;
  }
  bool ReadHistogram(int norder, double& V0norm, double& Vol_diagrams0, double& Vol_diagrams1, ostream& out)
  {
    string filename = "histogram."s+to_string(norder)+".dat";
    ifstream inp(filename);
    string comm;
    inp>>comm;
    if (inp and !inp.eof()){
      inp>>V0norm>>Vol_diagrams0>>Vol_diagrams1;
      getline (inp,comm);
      out<<"Read "<<filename<<" and set V0norm["<<norder<<"]="<<V0norm<<" Vol_dig0["<<norder<<"]="<<Vol_diagrams0<<" Vol_diag1["<<norder<<"]="<<Vol_diagrams1<<endl;
      return ReadStatus(inp);
    }else{
      out<<"Did not read "<<filename<<". Starting from scratch using V0norm["<<norder<<"]="<<V0norm<<" Vol_dig0["<<norder<<"]="<<Vol_diagrams0<<" Vol_diag1["<<norder<<"]="<<Vol_diagrams1<<endl;
      return false;
    }
  }
  bool WriteHistogram(int norder, double V0norm, double Vol_diagrams0, double Vol_diagrams1, ostream& out)
  {
    string filename = "histogram."s+to_string(norder)+".dat";
    ofstream outf(filename);
    outf.precision(16);
    outf << "# "<<V0norm<< " " << Vol_diagrams0 << " " << Vol_diagrams1 << endl;
    WriteStatus(outf);
    return true;
  }
};


class CumulativeInt1D2{
  // This class prepares cumulative integral, so that when
  // operator() is invoked, it returns  operator()(double t) -> Integrate[ h(x), {x,a,t}] within trapezoid rule
public:
  double a, b, dh, Iall;
  bl::Array<double,1> h;
  bl::Array<double,1> Ih;
  CumulativeInt1D2(const bl::Array<double,1>& _h_, double _a_, double _b_) : a(_a_), b(_b_), h(_h_.extent(0)), Ih(_h_.extent(0))
  {
    int Nbin = _h_.extent(0);
    dh = (b-a)/Nbin;
    for (int i=0; i<Nbin; ++i) h(i) = _h_(i);
    double cum = 0;
    Ih = 0;
    for (int i=0; i<Nbin; ++i){
      Ih(i) = cum;
      cum += h(i)*dh;
    }
    Iall = cum;
    //cout << "Ih(N-1)=" << Ih(Nbin-1) << endl;
  }
  double operator()(double t)
  {
    int Nbin = h.extent(0);
    int ii = static_cast<int>( (t-a)/(b-a)*Nbin );
    if (ii < 0) return 0;
    if (ii >= Nbin ) return Iall;
    return Ih(ii) + h(ii)*(t-a-ii*dh);
  }
};

bl::Array<double,1> NumericInt1DType2(const bl::Array<double,2>& gx, const bl::Array<double,3>& gx_off, double a, double b, double _two_, int i_first_time=1, int r=3)
{
  //  This routine numerically (exactly) integrates the 1D function of the form
  //   It[isk] = \int_a^b dt1 g[isk](t1)  \product_{ik\ne isk} \int_a^b dt2 g[ik](t2) g_off[isk,ik](t2-t1)
  //  in which g[isk] and g_off[isk,ik] are pice-wise step functions.
  //
  //   It[isk] = \sum_{j=0}^{Nbin} \int_{a+j*dh}^{a+(j+1)*dh} dt1 g[isk](j) \product_{ik\ne isk} ( \sum_{i1=0}^{Nbin} g[ik](i1) \int_{a+i1*dh-t1}^{a+(i1+1)*dh-t1} dx g_off[isk,ik](x) )
  //
  //  For pice-wise function g_off[isk,ik](x) we compute cumulative integral using trapezoid rule 
  //   I_ci(x) = \int_{a}^x  g_off(u) du
  //
  //  We first compute \int_{a+i1*dh-t1}^{a+(i1+1)*dh-t1} dx g_off[isk,ik](x)
  //  using this precomputed cumulative integrals.
  //  For variable t1 we use very dense mesh of Nr points, i.e., l=0...(Nr-1) and t1 = a + dh*j + dh* l/(Nr-1), and at the end use it with the Romberg method.
  //  We compute the following quantity in the inside loop:
  //                 \sum_{i1=0}^{Nbin} g[ik](i1) \int_{a+i1*dh-t1}^{a+(i1+1)*dh-t1} dx g_off[isk,ik](x)
  //  and than we evaluate the product \product_{ik\ne isk}. All these quantities are still functions of t1.
  // Finally we carry out the integral \int_{a+j*dh}^{a+(j+1)*dh} dt1  using Romberg method
  // and we need to sum over j and multiply with    g[isk](j) . 
  // 
  //int Nloops = gx.extent(0); // This is N=Nloops-1 in the above equation
  int Nbin = gx.extent(1);   // number of bins
  //int Noff = Nloops-1- i_first_time; // Number of off-diagonal functions Noff above
  int Nr = ipower(2,r)+1;                // How many points is used for romberg integration below. It needs to be exact for polynomial of correct degree.
  double dh = (b-a)/Nbin;
  bl::Array<double,1> Ps(Nr), Pt(Nr);
  bl::Array<double,1> t1(Nr);
  bl::Array<double,1> Integ(gx_off.extent(0));
  Integ=0;
  for (int isk=i_first_time; isk<gx_off.extent(0); isk++){ // variable j in the above equation
    double dsum = 0.0;
    for (int j=0; j<Nbin; ++j){                        // bins corresponding to variable j
      for (int l=0; l<Nr; ++l) t1(l) = a + dh*j + dh * l/(Nr-1.); // more precise mesh (for romberg integration) of variable r_j, i.e., the last integral in the above equation
      Pt(bl::Range::all()) = 1.0;                         // Will contain the product of functions \product_{i\ne j and i=1}^{N} .....
      // dsum += Int[ gx[isk](t1) Int[ h[isk,ik](t2-t1)*gx[ik](t2), {t2,a+i1*dh,a+(i1+1)*dh} ] , {t1,a+j*dh,a+(j+1)*dh}]
      for (int ik=i_first_time; ik<gx_off.extent(1); ++ik){     // Here we compute the product i\ne j, and we use i==ik and j==isk
	if (ik!=isk){
	  // This class calculates numerical integral of picewise function Integrate[gx_off[isk,ik](u),{u,a,x}] at an arbitrary point x. 
	  CumulativeInt1D2 I_ci(gx_off(isk,ik,bl::Range::all()), _two_*a, _two_*b); 
	  Ps(bl::Range::all()) = 0;                      // Will contain the product of functions \product_{i\ne j and i=1}^{N}
	  for (int i1=0; i1<Nbin; i1++)
	    for (int l=0; l<Nr; ++l)
	      // Int[ gx[isk](t1) Int[ h[isk,ik](t2-t1)*gx[ik](t2), {t2,a+i1*dh,a+(i1+1)*dh} ] , {t1,a+j*dh,a+(j+1)*dh}]
	      Ps(l) += ( I_ci( a + (i1+1.)*dh - t1(l) ) - I_ci( a + i1*dh - t1(l) ) ) * gx(ik,i1); 
	  Pt(bl::Range::all()) *= Ps(bl::Range::all()); // Product over all Ps functions. There are Noff products.
	  //cout << "isk=" << isk << " ik=" << ik << " Ps=" << Ps << " Pt=" << Pt << endl;
	}
      }
      double Pt_Int = romberg(Pt, dh);
      dsum += Pt_Int * gx(isk,j);
    }
    Integ(isk) = fabs(dsum);
  }
  return Integ;
}

class WmeassureWeight{
protected:
  //Meassure function used here has the form:
  //
  //   fm(x) = 1/Nxs * \sum_{i=0}^{Nxs-1} gt[i](x_i) \prod_{j\ne i} gt[j](x_j)*gt_off[i,j](x_j-x_i)
  //
  // Note that g[i](x) functions have range [a,b] with Nbin points
  // while gt_off[i,j](x) have range [two*a,two*b]
  // Both are computed self-consistently from t_hist, which is histogram of projections to variable x_i and to
  // convolution of variables x_j-x_i.
  //
  // Also notice that the functions are stored with symmetry taken into account, which allows us to use
  //     gt[i] = gt[indsym(i,i)] and
  //     gt_off[i,j]=gt[indsym(i,j)]
  // here indsym at the order 2 has the form
  // indsym=[[0,2,2]
  //         [3,1,4]
  //         [3,4,1]]
  // and at order 4 indsym is
  //indsym=[[0, 2, 2, 2, 2]
  //        [3, 1, 4, 4, 4]
  //        [3, 4, 1, 4, 4]
  //        [3, 4, 4, 1, 4]   
  //        [3, 4, 4, 4, 1]]
  //
public:
  bool self_consistent;
  double x_cutoff_a, x_cutoff_b, x_width;
  double pref_time;
  int Nbin, Nxs;
  bl::Array<double,2> gt;
  bl::Array<double,2> t_hist;
  const bool OFFDT=true;
  const bool KEEP_TAILS=false;
  const double _two_=2.0;
  double x_width2, norm_gauss;
  Array<int,2> indsym; int Nsym;
public:
  // 1 : w-acceptanece=0.01599, occurence=0.04447
  // 2 : w-acceptanece=0.01599, occurence=0.04447
  void set_indsym()
  {
    Nsym=5;
    indsym(0,0)=0;
    for (int i=1; i<Nxs; ++i){
      indsym(i,i)=1;
      indsym(0,i)=2;
      indsym(i,0)=3;
      for (int j=1; j<Nxs; ++j)
	if (i!=j) indsym(i,j)=4;
    }
  }
  WmeassureWeight(){};
  WmeassureWeight(double _x_width_, double _x_cutoff_a_, double _x_cutoff_b_, int _Nbin_, int _Nxs_):
    x_width(_x_width_), x_cutoff_a(_x_cutoff_a_), x_cutoff_b(_x_cutoff_b_), Nbin(_Nbin_), Nxs(_Nxs_), self_consistent(false)
  {
    //t_hist2.resize(Nxs,Nxs,Nbin);
    //t_hist2=0;
    //EF = kF*kF;
    pref_time = 1/(sqrt(pi)*x_width);
    
    double debug_x_width = 1.0;
    norm_gauss = 1/(sqrt(pi)*debug_x_width);
    x_width2 = debug_x_width*debug_x_width;

    
    indsym.resize(Nxs,Nxs);
    set_indsym();
    //cout<<"indsym["<<Nxs-1<<"]="<<indsym<<endl;
    t_hist.resize(Nsym,Nbin);
    t_hist=0;
  }
  void Init(double _x_width_, double _x_cutoff_a_, double _x_cutoff_b_, int _Nbin_, int _Nxs_)
  {
    x_width = _x_width_;
    x_cutoff_a = _x_cutoff_a_;
    x_cutoff_b = _x_cutoff_b_;
    Nbin = _Nbin_;
    Nxs = _Nxs_;
    self_consistent = false;
    //
    //t_hist2.resize(Nxs,Nxs,Nbin);
    //t_hist2=0;
    pref_time = 1/(sqrt(pi)*x_width);
    double debug_x_width = 1.0;
    norm_gauss = 1/(sqrt(pi)*debug_x_width);
    x_width2 = debug_x_width*debug_x_width;


    indsym.resize(Nxs,Nxs);
    set_indsym();
    //cout<<"indsym["<<Nxs-1<<"]="<<indsym<<endl;
    t_hist.resize(Nsym,Nbin);
    t_hist=0;
  }
  double trst(double x){ return fabs(x); }
  void WriteStatus(ofstream& out)
  {
    out << "# gt["<<Nxs<<","<<Nbin<<"]=" << endl;
    for (int i=0; i<Nbin; i++){
      //for (int it=0; it<Nxs; ++it)
      for (int it=0; it<Nsym; ++it) out << gt(it,i) << " ";
      out << endl;
    }
    /*
    if (OFFDT){
      out << "# gt_off["<<Nxs<<","<<Nxs<<","<<Nbin<<"]=" << endl;
      for (int i=0; i<Nbin; ++i){
	for (int its=0; its<Nxs; ++its)
	  for (int it=0; it<Nxs; ++it)
	    out << gt_off(its,it,i) << " ";
	out << endl;
      }
    }
    */
    out << "# t_hist["<<Nxs<<","<<Nbin<<"]=" << endl;
    for (int i=0; i<Nbin; ++i){
      //for (int its=0; its<Nxs; ++its)
      //for (int it=0; it<Nxs; ++it)
      //out << t_hist2(its,it,i) << " ";
      for (int it=0; it<Nsym; ++it) out << t_hist(it,i) << " ";
      out << endl;
    }
  }
  bool ReadStatus(ifstream& inp, ostream& cout){
    std::string line;
    std::getline (inp,line); // comment gt[Nxs,Nbin]
    //gt.resize(Nxs,Nbin);
    gt.resize(Nsym,Nbin);
    for (int i=0; i<Nbin; i++){
      //for (int it=0; it<Nxs; ++it){
      for (int it=0; it<Nsym; ++it){
	double x;
	inp >> x;
	gt(it,i) = x;
	if (!inp){
	  cout << "WARN: reading status file and could not read gt" << endl;
	  return false;
	}
      }
      std::getline (inp,line);
    }
    /*
    if (OFFDT){
      std::getline (inp,line); // comment gt_off[Nxs,Nbin]=
      gt_off.resize(Nxs,Nxs,Nbin);
      for (int i=0; i<Nbin; ++i){
	for (int its=0; its<Nxs; ++its){
	  for (int it=0; it<Nxs; ++it){
	    long double x;
	    inp >> x;
	    if (std::isnan(x)) x=0;
	    gt_off(its,it,i) = x;
	    if (!inp){
	      cout << "WARN: reading status file and could not read gt_off" << endl;
	      return false;
	    }
	  }
	}
	std::getline (inp,line);
      }
    }
    */
    std::getline (inp,line); // comment t_hist2[Nxs,Nbin]=
    t_hist.resize(Nsym,Nbin);
    //cout<<"comment line="<<line<<endl;
    for (int i=0; i<Nbin; ++i){
      //for (int its=0; its<Nxs; ++its){
      for (int it=0; it<Nsym; ++it){
	double x=0;
	try {// I had some issues reading numbers beyond range, like 1e-318, which would fail reading.
	  inp >> x;
	}catch (const std::ios_base::failure& e) {
	  std::cerr << "Caught a failbit exception: " << e.what() << "\n";
	  inp.clear(); // Clear the failbit
	  x=0;
	}
	if (std::isnan(x)) x=0;
	t_hist(it,i) = x;
	//cout<<"i="<<i<<"it="<<it<<" x="<<x<<endl;
	if (!inp){
	  cout << "WARN: reading status file and could not read t_hist2" << endl;
	  inp>>x;
	  return false;
	}
      }
      std::getline (inp,line);
    }
    self_consistent=true;
    cout << "mweight.status successfully read" << endl;
    return true;
  }
  double N0w(double x) const
  { // meassuring function for time distribution
    // exp(-x*x/(x_width*x_width))/(sqrt(pi)*x_width)
    return pref_time * exp(-sqr(x/x_width));
  }
  double Nsw(double x, int it) const{
    const double small=1e-16;
    if (x<x_cutoff_a or x>x_cutoff_b) return 0;
    int ip = static_cast<int>((x-x_cutoff_a)/(x_cutoff_b-x_cutoff_a)*Nbin);
    if (ip>=Nbin or ip<0) return 0;
    double res = gt(it,ip);
    return res > small ? res : small;
  }
  inline double Nw(double x, int it) const{
    return self_consistent ? Nsw(x,it) : N0w(x);
  }
  double operator()(const bl::Array<double,1>& xs) const
  { // Here we actually evaluate the function:
    // fm(x) = 1/Nxs * \sum_{i=0}^{Nxs-1} gt[i](x_i) \prod_{j\ne i} gt[j](x_j)*gt_off[i,j](x_j-x_i)
    double PQ_new = 1.0;
    if (! self_consistent){
      for (int ix=0; ix < Nxs; ix++){
	PQ_new *= N0w( xs(ix) );
	//cout<<" N0w("<<xs(ix)<<")="<<N0w(xs(ix))<<" PQ="<<PQ_new<<endl;
      }
    }else{
      for (int ix=0; ix < Nxs; ix++) PQ_new *= Nsw( xs(ix), indsym(ix,ix) );
      if (OFFDT and Nxs>0){
	double Psm = 0;
      	for (int its=0; its < Nxs; ++its){
	  double P = 1.0;
	  for (int it=0; it < Nxs; ++it){
	    if (it!=its){
	      double tc = xs(it) - xs(its);
	      double pp = 0;
	      if (tc>= x_cutoff_a*_two_ and tc< x_cutoff_b*_two_){
		int ip = static_cast<int>((tc-_two_*x_cutoff_a)/(_two_*(x_cutoff_b-x_cutoff_a))*Nbin);
		if (ip>=Nbin) ip=Nbin-1;
		//pp = gt_off(its,it,ip);
		pp = gt(indsym(its,it),ip);
	      }
	      P *= pp;
	    }
	  }
	  Psm += P;
	}
	PQ_new *= Psm/Nxs;
      }
    }
    return PQ_new;
  }
  void Recompute(bool SaveData=true){
    // Here we normalize the function
    //     fm(x) = 1/Nxs * \sum_{i=0}^{Nxs-1} gt[i](x_i) \prod_{j\ne i} gt[j](x_j)*gt_off[i,j](x_j-x_i)
    // First we take the collected histograms t_hist[i,j] and broaden it and make gt[i,j](x) \propto t_hist[i,j](x).
    // Next we numerically compute the integral and make sure that Integrate[fm(x),x]=1
    //
    if (Nxs>0){
      self_consistent=true; // from now on it is self-consistent
      // broadening
      gt.resize(Nsym,Nbin);
      gt=0;
      for (int it=0; it<Nsym; it++){ 
	gt(it,0) = trst( 0.5*(t_hist(it,0)+t_hist(it,1)) );
	for (int i=1; i<Nbin-1; i++)
	  gt(it,i) = trst( (t_hist(it,i-1)+t_hist(it,i)+t_hist(it,i+1))/3. );
	gt(it,Nbin-1) = trst( 0.5*(t_hist(it,Nbin-2)+t_hist(it,Nbin-1)) );
      }
      for (int it=0; it<2; it++){
	// First we compute rough normalization of diagonal parts by requiring that they are separately normalized, i.e.,
	//  Integrate[ gt[i](x_i), {xi,a,b}] = 1
	// Since gt[i] and picewise function, the (trapezoid) sum is exact integral.
	double intg1D = sum( gt(it,bl::Range::all()) ) * (x_cutoff_b-x_cutoff_a)/Nbin;
	gt(it,bl::Range::all()) *= 1./intg1D; 
	// If we want to keep tails so that the function does not exactly vanish in the tails.
	if (KEEP_TAILS){
	  for (int i=0; i<Nbin; i++){
	    double x = x_cutoff_a + (i+0.5)*(x_cutoff_b-x_cutoff_a)/Nbin;
	    double g0 = pref_time * exp(-sqr(x/x_width));
	    if (gt(it,i)<g0) gt(it,i)=g0;
	  }
	  double intg1D = sum( gt(it,bl::Range::all()) ) * (x_cutoff_b-x_cutoff_a)/Nbin;
	  gt(it,bl::Range::all()) *= 1./intg1D;
	}
      }
      if (OFFDT){
	// Here we normalize off-diagonal part of the function.
	// Keep in mind that the function is
	//     fm(x) = 1/Nxs * \sum_{i=0}^{Nxs-1} gt[i](x_i) \prod_{j\ne i} gt[j](x_j)*gt_off[i,j](x_j-x_i)
	// We will require
	//     1 = \int_a^b dx_i gt[i](x_i) \prod_{j\ne i} \int_a^b dx_j gt[j](x_j)*gt_off[i,j](x_j-x_i)
	//
	// We will produce gt1[i](x_i) == gt[i](x_i) and gt_off[i,j](x) from their symmetry equivalents.
	// The integral will need to compute the following
	//  nmc[iks] = \int_a^b dx_{ikx} gt1[iks](x_{iks})* \prod_{ik\ne iks} \int_a^b dx_{iks} gt_off[iks,ik](x_{ik}-x_{iks})*gt1[ik](x_{ik})
	// We here first prepare gt1[ik](t) = gt[indsym(it,it)](t) and
	//                       gt_off[iks,ik] = gt[indsym(iks,ik)](t)
	// Than we call NumericInt1DType2, which will evaluate the above integral
	Array<double,2> gt1(Nxs,Nbin);
	Array<double,3> gt_off(2,Nxs,Nbin);
	for (int it=0; it<Nxs; ++it) gt1(it,Range::all()) = gt(indsym(it,it),Range::all());
	for (int its=0; its<2; ++its){
	  for (int it=0; it<Nxs; ++it){
	    gt_off(its,it,Range::all()) = gt(indsym(its,it),Range::all());
	  }
	}
	// Evaluates numerically the integral:
	//  nmc[iks] = \int_a^b dx_{ikx} gt1[iks](x_{iks})* \prod_{ik\ne iks} \int_a^b dx_{iks} gt_off[iks,ik](x_{ik}-x_{iks})*gt1[ik](x_{ik})
	bl::Array<double,1> nmc = NumericInt1DType2(gt1, gt_off, x_cutoff_a, x_cutoff_b, _two_, 0);
	//cout<<"nmc="<<nmc<<" expecting nmc[1]=nmc[2]=...."<<endl;
	// Finally normalizing gt_off[iks,ik] by (Nxs-1)-th root, because there are (Nxs-1) gt_off functions in each product.
	double nrm=1;
	// Because of the special symmetry, we actually have only two independent normalization factors nmc[0], and nmc[1]=nmc[2]=nmc[3]...
	// hence normalization can be done over all components by normalizing
	// gt[indsym(0,1)], gt[indsym(1,0)] and gt[indsym(1,2)]. Notice that this relays on the current symmetry such that we have only
	// three independent off-diagonal functions. If more independ functions are introduced, this will need to change.
	// (We also have two diagonal functions).
	nrm = 1./std::pow(nmc(0),1./(Nxs-1.));
	gt(indsym(0,1),Range::all()) *= nrm;
	nrm = 1./std::pow(nmc(1),1./(Nxs-1.));
	gt(indsym(1,0),Range::all()) *= nrm;
	gt(indsym(1,2),Range::all()) *= nrm;
	/*
	for (int its=0; its<Nxs; its++){
	  double nrm = 1./std::pow(nmc(its),1./(Nxs-1.));
	  for (int it=0; it<Nxs; ++it)
	    gt_off(its,it,bl::Range::all()) *= nrm;
	}
	*/
	/*
	{
	  Array<double,2> gt1(Nxs,Nbin);
	  Array<double,3> gt_off(Nxs,Nxs,Nbin);
	  for (int its=0; its<Nxs; ++its){
	    gt1(its,Range::all()) = gt(indsym(its,its),Range::all());
	    for (int it=0; it<Nxs; ++it){
	      gt_off(its,it,Range::all()) = gt(indsym(its,it),Range::all());
	    }
	  }
	  bl::Array<double,1> nmc = NumericInt1DType2(gt1, gt_off, x_cutoff_a, x_cutoff_b, _two_, 0);
	  cout<<"checking nmc="<<nmc<<endl;
	}
	*/
      }
      if (SaveData){
	for (int it=0; it<Nsym; it++){
	  std::ofstream out( "xmeasure_weight."s+to_string(Nxs-1)+"."s+to_string(it) );
	  for (int i=0; i<Nbin; i++){
	    double x = x_cutoff_a + (i+0.5)*(x_cutoff_b-x_cutoff_a)/Nbin;
	    if (it>=2) x *= _two_;
	    out<< x << " "<< gt(it,i) << "  "<< N0w(x) <<endl;
	  }
	  out.close();
	}
      }
    }
  }
  void Add_to_histogram(double dk_hist, const bl::Array<double,1>& xs)
  {
    if (OFFDT){
      for (int its=0; its<Nxs; its++){
	for (int it=0; it<Nxs; it++){
	  double tc;
	  if (it == its){ // diagonal
	    tc = xs(it);
	    int ip = (tc - x_cutoff_a)/(x_cutoff_b-x_cutoff_a)*Nbin;
	    if (ip>=Nbin) ip=Nbin-1;
	    if (ip<0) ip=0;
	    t_hist(indsym(its,it),ip) += dk_hist;
	  }else{
	    tc = xs(it)-xs(its);
	    int ip = (tc - _two_*x_cutoff_a)/(_two_*(x_cutoff_b-x_cutoff_a))*Nbin;
	    if (ip>=Nbin) ip=Nbin-1;
	    if (ip<0) ip=0;
	    t_hist(indsym(its,it),ip) += dk_hist;
	  }
	}
      }
    }else{
      for (int it=0; it<Nxs; it++){
	double tc = xs(it);
	int ip = (tc - x_cutoff_a)/(x_cutoff_b-x_cutoff_a)*Nbin;
	if (ip>=Nbin) ip=Nbin-1;
	if (ip<0) ip=0;
	t_hist(indsym(it,it),ip) += dk_hist;
      }
    }
  }
  double Histogram_Sum()
  {
    return sum(t_hist);
  }
  double Normalize_K_histogram(){
    double dnrm = 2./(sum(t_hist(0,Range::all())) + sum(t_hist(1,Range::all())));
    t_hist *= dnrm;
    return dnrm;
  }
  bool ReadHistogram(int norder, double& V0norm, double& Vol_diagrams0, double& Vol_diagrams1, ostream& out)
  {
    string filename = "histogram."s+to_string(norder)+".dat";
    ifstream inp(filename);
    string comm;
    inp>>comm;
    if (inp and !inp.eof()){
      inp>>V0norm>>Vol_diagrams0>>Vol_diagrams1;
      getline (inp,comm);
      out<<"Read "<<filename<<" and set V0norm["<<norder<<"]="<<V0norm<<" Vol_dig0["<<norder<<"]="<<Vol_diagrams0<<" Vol_diag1["<<norder<<"]="<<Vol_diagrams1<<endl;
      return ReadStatus(inp,out);
    }else{
      out<<"Did not read "<<filename<<". Starting from scratch using V0norm["<<norder<<"]="<<V0norm<<" Vol_dig0["<<norder<<"]="<<Vol_diagrams0<<" Vol_diag1["<<norder<<"]="<<Vol_diagrams1<<endl;
      return false;
    }
  }
  bool WriteHistogram(int norder, double V0norm, double Vol_diagrams0, double Vol_diagrams1, ostream& out)
  {
    string filename = "histogram."s+to_string(norder)+".dat";
    ofstream outf(filename);
    outf.precision(16);
    outf << "# "<<V0norm<< " " << Vol_diagrams0 << " " << Vol_diagrams1 << endl;
    WriteStatus(outf);
    return true;
  }
};
