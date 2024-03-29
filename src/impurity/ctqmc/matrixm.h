// @Copyright 2007 Kristjan Haule
// 
#include "timer.h"
class MatrixM{
  function1D<dcomplex> sumt;
  function2D<dcomplex> sum1, sum2;
  double beta;
  const mesh1D* ptau;
  function2D<spline1D<double>* > vfun;
  const mesh1D* piom;
  function1D<double> d_ln, d_nl, d_lm, d_ml;
  double d_nn, ratio;
  function1D<double> Lt, Rt, Lt2, Rt2;
  function1D<double> temp;
  function1D<int> ipiv;
public:
  function2D<dcomplex> Gf;
  function2D<dcomplex> Ff;
private:
  function2D<dcomplex> sum2c;
  function2D<dcomplex> sum3;
  function1D<dcomplex> dexp0;
  function2D<dcomplex> dG;
  vector<function2D<dcomplex> > dMv;
  const function2D<int>* ptfl_index;
  deque<int> v2fl_index;
  int nomv; // maximum fermionic frequency for vertex computation
  function1D<double> M_ul;
public:
  //double tm1, tm2, tm3, tm4, tm5, tm6, tm7, tm8, tm9, tm10, tm11, tm12, tm13;
public:
  void copy_data(MatrixM& source);
  // small inlined functions to access members
  int fullsize(){return Lt.fullsize();}
  const function2D<dcomplex>& gf() const { return Gf;}
  double iom_size(){return piom->size();}
  
  void SetUp(int Nmax, const mesh1D& tau_, const function2D<spline1D<double>* >& fun_, const function2D<int>& tfl_index_, const mesh1D& iom, double beta_, deque<int>& v2fl_index_, int nomv_, int svd_lmax);
  
  double AddDetRatio(function2D<double>& MD, int Nk, double t_start, int btype_s, double t_end, int btype_e, const nIntervals& interval);
  double AddUpdateMatrix(function2D<double>& MD, int Nk, double t_start, int is, int btype_s, double t_end, int ie, int btype_e, const nIntervals& interval, int istep);
  double AddUpdateMatrixSlow(function2D<double>& MD, int Nk, double t_start, int is, int btype_s, double t_end, int ie, int btype_e, const nIntervals& interval, int istep);
  void AddUpdate_Gf(const nIntervals& interval, vector<function2D<dcomplex> >& Mv, int istep);
  
  double RemoveDetRatio(const function2D<double>& MD, int Nk, int is, int ie);
  double RemoveDetRatio(const function2D<double>& MD, int Nk, int is1, int ie1, int is2, int ie2);
  double RemoveUpdateMatrix(function2D<double>& MD, int Nk, int is, int ie, const nIntervals& interval);
  double RemoveUpdateMatrixSlow(function2D<double>& MD, int Nk, int is, int ie, const nIntervals& interval);
  
  void RemoveUpdate_Gf(const function2D<double>& MD, int Nk, int is, int ie, const nIntervals& interval, vector<function2D<dcomplex> >& Mv);
  
  double Move_start_DetRatio(const function2D<double>& MD, int Nk, double ts_old, double ts_new, int btype_s, int is_old, const nIntervals& interval);
  double Move_end_DetRatio(const function2D<double>& MD, int Nk, double te_old, double te_new, int btype_e, int ie_old, const nIntervals& interval);
  void Move_start_UpdateMatrix(function2D<double>& MD, int Nk, double ts_old, double ts_new, int btype_s, int is_old, int is_new, const nIntervals& interval);
  void Move_end_UpdateMatrix(function2D<double>& MD, int Nk, double te_old, double te_new, int btype_e, int ie_old, int ie_new, const nIntervals& interval);
  void MoveUpdate_Gf(const nIntervals& interval, vector<function2D<dcomplex> >& Mv, int istep);
  
  
  void CleanUpdateMatrix(function2D<double>& MD, int Nk, const nIntervals& interval, function2D<double>& mD);
  void CleanUpdateGf(function2D<double>& MD, int Nk, const nIntervals& interval, vector<function2D<dcomplex> >& Mv, int istep);
  
  double AddDetRatio(function2D<double>& MD, int Nk, int btype_s1, double t_start1, int btype_e1, double t_end1, 
		     int btype_s2, double t_start2, int btype_e2, double t_end2, const nIntervals& interval);

  
  double GlobalFlipDetRatio(const nIntervals& interval, const function2D<double>& MD, const function2D<spline1D<double>* >& Delta_flipped, /*function2D<double>& C,*/ function2D<double>& mD);
  void ComputeGtau(const function2D<double>& MD, const nIntervals& interval, function2D<double>& Gtau);
  void ComputeG_svd(function2D<double>& gl, double tsign, const SVDFunc& svd, const function2D<double>& MD, const nIntervals& interval);
  void ComputeFG_svd(vector<function2D<double> >& gl, double tsign, const SVDFunc& svd, const function2D<double>& MD, const nIntervals& interval);
  void ComputeFGF_svd(vector<function2D<double> >& gs, vector<function2D<double> >& ge, vector<function3D<double> >& ges, bool Vertex_subtract, double tsign, const SVDFunc& svd, const function2D<double>& MD, const nIntervals& interval);
  void Resize(int new_Nmax);
  void ReleaseSomeTempMemory();
  
  void ComputeG_svd_debug(function2D<double>& gl, function2D<double>& Gtau, double tsign, const SVDFunc& svd, const function2D<double>& MD, const nIntervals& interval);
  
private:
  void AddComputeGs(int Nk, double t_start, int btype_s, double t_end, int btype_e, const nIntervals& interval);
  void AddComputeGs(int Nk, double t_start, int btype_s, double t_end, int btype_e, const nIntervals& interval, function1D<double>& d_ln, function1D<double>& d_nl, double& dnn);
  friend void swapGf(MatrixM& m1, MatrixM& m2, vector<function2D<dcomplex> >& Mv);
  double Interp_Delta(int bs, int be, const intpar& p);
  double Interp_Antiperiodic_Delta(int bs, int be, double dt);
  void SetMatrix_(int Nk, const nIntervals& interval, function2D<double>& mD, const function2D<spline1D<double>* >& Delta_any);
};

inline void MatrixM::SetUp(int Nmax, const mesh1D& tau_, const function2D<spline1D<double>* >& fun_, const function2D<int>& tfl_index_, const mesh1D& iom, double beta_, deque<int>& v2fl_index_, int nomv_, int svd_lmax)
{
  beta = beta_;
  ptau = &tau_;
  ptfl_index = &tfl_index_;
  vfun = fun_;
  d_ln.resize(Nmax);
  d_nl.resize(Nmax);
  d_lm.resize(Nmax);
  d_ml.resize(Nmax);
  Lt.resize(Nmax);
  Rt.resize(Nmax);
  Lt2.resize(Nmax);
  Rt2.resize(Nmax);
  ipiv.resize(Nmax);
  int dim = ptfl_index->size_N();
  nomv = nomv_;  // maximum fermionic frequency for vertex computation
  if (! common::Qsvd){
    piom = &iom;
    Gf.resize(dim*dim,iom.size());
    if (common::QHB2) Ff.resize(Nmax,iom.size());
    sum1.resize(dim,iom.size());
    sum2.resize(dim,iom.size());
    if (common::QHB2) sum2c.resize(Nmax/2+1,iom.size());
    dG.resize(dim*dim,iom.size());
    dexp0.resize(iom.size());
    sum3.resize(dim*dim,iom.size());
    sumt.resize(iom.size());
    if (common::cmp_vertex){
      dMv.resize(dim*dim);
      for (unsigned i=0; i<dMv.size(); i++) dMv[i].resize(2*nomv,2*nomv);//(-nomv,nomv,-nomv,nomv);
    }
  }else{
    //if (common::cmp_vertex)
      M_ul.resize(svd_lmax);
  }
  temp.resize(Nmax);
  v2fl_index = v2fl_index_;
  //tm1=tm2=tm3=tm4=tm5=tm6=tm7=tm8=tm9=tm10=tm11=tm12=tm13=0;
}

inline void MatrixM::Resize(int new_Nmax)
{
  d_ln.resize(new_Nmax);
  d_nl.resize(new_Nmax);
  d_lm.resize(new_Nmax);
  d_ml.resize(new_Nmax);
  Lt.resize(new_Nmax);
  Rt.resize(new_Nmax);
  Lt2.resize(new_Nmax);
  Rt2.resize(new_Nmax);
  ipiv.resize(new_Nmax);
  temp.resize(new_Nmax);
  if (!common::Qsvd && common::QHB2){
    sum2c.resize(new_Nmax/2+1,piom->size());
    //int _N_ = Ff.size_N();
    //int _Nd_ = Ff.size_Nd();
    Ff.resize_preserve(new_Nmax,piom->size());
    //Ff.resize(_N_,_Nd_); // virtual resizing, does not allocate/deallocate memory
  }
}

void MatrixM::ReleaseSomeTempMemory()
{
  d_ln.ReleaseMemory();//~function1D();
  d_nl.ReleaseMemory();//~function1D();
  d_lm.ReleaseMemory();//~function1D();
  d_ml.ReleaseMemory();//~function1D();
  Lt.ReleaseMemory();//~function1D();
  Rt.ReleaseMemory();//~function1D();
  Lt2.ReleaseMemory();//~function1D();
  Rt2.ReleaseMemory();//~function1D();
  ipiv.ReleaseMemory();//~function1D();
  
  if (! common::Qsvd){
    Gf.ReleaseMemory();//~function2D();
    if (common::QHB2) Ff.ReleaseMemory();//~function2D();
    sum1.ReleaseMemory();//~function2D();
    sum2.ReleaseMemory();//~function2D();
    if (common::QHB2) sum2c.ReleaseMemory();//~function2D();
    dG.ReleaseMemory();//~function2D();
    dexp0.ReleaseMemory();//~function1D();
    sum3.ReleaseMemory();//~function2D();
    sumt.ReleaseMemory();//~function1D();
    if (common::cmp_vertex)
      for (unsigned i=0; i<dMv.size(); i++) dMv[i].ReleaseMemory();//~function2D();
  }else{
    //if (common::cmp_vertex)
    M_ul.ReleaseMemory();//~function1D();
  }
  temp.ReleaseMemory();//~function1D();
}


inline double Interp_Any_Delta(const spline1D<double>& Delta_any, const intpar& p, bool bs_equal_be)
{
  // Interpolates Delta(tau) and makes sure that Delta is causal
  double Delta_tau = Delta_any(p); // interpolation of Delta(tau) in point p.x==tau
  if (bs_equal_be && Delta_tau>-common::minDeltat)
    Delta_tau = -common::minDeltat; // Fix any causality problem

  return Delta_tau;
}

inline double MatrixM::Interp_Delta(int bs, int be, const intpar& p)
{ // Interpolates Delta stored in MatrixM class. This Interpolation is for speed only.
  // It could be enough to have only Interp_Antiperiodic_Delta. In the latter case, 
  // program would not be optimized for correlated lookups.
  return Interp_Any_Delta(*vfun[bs][be], p, bs==be);
}

inline double MatrixM::Interp_Antiperiodic_Delta(int bs, int be, double dt)
{  // Interpolates Delta antiperiodically (because Delta is fermionic)
  if (dt>0)
    return Interp_Delta(bs, be, ptau->Interp(dt));
  else
    return -Interp_Delta(bs, be, ptau->Interp(dt+beta));
}

inline void MatrixM::AddComputeGs(int Nk, double t_start, int btype_s, double t_end, int btype_e, const nIntervals& interval)
{
  tint pos1 = ptau->InitInterpLeft();
  tint pos2 = ptau->InitInterpLeft();
  for (int i=0; i<Nk; i++){
    int bs = interval.btype_s(i);
    int be = btype_e;
    if (interval.time_s(i)<t_end)
      d_ln[i] = -Interp_Delta(bs, be, ptau->InterpLeft(interval.time_s(i)-t_end+beta,pos1));
    else
      d_ln[i] = Interp_Delta(bs, be, ptau->InterpLeft(interval.time_s(i)-t_end,pos2));
  }
  
  pos1 = ptau->InitInterpRight();
  pos2 = ptau->InitInterpRight();
  for (int i=0; i<Nk; i++){
    int bs = btype_s; 
    int be = interval.btype_e(i);
    if (t_start<interval.time_e(i))
      d_nl[i] = -Interp_Delta(bs, be, ptau->InterpRight(t_start-interval.time_e(i)+beta,pos1));
    else
      d_nl[i] = Interp_Delta(bs, be, ptau->InterpRight(t_start-interval.time_e(i),pos2));
  }
  d_nn = Interp_Antiperiodic_Delta(btype_s, btype_e, t_start-t_end);
}

inline double MatrixM::AddDetRatio(function2D<double>& MD, int Nk, double t_start, int btype_s, double t_end, int btype_e, const nIntervals& interval)
{
  AddComputeGs(Nk, t_start, btype_s, t_end, btype_e, interval);
  double* __restrict__ _Lt = Lt.MemPt();
  double* __restrict__ _Rt = Rt.MemPt();
  const double* __restrict__ _d_ln = d_ln.MemPt();
  const double* __restrict__ _d_nl = d_nl.MemPt();
  int Nd0 = MD.fullsize_Nd();
  const double* __restrict__ _MD = MD.MemPt();
  for (int i=0; i<Nk; i++){
    double sumL=0;
    double sumR=0;
    for (int l=0; l<Nk; l++){// should be optimised by dgemv
      //sumL += MD(i,l)*_d_ln[l];
      //sumR += _d_nl[l]*MD(l,i);
      sumL += _MD[i*Nd0+l]*_d_ln[l];
      sumR += _d_nl[l]*_MD[l*Nd0+i];
    }
    _Lt[i] = sumL;
    _Rt[i] = sumR;
  }
  double qsum=0;
  for (int l=0; l<Nk; l++) qsum += _d_nl[l]*_Lt[l]; // should be optimized by ddot
  ratio = d_nn - qsum;
  return ratio;
}

inline double MatrixM::Move_start_DetRatio(const function2D<double>& MD, int Nk, double ts_old, double ts_new, int btype_s, int is_old, const nIntervals& interval)
{
  tint pos1 = ptau->InitInterpRight();
  tint pos2 = ptau->InitInterpRight();
  for (int i=0; i<Nk; i++){
    int bs = btype_s;
    int be = interval.btype_e(i);
    if (ts_old<interval.time_e(i))
      d_nl[i] = -Interp_Delta(bs, be, ptau->InterpRight(ts_old-interval.time_e(i)+beta,pos1));
    else
      d_nl[i] = Interp_Delta(bs, be, ptau->InterpRight(ts_old-interval.time_e(i),pos2));
  }
  pos1 = ptau->InitInterpRight();
  pos2 = ptau->InitInterpRight();
  for (int i=0; i<Nk; i++){
    int bs = btype_s;
    int be = interval.btype_e(i);
    if (ts_new<interval.time_e(i))
      d_ln[i] = -Interp_Delta(bs, be, ptau->InterpRight(ts_new-interval.time_e(i)+beta,pos1));
    else
      d_ln[i] = Interp_Delta(bs, be, ptau->InterpRight(ts_new-interval.time_e(i),pos2));
  }
  for (int i=0; i<Nk; i++) d_nl[i] = d_ln[i]-d_nl[i]; // new-old
  
  double sumR=0;
  for (int l=0; l<Nk; l++) sumR += d_nl[l]*MD(l,is_old);
  return 1. + sumR;
}

inline double MatrixM::Move_end_DetRatio(const function2D<double>& MD, int Nk, double te_old, double te_new, int btype_e, int ie_old, const nIntervals& interval)
{
  tint pos1 = ptau->InitInterpLeft();
  tint pos2 = ptau->InitInterpLeft();
  for (int i=0; i<Nk; i++){
    int bs = interval.btype_s(i);
    int be = btype_e;
    if (interval.time_s(i)<te_old)
      d_ln[i] = -Interp_Delta(bs, be, ptau->InterpLeft(interval.time_s(i)-te_old+beta,pos1));
    else
      d_ln[i] = Interp_Delta(bs, be, ptau->InterpLeft(interval.time_s(i)-te_old,pos2));
  }
  pos1 = ptau->InitInterpLeft();
  pos2 = ptau->InitInterpLeft();
  for (int i=0; i<Nk; i++){
    int bs = interval.btype_s(i);
    int be = btype_e;
    if (interval.time_s(i)<te_new)
      d_nl[i] = -Interp_Delta(bs, be, ptau->InterpLeft(interval.time_s(i)-te_new+beta,pos1));
    else
      d_nl[i] = Interp_Delta(bs, be, ptau->InterpLeft(interval.time_s(i)-te_new,pos2));
  }
  for (int i=0; i<Nk; i++) d_ln[i] = d_nl[i]-d_ln[i]; // new-old
  
  double sumL=0;
  for (int l=0; l<Nk; l++) sumL += MD(ie_old,l)*d_ln[l];
  return 1. + sumL;
}

inline void MatrixM::Move_start_UpdateMatrix(function2D<double>& MD, int Nk, double ts_old, double ts_new, int btype_s, int is_old, int is_new, const nIntervals& interval)
{
  double q = Move_start_DetRatio(MD, Nk, ts_old, ts_new, btype_s, is_old, interval);

  int bs = btype_s;
  if (!common::Qsvd){
    // Correction to the Green's function will be needed
    // This is because M^{new} = M^{old} + L_i R_j
    // and although M^{old} is unchanged, here and all corrections are in L and R,
    // it needs to be multiplied by the new time_s rather than with the old time_s
    for (int im=0; im<piom->size(); im++){
      double xs = ts_new*(*piom)[im]; 
      dexp0[im] = -(dcomplex(cos(xs),-sin(xs))-interval.exp_s(is_old)[im].conj())/beta;// dexp0[iom] = -(e^{-iom*ts_new}-e^{-iom*ts_old})/beta
    }
    sum3=0;
    for (int i=0; i<Nk; i++){
      int be = interval.btype_e(i);
      int ind = (*ptfl_index)[be][bs];
      double md = MD(i,is_old);
      for (int im=0; im<piom->size(); im++) sum3[ind][im] += interval.exp_e(i)[im]*md; // sum3[(be,bs)][iom] += e^{iom*te}M[te,ts_old]
    }
    for (int ind=0; ind<sum3.size_N(); ind++)
      for (int im=0; im<piom->size(); im++)
	dG[ind][im] = sum3[ind][im]*dexp0[im];  // dG[(be,bs)][iom] = -e^{iom*te}M[te,ts_old](e^{-iom*ts_new}-e^{-iom*ts_old})/beta

    if (common::QHB2){
      dcomplex* __restrict__ _Ff = Ff[is_old].MemPt();
      dcomplex* __restrict__ _dG = dG[0].MemPt();
      for (int im=0; im<piom->size(); im++) _Ff[im] += _dG[im]; // Ff[is_old][iom] += dG[iom]
    }
  
    if (common::cmp_vertex){
      for (int ind=0; ind<sum3.size_N();ind++){
	function2D<dcomplex>& _dMv = dMv[ind];
	for (int im1=0; im1<nomv; im1++){
	  dcomplex sm3 = sum3[ind][im1];
	  funProxy<dcomplex>& __dMv = _dMv[im1+nomv];
	  for (int im2=0; im2<nomv; im2++)  __dMv[im2+nomv] = sm3*dexp0[im2];  // _dMv(im1+nomv,im2+nomv) = sm3*dexp0[im2];
	  for (int im2=-nomv; im2<0; im2++) __dMv[im2+nomv] = sm3*dexp0[-im2-1].conj();//_dMv(im1+nomv,im2+nomv) = sm3*dexp0[-im2-1].conj();
	}
      }
    }
  }

  // Prepares Lt and Rt to update M and G
  Lt.resize(Nk); Rt.resize(Nk);
  for (int i=0; i<Nk; i++){
    double sumR=0;
    for (int l=0; l<Nk; l++){// should be optimised by dgemv
      sumR += d_nl[l]*MD(l,i);
    }
    Rt[i] = sumR/q;
  }
  for (int i=0; i<Nk; i++) Lt[i] = -MD(i,is_old);
  for (int i=0; i<Nk; i++)
    for (int j=0; j<Nk; j++)
      MD(i,j) += Lt[i]*Rt[j];


  // Shufle rows if time order changed because of move
  if (is_new!=is_old){
    double rt=0;
    for (int i=0; i<Nk; i++){ // remember the old row
      temp[i] = MD(i,is_old);
      rt = Rt[is_old];
    }
    if (is_new<is_old){// shuffle rows
      for (int j=is_old; j>is_new; j--){
	for (int i=0; i<Nk; i++) MD(i,j) = MD(i,j-1);
	Rt[j] = Rt[j-1];
      }
    }else{
      for (int j=is_old; j<is_new; j++){
	for (int i=0; i<Nk; i++) MD(i,j) = MD(i,j+1);
	Rt[j] = Rt[j+1];
      }
    }
    for (int i=0; i<Nk; i++) MD(i,is_new) = temp[i];
    Rt[is_new] = rt;
    
    if (!common::Qsvd && common::QHB2){
      for (int im=0; im<piom->size(); im++) sum1(0,im) = Ff(is_old,im);
      if (is_new<is_old){
	for (int j=is_old; j>is_new; j--)
	  for (int im=0; im<piom->size(); im++) Ff(j,im) = Ff(j-1,im);
      }else{
	for (int j=is_old; j<is_new; j++)
	  for (int im=0; im<piom->size(); im++) Ff(j,im) = Ff(j+1,im);
      }
      for (int im=0; im<piom->size(); im++) Ff(is_new,im) = sum1(0,im);
    }
  }
}

inline void MatrixM::Move_end_UpdateMatrix(function2D<double>& MD, int Nk, double te_old, double te_new, int btype_e, int ie_old, int ie_new, const nIntervals& interval)
{
  double q = Move_end_DetRatio(MD, Nk, te_old, te_new, btype_e, ie_old, interval);
  int be = btype_e;

  if (!common::Qsvd){
    // Correction to the Green's function will be needed
    // This is because M^{new} = M^{old} + L_i R_j
    // and although M^{old} is unchanged, here and all corrections are in L and R,
    // it needs to be multiplied by the new time_e rather than with the old time_e
    for (int im=0; im<piom->size(); im++){
      double xe = te_new*(*piom)[im];
      dexp0[im] = -(dcomplex(cos(xe),sin(xe))-interval.exp_e(ie_old)[im])/beta;
    }
    sum3=0;
    for (int i=0; i<Nk; i++){
      int bs = interval.btype_s(i);
      int ind = (*ptfl_index)[be][bs];
      double md = MD(ie_old,i);
      const funProxy<dcomplex>& exp_s = interval.exp_s(i);
      funProxy<dcomplex>& _sum3 = sum3[ind];
      for (int im=0; im<piom->size(); im++) _sum3[im] += md*exp_s[im].conj();
    }
  
    for (int ind=0; ind<sum3.size_N(); ind++){
      funProxy<dcomplex>& _dG = dG[ind];
      funProxy<dcomplex>& _sum3 = sum3[ind];
      for (int im=0; im<piom->size(); im++) _dG[im] = dexp0[im]*_sum3[im];
    }

    if (common::QHB2){
      for (int i=0; i<Nk; i++){
	double md = MD(ie_old,i);
	const funProxy<dcomplex>& exp_s = interval.exp_s(i);
	for (int im=0; im<piom->size(); im++)
	  Ff(i,im) += dexp0[im]*md*exp_s[im].conj();
      }
    }
  
    if (common::cmp_vertex){
      for (int ind=0; ind<sum3.size_N(); ind++){
	funProxy<dcomplex>& _sum3 = sum3[ind];
	function2D<dcomplex>& _dMv = dMv[ind];
	for (int im1=0; im1<nomv; im1++){
	  dcomplex dexp01 = dexp0[im1];
	  funProxy<dcomplex>& __dMv = _dMv[im1+nomv];
	  for (int im2=0; im2<nomv; im2++)  __dMv[im2+nomv] = dexp01*_sum3[im2]; //_dMv(im1+nomv,im2+nomv) = dexp01*_sum3[im2];
	  for (int im2=-nomv; im2<0; im2++) __dMv[im2+nomv] = dexp01*_sum3[-im2-1].conj(); //_dMv(im1+nomv,im2+nomv) = dexp01*_sum3[-im2-1].conj();
	}
      }
    }
  }
  
  // Prepares Lt and Rt to update M and G
  Lt.resize(Nk); Rt.resize(Nk);
  for (int i=0; i<Nk; i++){
    double sumL=0;
    for (int l=0; l<Nk; l++){// should be optimised by dgemv
      sumL += MD(i,l)*d_ln[l];
    }
    Lt[i] = sumL/q;
  }
  for (int i=0; i<Nk; i++) Rt[i] = -MD(ie_old,i);
  for (int i=0; i<Nk; i++)
    for (int j=0; j<Nk; j++)
      MD(i,j) += Lt[i]*Rt[j];
  
  // Shufle columns if time order changed because of move
  if (ie_new!=ie_old){
    double lt = Lt[ie_old];
    for (int i=0; i<Nk; i++) temp[i] = MD(ie_old,i); // remember the old row
    if (ie_new<ie_old){// shuffle columns
      for (int j=ie_old; j>ie_new; j--){
	for (int i=0; i<Nk; i++) MD(j,i) = MD(j-1,i);
	Lt[j] = Lt[j-1];
      }
      for (int i=0; i<Nk; i++) MD(ie_new,i) = temp[i];
      Lt[ie_new] = lt;
    }else{
      for (int j=ie_old; j<ie_new; j++){
	for (int i=0; i<Nk; i++) MD(j,i) = MD(j+1,i);
	Lt[j] = Lt[j+1];
      }
      for (int i=0; i<Nk; i++) MD(ie_new,i) = temp[i];
      Lt[ie_new] = lt;
    }
  }
}

inline double MatrixM::AddUpdateMatrixSlow(function2D<double>& MD, int Nk, double t_start, int is, int btype_s, double t_end, int ie, int btype_e, const nIntervals& interval, int istep)
{
  AddDetRatio(MD, Nk, t_start, btype_s, t_end, btype_e, interval);
  Lt.resize(Nk+1); Rt.resize(Nk+1);
  
  for (int i=Nk; i>ie; i--) Lt[i] = Lt[i-1];
  Lt[ie] = -1;
  for (int i=Nk; i>is; i--) Rt[i] = Rt[i-1];
  Rt[is] = -1;
  double v = 1/ratio;
  for (int i=0; i<Nk+1; i++) Lt[i] *= v;
  
  for (int i=Nk; i>ie; i--)
    for (int j=Nk; j>is; j--)
      MD(i,j) = MD(i-1,j-1);
  
  for (int i=0; i<ie; i++)
    for (int j=Nk; j>is; j--)
      MD(i,j) = MD(i,j-1);
  
  for (int i=Nk; i>ie; i--)
    for (int j=0; j<is; j++)
      MD(i,j) = MD(i-1,j);
  
  for (int i=0; i<Nk+1; i++) MD(i,is)=0;
  for (int j=0; j<Nk+1; j++) MD(ie,j)=0;
  
  for (int i=0; i<Nk+1; i++)// should be optimized by dger
    for (int j=0; j<Nk+1; j++)
      MD(i,j) += Lt[i]*Rt[j];
  
  MD.resize(Nk+1,Nk+1);

  if (!common::Qsvd && common::QHB2){
    for (int j=Nk; j>is; j--)
      for (int im=0; im<piom->size(); im++)
	Ff(j,im) = Ff(j-1,im);
    
    for (int im=0; im<piom->size(); im++) Ff(is,im)=0.0;
    Ff.resize(Nk+1,Ff.size_Nd());
  }
  
  return ratio*(1-2*((is+ie)%2));
}
inline double MatrixM::AddUpdateMatrix(function2D<double>& MD, int Nk, double t_start, int is, int btype_s, double t_end, int ie, int btype_e, const nIntervals& interval, int istep)
{
  AddDetRatio(MD, Nk, t_start, btype_s, t_end, btype_e, interval);
  Lt.resize(Nk+1); Rt.resize(Nk+1);

  double* __restrict__ _Lt = Lt.MemPt();
  double* __restrict__ _Rt = Rt.MemPt();
  for (int i=Nk; i>ie; i--) _Lt[i] = _Lt[i-1];
  _Lt[ie] = -1;
  for (int i=Nk; i>is; i--) _Rt[i] = _Rt[i-1];
  _Rt[is] = -1;
  double v = 1/ratio;
  for (int i=0; i<Nk+1; i++) _Lt[i] *= v;
  
  // ** slow equivalent **
  // for (int i=Nk; i>ie; i--)
  //   for (int j=Nk; j>is; j--)
  //     MD(i,j) = MD(i-1,j-1);
  for (int i=Nk; i>ie; i--)
    if (is<Nk) memcpy(MD[i].MemPt()+(is+1), MD[i-1].MemPt()+is, (Nk-is)*sizeof(double));
  
  // ** slow equivalent **
  // for (int i=0; i<ie; i++)
  //   for (int j=Nk; j>is; j--)
  //     MD(i,j) = MD(i,j-1);
  for (int i=0; i<ie; i++){
    double* __restrict__ _MD_i=MD[i].MemPt();
    for (int j=Nk; j>is; j--) _MD_i[j] = _MD_i[j-1];
  }

  // ** slow equivalent **
  // for (int i=Nk; i>ie; i--)
  //   for (int j=0; j<is; j++)
  //     MD(i,j) = MD(i-1,j);
  for (int i=Nk; i>ie; i--) memcpy(MD[i].MemPt(), MD[i-1].MemPt(), sizeof(double)*is);

  for (int i=0; i<Nk+1; i++) MD(i,is)=0;
  for (int j=0; j<Nk+1; j++) MD(ie,j)=0;

  // ** slow equivalent **
  // for (int i=0; i<Nk+1; i++)
  //   for (int j=0; j<Nk+1; j++)
  //     MD(i,j) += Lt[i]*Rt[j];
  for (int i=0; i<Nk+1; i++){
    double _lt = Lt[i];
    double* __restrict__ _MD_i = MD[i].MemPt();
    const double* __restrict__ _Rt = Rt.MemPt();
    for (int j=0; j<Nk+1; j++) _MD_i[j] += _Rt[j]*_lt;
  }

  MD.resize(Nk+1,Nk+1);

  if (!common::Qsvd && common::QHB2){
    /// Qsvd !!! Qsvd
    for (int j=Nk; j>is; j--)
      for (int im=0; im<piom->size(); im++)
	Ff(j,im) = Ff(j-1,im);
    
    for (int im=0; im<piom->size(); im++) Ff(is,im)=0.0;
    Ff.resize(Nk+1,Ff.size_Nd());
  }

  return ratio*(1-2*((is+ie)%2));
}

inline double MatrixM::RemoveDetRatio(const function2D<double>& MD, int Nk, int is, int ie)
{
  ratio = MD(ie,is);
  return ratio*(1-2*((is+ie)%2));
}
inline double MatrixM::RemoveDetRatio(const function2D<double>& MD, int Nk, int is1, int ie1, int is2, int ie2)
{
  ratio = MD(ie1,is1)*MD(ie2,is2)-MD(ie2,is1)*MD(ie1,is2);
  return ratio*(1-2*((is1+ie1+is2+ie2)%2));
}
inline double MatrixM::RemoveUpdateMatrixSlow(function2D<double>& MD, int Nk, int is, int ie, const nIntervals& interval)
{
  RemoveDetRatio(MD,Nk,is,ie);
  
  for (int i=0; i<Nk; i++)
    for (int j=0; j<Nk; j++){
      if (i!=ie && j!=is)
	MD(i,j) -= MD(i,is)*MD(ie,j)/MD(ie,is);
    }

  for (int i=ie; i<Nk-1; i++)
    for (int j=is; j<Nk-1; j++)
      MD(i,j) = MD(i+1,j+1);
  
  for (int i=0; i<ie; i++)
    for (int j=is; j<Nk-1; j++)
      MD(i,j) = MD(i,j+1);
  
  for (int i=ie; i<Nk-1; i++)
    for (int j=0; j<is; j++)
      MD(i,j) = MD(i+1,j);

  MD.resize(Nk-1,Nk-1);

  if (!common::Qsvd && common::QHB2){
    for (int j=is; j<Nk-1; j++)
      for (int im=0; im<piom->size(); im++)
	Ff(j,im) = Ff(j+1,im);
    Ff.resize(Nk-1,Ff.size_Nd());
  }
  
  return ratio*(1-2*((is+ie)%2));
}


inline double MatrixM::RemoveUpdateMatrix(function2D<double>& MD, int Nk, int is, int ie, const nIntervals& interval)
{
  RemoveDetRatio(MD,Nk,is,ie);
  //// ** slow equivalent **
  // for (int i=0; i<Nk; i++)
  //   for (int j=0; j<Nk; j++){
  //     if (i!=ie && j!=is)
  // 	MD(i,j) -= MD(i,is)*MD(ie,j)/MD(ie,is);
  //   }
  double mc = MD(ie,is);
  for (int i=0; i<Nk; i++){
    if (i==ie) continue;
    double mdc = MD(i,is)/mc;
    const double* __restrict__ _MD_ie = MD[ie].MemPt();
    double* __restrict__       _MD_i  = MD[i].MemPt();
    for (int j=0; j<is; j++)    _MD_i[j] -= _MD_ie[j]*mdc;
    for (int j=is+1; j<Nk; j++) _MD_i[j] -= _MD_ie[j]*mdc;
  }
  // ** slow equivalent **
  // for (int i=ie; i<Nk-1; i++)
  //   for (int j=is; j<Nk-1; j++)
  //     MD(i,j) = MD(i+1,j+1);   
  for (int i=ie; i<Nk-1; i++)
    if(Nk-1>is) memcpy( MD[i].MemPt()+is, MD[i+1].MemPt()+(is+1), (Nk-1-is)*sizeof(double));
  
  // ** slow equivalent **
  // for (int i=0; i<ie; i++)
  //   for (int j=is; j<Nk-1; j++)
  //     MD(i,j) = MD(i,j+1);
  for (int i=0; i<ie; i++){
    double* __restrict__ _MD_i = MD[i].MemPt();
    for (int j=is; j<Nk-1; j++) _MD_i[j] = _MD_i[j+1];
  }

  // ** slow equivalent **
  // for (int i=ie; i<Nk-1; i++)
  //   for (int j=0; j<is; j++)
  //      MD(i,j) = MD(i+1,j);
  for (int i=ie; i<Nk-1; i++)
    if (is>0) memcpy( MD[i].MemPt(), MD[i+1].MemPt(), sizeof(double)*is);
  
  MD.resize(Nk-1,Nk-1);

  if (!common::Qsvd && common::QHB2){
    for (int j=is; j<Nk-1; j++)
      for (int im=0; im<piom->size(); im++)
	Ff(j,im) = Ff(j+1,im);
    Ff.resize(Nk-1,Ff.size_Nd());
  }
  
  return ratio*(1-2*((is+ie)%2));
}

inline void MatrixM::AddUpdate_Gf(const nIntervals& interval, vector<function2D<dcomplex> >& Mv, int istep)
{// This is currently the slowest routine. You should try to use intel-vml routines
  //if (common::Qsvd) return; // Nothing done when projecting to SVD functions.
  
  // Updates Green's function
  sum1=0; sum2=0;
  double sbta = 1/sqrt(beta);
  int nomega = piom->size();
  int dim = ptfl_index->size_N();
  for (int it=0; it<Lt.size(); it++){
    const dcomplex* __restrict__ _exp_e = interval.exp_e(it).MemPt();
    const dcomplex* __restrict__ _exp_s = interval.exp_s(it).MemPt();
    
    int be = interval.btype_e(it);
    int bs = interval.btype_s(it);
    dcomplex* __restrict__ _sum1 = sum1[be].MemPt();
    dcomplex* __restrict__ _sum2 = sum2[bs].MemPt();
    double lt = Lt[it]*sbta;
    double rt = Rt[it]*sbta;
    dcomplex* __restrict__ _sumt = sumt.MemPt();
    //xaxpy(nomega, lt, _exp_e, _sum1);
    for (int im=0; im<nomega; im++){
      _sum1[im] += _exp_e[im]*lt;
      dcomplex z = _exp_s[im]*rt;
      _sumt[im] = z;
      _sum2[im] += z;
    }
    
    if (common::QHB2){
      dcomplex* __restrict__ _sum2c = sum2c[it].MemPt();
      memcpy(_sum2c, _sumt, 2*sizeof(double)*nomega);
    }
  }
  for (int be=0; be<dim; be++){
    for (int bs=0; bs<dim; bs++){
      dcomplex* __restrict__ _sum1 = sum1[be].MemPt();
      dcomplex* __restrict__ _sum2 = sum2[bs].MemPt();
      int ind = (*ptfl_index)[be][bs];
      dcomplex* __restrict__ _Gf = Gf[ind].MemPt();// can use zdotc
      for (int im=0; im<nomega; im++) _Gf[im] -= _sum1[im]*_sum2[im].conj();
    }
  }
  
  if (common::QHB2){
    for (int it=0; it<Lt.size(); it++){
      dcomplex* __restrict__ _sum1 = sum1[0].MemPt();
      dcomplex* __restrict__ _sum2c = sum2c[it].MemPt();
      dcomplex* __restrict__ _Ff = Ff[it].MemPt();
      for (int im=0; im<nomega; im++) _Ff[im] -= _sum1[im]*_sum2c[im].conj(); // works only when baths are not multidimensional
    }
  }

  if (common::cmp_vertex){
    for (int be=0; be<dim; be++){
      for (int bs=0; bs<dim; bs++){
	const funProxy<dcomplex>& _sum1 = sum1[be];
	const funProxy<dcomplex>& _sum2 = sum2[bs];
	int ind = (*ptfl_index)[be][bs];
	int ind2 = v2fl_index[ind];
	function2D<dcomplex>& _Mv = Mv[ind2];
	for (int im1=0; im1<nomv; im1++){
	  dcomplex sm1 = _sum1[im1];
	  funProxy<dcomplex>& __Mv = _Mv[im1+nomv];
	  for (int im2=0; im2<nomv; im2++)  __Mv[im2+nomv] -= sm1*_sum2[im2].conj(); //_Mv(im1+nomv,im2+nomv) -= sm1*_sum2[im2];
	  for (int im2=-nomv; im2<0; im2++) __Mv[im2+nomv] -= sm1*_sum2[-im2-1]; //_Mv(im1+nomv,im2+nomv) -= sm1*_sum2[-im2-1].conj();
	}
      }
    }
  }

  if (common::QHB2 && istep!=-1){
    for (int im=0; im<piom->size(); im++){
      dcomplex csum=0;
      for (int i=0; i<Lt.size(); i++) csum += Ff(i,im);
      if (abs(Gf(0,im)-csum)>1e-6){
	cout<<"Green's function and F are different (add)"<<Gf(0,im)-csum<<endl;
      }
    }
  }
}

inline void MatrixM::MoveUpdate_Gf(const nIntervals& interval, vector<function2D<dcomplex> >& Mv, int istep)
{
  //if (common::Qsvd) return; // Nothing done when projecting to SVD functions.
  
  AddUpdate_Gf(interval, Mv, -1);
  
  for (int ind=0; ind<Gf.size_N(); ind++)
    for (int im=0; im<piom->size(); im++) Gf[ind][im] += dG[ind][im];
  
  if (common::cmp_vertex){
    for (int ind=0; ind<Gf.size_N(); ind++){
      int ind2 = v2fl_index[ind];
      function2D<dcomplex>& _Mv = Mv[ind2];
      const function2D<dcomplex>& _dMv = dMv[ind];
      for (int im1=0; im1<nomv; im1++){
	funProxy<dcomplex>& __Mv = _Mv[im1+nomv];
	const funProxy<dcomplex>& __dMv = _dMv[im1+nomv];
	
	for (int im2=0; im2<2*nomv; im2++) __Mv[im2] += __dMv[im2];
      }
    }
  }

  if (common::QHB2){
    for (int im=0; im<piom->size(); im++){
      dcomplex csum=0;
      for (int i=0; i<Lt.size(); i++) csum += Ff(i,im);
      if (abs(Gf(0,im)-csum)>1e-6) cout<<"Green's function and F are different (move)"<<Gf(0,im)-csum<<endl;
    }
  }
}

inline void MatrixM::RemoveUpdate_Gf(const function2D<double>& MD, int Nk, int is, int ie, const nIntervals& interval, vector<function2D<dcomplex> >& Mv)
{
  //if (common::Qsvd) return; // Nothing done when projecting to SVD functions.
  // Updates Green's function
  sum1=0; sum2=0;
  double sbta = 1/sqrt(beta);
  int nomega = piom->size();
  int dim = ptfl_index->size_N();
  for (int it=0; it<Nk; it++){
    //const funProxy<dcomplex>& exp_e = interval.exp_e(it);
    //const funProxy<dcomplex>& exp_s = interval.exp_s(it);
    const dcomplex* __restrict__ exp_e = interval.exp_e(it).MemPt();
    const dcomplex* __restrict__ exp_s = interval.exp_s(it).MemPt();
    int be = interval.btype_e(it);
    int bs = interval.btype_s(it);
    double me = MD(it,is)*sbta;
    double ms = MD(ie,it)*sbta;
    dcomplex* __restrict__ _sum1 = sum1[be].MemPt();
    dcomplex* __restrict__ _sum2 = sum2[bs].MemPt();
    dcomplex* __restrict__ _sumt = sumt.MemPt();
    for (int im=0; im<nomega; im++){
      _sum1[im] += exp_e[im]*me;
      dcomplex z = exp_s[im]*ms;
      _sumt[im]=z;
      _sum2[im] += z;
    }
    //for (int im=0; im<nomega; im++) _sum2[im] += sumt[im];
    if (common::QHB2){
      dcomplex* __restrict__ _sum2c = sum2c[it].MemPt();
      memcpy(sum2c[it].MemPt(),_sumt,2*sizeof(double)*nomega);
    }
  }
  for (int be=0; be<dim; be++){
    for (int bs=0; bs<dim; bs++){
      int ind = (*ptfl_index)[be][bs];  
      double mdi = 1/MD(ie,is);
      dcomplex* __restrict__ _sum1 = sum1[be].MemPt();
      dcomplex* __restrict__ _sum2 = sum2[bs].MemPt();
      dcomplex* __restrict__ _Gf = Gf[ind].MemPt();
      for (int im=0; im<piom->size(); im++)
	_Gf[im] += _sum1[im]*_sum2[im].conj()*mdi;
    }
  }
  if (common::QHB2){
    double mdi = 1/MD(ie,is);
    for (int it=0; it<Nk; it++){
      dcomplex* __restrict__ _sum1 = sum1[0].MemPt();
      dcomplex* __restrict__ _sum2c = sum2c[it].MemPt();
      dcomplex* __restrict__ _Ff = Ff[it].MemPt();
      for (int im=0; im<piom->size(); im++) _Ff[im] += _sum2c[im].conj()*_sum1[im]*mdi;// works only when baths are not multidimensional
    }
  }
  if (common::cmp_vertex){
    for (int be=0; be<ptfl_index->size_N(); be++){
      for (int bs=0; bs<ptfl_index->size_Nd(); bs++){
	const funProxy<dcomplex>& _sum1 = sum1[be];
	const funProxy<dcomplex>& _sum2 = sum2[bs];
	int ind = (*ptfl_index)[be][bs];
	int ind2 = v2fl_index[ind];
	double mdi = 1/MD(ie,is);
	function2D<dcomplex>& _Mv = Mv[ind2];
	for (int im1=0; im1<nomv; im1++){
	  funProxy<dcomplex>& __Mv = _Mv[im1+nomv];
	  dcomplex md1 = mdi*_sum1[im1];
	  for (int im2=0; im2<nomv; im2++)  __Mv[im2+nomv] += md1*_sum2[im2].conj();
	  for (int im2=-nomv; im2<0; im2++) __Mv[im2+nomv] += md1*_sum2[-im2-1];
	}
      }
    }
  }

  if (common::QHB2){
    for (int im=0; im<piom->size(); im++){
      dcomplex csum=0;
      for (int i=0; i<Nk; i++) csum += Ff(i,im);
      if (abs(Gf(0,im)-csum)>1e-6) cout<<"Green's function and F are different (remove)"<<Gf(0,im)<<" "<<csum<<endl;
    }
  }
  
}

void MatrixM::SetMatrix_(int Nk, const nIntervals& interval, function2D<double>& mD, const function2D<spline1D<double>* >& Delta_any)
{
  mD.resize(Nk,Nk);
  for (int i=0; i<Nk; i++){
    double t_start = interval.time_s(i);
    int bs = interval.btype_s(i);
    tint pos1 = ptau->InitInterpRight(), pos2 = ptau->InitInterpRight();
    for (int j=0; j<Nk; j++){
      double dt = t_start-interval.time_e(j);
      int be = interval.btype_e(j);
      if (dt<0){
	mD(i,j) = -Interp_Any_Delta(*Delta_any[bs][be], ptau->InterpRight(dt+beta,pos1), bs==be);
      }else{
	mD(i,j) =  Interp_Any_Delta(*Delta_any[bs][be], ptau->InterpRight(dt,pos2), bs==be);
      }
    }
  }
}

inline void MatrixM::CleanUpdateMatrix(function2D<double>& MD, int Nk, const nIntervals& interval, function2D<double>& mD)
{
  SetMatrix_(Nk, interval, mD, vfun);
  Inverse(mD, MD, ipiv);
}

inline void MatrixM::CleanUpdateGf(function2D<double>& MD, int Nk, const nIntervals& interval, vector<function2D<dcomplex> >& Mv, int istep)
{
  //if (common::Qsvd) return; // Nothing done when projecting to SVD functions.
  sum3=0;
  if (common::cmp_vertex)
    for (unsigned ind=0; ind<dMv.size(); ind++) dMv[ind]=0;

  if (common::QHB2){
    Ff.resize(Nk,piom->size());
    Ff=0.0;
  }
  for (int ie=0; ie<Nk; ie++){
    const funProxy<dcomplex>& exp_e = interval.exp_e(ie);
    int be = interval.btype_e(ie);
    for (int is=0; is<Nk; is++){
      const funProxy<dcomplex>& exp_s = interval.exp_s(is);
      int bs = interval.btype_s(is);      
      double md = -MD(ie,is)/beta;
      int ind = (*ptfl_index)[be][bs];
      dcomplex sm1;
      for (int im=0; im<piom->size(); im++){
	sm1 = exp_e[im]*md*conj(exp_s[im]);
	sum3[ind][im] += sm1;
	if (common::QHB2) Ff(is,im) += sm1;
      }
  
      if (common::cmp_vertex){
	function2D<dcomplex>& _dMv = dMv[ind];
	for (int im1=0; im1<nomv; im1++){
	  dcomplex exp_e_md = exp_e[im1]*md;
	  for (int im2=0; im2<nomv; im2++)  _dMv(im1+nomv,im2+nomv) += exp_e_md*conj(exp_s[im2]);
	  for (int im2=-nomv; im2<0; im2++) _dMv(im1+nomv,im2+nomv) += exp_e_md*exp_s[-im2-1];
	}
      }
    }
  }

  for (int ind=0; ind<Gf.size_N(); ind++)
    for (int im=0; im<piom->size(); im++) Gf[ind][im] = sum3[ind][im];
  
  if (common::cmp_vertex){
    for (int ind=0; ind<Gf.size_N(); ind++){
      int ind2 = v2fl_index[ind];
      function2D<dcomplex>& _Mv  = Mv[ind2];
      const function2D<dcomplex>& _dMv = dMv[ind];
      for (int im1=0; im1<nomv; im1++)
	for (int im2=-nomv; im2<nomv; im2++)  _Mv(im1+nomv,im2+nomv) = _dMv(im1+nomv,im2+nomv);
    }
  }
}

inline void MatrixM::AddComputeGs(int Nk, double t_start, int btype_s, double t_end, int btype_e, const nIntervals& interval,
				  function1D<double>& d_ln, function1D<double>& d_nl, double& d_nn)
{
  tint pos1 = ptau->InitInterpLeft();
  tint pos2 = ptau->InitInterpLeft();
  for (int i=0; i<Nk; i++){
    int bs = interval.btype_s(i);
    int be = btype_e;
    if (interval.time_s(i)<t_end)
      d_ln[i] = -Interp_Delta(bs, be, ptau->InterpLeft(interval.time_s(i)-t_end+beta,pos1));
    else
      d_ln[i] = Interp_Delta(bs, be, ptau->InterpLeft(interval.time_s(i)-t_end,pos2));
  }
  pos1 = ptau->InitInterpRight();
  pos2 = ptau->InitInterpRight();
  for (int i=0; i<Nk; i++){
    int bs = btype_s; 
    int be = interval.btype_e(i);
    if (t_start<interval.time_e(i))
      d_nl[i] = -Interp_Delta(bs, be, ptau->InterpRight(t_start-interval.time_e(i)+beta,pos1));
    else
      d_nl[i] = Interp_Delta(bs, be, ptau->InterpRight(t_start-interval.time_e(i),pos2));
  }
  d_nn = Interp_Antiperiodic_Delta(btype_s, btype_e, t_start-t_end);
}

inline double MatrixM::AddDetRatio(function2D<double>& MD, int Nk,
				   int btype_s1, double t_start1, int btype_e1, double t_end1, 
				   int btype_s2, double t_start2, int btype_e2, double t_end2, 
				   const nIntervals& interval)
{
  double d_nn;
  double d_mm;
  AddComputeGs(Nk, t_start1, btype_s1, t_end1, btype_e1, interval, d_ln, d_nl, d_nn);
  AddComputeGs(Nk, t_start2, btype_s2, t_end2, btype_e2, interval, d_lm, d_ml, d_mm);
  double d_nm = Interp_Antiperiodic_Delta(btype_s1, btype_e2, t_start1-t_end2);
  double d_mn = Interp_Antiperiodic_Delta(btype_s2, btype_e1, t_start2-t_end1);
  
  for (int i=0; i<Nk; i++){
    double sumL1=0;
    double sumR1=0;
    double sumL2=0;
    for (int l=0; l<Nk; l++){
      sumL1 += MD(i,l)*d_ln[l];
      sumR1 += d_nl[l]*MD(l,i);
      sumL2 += MD(i,l)*d_lm[l];
    }
    Lt[i] = sumL1;
    Rt[i] = sumR1;
    Lt2[i] = sumL2;
  }
  double sumL_mn=0;
  double sumR_nm=0;
  double sumL_nn=0;
  double sumL_mm=0;
  for (int l=0; l<Nk; l++){
    sumL_mn += d_ml[l]*Lt[l];
    sumR_nm += Rt[l]*d_lm[l];
    sumL_nn += d_nl[l]*Lt[l];
    sumL_mm += d_ml[l]*Lt2[l];
  }
  ratio = (d_nn - sumL_nn)*(d_mm - sumL_mm) - (sumL_mn - d_mn)*(sumR_nm - d_nm);
  return ratio;
}

double MatrixM::GlobalFlipDetRatio(const nIntervals& interval, const function2D<double>& MD, const function2D<spline1D<double>* >& Delta_flipped, /*function2D<double>& C,*/ function2D<double>& mD)
{
  int Nk = interval.size()/2;
  if (Nk<=0) return 1;
  SetMatrix_(Nk, interval, mD, Delta_flipped);
  if (mD.size_N()!=MD.size_N() || mD.size_Nd()!=MD.size_Nd()){
    cerr<<"WARNING: Sizes of mD and MD are different"<<endl;
  }
  //C.Product("N","N",mD,MD);
  //return Det(C);
  //double dd1 = Det(mD)*Det(MD);
  double dd2 = (logDet(mD)*logDet(MD)).dbl(); // ratio of determinants
  //if (fabs(dd1/dd2-1.0)>1e-5) cerr<<"The determinants computed in two different ways are different "<<dd1<<" "<<dd2<<" "<<endl;
  return dd2;
}


void MatrixM::ComputeGtau(const function2D<double>& MD, const nIntervals& interval, function2D<double>& Gtau)
{
  // Gtau=0;
  int Nk = interval.size()/2;
  int Gtsize = Gtau[0].size();
  for (int ie=0; ie<Nk; ie++){
    double te = interval.time_e(ie);
    int be = interval.btype_e(ie);
    for (int is=0; is<Nk; is++){
      double ts = interval.time_s(is);
      int bs = interval.btype_s(is);
      double md = -MD(ie,is)/beta;
      int idt;
      if (te-ts>0){
	idt = static_cast<int>((te-ts)/beta*Gtsize);
      }else{
	idt = static_cast<int>((beta+te-ts)/beta*Gtsize);
	md = -md;
      }
      if (idt>=Gtau[0].size()) idt = Gtau[0].size()-1;
      //if (idt<0) idt=0;
      int ind = (*ptfl_index)[be][bs];
      Gtau[ind][idt] += md;
    }
  }
}
void MatrixM::ComputeG_svd_debug(function2D<double>& gl, function2D<double>& Gtau, double tsign, const SVDFunc& svd, const function2D<double>& MD, const nIntervals& interval)
{ // Computes the Green's function coefficients g_l for the current distribution of kinks.
  // The formula is
  //    G[(be,bs)][tau] = -1/beta * \sum_{t_e,t_s} M[(be,bs)](t_e,t_s) * delta^{-}(tau-(t_e-t_s))
  //      where delta^{-}(tau) = \delta(tau) when tau>0 and is -\delta(tau+beta) when tau<0.
  //    When expanding in terms of u_l(\tau) basis functions, we get
  //    G[(be,bs)][l] = -1/beta * \sum_{t_e,t_s} M[(be,bs)](t_e,t_s) * u_l(t_e-t_s)   if tau>0
  //     and 
  //    G[(be,bs)][l] = -1/beta * \sum_{t_e,t_s} M[(be,bs)](t_e,t_s) * (-1)*u_l(t_e-t_s+beta)   if tau<0
  //
  int Nk = interval.size()/2;
  int Gtsize = Gtau[0].size();
  for (int ie=0; ie<Nk; ie++){
    double te = interval.time_e(ie);
    int be = interval.btype_e(ie);
    for (int is=0; is<Nk; is++){
      double ts = interval.time_s(is);
      int bs = interval.btype_s(is);
      int ind = (*ptfl_index)(be,bs);
      double md = -MD(ie,is) /beta *(tsign);
      double t = te-ts;
      intpar p;
      int idt;
      if (t>=0){
	p = svd.tau.Interp(t);
	idt = static_cast<int>(t/beta*Gtsize);
      }else{
	p = svd.tau.Interp(t+beta);
	idt = static_cast<int>((t+beta)/beta*Gtsize);
	md = -md;
      }
      //*** slower equivalent is ****//
      //for (int l=0; l<svd.lmax; l++) gl(ind,l) += md*svd.fU[l](p);
      svd.FastInterp(gl[ind], p, md);

      if (idt>=Gtsize) idt = Gtsize-1;
      //if (idt<0) idt=0;
      Gtau[ind][idt] += md;
    }
  }
}

void MatrixM::ComputeG_svd(function2D<double>& gl, double tsign, const SVDFunc& svd, const function2D<double>& MD, const nIntervals& interval)
{ // Computes the Green's function coefficients g_l for the current distribution of kinks.
  // The formula is
  //    G[(be,bs)][tau] = -1/beta * \sum_{t_e,t_s} M[(be,bs)](t_e,t_s) * delta^{-}(tau-(t_e-t_s))
  //      where delta^{-}(tau) = \delta(tau) when tau>0 and is -\delta(tau+beta) when tau<0.
  //    When expanding in terms of u_l(\tau) basis functions, we get
  //    G[(be,bs)][l] = -1/beta * \sum_{t_e,t_s} M[(be,bs)](t_e,t_s) * u_l(t_e-t_s)   if tau>0
  //     and 
  //    G[(be,bs)][l] = -1/beta * \sum_{t_e,t_s} M[(be,bs)](t_e,t_s) * (-1)*u_l(t_e-t_s+beta)   if tau<0
  //
  int Nk = interval.size()/2;
  for (int ie=0; ie<Nk; ie++){
    double te = interval.time_e(ie);
    int be = interval.btype_e(ie);
    for (int is=0; is<Nk; is++){
      double ts = interval.time_s(is);
      int bs = interval.btype_s(is);
      int ind = (*ptfl_index)(be,bs);
      double md = -MD(ie,is) /beta *(tsign);
      double t = te-ts;
      intpar p; 
      if (t>=0){
	p = svd.tau.Interp(t);
      }else{
	p = svd.tau.Interp(t+beta);
	md = -md;
      }
      //*** slower equivalent is ****//
      //for (int l=0; l<svd.lmax; l++) gl(ind,l) += md*svd.fU[l](p);
      svd.FastInterp(gl[ind], p, md);
    }
  }
}

void MatrixM::ComputeFG_svd(vector<function2D<double> >& gl, double tsign, const SVDFunc& svd, const function2D<double>& MD, const nIntervals& interval)
{ // This is very similar to the above Green's function calculation, except that it keeps also information about the starting kinks, i.e. t_s.
  // The formula is
  //    G[(be,bs)][t_s,tau] = -1/beta * \sum_{t_s} M[(be,bs)](t_e,t_s) * delta^{-}(tau-(t_e-t_s))
  //      where delta^{-}(tau) = \delta(tau) when tau>0 and is -\delta(tau+beta) when tau<0.
  //    When expanding in terms of u_l(\tau) basis functions, we get
  //    G[(be,bs)][t_s,l] = -1/beta * \sum_{t_e} M[(be,bs)](t_e,t_s) * u_l(t_e-t_s)   if tau>0
  //     and 
  //    G[(be,bs)][t_s,l] = -1/beta * \sum_{t_e} M[(be,bs)](t_e,t_s) * (-1)*u_l(t_e-t_s+beta)  if tau<0
  //
  int Nk = interval.size()/2;
  for (int ie=0; ie<Nk; ie++){
    double te = interval.time_e(ie);
    int be = interval.btype_e(ie);
    for (int is=0; is<Nk; is++){
      double ts = interval.time_s(is);
      int bs = interval.btype_s(is);
      int ind = (*ptfl_index)(be,bs);
      double md = -MD(ie,is) * tsign/beta;
      double t = te-ts;
      intpar p; 
      if (t>=0){
	p = svd.tau.Interp(t);
      }else{
	p = svd.tau.Interp(t+beta);
	md = -md;
      }
      //*** slower equivalent is ****//
      //for (int l=0; l<svd.lmax; l++) gl(ind,l) += md*svd.fU[l](p);
      
      svd.FastInterp( gl[ind][is], p, md);
      /*
      svd.FastInterp_( M_ul, p, md);
      gl[ind][is] += M_ul;
      */
    }
  }
}

void MatrixM::ComputeFGF_svd(vector<function2D<double> >& gs, vector<function2D<double> >& ge, vector<function3D<double> >& ges, bool Vertex_subtract, double tsign, const SVDFunc& svd, const function2D<double>& MD, const nIntervals& interval)
{ // This is very similar to the above two routines for the Green's function calculation, except that it
  // computes two arrays rather than one, in which one keeps also information about the starting kinks (t_s),
  // and the other keeps the information about the ending kinks (t_e).
  // The formula is
  //    Gs[(be,bs)][t_s,tau] = -1/beta * \sum_{t_e} M[(be,bs)](t_e,t_s) * delta^{-}(tau-(t_e-t_s))
  //    Ge[(be,bs)][t_e,tau] = -1/beta * \sum_{t_s} M[(be,bs)](t_e,t_s) * delta^{-}(tau-(t_e-t_s))
  //      where delta^{-}(tau) = \delta(tau) when tau>0 and is -\delta(tau+beta) when tau<0.
  //    When expanding in terms of u_l(\tau) basis functions, we get
  //    Gs[(be,bs)][t_s,l] = -1/beta * \sum_{t_e} M[(be,bs)](t_e,t_s) * u_l(t_e-t_s)   if tau>0
  //    Ge[(be,bs)][t_e,l] = -1/beta * \sum_{t_s} M[(be,bs)](t_e,t_s) * u_l(t_e-t_s)   if tau>0
  //     and 
  //    Gs[(be,bs)][t_s,l] = -1/beta * \sum_{t_e} M[(be,bs)](t_e,t_s) * (-1)*u_l(t_e-t_s+beta)  if tau<0
  //    Ge[(be,bs)][t_e,l] = -1/beta * \sum_{t_e} M[(be,bs)](t_e,t_s) * (-1)*u_l(t_e-t_s+beta)  if tau<0
  //
  double tsign_beta = tsign/beta;
  int Nk = interval.size()/2;
  for (int ie=0; ie<Nk; ie++){
    double te = interval.time_e(ie);
    int be = interval.btype_e(ie);
    for (int is=0; is<Nk; is++){
      double ts = interval.time_s(is);
      int bs = interval.btype_s(is);
      int ind = (*ptfl_index)(be,bs);
      int ind2 = v2fl_index[ind];
      double md = -MD(ie,is) * tsign_beta;
      double t = te-ts;
      intpar p; 
      if (t>=0){
	p = svd.tau.Interp(t);
      }else{
	p = svd.tau.Interp(t+beta);
	md = -md;
      }
      //*** slower equivalent is ****//
      //for (int l=0; l<svd.lmax; l++) gl(ind,l) += md*svd.fU[l](p);
      svd.FastInterp_( M_ul, p, md);
      gs[ind2][is] += M_ul;
      ge[ind2][ie] += M_ul;
      if (Vertex_subtract){
	for (int l=0; l<svd.lmax; l++) ges[ind2](ie,is,l) = M_ul[l];
      }
    }
  }
}


void swapGf(MatrixM& m1, MatrixM& m2, vector<function2D<dcomplex> >& Mv)
{
  //if (common::Qsvd) return; // Nothing done when projecting to SVD functions.
  static function2D<dcomplex> Gtmp;
  Gtmp = m1.Gf;
  m1.Gf = m2.Gf;
  m2.Gf = Gtmp;

  static function2D<dcomplex> Ftmp;
  Ftmp = m1.Ff;
  m1.Ff = m2.Ff;
  m2.Ff = Ftmp;
 
  if (common::cmp_vertex){
    for (int ind=0; ind<m1.Gf.size_N(); ind++){
      int inda = m1.v2fl_index[ind];
      int indb = m2.v2fl_index[ind];
      //swap Mv[a] with Mv[b]
      m1.dMv[0] = Mv[inda];
      Mv[inda] = Mv[indb];
      Mv[indb] = m1.dMv[0];
    }
  }
}
