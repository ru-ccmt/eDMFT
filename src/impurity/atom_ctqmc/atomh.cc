// @Copyright 2007 Kristjan Haule
// 
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <deque>
#include <sstream>
#include "util.h"
#include "complex.h"
#include "blas.h"
#include "function.h"

using namespace std;

vector<double> ComputeCoulombRatios(int l)
{
  // Ratio between F2,F4,F6 and J! At the end of the day, we want to use U and J only!
  vector<double> FkoJ(l+1);
  FkoJ[0]=0;//This is not defined
  if (l==1) FkoJ[1] = 5.;
  if (l==2) { FkoJ[1]=14./(1+0.625); FkoJ[2]=0.625*FkoJ[1];}
  if (l==3) { FkoJ[1]=6435./(286+195*0.668+250*0.494); FkoJ[2] = 0.668*FkoJ[1]; FkoJ[3]=0.494*FkoJ[1];}
  return FkoJ;
}

// Just old good !
int Binomial(int n, int m)
{
  double Mf = 1;
  for (int i=2; i<=m; i++) Mf *= i;
  double r = 1;
  for (int i=n; i>=n-m+1; i--) r*=i;
  return static_cast<int>(r/Mf);
}
double dFactorial(int j)
{
  if (j<0) {cerr<<"Factorial defined only for non-negative numbers!"<<endl;}
  if (j<=1) return 1;
  double x=1;
  for (int i=2; i<=j; i++) x*=i;
  return x;
}
double dFactorial(double x)
{
  double r=1;
  double y=x;
  while(y>1.0){
    r *= y;
    y--;
  }
  if (fabs(y-0.5)<1e-10) r/=M_2_SQRTPI;
  return r;
}
template <class T>
class svector : public std::vector<T>{
  int start_, stop_;
public:
  svector(int _start, int _stop) : std::vector<T>(_stop-_start), start_(_start), stop_(_stop){};
  T operator[](int i) const{return std::vector<T>::operator[](i-start_);}
  T& operator[](int i){return std::vector<T>::operator[](i-start_);}
  int start(){return start_;}
  int stop(){return stop_;}
};


class ClebschGordan{
  int nn, jn;
  map<int,double> stored;
public:
  int mone(int i) const {return 1-2*(i&1);}
  int encode(double j, double m)
  { return static_cast<int>(2*j*(j+1)+m+1e-6);}
  void decode(int i, double& j, double& m)
  {
    double kk = 0.25*(-1+sqrt(1.+8*(i+1)));
    j = 0.5*floor(2*(kk-0.0001));
    m = i - 2*j*(j+1);
  }
  int encode(double j, double m, double j1, double m1, double j2, double m2)
  {
    int i = encode(j,m);
    int i1 = encode(j1,m1);
    int i2 = encode(j2,m2);
    int ii = i + nn*(i1+nn*i2);
    return encode(j2,m2) + nn*(encode(j1,m1)+nn*encode(j,m));
  }
  void decode(int i, double& j, double& m, double& j1, double& m1, double& j2, double& m2)
  {
    int j2m2 = i%nn;
    int i1 = i/nn;
    int j1m1 = i1%nn;
    int jm = i1/nn;
    decode(jm,j,m);
    decode(j1m1, j1, m1);
    decode(j2m2, j2, m2);
  }
public:
  ClebschGordan() : nn(1000), jn(22) {};
  double CG(double j, double m, double j1, double m1, double j2, double m2) const
  {
    if (m1+m2 != m) return 0;
    if (fabs(m2)>j2 || fabs(m1)>j1 || fabs(m)>j) return 0;
    int tmin = static_cast<int>(max(max(0.0,j2-j-m1),j1-j+m2) + 1e-14);
    int tmax = static_cast<int>(min(min(j1+j2-j,j1-m1),j2+m2) + 1e-14);
    double sum=0;
    for (int t=tmin; t<=tmax; t++){
      double v1 = sqrt((2*j+1)*dFactorial(j1+m1)*dFactorial(j1-m1)*dFactorial(j2+m2)*dFactorial(j2-m2)*dFactorial(j+m)*dFactorial(j-m));
      double v2 = dFactorial(t)*dFactorial(j1+j2-j-t)*dFactorial(j1-m1-t)*dFactorial(j2+m2-t)*dFactorial(j-j2+m1+t)*dFactorial(j-j1-m2+t);
      sum += mone(t)*v1/v2;
    }
    return sum*Delta(j1,j2,j);
  }
  double CGs(double j, double m, double j1, double m1, double j2, double m2)
  {
    if (j>jn || j1>jn || j2>jn) return CG(j,m,j1,m1,j2,m2);
    int ii = encode(j,m,j1,m1,j2,m2);
    map<int,double>::iterator il;
    il = stored.find(ii);
    if (il!=stored.end()) return il->second;
    double v = CG(j,m,j1,m1,j2,m2);
    stored[ii]=v;
    return v;
  }
  double Delta(double j1, double j2, double j) const
  {return sqrt(dFactorial(j1+j2-j)*dFactorial(j1-j2+j)*dFactorial(-j1+j2+j)/dFactorial(j1+j2+j+1));}
  double f6j(double j1, double j2, double j3, double m1, double m2, double m3)
  {
    int tmin = static_cast<int>(max(max(max(j1+j2+j3,j1+m2+m3),m1+j2+m3),m1+m2+j3)+1e-14);
    int tmax = static_cast<int>(min(min(j1+j2+m1+m2,j1+j3+m1+m3),j2+j3+m2+m3)+1e-14);
    double sum=0;
    for (int t=tmin; t<=tmax; t++){
      double v1 = dFactorial(t-j1-j2-j3)*dFactorial(t-j1-m2-m3)*dFactorial(t-m1-j2-m3)*dFactorial(t-m1-m2-j3);
      double v2 = dFactorial(j1+j2+m1+m2-t)*dFactorial(j1+j3+m1+m3-t)*dFactorial(j2+j3+m2+m3-t);
      sum += mone(t)*dFactorial(t+1)/(v1*v2);
    }
    return Delta(j1,j2,j3)*Delta(j1,m2,m3)*Delta(m1,m2,j3)*sum;
  }
  double f3j(double j1, double m1, double j2, double m2, double j3, double m3) const
  {
    if (m1+m2+m3!=0) return 0;
    if (abs(j1-j2)-1e-14>j3 || j3>j1+j2+1e-14) return 0;
    if (abs(m1)>j1 || abs(m2)>j2 || abs(m3)>j3) return 0;
    int tmin = static_cast<int>(max(max(0.0,j2-j3-m1),j1-j3+m2)+1e-14);
    int tmax = static_cast<int>(min(min(j1+j2-j3,j1-m1),j2+m2)+1e-14);
    double sum=0;
    for (int t=tmin; t<=tmax; t++){
      double v1 = dFactorial(j3-j2+m1+t)*dFactorial(j3-j1-m2+t);
      double v2 = dFactorial(j1+j2-j3-t)*dFactorial(j1-m1-t)*dFactorial(j2+m2-t);
      sum += mone(t)/(dFactorial(t)*v1*v2);
    }
    double dn = dFactorial(j1+m1)*dFactorial(j1-m1)*dFactorial(j2+m2)*dFactorial(j2-m2)*dFactorial(j3+m3)*dFactorial(j3-m3);
    return mone(static_cast<int>(j1-j2-m3))*Delta(j1,j2,j3)*sqrt(dn)*sum;
  }
  double Gaunt(int l1, int m1, int l2, int m2, int l3, int m3) const
  {
    // Calculates <Y_{l1m1}|Y_{l2m2}|Y_{l3m3}>
    if (l1<0 || l2<0 || l3<0) {cerr<<"Quantum number l must be non-negative!"<<endl;}
    //    cout<<sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/(4*M_PI))<<" "<<f3j(l1,0,l2,0,l3,0)<<" "<<f3j(l1,-m1,l2,m2,l3,m3)<<endl;
    return mone(m1)*sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/(4*M_PI))*f3j(l1,0,l2,0,l3,0)*f3j(l1,-m1,l2,m2,l3,m3);
  }
};

class gaunt_ck{
  int l;
  vector<double> ck;
public:
  gaunt_ck(int l_, const ClebschGordan& cg);
  const double& operator()(int m1, int m2, int k) const{
    if (abs(m1)>l || abs(m2)>l || k<0 || k>2*l) { cerr<<"Gaunt coefficient does not exist!"<<endl; return ck[0];}
    int i =m1+l;
    int j =m2+l;
    int kp = k/2;
    int offset = kp + j*(l+1) + i*(l+1)*(2*l+1);
    return ck[offset];
  }
private:
  double& operator()(int m1, int m2, int k){
    if (abs(m1)>l || abs(m2)>l || k<0 || k>2*l) { cerr<<"Gaunt coefficient does not exist!"<<endl; return ck[0];}
    int i =m1+l;
    int j =m2+l;
    int kp = k/2;
    int offset = kp + j*(l+1) + i*(l+1)*(2*l+1);
    return ck[offset];
  }
};
gaunt_ck::gaunt_ck(int l_, const ClebschGordan& cg) : l(l_), ck((2*l+1)*(2*l+1)*(l+1))
{
  for (int m1=-l; m1<=l; m1++){
    for (int m2=-l; m2<=l; m2++){
      for (int k=0; k<=2*l; k+=2){
	double c = cg.Gaunt(l,m1,k,m1-m2,l,m2)*sqrt(4*M_PI/(2*k+1));
	if (abs(c)<1e-10) c=0;
	operator()(m1,m2,k) = c;
      }
    }
  }
}


void CoulombU0(function2D<function2D<double> >& Uf, int l, const vector<pair<int,double> >& ml_ms, const vector<pair<double,double> >& one_electron_base, const ClebschGordan& cg, const gaunt_ck& gck, const vector<double>& FkoJ)
{
  // We have U(m1,m2,m3,m4) psi_m1^+ * psi_m2^+ * psi_m3 * psi_m4
  vector<function2D<function2D<double> > > Uk(2);
  int baths = ml_ms.size();

  cout<<"ml_ms="<<endl;
  for (int i=0; i<baths; i++){
    cout<<i<<" "<<ml_ms[i].first<<" "<<ml_ms[i].second<<endl;
  }  
  for (int k=0; k<2; k++){
    Uk[k].resize(baths,baths);
    for (int i=0; i<baths; i++)
      for (int j=0; j<baths; j++)
	Uk[k](i,j).resize(baths,baths);
  }
  for (int i=0; i<baths; i++){
    int m1 = ml_ms[i].first;
    int s1 = round(ml_ms[i].second*2);
    for (int j=0; j<baths; j++){
      int m2 = ml_ms[j].first;
      int s2 = round(ml_ms[j].second*2);
      for (int a=0; a<baths; a++){
	int m3 = ml_ms[a].first;
	int s3 = round(ml_ms[a].second*2);
	for (int b=0; b<baths; b++){
	  int m4 = ml_ms[b].first;
	  int s4 = round(ml_ms[b].second*2);
	  if (m4-m1!=m2-m3) continue;
	  if (s1!=s4 || s2!=s3) continue;
	  //if (i==0 && a==0 && j==0) cout<<" b="<<b<<" s1="<<s1<<" s2="<<s2<<" s3="<<s3<<" s4="<<s4<<endl;
	    
	  //Uk[0](i,a)(j,b) = gck(m4,m1,0)*gck(m2,m3,0);
	  Uk[0](i,j)(a,b) = gck(m4,m1,0)*gck(m2,m3,0);
	  double sum=0;
	  for (int k=1; k<=l; k++) sum += gck(m4,m1,2*k)*gck(m2,m3,2*k)*FkoJ[k];
	  //Uk[1](i,a)(j,b) = sum;
	  Uk[1](i,j)(a,b) = sum;
	}
      }
    }
  }
  /*
  cout<<"U0 before transform"<<endl;
  for (int i=0; i<baths; i++){
    for (int j=0; j<baths; j++){
      for (int k=0; k<baths; k++){
	for (int l=0; l<baths; l++){
	  if (fabs(Uk[1](i,j)(k,l))>1e-6){
	    cout<<setw(2)<<i<<" "<<setw(2)<<j<<" "<<setw(2)<<k<<" "<<setw(2)<<l<<" "<<setw(20)<<Uk[1](i,j)(k,l)<<endl;
	  }
	}
      }
    }
  }
  */
  // U(m1,m3)(m2,m4)
  // Ts ((j,mj),(ml,ms))
  // Tst((ml,ms),(j,mj))
  function2D<double> Ts(baths,baths), Tst(baths,baths);
  for (int i=0; i<one_electron_base.size(); i++){
    for (int j=0; j<ml_ms.size(); j++){
      Ts(i,j) = cg.CG(one_electron_base[i].first,one_electron_base[i].second,l,ml_ms[j].first,0.5,ml_ms[j].second);
      Tst(j,i) = Ts(i,j);
    }
  }
  
  function2D<double> R(baths,baths), tmp(baths,baths);
  function2D<function2D<double> > UR(baths,baths);
  for (int i=0; i<baths; i++){
    for (int j=0; j<baths; j++){
      UR(i,j).resize(baths,baths);
      UR(i,j) = 0.0;
    }
  }
  
  for (int ms1=0; ms1<baths; ms1++){
    for (int ms2=0; ms2<baths; ms2++){
      //R = Ts * Uk[1](ms1,ms2) * Ts^+;
      tmp.MProduct(Ts,Uk[1](ms1,ms2));
      R.MProduct(tmp,Tst);
      for (int j3=0; j3<baths; j3++)
	for (int j4=0; j4<baths; j4++)
	  UR(j3,j4)(ms1,ms2) += R(j3,j4);
    }
  }
  /*
  cout<<"U1 after first transform"<<endl;
  for (int ms1=0; ms1<baths; ms1++){
    for (int ms2=0; ms2<baths; ms2++){
      for (int j3=0; j3<baths; j3++){
	for (int j4=0; j4<baths; j4++){
	  if (fabs(UR(j3,j4)(ms1,ms2))>1e-6){
	    cout<<setw(2)<<ms1<<" "<<setw(2)<<ms2<<" "<<setw(2)<<j3<<" "<<setw(2)<<j4<<" "<<setw(20)<<UR(j3,j4)(ms1,ms2)<<endl;
	  }
	}
      }
    }
  }
  */
  Uf.resize(baths,baths);
  for (int i=0; i<baths; i++)
    for (int j=0; j<baths; j++){
      Uf(i,j).resize(baths,baths);
      Uf(i,j) = 0.0;
    }
  
  function2D<double> Q(baths,baths);
  for (int j3=0; j3<baths; j3++){
    for (int j4=0; j4<baths; j4++){
      //Q = Ts * UR(j3,j4) * Ts^+;
      tmp.MProduct(Ts, UR(j3,j4));
      Q.MProduct(tmp,Tst);
      for (int j1=0; j1<baths; j1++)
	for (int j2=0; j2<baths; j2++)
	  Uf(j1,j2)(j3,j4) += Q(j1,j2);
    }
  }
  /*
  cout<<"U2 after second transform"<<endl;
  for (int j1=0; j1<baths; j1++){
    for (int j2=0; j2<baths; j2++){
      for (int j3=0; j3<baths; j3++){
	for (int j4=0; j4<baths; j4++){
	  if (fabs(Uf(j1,j2)(j3,j4))>1e-6){
	    cout<<setw(2)<<j1<<" "<<setw(2)<<j2<<" "<<setw(2)<<j3<<" "<<setw(2)<<j4<<" "<<setw(20)<<Uf(j1,j2)(j3,j4)<<endl;
	  }
	}
      }
    }
  }
  */
}

class operateLS{
  int l, baths, l2p1;
  vector<int> mask, mask_u, mask_d;
  vector<int> lz_u, lz_d, ind_lz_u, ind_lz_d;
  vector<int> lz, sz;
  vector<double> lp_fact, lm_fact;
  void startSort(vector<int>& base)
  {
    vector<int> Nn(baths+1);
    vector<int> Nmax(baths+1);
    for (int i=0; i<=baths; i++) Nn[i]=0;
    Nmax[0] = 0;
    for (int i=1; i<=baths; i++) Nmax[i] = Nmax[i-1] + Binomial(baths,i-1);
 
    for (int i=0; i<base.size(); i++){
      int n = Nel(i);
      base[Nmax[n]+Nn[n]] = i;
      Nn[n]++;
    }
  }
public:
  operateLS(vector<int>& base, vector<pair<int,double> > ml_ms, int l_) :
    l(l_), baths(2*(2*l+1)), l2p1(2*l+1), mask(baths), mask_u(l2p1), mask_d(l2p1), lz_u(l2p1), lz_d(l2p1),
    ind_lz_u(l2p1), ind_lz_d(l2p1), lz(baths), sz(baths), lp_fact(l2p1), lm_fact(l2p1)
  {
    for (int i=0; i<baths; i++) mask[i] = 1<<i;
    int iu=0, id=0, ib=0;
    for (int i=0; i<baths; i++){
      if (ml_ms[i].second>0) {lz_u[iu] = ml_ms[i].first; mask_u[iu] = mask[i]; iu++; lz[ib] = ml_ms[i].first; sz[ib] = static_cast<int>(2*ml_ms[i].second); ib++;}
      if (ml_ms[i].second<0) {lz_d[id] = ml_ms[i].first; mask_d[id] = mask[i]; id++; lz[ib] = ml_ms[i].first; sz[ib] = static_cast<int>(2*ml_ms[i].second); ib++;}
    }
    for (int i=0; i<iu; i++){
      ind_lz_u[lz_u[i]+l]=i;
      ind_lz_d[lz_d[i]+l]=i;
    }
    for (int i=0; i<l2p1; i++){
      int m = i-l;
      lp_fact[i] = sqrt(l*(l+1)-m*(m+1)+0.0);
      lm_fact[i] = sqrt(l*(l+1)-m*(m-1)+0.0);
    }
    startSort(base);
  }
  int estate(int state, int lz) const{
    bool lup = state & mask_u[ind_lz_u[lz+l]];
    bool ldo = state & mask_d[ind_lz_d[lz+l]];
    if (lup && ldo) return 2;
    if (lup) return 1;
    if (ldo) return -1;
    return 0;
  }
  int Nel(int state) const
  {
    int n=0;
    for (int k=0; k<baths; k++) if (mask[k]&state) n++;
    return n;
  }
  double Sz(int state) const
  {
    int nu=0, nd=0;
    for (int ilz=0; ilz<l2p1; ilz++){
      if (mask_u[ilz]&state) nu++;
      if (mask_d[ilz]&state) nd++;
    }
    return 0.5*(nu-nd);
  }
  int Lz(int state) const
  {
    int tLz=0;
    for (int ilz=0; ilz<l2p1; ilz++){
      if (mask_u[ilz] & state) tLz += lz_u[ilz];
      if (mask_d[ilz] & state) tLz += lz_d[ilz];
    }
    return tLz;
  }
  double Jz(int state) const
  { return Lz(state)+Sz(state);}
  void Sp(int state, list<int>& sts)
  {
    sts.clear();
    for (int ilz=0; ilz<l2p1; ilz++){
      if (mask_d[ilz]&state && !(mask_u[ilz]&state)){
	int newstate=state^mask_d[ilz]^mask_u[ilz];
	sts.push_back(newstate);
      }
    }
  }
  int sign(int state, int mask_min, int mask_max) const
  {
    int mask = mask_min<<1;
    int n=0;
    while (mask<mask_max){
      if (mask&state) n++;
      mask = mask<<1;
    }
    return 1-2*(n&1);
  }
  int sign(int state, int mask) const
  {
    int m = 1;
    int n=0;
    while (m<mask){
      if (m&state) n++;
      m = m<<1;
    }
    return 1-2*(n%2);
  }
  void S2(int state, list<pair<int,double> >& sts) const
  {
    sts.clear();
    // diagonal part
    double dd=0;
    for (int ilz=0; ilz<l2p1; ilz++){
      int up = (mask_u[ilz]&state) ? 1 : 0;
      int dn = (mask_d[ilz]&state) ? 1 : 0;
      // if only up or only down in certain lz
      if (up+dn==1) dd += 0.5;
    }
    // Sz^2
    double fct = sqr(Sz(state)) + dd;
    // store diagonal
    sts.push_back(make_pair(state,fct));
    // off diagonal
    for (int ilz=0; ilz<l2p1; ilz++){
      int im1 = mask_u[ilz];
      int im2 = mask_d[ilz];
      bool ib1 = state & im1;
      bool ib2 = state & im2;
      if (ib1 && !ib2){// S^-_i gives nonzero
	int isig = sign(state,min(im1,im2),max(im1,im2));
	int istate = state^im1^im2;
	for (int jlz=0; jlz<l2p1; jlz++){
	  if (ilz==jlz) continue;
	  int jm1 = mask_d[jlz];
	  int jm2 = mask_u[jlz];
	  bool jb1 = state & jm1;
	  bool jb2 = state & jm2;
	  if (jb1 && !jb2){//S^+_j gives nonzero
	    int jsig = sign(istate,min(jm1,jm2),max(jm1,jm2));
	    int jstate = istate^jm1^jm2;
	    sts.push_back(make_pair(jstate, isig*jsig));
	  }
	}
      }
    }
  }
  void L2_off_diag(int state, list<pair<int,double> >& sts,
		   const vector<int>& maski, const vector<int>& lzi,
		   const vector<int>& maskj, const vector<int>& lzj, bool cont) const
  {
    for (int ilz=1; ilz<l2p1; ilz++){
      int im1 = maski[ilz];
      int im2 = maski[ilz-1];
      bool ib1 = state & im1;
      bool ib2 = state & im2;
      if (ib1 && !ib2){// L^-_i gives nonzero
	double m = lzi[ilz];
	double ifct = sqrt(l*(l+1)-m*(m-1))*sign(state,min(im1,im2),max(im1,im2));
	int istate = state^im1^im2;
	for (int jlz=0; jlz<l2p1-1; jlz++){
	  int jm1 = maskj[jlz];
	  int jm2 = maskj[jlz+1];
	  if (cont && (ilz==jlz || im2==jm2)) continue;
	  bool jb1 = state & jm1;
	  bool jb2 = state & jm2;
	  if (jb1 && !jb2){//L^+_j gives nonzero
	    double m = lzj[jlz];
	    double jfct = sqrt(l*(l+1)-m*(m+1))*sign(istate,min(jm1,jm2),max(jm1,jm2));
	    int jstate = istate^jm1^jm2;
	    sts.push_back(make_pair(jstate, ifct*jfct));
	  }
	}
      }
    }
  }
  void L2(int state, list<pair<int,double> >& sts) const
  {
    sts.clear();
    // diagonal part
    double dd=0;
    for (int ilz=1; ilz<l2p1; ilz++){
      //l^+_i*l^-_i
      if ((mask_u[ilz]&state)&&(!(mask_u[ilz-1]&state)))
	dd += 0.5*(l*(l+1)-lz_u[ilz]*(lz_u[ilz]-1));
      if ((mask_d[ilz]&state)&&(!(mask_d[ilz-1]&state)))
	dd += 0.5*(l*(l+1)-lz_d[ilz]*(lz_d[ilz]-1));
    }
    for (int ilz=0; ilz<l2p1-1; ilz++){
      //l^-_i*l^+_i
      if ((mask_u[ilz]&state)&&(!(mask_u[ilz+1]&state)))
	dd += 0.5*(l*(l+1)-lz_u[ilz]*(lz_u[ilz]+1));
      if ((mask_d[ilz]&state)&&(!(mask_d[ilz+1]&state)))
	dd += 0.5*(l*(l+1)-lz_d[ilz]*(lz_d[ilz]+1));
    }
    // Lz^2
    double fct = sqr(Lz(state)) + dd;
    // store diagonal
    sts.push_back(make_pair(state,fct));
    
    // off diagonal
    L2_off_diag(state, sts, mask_u, lz_u, mask_u, lz_u, true);
    L2_off_diag(state, sts, mask_d, lz_d, mask_d, lz_d, true);
    L2_off_diag(state, sts, mask_u, lz_u, mask_d, lz_d, false);
    L2_off_diag(state, sts, mask_d, lz_d, mask_u, lz_u, false);
  }
  int N_el_before(int state, int i) const
  {
    int n=0;
    for (int q=0; q<i; q++) if (mask[q]&state) n++;
    return n;
  }
  void CoulombU(int state, const gaunt_ck& gck, const vector<double>& FkoJ,
		list<pair<int,vector<double> > >& sts) const
  {
    sts.clear();
    //    vector<double> Uk(l+1);
    vector<double> Uk(2); // very new change july 2006
    int ni=-1;
    for (int i=0; i<baths; i++){
      if (!(mask[i]&state)) continue;
      ni++;
      int state1 = state^mask[i];
      int m1 = lz[i];
      int s1 = sz[i];
      int nj=-1;
      for (int j=0; j<baths; j++){
	if (!(mask[j]&state1)) continue;
	nj++;
	// here we have: mask[i]&state && mask[j]&state
	int state2 = state1^mask[j];
	int m2 = lz[j];
	int s2 = sz[j];
	for (int a=0; a<baths; a++){
	  if (mask[a]&state2 || sz[a]!=s2) continue;
	  int na = N_el_before(state2,a);
	  int state3 = state2^mask[a];
	  int m3 = lz[a];
	  int s3 = sz[a];
	  for (int b=0; b<baths; b++){
	    if (mask[b]&state3 || sz[b]!=s1) continue;
	    int m4 = lz[b];
	    if (m4-m1!=m2-m3) continue;
	    int nb = N_el_before(state3,b);
	    int state4 = state3^mask[b];
	    int s4 = sz[b];
	    int sign = 1-2*((ni+nj+na+nb)&1);
	    Uk[0] = sign*gck(m4,m1,0)*gck(m2,m3,0);
	    double sum=0;
	    for (int k=1; k<=l; k++) sum += gck(m4,m1,2*k)*gck(m2,m3,2*k)*FkoJ[k];
	    Uk[1] = sign*sum;
	    sts.push_back(make_pair(state4, Uk));
	  }
	}
      }
    }
  }
  void Lp_to_N(int Ln, int state, map<int,double>& sts0, bool normalize=false){
    sts0.clear();
    sts0[state] = 1.0;
    for (int n=0; n<Ln; n++){
      map<int,double> sts;
      for (map<int,double>::const_iterator it=sts0.begin(); it!=sts0.end(); it++){
	int state = it->first;
	double fact = it->second;
	for (int ilz=0; ilz<l2p1-1; ilz++){
	  if (mask_u[ilz]&state && !(mask_u[ilz+1]&state))
	    sts[state^mask_u[ilz]^mask_u[ilz+1]] += fact*lp_fact[ilz];
	  if (mask_d[ilz]&state && !(mask_d[ilz+1]&state))
	    sts[state^mask_d[ilz]^mask_d[ilz+1]] += fact*lp_fact[ilz];
	}
      }
      sts0=sts;
    }
    if (normalize){
      double sum=0;
      for (map<int,double>::const_iterator l=sts0.begin(); l!=sts0.end(); l++) sum += sqr(l->second);
      double nrm = 1/sqrt(sum);
      for (map<int,double>::iterator l=sts0.begin(); l!=sts0.end(); l++) l->second *= nrm;
    }
  }
  void Ln_Sm(int Ln, int Sm, int state, map<int,double>& sts0) const
  {
    sts0.clear();
    sts0[state] = 1.0;
    // if Ln>0
    for (int n=0; n<Ln; n++){
      map<int,double> sts;
      for (map<int,double>::const_iterator it=sts0.begin(); it!=sts0.end(); it++){
	int state = it->first;
	double fact = it->second;
	for (int ilz=0; ilz<l2p1-1; ilz++){
	  if (mask_u[ilz]&state && !(mask_u[ilz+1]&state))
	    sts[state^mask_u[ilz]^mask_u[ilz+1]] += fact*lp_fact[ilz];
	  if (mask_d[ilz]&state && !(mask_d[ilz+1]&state))
	    sts[state^mask_d[ilz]^mask_d[ilz+1]] += fact*lp_fact[ilz];
	}
      }
      sts0=sts;
    }
    for (int n=0; n>Ln; n--){
      map<int,double> sts;
      for (map<int,double>::const_iterator it=sts0.begin(); it!=sts0.end(); it++){
	int state = it->first;
	double fact = it->second;
	for (int ilz=1; ilz<l2p1; ilz++){
	  if (mask_u[ilz]&state && !(mask_u[ilz-1]&state))
	    sts[state^mask_u[ilz]^mask_u[ilz-1]] += fact*lm_fact[ilz];
	  if (mask_d[ilz]&state && !(mask_d[ilz-1]&state))
	    sts[state^mask_d[ilz]^mask_d[ilz-1]] += fact*lm_fact[ilz];
	}
      }
      sts0=sts;
    }
    for (int n=0; n<Sm; n++){
      map<int,double> sts;
      for (map<int,double>::const_iterator it=sts0.begin(); it!=sts0.end(); it++){
	int state = it->first;
	double fact = it->second;
	for (int ilz=0; ilz<l2p1; ilz++)
	  if (mask_d[ilz]&state && !(mask_u[ilz]&state))
	    sts[state^mask_d[ilz]^mask_u[ilz]] += fact*sign(state,mask_d[ilz],mask_u[ilz]);// pazi, privzel da so down najprej!
      }
      sts0=sts;
    }
    for (int n=0; n>Sm; n--){
      map<int,double> sts;
      for (map<int,double>::const_iterator it=sts0.begin(); it!=sts0.end(); it++){
	int state = it->first;
	double fact = it->second;
	for (int ilz=0; ilz<l2p1; ilz++)
	  if (mask_u[ilz]&state && !(mask_d[ilz]&state))
	    sts[state^mask_u[ilz]^mask_d[ilz]] += fact*sign(state,mask_d[ilz],mask_u[ilz]);// pazi, privzel da so down najprej!
      }
      sts0=sts;
    }
  }
  void J2(int state, map<int,double>& sts0){
    map<int,double> sts_pm;
    Ln_Sm(1,-1,state,sts_pm);
    map<int,double> sts_mp;
    Ln_Sm(-1,1,state,sts_mp);
    list<pair<int,double> > sts_l2;
    L2(state,sts_l2);
    list<pair<int,double> > sts_s2;
    S2(state,sts_s2);
    sts0.clear();
    sts0[state] = 2*Lz(state)*Sz(state);
    for (map<int,double>::const_iterator i=sts_pm.begin(); i!=sts_pm.end(); i++) sts0[i->first] += i->second;
    for (map<int,double>::const_iterator i=sts_mp.begin(); i!=sts_mp.end(); i++) sts0[i->first] += i->second;
    for (list<pair<int,double> >::const_iterator i=sts_l2.begin(); i!=sts_l2.end(); i++) sts0[i->first] += i->second;
    for (list<pair<int,double> >::const_iterator i=sts_s2.begin(); i!=sts_s2.end(); i++) sts0[i->first] += i->second;
  }
  void li_dot_si(int state, map<int,double>& sts0) const
  {
    sts0.clear();
    double diag=0;
    for (int i=0; i<baths; i++) if (mask[i]&state) diag += lz[i]*0.5*sz[i];
    sts0[state] = diag;
    // ****** l_- s_+
    for (int ilz=1; ilz<l2p1; ilz++){
      //      cout<<ilz<<" "<<!!(mask_d[ilz]&state)<<" "<<!(mask_u[ilz-1]&state)<<endl;
      if (mask_d[ilz]&state && !(mask_u[ilz-1]&state))
	sts0[state^mask_d[ilz]^mask_u[ilz-1]] += 0.5*lm_fact[ilz]*sign(state,mask_d[ilz],mask_u[ilz-1]);
    }
    // ****** l_+ s_-
    for (int ilz=0; ilz<l2p1-1; ilz++){
      //      cout<<ilz<<" "<<!!(mask_u[ilz]&state)<<" "<<!(mask_d[ilz+1]&state)<<endl;
      if (mask_u[ilz]&state && !(mask_d[ilz+1]&state))
	sts0[state^mask_u[ilz]^mask_d[ilz+1]] += 0.5*lp_fact[ilz]*sign(state,mask_d[ilz+1],mask_u[ilz]);
    }
  }
  void Fp(int state, list<pair<int,double> >& sts0, double j, double jz, int l) const
  {

    ClebschGordan cg;
    sts0.clear();
    int lz;
    double sz;
    sz = 0.5;
    lz = static_cast<int>(round(jz-0.5));
    if (abs(lz)<=l){
      if (!(mask_u[l+lz]&state))
	sts0.push_back(make_pair(state^mask_u[l+lz], sign(state,mask_u[l+lz])*cg.CG(j,jz,l,lz,0.5,sz)));
    }
    sz = -0.5;
    lz = static_cast<int>(round(jz+0.5));
    if (abs(lz)<=l){
      if (!(mask_d[l+lz]&state))
	sts0.push_back(make_pair(state^mask_d[l+lz], sign(state,mask_d[l+lz])*cg.CG(j,jz,l,lz,0.5,sz)));
    }
  }
  string print(int state) const
  {
    stringstream stream;
    for (int lz=-l; lz<=l; lz++)
      stream<<setw(3)<<estate(state,lz)<<" ";
    return stream.str();
  }
  string printn(int state) const
  {
    stringstream stream;
    for (int lz=-l; lz<=l; lz++){
      int e = estate(state,lz);
      char ch;
      switch(e){
      case  0: ch='0'; break;
      case  2: ch='2'; break;
      case  1: ch='u'; break;
      case -1: ch='d'; break;
      }
      stream<<ch;
    }
    return stream.str();
  }
};



bool Project(function1D<double>& psi0, const function2D<double>& L2, int l0, int lmax)
{
  bool found = true;
  function1D<double> psi1(psi0.size());
  for (int l1=0; l1<lmax; l1++){
    if (l1==l0) continue;
    psi1.Product(L2,psi0);
    for (int i=0; i<psi1.size(); i++) psi1[i] = (psi1[i]-l1*(l1+1)*psi0[i])/(l0*(l0+1)-l1*(l1+1));
    psi0 = psi1;
  }
  double sum=0;
  for (int i=0; i<psi0.size(); i++) sum += sqr(psi0[i]);
  double renorm = 1/sqrt(sum);
  if (renorm>1e6) { renorm=0; found = false;}
  for (int i=0; i<psi0.size(); i++) psi0[i] *= renorm;
  //  cout<<renorm<<endl;
  // checking
  psi1.Product(L2,psi0);
  sum=0;
  for (int i=0; i<psi0.size(); i++) sum += fabs(psi1[i]-l0*(l0+1)*psi0[i]);
  if (found) found = sum<1e-6;
  return found;
}

void ConstructReducedBaseNLzSz(int n, int Lz, double Sz, const vector<int>& base, const operateLS& op,
			       vector<int>& reduced, vector<int>& red_base, vector<int>& red_ind)
{
  int kk=0;
  for (int i=0; i<base.size(); i++){
    int state = base[i];
    if (op.Nel(state)==n && op.Sz(state)==Sz && op.Lz(state)==Lz) {
      reduced.push_back(i);
      red_base.push_back(state);
      red_ind[state] = kk;
      kk++;
    }
  }
}
void ConstructReducedBaseNLzSz_1(int n, int Lz, double Sz, const vector<int>& base,
				 const vector<int>& bN, const vector<int>& bLz, const vector<int>& bSz,
				 vector<int>& reduced, vector<int>& red_ind)
{
  list<int> treduced;
  int kk=0;
  int iSz = static_cast<int>(round(Sz-0.25));
  for (int i=0; i<base.size(); i++){
    int state = base[i];
    if (bN[i]==n && bLz[i]==Lz && bSz[i]==iSz){
      treduced.push_back(i);
      red_ind[state] = kk;
      kk++;
    }
  }
  reduced.resize(treduced.size());
  int il=0;
  for (list<int>::iterator l=treduced.begin(); l!=treduced.end(); l++) reduced[il++] = *l;
}
void ConstructReducedBaseNJz(int n, double Jz, const vector<int>& base, const operateLS& op,
			     vector<int>& reduced, vector<int>& red_base, vector<int>& red_ind)
{
  int kk=0;
  for (int i=0; i<base.size(); i++){
    int state = base[i];
    if (op.Nel(state)==n && fabs(op.Jz(state)-Jz)<1e-10) {
      reduced.push_back(i);
      red_base.push_back(state);
      red_ind[state] = kk;
      kk++;
    }
  }
}
void ConstructReducedBaseNJz_1(int n, double Jz, const vector<int>& base, const operateLS& op,
			       vector<int>& reduced, vector<int>& red_ind)
{
  int kk=0;
  for (int i=0; i<base.size(); i++){
    int state = base[i];
    if (op.Nel(state)==n && fabs(op.Jz(state)-Jz)<1e-10) {
      reduced.push_back(i);
      red_ind[state] = kk;
      kk++;
    }
  }
}


void Find_LS_Base(const vector<int>& red_base, const vector<int>& red_ind, const operateLS& op, vector<pair<double,double> >& LSv,
		  function2D<double>& LSeigenv, vector<vector<int> >& LSk)
{
  // First, construct matrix S2 and L2 in the reduced base where Lz,Sz and n are good quantum numbers
  function2D<double> S2(red_base.size(),red_base.size());
  S2=0;
  list<pair<int,double> > lsp;
  for (int i=0; i<red_base.size(); i++){
    op.S2(red_base[i],lsp);
    for (list<pair<int,double> >::const_iterator l=lsp.begin(); l!=lsp.end(); l++) S2(red_ind[l->first],i)=l->second;
  }
  function2D<double> L2(red_base.size(),red_base.size());
  L2=0;
  for (int i=0; i<red_base.size(); i++){
    op.L2(red_base[i],lsp);
    for (list<pair<int,double> >::const_iterator l=lsp.begin(); l!=lsp.end(); l++) L2(red_ind[l->first],i)=l->second;
  }
  // Save them for testing purposes
  function2D<double> S2_orig(red_base.size(),red_base.size()), L2_orig(red_base.size(),red_base.size());
  S2_orig = S2;
  L2_orig = L2;
  // Diagonalize S2 and find eigensystem
  function1D<double> Seig(S2.size_N());
  function1D<double> work(4*S2.size_N());
  xsyev(S2.size_N(), S2.MemPt(), S2.fullsize_Nd(), Seig.MemPt(), work.MemPt(), work.size());
  // saves start and end of blocks that correspond to certain S2 eigenvalue
  list<pair<int,int> > start_end;
  double sw = 0.5*round(sqrt(1+4*Seig[0])-1);
  int start=0;
  for(int i=0; i<Seig.size(); i++){
    double ss = 0.5*round(sqrt(1+4*Seig[i])-1);
    if (fabs(ss-sw)>1e-12){
      if (i>start) start_end.push_back(make_pair(start,i));
      start=i; sw=ss;
    }
  }
  if (Seig.size()>start) start_end.push_back(make_pair(start,Seig.size()));
  // Now, calculate L2 in the eigenbase of S2. Must be block-diagonal since L2 and S2 commute
  function2D<double> Ln(red_base.size(), red_base.size());
  Ln.Product(L2,S2);
  L2.MProduct(S2,Ln);
  // Goes over all blocks of S2 (over all S2 values) and diagonalizes L2
  for (list<pair<int,int> >::const_iterator l=start_end.begin(); l!=start_end.end(); l++){
    int size = l->second-l->first;
    int offset = l->first;
    // The quantum number S in this block
    double St = 0.5*round(sqrt(1+4*Seig[offset])-1);
    // Takes out one block of L2
    function2D<double> Lnew(size,size);
    for (int i=0; i<size; i++)
      for (int j=0; j<size; j++) Lnew(i,j) = L2(offset+i,offset+j);
    // Solves for eigensystem of a block of L2
    function1D<double> Leig(size);
    xsyev(size, Lnew.MemPt(), Lnew.fullsize_Nd(), Leig.MemPt(), work.MemPt(), work.size());
    // Constructs eigenvectors in terms of direct-base functions (functions in red_base).
    // Also saves values of S and L for each eigenvector
    int kappa=0, L_last=-1;// index that distinguishes states with the same Lz,Sz,L,S
    for (int r=0; r<size; r++){
      for (int k=0; k<S2.size_Nd(); k++){
	double sum=0;
	for (int j=0; j<size; j++) sum += Lnew(r,j)*S2(j+offset,k);
	LSeigenv(r+offset,k) = sum;
      }
      // L for this eigenvector
      double Lt = 0.5*round(sqrt(1+4*Leig[r])-1);
      LSv[r+offset] = make_pair(Lt,St);
      vector<int> lsk(3);
      lsk[0] = static_cast<int>(Lt);
      lsk[1] = static_cast<int>(round(St-0.25));
      if (static_cast<int>(round(Lt))==L_last) kappa++;
      else{kappa=0; L_last=static_cast<int>(round(Lt));}
      lsk[2] = kappa;
      LSk[r+offset] = lsk;
    }
  }
  //Testing only below this point
  /*
  function2D<double> neki(red_base.size(),red_base.size());
  function2D<double> neki1(red_base.size(),red_base.size());
  neki.Product(S2_orig,LSeigenv);
  neki1.MProduct(LSeigenv,neki);
  for (int i=0; i<neki1.size_N(); i++){
    for (int j=0; j<neki1.size_Nd(); j++){
      if (fabs(neki1(i,j))<1e-8) neki1(i,j)=0;
      cout<<setw(4)<<neki1(i,j)<<" ";
    }
    cout<<endl;
  }
  */
}

void Calculate_Coulomb_Int(int l, const vector<int>& red_base, const vector<int>& red_ind, const operateLS& op,
			   const vector<pair<double,double> >& LSv, function2D<double>& LSeigenv, function1D<double>& LSU)
{
  // Ratio between F2,F4,F6 and J! At the end of the day, we want to use U and J only!
  vector<double> FkoJ = ComputeCoulombRatios(l);
  
//   FkoJ[0]=0;//This is not defined
//   if (l==1) FkoJ[1] = 5.;
//   if (l==2) { FkoJ[1]=14./(1+0.625); FkoJ[2]=0.625*FkoJ[1];}
//   if (l==3) { FkoJ[1]=6435./(286+195*0.668+250*0.494); FkoJ[2] = 0.668*FkoJ[1]; FkoJ[3]=0.494*FkoJ[1];}
  // Gaunt coefficients are calculated
  ClebschGordan cg;
  gaunt_ck gck(l, cg);
  // Matrix of Coulomb U (U[0]=U U[1]=J) in the reduced base (Lz,Sz)
  vector<function2D<double> > Us(2);
  for (int i=0; i<2; i++) { Us[i].resize(red_base.size(), red_base.size());   Us[i]=0; }
  
  list<pair<int,vector<double> > > sts;
  for (int i=0; i<red_base.size(); i++){
    op.CoulombU(red_base[i], gck, FkoJ, sts);
    for (list<pair<int,vector<double> > >::const_iterator r=sts.begin(); r!=sts.end(); r++){
      int newstate = r->first;
      int j = red_ind[newstate];
      for (int k=0; k<2; k++) Us[k](j,i) += 0.5*r->second[k]; // One-half prefactor because every term is caunted twice 
    }
  }
  // Calculates block diagonal form for U (in the base with good L,S,Lz,Sz)
  function2D<double> neki(red_base.size(),red_base.size()), Us_d(red_base.size(),red_base.size());
  neki.Product(Us[1],LSeigenv);
  Us_d.MProduct(LSeigenv,neki);
  // Find blocks of good L and S
  list<pair<int,int> > start_end;
  int start=0;
  double L0 = LSv[0].first;
  double S0 = LSv[0].second;
  for(int i=0; i<LSv.size(); i++){
    double L1 = LSv[i].first;
    double S1 = LSv[i].second;
    if (fabs(L1-L0)<1e-10 && fabs(S1-S0)<1e-10) continue;
    if (i>start) start_end.push_back(make_pair(start,i));
    start=i; L0 = L1; S0 = S1;
  }
  if (LSv.size()>start) start_end.push_back(make_pair(start,LSv.size()));
  // Goes over all blocks of good L and S and diagonalizes Coulomb U
  function1D<double> work(4*LSv.size());
  function2D<double> LSeigenv1(red_base.size(), red_base.size());
  LSeigenv1 = LSeigenv;
  for (list<pair<int,int> >::const_iterator l=start_end.begin(); l!=start_end.end(); l++){
    int size = l->second-l->first;
    int offset = l->first;
    if (size<=1) { LSU[offset] = Us_d(offset,offset); continue;}
    // Takes out the block
    function2D<double> Unew(size,size);
    for (int i=0; i<size; i++)
      for (int j=0; j<size; j++) Unew(i,j) = Us_d(offset+i,offset+j);
    // Solves for eigensystem of the block
    function1D<double> Ueig(size);
    xsyev(size, Unew.MemPt(), Unew.fullsize_Nd(), Ueig.MemPt(), work.MemPt(), work.size());
    // Constructs eigenvectors in terms of direct-base functions (functions in red_base).
    for (int r=0; r<size; r++){
      for (int k=0; k<LSeigenv.size_Nd(); k++){
	double sum=0;
	for (int j=0; j<size; j++) sum += Unew(r,j)*LSeigenv(j+offset,k);
	LSeigenv1(r+offset,k) = sum;
      }
      // U for this eigenvector
      LSU[r+offset] = Ueig[r];
    }
  }
  LSeigenv = LSeigenv1;
  
  // testing only
  /*
  function2D<double> temp(red_base.size(), red_base.size());
  temp.Product(Us[1],LSeigenv);
  Us_d.MProduct(LSeigenv,temp);
  for (int i=0; i<red_base.size(); i++){
    for (int j=0; j<red_base.size(); j++){
      if (fabs(Us_d(i,j))<1e-10) Us_d(i,j)=0;
      cout<<setw(11)<<Us_d(i,j)<<" ";
    }
    cout<<endl;
  }
  cout<<endl;
  */
}

void FindIndexToLargerBase(const vector<int>& base, const vector<int>& reducedn, const vector<int>& reducedm, vector<int>& ind_m_n)
{
  int k0=0;
  for (int i=0; i<reducedm.size(); i++){
    bool found=false;
    int k;
    for (k=k0; k<reducedn.size(); k++) if (reducedm[i] == reducedn[k]) {found=true;break;}
    if (!found) for (int k=0; k<k0; k++) if (reducedm[i] == reducedn[k]) {found=true;break;}
    if (!found){cerr<<"Can't find connection between base LzSz and corresponding Jz!"<<endl;}
    k0=k;
    ind_m_n[i]=k0;
  }
}

class LSJ_kappa{
  int Ln, Sn, Jn, kn;
  function2D<int> kappas;
  vector<vector<double> > LSJk;
  vector<int> index;
  int Lmax_, Smax_, Jmax_, kmax_;
public:
  void CreateIndex(int n, double Jz, const vector<vector<int> >& LSk, int& LSJ_size)
  {
    Lmax_=0; Smax_=0; Jmax_=0; kmax_=0;
    for (int i=0; i<LSk.size(); i++){
      if (LSk[i][0]>Lmax_) Lmax_ = LSk[i][0];
      if (LSk[i][1]>Smax_) Smax_ = LSk[i][1];
      if (LSk[i][0]+LSk[i][1]>Jmax_) Jmax_ = LSk[i][0]+LSk[i][1];
      if (LSk[i][2]+1>kmax_) kmax_ = LSk[i][2]+1;
    }
    Ln = Lmax_+1;
    Sn = Smax_+1;
    Jn = Jmax_+1;
    kn = kmax_;
    
    kappas.resize(Ln,Sn);
    kappas=0;
    for (int i=0; i<LSk.size(); i++) kappas(LSk[i][0],LSk[i][1]) = LSk[i][2]+1;

    index.resize(Ln*Sn*Jn*kn);
    for (int i=0; i<index.size(); i++){index[i]=-1;}
    vector<double> t(4);
    int nn=0;
    LSJk.clear();
    for (int S=0; S<=Smax_; S++){
      double St = S+(n%2)/2.;
      for (int L=0; L<=Lmax_; L++){
	if (kappas(L,S)==0) continue;
	nn+=kappas(L,S);
	double Jt = L+St;
	while(Jt>=(fabs(L-St)-1e-6)){
	  if (Jt<abs(Jz)) { Jt--; continue;}
	  for (int k=0; k<kappas(L,S); k++){
	    t[0] = L; t[1] = St; t[2] = Jt; t[3] = k;
	    LSJk.push_back(t);
	  }
	  Jt--;
	}
      }
    }
    if (nn!=LSk.size()){cerr<<"Sized are not correct!!"<<endl;}
    //    cout<<" Total="<<nn<<" "<<LSJk.size()<<endl;
    LSJ_size = LSJk.size();
    // Creating index array
    for (int i=0; i<LSJk.size(); i++) ind(static_cast<int>(LSJk[i][0]),LSJk[i][1],LSJk[i][2],static_cast<int>(LSJk[i][3])) = i;
  }
  int& ind(int L, double S, double J, int kappa){
    int iS = static_cast<int>(round(S-0.25));
    int iJ = static_cast<int>(round(J-0.25));
    int ii = kappa+kn*(iJ+Jn*(iS+Sn*L));
    if (L>=Ln || iS>=Sn || iJ>=Jn || kappa>=kn || L<0 || iS<0 || iJ<0 || kappa<0) {cerr<<"Index in LSJ_kappa out of range "<<iJ<<"-"<<Jn<<" "<<iS<<"-"<<Sn<<" "<<kappa<<"-"<<kn<<endl; return index[0];}
    return index[ii];
  }
  int ind(int L, double S, double J, int kappa) const{
    int iS = static_cast<int>(round(S-0.25));
    int iJ = static_cast<int>(round(J-0.25));
    int ii = kappa+kn*(iJ+Jn*(iS+Sn*L));
    if (L>=Ln || iS>=Sn || iJ>=Jn || kappa>=kn || L<0 || iS<0 || iJ<0 || kappa<0 || index[ii]==-1) {cerr<<"Index in LSJ_kappa out of range "<<iJ<<"."<<Jn<<endl; return 0;}
    return index[kappa+kn*(iJ+Jn*(iS+Sn*L))];
  }
  int nkappa(int L, double S){
    int iS = static_cast<int>(round(S-0.25));
    return kappas(L,iS);
  }
  const vector<double>& LSJkp(int i) const{
    if (i>=LSJk.size()){cerr<<"Out of range in LSJkp!"<<endl; exit(1); return LSJk[0];}
    return LSJk[i];
  }
  int Lmax(){return Lmax_;}
  int Smax(){return Smax_;}
  int Jmax(){return Jmax_;}
  int kmax(){return kmax_;}
};

class Ecmp{
  const vector<double>& E;
public:
  Ecmp(const vector<double>& E_) : E(E_){};
  bool operator()(int a, int b){
    if(a<0||a>=E.size() || b<0||b>=E.size()){
      cerr<<"Accesing elements out of range in sorting!!"<<endl; return false;
    }
    return E[a]<E[b];
  }
};
template<class T>
class Vcmp{
  const T& V;
public:
  Vcmp(const T& V_) : V(V_){};
  bool operator()(int a, int b){
    if(a<0||a>=V.size() || b<0||b>=V.size()){
      cerr<<"Accesing elements out of range in sorting!!"<<endl; return false;
    }
    return abs(V[a])>abs(V[b]);
  }
};


void PrintLSstates(const vector<pair<double,double> >& LSv, const vector<vector<int> >& LSk,  const function1D<double>& ULS)
{
  for (int i=0; i<LSv.size(); i++){
    cout<<" S="<<setw(6)<<left<<LSv[i].second<<right<<" L="<<setw(6)<<left<<LSv[i].first<<right<<" J="<<setw(10)<<left<<ULS[i]<<right<<" ";
    cout<<"    "<<" ~S="<<setw(2)<<left<<LSk[i][1]<<" kappa="<<setw(2)<<LSk[i][2];
    cout<<endl;
  }
}

inline int ind_Jz_i(int iJz, int i)
{
  static int iJzmax=100;
  return iJz + iJzmax*i;
}
inline int ind_Jz_i(pair<int,int>& piJz)
{ return ind_Jz_i(piJz.first,piJz.second);}
pair<int,int> inline Jz_i(int ind)
{
  static int iJzmax=100;
  return make_pair<int,int>(ind%iJzmax,ind/iJzmax);
}

class Eigenbase{
public:
  int n, l;
  double Jz;
  int iJz;
private:
  double Smax, Jmax;
  int Lmax, nJmax;
  vector<int> reduced, red_base;
  vector<int> red_ind;
  vector<pair<double,double> > LSv;
  vector<vector<int> > LSk;
  function2D<double> T_NLS;         // Transformation from direct base to (N,L,S,Lz,Sz)
  function1D<double> ULS;
  vector<int> red_indn, reducedn;                                                    // This is really needed
  LSJ_kappa LSJk;
  function2D<double> SpinOrbit;    // Spin-orbit interaction in (N,L,S,J,Jz)
  function1D<double> ULSJk;        // Coulomb interaction in (N,L,S,J,Jz)
  function2D<double> T_NLzSz_NLSJ; // Transformation (N,L,S,Lz,Sz) -> (N,L,S,J,Jz)
  function2D<double> T_NLSJ_Eig;   // Transformation (N,L,S,J,Jz) -> atom eigenbase
  function2D<double> T_NLzSz_Eig;   // Transformation (N,L,S,Lz,Sz) -> atom eigenbase // This is really needed
  friend void Compute_Fp(double j, double jz, const Eigenbase& E1, const Eigenbase& E2, const vector<int>& base, const operateLS& op, function2D<double>& mFp);
  //friend void Compute_Gz(const Eigenbase& Ei, const vector<int>& base, const operateLS& op, function2D<double>& Gz);
  friend void Compute_Moment(const Eigenbase& E1, const vector<int>& base, const operateLS& op, function2D<double>& Mu);
public:
  vector<int> Jkappa;
  vector<double> Jt_ein, eEner;    // J and E in eigenbase
  vector<vector<int> > JJzLSk;
  //  map<int,pair<int,double> > index0;
  //  vector<int> index1;
  vector<int> index2;
  vector<double> E_LS;
public:
  Eigenbase(int n_, int l_) : n(n_), l(l_) {};
  Eigenbase() : n(-1), l(-1) {};
  void Set_nl(int n_, int l_){n=n_;l=l_;}
  int size() const{return reducedn.size();}
  void CreateLSbase(const vector<int>& base, const operateLS& op)
  {
    int Lz = 0;         // the smallest Lz
    double Sz=0.5*(n&1);// the smallest allowed Sz
  
    red_ind.resize(base.size());
    ConstructReducedBaseNLzSz(n, Lz, Sz, base, op, reduced, red_base, red_ind);
    
    LSv.resize(red_base.size());
    LSk.resize(red_base.size());
  
    T_NLS.resize(red_base.size(), red_base.size());
    Find_LS_Base(red_base, red_ind, op, LSv, T_NLS, LSk);

    // Calculates Coulomb interaction
    ULS.resize(red_base.size());
    Calculate_Coulomb_Int(l, red_base, red_ind, op, LSv, T_NLS, ULS);
  }
  void PrintLS_states()
  {PrintLSstates(LSv, LSk, ULS);}
  void CreateNLSJbase_only(double Jz_, const vector<int>& base, const operateLS& op, const vector<int>& bN, const vector<int>& bLz, const vector<int>& bSz)
  {
    Jz = Jz_;
    reducedn.clear();
    int LSJ_size;
    LSJk.CreateIndex(n, Jz, LSk, LSJ_size); // Creates order for new LSJ base
    if (LSJ_size==0) return;
    red_indn.resize(base.size());           // index array for the base (N,Jz)
    
    Lmax = LSJk.Lmax();
    Smax = LSJk.Smax() + 0.5*(n%2);
    Jmax = LSJk.Jmax() + 0.5*(n%2);
    nJmax = LSJk.Jmax()+1;
    // And creates direct nJz base to be used for calculation of matrix elements
    ConstructReducedBaseNJz_1(n, Jz, base, op, reducedn, red_indn);
    if (reducedn.size()!=LSJ_size){cerr<<"Something wrong in constructing base LSJ_Jz"<<endl;}
  }
  void CreateNLSJbase(double Jz_, int iJz_, const vector<int>& base, const operateLS& op, const vector<int>& bN, const vector<int>& bLz, const vector<int>& bSz)
  {
    Jz = Jz_;
    iJz = iJz_;
    reducedn.clear();
    int LSJ_size;
    LSJk.CreateIndex(n, Jz, LSk, LSJ_size); // Creates order for new LSJ base
    if (LSJ_size==0) return;
    red_indn.resize(base.size());           // index array for the base (N,Jz)
    
    Lmax = LSJk.Lmax();
    Smax = LSJk.Smax() + 0.5*(n%2);
    Jmax = LSJk.Jmax() + 0.5*(n%2);
    nJmax = LSJk.Jmax()+1;
    // And creates direct nJz base to be used for calculation of matrix elements
    ConstructReducedBaseNJz_1(n, Jz, base, op, reducedn, red_indn);
    if (reducedn.size()!=LSJ_size){cerr<<"Something wrong in constructing base LSJ_Jz"<<endl;}
    
    // Finally, lets go to the LSJ base with good qantum numbers (n,L,S,J,Jz)
    ClebschGordan cg;
    vector<int> red_ind(base.size());
    T_NLzSz_NLSJ.resize(reducedn.size(), reducedn.size());
    T_NLzSz_NLSJ=0;
    int nall0=0;
    for (int Lz = -Lmax; Lz<=Lmax; Lz++){
      double Sz = Jz-Lz;
      if (fabs(Sz)>Smax) continue;
      //Now, constructing base for this particular Lz,Sz block
      vector<int> reducedm;
      for (int i=0; i<red_ind.size(); i++) red_ind[i] = -1;
      ConstructReducedBaseNLzSz_1(n, Lz, Sz, base, bN, bLz, bSz, reducedm, red_ind);
      if (reducedm.size()==0) continue; // This block is empty
      // And index to the larger n,Jz - base.
      vector<int> ind_m_n(reducedm.size());
      FindIndexToLargerBase(base, reducedn, reducedm, ind_m_n);
      // How many S^+ and L^+ needs to be performed to get to this particular Lz,Sz?
      int Lz_steps = Lz;
      int Sz_steps = static_cast<int>(round(Sz-(n%2)/2.));
      // Start performing operation (L^+)^{Lz_steps}(S^+)^{Sz_steps}
      function2D<double> Lpn(reduced.size(),reducedm.size());
      Lpn=0;
      for (int i=0; i<reduced.size(); i++){
	map<int,double> sts;
	op.Ln_Sm(Lz_steps, Sz_steps, base[reduced[i]], sts);
	for (map<int,double>::const_iterator l=sts.begin(); l!=sts.end(); l++){
	  int j = red_ind[l->first];
	  if (j<0) {cerr<<"Out of reduced base!"<<endl; exit(1);}
	  Lpn(i,j) = l->second;
	}
      }
      // These are eigenvectors with good L and S in the particular Lz,Sz sector
      function2D<double> T_NLS1(reduced.size(), reducedm.size());
      T_NLS1.MProduct(T_NLS,Lpn);
      // Normalizing new eigenvectors and store
      for (int i=0; i<reduced.size(); i++){
	int L = LSk[i][0];
	double S = LSk[i][1]+0.5*(n%2);
	if (L<abs(Lz)) continue;
	if (S<fabs(Sz)) continue;
	int kappa = LSk[i][2];
	double sum=0;
	for (int j=0; j<T_NLS1.size_Nd(); j++) sum += sqr(T_NLS1(i,j));
	double nrm = 1/sqrt(sum);
	if(nrm>1e10) {nrm=0; cerr<<" Did not expect null eigenvector in this sector!"<<endl;}
	for (int j=0; j<T_NLS1.size_Nd(); j++) T_NLS1(i,j) *= nrm;
	// Storing eigenvectors
	double Jt=fabs(L-S);
	while(Jt<=L+S+1e-6){
	  if (Jt<fabs(Jz)){Jt++; continue;}
	  int ii = LSJk.ind(L, S, Jt, kappa);
	  for (int j=0; j<T_NLS1.size_Nd(); j++)
	    T_NLzSz_NLSJ(ii,ind_m_n[j]) += T_NLS1(i,j)*cg.CG(Jt, Jz, L , Lz, S, Sz);
	  Jt++;
	}
	nall0++;
      }
    }
    if (nall0!=reducedn.size()){cerr<<"Something wrong in constructing J2 states!"<<endl;}
    //    cout<<"nall="<<nall0<<endl;
  }
  void Create_Coulomb_and_SO_NLSJ(const vector<int>& base, const operateLS& op)
  {
    // Spin orbit interaction can be written only in this enlarged (reducedn) base, where all possible L and S states are included
    SpinOrbit.resize(reducedn.size(),reducedn.size());
    SpinOrbit=0;
    for (int i=0; i<reducedn.size(); i++){
      map<int,double> sts;
      op.li_dot_si(base[reducedn[i]], sts);
      for (map<int,double>::const_iterator l=sts.begin(); l!=sts.end(); l++){
	int j = red_indn[l->first];
	if (j<0) {cerr<<"Out of reduced base!"<<endl; exit(1);}
	SpinOrbit(i,j) += l->second;
      }
    }
    // We need spin-orbit interaction in the LSJ base too. It has some offdiagonal elements
    function2D<double> temp(reducedn.size(), reducedn.size());
    temp.MProduct(T_NLzSz_NLSJ,SpinOrbit);
    SpinOrbit.Product(temp,T_NLzSz_NLSJ);

    // Find Coulomb interaction in new LSJk base
    ULSJk.resize(reducedn.size());
    for (int i=0; i<LSk.size(); i++){
      int L = LSk[i][0];
      double S = LSk[i][1]+0.5*(n%2);
      int kappa = LSk[i][2];
      double Jt=fabs(L-S);
      while(Jt<=L+S+1e-6){
	if (Jt<abs(Jz)) { Jt++; continue;}
	int ii = LSJk.ind(L, S, Jt, kappa);
	ULSJk[ii] = ULS[i];
	Jt++;
      }
    }
  }
  void ConstructEigenbase(double c_spinorb, double J_coulomb, const vector<int>& base, const operateLS& op)
  {// c_spinorb+U_coulomb gives H
    function2D<double> Energies(reducedn.size(), reducedn.size());
    for (int i=0; i<reducedn.size(); i++){
      for (int j=0; j<reducedn.size(); j++)
	Energies(i,j) = c_spinorb*SpinOrbit(i,j);
      Energies(i,i) += J_coulomb*ULSJk[i];
    }
    
    // Transformation from NSLJk base to atom eigenbase 
    T_NLSJ_Eig.resize(reducedn.size(),reducedn.size());
    T_NLSJ_Eig=0;
    // Since the matrix is block-diagonal in total spin J, we diagonalize each block separately
    Jt_ein.resize(reducedn.size());
    eEner.resize(reducedn.size());
    int offset=0;
    double Jt = (n%2)/2.;
    do{
      vector<int> tind;
      for (int i=0; i<reducedn.size(); i++)
	if (fabs(LSJk.LSJkp(i)[2]-Jt)<1e-6) tind.push_back(i);

      if (tind.size()==0){Jt++; continue;}
	
      function2D<double> Ener(tind.size(),tind.size());
      for (int i=0; i<tind.size(); i++)
	for (int j=0; j<tind.size(); j++)
	  Ener(i,j) = Energies(tind[i],tind[j]);
    
      function1D<double> eig(tind.size());
      function1D<double> work(tind.size()*5);
      xsyev(Ener.size_N(), Ener.MemPt(), Ener.fullsize_Nd(), eig.MemPt(), work.MemPt(), work.size());

      for (int i=0; i<tind.size(); i++)
	for (int j=0; j<tind.size(); j++)
	  T_NLSJ_Eig(i+offset, tind[j]) = Ener(i,j);
    
      for (int i=0; i<tind.size(); i++) eEner[offset+i] = eig[i];
      
      for (int i=0; i<tind.size(); i++) Jt_ein[offset+i]=Jt;
      offset+=tind.size();
      Jt++;
    }while(Jt-1e-6<=Jmax);
 
    // Sorts energies
    Ecmp ecmp(eEner);
    vector<int> eind(reducedn.size());
    for (int i=0; i<reducedn.size(); i++) eind[i] = i;
    sort(eind.begin(),eind.end(),ecmp);
    // Sorts eigenvectors and the rest of the things...
    {
      function2D<double> t_T(reducedn.size(),reducedn.size());
      vector<double> t_Jt(reducedn.size());
      vector<double> t_ene(reducedn.size());
      for (int i=0; i<reducedn.size(); i++){
	for (int j=0; j<reducedn.size(); j++) t_T(i,j) = T_NLSJ_Eig(eind[i],j);
	t_Jt[i] = Jt_ein[eind[i]];
	t_ene[i] = eEner[eind[i]];
      }
      T_NLSJ_Eig = t_T;
      Jt_ein = t_Jt;
      eEner = t_ene;
    }
    // Creates direct Transformation from (Lz,Sz) to eigenbase
    T_NLzSz_Eig.resize(reducedn.size(),reducedn.size());
    T_NLzSz_Eig.MProduct(T_NLSJ_Eig,T_NLzSz_NLSJ);
    // Each state can be distinguished by the dominant JJzLS quantum number. Create it.
    JJzLSk.resize(reducedn.size());
    for (int i=0; i<reducedn.size(); i++){
      vector<int> st(5);
      st[0] = static_cast<int>(round(Jt_ein[i]-0.25)); // J
      st[1] = static_cast<int>(round(Jt_ein[i]+Jz));   // Jz
      // find dominant state
      vector<int> vind(reducedn.size());
      for (int j=0; j<reducedn.size(); j++) vind[j] = j;
      Vcmp<functionb<double> > vcmp(T_NLSJ_Eig[i]);
      sort(vind.begin(),vind.end(),vcmp);
      int i0 = vind[0];
      st[2] = static_cast<int>(round(LSJk.LSJkp(i0)[0]));      // L
      st[3] = static_cast<int>(round(LSJk.LSJkp(i0)[1]-0.25)); // S
      st[4] = static_cast<int>(round(LSJk.LSJkp(i0)[3]));      // kappa
      if (LSJk.LSJkp(i0)[2]!=Jt_ein[i]) cout<<"ERROR: Strange! Seems the dominant character is not equal to the good quantum number J"<<endl;
      JJzLSk[i] = st;
    }
    // Creates index array Jkappa which gives any eigenstate unique number within states of the same J
    // in another words: states with the same J can be distinguished by themself with unique number kappa -> (J,kappa) is unique index
    Jkappa.resize(reducedn.size());
    //    index1.resize(reducedn.size());
    index2.resize(reducedn.size());
    map<int,int> countJ;
    for (int i=0; i<reducedn.size(); i++){
      int iJ = static_cast<int>(round(Jt_ein[i]-0.25));
      int local_iJz = static_cast<int>(round(Jt_ein[i]+Jz));
      Jkappa[i] = countJ[iJ]++;
      //      int indJkappaJz = indJkappa(iJ,Jkappa[i],local_iJz);
      //      index0.push_back(make_pair(indJkappa(iJ,Jkappa[i]),eEner[i]));
      //      index0[indJkappaJz] = make_pair(i,eEner[i]);
      //      index1[i] = indJkappaJz;
      index2[i] = ind_Jz_i(iJz,i);
    }


    /// POSKUS!!!
    // Transforming SpinOrbit Energies to Eigenbase. It can be usefull to subtract the average spin-orbit interaction
    // cause it is usually contained in impurity levels.
    function2D<double> temp(reducedn.size(), reducedn.size());
    temp.MProduct(T_NLSJ_Eig,SpinOrbit);
    SpinOrbit.Product(temp,T_NLSJ_Eig);
    E_LS.resize(reducedn.size());
    for (int i=0; i<reducedn.size(); i++) E_LS[i] = c_spinorb*SpinOrbit(i,i);
//     cout<<"Testing eigeneneries!!"<<endl;
//     for (int i=0; i<reducedn.size(); i++){
//       for (int j=0; j<reducedn.size(); j++){
// 	if (fabs(Energies(i,j))<1e-6) SpinOrbit(i,j)=0;
// 	cout<<setw(11)<<SpinOrbit(i,j)<<" ";
//       }
//       cout<<endl;
//     }
  }
  void PrintAtomEnergies()
  {
    // Prints energies and the dominant character of the eigenvectors
    for (int i=0; i<Jt_ein.size(); i++) {
      vector<int> vind(reducedn.size());
      for (int j=0; j<reducedn.size(); j++) vind[j] = j;
      Vcmp<functionb<double> > vcmp(T_NLSJ_Eig[i]);
      sort(vind.begin(),vind.end(),vcmp);
      int i0 = vind[0], i1;
      double P0 = sqr(T_NLSJ_Eig[i][i0]), P1;
      bool PrintP1=false;
      if (P0<0.999999){
	i1 = vind[1];
	P1 = sqr(T_NLSJ_Eig[i][i1]);
	PrintP1=true;
      }
      cout<<setw(3)<<i<<" J="<<setw(11)<<Jt_ein[i]<<" E="<<setw(11)<<eEner[i]<<right<<" Jkappa="<<setw(3)<<Jkappa[i];
      cout<<" P0="<<setw(11)<<P0<<" (LSJk)=("<<setw(2)<<LSJk.LSJkp(i0)[0]<<","<<setw(2)<<LSJk.LSJkp(i0)[1]<<","<<setw(4)<<LSJk.LSJkp(i0)[2]<<","<<setw(2)<<LSJk.LSJkp(i0)[3]<<")  |  ";
      if (PrintP1) cout<<" P1="<<setw(11)<<P1<<" (LSJk)=("<<setw(2)<<LSJk.LSJkp(i1)[0]<<","<<setw(2)<<LSJk.LSJkp(i1)[1]<<","<<setw(4)<<LSJk.LSJkp(i1)[2]<<","<<setw(2)<<LSJk.LSJkp(i1)[3]<<")";
      cout<<left<<endl;;    
    }
  }
  void FindEquivalent(map<double,list<vector<int> > >& states)
  {
    for (int i=0; i<reducedn.size(); i++)
      if (states.find(eEner[i])!=states.end())
	states[eEner[i]].push_back(JJzLSk[i]);
      else{
	bool found=false;
	for (map<double,list<vector<int> > >::iterator j=states.begin(); j!=states.end(); j++){
	  if (fabs(j->first-eEner[i])<1e-6){
	    j->second.push_back(JJzLSk[i]);
	    found=true;
	    break;
	  }
	}
	if (!found) states[eEner[i]].push_back(JJzLSk[i]);
      }
  }
};

//////////////////////////////// Small utility functions below this point ////////////////////////////////////////////////
void Compute_Fp(double j, double jz, const Eigenbase& E1, const Eigenbase& E2, const vector<int>& base, const operateLS& op, function2D<double>& mFp)
{
  if (E2.n!=E1.n+1 || E2.Jz!=E1.Jz+jz) cerr<<"Not correct eigenbasis in Compute_Fp!"<<endl;
  function2D<double> Fp(E2.reducedn.size(), E1.reducedn.size());
  Fp = 0;
  list<pair<int,double> > sts;
  for (int i=0; i<E1.reducedn.size(); i++){
    int state = base[E1.reducedn[i]];
    op.Fp(state, sts, j, jz, E1.l);
    for (list<pair<int,double> >::const_iterator jt=sts.begin(); jt!=sts.end(); jt++){
      int ii = E2.red_indn[jt->first];
      if (ii<0 || ii>=E2.reducedn.size()) {cerr<<"Out of range in Compute_Fp"<<endl; ii=0;}
      Fp(ii,i) += jt->second;
    }
  }

  function2D<double> temp1(E2.reducedn.size(), E1.reducedn.size());
  temp1.Product(Fp,E1.T_NLzSz_Eig);
  mFp.MProduct(E2.T_NLzSz_Eig,temp1);
}


//////////////////////////////// Small utility functions below this point ////////////////////////////////////////////////
void Compute_Moment(const Eigenbase& E1, const vector<int>& base, const operateLS& op, function2D<double>& Mu)
{
  function2D<double> Mt(E1.reducedn.size(), E1.reducedn.size());
  Mt = 0;
  for (int i=0; i<E1.reducedn.size(); i++){
    int state = base[E1.reducedn[i]];
    double m0 = op.Lz(state) + 2*op.Sz(state);
    Mt(i,i) = m0;
  }

  function2D<double> temp1(E1.reducedn.size(), E1.reducedn.size());
  temp1.Product(Mt,E1.T_NLzSz_Eig);
  Mu.MProduct(E1.T_NLzSz_Eig,temp1);
}



// void Compute_Gz(const Eigenbase& Ei, const vector<int>& base, const operateLS& op, function2D<double>& Gz)
// {
//   function2D<double> Gz0(Ei.reducedn.size(), Ei.reducedn.size());
//   Gz0 = 0;
//   for (int i=0; i<Ei.reducedn.size(); i++){
//     int state = base[Ei.reducedn[i]];
//     Gz0(i,i) = 0.5*(op.Lz(state)+2*op.Sz(state));
//   }

//   function2D<double> temp1(Ei.reducedn.size(), Ei.reducedn.size());
//   temp1.Product(Gz0,Ei.T_NLzSz_Eig);
//   Gz.MProduct(Ei.T_NLzSz_Eig,temp1);
// }

// void SimpleBubbles(const Eigenbase& E1, const Eigenbase& E2, const function2D<double>& Fp, map<int, map<int,double> >& bubbles)
// {
//   for (int i=0; i<E1.size(); i++){
//     int i0 = ind_Jz_i(E1.iJz,i);
//     for (int j=0; j<E2.size(); j++){
//       int j0 = ind_Jz_i(E2.iJz,j);
//       if (fabs(Fp(j,i))<1e-3) continue;
//       bubbles[i0][j0] += sqr(Fp(j,i));
//     }
//   }
// }

// class OCAd{
// public:
//   vector<int> states;
//   int ib1, ib2;
//   double f;
//   OCAd(int i0, int i1, int i2, int i3, int ib1_, int ib2_, double f_) : states(4), ib1(ib1_), ib2(ib2_), f(f_)
//   { states[0] = i0;    states[1] = i1;    states[2] = i2;    states[3] = i3;  }
//   bool operator==(const OCAd& oc)
//   {
//     bool eq  = (states[0]==oc.states[0] && states[1]==oc.states[1] && states[2]==oc.states[2] && states[3]==oc.states[3] && ib1==oc.ib1 && ib2==oc.ib2);
//     bool sim = (states[0]==oc.states[0] && states[1]==oc.states[3] && states[2]==oc.states[2] && states[3]==oc.states[1] && ib1==oc.ib2 && ib2==oc.ib1);
//     return eq || sim;
//   }
// };

// bool cmpOCA(const OCAd& a, const OCAd& b)
// {  return fabs(a.f)>fabs(b.f);}

// void OcaDiagrams(const Eigenbase& E0, const Eigenbase& E1, const Eigenbase& E2, const Eigenbase& E3,
// 		 const function2D<double>& Fp1, const function2D<double>& Fp2, const function2D<double>& Fp3, const function2D<double>& Fp4,
// 		 list<OCAd>& diagrams, int ib1, int ib2, int in, svector<double> E_ground_state, double E_oca)
// {
//   for (int i0=0; i0<E0.size(); i0++){
//     int il0 = ind_Jz_i(E0.iJz,i0);
//     if (fabs(E0.eEner[i0]-E_ground_state[in])>E_oca) continue;
//     for (int i1=0; i1<E1.size(); i1++){
//       int il1 = ind_Jz_i(E1.iJz,i1);
//       if (fabs(E1.eEner[i1]-E_ground_state[in+1])>E_oca) continue;
//       for (int i2=0; i2<E2.size(); i2++){
// 	int il2 = ind_Jz_i(E2.iJz,i2);
// 	if (fabs(E2.eEner[i2]-E_ground_state[in+2])>E_oca) continue;
// 	for (int i3=0; i3<E3.size(); i3++){
// 	  int il3 = ind_Jz_i(E3.iJz,i3);
// 	  if (fabs(E3.eEner[i3]-E_ground_state[in+1])>E_oca) continue;
// 	  double v1 = Fp1(i1,i0);
// 	  double v2 = Fp2(i2,i1);
// 	  double v3 = Fp3(i2,i3);
// 	  double v4 = Fp4(i3,i0);
// 	  double v = v1*v2*v3*v4;
// 	  if (fabs(v)<1e-3) continue;
// 	  OCAd oca(il0,il1,il2,il3,ib1,ib2,v);
// 	  diagrams.push_back(oca);
// 	}
//       }
//     }
//   }
// }

// void AddDiagram(const OCAd& ocad, list<OCAd>& list_oca)
// {
//   for (list<OCAd>::iterator p=list_oca.begin(); p!=list_oca.end(); p++)
//     if ((*p)==ocad) {p->f += ocad.f; return;}
//   list_oca.push_back(ocad);
// }
// void CompressOCA(const list<OCAd>& diagrams, map<int,int>& index0, map<int,int>& index1, map<int,int>& index2,
// 		 const vector<int>& cmp_one_electron_1, list<OCAd>& list_oca)
// {
//   for (list<OCAd>::const_iterator it0=diagrams.begin(); it0!=diagrams.end(); it0++){
//     int i0 = index0[it0->states[0]];
//     int i1 = index1[it0->states[1]];
//     int i2 = index2[it0->states[2]];
//     int i3 = index1[it0->states[3]];
//     int ib1 = cmp_one_electron_1[it0->ib1];
//     int ib2 = cmp_one_electron_1[it0->ib2];
//     OCAd ocad(i0,i1,i2,i3,ib1,ib2,it0->f);
//     AddDiagram(ocad, list_oca);
//   }
// }

// void ChiBubbles(const Eigenbase& Eb, const function2D<double>& Gz, map<int, map<int,double> >& Chi_bubbles)
// {
//   for (int i=0; i<Eb.size(); i++){
//     int i0 = ind_Jz_i(Eb.iJz,i);
//     for (int j=0; j<Eb.size(); j++){
//       int j0 = ind_Jz_i(Eb.iJz,j);
//       if (fabs(Gz(j,i))<1e-3) continue;
//       Chi_bubbles[i0][j0] += sqr(Gz(j,i));
//     }
//   }
// }

void Set_fmfp(const Eigenbase& E1, const Eigenbase& E2, const function2D<double>& Fp, map<int,double>& fpfm, map<int,double>& fmfp)
{
  for (int i=0; i<E1.size(); i++){
    int i0 = ind_Jz_i(E1.iJz,i);
    double sum=0;
    for (int j=0; j<E2.size(); j++) sum += sqr(Fp(j,i));
    fmfp[i0] += sum;
  }
  for (int i=0; i<E2.size(); i++){
    int i0 = ind_Jz_i(E2.iJz,i);
    double sum=0;
    for (int j=0; j<E1.size(); j++) sum += sqr(Fp(i,j));
    fpfm[i0] += sum; 
  }
}


// void CompressBubbles(map<int, map<int,double> >& bubbles, map<int,int>& indexn, map<int,int>& indexp, vector<vector<double> >& nca)
// {
//   for (map<int,map<int,double> >::iterator it0=bubbles.begin(); it0!=bubbles.end(); it0++){
//     int i0 = it0->first;
//     int i = indexn[i0];
//     if (i>=nca.size()) cerr<<"First index is out of range in CompessBubbles!"<<endl;
//     for (map<int,double>::iterator jt0=bubbles[i0].begin(); jt0!=bubbles[i0].end(); jt0++){
//       int j0 = jt0->first;
//       int j = indexp[j0];
//       if (j>=nca[i].size()) cerr<<"Second index is out of range in CompessBubbles!"<<endl;
//       nca[i][j] += jt0->second;
//     }
//   }
// }

// void CompressBubbles_Gz(map<int, map<int,double> >& bubbles, map<int,int>& indexn, vector<vector<double> >& nca)
// {
//   for (map<int,map<int,double> >::iterator it0=bubbles.begin(); it0!=bubbles.end(); it0++){
//     int i0 = it0->first;
//     int i = indexn[i0];
//     if (i>=nca.size()) cerr<<"First index is out of range in CompessBubbles_Gz!"<<endl;
//     for (map<int,double>::iterator jt0=bubbles[i0].begin(); jt0!=bubbles[i0].end(); jt0++){
//       int j0 = jt0->first;
//       int j = indexn[j0];
//       if (j>=nca[i].size()) cerr<<"Second index is out of range in CompessBubbles_Gz!"<<endl;
//       nca[i][j] += jt0->second;
//     }
//   }
// }

// bool CmpStates_(const pair<int,double>& a, const pair<int,double>& b)
// {return fabs(a.second)>fabs(b.second); }

// class CmpE{
//   map<int, vector<double> >& Es;
// public:
//   CmpE(map<int, vector<double> >& Es_) : Es(Es_) {};
//   bool operator()(int a, int b){
//     pair<int,int> pa = Jz_i(a), pb = Jz_i(b);
//     return Es[pa.first][pa.second]<Es[pb.first][pb.second];
//   }
// };

class TiJz{
public:
  int n;
  double Jmax;
  int iJmaxi;
  
  TiJz(int n_, int l) : n(n_)
  {
    //    Jmax = l*n+0.5*(n%2) + 2;
    Jmax = (l+0.5)*n;
    iJmaxi = static_cast<int>(2*(Jmax+2*l+1));
  }
  int give_iJz(double Jz)
  {return static_cast<int>(2*Jz) + iJmaxi;}
  double give_Jz(int iJz)
  { return (iJz-iJmaxi)/2.;}
  int iJz_min(int in)
  { return  (n+in+1)%2 + iJmaxi; }
  double J_max(){return Jmax;}
};

void RememberParams (int argc, char *argv[]){
  ofstream param ("history.atom", ios::app);
  if (!param) cerr<<" Didn't suceeded to open params file!"<<"history.atom"<<endl;
  for (int i=0; i<argc; i++) param << argv[i] << " ";
  param << endl;
}



int main(int argc, char *argv[], char *env[])
{
  //Sergej's program // Jur. Chem. Phys. 40, 3428, (1964)
  const double spinorb_default=-100;
  double J_coulomb = 0.0;//0.680290;//0.604
  double c_spinorb = spinorb_default;//0.367864976619714;//0.272116;//0.316
  int ns=0;
  int ne=2;
  int l=3;
  bool susceptibility = false;
  bool paramagnetic = false;
  //Energy window of excited states to keep!
  double Ekeep=100;
  double Ekeepc=100;
  int Nmax=1000;
  int Nmaxc=100;
  bool display_help = false;
  int i=0;
  string Impq;
  string Eimp;
  bool PrintMoment=false;
  bool ImpTot=false;
  bool QHB2=false;
  if (argc<=1) {display_help = true; i=-1;}
  while (++i<argc || display_help){
    std::string str(argv[i]);
    if (str=="-J" && i<argc-1) J_coulomb = atof(argv[++i]);
    if (str=="-cx" && i<argc-1) c_spinorb = atof(argv[++i]);
    if (str=="-ns" && i<argc-1) ns = atoi(argv[++i]);
    if (str=="-ne" && i<argc-1) ne = atoi(argv[++i]);
    if (str=="-l" && i<argc-1) l = atoi(argv[++i]);
    if (str=="-Ekeep" && i<argc-1) Ekeep = atof(argv[++i]);
    if (str=="-Ekeepc" && i<argc-1) Ekeepc = atof(argv[++i]);
    if (str=="-Nmax" && i<argc-1) Nmax = atoi(argv[++i]);
    if (str=="-Nmaxc" && i<argc-1) Nmaxc = atoi(argv[++i]);
    if (str=="-para") paramagnetic = true;
    if (str=="-Impq") Impq = argv[++i];
    if (str=="-Eimp") Eimp = argv[++i];
    if (str=="-pm") PrintMoment=true;
    if (str=="-ImpTot") ImpTot=true;
    if (str=="-HB2") QHB2=true;
    if (str=="-h" || str=="--help" || display_help){
      std::clog<<"**      Exact diagonalization of the atom           **\n";
      std::clog<<"**                                                  **\n";
      std::clog<<"**      Copyright Kristjan Haule, 26.09.2004        **\n";
      std::clog<<"******************************************************\n";
      std::clog<<"\n";
      std::clog<<"atom [Options]\n" ;
      std::clog<<"Options:   -J      Slater integrals F2=J*11.9219653179191 ("<<J_coulomb<<")\n";
      std::clog<<"           -cx     Spin-orbit coupling ("<<0.0<<")\n";
      std::clog<<"           -ns     Lowest occupancy of the atom ("<<ns<<")\n";
      std::clog<<"           -ne     Highest occupancy of the atom ("<<ne<<")\n";
      std::clog<<"           -l      Orbital angular momentum of the shel ("<<l<<")\n";
      std::clog<<"           -Ekeep  Energy window treated exactly ("<<Ekeep<<")\n";
      std::clog<<"           -Ekeepc Energy window treated exactly in core ("<<Ekeepc<<")\n";
      std::clog<<"           -Nmax   Maximum number of states kept in superstate ("<<Nmax<<")\n";
      std::clog<<"           -Nmaxc  Maximum number of states kept in superstate for core ("<<Nmaxc<<")\n";
      std::clog<<"           -para   Runs in paramagnetic mode ("<<paramagnetic<<")\n";
      std::clog<<"           -Impq   List of index for equivalent orbitals ("<<Impq<<")\n";
      std::clog<<"           -Eimp   List of impurity levels ("<<Eimp<<")\n";
      std::clog<<"           -pm     Prints moment instead of Sz ("<<PrintMoment<<")\n";
      std::clog<<"           -ImpTot Take full Eimp and not just splitting of Eimp - do not subtract Eimp[0] ("<<ImpTot<<")\n";
      std::clog<<"           -HB2    If we want to compute self-energy from Swinger-like equation and two particle response ("<<QHB2<<")\n";
      std::clog<<"*****************************************************\n";
      return 0;
    }
  }

  RememberParams (argc, argv);
  if (Nmaxc>Nmax) Nmaxc=Nmax;
  if (Ekeepc>Ekeep) Ekeepc = Ekeep;
  
  ClebschGordan cg;
  int dim = 2*(2*l+1);
  int global_l = l;
  // Creates basic order for encoding states
  vector<pair<int,double> > ml_ms;
  vector<pair<double,double> > j_mj;
  for (double ms = -0.5; ms<=0.5; ms++){
    for (int ml=-l; ml<=l; ml++){
      ml_ms.push_back(make_pair(ml,ms));
    }
  }
  for (double j=abs(l-0.5); j<=l+0.5; j++){
    for (double mj=-j; mj<=j; mj++){
      j_mj.push_back(make_pair(j,mj));
    }
  }
  if (dim != ml_ms.size() || dim != j_mj.size()){ cerr<<"Something wrong in dimensions lm,jmj"<<endl; exit(1);}

  
  deque<int> Impqs;
  if (Impq!=""){
    istringstream stream(Impq);
    stream.ignore(1);
    for (int i=0; i<200; i++){
      int a;
      stream>>a;
      if (!stream.good()) break;
      stream.ignore(1);
      Impqs.push_back(a);
    }
  }else{
    double j=0;
    int curind=-1;
    for (int i=0; i<j_mj.size(); i++)
      if (j_mj[i].first == j){
	Impqs.push_back(curind);
      }else{
	j=j_mj[i].first;
	curind += 1;
	Impqs.push_back(curind);
      }
  }
  deque<double> dE;
  if (Eimp!=""){
    istringstream stream(Eimp);
    stream.ignore(1);
    for (int i=0; i<200; i++){
      double a;
      stream>>a;
      if (!stream.good()) break;
      stream.ignore(1);
      dE.push_back(a);
    }
  }else{
    int nsize=0;
    for (int i=0; i<Impqs.size(); i++)
      if (Impqs[i]>nsize) nsize=Impqs[i];
    nsize++;
    dE.resize(nsize);
    for (int i=0; i<dE.size(); i++) dE[i]=0;
  }
  if (!ImpTot){
    double dE0 = dE[0];
    for (int i=0; i<dE.size(); i++) dE[i] -= dE0;
  }
  

  cout<<"Impqs="<<endl;
  for (int i=0; i<Impqs.size(); i++){
    cout<<i<<" "<<Impqs[i]<<endl;
  }
  cout<<"dE="<<endl;
  for (int i=0; i<dE.size(); i++){
    cout<<i<<" "<<dE[i]<<endl;
  }

  /*
    This code allows to input c_spinorb, which might be different than the splitting of the impurity levels.
    By default we take c_spinorb to be compatible with the splitting of the 5/2 and 7/2 impurity level
   */
  if (c_spinorb==spinorb_default){
    double E_l_m_1_2=0;
    for (int k=0; k<2*l; k++) E_l_m_1_2 += dE[Impqs[k]];
    E_l_m_1_2 *= 1./(2*l);
    double E_l_p_1_2=0;
    for (int k=2*l; k<2*(2*l+1); k++) E_l_p_1_2 += dE[Impqs[k]];
    E_l_p_1_2 *= 1./(2*(l+1.));
    c_spinorb = (E_l_p_1_2-E_l_m_1_2)/(l+0.5);
    cout<<"c_spinorb changed to "<<c_spinorb<<endl;
  }

  // ********** This is the true start!!!
  int Nbase = 1<<dim;
  vector<int> base(Nbase);
  // Creates direct base with lzi and szi as good quantum numbers
  // This class can perform certain simple operations on a direct state such as L^+, S^+, S^2, L^2
  operateLS op(base,ml_ms,l);

  vector<int>  bN(Nbase), bLz(Nbase), bSz(Nbase);
  for (int i=0; i<base.size(); i++){
    bN[i] = op.Nel(base[i]);
    bLz[i] = op.Lz(base[i]);
    bSz[i] = static_cast<int>(round(op.Sz(base[i])-0.25));
  }
  
  
  TiJz tiJz(ns+1+2*((ne-ns)/2+1),l); // How to convert Jz to array index

  int nstart = -1;
  int nstop = ne-ns+2;
  
  //  int nstart=0, nstop=3;
  //  if (add_core){ nstart=-1; nstop=4;}
  if (nstart+ns<0) nstart=-ns;
  if (ns+nstop-1>2*(2*l+1)) nstop=1-ns+2*(2*l+1);

  
  svector<int> Np(nstart, nstop);
  for (int i=nstart; i<nstop; i++) Np[i] = ns+i;

  deque<vector<double> > Es, E_LSs; // atom energy and Spin-Orbit energy
  
  ///////////////////// Creates One electron basis for the bath function /////////////////////////////////
  vector<pair<double,double> > one_electron_base(2*(2*l+1));
  vector<int> equiv(one_electron_base.size());
  vector<int> partner(one_electron_base.size());
  function1D<int> global_flip(one_electron_base.size());
  map<int,int> deg;
  {
    int kk=0;
    for (int ic=-1; ic<=2; ic+=2){
      double jc=l+0.5*ic;
      if (jc<0) continue;
      for (int j2z=-static_cast<int>(2*jc); j2z<=static_cast<int>(2*jc); j2z+=2){
	double jz = j2z/2.;
	one_electron_base[kk] = make_pair(jc,jz);
	equiv[kk] = Impqs[kk];
	deg[Impqs[kk]]++;
	kk++;
      }
    }
    global_flip=0;
    kk=0;
    int ii=0;
    for (int ic=-1; ic<=2; ic+=2){
      double jc=l+0.5*ic;
      if (jc<0) continue;
      int j2z=-static_cast<int>(2*jc);
      for (; j2z<=0; j2z+=2){
	global_flip[kk]=ii;
	kk++;
	ii++;
      }
      for (; j2z<=static_cast<int>(2*jc); j2z+=2){
	ii--;
	global_flip[kk]=ii;
	kk++;
      }
      for(int i=0; i<kk; i++) ii = max(ii,global_flip[i]);
      ii++;
    }
    cout<<"global flip:"<<endl;
    for(int i=0; i<global_flip.size(); i++){
      cout<<i<<" "<<"global_flip="<<global_flip[i]<<endl;
    }
  }

  /*
  vector<function2D<double> > Uc(one_electron_base.size());
  if (QHB2){
    ClebschGordan cg;
    gaunt_ck gck(l, cg);
    vector<double> FkoJ = ComputeCoulombRatios(l);
    function2D<function2D<double> > Uf;
    CoulombU0(Uf, l, ml_ms, one_electron_base, cg, gck, FkoJ);

    for (int j1=0; j1<Uc.size(); j1++) Uc[j1].resize(one_electron_base.size(),one_electron_base.size());
    
    for (int j1=0; j1<Uf.size_N(); j1++){
      for (int j2=0; j2<Uf.size_Nd(); j2++){
	for (int j3=0; j3<Uf.size_N(); j3++){
	  double u = J_coulomb*(Uf(j1,j2)(j3,j1)-Uf(j1,j2)(j1,j3));
	  if (fabs(u)<1e-10) u=0;
	  Uc[j1](j2,j3) = u;
	  //cout<<u<<" ";
	}
	//cout<<endl;
      }
    }
  }
  */
  ///////////////////// Data structure to store results /////////////////////////////////
  svector<vector<map<int, double> > > fpfm(nstart,nstop), fmfp(nstart,nstop);
  for (int in=fpfm.start(); in<fpfm.stop(); in++){
    fpfm[in].resize(one_electron_base.size());
    fmfp[in].resize(one_electron_base.size());
  }
  int instart=nstart, instop=nstop;

  ///////////////////// Goes over all atomic states and diagonalizes blocks of H /////////////////////////////////
  {
    // Basic states for all necessary occupations and spins are created
    svector<vector<Eigenbase> > Ebase(nstart,nstop);
    svector<vector<vector<double> > > Moment(nstart,nstop);
    svector<vector<int> > state(nstart,nstop);
    deque<pair<int,int> > istate;
    deque<int> isize;
    int iq=0;
    for (int in=nstart; in<nstop; in++){
      int tn = Np[in];
      //      double tJmax = tn*l + 0.5*(tn%2);
      double tJmax = tiJz.Jmax + 0.5*((in+1)%2);
      int itJmax = tiJz.give_iJz(tJmax);//+1;
      Ebase[in].resize(itJmax);
      Moment[in].resize(itJmax);
      state[in].resize(itJmax);
      double Jz = -tJmax;
      int tot_size=0;
      do{
	int iJz = tiJz.give_iJz(Jz);
	if (abs(static_cast<int>(2*Jz))%2 != tn%2) continue;// half-integer spin only for odd number of electrons and vice versa
	if (fabs(Jz)>tn*l+0.5) continue;        // total J can not be bigger than that
	
	cout<<"iJz = "<<iJz<<" "<<itJmax<<endl;
	Ebase[in][iJz].Set_nl(Np[in], l);
	Ebase[in][iJz].CreateLSbase(base, op);
      
	Ebase[in][iJz].CreateNLSJbase(Jz, iJz, base, op, bN, bLz, bSz);
	if (Ebase[in][iJz].size()==0) continue;
	Ebase[in][iJz].Create_Coulomb_and_SO_NLSJ(base, op);
	Ebase[in][iJz].ConstructEigenbase(c_spinorb, J_coulomb, base, op);
	cout<<"--------------n="<<Ebase[in][iJz].n<<" Jz="<<Ebase[in][iJz].Jz<<"---------------"<<endl;
	Ebase[in][iJz].PrintAtomEnergies();


	function2D<double> Mu(Ebase[in][iJz].size(), Ebase[in][iJz].size());
	Compute_Moment(Ebase[in][iJz], base, op, Mu);
	cout<<"--------------mu(0,0)="<<Mu(0,0)<<endl;
	Moment[in][iJz].resize(Ebase[in][iJz].size());
	for (int ii=0; ii<Moment[in][iJz].size(); ii++) Moment[in][iJz][ii] = Mu(ii,ii);
	
	
	// To save for later use
	//	indexes[in][iJz] = Ebase[in][iJz].index2;
	//	JJzLSks[in][iJz] = Ebase[in][iJz].JJzLSk;
	
	Es.push_back(Ebase[in][iJz].eEner);
	E_LSs.push_back(Ebase[in][iJz].E_LS);
	
	istate.push_back(make_pair(in,iJz));
	state[in][iJz] = iq;
	isize.push_back(Ebase[in][iJz].size());
	iq++;
	
	tot_size += Ebase[in][iJz].size();
      } while(++Jz<=tJmax);
      cout<<"Total size = "<<tot_size<<endl;
    }
    
//     for (int i=0; i<istate.size(); i++){
//       int N = Np[istate[i].first];
//       for (int l=0; l<isize[i]; l++){
// 	Es[i][l] += F0*0.5*N*(N-1.);
//       }
//     }
    
    ///////////////////// Computes all possible bubbles between atomic states, also F^+F and FF^+ are computed ///////////////////
    function2D<int> Fp_index(istate.size(),one_electron_base.size());
    Fp_index = -1;
    vector<vector<function2D<double> > > Fp_matrix(istate.size());
    for (int i=0; i<Fp_matrix.size(); i++) Fp_matrix[i].resize(one_electron_base.size());
    
    for (int in=nstart; in<nstop-1; in++){
      double Jmax = tiJz.Jmax + 0.5*((in+1)%2);// + (in>0 ? in+2 : 0);
      double Jz = -Jmax;
      do{
	int iJzn = tiJz.give_iJz(Jz);
	if (iJzn>=Ebase[in].size() || Ebase[in][iJzn].size()==0) continue;
	for (int ij=0; ij<one_electron_base.size(); ij++){
	  double j = one_electron_base[ij].first;
	  double jz = one_electron_base[ij].second;
	  cout<<"Constructing diagrams for in="<<in<<" Jz="<<Jz<<" j="<<j<<" jz="<<jz<<" ---------------"<<endl;
	  int iJzp = tiJz.give_iJz(Jz + jz);
	  if (iJzp>=Ebase[in+1].size() || Ebase[in+1][iJzp].size()==0) continue;

	  int ist = state[in][iJzn];
	  int jst = state[in+1][iJzp];
	  Fp_index[ist][ij] = jst;
	  Fp_matrix[ist][ij].resize(isize[jst],isize[ist]);
	  
	  function2D<double> Fp(Ebase[in+1][iJzp].size(), Ebase[in][iJzn].size());
	  Compute_Fp(j, jz, Ebase[in][iJzn], Ebase[in+1][iJzp], base, op, Fp); // This computes (F^{j,jz,+})_{n+1,n}
	  Fp_matrix[ist][ij] = Fp;
	  Set_fmfp(Ebase[in][iJzn], Ebase[in+1][iJzp], Fp, fpfm[in+1][ij], fmfp[in][ij]);// Stores F.F^+ and F^+.F
	}
      } while(++Jz<=Jmax);
    }
	 
    // high-frequency moments
    vector<vector<function2D<double> > > FpFm(2*deg.size());
    vector<vector<function2D<double> > > FmFp(2*deg.size());
    for (int ib=0; ib<FpFm.size(); ib++){
      FpFm[ib].resize(istate.size());
      FmFp[ib].resize(istate.size());
      for (int l=0; l<FpFm[ib].size(); l++){
	FpFm[ib][l].resize(isize[l],isize[l]);
	FmFp[ib][l].resize(isize[l],isize[l]);
	FpFm[ib][l]=0;
	FmFp[ib][l]=0;
      }
    }
    int max_size=0;
    for (int i=0; i<isize.size(); i++) if (isize[i]>max_size) max_size = isize[i];    
    function2D<double> res3(max_size,max_size);
    
    for (int ist=0; ist<istate.size(); ist++){
      cout<<"ist="<<ist<<endl;
      for (int ij=0; ij<one_electron_base.size(); ij++){
	double j = one_electron_base[ij].first;
	int jst = Fp_index[ist][ij];//ist:N, jst:N+1
	if (jst<0) continue;
	function2D<double>& fp = Fp_matrix[ist][ij];
	
 	res3.Product("T", "N", fp, fp);
	FmFp[equiv[ij]][ist] += res3;

	res3.Product("N", "T", fp, fp);
	FpFm[equiv[ij]][jst] += res3;

      }
    }

    
    ///////////////////////////////////////////////////////////////
    ///// Decides how to cut some excited states   ////////////////
    ///// The parameters used are: Nmax and Ekeep  ////////////////
    ///////////////////////////////////////////////////////////////
    // first remove core
    int tin=-1;
    int istart = -1;
    do{
      tin = istate[++istart].first;
    }while(tin<0);
    tin = ne-ns+1;
    int istop = istate.size();
    do{
      tin = istate[--istop].first;
    }while(tin>ne-ns);
    istop++;
    
    function1D<int> nisize(istop);//(isize.size());
    for (int i=istart; i<istop; i++) nisize[i] = isize[i];
    for (int in = nstart; in<nstop; in++){
      double Egs = 1e10;
      for (int i=istart; i<istop; i++)
	if (istate[i].first==in && Es[i][0]<Egs) Egs = Es[i][0];
      cout<<"Egs["<<in<<"]= "<<Egs<<endl;

      double iEkeep = (in<0||in>2) ? Ekeepc : Ekeep;
      for (int i=istart; i<istop; i++){
	if (istate[i].first!=in) continue;
	cout<<i<<"  ";
	for (int m=0; m<isize[i]; m++){
	  if (fabs(Es[i][m]-Egs)>iEkeep){
	    nisize[i]=m;
	    break;
	  }
	  cout<<setw(12)<<Es[i][m]-Egs<<" ";
	}
	cout<<"  "<<nisize[i]<<endl;
      }
    }
    for (int i=istart; i<istop; i++){
      int iNmax = Nmax;
      int in = istate[i].first;
      if (in<0 || in>2) iNmax = Nmaxc;
      nisize[i] = min(nisize[i],iNmax);
    }
    deque<int> nindex;
    {
      int l=0;
      for (int i=istart; i<istop; i++){
	if (nisize[i]==0) continue;
	nindex.push_back(i);
      }
    }
    function1D<int> inindex(istate.size());//(istop-istart);
    inindex=-1;
    for (int s=0; s<nindex.size(); s++) inindex[nindex[s]]=s;

    cout<<"inindex:"<<endl;
    for (int i=0; i<inindex.size(); i++) cout<<i<<" "<<inindex[i]<<endl;
					
    for (int i=0; i<nindex.size(); i++){
      cout<<setw(3)<<i<<" "<<setw(3)<<nindex[i]<<" "<<setw(5)<<nisize[nindex[i]]<<endl;
    }
    
    
    
    map<int,double> epsk;
    for (int b=0; b<one_electron_base.size(); b++){
      double j = one_electron_base[b].first;
      double ls = c_spinorb*0.5*(j*(j+1)-global_l*(global_l+1)-0.5*1.5);
      //cout<<"adding to "<<Impqs[b]<<" value "<<ls<<" "<<one_electron_base[b].first<<" "<<one_electron_base[b].second<<endl;
      epsk[Impqs[b]] += ls;
    }
    for (int i=0; i<epsk.size(); i++) epsk[i]/=deg[i];

    
    cout<<"epsk_SO=";
    for (int i=0; i<epsk.size(); i++) cout<<epsk[i]<<" ";
    cout<<endl;
    
    // epsk now contains only spin-orbit.
    // Spin-orbit is in both : should be subtracted from dE!
    for (int i=0; i<dE.size(); i++) dE[i] -= epsk[i];
    
    
    // if c_spinorb is very large, we just want to keep 5/2's and remove 7/2's. In this case we should scale back
    // the E_{5/2} energies.
    if (c_spinorb>=1000.)
      for (int i=0; i<dE.size(); i++) dE[i] += c_spinorb*2;
    
    
    // epsk should be corrected for crystal fields
    for (int i=0; i<epsk.size(); i++) epsk[i] += dE[i];
    
    cout<<"epsk=";
    for (int i=0; i<epsk.size(); i++) cout<<epsk[i]<<" ";
    cout<<endl<<"dE=";
    for (int i=0; i<epsk.size(); i++) cout<<dE[i]<<" ";
    cout<<endl;

    
    cout<<"Correcting energies"<<endl;
    for (int i=0; i<istate.size(); i++){
      int in = istate[i].first;
      int N = Np[istate[i].first];
      for (int l=0; l<isize[i]; l++){
 	cout<<setw(3)<<i+1<<" "<<setw(3)<<l<<") "<<Es[i][l]<<" ";
	double Ntot = 0;
	double dEtot=0;
	for (int ik=0; ik<deg.size(); ik++){
	  double pm = FpFm[ik][i](l,l);
	  if (in<=0) pm = -FmFp[ik][i](l,l) + deg[ik];
	  Ntot += pm;
	  
	  dEtot += pm*dE[ik];
	  if (fabs(pm)<1e-10) pm=0;
	  cout<<pm<<" ";	  
	}	
	cout<<" "<<Ntot<<" "<<dEtot<<endl;
	Es[i][l] += dEtot;
      }
    }

    

    max_size=0;
    for (int i=istart; i<istop; i++)
      if (nisize[i]>max_size) max_size = nisize[i];

//     map<int,double> epsk;
//     for (int b=0; b<one_electron_base.size(); b++){
//       double j = one_electron_base[b].first;
//       double ls = c_spinorb*0.5*(j*(j+1)-global_l*(global_l+1)-0.5*1.5);
      
//       cout<<"adding to "<<Impqs[b]<<" value "<<ls<<" "<<one_electron_base[b].first<<" "<<one_electron_base[b].second<<endl;
      
//       epsk[Impqs[b]] += ls;
//     }
//     for (int i=0; i<epsk.size(); i++) epsk[i]/=deg[i];
//     cout<<"EPSK="<<endl;
//     for (int i=0; i<epsk.size(); i++) {
//       epsk[i] += dE[i];
//       cout<<epsk[i]<<" "<<dE[i]<<endl;
//     }


    if (c_spinorb>=1000.) one_electron_base.resize(6);
    int unique=0;
    for (int ij=0; ij<one_electron_base.size(); ij++) if (Impqs[ij]>=unique) unique++;
    
    
    cout<<"Printing"<<endl;
    ofstream gout("actqmc.cix"); gout.precision(12);
    gout<<"# Cix file for cluster DMFT with CTQMC: J="<<J_coulomb<<" cx="<<c_spinorb<<" Ekeep="<<Ekeep<<" Ekeepc="<<Ekeepc<<" Nmax="<<Nmax<<" Nmaxc="<<Nmaxc<<endl;
    gout<<"# cluster_size, number of states, number of baths, maximum_matrix_size "<<endl;
    
    int nbaths=Impqs.size();
    gout<<1<<" "<<nindex.size()<<" "<<one_electron_base.size()<<" "<<max_size<<endl;
    
    gout<<"# baths, dimension, symmetry"<<endl;
    //for (int ij=0; ij<Impqs.size(); ij++){
    for (int ij=0; ij<one_electron_base.size(); ij++){
      gout<<left<<setw(3)<<ij<<"  1 "<<setw(4)<<Impqs[ij]<<" "<<setw(4)<<global_flip[ij]<<endl;
    }

    
    gout<<right;
    gout<<"# cluster energies for non-equivalent baths, eps[k]"<<endl;

    for (int ig=0; ig<epsk.size(); ig++) gout<<epsk[ig]<<" ";
    gout<<endl;
    
    gout<<setw(2)<<"#"<<" "<<setw(3)<<"N"<<" "<<setw(3)<<"K"<<" "<<setw(7)<<"Jz"<<" "<<setw(3)<<"size"<<endl;
    for (int ia=0; ia<nindex.size(); ia++){
      int i = nindex[ia];
      int in = istate[i].first;
      int iJz = istate[i].second;
      double Jz = Ebase[in][iJz].Jz;
      double n = Ebase[in][iJz].n;
      double moment = Moment[in][iJz][0];
      
      gout<<setw(2)<<inindex[i]+1<<" "<<setw(3)<<n<<" "<<setw(3)<<0<<" ";
      if (!PrintMoment)
	gout<<setw(4)<<Jz<<" ";
      else{
	if (fabs(moment)<1e-10) moment=0;
	int defaultp = gout.precision();
	gout.precision(4);
	gout<<setw(7)<<moment<<" ";
	gout.precision(defaultp);
      }
      gout<<setw(3)<<nisize[i]<<"    ";
      
      for (int ij=0; ij<one_electron_base.size(); ij++){
	int jst = Fp_index[i][ij];
	int fi;
	if (jst>=istop || jst<istart) fi=0;
	else fi = inindex[jst]+1;
	gout<<setw(3)<<fi<<" ";
      }

      for (int l=0; l<nisize[i]; l++){
	double Eatom_LS = Es[i][l];
	if (fabs(Eatom_LS)<1e-10) Eatom_LS=0;
	gout<<setw(18)<<Eatom_LS<<" ";
      }
      for (int l=0; l<nisize[i]; l++) gout<<setw(18)<<Ebase[in][iJz].Jt_ein[l]<<" ";
      gout<<endl;
    }
    gout<<"# matrix elements"<<endl;
    for (int ia=0; ia<nindex.size(); ia++){
      int i = nindex[ia];
      int in = istate[i].first;
      int iJz = istate[i].second;
      double Jz = Ebase[in][iJz].Jz;
      double n = Ebase[in][iJz].n;
      for (int ij=0; ij<one_electron_base.size(); ij++){
	gout<<setw(2)<<inindex[i]+1<<" ";
	int jst = Fp_index[i][ij];
	int fi;
	if (jst>=istop || jst<istart) fi=0;
	else fi = inindex[jst]+1;
	gout<<setw(3)<<fi<<"  ";
	if (fi==0 || nisize[i]==0 || nisize[jst]==0) { gout<<setw(3)<<0<<" "<<setw(3)<<0<<"  "<<endl;; continue; }
	function2D<double>& fp = Fp_matrix[i][ij];
	gout<<setw(3)<<nisize[i]<<" "<<setw(3)<<nisize[jst]<<"  ";
	for (int it=0; it<nisize[i]; it++)
	  for (int jt=0; jt<nisize[jst]; jt++){
	    double ff = fp(jt,it);
	    if (fabs(ff)<1e-10) ff=0;
	    gout<<setw(20)<<ff<<" ";
	  }
	gout<<endl;
      }
    }

    
    if (QHB2){
      gout<<"HB2"<<endl;
      gout<<"# UCoulomb : (jm1) (jm2) (jm3) (jm4)  Uc[jm1,jm2,jm3,jm4]"<<endl;
      gout.precision(8);
      ClebschGordan cg;
      gaunt_ck gck(l, cg);
      vector<double> FkoJ = ComputeCoulombRatios(l);
      function2D<function2D<double> > Uf;
      CoulombU0(Uf, l, ml_ms, one_electron_base, cg, gck, FkoJ);
      for (int j1=0; j1<Uf.size_N(); j1++){
	for (int j2=0; j2<Uf.size_Nd(); j2++){
	  for (int j3=0; j3<Uf.size_N(); j3++){
	    for (int j4=0; j4<Uf.size_Nd(); j4++){
	      double u = J_coulomb*Uf(j1,j2)(j3,j4);
	      if (fabs(u)>1e-6)
		gout<<setw(3)<<j1<<" "<<setw(3)<<j2<<" "<<setw(3)<<j3<<" "<<setw(3)<<j4<<" "<<setw(15)<<u<<endl;
	    }
	  }
	}
      }
      gout.precision(12);
    }else{
      gout<<"HB1"<<endl;
    }
    
    gout<<"# number of operators needed"<<endl;
    gout<<"1"<<endl;
    gout<<"# high frequency moment matrix elements, called Op1 above"<<endl;
    for (int ia=0; ia<nindex.size(); ia++){
      int i = nindex[ia];
      int in = istate[i].first;
      int wisize = nisize[i];
      //for (int ik=0; ik<deg.size(); ik++){
      for (int ik=0; ik<unique; ik++){
	gout<<setw(2)<<inindex[i]+1<<" ";
	function2D<double>& MEpm = FpFm[ik][i];
	function2D<double>& MEmp = FmFp[ik][i];
	gout<<setw(3)<<wisize<<" "<<setw(3)<<wisize<<"  ";
	for (int it=0; it<wisize; it++)
	  for (int jt=0; jt<wisize; jt++){
	    double pm = MEpm(it,jt);
	    if (in<=0 && ik<deg.size()){
	      pm = -MEmp(it,jt);
	      if (it==jt) pm += deg[ik];
	    }
	    gout<<setw(20)<<(fabs(pm)<1e-10 ? 0.0 : pm)<<" ";
	  }
	gout<<endl;
      }
      //gout<<setw(2)<<inindex[i]+1<<" "<<setw(3)<<wisize<<" "<<setw(3)<<wisize<<"  ";
      //for (int it=0; it<wisize; it++) for (int jt=0; jt<wisize; jt++) gout<<setw(20)<<0<<" ";
      //gout<<endl;
    }
    
    for (int i=0; i<istate.size(); i++) if (isize[i]>max_size) max_size = isize[i];

    gout<<"# Data for HB1"<<endl;
    gout<<1<<" "<<istate.size()<<" "<<one_electron_base.size()<<" "<<max_size<<endl;
    gout<<setw(2)<<"#"<<" "<<setw(3)<<"ind"<<" "<<setw(3)<<"N"<<" "<<setw(3)<<"K"<<" "<<setw(4)<<"Jz"<<" "<<setw(3)<<"size"<<endl;
    for (int i=0; i<istate.size(); i++){
      int ia = inindex[i];
      int in = istate[i].first;
      int iJz = istate[i].second;
      double Jz = Ebase[in][iJz].Jz;
      double n = Ebase[in][iJz].n;
      
      gout<<setw(2)<<i+1<<" "<<setw(3)<<ia+1<<" "<<setw(3)<<n<<" "<<setw(3)<<0<<" "<<setw(4)<<Jz<<" "<<setw(3)<<isize[i]<<"    ";
      for (int ij=0; ij<one_electron_base.size(); ij++){
	int jst = Fp_index[i][ij];
	int fi;
	if (jst>=istate.size() || jst<0) fi=0;
	else fi = jst/*-istart*/+1;
	gout<<setw(3)<<fi<<" ";
      }

      for (int l=0; l<isize[i]; l++){
	double Eatom_LS = Es[i][l];
	if (fabs(Eatom_LS)<1e-10) Eatom_LS=0;
	gout<<setw(18)<<Eatom_LS<<" ";
      }
      for (int l=0; l<isize[i]; l++) gout<<setw(18)<<Ebase[in][iJz].Jt_ein[l]<<" ";
      gout<<endl;
    }
    gout<<"# matrix elements"<<endl;
    for (int i=0; i<istate.size(); i++){
      int in = istate[i].first;
      int iJz = istate[i].second;
      double Jz = Ebase[in][iJz].Jz;
      double n = Ebase[in][iJz].n;
      for (int ij=0; ij<one_electron_base.size(); ij++){
	gout<<setw(2)<<i/*-istart*/+1<<" ";
	int jst = Fp_index[i][ij];
	int fi;
	if (jst>=istate.size() || jst<0) fi=0;
	else fi = jst/*-istart*/+1;
	gout<<setw(3)<<fi<<"  ";
	if (fi==0) { gout<<setw(3)<<0<<" "<<setw(3)<<0<<"  "<<endl;; continue; }
	function2D<double>& fp = Fp_matrix[i][ij];
	gout<<setw(3)<<isize[i]<<" "<<setw(3)<<isize[jst]<<"  ";
	for (int it=0; it<isize[i]; it++)
	  for (int jt=0; jt<isize[jst]; jt++){
	    double ff = fp(jt,it);
	    if (fabs(ff)<1e-10) ff=0;
	    gout<<setw(20)<<ff<<" ";
	  }
	gout<<endl;
      }
    }
    
  }
    
  return 0;
}
