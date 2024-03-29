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
#include "util.h"
#include "complex.h"
#include "blas.h"
#include "function.h"
#include "parser.h"
using namespace std;

namespace Parameters{
  Par<int>    nf              ("n",          "Valencies of the atom to be considered");
  Par<int>    l               ("l",      3,  "Orbital angular momentum of the shel");
  Par<double> J_coulomb       ("J",    0.0,  "Slatter integrals F2=J*11.9219653179191");
  Par<double> c_spinorb       ("cx",   0.0,  "Spin-orbit coupling");
  Par<double> E_exact         ("Ex",   0.5,  "Energy window treated exactly");
  Par<double> E_approx        ("Ep",   3.0,  "Energy window treated approximately");
  Par<double> E_oca           ("Eoca", 0.2,  "Energy window for OCA diagrams");
  Par<int>    qOCA            ("qOCA",   0,  "OCA diagrams are computed");
  Par<int>    qatom           ("qatom",  0,  "Prints full atomic energis rather than E-E_{ls}");
  Par<int>    susceptibility  ("qsusc",  0,  "Diagrams for susceptibility are computed");
  Par<int>    paramagnetic    ("para",   1,  "Runs in paramagnetic mode");
  Par<double> mOCA            ("mOCA", 1e-3, "Matrix element for OCA should be greater than that");
  Par<int>    kOCA            ("kOCA", -1,   "If set to >0, only diagrams which have kOCA pseudoparticle should be kept");
  Par<int>    Ncentral        ("Ncentral", -1, "OCA diagrams are selected such that central occupancy is in Ncentral");
  //  Par<int>    add_core        ("ncore",  0,  "Core diagrams are not added to the cix file");
}

using namespace std;

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
    //int ii = i + nn*(i1+nn*i2);
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
  double CG(double j, double m, double j1, double m1, double j2, double m2)
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
  vector<double> FkoJ(l+1);
  FkoJ[0]=0;//This is not defined
  if (l==1) FkoJ[1] = 5.;
  if (l==2) { FkoJ[1]=14./(1+0.625); FkoJ[2]=0.625*FkoJ[1];}
  if (l==3) { FkoJ[1]=6435./(286+195*0.668+250*0.494); FkoJ[2] = 0.668*FkoJ[1]; FkoJ[3]=0.494*FkoJ[1];}
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
  static int iJzmax=400;
  return iJz + iJzmax*i;
}
inline int ind_Jz_i(pair<int,int>& piJz)
{ return ind_Jz_i(piJz.first,piJz.second);}
pair<int,int> inline Jz_i(int ind)
{
  static int iJzmax=400;
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
  vector<int> Jkappa;
  friend void Compute_Fp(double j, double jz, const Eigenbase& E1, const Eigenbase& E2, const vector<int>& base, const operateLS& op, function2D<double>& mFp);
  friend void Compute_Gz(const Eigenbase& Ei, const vector<int>& base, const operateLS& op, function2D<double>& Gz);
public:
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
  {
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

void Compute_Gz(const Eigenbase& Ei, const vector<int>& base, const operateLS& op, function2D<double>& Gz)
{
  function2D<double> Gz0(Ei.reducedn.size(), Ei.reducedn.size());
  Gz0 = 0;
  for (int i=0; i<Ei.reducedn.size(); i++){
    int state = base[Ei.reducedn[i]];
    Gz0(i,i) = 0.5*(op.Lz(state)+2*op.Sz(state));
  }

  function2D<double> temp1(Ei.reducedn.size(), Ei.reducedn.size());
  temp1.Product(Gz0,Ei.T_NLzSz_Eig);
  Gz.MProduct(Ei.T_NLzSz_Eig,temp1);
}

void SimpleBubbles(const Eigenbase& E1, const Eigenbase& E2, const function2D<double>& Fp, map<int, map<int,double> >& bubbles)
{
  for (int i=0; i<E1.size(); i++){
    int i0 = ind_Jz_i(E1.iJz,i);
    for (int j=0; j<E2.size(); j++){
      int j0 = ind_Jz_i(E2.iJz,j);
      if (fabs(Fp(j,i))<1e-3) continue;
      bubbles[i0][j0] += sqr(Fp(j,i));
    }
  }
}

class OCAd{
public:
  vector<int> states;
  int ib1, ib2;
  double f;
  OCAd(int i0, int i1, int i2, int i3, int ib1_, int ib2_, double f_) : states(4), ib1(ib1_), ib2(ib2_), f(f_)
  { states[0] = i0;    states[1] = i1;    states[2] = i2;    states[3] = i3;  }
  bool operator==(const OCAd& oc)
  {
    bool eq  = (states[0]==oc.states[0] && states[1]==oc.states[1] && states[2]==oc.states[2] && states[3]==oc.states[3] && ib1==oc.ib1 && ib2==oc.ib2);
    bool sim = (states[0]==oc.states[0] && states[1]==oc.states[3] && states[2]==oc.states[2] && states[3]==oc.states[1] && ib1==oc.ib2 && ib2==oc.ib1);
    return eq || sim;
  }
};

bool cmpOCA(const OCAd& a, const OCAd& b)
{  return fabs(a.f)>fabs(b.f);}

void OcaDiagrams(const Eigenbase& E0, const Eigenbase& E1, const Eigenbase& E2, const Eigenbase& E3,
		 const function2D<double>& Fp1, const function2D<double>& Fp2, const function2D<double>& Fp3, const function2D<double>& Fp4,
		 list<OCAd>& diagrams, int ib1, int ib2, int in, svector<double> E_ground_state, double E_oca, double mOCA)
{
  for (int i0=0; i0<E0.size(); i0++){
    int il0 = ind_Jz_i(E0.iJz,i0);
    if (fabs(E0.eEner[i0]-E_ground_state[in])>E_oca) continue;
    for (int i1=0; i1<E1.size(); i1++){
      int il1 = ind_Jz_i(E1.iJz,i1);
      if (fabs(E1.eEner[i1]-E_ground_state[in+1])>E_oca) continue;
      for (int i2=0; i2<E2.size(); i2++){
	int il2 = ind_Jz_i(E2.iJz,i2);
	if (fabs(E2.eEner[i2]-E_ground_state[in+2])>E_oca) continue;
	for (int i3=0; i3<E3.size(); i3++){
	  int il3 = ind_Jz_i(E3.iJz,i3);
	  if (fabs(E3.eEner[i3]-E_ground_state[in+1])>E_oca) continue;
	  double v1 = Fp1(i1,i0);
	  double v2 = Fp2(i2,i1);
	  double v3 = Fp3(i2,i3);
	  double v4 = Fp4(i3,i0);
	  double v = v1*v2*v3*v4;
	  if (fabs(v)<0.5*mOCA) continue;
	  OCAd oca(il0,il1,il2,il3,ib1,ib2,v);
	  diagrams.push_back(oca);
	}
      }
    }
  }
}

void AddDiagram(const OCAd& ocad, list<OCAd>& list_oca)
{
  for (list<OCAd>::iterator p=list_oca.begin(); p!=list_oca.end(); p++)
    if ((*p)==ocad) {p->f += ocad.f; return;}
  list_oca.push_back(ocad);
}
void CompressOCA(const list<OCAd>& diagrams, map<int,int>& index0, map<int,int>& index1, map<int,int>& index2,
		 const vector<int>& cmp_one_electron_1, list<OCAd>& list_oca)
{
  for (list<OCAd>::const_iterator it0=diagrams.begin(); it0!=diagrams.end(); it0++){
    int i0 = index0[it0->states[0]];
    int i1 = index1[it0->states[1]];
    int i2 = index2[it0->states[2]];
    int i3 = index1[it0->states[3]];
    int ib1 = cmp_one_electron_1[it0->ib1];
    int ib2 = cmp_one_electron_1[it0->ib2];
    OCAd ocad(i0,i1,i2,i3,ib1,ib2,it0->f);
    AddDiagram(ocad, list_oca);
  }
}

void ChiBubbles(const Eigenbase& Eb, const function2D<double>& Gz, map<int, map<int,double> >& Chi_bubbles)
{
  for (int i=0; i<Eb.size(); i++){
    int i0 = ind_Jz_i(Eb.iJz,i);
    for (int j=0; j<Eb.size(); j++){
      int j0 = ind_Jz_i(Eb.iJz,j);
      if (fabs(Gz(j,i))<1e-3) continue;
      Chi_bubbles[i0][j0] += sqr(Gz(j,i));
    }
  }
}

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


void CompressBubbles(map<int, map<int,double> >& bubbles, map<int,int>& indexn, map<int,int>& indexp, vector<vector<double> >& nca)
{
  for (map<int,map<int,double> >::iterator it0=bubbles.begin(); it0!=bubbles.end(); it0++){
    int i0 = it0->first;
    int i = indexn[i0];
    if (i>=nca.size()) cerr<<"First index is out of range in CompessBubbles!"<<endl;
    for (map<int,double>::iterator jt0=bubbles[i0].begin(); jt0!=bubbles[i0].end(); jt0++){
      int j0 = jt0->first;
      int j = indexp[j0];
      if (j>=nca[i].size()) cerr<<"Second index is out of range in CompessBubbles!"<<endl;
      nca[i][j] += jt0->second;
    }
  }
}

void CompressBubbles_Gz(map<int, map<int,double> >& bubbles, map<int,int>& indexn, vector<vector<double> >& nca)
{
  for (map<int,map<int,double> >::iterator it0=bubbles.begin(); it0!=bubbles.end(); it0++){
    int i0 = it0->first;
    int i = indexn[i0];
    if (i>=nca.size()) cerr<<"First index is out of range in CompessBubbles_Gz!"<<endl;
    for (map<int,double>::iterator jt0=bubbles[i0].begin(); jt0!=bubbles[i0].end(); jt0++){
      int j0 = jt0->first;
      int j = indexn[j0];
      if (j>=nca[i].size()) cerr<<"Second index is out of range in CompessBubbles_Gz!"<<endl;
      nca[i][j] += jt0->second;
    }
  }
}

bool CmpStates_(const pair<int,double>& a, const pair<int,double>& b)
{return fabs(a.second)>fabs(b.second); }

class CmpE{
  map<int, vector<double> >& Es;
public:
  CmpE(map<int, vector<double> >& Es_) : Es(Es_) {};
  bool operator()(int a, int b){
    pair<int,int> pa = Jz_i(a), pb = Jz_i(b);
    double Ea = Es[pa.first][pa.second];
    double Eb = Es[pb.first][pb.second];
    return Ea<Eb;
  }
};

class TiJz{
public:
  int n;
  double Jmax;
  int iJmaxi;
  
  TiJz(int n_, int l) : n(n_)
  {
    Jmax = l*n+0.5*(n%2);
    iJmaxi = static_cast<int>(2*(Jmax+2*l+1));
  }
  int give_iJz(double Jz)
  {return static_cast<int>(2*Jz) + iJmaxi;}
  double give_Jz(int iJz)
  { return (iJz-iJmaxi)/2.;}
  int iJz_min(int Np)
  { return  (Np)%2 + iJmaxi; }
  double J_max(){return Jmax;}
};

void RememberParams (int argc, char *argv[]){
  ofstream param ("history.atom", ios::app);
  if (!param) cerr<<" Didn't suceeded to open params file!"<<"history.atom"<<endl;
  for (int i=0; i<argc; i++) param << argv[i] << " ";
  param << endl;
}

int main(int argc, char *argv[])
{
  setvbuf (stdout, NULL, _IONBF, BUFSIZ);
  using namespace Parameters;
  DynArray arr(100, &nf, &l, &J_coulomb, &c_spinorb, &qOCA, &E_oca, &E_exact, &E_approx, &qatom, &susceptibility, &paramagnetic, &Ncentral, &mOCA, &kOCA, NULL);

  if (argc<2) {
    arr.printl(clog);
    return 0;
  }

  arr.ParsComLine(argc, argv);
  arr.print(clog);
  
// int main(int argc, char *argv[], char *env[])
// {
//   //Sergej's program // Jur. Chem. Phys. 40, 3428, (1964)
//   double J_coulomb = 0.60;//0.680290;//0.604
//   double c_spinorb = 0.31;//0.367864976619714;//0.272116;//0.316
//   int n=1;
//   int l=3;
//   bool susceptibility = false;
//   bool qOCA = false;
//   bool add_core = true;
//   bool paramagnetic = false;
//   //Energy window of excited states to keep!
//   double E_exact=0.5, E_approx=3.; //long: (0.5,3.0), long_long: (0.8,5), original: (1e-4,3)
//   double E_oca = 1e-3;
//   bool display_help = false;
//   bool qatom = false;
//   int i=0;

    
//   if (argc<=1) {display_help = true; i=-1;}
//   while (++i<argc || display_help){
//     std::string str(argv[i]);
//     if (str=="-J" && i<argc-1) J_coulomb = atof(argv[++i]);
//     if (str=="-cx" && i<argc-1) c_spinorb = atof(argv[++i]);
//     if (str=="-n" && i<argc-1) n = atoi(argv[++i]);
//     if (str=="-l" && i<argc-1) l = atoi(argv[++i]);
//     if (str=="-Ex" && i<argc-1) E_exact = atof(argv[++i]);
//     if (str=="-Ep" && i<argc-1) E_approx = atof(argv[++i]);
//     if (str=="-Eoca" && i<argc-1) E_oca = atof(argv[++i]);
//     if (str=="-qsusc") susceptibility = true;
//     if (str=="-ncore") add_core = false;
//     if (str=="-para") paramagnetic = true;
//     if (str=="-qOCA") qOCA = true;
//     if (str=="-qatom") qatom = true;
//     if (str=="-h" || str=="--help" || display_help){
//       std::clog<<"**      Exact diagonalization of the atom           **\n";
//       std::clog<<"**                                                  **\n";
//       std::clog<<"**      Copyright Kristjan Haule, 26.09.2004        **\n";
//       std::clog<<"******************************************************\n";
//       std::clog<<"\n";
//       std::clog<<"atom [Options]\n" ;
//       std::clog<<"Options:   -J     Slatter integrals F2=J*11.9219653179191 ("<<J_coulomb<<")\n";
//       std::clog<<"           -cx    Spin-orbit coupling ("<<c_spinorb<<")\n";
//       std::clog<<"           -n     Average occupancy of the atom ("<<n<<")\n";
//       std::clog<<"           -l     Orbital angular momentum of the shel ("<<l<<")\n";
//       std::clog<<"           -Ex    Energy window treated exactly ("<<E_exact<<")\n";
//       std::clog<<"           -Ep    Energy window treated approximately ("<<E_approx<<")\n";
//       std::clog<<"           -Eoca  Energy window for OCA diagrams ("<<E_oca<<")\n";
//       std::clog<<"           -susc  Diagrams for susceptibility are computed ("<<susceptibility<<")\n";
//       std::clog<<"           -ncore Core diagrams are not added to the cix file ("<<add_core<<")\n";
//       std::clog<<"           -para  Runs in paramagnetic mode ("<<paramagnetic<<")\n";
//       std::clog<<"           -qOCA  OCA diagrams are computes ("<<qOCA<<")\n";
//       std::clog<<"           -qatom Prints full atomic energis rather than E-E_{ls} ("<<qatom<<")\n";
//       std::clog<<"*****************************************************\n";
//       return 0;
//     }
//   }


  RememberParams (argc, argv);
  
  deque<int> Np;
  while (nf.next()) Np.push_back(nf);
  cout<<"n size ="<<Np.size()<<endl;
  for (int i=0; i<Np.size(); i++) cout<<"n: "<<Np[i]<<endl;

  deque<int> qNcentral;
  if (Ncentral>=0){
    while (Ncentral.next()) qNcentral.push_back(Ncentral);
    cout<<"Ncentral = ";
    for (int i=0; i<qNcentral.size(); i++) cout<<qNcentral[i]<<" ";
    cout<<endl;
  }
  deque<int> qkOCA;
  if (kOCA>0){
    while (kOCA.next()) qkOCA.push_back(kOCA);
    cout<<"kOCA = ";
    for (int i=0; i<qkOCA.size(); i++) cout<<qkOCA[i]<<" ";
    cout<<endl;
  }
  
  deque<double> e_approx;
  while (E_approx.next()) e_approx.push_back(E_approx);
  vector<double> _E_approx(Np.size());
  for (int i=0; i<min(e_approx.size(),Np.size()); i++) _E_approx[i] = e_approx[i];
  for (int i=e_approx.size(); i<Np.size(); i++) _E_approx[i] = e_approx[e_approx.size()-1];
  
  deque<double> e_exact;
  while (E_exact.next()) e_exact.push_back(E_exact);
  vector<double> _E_exact(Np.size());
  if (paramagnetic){
    for (int i=0; i<Np.size(); i++) _E_exact[i]=-1; // for paramagnetic state, we can always combine states with different Jz and same J
  }else{
    for (int i=0; i<min(e_exact.size(),Np.size()); i++) _E_exact[i] = e_exact[i];
    for (int i=e_exact.size(); i<Np.size(); i++) _E_exact[i] = e_exact[e_exact.size()-1];
  }

  for (int i=0; i<_E_exact.size(); i++){
    cout<<"Eexact ="<<_E_exact[i]<<endl;
    cout<<"Eapprox="<<_E_approx[i]<<endl;
  }

  ClebschGordan cg;
  int dim = 2*(2*l+1);
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
  

  // THIS CAN BE PROBLEMATIC. CHECK THAT IT ALWAYS WORKS
  TiJz tiJz(Np[Np.size()/2],l); // How to convert Jz to array index
  
  int nstart=0, nstop=Np.size();
  //  if (add_core){ nstart=-1; nstop=4;}
  
  //  svector<double> _E_exact(nstart,nstop), _E_approx(nstart,nstop);
  //  for (int in=nstart; in<0; in++) {_E_exact[in]=-1.;     _E_approx[in]=2;}// in core, states are not treated exactly
  //  for (int in=0; in<3; in++)      {_E_exact[in]=E_exact; _E_approx[in]=E_approx;}
  //  for (int in=3; in<nstop; in++)  {_E_exact[in]=-1;      _E_approx[in]=2;}
  
  //  svector<int> Np(nstart, nstop);
  //  for (int i=nstart; i<nstop; i++) Np[i] = n+i-1;

  svector<map<int, vector<int> > > indexes(nstart,nstop);
  svector<map<int, vector<double> > > Es(nstart,nstop), E_LSs(nstart,nstop); // atom energy and Spin-Orbit energy
  svector<map<int, vector<vector<int> > > > JJzLSks(nstart,nstop);

  ///////////////////// Creates One electron basis for the bath function /////////////////////////////////
  vector<pair<double,double> > one_electron_base(2*(2*l+1));
  {
    int kk=0;
    for (int sz=-1; sz<=1; sz+=2){
      double j = l-0.5*sz;
      int j2 = static_cast<int>(2*j);
      for (double j2z=-j2; j2z<=j2; j2z+=2){
  	double jz = j2z/2.;
  	one_electron_base[kk++] = make_pair(j,jz);
      }
    }
  }
  
  ///////////////////// Data structure to store results /////////////////////////////////
  svector<vector<map<int, map<int,double> > > > bubbles(nstart,nstop); // largest data structure containing all possible bubbles between states
  svector<vector<map<int, double> > > fpfm(nstart,nstop), fmfp(nstart,nstop);
  for (int in=fpfm.start(); in<fpfm.stop(); in++){
    fpfm[in].resize(one_electron_base.size());
    fmfp[in].resize(one_electron_base.size());
    bubbles[in].resize(one_electron_base.size());
  }
  svector<map<int, map<int,double> > > Chi_bubbles(nstart,nstop); 
  int instart=nstart, instop=nstop;
  svector<list<OCAd> > OCAdiag(instart,instop-2);

  ///////////////////// Goes over all atomic states and diagonalizes blocks of H /////////////////////////////////
  {
    // Basic states for all necessary occupations and spins are created
    svector<vector<Eigenbase> > Ebase(nstart,nstop);
    for (int in=nstart; in<nstop; in++) {
      int tn = Np[in];
      double tJmax = tn*l + 0.5*(tn%2);
      int itJmax = tiJz.give_iJz(tJmax)+1;
      Ebase[in].resize(itJmax);
      double Jz = -tJmax;
      int tot_size=0;
      do{
	int iJz = tiJz.give_iJz(Jz);//static_cast<int>(2*Jz) + iJmaxi;
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

	// To save for later use
	indexes[in][iJz] = Ebase[in][iJz].index2;
	Es[in][iJz]      = Ebase[in][iJz].eEner;
	JJzLSks[in][iJz] = Ebase[in][iJz].JJzLSk;
	E_LSs[in][iJz]   = Ebase[in][iJz].E_LS;

	tot_size += Ebase[in][iJz].size();
      } while(++Jz<=tJmax);
      cout<<"Total size = "<<tot_size<<endl;
    }
    
  ///////////////////// Computes all possible bubbles between atomic states, also F^+F and FF^+ are computed ///////////////////
    for (int in=nstart; in<nstop-1; in++){
      double Jmax = l*Np[in] + 0.5*(Np[in]%2);
      double Jz = -Jmax;
      do{
	int iJzn = tiJz.give_iJz(Jz);
	if (iJzn>=Ebase[in].size() || Ebase[in][iJzn].size()==0) continue;
	for (int ij=0; ij<one_electron_base.size(); ij++){
	  double j = one_electron_base[ij].first;
	  double jz = one_electron_base[ij].second;
	  cout<<"Constructing diagrams for Jz="<<Jz<<" j="<<j<<" jz="<<jz<<" ---------------"<<endl;
	  int iJzp = tiJz.give_iJz(Jz + jz);
	  if (iJzp>=Ebase[in+1].size() || Ebase[in+1][iJzp].size()==0) continue;
	  
	  function2D<double> Fp(Ebase[in+1][iJzp].size(), Ebase[in][iJzn].size());
	  Compute_Fp(j, jz, Ebase[in][iJzn], Ebase[in+1][iJzp], base, op, Fp); // This computes (F^{j,jz,+})_{n+1,n}
	  SimpleBubbles(Ebase[in][iJzn], Ebase[in+1][iJzp], Fp, bubbles[in][ij]); // Enumerates bubbles according to the "global" index
	  Set_fmfp(Ebase[in][iJzn], Ebase[in+1][iJzp], Fp, fpfm[in+1][ij], fmfp[in][ij]);// Stores F.F^+ and F^+.F
	}
      } while(++Jz<=Jmax);
    }
    
    ///////////////// Second order diagrams ///////////////////////////////////
    svector<double> E_ground_state(nstart,nstop);
    for (int in=nstart; in<nstop; in++) E_ground_state[in] = Es[in][tiJz.iJz_min(Np[in])][0];

    if (qOCA){
      for (int in=instart; in<instop-2; in++){
	for (int iJz0=0; iJz0<Ebase[in].size(); iJz0++){// state with the lowest occupancy |0>
	  if (Ebase[in][iJz0].size()==0) continue;
	  double Jz0 = tiJz.give_Jz(iJz0);

	  if (fabs(Ebase[in][iJz0].eEner[0]-E_ground_state[in])>E_oca) continue; // only between ground states
      
	  for (int ij1=0; ij1<one_electron_base.size(); ij1++){
	    double j1 = one_electron_base[ij1].first;
	    double jz1 = one_electron_base[ij1].second;
	    double Jz1 = Jz0 + jz1;
	    int iJz1 = tiJz.give_iJz(Jz1);
	  
	    cout<<"Constructing OCA diagrams for Jz="<<Jz0<<" j1="<<j1<<" jz1="<<jz1<<" ---------------"<<endl;
	    if (iJz1>=Ebase[in+1].size() || Ebase[in+1][iJz1].size()==0) continue; // base for state |1> empty?
	  
	    if (fabs(Ebase[in+1][iJz1].eEner[0]-E_ground_state[in+1])>E_oca) continue; // only between ground states
	
	    function2D<double> Fp1(Ebase[in+1][iJz1].size(), Ebase[in][iJz0].size());
	    Compute_Fp(j1, jz1, Ebase[in][iJz0], Ebase[in+1][iJz1], base, op, Fp1); // This computes (F^{j1,jz1,+})_{1,0}
	  
	    for (int ij2=0; ij2<one_electron_base.size(); ij2++){
	      double j2 = one_electron_base[ij2].first;
	      double jz2 = one_electron_base[ij2].second;
	      double Jz2 = Jz1 + jz2;
	      int iJz2 = tiJz.give_iJz(Jz2);
	      double Jz3 = Jz2 - jz1;
	      int iJz3 = tiJz.give_iJz(Jz3);


	      if (iJz2>=Ebase[in+2].size() || Ebase[in+2][iJz2].size()==0) continue;
	      if (iJz3>=Ebase[in+1].size() || Ebase[in+1][iJz3].size()==0) continue;

	      if (fabs(Ebase[in+2][iJz2].eEner[0]-E_ground_state[in+2])>E_oca) continue; // only between ground states
	      if (fabs(Ebase[in+1][iJz3].eEner[0]-E_ground_state[in+1])>E_oca) continue; // only between ground states

	      bool QNcentral = qNcentral.size()==0; // if Ncentral is -1 in the command line, we do nothing.
	      for (int iq=0;  iq<qNcentral.size(); iq++)
		if (qNcentral[iq]==Np[in+1]) QNcentral = true;
	      if (!QNcentral) continue;

	      function2D<double> Fp2(Ebase[in+2][iJz2].size(), Ebase[in+1][iJz1].size());
	      Compute_Fp(j2, jz2, Ebase[in+1][iJz1], Ebase[in+2][iJz2], base, op, Fp2); // This computes (F^{j2,jz2,+})_{2,1}
	      function2D<double> Fp3(Ebase[in+2][iJz2].size(), Ebase[in+1][iJz3].size());
	      Compute_Fp(j1, jz1, Ebase[in+1][iJz3], Ebase[in+2][iJz2], base, op, Fp3); // This computes (F^{j1,jz1,+})_{2,1}
	      function2D<double> Fp4(Ebase[in+1][iJz3].size(), Ebase[in][iJz0].size());
	      Compute_Fp(j2, jz2, Ebase[in][iJz0], Ebase[in+1][iJz3], base, op, Fp4); // This computes (F^{j2,jz2,+})_{1,0}
	    
	      OcaDiagrams(Ebase[in][iJz0],Ebase[in+1][iJz1],Ebase[in+2][iJz2],Ebase[in+1][iJz3], Fp1, Fp2, Fp3, Fp4, OCAdiag[in], ij1, ij2, in, E_ground_state, E_oca, mOCA);
	    }
	  }
	}
      }
    }
    ///////////////// Second order diagrams ///////////////////////////////////

    {
      for (int in=nstart; in<nstop; in++){
	int tn = Np[in];
	double tJmax = tn*l + 0.5*(tn%2);
	double Jz = -tJmax;
	do{
	  int iJz = tiJz.give_iJz(Jz);
	  if (Ebase[in].size()>iJz && Ebase[in][iJz].size()!=0){
	    clog<<" Constructing G "<<in<<" "<<Jz<<endl;
	    function2D<double> Gz(Ebase[in][iJz].size(), Ebase[in][iJz].size());
	    Compute_Gz(Ebase[in][iJz], base, op, Gz);
	    ChiBubbles(Ebase[in][iJz], Gz, Chi_bubbles[in]);
	  }
	}while(++Jz<=tJmax);
      }
    }
  }

  ///////////////////// Here we compress states and bubbles. Ground states and some ////////////////////
  ///////////////////// excited states are kept. First we create index array indexp ////////////////////
  ///////////////////// which stores index for all ground states                    ////////////////////
  // This part of the code stores to indexp all ground states (states with the same ground-state energy
  // but different Jz).
  svector<deque<int> > indexp(nstart,nstop);
  for (int in=nstart; in<nstop; in++){
    int iJz_min = tiJz.iJz_min(Np[in]);
    cout<<" ------------------ in="<<in<<" -----------------"<<endl;
    if (Es[in].find(iJz_min)==Es[in].end()) cerr<<"This is not right!"<<endl;
    for (int i = 0; i< Es[in][iJz_min].size(); i++) cout<<"Energy"<<i<<" ="<<Es[in][iJz_min][i]<<endl;
  
    double E_ground_state = Es[in][iJz_min][0];
    // All the ground states are kept first
    for (map<int,vector<int> >::iterator p=indexes[in].begin(); p!=indexes[in].end(); p++){// over Jz
      int iJz = p->first;
      double Jz = tiJz.give_Jz(iJz);
      if (!Es[in][iJz].size()) continue;
      if (fabs(Es[in][iJz][0]-E_ground_state)<1e-4){
  	int vind = (p->second)[0]; // index to the ground-state
  	indexp[in].push_back(vind);
  	clog<<"Jz="<<Jz<<" Energy="<<Es[in][iJz][0]<<" index="<<vind<<endl;
      }
    }
  }

  //////////////// Here we add to the indexp all states which are coupled to  ////////////////////////
  //////////////// any ground state through a bubble -                        ////////////////////////
  ///////////////  all states which contribute to atomic ground state         ////////////////////////
  // This part of the code adds to the indexp all states that have nonzero bubbles with
  // any ground state - only those states will contribute to atomic Green's function at low
  // temperature
  svector<deque<int> > indexp_orig(indexp);
  for (int in=nstart; in<nstop-1; in++){
     
    for (int i=0; i<indexp_orig[in].size(); i++){
      int igs = indexp_orig[in][i];
      for (int ij=0; ij<one_electron_base.size(); ij++){
 	for (map<int,double>::iterator it = bubbles[in][ij][igs].begin(); it!=bubbles[in][ij][igs].end(); it++){
 	  bool not_found = find(indexp[in+1].begin(),indexp[in+1].end(),it->first)==indexp[in+1].end();
 	  if (fabs(it->second)>1e-3 && not_found)
 	    indexp[in+1].push_back(it->first);
 	}
      }
    }
    for (int i=0; i<indexp_orig[in+1].size(); i++){
      int igs = indexp_orig[in+1][i];
      for (int ij=0; ij<one_electron_base.size(); ij++){
 	for (map<int, map<int,double> >::iterator it = bubbles[in][ij].begin(); it!=bubbles[in][ij].end(); it++){
 	  map<int,double>& bub = it->second;
 	  //	  for (map<int,double>::iterator iw=bub.begin(); iw!=bub.end(); iw++) clog<<iw->first<<" "<<iw->second<<endl;
 	  if (bub.find(igs)==bub.end()) continue; // it is not coupled to the ground state
 	  int j = it->first;
 	  bool not_found = find(indexp[in].begin(),indexp[in].end(),j)==indexp[in].end();
 	  if (fabs(bub[igs])>1e-3 && not_found)	indexp[in].push_back(j); // add the guy if it is connected to the ground state and is not yet in
 	}
      }
    }
  }
  // States in indexp are sorted according to atomic energy
  // HERE IS A BUG
  for (int in=nstart; in<nstop; in++){
    cout<<"******** in="<<in<<" *********"<<endl;
    CmpE cmp(Es[in]);
    sort(indexp[in].begin(), indexp[in].end(), cmp);

    cout<<"indexp="<<endl;
    for (int i=0; i<indexp[in].size(); i++){
      cout<<setw(3)<<i<<" "<<setw(5)<<indexp[in][i]<<endl;
    }
  }
 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///// Here we create another index array (indexf) which is array of pseudoparticles:
  ///// indexf[in]={{q0,q1,q2},{q3,q4,...},...}; here "in" is occupation while {q0,q1,q2} are atomic states
  ///// represented with first pseudoparticle (they are usually degenerate)
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///// Here we create pseudoparticles using indexp array
  svector<deque<deque<int> > > indexf(nstart,nstop);
  for (int in=nstart; in<nstop; in++){
    pair<int,int> p_gs = Jz_i(indexp[in][0]);
    double Egs = Es[in][p_gs.first][p_gs.second];
    double Et;
    int i;
    for (i=0; i<indexp[in].size(); i++){// States in the range [Egs..Egs+E_exact] treat exactly: Each state has its own pseudoparticle
      pair<int,int> p = Jz_i(indexp[in][i]);
      Et = Es[in][p.first][p.second];
      if ((Et-Egs)>_E_exact[in]) break;
      indexf[in].push_back(deque<int>(1, indexp[in][i]));// queue has only one element  
    }
    deque<int> tque;
    int j;
    for (j=i; j<indexp[in].size(); j++){// States in the range [Egs+E_exact...Egs+E_approx] treat approximately: States with same Jz correspond to the same pseudoparticle
      pair<int,int> p = Jz_i(indexp[in][j]);
      double E = Es[in][p.first][p.second];
      if (fabs(E-Et)<1e-4) {tque.push_back(indexp[in][j]); continue;}// queue containe elements which have the same energy
      Et = E;                 // Energy changes, update energy
      indexf[in].push_back(tque); // Store previous set of states
      tque.clear();           // And prepare new queue
      tque.push_back(indexp[in][j]);// Stores the new state
      if ((E-Egs)>_E_approx[in]) break;  // If state is to high in energy, exits
    }
    for (int k=j+1; k<indexp[in].size(); k++){// The rest of the states are punched together into one state
      tque.push_back(indexp[in][k]);
    }
    if (tque.size()) indexf[in].push_back(tque); 
  }
  //// States contained in indexp are remembered in an STL set 
  svector<set<int> > already_sorted_states(nstart,nstop);
  for (int in=nstart; in<nstop; in++)
    for (int i=0; i<indexp[in].size(); i++) already_sorted_states[in].insert(indexp[in][i]);
   
  /// The rest of the states not in indexp are added to indexf array (containing pseudoparticles)
  for (int in=nstart; in<nstop; in++){
    deque<int> tque;
    for (map<int,vector<int> >::iterator p=indexes[in].begin(); p!=indexes[in].end(); p++){// over Jz
      int iJz = p->first;
      double Jz = tiJz.give_Jz(iJz);//(iJz-iJmaxi)/2.;
      vector<int>& tindex2 = p->second;
      for (int j=0; j<tindex2.size(); j++){
   	if (already_sorted_states[in].find(tindex2[j])==already_sorted_states[in].end())// this particular state is not yet in the indexf
   	  tque.push_back(tindex2[j]);
      }
    }
    if (tque.size()) indexf[in].push_back(tque); 
  }


  /////////////// The order in which states of different occupation are printed //////////////////
  /////////////// first the valence states between 0-2 and latter core states   //////////////////
  /////////////// with in=-1 and in=3                                           //////////////////
  vector<int> nord(nstop-nstart);
  for (int i=0; i<nord.size(); i++) nord[i]=i;
  //  if (add_core) {nord[3]=3;nord[4]=-1;}
  ////////////// Pseudoparticles are enumerated such that first come valence   ////////////////////
  ////////////// states and later core states                                  ////////////////////
  svector<vector<int> > pseudo(nstart,nstop);
  for (int in=nstart; in<nstop; in++) pseudo[in].resize(indexf[in].size()); // number of pseudoparticles
  ///////////// Also inverse of indexf is created which is used to re-enumerate ///////////////////
  ///////////// and compress bubbles and states (energy, degeneracy,....)       ///////////////////
  svector<map<int,int> > indexf_1(nstart,nstop);
  int pseudo_number=0;
  for (int iin=0; iin<nstop-nstart; iin++){
    int in = nord[iin];
    for (int i=0; i<indexf[in].size(); i++){
      for (int j=0; j<indexf[in][i].size(); j++){
	int indf = indexf[in][i][j];
	indexf_1[in][indf] = pseudo_number; // the pseudoparticle number is i
      }
      pseudo[in][i]=pseudo_number;
      pseudo_number++;
    }
  }

  /////////// Printing of pseudoparticles ///////////////////////////////////////////////////////
  cout<<"-------------------------------------------------------------------------------------"<<endl;
  cout<<"-------------------------------------------------------------------------------------"<<endl;
  for (int inn=0; inn<nstop-nstart; inn++){
    int in = nord[inn];
    cout<<"in="<<in<<endl;
    for (int i=0; i<indexf[in].size(); i++){
      cout<<"------------- New pseudoparticle ---------------------------------------------------------"<<endl;
      for (int j=0; j<indexf[in][i].size(); j++){
   	pair<int,int> p = Jz_i(indexf[in][i][j]);
	int ii = indexf[in][i][j];
	int jj = indexf_1[in][ii];
   	cout<<setw(4)<<jj<<" "<<setw(20)<<Es[in][p.first][p.second]<<" "<<setw(5)<<p.second<<" "<<setw(12)<<tiJz.give_Jz(p.first)<<" "<<setw(10)<<indexf[in][i][j]<<endl;
      }
    }
  }
  for (int in=nstart; in<nstop; in++){
    int isum1=0;
    for (int i=0; i<indexf[in].size(); i++) isum1 += indexf[in][i].size();
    int isum2=0;
    for (map<int,vector<int> >::iterator p=indexes[in].begin(); p!=indexes[in].end(); p++){
      int iJz = p->first;
      isum2 += indexes[in][iJz].size();
    }
    cout<<"Number of all states="<<isum1<<","<<isum2<<"  while the number of pseudoparticles is "<<indexf[in].size()<<endl;
  }
  /////////// Printing of pseudoparticles ///////////////////////////////////////////////////////


  deque<deque<int> > cmp_one_electron;
  vector<int> cmp_one_electron_1(one_electron_base.size());
  {
    int No = one_electron_base.size();
    double j  = one_electron_base[No-1].first;
    double jz = one_electron_base[No-1].second;
    deque<int> tque;
    for (int i=No-1; i>=0; i--){
      double tj = one_electron_base[i].first;
      double tjz = one_electron_base[i].second;
      if (fabs(tj-j)<1e-4 && (paramagnetic || fabs(tjz-jz)<1e-4)) {tque.push_back(i); continue;}
      j = tj; jz = tjz;
      cmp_one_electron.push_back(tque);
      tque.clear();
      tque.push_back(i);
    }
    cmp_one_electron.push_back(tque);
    for (int i=0; i<cmp_one_electron.size(); i++)
      for (int j=0; j<cmp_one_electron[i].size(); j++)
	cmp_one_electron_1[cmp_one_electron[i][j]] = i;
    
//     cout<<"One electron base: "<<endl;
//     for (int i=0; i<cmp_one_electron_1.size(); i++)
//       cout<<i<<" "<<cmp_one_electron_1[i]<<endl;
  }
  /////////////////////// Here we compress bubbles and states using ////////////////////////////
  ////////////////////// above created indexf and indexf_1 array    ////////////////////////////
  
  vector<vector<vector<double> > > nca(cmp_one_electron.size());
  for (int ij=0; ij<cmp_one_electron.size(); ij++){
    nca[ij].resize(pseudo_number);
    for (int i=0; i<pseudo_number; i++) {
      nca[ij][i].resize(pseudo_number);
      for (int j=0; j<pseudo_number; j++) nca[ij][i][j]=0;
    }
  }
  vector<vector<double> > nca_Gz(pseudo_number);
  for (int i=0; i<pseudo_number; i++){
    nca_Gz[i].resize(pseudo_number);
    for (int j=0; j<pseudo_number; j++){
      nca_Gz[i][j]=0;
    }
  }
  
  for (int in=nstart; in<nstop-1; in++)
    for (int ij=0; ij<one_electron_base.size(); ij++)
      CompressBubbles(bubbles[in][ij], indexf_1[in], indexf_1[in+1], nca[cmp_one_electron_1[ij]]);

  for (int in=0; in<3; in++)
    CompressBubbles_Gz(Chi_bubbles[in], indexf_1[in], nca_Gz);

  list<OCAd> list_oca;
  if (qOCA){
    for (int in=instart; in<instop-2; in++)
      CompressOCA(OCAdiag[in], indexf_1[in], indexf_1[in+1], indexf_1[in+2], cmp_one_electron_1,list_oca);

    cout<<"OCA diagrams :"<<endl;
    for (list<OCAd>::const_iterator p=list_oca.begin(); p!=list_oca.end(); p++){
      cout<<setw(3)<<p->states[0]<<" "<<setw(3)<<p->states[1]<<" "<<setw(3)<<p->states[2]<<" "<<setw(3)<<p->states[3]<<"    ";
      cout<<setw(3)<<p->ib1<<" "<<setw(3)<<p->ib2<<"    "<<setw(20)<<-p->f<<endl;
    }
    cout<<"-----------------"<<endl;
  }


  
  //////////////////// Checking the constraint (FF^+)_{ii} which should be      ///////////////
  /////////////////// sum over all states withing pseudoparticle (1-occupancy)  ///////////////
  for (int ij=0; ij<cmp_one_electron.size(); ij++){
    for (int i=0; i<pseudo_number; i++){
      double sum=0;
      for (int j=0; j<pseudo_number; j++) sum += nca[ij][i][j];
      cout<<setw(5)<<sum<<" ";
    }
    cout<<endl;
  }

  // For the lowest N, one need to use 1-F F^+ to computer  F^+F 
  for (int ij=0; ij<one_electron_base.size(); ij++){
    for (map<int,vector<int> >::iterator p=indexes[nstart].begin(); p!=indexes[nstart].end(); p++){// over Jz
      int iJz = p->first;
      double Jz = tiJz.give_Jz(iJz);
      //if (!Es[nstart][iJz].size()) continue;
      vector<int>& tindex2 = p->second;
      for (int j=0; j<tindex2.size(); j++){
	int ikk = tindex2[j];
	fpfm[nstart][ij][ikk] = 1-fmfp[nstart][ij][ikk];
      }
    }
  }

  
  
  {
    cout<<"Gz: "<<endl;
    for (int i=0; i<pseudo_number; i++){
      double sum=0;
      for (int j=0; j<pseudo_number; j++) sum += nca_Gz[j][i];
      cout<<setw(5)<<i<<" "<<sum<<" "<<endl;
    }
  }
  /////////////// Calculation of Energy, degeneracy, SO-energy and occupation   ///////////////
  /////////////// for all new created pseudoparticles                           ///////////////
  svector<vector<double> > Energ(nstart,nstop), Deg(nstart,nstop), E_LS(nstart,nstop);
  svector<vector<vector<double> > > vfpfm(nstart,nstop);
  for (int in=nstart; in<nstop; in++){
    Energ[in].resize(pseudo[in].size());
    Deg[in].resize(pseudo[in].size());
    E_LS[in].resize(pseudo[in].size());
    for (int i=0; i<indexf[in].size(); i++){
      double averE=0;
      double deg=0;
      double aver_E_LS=0;
      for (int j=0; j<indexf[in][i].size(); j++){
	pair<int,int> p = Jz_i(indexf[in][i][j]);
	averE += Es[in][p.first][p.second];
	aver_E_LS += E_LSs[in][p.first][p.second];
	deg++;
      }
      averE /= indexf[in][i].size();
      aver_E_LS /= indexf[in][i].size();
      Energ[in][i] = averE;
      E_LS[in][i] = aver_E_LS;
      Deg[in][i] = deg;
    }
    
    vfpfm[in].resize(cmp_one_electron.size());
    for (int ij=0; ij<vfpfm[in].size(); ij++){
      vfpfm[in][ij].resize(pseudo[in].size());
      for (int i=0; i<vfpfm[in][ij].size(); i++) vfpfm[in][ij][i]=0;
    }
	
    for (int ij=0; ij<one_electron_base.size(); ij++){
      for (int i=0; i<indexf[in].size(); i++){
	for (int j=0; j<indexf[in][i].size(); j++){
	  int ind = indexf[in][i][j];
	  if (fpfm[in][ij].find(ind)!=fpfm[in][ij].end()){
	    vfpfm[in][cmp_one_electron_1[ij]][i] += fpfm[in][ij][ind];
	  }
	}
      }
    }
      
    for (int ij=0; ij<vfpfm[in].size(); ij++)
      for (int i=0; i<vfpfm[in][ij].size(); i++) if(fabs(vfpfm[in][ij][i])<1e-6) vfpfm[in][ij][i]=0;
  }


  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////// Bubbles are finally reordered in a way convenient for printing below  //////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  vector<vector<deque<pair<int,double> > > > ncab(cmp_one_electron.size()), ncaf(cmp_one_electron.size());
  for (int ij=0; ij<cmp_one_electron.size(); ij++){
    ncaf[ij].resize(nca[ij].size());
    ncab[ij].resize(nca[ij].size());
    for (int i=0; i<nca[ij].size(); i++){
      deque<pair<int,double> > tncab;
      for (int j=0; j<nca[ij][i].size(); j++){
	double value = nca[ij][i][j];
	if (fabs(value)>1e-5) tncab.push_back(make_pair(j,value));
      }
      //      sort(tncab.begin(),tncab.end(),CmpStates_);
      ncab[ij][i]=tncab;
      deque<pair<int,double> > tncaf;
      for (int j=0; j<nca[ij][i].size(); j++){
	double value = nca[ij][j][i];
	if (fabs(value)>1e-5) tncaf.push_back(make_pair(j,value));
      }
      //      sort(tncaf.begin(),tncaf.end(),CmpStates_);
      ncaf[ij][i]=tncaf;
    }
  }
  vector<deque<pair<int,double> > > nca_Sz(nca_Gz.size());
  for (int i=0; i<nca_Gz.size(); i++){
    deque<pair<int,double> > tnca;
    for (int j=0; j<nca_Gz[i].size(); j++){
      double value = nca_Gz[i][j];
      if (fabs(value)>1e-5) tnca.push_back(make_pair(j,value));
    }
    //      sort(tncab.begin(),tncab.end(),CmpStates_);
    nca_Sz[i]=tnca;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////// Printing of cix file  ///////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ofstream cix("out.cix");
  cix<<"# Input file for impurity solver. J="<<J_coulomb<<" c="<<c_spinorb<<" Ex=[";
  for (int i=0; i<_E_exact.size()-1; i++) cix<<_E_exact[i]<<",";
  cix<<_E_exact[_E_exact.size()-1]<<"] ";
  cix<<" Ep=[";
  for (int i=0; i<_E_approx.size()-1; i++) cix<<_E_approx[i]<<",";
  cix<<_E_approx[_E_approx.size()-1]<<"] ";
  
  cix<<" Eoca="<<E_oca<<" l="<<l<<" para="<<paramagnetic<<endl;
  cix<<cmp_one_electron.size()<<" ";
  for (int i=0; i<cmp_one_electron.size(); i++) cix<<cmp_one_electron[i].size()<<" ";

  int npseudo=0;
  for (int in=nstart; in<nstop; in++) npseudo += pseudo[in].size();
  
  //  cix<<pseudo[0].size()+pseudo[1].size()+pseudo[2].size()<<" ";
  //  if (nstop-nstart>3) cix<<pseudo[3].size()+pseudo[-1].size()<<" "; else cix<<0<<" ";
  cix<<npseudo<<" "<<0<<" ";
  
  cix<<" # Number of baths it's degeneracy an number of local valence and local core states"<<endl;
  cix<<right<<setw(4)<<"#"<<" ";
  for (int ij=0; ij<cmp_one_electron.size(); ij++) cix<<setw(11)<<"N"<<ij;
  cix<<setw(2)<<"Mtot "<<setw(5)<<"deg      Eatom-E_SO ";
  for (int ij=0; ij<cmp_one_electron.size(); ij++) cix<<"#b  ";
  for (int ij=0; ij<cmp_one_electron.size(); ij++) cix<<"#f  ";
  cix <<endl;
  
  for (int iin=0; iin<nstop-nstart; iin++){   // over all occupancies
    int in = nord[iin];
    for (int i=0; i<pseudo[in].size(); i++){  // over all pseudoparticles
      int ii = pseudo[in][i];
      cix<<right<<setw(4)<<ii<<" ";
      for (int ij=0; ij<cmp_one_electron.size(); ij++)	cix<<setw(11)<<vfpfm[in][ij][i]/Deg[in][i]<<" ";
      double Eatom_LS = Energ[in][i]-(qatom ? 0 : E_LS[in][i]); if (fabs(Eatom_LS)<1e-10) Eatom_LS=0;
      cix<<setw(2)<<Np[in]<<" "<<setw(5)<<Deg[in][i]<<" "<<setw(11)<<Eatom_LS<<" ";
      for (int ij=0; ij<cmp_one_electron.size(); ij++)	cix<<setw(3)<<ncab[ij][ii].size()<<" ";
      for (int ij=0; ij<cmp_one_electron.size(); ij++)	cix<<setw(3)<<ncaf[ij][ii].size()<<" ";
      for (int ij=0; ij<cmp_one_electron.size(); ij++)
	for (int j=0; j<ncab[ij][ii].size(); j++) cix<<right<<setw(11)<<ncab[ij][ii][j].second<<" x "<<left<<setw(3)<<ncab[ij][ii][j].first<<" ";
      for (int ij=0; ij<cmp_one_electron.size(); ij++)
	for (int j=0; j<ncaf[ij][ii].size(); j++) cix<<right<<setw(11)<<ncaf[ij][ii][j].second<<" x "<<left<<setw(3)<<ncaf[ij][ii][j].first<<" ";
      cix<<endl;
    }
  }

  list_oca.sort(cmpOCA);
  if (qOCA){
    cix<<"# OCA diagrams, information is (pp0,pp1,pp2,pp3) (b1,b2) fact , where pp is pseudoparticle and b is bath"<<endl;
    for (list<OCAd>::const_iterator p=list_oca.begin(); p!=list_oca.end(); p++){

      bool QkOCA = qkOCA.size()==0; // if kOCA is -1 in the command line, we do nothing.
      for (int iq=0; iq<qkOCA.size(); iq++)
	if (qkOCA[iq]==p->states[1] || qkOCA[iq]==p->states[3]) QkOCA = true;
      if (!QkOCA) continue;
		
      cix<<setw(3)<<p->states[0]<<" "<<setw(3)<<p->states[1]<<" "<<setw(3)<<p->states[2]<<" "<<setw(3)<<p->states[3]<<"    ";
      cix<<setw(3)<<p->ib1<<" "<<setw(3)<<p->ib2<<"    "<<setw(20)<<-p->f<<endl;
    }
  }
  if (susceptibility){
    // Printing susceptibility
    cix<<"# Diagrams for calculating magnetic susceptibility"<<endl;
    for (int iin=0; iin<nstop-nstart; iin++){   // over all occupancies
      int in = nord[iin];
      for (int i=0; i<pseudo[in].size(); i++){  // over all pseudoparticles
	int ii = pseudo[in][i];
	cix<<right<<setw(4)<<ii<<" ";
	
	double sum=0;
	for (int j=0; j<nca_Sz[ii].size(); j++) sum += nca_Sz[ii][j].second;
	sum *=12/Deg[in][i];
	cix<<setw(12)<<sqrt(sum)<<"   ";
	
	cix<<setw(3)<<nca_Sz[ii].size()<<" ";
	for (int j=0; j<nca_Sz[ii].size(); j++) cix<<right<<setw(11)<<nca_Sz[ii][j].second<<" x "<<left<<setw(3)<<nca_Sz[ii][j].first<<" ";
	cix<<endl;
      }
    }
  }

  return 0;
}
