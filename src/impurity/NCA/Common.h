// @Copyright 2007 Kristjan Haule
// 
#ifndef _COMMON_
#define _COMMON_
#include "zeroin.h"
#include "average.h"
#include <map>
#include <vector>
#ifdef _STRSTREAM
#include <strstream>
#endif

using namespace std;
typedef vector<int>::size_type vint;

int Binomial(int n, int m)
{
  int Mf = 1;
  for (int i=2; i<=m; i++) Mf *= i;
  int r = 1;
  for (int i=n; i>=n-m+1; i--) r*=i;
  return r/Mf;
}

//Common constants and variables
class common{
public:
  static double U;
  static double T;
  //  static double J;
  static int baths;
  static int Na, Nc;
  static function1D<int> Ns;
  static function2D<double> Ms;
  static function1D<int> Mtot;
  static function1D<int> deg;
  static function1D<double> sJc;
  static vector<vector<map<int,double> > > sncab;	   // index for hole diagrams
  static vector<vector<map<int,double> > > sncaf;	   // index for particle diagrams
  static vector<map<int,double> > suscb;                   // index for susceptibility
  static function2D<int> ncab;     // index for hole diagrams
  static function2D<int> ncaf;     // index for particle diagrams
  static function2D<double> prefactb; // prefactor for hole digrams
  static function2D<double> prefactf; // prefactor for particle diagrams
  static function2D<double> prefactG; // prefactor to calculate local Green's function
  static function1D<double> Ed;
  static function1D<double> Sinfty;
  static function1D<double> nalpha;
  static function1D<double> miss_nd;
  static function2D<double> moment;
  static double beta;
  static double delta;
  static double Q;
  static double Q0;
  static double nd;
  static double nd0;
  static double lambda0;
  static string outdir;
  static int totDeg;
  static function1D<string> Eds;
  static int N_ac;
  static double dom_ac;
  static int acore, pcore;
  static bool SubtractLorentz;
  static double LorentzMaxRatio;
  static double SearchLorentz;
  static int FirstLorentz;
  static int LastLorentz;
  static double dlmin;
  static bool renorm_core, renorm;
  static bool cmp_susc;
  static double Fimp, Epot, TrLogGimp;
  static void SetParameters(Par<double>& Ed_, double U_, /*double J_, */double T_, double Q0_, const string& outdir_, int N_ac_,
			    double dom_ac_, int acore_, int pcore_, bool SubtractLorentz_, double SearchLorentz_,
			    double LorentzMaxRatio_, int FirstLorentz_, int LastLorentz_,
			    bool renorm_core_, bool renorm_)
  {
    dlmin = 2.0;
    LorentzMaxRatio = LorentzMaxRatio_;
    SearchLorentz = SearchLorentz_;
    SubtractLorentz=SubtractLorentz_;
    FirstLorentz=FirstLorentz_; // First pseudoparticle which could be augmented with lorentz
    LastLorentz=LastLorentz_; // Last pseudoparticle which could be augmented with lorentz
    Ed.resize(baths);
    int i=0;
    while (Ed_.next() && i<baths) {
      Ed[i] = Ed_;
      i++;
    }
    for (int j=i; j<baths; j++) Ed[j]=Ed[i-1];
    T = T_; 
    U = U_;
    //    J = J_;
    beta=1/T_;
    Q0 = Q0_;
    outdir = outdir_;
    Eds.resize(baths);
    for (int i=0; i<baths; i++){
      stringstream t; t<<"E"<<i;
      Eds[i] = t.str();
    }
    nalpha.resize(baths);
    miss_nd.resize(baths);
    for (int i=0; i<baths; i++) miss_nd[i]=0;
    N_ac = N_ac_;
    dom_ac = dom_ac_;
    acore = acore_;
    pcore = pcore_;
    renorm_core=renorm_core_;
    renorm=renorm_;
    moment.resize(baths,2);
    Fimp=Epot=TrLogGimp=0.0;
  }
  static void ParsInputFile(const string& filename);
  static void PrintParsedData(ostream& stream);
  static ostream& printHead(ostream& stream);
};

class sLorentz{
public:
  double x0, gamma, P;
  bool exist;
  sLorentz() : x0(0), gamma(1), P(0), exist(false){};
  void Set(double zero, double eps, double a, double p, double q, double r)
  {
    exist = true;
    //double A = (sqr(1-p)+sqr(q))/a-2*eps*q*r/sqr(a)+sqr(eps*r/a)/a;
    double A = (sqr(1-p)+sqr(q))/a-2*eps*q*(r/a)/a+sqr(eps*r/a)/a;
    double B = eps*q/a-sqr(eps)/a*(r/a)/2;
    double C = sqr(eps)/a;
    double b2 = C/A-sqr(B/A);
    x0 = -B/A;
    gamma = (b2>0)? sqrt(b2) : sqrt(abs(C/A));
    if (gamma==0) {
      exist=false; P=0;
      return;
    }
    //cout<<"a="<<a<<" A="<<A<<" B="<<B<<" C="<<C<<" b2="<<b2<<" gamma="<<gamma<<endl;
    P = 1/(A*gamma);
    x0 += zero;
  }
  void SetFalse(){exist=false; P=0;}
private:
  double IntgA(double om0, double om1, double A0, double A1, double omega, double x0) const
  {
    if (!exist) return 0;
    if (fabs(om1-om0)*100<gamma) return P*gamma*0.5*(A0+A1)*(om1-om0)/(sqr(0.5*(om0+om1)+omega-x0)+sqr(gamma));
    double c0 = om0 + omega - x0;
    double c1 = om1 + omega - x0;
    double dA = (A1-A0)/(om1-om0);
    if (abs(c0)>100*gamma && abs(c1)>100*gamma && c0*c1>0) return P*gamma*( (A0-dA*c0)*(1/c0-1/c1)+dA*log(c1/c0)+0.5*dA*(sqr(gamma/c1)-sqr(gamma/c0)) ); ///// HERE WAS A BUG!! Corrected Dec/6/2013.
    if (abs(c0)>100*gamma && abs(c1)>100*gamma && c1-c0>199.9*gamma) return P*( (A0-dA*c0)*(M_PI+gamma*(1/c0-1/c1))+dA*gamma*log(abs(c1/c0))+0.5*dA*gamma*(sqr(gamma/c1)-sqr(gamma/c0)) ); ///// HERE WAS A BUG!! Corrected Dec/6/2013.
    //if (abs(c0)>1 && abs(c1)>1) return P*gamma*(c1-c0)*0.5*(A1+A0)/(c1*c0); ///// HERE WAS A BUG!! Corrected Dec/6/2013.
    double a0 = c0/gamma;
    double a1 = c1/gamma;
    double R;
    if (fabs(gamma)<1e-30){
      R= P*gamma*((A0-dA*c0)*(1/c0-1/c1)+dA*log(fabs(c1/c0)));
    }else{
      R = P*((A0-dA*c0)*(atan(a1)-atan(a0))+0.5*gamma*dA*log((1+sqr(a1))/(1+sqr(a0))));
    }
    if (std::isnan(R) || std::isinf(R)){
      cerr<<"R is nan or inf "<<R<<" "<<om0<<" "<<om1<<" "<<A0<<" "<<A1<<" "<<omega<<" "<<x0<<" "<<c0<<" "<<c1<<endl;
      cerr<<"to "<<(1+sqr(a1))<<" "<<(1+sqr(a0))<<" a0="<<a0<<" a1="<<a1<<" gamma="<<gamma<<" c0="<<c0<<" c1="<<c1<<" "<<atan(a1)-atan(a0)<<" "<<(A0-dA*c0)<<" "<<log((1+sqr(a1))/(1+sqr(a0)))<<endl;
    }
    return R;
  }
public:
  double IntgAp(double om0, double om1, double A0, double A1, double omega)const{
    return  IntgA(om0, om1, A0, A1, omega, x0);}
  double IntgAm(double om0, double om1, double A0, double A1, double omega)const{
    return  IntgA(om0, om1, A0, A1, -omega, -x0);}
  double IntgApLL(const sLorentz& l, double omega) const
  {
    return P*l.P*M_PI*(gamma+l.gamma)/(sqr(gamma+l.gamma)+sqr(x0-l.x0-omega));
  }
  double V(double x){ return P*gamma/(sqr(x-x0)+sqr(gamma));}
  friend ostream& operator<<(ostream& stream, const sLorentz& s);
};
ostream& operator<<(ostream& stream, const sLorentz& s)
{
  if (s.exist)
    stream<<setw(15)<<s.x0<<" "<<setw(15)<<s.gamma<<" "<<setw(15)<<s.P<<" ";
  return stream;
}

// Auxiliary self-energies and spectral functions
class Auxiliary{
  const int Na, Nc, baths;
  mesh1D om;
  function1D<double> fe;
  function1D<double> fedh;
  function1D<double> logo;
  function2D<double> Sigt;
  function2D<double> Sigtn;
  function2D<dcomplex> Sigc;
  function2D<dcomplex> Sigcore;
  function2D<double> Gt;
  function2D<double> Gp;
  function2D<double> Gm;
  vector<function2D<double> > aAc;
  function1D<double> Acx;
  function1D<double> Acy;
  function2D<double> Acp, Acm;
  AvFun<double> aF;
  function1D<double> Energy;
  function1D<double> Probability;
  mesh1D oml;
  function2D<dcomplex> Deltam_ac, Deltap_ac;
  function1D<double> mom_Deltam_ac, mom_Deltap_ac;
  function1D<dcomplex> Sigtmp;
  int mpos, m0, m1;
  function2D<double> GtA1, GtA2;
  vector<sLorentz> lorentzm, lorentzp;
public:
  Auxiliary (int Na_, int Nc_, int baths_) : Na(Na_), Nc(Nc_), baths(baths_), aAc(2*baths), mom_Deltam_ac(baths), mom_Deltap_ac(baths),
					     lorentzm(Na), lorentzp(Na),Probability(Na){};
  bool ReadSelfEnergy(const string& filename, const Par<double>& Ed, const Par<double>& T, const Par<double>& U, const mesh1D& ph_omd, const function2D<double>& ph_Ac);
  void KramarsKronig();
  double DeterminSpectralFunctions(double StartLambda, double EndLambda, double dLamdba, int followPeak);
  void PrintOutMeanQ(double StartLambda, double EndLambda);
  void PrintNorm(ostream& stream);
  void Print(int l, string dir);
  void Printn(int l);
  void SetSignToZero(){Sigtn=0.0;Sigcore=0.0;}
  
  void SetUpAverageAc(const mesh1D& omd, const mesh1D& momd, const function2D<double>& Ack, const function1D<double>& fed);
  void CalcSigmab(const mesh1D& omd);
  void CalcSigmaf(const mesh1D& omd);
  
  double Difference();
  double DeterminSelfEnergies(double alpha, int CmpDiff);
  const mesh1D& omega() const {return om;}
  double ferm(int i) const {return fe[i];}
  const function2D<double>& _Gp() const {return Gp;}
  const function2D<double>& _Gm() const {return Gm;}
  
  void PrintSign();
  
  double Q(double lambda);
  double operator()(double lambda);
  double minEnergy;
  void PrintCore(const string& filename);
  const function1D<double>& Energ() const{return Energy;}
  const vector<sLorentz>& Lorentzm()const{return lorentzm;}
  const vector<sLorentz>& Lorentzp()const{return lorentzp;}
  
  void CreateSigma000(const mesh1D& omd, const function2D<double>& Ac);
private:
  void Print_aAc(int l);
  void Print_Qux(int l);
  void Print_Sign(int l, int st, int en);
  void PrintOutMeanQ(int M, double StartLambda, double EndLambda);
};

// Physical electron spectral function and suscpetibility
// Physical observables
class Physical{
public:
  const int Na, Nc, baths;
  mesh1D omd;
  function2D<dcomplex> G00;
  function2D<double> A00;
  function1D<double> C00;
  function1D<dcomplex> Chi;
  function2D<double> A00c;
  function2D<dcomplex> Sig;
private:
  mesh1D momd;
  function1D<double> fed;
  function1D<double> logod;
  function1D<double> th;
  function2D<double> Ac;
  function2D<dcomplex> Delta0;
  vector<AvFun<double> > aF;
  function2D<double> Gtx;
  function2D<double> Cmp;
  function1D<double> tG;
  function1D<bool> Pexists;
public:
  Physical(int Na_, int Nc_, int baths_);
  bool ReadBathFunction(const string& filename, bool spectra);
  void CalculateA00(const mesh1D& omega, const function2D<double>& Gp, const function2D<double>& Gm,
		    const function1D<double>& Energy, const vector<sLorentz>& lorentzm, const vector<sLorentz>& lorentzp);
  void KramarsKronig();
  void DeterminG00(double alpha,ostream& loging);
  double Difference();
  void Print(int l, string dir);
  void Print0(const string& filename);
  
  const mesh1D& omega() const {return omd;}
  const mesh1D& momega() const {return momd;}
  const function1D<double>& fe() const {return fed;}
  const function2D<double>& Ac0() const {return Ac;}
  
  void PrintA00(ostream& out);
  void CalcSelfEnergy();
  void MissingDoping(double start);
private:
  void CalculateProducts(double u, double fu, const mesh1D& om, const function2D<double>& Gm);
  
  bool ReadBeginning(const string& filename, istream& input, int& n, int& m, bool& begincomment, double& center);
};

void AverageFunction(const mesh1D& omx, double u, const mesh1D& eps, AvFun<double>& aF, functionb<double>& aAc)
{
  apar ap;
  cintpar pi;
  tint position = omx.InitInterpLeft();
  InterpLeft(eps[0]-u, omx, position, pi);
  aF.InterpolateFirst(pi);
  InterpLeft(eps[1]-u, omx, position, pi);
  ap.SetUpCsFirst(u, eps);
  aAc[0] = aF.InterpolateNext(pi, ap) * eps.Dh(0);
  for (int j=1; j<eps.size()-1; j++){
    InterpLeft(eps[j+1]-u, omx, position, pi);
    ap.SetUpCs(u, j, eps, omx.Dh(pi.i));
    aAc[j] = aF.InterpolateNext(pi, ap) * eps.Dh(j);
  }
  ap.SetUpCsLast(u, eps);
  aAc[eps.size()-1] = aF.InterpolateLast(ap) * eps.Dh(eps.size()-1);
}

inline double product(const double* A, const double* G, int size)
{
  double sum = 0;
  for (int i=0; i<size; i++) sum += A[i]*G[i];
  return sum;
}

void Auxiliary::SetUpAverageAc(const mesh1D& omd, const mesh1D& momd, const function2D<double>& Ack, const function1D<double>& fed)
{
  int m = om.find_(0.0)+1;
  
  Acx.resize(omd.size());

  for (int b=0; b<baths; b++){
    aAc[b].resize(om.size(),om.size());

    for (int i=0; i<omd.size(); i++) Acx[i] = Ack[b][i]*(1-fed[i]);
    aF.SetUp(Acx,omd);
    for (int i=0; i<m; i++) AverageFunction(omd,om[i],om,aF,aAc[b][i]);

    for (int i=0; i<omd.size(); i++) Acx[i] = Ack[b][i]*fed[i];
    aF.SetUp(Acx,omd);
    for (int i=m; i<om.size(); i++) AverageFunction(omd,om[i],om,aF,aAc[b][i]);
  
    aAc[baths+b].resize(om.size(),om.size());
  
    for (int i=0; i<momd.size(); i++) Acx[momd.size()-i-1] = Ack[b][i]*fed[i];
    aF.SetUp(Acx,momd);
    for (int i=0; i<m; i++) AverageFunction(momd,om[i],om,aF,aAc[baths+b][i]); 

    for (int i=0; i<momd.size(); i++) Acx[momd.size()-i-1] = Ack[b][i]*(1-fed[i]);
    aF.SetUp(Acx,momd);
    for (int i=m; i<om.size(); i++) AverageFunction(momd,om[i],om,aF,aAc[baths+b][i]);
  }

  // For core part need Delta in more extended range
  Acy.resize(omd.size());
  oml.resize(omd.size()+2*common::N_ac);
  for (int i=0; i<common::N_ac; i++) oml[i] = omd[0]-(common::N_ac-i)*common::dom_ac;
  for (int i=0; i<omd.size(); i++)  oml[i+common::N_ac]            = omd[i];
  for (int i=0; i<common::N_ac; i++) oml[omd.size()+common::N_ac+i] = omd.last()+(i+1)*common::dom_ac;
  oml.SetUp(omd.dcenter());
  
  Deltam_ac.resize(baths,oml.size());
  Deltap_ac.resize(baths,oml.size());
  Acp.resize(baths,omd.size());
  Acm.resize(baths,omd.size());
  for (int b=0; b<baths; b++){
    for (int i=0; i<omd.size(); i++){
      Acm(b,i) = Ack[b][i]*fed[i];
      Acp(b,i) = Ack[b][i]*(1-fed[i]);
    }
    int ofst=0;
    #pragma omp parallel for
    for (int i=0; i<common::N_ac; i++){
      double Deltar = ::KramarsKronig(Acm[b], omd, oml[i], 0, 0.0);
      Deltam_ac[b][i] = dcomplex(-M_PI*Deltar,0.0);
      Deltar = ::KramarsKronig(Acp[b], omd, oml[i], 0, 0.0);
      Deltap_ac[b][i] = dcomplex(-M_PI*Deltar,0.0);
    }
    ofst=common::N_ac;
    #pragma omp parallel for
    for (int i=0; i<omd.size(); i++){
      double Deltar = ::KramarsKronig(Acm[b], omd, omd[i], i, Acm[b][i]);
      Deltam_ac[b][ofst+i] = dcomplex(-M_PI*Deltar,-M_PI*Acm[b][i]);
      Deltar = ::KramarsKronig(Acp[b], omd, omd[i], i, Acp[b][i]);
      Deltap_ac[b][ofst+i] = dcomplex(-M_PI*Deltar,-M_PI*Acp[b][i]);
    }
    ofst=common::N_ac+omd.size();
    #pragma omp parallel for
    for (int i=0; i<common::N_ac; i++){
      double Deltar = ::KramarsKronig(Acm[b], omd, oml[omd.size()+common::N_ac+i], omd.size()-1, 0.0);
      Deltam_ac[b][ofst+i] = dcomplex(-M_PI*Deltar, 0.0);
      Deltar = ::KramarsKronig(Acp[b], omd, oml[omd.size()+common::N_ac+i], omd.size()-1, 0.0);
      Deltap_ac[b][ofst+i] = dcomplex(-M_PI*Deltar, 0.0);
    }
    double summ=0;
    for (int i=0; i<omd.size(); i++) summ += Acm[b][i]*omd.Dh(i);
    double sump=0;
    for (int i=0; i<omd.size(); i++) sump += Acp[b][i]*omd.Dh(i);
    mom_Deltam_ac[b] = summ;
    mom_Deltap_ac[b] = sump;
  }
}

void Auxiliary::CalcSigmab(const mesh1D& omd)
{
  for (int b=0; b<baths; b++){
    GtA1.Product(Gm,aAc[b],0,mpos);               // Gm[f,eps]*Acfm[x,eps]
    GtA2.Product(Gp,aAc[b],mpos,aAc[b].size_N()); // Gp[f,eps]*Acfp[x,eps]

    if (common::SubtractLorentz){
      #pragma omp parallel for
      for (int j=0; j<Na; j++){
	if (lorentzm[j].exist){
	  tint pos0=omd.size()-2, pos1=omd.size()-2;
	  double dlmin_x0 = -common::dlmin + lorentzm[j].x0;
	  double dlmin_x1 =  common::dlmin + lorentzm[j].x0;
	  for (int i=0; i<mpos; i++){
	    int k0 = omd._find(dlmin_x0 - om[i], 0, pos0);
	    int k1 = omd._find(dlmin_x1 - om[i], 0, pos1);
	    double sum=0;
	    for (int k=k0; k<k1; k++)
	      sum += lorentzm[j].IntgAp(omd[k], omd[k+1], Acp(b,k), Acp(b,k+1), om[i]);
	    GtA1(j,i) += sum;
	  }
	}
	if (lorentzp[j].exist){
	  tint pos0=omd.size()-2, pos1=omd.size()-2;
	  double dlmin_x0 = -common::dlmin + lorentzp[j].x0;
	  double dlmin_x1 =  common::dlmin + lorentzp[j].x0;
	  for (int i=mpos; i<om.size(); i++){
	    int k0 = omd._find(dlmin_x0 - om[i], 0, pos0);
	    int k1 = omd._find(dlmin_x1 - om[i], 0, pos1);
	    double sum = 0;
	    for (int k=k0; k<k1; k++)
	      sum += lorentzp[j].IntgAp(omd[k], omd[k+1], Acm(b,k), Acm(b,k+1), om[i]);
	    GtA2(j,i-mpos) += sum;
	  }
	}
      }
    }
    
    #pragma omp parallel for
    for (int j=0; j<Na; j++){
      for (map<int,double>::const_iterator l=common::sncab[j][b].begin(); l!=common::sncab[j][b].end(); l++){
	int ind = l->first;
	if (ind>=0 && ind<Na){
	  double prf = l->second/static_cast<double>(common::deg[j]);
	  for (int i=0; i<mpos; i++)         Sigtn(j,i) += prf * GtA1(ind,i)/fe[i];
	  for (int i=mpos; i<om.size(); i++) Sigtn(j,i) += prf * GtA2(ind,i-mpos)/(1-fe[i]);
	}
      }
    }
  }
  if (!common::acore) return;
  
  for (int b=0; b<baths; b++){
    for (int j=0; j<Na; j++){
      for (map<int,double>::const_iterator l=common::sncab[j][b].begin(); l!=common::sncab[j][b].end(); l++){
	int ind = l->first;
	if (ind>=Na && ind<Na+Nc){
	  double prf = l->second/static_cast<double>(common::deg[j]);
	  tint position = oml.InitInterpRight();
	  for (int i=0; i<om.size(); i++){
	    double x = Energy[ind]-common::lambda0-om[i];
	    dcomplex Delta=0;
	    if (x>oml.last()) Delta = mom_Deltam_ac[b]/x;
	    else Delta = Deltam_ac[b](oml.InterpRight(x, position));
	    Sigcore[j][i] += prf*Delta;
	  }
	}
      }
    }
  }
}



void Auxiliary::CalcSigmaf(const mesh1D& omd)
{
  for (int b=0; b<baths; b++){
    GtA1.Product(Gm,aAc[baths+b],0,mpos);
    GtA2.Product(Gp,aAc[baths+b],mpos,aAc[baths+b].size_N());

    if (common::SubtractLorentz){
      #pragma omp parallel for
      for (int j=0; j<Na; j++){
	if (lorentzm[j].exist){
	  tint pos0=0, pos1=0;
	  double dlmin_x0 = -common::dlmin - lorentzm[j].x0;
	  double dlmin_x1 =  common::dlmin - lorentzm[j].x0;
	  for (int i=0; i<mpos; i++){
	    int k0 = omd.find_(dlmin_x0 + om[i], pos0);
	    int k1 = omd.find_(dlmin_x1 + om[i], pos1);
	    double sum = 0;
	    //for (int k=0; k<omd.size()-1; k++)
	    for (int k=k0; k<k1; k++)
	      sum += lorentzm[j].IntgAm(omd[k], omd[k+1], Acm(b,k), Acm(b,k+1), om[i]);
	    GtA1(j,i) += sum;
	  }
	}
	if (lorentzp[j].exist){
	  tint pos0=0, pos1=0;
	  double dlmin_x0 = -common::dlmin - lorentzp[j].x0;
	  double dlmin_x1 =  common::dlmin - lorentzp[j].x0;
	  for (int i=mpos; i<om.size(); i++){
	    int k0 = omd.find_(dlmin_x0 + om[i], pos0);
	    int k1 = omd.find_(dlmin_x1 + om[i], pos1);
	    double sum = 0;
	    //	    for (int k=0; k<omd.size()-1; k++)
	    for (int k=k0; k<k1; k++)
	      sum += lorentzp[j].IntgAm(omd[k], omd[k+1], Acp(b,k), Acp(b,k+1), om[i]);
	    GtA2(j,i-mpos) += sum;
	  }
	}
      }
    }

    #pragma omp parallel for
    for (int j=0; j<Na; j++){
      for (map<int,double>::const_iterator l=common::sncaf[j][b].begin(); l!=common::sncaf[j][b].end(); l++){
	int ind = l->first;
	if (ind>=0 && ind<Na){
	  double prf = l->second/static_cast<double>(common::deg[j]);
	  for (int i=0; i<mpos; i++)         Sigtn(j,i) += prf * GtA1(ind,i)/fe[i];
	  for (int i=mpos; i<om.size(); i++) Sigtn(j,i) += prf * GtA2(ind,i-mpos)/(1-fe[i]);
	}
      }
    }
  }
  if (!common::acore) return;
  
  for (int b=0; b<baths; b++){
    for (int j=0; j<Na; j++){
      for (map<int,double>::const_iterator l=common::sncaf[j][b].begin(); l!=common::sncaf[j][b].end(); l++){
	int ind = l->first;
	if (ind>=Na && ind<Na+Nc){
	  double prf = l->second/static_cast<double>(common::deg[j]);
	  tint position = oml.InitInterpLeft();
	  for (int i=0; i<om.size(); i++){
	    double x = om[i]-Energy[ind]+common::lambda0;
	    dcomplex Delta=0;
	    if (x<om[0]) Delta = mom_Deltap_ac[b]/x;
	    else Delta = Deltap_ac[b](oml.InterpLeft(x, position));
	    Sigcore[j][i] += prf*Delta;
	  }
	}
      }
    }
  }
}

void Auxiliary::CreateSigma000(const mesh1D& omd, const function2D<double>& Ac)
{// If inteligence guess for the pseudo-particles self-energy is not found,
 // it creates a guess using atomic type of approximation.
  Sigt=0;
  for (int b=0; b<baths; b++){
    for (int j=0; j<Na; j++){
      for (map<int,double>::const_iterator l=common::sncab[j][b].begin(); l!=common::sncab[j][b].end(); l++){
	int ind = l->first;
	if (ind>=0 && ind<Na){
	  double Em = Energy[ind]-minEnergy;
	  double prf = l->second/static_cast<double>(common::deg[j]);
	  tint pos = omd.InitInterpRight();
	  for (int i=0; i<om.size(); i++){
	    double ff;
	    if (om[i]>0) ff = ferm_f((Em-om[i])/common::T)/(1-fe[i]);
	    else{
	      double eom = exp(om[i]/common::T);
	      ff = (eom+1.)/(eom+exp(Em/common::T));
	    }
	    Sigt(j,i) += -M_PI*prf*ff*Ac[b](omd.InterpRight(Em-om[i],pos));
	  }
	}
      }
      for (map<int,double>::const_iterator l=common::sncaf[j][b].begin(); l!=common::sncaf[j][b].end(); l++){
	int ind = l->first;
	if (ind>=0 && ind<Na){
	  double Em = Energy[ind]-minEnergy;
	  double prf = l->second/static_cast<double>(common::deg[j]);
	  tint pos = omd.InitInterpLeft();
	  for (int i=0; i<om.size(); i++){
	    double ff;
	    if (om[i]>0) ff = ferm_f((Em-om[i])/common::T)/(1-fe[i]);
	    else{
	      double eom = exp(om[i]/common::T);
	      ff = (eom+1.)/(eom+exp(Em/common::T));
	    }
	    Sigt(j,i) += -M_PI*prf*ff*Ac[b](omd.InterpLeft(om[i]-Em,pos));
	  }
	}
      }
    }
  }
  KramarsKronig();
}

inline ostream& common::printHead(ostream& stream)
{
  stream<<"# ";
  stream<<" nb="<<baths<<" ";

  //stream<<" T="<<T<<" ntot="<<nd<<" U="<<U<<" lambda0="<<lambda0<<"  ";
  stream<<" T="<<T<<" ntot="<<nd<<" U="<<U<<" dFimpG="<<Fimp-TrLogGimp<<" Fimp="<<Fimp<<" Epot="<<Epot<<" TrLogGimp="<<TrLogGimp<<" lambda0="<<lambda0<<"  ";

  stream<<" Ns=[";
  for (int i=0; i<baths-1; i++) stream<<Ns[i]<<",";
  stream<<Ns[baths-1]<<"]  ";
  
  stream<<" Eimp=[";
  for (int i=0; i<baths-1; i++) stream<<Ed[i]<<",";
  stream<<Ed[baths-1]<<"]  ";

  stream<<" nf=[";
  for (int i=0; i<baths-1; i++) stream<<nalpha[i]<<",";
  stream<<nalpha[baths-1]<<"]  ";
  
  stream<<" md=[";
  for (int i=0; i<baths-1; i++) stream<<miss_nd[i]<<",";
  stream<<miss_nd[baths-1]<<"]  ";

  stream<<" moment=[";
  for (int i=0; i<baths-1; i++) stream<<"["<<moment[i][0]<<","<<moment[i][1]<<"],";
  stream<<"["<<moment[baths-1][0]<<","<<moment[baths-1][1]<<"]]  ";
  
  if (Sinfty.size()>0){
    double aS=0; for (int i=0; i<baths; i++) aS += Sinfty[i]; aS/=baths;
    stream<<" aSinfty="<<aS<<" ";
    stream<<" Sinfty=(";
    for (int i=0; i<baths-1; i++)stream<<Sinfty[i]<<",";
    stream<<Sinfty[baths-1]<<") ";
  }

  return stream;
}

void RememberParams (int argc, char *argv[]){
  ofstream param ((common::outdir+"/history.nca").c_str(), ios::app);
  if (!param) cerr<<" Didn't suceeded to open params file!"<<(common::outdir+"/history.nca")<<endl;
  for (int i=0; i<argc; i++) param << argv[i] << " ";
  param << endl;
}

template <class T>
bool ReadValue(T& a, const std::string& variable, const std::string& str){
  std::string::size_type pos = str.find(variable);
  if (pos < std::string::npos){
    std::string::size_type poseq = str.find("=",pos);
    if (poseq<std::string::npos){
      std::istringstream streambuff(std::string(str,poseq+1));
      streambuff >> a;
    }
    return true;
  }
  return false;
}

bool Auxiliary::ReadSelfEnergy(const string& filename, const Par<double>& Ed, const Par<double>& T, const Par<double>& U, const mesh1D& ph_omd, const function2D<double>& ph_Ac){
  ifstream inputf(filename.c_str());
  istream input(inputf.rdbuf());
  input.seekg(0,ios::beg);
  
  if (!input) {
    cerr << "Can't open input file: " << filename << endl;
    return false;
  }
  // Is the input file started with comment?
  bool begincomment = false;
  int n = 0;
  string str;
  const double SpecNumber = -100000; 
  double T_ = SpecNumber, U_ = SpecNumber;
  function1D<double> Ed_(baths);
  Ed_ = SpecNumber;
  double center = 0;
  getline(input,str);
  if (str.find('#')<string::npos){
    begincomment = true;
    for (int i=0; i<baths; i++) ReadValue(Ed_[i], common::Eds[i], str);
    ReadValue(T_, "T", str);
    ReadValue(U_, "U", str);
    if (!ReadValue(center, "peakposition", str)) center=0;
  } else n++;
  
  if (!Ed.IsSet() && Ed_[0]!=SpecNumber) for (int i=0; i<baths; i++) common::Ed[i] = Ed_[i];
  if (!T.IsSet()  && T_!=SpecNumber) common::T = T_;
  if (!U.IsSet() && U_!=SpecNumber) common::U = U_;
  common::beta = 1./common::T;

  Energy.resize(Na+Nc);
  minEnergy=0;
  // Calculates auxiliary Energies
  for (int i=0; i<Na+Nc; i++){
    Energy[i] = 0;
    for (int j=0; j<baths; j++)  Energy[i] += common::Ed[j]*common::Ms[i][j];
    //    Energy[i] += 0.5*common::Mtot[i]*(common::Mtot[i]-1)*(common::U-0.5*common::J);
    //    Energy[i] += common::J*common::sJc[i];
    Energy[i] += 0.5*common::Mtot[i]*(common::Mtot[i]-1)*common::U;
    Energy[i] += common::sJc[i];
    if (Energy[i]<minEnergy) minEnergy = Energy[i];
  }

  clog<<"************* Parameters ****************"<<endl;
  clog<<"  		U  = "<<common::U<<endl;
  for (int i=0; i<baths; i++)
    clog<<"  		Ed"<<i<<" = "<<common::Ed[i]<<endl;
  clog<<"  		T  = "<<common::T<<endl;
  for (int i=0; i<baths; i++)
    clog<<"  		N"<<i<<" = "<<common::Ns[i]<<endl;
  for (int i=0; i<Na+Nc; i++){
    if (i<Na) clog<<" valence state"<<setw(2)<<left<<i<<right<<" = ";
    else clog<<" core    state"<<i<<" = ";
    for (int j=0; j<baths; j++) clog<<setw(2)<<common::Ms[i][j];
    clog<<" with Energy"<<setw(2)<<left<<i<<right<<" = "<<Energy[i]<<endl;
  }
  clog<<"*****************************************"<<endl;
  
  // Computes the number of columns in file
  if (!input) {
    cerr << "ERROR: Wrong file format for Sigm" << endl;
    return false;
  }
  getline(input,str);  n++;
#ifdef _STRSTREAM
  strstream oneline;
  oneline << str <<ends;
#else
  istringstream oneline(str);
#endif
  int m=0; double t;
  while (oneline){oneline>>t; m++;}
  m--;
  while (input){ getline(input,str); n++;}
  n--;
  
  clog << filename << ": Number of entries: "<< n <<endl;
  clog << filename << ": Number of columns: "<< m <<endl;
  clog << filename << ": Peak-position "<< center <<endl;

  bool CreateDefault = false;
  if (m<2*Na+1){
    //cerr<<"ERROR: Not enough columns is input Sigma file. Exiting!"<<endl;
    clog<<"WARRNING: Not enough columns is input self-energy for pseudoparticles.... Creating default!"<<endl;
    CreateDefault = true;
  }
  inputf.seekg(0,ios::beg);
  //  clog<<"Premaknil na "<< inputf.tellg()<<endl;
  if (begincomment) inputf.ignore(10000,'\n');
  if (!inputf){ cerr<<"Reopening didn't suceeded!"<<endl; return false;}
  
  om.resize(n);
  Sigt.resize(Na,n);
  Sigc.resize(Na,n);
  
  int l=0;
  double omega;
  while (inputf>>omega && l<n){
    om[l] = omega;
    if (!CreateDefault){
      for (int i=0; i<Na; i++){
	double Sr, St;
	inputf>>Sr;
	inputf>>St;
	Sigc(i,l) = dcomplex(Sr,-St);
	Sigt(i,l) = -St;
      }
    }
    getline(inputf, str);
    l++;
  }
  inputf.close();
  if (l<n) cerr<<"Something wrong by reading file "<<filename<<endl;
  om.SetUp(center);
  mpos = om.find_(0.0)+1;
  m0 = om.find_(-common::SearchLorentz);
  m1 = om.find_(common::SearchLorentz)+1;

  GtA1.resize(Na,mpos);
  GtA2.resize(Na,om.size()-mpos);
  Sigcore.resize(Na,om.size());
  Sigtn.resize(Na,om.size());
  Gt.resize(Na,om.size());
  Gp.resize(Na,om.size());
  Gm.resize(Na,om.size());
  fe.CalcFermOnMesh(common::beta, om);
  logo.CalcLogOnMesh(om);
  fedh.resize(om.size());
  for (int i=0; i<om.size(); i++) fedh[i] = fe[i]*om.Dh(i);
  

  if (CreateDefault){
    CreateSigma000(ph_omd, ph_Ac);
  }else{
    for (int j=0; j<Na; j++){
      for (int i=0; i<om.size(); i++)
	Sigc(j,i) = dcomplex(Sigc(j,i).real(), Sigc(j,i).imag()*(1-fe[i]));
    }
  }
  return true;
}

void Auxiliary::KramarsKronig()
{
  for (int l=0; l<Na; l++){
    for (int i=0; i<om.size(); i++) Sigc(l,i).imag() = Sigt(l,i)*(1-fe[i]);
    Sigc[l].KramarsKronig(om, logo);
  }
}

double Lambda(double E, const functionb<dcomplex>& Sigc, const functionb<double>& Sigx, const mesh1D& om)
{
  // looking for lambda such that \widetilde{G} has maximum at zero frequency.
  // Sufficient condition is that the derivative of 1/\widetilde{G} is zero at zero frequency.
  // One gets a quadratic equation for lambda and thus two roots. Then one chooses the root that maximizes \widetilde{G}.
  // If no root exists, than we take lambda that minimizes linear coeficient in the expansion of 1/\widetilde{G}.
  // The latter equation is linear and one always gets unique solution.
  intpar p = om.Interp(0.0); int i=p.i;
  dcomplex cs = -E-Sigc(p);
  dcomplex ds = (Sigc[i+1]-Sigc[i])*om.Delta(i);
  double cr = cs.real();
  double ci = cs.imag();
  double dcr = 1-ds.real();
  double dci = -ds.imag();
  double dSigx   = (Sigx[i+1]-Sigx[i])*om.Delta(i);
  double x = Sigx[i]/dSigx;
  double determinant2 = x*(x*dcr*dcr+2*ci*dci)-ci*ci;
  // Minimum can not be at zero. Try to find lambda that minimizes the linear coefficient in the expansion of 1/G
  // If 1/G = a + b omega + c omega^2 +... and the below determinant is smaller than zero, coefficient b can not be
  // set to zero. Than return lambda that gives the smallest b.
  if (determinant2<=0) return dcr*x-cr;
  double d2 = -sqrt(determinant2);
  double d1 = -cr + dcr*x;
  double v1 = 1/(sqr(ci)+sqr(cr+d1+d2));
  double v2 = 1/(sqr(ci)+sqr(cr+d1-d2));
  cout<<"Lambda="<<d1+d2<<" "<<d1-d2<<" "<<v1<<" "<<v2<<endl;
  if (fabs(v1)>fabs(v2)) return d1+d2;
  else return d1-d2;
}

double Auxiliary::Q(double lambda)
{
  double sumQ=0;
  for (int j=0; j<Na; j++){
    double mune = -Energy[j]+lambda;

    sLorentz lorentz;
    if (common::SubtractLorentz && j>=common::FirstLorentz && j<=common::LastLorentz){
      double v0 = om[m0]+mune-Sigc(j,m0).real(), v=v0;
      int ii=0;
      for (ii=m0+1; ii<m1; ii++) {
	v = om[ii]+mune-Sigc(j,ii).real();
	if (sign(v)*sign(v0)<0) break;
      }
      double denom = om[ii]-om[ii-1]-Sigc(j,ii).real()+Sigc(j,ii-1).real();
      if (denom==0) cout<<"denom="<<denom<<endl;
      if (sign(v)*sign(v0)<0 && denom!=0){
	double zero = om[ii-1]-(om[ii]-om[ii-1])*(om[ii-1]+mune-Sigc(j,ii-1).real())/(om[ii]-om[ii-1]-Sigc(j,ii).real()+Sigc(j,ii-1).real());
	intpar ip(ii-1,(zero-om[ii-1])/(om[ii]-om[ii-1]));
	double dom = om[ii]-om[ii-1];
	dcomplex  Sc = Sigc[j](ip);
	double ratio = abs(Sc.imag()/dom);
	if (ratio<common::LorentzMaxRatio){
	  double    Sm = Sigt[j](ip)*fe(ip);
	  dcomplex dSc = (Sigc[j][ii]-Sigc[j][ii-1])/dom; //(om[ii]-om[ii-1]);
	  double   dSm = (Sigt[j][ii]*fe[ii]-Sigt[j][ii-1]*fe[ii-1])/dom; //(om[ii]-om[ii-1]);
	  double   Sc_im = Sc.imag();
	  if (fabs(Sc_im)<1e-20) Sc_im=-1e-20;
	  if (fabs(Sm)<1e-20) Sm=-1e-20;
	  if (fabs(Sc_im)>=1e-20 && fabs(Sm)>=1e-20){ 
	    lorentz.Set(zero, Sc_im, Sm, dSc.real(), dSc.imag(), dSm);
	    //cout<<"QFound zero "<<setw(2)<<left<<j<<right<<setw(10)<<zero<<" "<<lorentz<<endl;//setw(15)<<Sc<<" "<<setw(15)<<-St<<" ";
	  }
	}
      }
    }
    double sum=0, v;
    for (int i=0; i<om.size(); i++){
      v = fedh[i]*Sigt(j,i)/(sqr(om[i]+mune-Sigc(j,i).real())+sqr(Sigc(j,i).imag()));
      if (lorentz.exist) v -= om.Dh(i)*lorentz.V(om[i]);
      sum -= v;
    }
    sum -= lorentz.P*M_PI;
    sumQ += sum*common::deg[j];
  }
  return (sumQ/M_PI);
}

inline double Auxiliary::operator()(double lambda)
{
  double Q_ = Q(lambda);
  return Q_-common::Q0;
}

void Auxiliary::PrintOutMeanQ(double StartLambda, double EndLambda)
{
  double a0 = StartLambda;
  int M = 100;
  double da0 = (EndLambda-StartLambda)/M;
  cout.precision(16);
  for (int i=0; i<M; i++){
    cout << a0 << setw(25) << operator()(a0) << endl;
    a0 += da0;
  }
}

double Auxiliary::DeterminSpectralFunctions(double StartLambda, double EndLambda, double dLambda, int followPeak)
{
  double lambda0;
  if (followPeak>=0 && followPeak<Na)
    lambda0 = Lambda(Energy[followPeak], Sigc[followPeak], Sigt[followPeak], om);
  else if (followPeak==-2){
    lambda0 = minEnergy;
  }else{
    double a0 = StartLambda, b0 = 0;
    int sign=0, nn=0;
    while (!sign && nn++<100){
      double pQ = operator()(a0);
      while (!sign && a0<=b0) {
	double sQ = operator()(a0+dLambda);
	sign = pQ*sQ<0;
	pQ = sQ;
	if (!sign) a0 += dLambda;
      }
      if (!sign) dLambda /= 2.0;
    }
    
    if (nn>=100) {
      cerr << "Can't find root for <Q>" << endl;
      PrintOutMeanQ(StartLambda, EndLambda);
      exit(1);
    }
    
    // loking for zero (lambda0)
    lambda0 = zeroin(a0, a0+dLambda, *this, 1e-15*common::Q0);
  }

  common::lambda0 = lambda0;
  clog << setprecision(16) << "; lambda = "<<lambda0<<"  "<<lambda0-minEnergy<<endl;

  double sumQ = 0, sumnd=0;
  function1D<double> dQ(Na);
  for (int j=0; j<Na; j++){
    double mune = -Energy[j]+lambda0;

    if (common::SubtractLorentz && j>=common::FirstLorentz && j<=common::LastLorentz){
      double v = om[m0]+mune-Sigc(j,m0).real(), v0=v;
      int ii=0;
      for (ii=m0+1; ii<m1; ii++) {
	v = om[ii]+mune-Sigc(j,ii).real();
	if (sign(v)*sign(v0)<0) break;
      }
      bool found = false;
      double denom = om[ii]-om[ii-1]-Sigc(j,ii).real()+Sigc(j,ii-1).real();
      if (sign(v)*sign(v0)<0 && denom!=0){
	double zero = om[ii-1]-(om[ii]-om[ii-1])*(om[ii-1]+mune-Sigc(j,ii-1).real())/(om[ii]-om[ii-1]-Sigc(j,ii).real()+Sigc(j,ii-1).real());
	intpar ip(ii-1,(zero-om[ii-1])/(om[ii]-om[ii-1]));
	double dom = om[ii]-om[ii-1];
	dcomplex  Sc = Sigc[j](ip);
	double ratio = abs(Sc.imag()/dom);
	//clog<<"ps"<<j<<" ratio="<<ratio<<endl;
	if (ratio<common::LorentzMaxRatio){
	  double    Sm = Sigt[j](ip)*fe(ip);
	  dcomplex dSc = (Sigc[j][ii]-Sigc[j][ii-1])/dom;
	  double   dSm = (Sigt[j][ii]*fe[ii]-Sigt[j][ii-1]*fe[ii-1])/dom;
	  double   Sc_im = Sc.imag();
	  if (fabs(Sc_im)<1e-20) Sc_im=-1e-20;
	  if (fabs(Sm)<1e-20) Sm=-1e-20;
	  if (fabs(Sc_im)>=1e-20 && fabs(Sm)>=1e-20){
	    found = true;
	    lorentzm[j].Set(zero, Sc_im, Sm,        dSc.real(), dSc.imag(), dSm);
	    lorentzp[j].Set(zero, Sc_im, Sc_im, dSc.real(), dSc.imag(), dSc.imag());
	    //cout<<"Sc.im="<<Sc.imag()<<" Sm="<<Sm<<" dSc.r="<<dSc.real()<<" dSc.i="<<dSc.imag()<<" dSm="<<dSm<<endl;
	    //cout<<"zero="<<zero<<" ratio="<<ratio<<" Sm="<<Sm<<" dSm="<<dSm<<" Sc_im="<<Sc_im<<endl;
	    cout<<"Found lorentz at "<<setw(4)<<left<<j<<right<<setw(10)<<zero<<" lm="<<lorentzm[j]<<" lp="<<lorentzp[j]<<" r-"<<setw(15)<<ratio<<endl;
	  }
	}
      }
      if (!found){
	lorentzp[j].SetFalse();
	lorentzm[j].SetFalse();
      }
    }
  }
//   // We want to make sure that only one integer occupacition is treated with lorentz
//   // because we did not yet implement Lorentz*Lorentz
//   int MaxMtot=0;
//   for (int i=0; i<Na; i++) if (MaxMtot<common::Mtot[i]) MaxMtot = common::Mtot[i];
//   function1D<int> lorex(MaxMtot+1);lorex=0;
//   for (int j=0; j<Na; j++) if (lorentzm[j].exist ||lorentzp[j].exist) lorex[common::Mtot[j]]++;
//   int imaxLorentz=0;
//   for (int i=0; i<=MaxMtot; i++) if (lorex[i]>lorex[imaxLorentz]) imaxLorentz=i;
//   for (int i=0; i<Na; i++){
//     if (lorentzm[i].exist && common::Mtot[i]!=imaxLorentz) { cout<<"Lorentzm for "<<i<<" not accepted!"<<endl; lorentzm[i].SetFalse();}
//     if (lorentzp[i].exist && common::Mtot[i]!=imaxLorentz) { cout<<"Lorentzp for "<<i<<" not accepted!"<<endl; lorentzp[i].SetFalse();}
//   }

  for (int j=0; j<Na; j++){
    double mune = -Energy[j]+lambda0;
    dQ[j]=0;
    for (int i=0; i<om.size(); i++){
      Gt(j,i) = Sigt(j,i)/(sqr(om[i]+mune-Sigc(j,i).real())+sqr(Sigc(j,i).imag()));
      Gm(j,i) = fe[i]*Gt(j,i);
      Gp(j,i) = (1-fe[i])*Gt(j,i);
      if (lorentzm[j].exist) Gm(j,i) -= lorentzm[j].V(om[i]);
      if (lorentzp[j].exist) Gp(j,i) -= lorentzp[j].V(om[i]);
      dQ[j] -= Gm(j,i)*om.Dh(i);
    }
    dQ[j] -= lorentzm[j].P*M_PI;
    dQ[j] *= common::deg[j]/M_PI;
    sumQ += dQ[j];
    sumnd += dQ[j]*common::Mtot[j];
  }
  clog<<"       Q = "<<sumQ<<endl;
  for (int j=0; j<Na; j++){
    Probability[j] = dQ[j]/sumQ;
    clog<<setprecision(16)<<"       n"<<j<<"="<<dQ[j]/sumQ<<endl;
  }

  for (int b=0; b<baths; b++){
    common::nalpha[b]=0;
    for (int j=0; j<Na; j++) common::nalpha[b] += dQ[j]*common::Ms[j][b];
    common::nalpha[b]/=sumQ;
  }
  common::Q = sumQ;
  
  common::Fimp = common::lambda0-common::T * ::log(common::Q);
  double Epot=0;
  for (int j=0; j<Na; j++) Epot += Probability[j]*Energy[j];
  double dEpot=0;
  for (int b=0; b<baths; b++) dEpot += common::Ed[b]*common::nalpha[b];
  common::Epot = Epot-dEpot;
  clog<<"        Fimp="<<common::Fimp<<" Epot="<<common::Epot<<" Epot+OneP="<<Epot<<endl;
  
  //  if (fabs(sumQ-common::Q0)>1e-10) cerr<<"Something wrong with Q "<<sumQ<<"!"<<endl;
  clog<<"        Q is here equal to "<<sumQ<<endl;
  return sumnd/sumQ;
}

void Auxiliary::Print(int l, string dir="")
{
  string filename;
  if (l<0) filename = common::outdir+"/Sigma"+dir;
  else filename = NameOfFile(common::outdir+"/Sigma", l);
  ofstream out1(filename.c_str());  out1.precision(16);
  common::printHead(out1)<<" peakposition="<<om.dcenter()<<endl;
  for (int i=0; i<om.size(); i++){
    out1<<setw(25)<<om[i];
    for (int j=0; j<Na; j++) out1<<setw(25)<<Sigc(j,i).real()<<" "<<setw(25)<<-Sigt(j,i);
    out1<<endl;
  }
  if (l<0) filename = common::outdir+"/Spec"+dir;
  else filename = NameOfFile(common::outdir+dir+"/Spec", l);
  ofstream out2(filename.c_str());  out2.precision(16);
  common::printHead(out2)<<" peakposition="<<om.dcenter()<<endl;
  for (int i=0; i<om.size(); i++){
    out2<<setw(25)<<om[i];
    for (int j=0; j<Na; j++) out2<<setw(25)<<-Gt(j,i);
    for (int j=0; j<Na; j++) out2<<setw(25)<<-Gp(j,i);
    for (int j=0; j<Na; j++) out2<<setw(25)<<-Gm(j,i);
    out2<<endl;
  }
}

void Auxiliary::Printn(int l)
{
  string filename;
  filename = NameOfFile(common::outdir+"/nSigma", l);
  ofstream out1(filename.c_str());  out1.precision(16);
  common::printHead(out1)<<" peakposition="<<om.dcenter()<<endl;
  for (int i=0; i<om.size(); i++){
    out1<<setw(25)<<om[i];
    for (int j=0; j<Na; j++) out1<<setw(25)<<-Sigtn(j,i);
    out1<<endl;
  }
}

Physical::Physical(int Na_, int Nc_, int baths_) : Na(Na_), Nc(Nc_), baths(baths_), aF(Na)
{
  Pexists.resize(Na);
  for (int j=0; j<Na; j++){
    Pexists[j]=false;
    for (int b=0; b<baths; b++){
      for (map<int,double>::const_iterator l=common::sncab[j][b].begin(); l!=common::sncab[j][b].end(); l++){
	if (l->first >=0 && l->first < Na){
	  Pexists[j]=true;
	  break;
	}
      }
    }
    if (!Pexists[j] && common::cmp_susc){
      for (map<int,double>::const_iterator l=common::suscb[j].begin(); l!=common::suscb[j].end(); l++)
	if (l->first >=0 && l->first < Na){
	  Pexists[j]=true;
	  break;
	}
    }
  }
}

bool Physical::ReadBeginning(const string& filename, istream& input, int& n, int& m, bool& begincomment, double& center)
{
  if (!input) {
    cerr << "Can't open input file: " << filename << endl;
    return false;
  }
  // Is the input file started with comment?
  begincomment = false;
  n = 0;
  string str;
  getline(input,str);
  if (str.find('#')<string::npos){
    begincomment = true;
    if (!ReadValue(center, "peakposition", str)) center=0;
  } else n++;
  // Computes the number of columns in file
  if (!input) {
    cerr << "ERROR: Wrong file format for Sigm" << endl;
    return false;
  }
  getline(input,str);  n++;
  stringstream oneline;
  oneline << str << ends;
  m=0; double t;
  while (oneline){oneline>>t; m++;}
  m--;
  while (input){ getline(input,str); n++;}
  n--;

  clog << filename << ": Number of entries: "<< n <<endl;
  clog << filename << ": Number of columns: "<< m <<endl;
  clog << filename << ": Peak-position "<< center <<endl;
  
  input.seekg(0, ios::beg);
  input.clear();
  if (begincomment) getline(input, str);
  
  return true;
}

bool Physical::ReadBathFunction(const string& filename, bool spectra=true) // spectra=true: only spectral function will be read not the retarded quantity
{
  ifstream inputf(filename.c_str());
  istream input(inputf.rdbuf());
  input.seekg(0,ios::beg);

  if (!input) {
    cerr << "Can't open input file: " << filename << endl;
    return false;
  }
  // Is the input file started with comment?
  bool begincomment = false;
  int n = 0;
  string str;
  double center=0;
  getline(input,str);
  if (str.find('#')<string::npos){
    begincomment = true;
    if (!ReadValue(center, "peakposition", str)) center=0;
  } else n++;
  // Computes the number of columns in file
  if (!input) {
    cerr << "ERROR: Wrong file format for " << filename << endl;
    return false;
  }
  getline(input,str);  n++;
#ifdef _STRSTREAM
  strstream oneline;
  oneline << str <<ends;
#else
  istringstream oneline(str);
#endif
  int m=0; double t;
  while (oneline){oneline>>t; m++;}
  m--;
  while (input){ getline(input,str); n++;}
  n--;

  clog << filename << ": Number of entries: "<< n <<endl;
  clog << filename << ": Number of columns: "<< m <<endl;
  clog << filename << ": Peak-position "<< center <<endl;

  int number_cols = baths+1;
  if (!spectra) number_cols = 2*baths+1;
  
  if (m<number_cols){
    cerr<<"ERROR: Not enough columns in bath input file! Exiting..."<<endl;
    return false;
  }
  inputf.seekg(0, ios::beg);
  clog<<"Premaknil na "<< inputf.tellg()<<endl;
  if (begincomment) inputf.ignore(1000,'\n');
  if (!inputf){ cerr<<"Reopening didn't suceeded!"<<endl; return false;}
  
  omd.resize(n);
  momd.resize(n);
  G00.resize(baths,n);
  A00.resize(baths,n);
  A00c.resize(baths,n);
  Sig.resize(baths,n);
  Ac.resize(baths,n);
  Delta0.resize(baths,n);
  if (common::cmp_susc){
    C00.resize(n);
    Chi.resize(n);
  }
  int l=0;
  double omega;
  while (inputf>>omega && l<n){
    omd[l] = omega;
    if (spectra)
      for (int j=0; j<baths; j++) inputf>>Ac(j,l);
    else{
      for (int j=0; j<baths; j++) {
	double dr, di;
	inputf>>dr; inputf>>di;
	Ac(j,l) = -di/M_PI;
	Delta0(j,l) = dcomplex(dr,di);
      }
    }
    getline(inputf, str);
    momd[n-l-1] = -omd[l];
    l++;
  }
  inputf.close();
  if (l<n) cerr<<"Something wrong by reading file "<<filename<<endl;
  omd.SetUp(center);
  momd.SetUp(-center);

  fed.CalcFermOnMesh(common::beta, omd);
  th.CalcTanhOnMesh(common::beta, omd);
  logod.CalcLogOnMesh(omd);

  if (spectra){
    for (int b=0; b<baths; b++){
      for (int i=0; i<omd.size(); i++){
	double Deltar = ::KramarsKronig(Ac[b], omd, omd[i], i, Ac[b][i]);
	Delta0(b,i) = dcomplex(-M_PI*Deltar,-M_PI*Ac[b][i]);
      }
    }
  }
  return true;
}

void Physical::CalculateProducts(double u, double fu, const mesh1D& om, const function2D<double>& Gm)
{
  apar ap;
  cintpar pi;
  tint position = om.InitInterpLeft();
  InterpLeft(om[0]-u, om, position, pi);
  
  #pragma omp parallel for
  for (int i=0; i<Na; i++) if (Pexists[i]) aF[i].InterpolateFirst(pi);
  
  InterpLeft(om[1]-u, om, position, pi);
  ap.SetUpCsFirst(u, om);
  
  #pragma omp parallel for
  for (int i=0; i<Na; i++) if (Pexists[i]) Gtx(i,0) = aF[i].InterpolateNext(pi, ap) * om.Dh(0);
  
  for (int j=1; j<om.size()-1; j++){
    InterpLeft(om[j+1]-u, om, position, pi);
    ap.SetUpCs(u, j, om, om.Dh(pi.i+1));

    //#pragma omp parallel for
    for (int i=0; i<Na; i++) if (Pexists[i]) Gtx(i,j) = aF[i].InterpolateNext(pi, ap) * om.Dh(j);
  }
  ap.SetUpCsLast(u, om);
  
  #pragma omp parallel for
  for (int i=0; i<Na; i++) if (Pexists[i]) Gtx(i,om.size()-1) = aF[i].InterpolateLast(ap) * om.Dh(om.size()-1);

  Cmp.resize(Na,Na);

  #pragma omp parallel for
  for (int i=0; i<Na; i++){
    for (int b=0; b<baths; b++){
      for (map<int,double>::const_iterator l=common::sncab[i][b].begin(); l!=common::sncab[i][b].end(); l++){
	int ind = l->first;
	if (ind>=0 && ind<Na) Cmp(i,ind) = product(Gtx[i].MemPt(),Gm[ind].MemPt(),om.size())/fu;
      }
    }
    if (common::cmp_susc){
      for (map<int,double>::const_iterator l=common::suscb[i].begin(); l!=common::suscb[i].end(); l++){
	int ind = l->first;
	if (ind>=0 && ind<Na) Cmp(i,ind) = product(Gtx[i].MemPt(),Gm[ind].MemPt(),om.size())/fu;
      }
    }
  }
}

void Physical::CalculateA00(const mesh1D& omega, const function2D<double>& Gp, const function2D<double>& Gm,
			    const function1D<double>& Energy,
			    const vector<sLorentz>& lorentzm, const vector<sLorentz>& lorentzp)
{
  int m = omd.find_(0.0)+1;

  Gtx.resize(Na, omega.size());
  #pragma omp parallel for
  for (int i=0; i<Na; i++) if (Pexists[i]) aF[i].SetUp(Gp[i],omega);
  
  for (int i=0; i<m; i++){

    CalculateProducts(omd[i], fed[i], omega, Gm);

    #pragma omp parallel for
    for (int b=0; b<baths; b++){
      double sum=0;
      for (int j=0; j<Na; j++)
	for (map<int,double>::const_iterator l=common::sncab[j][b].begin(); l!=common::sncab[j][b].end(); l++){
	  int ind = l->first;
	  double prf = l->second/common::Ns[b];
	  if (ind>=0 && ind<Na) sum += prf*Cmp(j,ind);
	}
      A00(b,i) = sum/(M_PI*M_PI*common::Q);
    }
    if (common::cmp_susc){
      double sum=0;
      for (int j=0; j<Na; j++)
	for (map<int,double>::const_iterator l=common::suscb[j].begin(); l!=common::suscb[j].end(); l++){
	  int ind = l->first;
	  double prf = l->second;
	  if (ind>=0 && ind<Na)	sum += prf*Cmp(j,ind);
	}
      C00[i] = sum*th[i]/(M_PI*common::Q);
    }
  }

  #pragma omp parallel for
  for (int i=0; i<Na; i++) if (Pexists[i]) aF[i].SetUp(Gm[i],omega);
  
  for (int i=m; i<omd.size(); i++){
    
    CalculateProducts(omd[i], (1-fed[i]), omega, Gp);
    
    #pragma omp parallel for
    for (int b=0; b<baths; b++){
      double sum=0;
      for (int j=0; j<Na; j++)
	for (map<int,double>::const_iterator l=common::sncab[j][b].begin(); l!=common::sncab[j][b].end(); l++){
	  int ind = l->first;
	  double prf = l->second/common::Ns[b];
	  if (ind>=0 && ind<Na)	sum += prf*Cmp(j,ind);
	}
      A00(b,i) = sum/(M_PI*M_PI*common::Q);
    }
    if (common::cmp_susc){
      double sum=0;
      for (int j=0; j<Na; j++)
	for (map<int,double>::const_iterator l=common::suscb[j].begin(); l!=common::suscb[j].end(); l++){
	  int ind = l->first;
	  double prf = l->second;
	  if (ind>=0 && ind<Na)	sum += prf*Cmp(j,ind);
	}
      C00[i] = sum*th[i]/(M_PI*common::Q);
    }
  }

  if (common::SubtractLorentz){
    for (int b=0; b<baths; b++){
      //cout<<"Starting parallel part"<<endl;
      
      double* A00_private = new double[omd.size()];
      for (int s=0; s<omd.size(); s++) A00_private[s]=0.0;
      
      for (int i=0; i<Na; i++){
	for (map<int,double>::const_iterator l=common::sncab[i][b].begin(); l!=common::sncab[i][b].end(); l++){
	  int ind = l->first;
	  if (ind>=0 && ind<Na){
	    double prf = (l->second/common::Ns[b])/(M_PI*M_PI)/common::Q;
	    if (lorentzm[ind].exist){
              #pragma omp parallel for
	      for (int j=0; j<m; j++){
		double sum=0;
		for (int k=0; k<omega.size()-1; k++)
		  sum += lorentzm[ind].IntgAp(omega[k], omega[k+1], Gp(i,k), Gp(i,k+1), omd[j]);
		//A00(b,j) += sum*prf/fed[j];
		A00_private[j] += sum*prf/fed[j];
	      }
	    }
	    if (lorentzp[ind].exist){
              #pragma omp parallel for
	      for (int j=m; j<omd.size(); j++){
		double sum=0;
		for (int k=0; k<omega.size()-1; k++)
		  sum += lorentzp[ind].IntgAp(omega[k], omega[k+1], Gm(i,k), Gm(i,k+1), omd[j]);
		//A00(b,j) += sum*prf/(1-fed[j]);	       
		A00_private[j] += sum*prf/(1-fed[j]);	       
	      }
	    }
	    if (lorentzp[i].exist){
              #pragma omp parallel for
	      for (int j=0; j<m; j++){
		double sum=0;
		for (int k=0; k<omega.size()-1; k++)
		  sum += lorentzp[i].IntgAp(omega[k], omega[k+1], Gm(ind,k), Gm(ind,k+1), -omd[j]);
		//A00(b,j) += sum*prf/fed[j];
		A00_private[j] += sum*prf/fed[j];
	      }
	    }
	    if (lorentzm[i].exist){
              #pragma omp parallel for
	      for (int j=m; j<omd.size(); j++){
		double sum=0;
		for (int k=0; k<omega.size()-1; k++)
		  sum += lorentzm[i].IntgAp(omega[k], omega[k+1], Gp(ind,k), Gp(ind,k+1), -omd[j]);
		//A00(b,j) += sum*prf/(1-fed[j]);
		A00_private[j] += sum*prf/(1-fed[j]);
	      }
	    }
	    if (lorentzm[ind].exist && lorentzp[i].exist)
              #pragma omp parallel for
	      for (int j=0; j<m; j++){
		//A00(b,j) += lorentzm[ind].IntgApLL(lorentzp[i], omd[j]) * prf/fed[j];
		A00_private[j] += lorentzm[ind].IntgApLL(lorentzp[i], omd[j]) * prf/fed[j];
	      }
	    if (lorentzp[ind].exist && lorentzm[i].exist)
              #pragma omp parallel for
	      for (int j=m; j<omd.size(); j++){
		//A00(b,j) += lorentzp[ind].IntgApLL(lorentzm[i], omd[j]) * prf/(1-fed[j]);
		A00_private[j] += lorentzp[ind].IntgApLL(lorentzm[i], omd[j]) * prf/(1-fed[j]);
	      }
	  }
	}
      }
      for (int s=0; s<omd.size(); s++) A00(b,s) += A00_private[s];
      delete[] A00_private;
      //cout<<"Just ended parallel part"<<endl;
    }

    if (common::cmp_susc){
      for (int i=0; i<Na; i++){
	for (map<int,double>::const_iterator l=common::suscb[i].begin(); l!=common::suscb[i].end(); l++){
	  int ind = l->first;
	  if (ind>=0 && ind<Na){
	    double prf = (l->second)/(M_PI*common::Q);
	    if (lorentzm[ind].exist){
	      for (int j=0; j<m; j++){
		double sum=0;
		for (int k=0; k<omega.size()-1; k++)
		  sum += lorentzm[ind].IntgAp(omega[k], omega[k+1], Gp(i,k), Gp(i,k+1), omd[j]);
		C00[j] += sum*prf*th[j]/fed[j];
	      }
	    }
	    if (lorentzp[ind].exist){
	      for (int j=m; j<omd.size(); j++){
		double sum=0;
		for (int k=0; k<omega.size()-1; k++)
		  sum += lorentzp[ind].IntgAp(omega[k], omega[k+1], Gm(i,k), Gm(i,k+1), omd[j]);
		C00[j] += sum*prf*th[j]/(1-fed[j]);
	      }
	    }
	    if (lorentzp[i].exist){
	      for (int j=0; j<m; j++){
		double sum=0;
		for (int k=0; k<omega.size()-1; k++)
		  sum += lorentzp[i].IntgAp(omega[k], omega[k+1], Gm(ind,k), Gm(ind,k+1), -omd[j]);
		C00[j] += sum*prf*th[j]/fed[j];
	      }
	    }
	    if (lorentzm[i].exist){
	      for (int j=m; j<omd.size(); j++){
		double sum=0;
		for (int k=0; k<omega.size()-1; k++)
		  sum += lorentzm[i].IntgAp(omega[k], omega[k+1], Gp(ind,k), Gp(ind,k+1), -omd[j]);
		C00[j] += sum*prf*th[j]/(1-fed[j]);
	      }
	    }
	    if (lorentzm[ind].exist && lorentzp[i].exist)
	      for (int j=0; j<m; j++)
		C00[j] += lorentzm[ind].IntgApLL(lorentzp[i], omd[j]) * prf * th[j]/fed[j];
	    
	    if (lorentzp[ind].exist && lorentzm[i].exist)
	      for (int j=m; j<omd.size(); j++)
		C00[j] += lorentzp[ind].IntgApLL(lorentzm[i], omd[j]) * prf * th[j]/(1-fed[j]);
	  }
	}
      }
    }
  }

  
  if (common::pcore){
    // core stuff
    for (int b=0; b<baths; b++){
      for (int i=0; i<omd.size(); i++){
	double sum1=0;
	for (int j=0; j<Na; j++){
	  for (map<int,double>::const_iterator l=common::sncab[j][b].begin(); l!=common::sncab[j][b].end(); l++){
	    if (l->first >= Na){
	      int ind = l->first;
	      double x = Energy[ind]-common::lambda0-omd[i];
	      double prf = l->second/common::Ns[b];
	      sum1 -= prf*Gm[j](omega.Interp(x))/common::Q/M_PI;
	    }
	  }
	}
	double sum2=0;
	for (int j=Na; j<Na+Nc; j++){
	  for (map<int,double>::const_iterator l=common::sncab[j][b].begin(); l!=common::sncab[j][b].end(); l++){
	    if (l->first >= 0 && l->first<Na){
	      int ind = l->first;
	      double x = Energy[j]-common::lambda0+omd[i];
	      double prf = l->second/common::Ns[b];
	      sum2 -= prf*Gm[ind](omega.Interp(x))/common::Q/M_PI;
	    }
	  }
	}
	A00c(b,i) = sum1+sum2;
      }
    }
    // Checking doping!
    for (int b=0; b<baths; b++){
      double suma = 0, sumc = 0;
      for (int i=0; i<omd.size(); i++) {
	suma += A00(b,i)*fed[i]*omd.Dh(i);
	sumc += A00c(b,i)*fed[i]*omd.Dh(i);
      }
      suma *= common::Ns[b];
      sumc *= common::Ns[b];
      double miss_nd = common::nalpha[b]-(suma+sumc);
      double core_fact = 1.;
      if (sumc!=0 && common::renorm_core){
	core_fact = (common::nalpha[b]-suma)/sumc;
	if (core_fact<0) core_fact=0;
	if (core_fact>10) core_fact = 10;
	cout<<b<<" : "<<miss_nd<<" renormaliziang core part by "<<core_fact<<endl;
      }
      for (int i=0; i<omd.size(); i++) A00(b,i) += A00c(b,i)*core_fact;

      if (common::renorm){
	double suml=0, sumr=0;
	for (int i=0; i<omd.size(); i++){
	  suml += A00(b,i)*fed[i]*omd.Dh(i);
	  sumr += A00(b,i)*(1-fed[i])*omd.Dh(i);
	}
	int izero = omd.find_(0.0);
	double ml1=0, mr1=0;
	for (int i=0; i<izero; i++) {
	  ml1 += omd[i]*A00(b,i)*fed[i]*omd.Dh(i);
	  mr1 += omd[i]*A00(b,i)*(1-fed[i])*omd.Dh(i);
	}
	double ml2=0, mr2=0;
	for (int i=izero+1; i<omd.size(); i++) {
	  ml2 += omd[i]*A00(b,i)*fed[i]*omd.Dh(i);
	  mr2 += omd[i]*A00(b,i)*(1-fed[i])*omd.Dh(i);
	}
	double n0 = common::nalpha[b]/common::Ns[b];
	double C = (-ml2 + ml2*n0 + mr2*n0 - mr2*suml + ml2*sumr)/(ml1*mr2-ml2*mr1);
	double D = (ml1 - ml1*n0 - mr1*n0 + mr1*suml - ml1*sumr)/(ml1*mr2-ml2*mr1);

	if (1+C*omd[0]<0) C = -1/omd[0];
	if (1+D*omd.last()<0)  D = -1/omd.last();
	
	for (int i=0; i<izero; i++) A00(b,i) *= (1+C*omd[i]);
	for (int i=izero+1; i<omd.size(); i++) A00(b,i) *= (1+D*omd[i]);
      
	cout<<"Renormalizing A["<<b<<"] by "<<C<<", "<<D<<"at negative and positive frequency"<<endl;
      }
    }
  }
  
//   ofstream out("Aloc.imp"); out.precision(16);
//   for (int i=0; i<omd.size(); i++){
//     out<<setw(25)<<omd[i]<<" ";
//     for (int b=0; b<baths; b++) out<<setw(25)<<A00(b,i)<<" ";
//     out<<endl;
//   }
}

inline void Physical::KramarsKronig()
{
  for (int b=0; b<baths; b++) G00[b].KramarsKronig(omd,logod);
}

void Physical::CalcSelfEnergy()
{
  for (int b=0; b<baths; b++){
    for (int i=0; i<omd.size(); i++){
      //double Deltar = ::KramarsKronig(Ac[b], omd, omd[i], i, Ac[b][i]);
      //dcomplex Delta(-M_PI*Deltar,-M_PI*Ac[b][i]);
      Sig[b][i] = omd[i]-common::Ed[b]-Delta0[b][i]-1/G00[b][i];
      if (Sig[b][i].imag()>0) Sig[b][i].imag()=0.0;
    }
  }
  if (common::cmp_susc){
    for (int i=0; i<omd.size(); i++)
      Chi[i] = dcomplex(::KramarsKronig(C00, omd, omd[i], i, C00[i]),C00[i]);
  }
}

void Physical::Print(int n, string dir="")
{
  string filename;
  if (n<0) filename = common::outdir+"/A00"+dir;
  else filename = common::outdir+NameOfFile("/A00",n,3);
  ofstream out(filename.c_str()); out.precision(16);
  common::printHead(out)<<" peakposition=" << omd.dcenter()<<endl;
  for (int i=0; i<omd.size(); i++){
    out <<setw(25)<<omd[i];
    for (int b=0; b<baths; b++)
      out<<setw(25)<<A00[b][i]<<setw(25)<<G00[b][i]<<setw(25)<<-Sig[b][i];
    out<<endl;
  }

  if (n<0) filename = common::outdir+"/Susc"+dir;
  else filename = common::outdir+NameOfFile("/Susc",n,3);
  ofstream outs(filename.c_str()); outs.precision(16);
  common::printHead(outs)<<" peakposition=" << omd.dcenter()<<endl;
  for (int i=0; i<omd.size(); i++)
    outs <<setw(25)<<omd[i]<<setw(25)<<Chi[i]<<endl;
}

void Physical::Print0(const string& filename)
{
  ofstream out(filename.c_str()); out.precision(16);
  common::printHead(out)<<" peakposition=" << omd.dcenter()<<endl;
  for (int i=0; i<omd.size(); i++){
    out <<setw(25)<<omd[i];
    for (int b=0; b<baths; b++) out<<setw(25)<<A00[b][i];
    for (int b=0; b<baths; b++) out<<setw(25)<<G00[b][i];
    for (int b=0; b<baths; b++) out<<setw(25)<<-Sig[b][i];
    out<<endl;
  }
}

double Auxiliary::DeterminSelfEnergies(double alpha,int CmpDiff){
  double beta=1-alpha;
  Sigtmp.resize(om.size());
  if (CmpDiff<0) CmpDiff = Na;
  double diff=0, norm=0;
  for (int j=0; j<Na; j++){
    for (int i=0; i<om.size(); i++) if (Sigtn(j,i)>0) Sigtn(j,i)=0;
    for (int i=0; i<om.size(); i++) Sigtmp[i].imag() = Sigtn(j,i)*(1-fe[i]);
    Sigtmp.KramarsKronig(om, logo);
    for (int i=0; i<om.size(); i++){
      dcomplex Sigcn = Sigtmp[i] + Sigcore(j,i);
      Sigtn(j,i) += Sigcore(j,i).imag();

      if (j<CmpDiff){
	diff += fabs(Sigtn(j,i)-Sigt(j,i));
	norm += fabs(Sigt(j,i));
      }
      
      Sigt(j,i) = beta*Sigt(j,i)+alpha*Sigtn(j,i);
      Sigc(j,i) = beta*Sigc(j,i)+alpha*Sigcn;
    }
  }
  return diff/norm;
}

void Physical::DeterminG00(double alpha,ostream& loging)
{
  double beta=1-alpha;
  double alphapi=-alpha*M_PI;
  for (int b=0; b<baths; b++){
    for (int j=0; j<omd.size(); j++)
      G00[b][j].imag()=beta*G00[b][j].imag()+alphapi*A00[b][j];
    G00[b].KramarsKronig(omd,logod);
  }
  
  common::TrLogGimp=0.0;
  for (int b=0; b<baths; b++){
    double Ndf=0.0;
    double dsum=0;
    for (int j=0; j<omd.size(); j++){
      dsum += -log(-G00[b][j]).imag()*fed[j]*omd.Dh(j)/M_PI;
      Ndf += -G00[b][j].imag()*fed[j]*omd.Dh(j)/M_PI;
    }
    common::TrLogGimp += dsum*common::Ns[b];
    Ndf *= common::Ns[b];
    loging<<"Expected density:"<<common::nalpha[b]<<" numerical density:"<<Ndf<<endl;
  }
  loging<<"TrLogGimp="<<common::TrLogGimp<<endl;
}

void Auxiliary::PrintNorm(ostream& stream)
{
  stream<<"    Norm of Spectral functions: "<<endl<<"   ";
  stream.setf(ios::fixed);
  for (int i=0; i<Na; i++){
    double sum=0;
    for (int j=0; j<om.size(); j++)
      sum += Gp(i,j)*om.Dh(j);
    sum += lorentzp[i].P*M_PI;
    
    sum/=-M_PI;
    double norm0=1;
    stream<<setprecision(4)<<" ";
    
    if (fabs(sum-norm0)<1e-2)
      stream<<COLOR(GREEN,setw(2)<<i<<":"<<setw(8)<<sum)<<" ";
    else if (fabs(sum-norm0)<1e-1)
      stream<<COLOR(YELLOW,setw(2)<<i<<":"<<setw(8)<<sum)<<" ";
    else 
      stream<<COLOR(PURPLE,setw(2)<<i<<":"<<setw(8)<<sum)<<" ";
    if ((i+1)%6==0) stream<<endl<<"   ";
  }
  stream<<endl;
  for (int b=0; b<baths; b++){
    stream<<setprecision(4)<<" "<<COLOR(BLUE,setw(2)<<b<<":"<<setw(8)<<common::nalpha[b])<<" ";
  }
  stream<<endl;
  stream.unsetf(ios::fixed);
}

void Physical::PrintA00(ostream& out)
{
  out.precision(16);
  common::printHead(out)<<" peakposition=" << omd.dcenter()<<endl;
  for (int i=0; i<omd.size(); i++){
    out<<setw(25)<<omd[i];
    for (int b=0; b<baths; b++)
      out<<setw(25)<<A00[i];
    out<<endl;
  }
}

double Auxiliary::Difference(){
  double diff=0, norm=0;
  for (int j=0; j<Na; j++){
    for (int i=0; i<om.size(); i++){
      diff += fabs(Sigtn(j,i)-Sigt(j,i));
      norm += 0.5*fabs(Sigtn(j,i)+Sigtn(j,i));
    }
  }
  return diff/norm;
}

/******************* Used only for debugging **********************/
void Auxiliary::PrintSign()
{
  for (int i=0; i<Na; i++){
    ofstream out(NameOfFile("Sign",i,2).c_str());
    out.precision(16);
    for (int j=0; j<om.size(); j++)
      out<<setw(25)<<om[j]<<setw(25)<<-Sigtn[i][j]<<endl;
  }
}

void Auxiliary::Print_aAc(int l)
{
  for (int i=0; i<aAc[0].size_N(); i++){
    ofstream out(NameOfFile_("aAc",l,i,1,3).c_str());
    out.precision(16);
    for (int j=0; j<aAc[0].size_Nd(); j++){
      out<<setw(25)<<om[j]<<setw(25)<<aAc[0][i][j]/om.Dh(j)<<endl;
    }
  }
}

/******************* New things ******************************/
void common::ParsInputFile(const string& filename)
{
  ifstream input(filename.c_str());
  string line;
  getline(input,line);
  input>>baths;
  Ns.resize(baths);
  for (int i=0; i<baths; i++) input>>Ns[i];
  input>>Na;
  input>>Nc;
  getline(input,line); getline(input,line);
  if (!input){ cerr<<filename<<" file not recognized. Error in first 3 lines!"<<endl; exit(1);}
  deg.resize(Na+Nc);
  Ms.resize(Na+Nc,baths);
  Mtot.resize(Na+Nc);
  sJc.resize(Na+Nc);
  ncab.resize(Na+Nc, baths);
  ncaf.resize(Na+Nc, baths);
  prefactb.resize(Na+Nc, baths);
  prefactf.resize(Na+Nc, baths);
  prefactG.resize(Na+Nc, baths);
  sncab.resize(Na+Nc);
  sncaf.resize(Na+Nc);
  for (int i=0; i<Na+Nc; i++) sncab[i].resize(baths);
  for (int i=0; i<Na+Nc; i++) sncaf[i].resize(baths);
  
  vector<int> Nncab(baths), Nncaf(baths);
  for (int i=0; i<Na+Nc; i++){
    getline(input, line);
    if (!input){ cerr<<filename<<" file not recognized. Error in line number "<<i+3<<endl; exit(1);}
    stringstream thisline(line);
    int lc;
    thisline>>lc;
    for (int j=0; j<baths; j++) thisline>>Ms[i][j];
    thisline>>Mtot[i]>>deg[i]>>sJc[i];
    for (int j=0; j<baths; j++) thisline>>Nncab[j];
    for (int j=0; j<baths; j++) thisline>>Nncaf[j];

    string cross; double fct; int ind;
    for (int j=0; j<baths; j++){
      for (int k=0; k<Nncab[j]; k++){
	thisline>>fct>>cross>>ind;
	sncab[i][j][ind]=fct;
      }
    }
    for (int j=0; j<baths; j++){
      for (int k=0; k<Nncaf[j]; k++){
	thisline>>fct>>cross>>ind;
	sncaf[i][j][ind]=fct;
      }
    }
    if (!input){ cerr<<filename<<" file not recognized. Error in line number "<<i+3<<endl; exit(1);}
  }
  getline(input, line);// comment
  cmp_susc = false;
  if (input){
    suscb.resize(Na);
    for (int i=0; i<Na; i++){
      getline(input, line);
      if (!input) goto exit_loop;
      stringstream thisline(line);
      int lc;
      thisline>>lc;
      int ndiagram;
      thisline>>ndiagram;
      string cross; double fct; int ind;
      for (int j=0; j<ndiagram; j++){
	thisline>>fct>>cross>>ind;
	suscb[i][ind]=fct;
      }
    }
    cmp_susc = true;
  }
 exit_loop:
  
  PrintParsedData(cout);
  totDeg = 0;
  for (int i=0; i<Na; i++) totDeg += deg[i];
}

void common::PrintParsedData(ostream& stream)
{
  stream<<baths<<" ";
  for (int i=0; i<baths; i++) stream<<Ns[i]<<" ";
  stream<<Na<<" "<<Nc<<endl;
  for (int i=0; i<Na+Nc; i++){
    stream<<setw(3)<<i<<" ";
    if (i<Na) stream<<"v ";
    else stream<<"c ";
    for (int j=0; j<baths; j++) stream<<setw(10)<<Ms[i][j];
    stream<<setw(4)<<Mtot[i]<<setw(5)<<deg[i]<<setw(6)<<sJc[i];
    
    for (int b=0; b<baths; b++) stream<<setw(2)<<sncab[i][b].size()<<" ";
    for (int b=0; b<baths; b++) stream<<setw(2)<<sncaf[i][b].size()<<" ";
    
    for (int b=0; b<baths; b++)
      for (map<int,double>::const_iterator l=sncab[i][b].begin(); l!=sncab[i][b].end(); l++)
	stream<<setw(6)<<l->second<<" x "<<setw(4)<<left<<l->first<<right;

    for (int b=0; b<baths; b++)
      for (map<int,double>::const_iterator l=sncaf[i][b].begin(); l!=sncaf[i][b].end(); l++)
	stream<<setw(6)<<l->second<<" x "<<setw(4)<<left<<l->first<<right;

    stream<<endl;
  }
  if (!cmp_susc) return;
  stream<<"Susceptibility digrams:"<<endl;
  for (int i=0; i<Na; i++){
    stream<<setw(3)<<i<<" ";
    for (map<int,double>::const_iterator l=suscb[i].begin(); l!=suscb[i].end(); l++)
      stream<<setw(6)<<l->second<<" x "<<setw(4)<<left<<l->first<<right;
    stream<<endl;
  }
}

void print(std::ostream& stream, const mesh1D& om, const function2D<dcomplex>& f, int width=20)
{
  if (om.size()!=f.size_Nd()) std::cerr<<"Can't print objectc of different size!"<<std::endl;
  for (int i=0; i<om.size(); i++){
    stream <<std::setw(width)<<om[i];
    for (int j=0; j<f.size_N(); j++) stream<<std::setw(width)<<f(j,i);
    stream<<std::endl;
  }
}

void Physical::MissingDoping(double start)
{
  cout<<"Missing doping : ";
  for (int b=0; b<baths; b++){
    double sum = 0;
    for (int i=0; i<omd.size(); i++) {
      if (omd[i]>start) sum += G00[b][i].imag()*fed[i]*omd.Dh(i);
    }
    sum *= -common::Ns[b]/M_PI;
    common::miss_nd[b] = common::nalpha[b]-sum;
    cout<<b<<" : "<<common::miss_nd[b]<<" ";
  }
  cout<<endl;

  common::Sinfty.resize(baths);
  for (int b=0; b<baths; b++){
    double sum0 = 0, sum1 = 0;
    for (int i=0; i<omd.size(); i++) {
      sum0 += A00(b,i)*omd.Dh(i);
      sum1 += A00(b,i)*omd[i]*omd.Dh(i);
    }
    common::moment[b][0] = sum0;
    common::moment[b][1] = sum1;
    common::Sinfty[b] = sum1/sum0-common::Ed[b];
  }
}


void Auxiliary::PrintCore(const string& filename)
{
  ofstream out(filename.c_str());
  for (int i=0; i<om.size(); i++){
    out<<setw(20)<<om[i]<<" ";
    for (int j=0; j<Na; j++){
      out<<setw(20)<<Sigcore[j][i]<<" ";
    }
    out<<endl;
  }
}


#endif
