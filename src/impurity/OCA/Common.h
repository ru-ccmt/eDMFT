// @Copyright 2007 Kristjan Haule
// 
#ifndef _COMMON_
#define _COMMON_
#include "zeroin.h"
#include "average.h"
#include "definitions.h"
#include <list>
#include <map>
#include <deque>

using namespace std;



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
  static deque<OCAd> ocaPhi;
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
  static bool renorm_core, renorm;
  static double Fimp, Epot, TrLogGimp;
  static void SetParameters(Par<double>& Ed_, double U_, /*double J_, */double T_, double Q0_, const string& outdir_, int N_ac_,
			    double dom_ac_, int acore_, int pcore_, bool renorm_core_, bool renorm_)
  {
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
  static void ParsInputFile(const string& filename, ofstream& log);
  static void PrintParsedData(ostream& stream);
  static ostream& printHead(ostream& stream);
};

// Auxiliary self-energies and spectral functions
class Auxiliary{
  const int Na, Nc, baths;
  mesh1D om;
  function1D<double> fe;
  function1D<double> fedh;
  function1D<double> logo;
  function2D<double> Sigt;
  function2D<double> Sigtn_;
  function2D<dcomplex> Sigc;
  function2D<Gmpr> G;
  function1D<double> Energy;
  //  int StrFreq, EndFreq;
  function1D<double> Probability;
  function2D<dcomplex> Sigcore;
  function2D<double> Gt;
  function2D<double> Gp;
  function2D<double> Gm;
  vector<function2D<double> > aAc;
  function1D<double> Acx;
  function1D<double> Acy;
  function2D<double> Acp, Acm;
  AvFun<double> aF;
  function2D<dcomplex> Deltam_ac, Deltap_ac;
  function1D<double> mom_Deltam_ac, mom_Deltap_ac;
  mesh1D oml;
  function2D<double> GtA1, GtA2;
  int mpos;
  
public:
  Auxiliary (int Na_, int Nc_, int baths_) : Na(Na_), Nc(Nc_), baths(baths_), aAc(2*baths), mom_Deltam_ac(baths), mom_Deltap_ac(baths), Probability(Na){}
  bool ReadSelfEnergy(const string& filename, const Par<double>& Ed, const Par<double>& T, const Par<double>& U, const mesh1D& ph_omd, const function2D<double>& ph_Ac, ofstream& log);
  void KramarsKronig();
  double DeterminSpectralFunctions(double StartLambda, double EndLambda, double dLamdba, int followPeak, ofstream& log);
  void PrintOutMeanQ(double StartLambda, double EndLambda);
  void PrintNorm(ostream& stream);
  void Print(int l);
  void Print(string filename, bool PrtSF);
  void SetSignToZero(){Sigtn_=0.0;Sigcore=0.0;}
  
  void Calc_OCA_SelfEnergy(int i, const mesh1D& omd, const function2D<Amp>& pAc, const mesh1D& momd, const function2D<Amp>& mAc);
  void Calc_NCA_Sigmab(const mesh1D& omd);
  void Calc_NCA_Sigmaf(const mesh1D& omd);
  function1D<double> Calc_OCA_A00(double omega, const mesh1D& momd, const function2D<Amp>& mAc);

  void SetUpAverageAc(const mesh1D& omd, const mesh1D& momd, const function2D<double>& Ack, const function1D<double>& fed);
  
  double Difference();
  void DeterminSelfEnergies(double alpha);
  const mesh1D& omega() const {return om;}
  double ferm(int i) const {return fe[i];}
  
  double Q(double lambda);
  double operator()(double lambda);
  //  void Start(int StrFreq_) {StrFreq = StrFreq_;}
  //  void Stop (int EndFreq_) {EndFreq = EndFreq_;}
  double minEnergy;
  
  const function2D<double>& _Gp() const {return Gp;}
  const function2D<double>& _Gm() const {return Gm;}
  const function1D<double>& Energ() const{return Energy;}
  void PrintSign();

  double* Sign(){return Sigtn_.MemPt();} // for MPI-calls
  
  void CreateSigma000(const mesh1D& omd, const function2D<double>& Ac);
  
private:
  void PrintOutMeanQ(int M, double StartLambda, double EndLambda);
  void MakeOneHsList(const deque<OCAd>& oca_diag, list<pair<int,int> >& hslist, int pp1, int pp2);
  void MakeTwoHsList(const deque<OCAd>& oca_diag, list<pair<int,int> >& hslist1, list<pair<int,int> >& hslist2);
  int MakeCombFunc(double omega, const mesh1D& omd, const function2D<Amp>& Ac, int j, const list<pair<int,int> >& xs,
		   map<pair<int,int>,int>& index_1, mesh1D& eps, function1D<Hxy>* Hx);
  double CalcUNCAVertex(double omega, int ipn, double center_inside, map<pair<int,int>,int>& index_1, 
			const mesh1D& eps, const mesh1D& eps_inside, const function1D<Hxy>* Hx, 
			const diag& d, int b1, int b2);
//   double CalcNCABubble(int ip, const list<pair<int,int> >& hslist, map<pair<int,int>,int>& index_1, const mesh1D& eps,
// 		       const function1D<Hxy>* Hx, const function2D<double>& pre_nca_x, function2D<int>& nca_xcopy);
//   double CalcNCABubbleInside(const mesh1D& eps, const function1D<Hxy>& Hx);
  void MakeLsList(const deque<OCAd>& oca_diag, list<pair<int,int> >& lslist_fb, list<pair<int,int> >& lslist_af);
  void MakeCombFunc(double omega, const list<pair<int,int> >& lslist, map<pair<int,int>,int>& index_1, mesh1D& eps, function1D<Lxy>* Lx);
};

// Physical electron spectral function and suscpetibility
// Physical observables
class Physical{
  const int Na, Nc, baths;
public:  
  mesh1D omd;
  function2D<dcomplex> G00;
  function2D<double> A00;
  function2D<double> A00c;
  function2D<dcomplex> Sig;
private:
  mesh1D momd;
  function1D<double> fed;
  function1D<double> logod;
  function2D<double> Ac;
  function2D<dcomplex> Delta0;
  function2D<Amp> pAc;
  function2D<Amp> mAc;
  //  int StrFreq, EndFreq;
  function1D<bool> Pexists;
  vector<AvFun<double> > aF;
  function2D<double> Gtx;
  function2D<double> Cmp;
public:
  Physical(int Na_, int Nc_, int baths_);
  bool ReadBathFunction(const string& filename, ofstream& log, bool spectra);
  void KramarsKronig();
  void DeterminG00(ofstream& loging);
  double Difference();
  void Print(const string& filename);
  
  const mesh1D& omega() const {return omd;}
  const mesh1D& momega() const {return momd;}
  const function2D<double>& Ac0() const {return Ac;}
  const function1D<double>& fe() const {return fed;}
  const function2D<Amp>& pA() const {return pAc;}
  const function2D<Amp>& mA() const {return mAc;}
  
  void PrintA00();
  void CalcSelfEnergy();
  void SetA00(int i, const function1D<double>& a00);
  
  //  void Start(int StrFreq_) {StrFreq = StrFreq_;}
  //  void Stop (int EndFreq_) {EndFreq = EndFreq_;}
  void MissingDoping(double start, ofstream& log);
  void Calculate_NCA_A00(const mesh1D& omega, const function2D<double>& Gp, const function2D<double>& Gm, const function1D<double>& Energy);
  void AddCore(const mesh1D& omega, const function2D<double>& Gp, const function2D<double>& Gm, const function1D<double>& Energy, ofstream& log);
  void ReadAlocImp(const string& filename);
  double* P_A00(){return A00.MemPt();}
private:
  bool ReadBeginning(const string& filename, istream& input, int& n, int& m, bool& begincomment, double& center, ofstream& log);
  void CalculateProducts(double u, double fu, const mesh1D& om, const function2D<double>& Gm);
};

bool Auxiliary::ReadSelfEnergy(const string& filename, const Par<double>& Ed, const Par<double>& T, const Par<double>& U, const mesh1D& ph_omd, const function2D<double>& ph_Ac, ofstream& log){
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
  //  for (int i=0; i<baths; i++) {Ed.next(); Ed_[i] = Ed;}
  //  Ed_[0] = SpecNumber;
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
    Energy[i] += 0.5*common::Mtot[i]*(common::Mtot[i]-1)*(common::U/*-0.5*common::J*/);
    Energy[i] += /*common::J*/common::sJc[i];
    if (Energy[i]<minEnergy) minEnergy = Energy[i];
  }

  log<<"************* Parameters ****************"<<endl;
  log<<"  		U  = "<<common::U<<endl;
  for (int i=0; i<baths; i++)
    log<<"  		Ed"<<i<<" = "<<common::Ed[i]<<endl;
  log<<"  		T  = "<<common::T<<endl;
  for (int i=0; i<baths; i++)
    log<<"  		N"<<i<<" = "<<common::Ns[i]<<endl;
  for (int i=0; i<Na+Nc; i++){
    if (i<Na) log<<" valence state"<<setw(2)<<left<<i<<right<<" = ";
    else log<<" core    state"<<i<<" = ";
    for (int j=0; j<baths; j++) log<<setw(10)<<common::Ms[i][j]<<" ";
    log<<" with Energy"<<setw(2)<<left<<i<<right<<" = "<<Energy[i]<<endl;
  }
  log<<"*****************************************"<<endl;
  
  // Computes the number of columns in file
  if (!input) {
    cerr << "ERROR: Wrong file format for Sigm" << endl;
    return false;
  }
  getline(input,str);  n++;
  istringstream oneline(str);
  int m=0; double t;
  while (oneline){oneline>>t; m++;}
  m--;
  while (input){ getline(input,str); n++;}
  n--;
  
  log << filename << ": Number of entries: "<< n <<endl;
  log << filename << ": Number of columns: "<< m <<endl;
  log << filename << ": Peak-position "<< center <<endl;

  bool CreateDefault = false;
  if (m<2*Na+1){
    //    cerr<<"ERROR: Not enough columns is input Sigma file. Exiting!"<<endl;
    //    return false;
    clog<<"WARRNING: Not enough columns is input self-energy for pseudoparticles.... Creating default!"<<endl;
    CreateDefault = true;
  }
  inputf.seekg(0,ios::beg);
  //  log<<"Premaknil na "<< inputf.tellg()<<endl;
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

  GtA1.resize(Na,mpos);
  GtA2.resize(Na,om.size()-mpos);
  Sigcore.resize(Na,om.size());
  Gp.resize(Na,om.size());
  Gm.resize(Na,om.size());
  
  Sigc.resize(Na,om.size());
  Sigtn_.resize(om.size(),Na);
  G.resize(Na,om.size());
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
    double sum=0;
    for (int i=0; i<om.size(); i++)
      sum -= Sigt(j,i)*fedh[i]/(sqr(om[i]+mune-Sigc(j,i).real())+sqr(Sigc(j,i).imag()));
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

double Auxiliary::DeterminSpectralFunctions(double StartLambda, double EndLambda, double dLambda, int followPeak, ofstream& log)
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
  log << setprecision(16) << "; lambda = "<<lambda0<<endl;

  double sumQ = 0, sumnd=0;
  function1D<double> dQ(Na);
  function1D<double> Gi(om.size());
  for (int j=0; j<Na; j++){
    double mune = -Energy[j]+lambda0;
    dQ[j]=0;
    for (int i=0; i<om.size(); i++){
      double Gt = Sigt(j,i)/(sqr(om[i]+mune-Sigc(j,i).real())+sqr(Sigc(j,i).imag()));
      G(j,i).p = fe[i]*Gt;
      G(j,i).m = (1-fe[i])*Gt;

      Gm(j,i) = G(j,i).p;
      Gp(j,i) = G(j,i).m;
      
      Gi[i] = G(j,i).m;
      dQ[j] -= Gt*fedh[i];
    }
    for (int i=0; i<om.size(); i++)
      G(j,i).r = ::KramarsKronig(Gi, om, om[i], i, Gi[i]);
    
    dQ[j] *= common::deg[j]/M_PI;
    sumQ += dQ[j];
    sumnd += dQ[j]*common::Mtot[j];
  }
  log<<"       Q = "<<sumQ<<endl;
  for (int j=0; j<Na; j++){
    Probability[j] = dQ[j]/sumQ;
    log<<setprecision(16)<<"       n"<<j<<"="<<dQ[j]/sumQ<<endl;
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
  log<<"        Fimp="<<common::Fimp<<" Epot="<<common::Epot<<" Epot+OneP="<<Epot<<endl;
  log<<"        Q is here equal to "<<sumQ<<endl;
  return sumnd/sumQ;
}

void Auxiliary::Print(string filename, bool PrtSF=true)
{
  ofstream out1(filename.c_str());  out1.precision(16);
  common::printHead(out1)<<" peakposition="<<om.dcenter()<<endl;
  //  for (int i=StrFreq; i<EndFreq; i++){
  for (int i=0; i<om.size(); i++){
    out1<<setw(25)<<om[i];
    for (int j=0; j<Na; j++) out1<<setw(25)<<Sigc(j,i).real()<<" "<<setw(25)<<-Sigt(j,i);
    out1<<endl;
  }

  if (PrtSF){
    char buffer[200]; strcpy(buffer, filename.c_str());
    string filename1 = string(buffer) + "~";
    ofstream out2(filename1.c_str());  out2.precision(16);
    common::printHead(out2)<<" peakposition="<<om.dcenter()<<endl;
    //    for (int i=StrFreq; i<EndFreq; i++){
    for (int i=0; i<om.size(); i++){
      out2<<setw(25)<<om[i];
      for (int j=0; j<Na; j++) out2<<setw(25)<<-(G(j,i).p+G(j,i).m);
      for (int j=0; j<Na; j++) out2<<setw(25)<<-Gp(j,i);
      for (int j=0; j<Na; j++) out2<<setw(25)<<-Gm(j,i);
      out2<<endl;
    }
  }
}
void Auxiliary::Print(int l)
{ Print(NameOfFile(common::outdir+"/Sigma", l));}


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
    double Deltar;
    int ofst=0;
    for (int i=0; i<common::N_ac; i++){
      Deltar = ::KramarsKronig(Acm[b], omd, oml[i], 0, 0.0);
      Deltam_ac[b][i] = dcomplex(-M_PI*Deltar,0.0);
      Deltar = ::KramarsKronig(Acp[b], omd, oml[i], 0, 0.0);
      Deltap_ac[b][i] = dcomplex(-M_PI*Deltar,0.0);
    }
    ofst=common::N_ac;
    for (int i=0; i<omd.size(); i++){
      Deltar = ::KramarsKronig(Acm[b], omd, omd[i], i, Acm[b][i]);
      Deltam_ac[b][ofst+i] = dcomplex(-M_PI*Deltar,-M_PI*Acm[b][i]);
      Deltar = ::KramarsKronig(Acp[b], omd, omd[i], i, Acp[b][i]);
      Deltap_ac[b][ofst+i] = dcomplex(-M_PI*Deltar,-M_PI*Acp[b][i]);
    }
    ofst=common::N_ac+omd.size();
    for (int i=0; i<common::N_ac; i++){
      Deltar = ::KramarsKronig(Acm[b], omd, oml[omd.size()+common::N_ac+i], omd.size()-1, 0.0);
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

void Auxiliary::MakeOneHsList(const deque<OCAd>& oca_diag, list<pair<int,int> >& hslist, int pp1, int pp2)
{
  for (deque<OCAd>::const_iterator il=oca_diag.begin(); il!=oca_diag.end(); il++){
    pair<int,int> p0(il->states[pp1],il->ib1);
    pair<int,int> p2(il->states[pp2],il->ib2);
    if (find(hslist.begin(),hslist.end(), p0) == hslist.end()) hslist.push_back(p0);
    if (find(hslist.begin(),hslist.end(), p2) == hslist.end()) hslist.push_back(p2);
  }
}

void Auxiliary::MakeTwoHsList(const deque<OCAd>& oca_diag, list<pair<int,int> >& hslist1, list<pair<int,int> >& hslist2)
{
  for (deque<OCAd>::const_iterator il=oca_diag.begin(); il!=oca_diag.end(); il++){
    pair<int,int> p2(il->states[2],il->ib2);
    pair<int,int> p0(il->states[0],il->ib1);
    if (find(hslist1.begin(),hslist1.end(), p2) == hslist1.end()) hslist1.push_back(p2);
    if (find(hslist2.begin(),hslist2.end(), p0) == hslist2.end()) hslist2.push_back(p0);
    if (il->ib2==il->ib1) continue;
    pair<int,int> p2n(il->states[2],il->ib1);
    pair<int,int> p0n(il->states[0],il->ib2);
    if (find(hslist1.begin(),hslist1.end(), p2n) == hslist1.end()) hslist1.push_back(p2n);
    if (find(hslist2.begin(),hslist2.end(), p0n) == hslist2.end()) hslist2.push_back(p0n);
  }
}

int Auxiliary::MakeCombFunc(double omega, const mesh1D& omd, const function2D<Amp>& Ac, int j,
			    const list<pair<int,int> >& hslist,
			    map<pair<int,int>,int>& index_1, mesh1D& eps, function1D<Hxy>* Hx)
{
  for (list<pair<int,int> >::const_iterator i=hslist.begin(); i!=hslist.end(); i++,j++){
    index_1[*i]=j;
    int ptn = i->first;
    int bth = i->second;
    MakeCombinedFunctions(eps, Hx[j], om, G[ptn], 0.0, omd, Ac[bth], omega);
#ifdef HEAVY_DEBUG    
    clog<<"1 Making combined functions for "<<ptn<<" "<<bth<<" at position "<<j<<endl;
#endif    
  }
  return j;
}

double Auxiliary::CalcUNCAVertex(double omega, int ipn, double center_inside, map<pair<int,int>,int>& index_1, 
				 const mesh1D& eps, const mesh1D& eps_inside, const function1D<Hxy>* Hx, 
				 const diag& d, int b1, int b2)
  
{
  double sum=0;
  for (int i=0; i<eps.size(); i++){
    double x[3] = {0.0, omega, omega-eps[i]};
    double peak[3] = {om.dcenter()+x[0], center_inside+x[1], om.dcenter()+x[2]};
    int j = index_1[make_pair(d[2],b2)];
    Hxy hg = Integrate2c1(eps_inside, Hx[j], om, G[d[1]], x[0], x[1], x[2], peak[0], peak[1], peak[2]);
    j = index_1[make_pair(d[0],b1)];
    sum += Product(Hx[j][i], hg)*eps.Dh(i);
  }
#ifdef _DEBUG2
  clog<<"2 Calculating Vertex for "<<ipn<<" "<<d[0]<<" "<<d[1]<<" "<<d[2]<<" baths: "<<b1<<" "<<b2<<" res="<<sum<<endl;
  //  clog<<" index: "<<index_1[make_pair(d[0],b1)]<<" "<<index_1[make_pair(d[2],b2)]<<endl;
#endif	
  return sum;
}

void Auxiliary::Calc_OCA_SelfEnergy(int iom, const mesh1D& omd, const function2D<Amp>& pAc, const mesh1D& momd, const function2D<Amp>& mAc)
{
#ifdef _DEBUG1
  static ofstream outp("mvertex.txt"); outp.precision(16);
#endif	

  double omega = om[iom];
  
  function1D<double> Sigtv(Na);
  for (int i=0; i<Na; i++) Sigtv[i]=0;

  {
    // back-back diagrams
    list<pair<int,int> > fs;
    // Prepares a list of pairs defining the combination of functions that are needed
    // for example G[i]*A_c[b] with any i or b. 
    MakeOneHsList(common::ocaPhi, fs, 1, 3);

    function1D<Hxy>* Hf = new function1D<Hxy>[fs.size()];
    map<pair<int,int>,int> index_1;
    mesh1D eps;
    // Actually makes this combination of functions defined above
    MakeCombFunc(omega, omd, pAc, 0, fs, index_1, eps, Hf);
    // And finally calculates the integrals
    for (deque<OCAd>::const_iterator il=common::ocaPhi.begin(); il!=common::ocaPhi.end(); il++){
      int i = il->states[0];
      diag d(il->states[1],il->states[2],il->states[3]);
      double fact = il->f/common::deg[i];
      double tt = fact*CalcUNCAVertex(omega, i, omd.dcenter(), index_1, eps, eps, Hf, d, il->ib1, il->ib2);
      Sigtv[i] += tt;
    }
    delete[] Hf;
  }
  {
    // forward-forward diagrams
    list<pair<int,int> > fs;
    MakeOneHsList(common::ocaPhi, fs, 3, 1);
    
    function1D<Hxy>* Hf = new function1D<Hxy>[fs.size()];
    map<pair<int,int>,int> index_1;
    mesh1D meps;
    MakeCombFunc(omega, momd, mAc, 0, fs, index_1, meps, Hf);
    for (deque<OCAd>::const_iterator il=common::ocaPhi.begin(); il!=common::ocaPhi.end(); il++){
      int i = il->states[2];
      diag d(il->states[3],il->states[0],il->states[1]);
      double fact = il->f/common::deg[i];
      double tt = fact*CalcUNCAVertex(omega, i, momd.dcenter(), index_1, meps, meps, Hf, d, il->ib1, il->ib2);
      Sigtv[i] += tt;
    }
    delete[] Hf;
  }
  {
    // back-forward diagrams
    list<pair<int,int> > bs;
    list<pair<int,int> > as;
    MakeTwoHsList(common::ocaPhi, as, bs);
    
    function1D<Hxy>* Hba = new function1D<Hxy>[bs.size()+as.size()];
    map<pair<int,int>,int> index_1;
    mesh1D eps, meps;
    int j = MakeCombFunc(omega, momd, mAc, 0, bs, index_1, meps, Hba);
    MakeCombFunc(omega, omd, pAc, j, as, index_1, eps, Hba);
    
    for (deque<OCAd>::const_iterator il=common::ocaPhi.begin(); il!=common::ocaPhi.end(); il++){
      int i = il->states[1];
      diag d(il->states[2],il->states[3],il->states[0]);
      double fact = il->f/common::deg[i];
      double val0 = fact*CalcUNCAVertex(omega, i, momd.dcenter(), index_1, eps, meps, Hba, d, il->ib2, il->ib1);
      Sigtv[i] += val0;
      if (il->states[1]==il->states[3] && il->ib2==il->ib1){
	Sigtv[i] += val0;
      }else{
	int i1 = il->states[3];
	diag d(il->states[2],il->states[1],il->states[0]);
	double fact1 = il->f/common::deg[i1];
	double val1 = fact1*CalcUNCAVertex(omega, i1, momd.dcenter(), index_1, eps, meps, Hba, d, il->ib1, il->ib2);
	Sigtv[i1] += val1;
      }
    }
    delete[] Hba;
  }

#ifdef _DEBUG1
  outp<<omega<<" ";
  for (int i=0; i<Na; i++) outp<<-Sigtn_[iom][i]<<" "<<-Sigtv[i]<<" ";
  outp<<endl;
  clog<<"-------------------------------------------"<<endl;
#endif
  
  for (int i=0; i<Na; i++) Sigtn_[iom][i] += Sigtv[i];
}


void Auxiliary::DeterminSelfEnergies(double alpha){
  double beta=1-alpha;
  for (int j=0; j<Na; j++)
    //    for (int i=StrFreq; i<EndFreq; i++)
    for (int i=0; i<om.size(); i++)
      Sigt(j,i) = beta*Sigt(j,i)+alpha*Sigtn_(i,j);
}

void Auxiliary::PrintNorm(ostream& stream)
{
  stream<<"    Norm of Spectral functions: "<<endl<<"   ";
  stream.setf(ios::fixed);
  for (int i=0; i<Na; i++){
    double sum=0;
    for (int j=0; j<om.size(); j++)
      sum += G(i,j).m*om.Dh(j);

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
  }
}

void Physical::SetA00(int i, const function1D<double>& a00)
{
  for (int j=0; j<baths; j++) A00[i][j] += a00[j];
}

bool Physical::ReadBeginning(const string& filename, istream& input, int& n, int& m, bool& begincomment, double& center, ofstream& log)
{
  if (!input) {
    log << "Can't open input file: " << filename << endl;
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
    log << "ERROR: Wrong file format for Sigm" << endl;
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

  log << filename << ": Number of entries: "<< n <<endl;
  log << filename << ": Number of columns: "<< m <<endl;
  log << filename << ": Peak-position "<< center <<endl;
  
  input.seekg(0, ios::beg);
  input.clear();
  if (begincomment) getline(input, str);
  
  return true;
}

bool Physical::ReadBathFunction(const string& filename, ofstream& log, bool spectra=true)// spectra=true: only spectral function will be read not the retarded quantity
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
  istringstream oneline(str);
  int m=0; double t;
  while (oneline){oneline>>t; m++;}
  m--;
  while (input){ getline(input,str); n++;}
  n--;

  log << filename << ": Number of entries: "<< n <<endl;
  log << filename << ": Number of columns: "<< m <<endl;
  log << filename << ": Peak-position "<< center <<endl;

  int number_cols = baths+1;
  if (!spectra) number_cols = 2*baths+1;
  
  if (m<number_cols){
    cerr<<"ERROR: Not enough columns in bath input file! Exiting..."<<endl;
    return false;
  }
  inputf.seekg(0, ios::beg);
  log<<"Premaknil na "<< inputf.tellg()<<endl;
  if (begincomment) inputf.ignore(10000,'\n');
  if (!inputf){ cerr<<"Reopening didn't suceeded!"<<endl; return false;}
  
  omd.resize(n);
  momd.resize(n);
  G00.resize(baths,n);
  A00.resize(n,baths);
  A00c.resize(n,baths);
  Sig.resize(baths,n);
  Ac.resize(baths,n);
  pAc.resize(baths,n);
  mAc.resize(baths,n);
  Delta0.resize(baths,n);
  
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

  for (int b=0; b<baths; b++){
    for (int i=0; i<omd.size(); i++){
      pAc(b,i).p =    fed[i] *Ac(b,i);
      pAc(b,i).m = (1-fed[i])*Ac(b,i);
      mAc(b,n-i-1).p = pAc(b,i).m;
      mAc(b,n-i-1).m = pAc(b,i).p;
    }
  }
  
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

void Physical::CalcSelfEnergy()
{
  for (int b=0; b<baths; b++){
    for (int i=0; i<omd.size(); i++){
      //      double Deltar = ::KramarsKronig(Ac[b], omd, omd[i], i, Ac[b][i]);
      //      dcomplex Delta(-M_PI*Deltar,-M_PI*Ac[b][i]);
      Sig[b][i] = omd[i]-common::Ed[b]-Delta0[b][i]-1/G00[b][i];
      if (Sig[b][i].imag()>0) Sig[b][i].imag()=0.0;
    }
  }
}

void Physical::Print(const string& filename)
{
  ofstream out(filename.c_str()); out.precision(16);
  common::printHead(out)<<" peakposition=" << omd.dcenter()<<endl;
  for (int i=0; i<omd.size(); i++){
    out <<setw(25)<<omd[i];
    for (int b=0; b<baths; b++) out<<setw(25)<<A00[i][b];
    for (int b=0; b<baths; b++) out<<setw(25)<<G00[b][i];
    for (int b=0; b<baths; b++) out<<setw(25)<<-Sig[b][i];
    out<<endl;
  }
}


void Physical::DeterminG00(ofstream& loging)
{
  for (int b=0; b<baths; b++){
    for (int j=0; j<omd.size(); j++) G00[b][j].imag()= -A00[j][b]*M_PI;
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

void Physical::PrintA00()
{
  ofstream out("A00.current");
  out.precision(16);
  common::printHead(out)<<" peakposition=" << omd.dcenter()<<endl;
  for (int i=0; i<omd.size(); i++){
    out<<setw(25)<<omd[i];
    for (int b=0; b<baths; b++)
      out<<setw(25)<<A00[i][b];
    out<<endl;
  }
}

double Auxiliary::Difference(){
  double diff=0, norm=0;
  for (int j=0; j<Na; j++){
    for (int i=0; i<om.size(); i++){
      diff += fabs(Sigtn_(i,j)-Sigt(j,i));
      norm += fabs(Sigtn_(i,j));
    }
  }
  return diff/norm;
}

ostream& common::printHead(ostream& stream)
{
  stream<<"# ";
  stream<<" nb="<<baths<<" ";

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

void common::ParsInputFile(const string& filename, ofstream& log)
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
  getline(input, line); // comment
  while (input){
    OCAd diag;
    input>>diag;
    getline(input,line);
    if (!input) break;
    ocaPhi.push_back(diag);
  }
  
  PrintParsedData(log);
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
    for (int j=0; j<baths; j++) stream<<setw(12)<<Ms[i][j];
    stream<<setw(4)<<Mtot[i]<<setw(4)<<deg[i]<<setw(12)<<sJc[i]<<" ";
    
    for (int b=0; b<baths; b++) stream<<setw(2)<<sncab[i][b].size()<<" ";
    for (int b=0; b<baths; b++) stream<<setw(2)<<sncaf[i][b].size()<<" ";
    
    for (int b=0; b<baths; b++)
      for (map<int,double>::const_iterator l=sncab[i][b].begin(); l!=sncab[i][b].end(); l++)
	stream<<right<<setw(11)<<l->second<<" x "<<left<<setw(3)<<l->first<<right<<" ";

    for (int b=0; b<baths; b++)
      for (map<int,double>::const_iterator l=sncaf[i][b].begin(); l!=sncaf[i][b].end(); l++)
	stream<<right<<setw(11)<<l->second<<" x "<<left<<setw(3)<<l->first<<right<<" ";

    stream<<endl;
  }
  stream<<"# OCA diagrams: "<<ocaPhi.size()<<endl;
  for (int i=0; i<ocaPhi.size(); i++) stream<<ocaPhi[i]<<endl;
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


// double Auxiliary::CalcNCABubbleInside(const mesh1D& eps, const function1D<Hxy>& Hx)
// {
//   double sum=0;
//   for (int l=0; l<eps.size(); l++) sum += Hx[l].bh*eps.Dh(l);
//   return sum;
// }
// double Auxiliary::CalcNCABubble(int ip, const list<pair<int,int> >& hslist, map<pair<int,int>,int>& index_1,
// 				const mesh1D& eps, const function1D<Hxy>* Hx,
// 				const function2D<double>& pre_nca_x, function2D<int>& nca_xcopy)
// {
//   double sumb=0;
//   for (int ib=0; ib<baths; ib++){
//     pair<int,int> p0(nca_xcopy[ip][ib],ib);
//     if (nca_xcopy[ip][ib]>=0 && find(hslist.begin(),hslist.end(),p0)!=hslist.end()){
// #ifdef HEAVY_DEBUG    
//       clog<<"3 Calculating bubble for "<<ip<<" "<<nca_xcopy[ip][ib]<<" with index "<<index_1[p0]<<endl;
// #endif	
//       nca_xcopy[ip][ib]=-1;
//       int j=index_1[p0];
//       sumb += CalcNCABubbleInside(eps,Hx[j]) * pre_nca_x(ip,ib);
//     }
//   }
//   return sumb;
// }

void Auxiliary::Calc_NCA_Sigmab(const mesh1D& omd)
{
  for (int b=0; b<baths; b++){
    GtA1.Product(Gm,aAc[b],0,mpos);
    GtA2.Product(Gp,aAc[b],mpos,aAc[b].size_N());

    for (int j=0; j<Na; j++){
      for (map<int,double>::const_iterator l=common::sncab[j][b].begin(); l!=common::sncab[j][b].end(); l++){
	int ind = l->first;
	if (ind>=0 && ind<Na){
	  double prf = l->second/static_cast<double>(common::deg[j]);
	  for (int i=0; i<mpos; i++)         Sigtn_(i,j) += prf * GtA1(ind,i)/fe[i];
	  for (int i=mpos; i<om.size(); i++) Sigtn_(i,j) += prf * GtA2(ind,i-mpos)/(1-fe[i]);
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

void Auxiliary::Calc_NCA_Sigmaf(const mesh1D& omd)
{
  for (int b=0; b<baths; b++){
    GtA1.Product(Gm,aAc[baths+b],0,mpos);
    GtA2.Product(Gp,aAc[baths+b],mpos,aAc[baths+b].size_N());

    for (int j=0; j<Na; j++){
      for (map<int,double>::const_iterator l=common::sncaf[j][b].begin(); l!=common::sncaf[j][b].end(); l++){
	int ind = l->first;
	if (ind>=0 && ind<Na){
	  double prf = l->second/static_cast<double>(common::deg[j]);
	  for (int i=0; i<mpos; i++)         Sigtn_(i,j) += prf * GtA1(ind,i)/fe[i];
	  for (int i=mpos; i<om.size(); i++) Sigtn_(i,j) += prf * GtA2(ind,i-mpos)/(1-fe[i]);
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

void Auxiliary::MakeLsList(const deque<OCAd>& oca_diag, list<pair<int,int> >& lslist_fb, list<pair<int,int> >& lslist_af)
{
  for (deque<OCAd>::const_iterator il=oca_diag.begin(); il!=oca_diag.end(); il++){
    pair<int,int> p0(il->states[1],il->states[0]);
    pair<int,int> p2(il->states[2],il->states[3]);
    if (find(lslist_fb.begin(),lslist_fb.end(), p0) == lslist_fb.end()) lslist_fb.push_back(p0);
    if (find(lslist_af.begin(),lslist_af.end(), p2) == lslist_af.end()) lslist_af.push_back(p2);
    if (il->states[1]!=il->states[3]){
      pair<int,int> p0(il->states[3],il->states[0]);
      pair<int,int> p2(il->states[2],il->states[1]);
      if (find(lslist_fb.begin(),lslist_fb.end(), p0) == lslist_fb.end()) lslist_fb.push_back(p0);
      if (find(lslist_af.begin(),lslist_af.end(), p2) == lslist_af.end()) lslist_af.push_back(p2);
    }
  }
}
void Auxiliary::MakeCombFunc(double omega, const list<pair<int,int> >& lslist,
			    map<pair<int,int>,int>& index_1, mesh1D& eps, function1D<Lxy>* Lx)
{
  int j=0;
  for (list<pair<int,int> >::const_iterator i=lslist.begin(); i!=lslist.end(); i++,j++){
    index_1[*i]=j;
    int pn1 = i->first;
    int pn2 = i->second;
    funProxy<Gmpr1>* G1 = reinterpret_cast<funProxy<Gmpr1>* >(&(G[pn1]));
    funProxy<Gmpr2>* G2 = reinterpret_cast<funProxy<Gmpr2>* >(&(G[pn2]));
    MakeCombinedFunctions(eps, Lx[j], om, *G1, 0.0, om, *G2, omega);
#ifdef HEAVY_DEBUG    
    clog<<"1 Making combined functions for "<<pn1<<" "<<pn2<<" at position "<<j<<endl;
#endif
  }
}


double CalcA00_vertex(double omega, const mesh1D& om, const mesh1D& momd, const function2D<Amp>& mAc,
		      const mesh1D& eps1, const mesh1D& eps2,
		      map<pair<int,int>,int>& indexfb_1, map<pair<int,int>,int>& indexaf_1,
		      function1D<Lxy>* L1, function1D<Lxy>* L2,
		      int i0, int i1, int i2, int i3, int b2)
{
//          / | \
//      i1 /  |  \ i2
//        /  /|\  \
//       /    |    \
//  b1 - eps2 |b2 eps1- b1
//       \    |    /
//        \   |   /
//      i0 \  |  / i3
//          \ | /
//
//       L2        L1
  double sum=0;
  for (int i=0; i<eps1.size(); i++){
    double x[3] = {0.0, omega, eps1[i]};
    double peak[3] = {om.dcenter()+x[0], om.dcenter()+x[1], momd.dcenter()+x[2]};
    int j = indexfb_1[make_pair(i1,i0)];
    Lxy Lr = Integrate2c1(eps2, L2[j], momd, mAc[b2], x[0], x[1], x[2], peak[0], peak[1], peak[2]);
    int l = indexaf_1[make_pair(i2,i3)];
    sum += Product(Lr, L1[l][i])*eps1.Dh(i);
  }
  return sum;
}

function1D<double> Auxiliary::Calc_OCA_A00(double omega, const mesh1D& momd, const function2D<Amp>& mAc)
{
  function1D<double> vertex(baths);
  for (int b=0; b<baths; b++) vertex[b]=0;
  
  list<pair<int,int> > lslist_fb, lslist_af; // list of int-pars that represent left (right) half-circle
  
  MakeLsList(common::ocaPhi, lslist_fb, lslist_af);
  
  mesh1D eps2, eps1;
  map<pair<int,int>,int> indexfb_1, indexaf_1;// index_1 gives index of pair in the above lists
  function1D<Lxy>* L2 = new function1D<Lxy>[lslist_fb.size()];//Lxy has 6 doubles that correspond to ({+-,-+},-R,+R,R-,R+,RR)
  function1D<Lxy>* L1 = new function1D<Lxy>[lslist_af.size()];
  MakeCombFunc(omega, lslist_fb, indexfb_1, eps2, L2);
  MakeCombFunc(omega, lslist_af, indexaf_1, eps1, L1);
  
  for (deque<OCAd>::const_iterator il=common::ocaPhi.begin(); il!=common::ocaPhi.end(); il++){
    
    int b = il->ib1;
    double v0 = CalcA00_vertex(omega, om, momd, mAc, eps1, eps2, indexfb_1, indexaf_1, L1, L2,
			       il->states[0], il->states[1], il->states[2], il->states[3], il->ib2);
    
    double fact = il->f/common::Ns[b];
    vertex[b] += fact*v0/(M_PI*M_PI*common::Q);

    if (il->states[1]==il->states[3] && il->ib1==il->ib2){
      vertex[b] += fact*v0/(M_PI*M_PI*common::Q);
    } else{
      b = il->ib2;
      v0 = CalcA00_vertex(omega, om, momd, mAc, eps1, eps2, indexfb_1, indexaf_1, L1, L2,
			  il->states[0], il->states[3], il->states[2], il->states[1], il->ib1);
      fact = il->f/common::Ns[b];
      
      vertex[b] += fact*v0/(M_PI*M_PI*common::Q);
    }
    
  }
  
  delete[] L2;
  delete[] L1;
  
  //  static ofstream outp("cvertex.txt"); outp.precision(16);
  //  outp<<omega<<" ";
  //  for (int b=0; b<baths; b++) outp<<vertex[b]<<" ";
  //  outp<<endl;

  return vertex;
}

void Physical::MissingDoping(double start, ofstream& log)
{
  log<<"Missing doping : ";
  for (int b=0; b<baths; b++){
    double sum = 0;
    for (int i=0; i<omd.size(); i++) {
      if (omd[i]>start) sum += G00[b][i].imag()*fed[i]*omd.Dh(i);
    }
    sum *= -common::Ns[b]/M_PI;
    common::miss_nd[b] = common::nalpha[b]-sum;
    log<<b<<" : "<<common::miss_nd[b]<<" ";
  }
  log<<endl;

  common::Sinfty.resize(baths);
  for (int b=0; b<baths; b++){
    double sum0 = 0, sum1 = 0;
    for (int i=0; i<omd.size(); i++) {
      sum0 += A00(i,b)*omd.Dh(i);
      sum1 += A00(i,b)*omd[i]*omd.Dh(i);
    }
    common::moment[b][0] = sum0;
    common::moment[b][1] = sum1;
    common::Sinfty[b] = sum1/sum0-common::Ed[b];
  }
}

void Auxiliary::PrintSign()
{
  ofstream outp("Sign.dat");
  for (int iom=0; iom<om.size(); iom++){
    outp<<om[iom]<<" ";
    for (int i=0; i<Na; i++) outp<<-Sigtn_[iom][i]<<" ";
    outp<<endl;
  }
}

void Physical::CalculateProducts(double u, double fu, const mesh1D& om, const function2D<double>& Gm)
{
  apar ap;
  cintpar pi;
  tint position = om.InitInterpLeft();
  InterpLeft(om[0]-u, om, position, pi);
  for (int i=0; i<Na; i++) if (Pexists[i]) aF[i].InterpolateFirst(pi);
  InterpLeft(om[1]-u, om, position, pi);
  ap.SetUpCsFirst(u, om);
  for (int i=0; i<Na; i++) if (Pexists[i]) Gtx(i,0) = aF[i].InterpolateNext(pi, ap) * om.Dh(0);
  for (int j=1; j<om.size()-1; j++){
    InterpLeft(om[j+1]-u, om, position, pi);
    ap.SetUpCs(u, j, om, om.Dh(pi.i+1));
    for (int i=0; i<Na; i++) if (Pexists[i]) Gtx(i,j) = aF[i].InterpolateNext(pi, ap) * om.Dh(j);
  }
  ap.SetUpCsLast(u, om);
  for (int i=0; i<Na; i++) if (Pexists[i]) Gtx(i,om.size()-1) = aF[i].InterpolateLast(ap) * om.Dh(om.size()-1);

  Cmp.resize(Na,Na);
  for (int i=0; i<Na; i++){
    for (int b=0; b<baths; b++){
      for (map<int,double>::const_iterator l=common::sncab[i][b].begin(); l!=common::sncab[i][b].end(); l++){
	int ind = l->first;
	if (ind>=0 && ind<Na) Cmp(i,ind) = product(Gtx[i].MemPt(),Gm[ind].MemPt(),om.size())/fu;
      }
    }
  }
}

void Physical::Calculate_NCA_A00(const mesh1D& omega, const function2D<double>& Gp, const function2D<double>& Gm,
				 const function1D<double>& Energy)
{
  int m = omd.find_(0.0)+1;

  Gtx.resize(Na, omega.size());
  
  for (int i=0; i<Na; i++) if (Pexists[i]) aF[i].SetUp(Gp[i],omega);
  for (int i=0; i<m; i++){
    CalculateProducts(omd[i], fed[i], omega, Gm);
    for (int b=0; b<baths; b++){
      double sum=0;
      for (int j=0; j<Na; j++)
	for (map<int,double>::const_iterator l=common::sncab[j][b].begin(); l!=common::sncab[j][b].end(); l++){
	  int ind = l->first;
	  double prf = l->second/common::Ns[b];
	  if (ind>=0 && ind<Na) sum += prf*Cmp(j,ind);
	}
      A00(i,b) = sum/(M_PI*M_PI*common::Q);
    }
  }

  for (int i=0; i<Na; i++) if (Pexists[i]) aF[i].SetUp(Gm[i],omega);
  for (int i=m; i<omd.size(); i++){
    CalculateProducts(omd[i], (1-fed[i]), omega, Gp);
    for (int b=0; b<baths; b++){
      double sum=0;
      for (int j=0; j<Na; j++)
	for (map<int,double>::const_iterator l=common::sncab[j][b].begin(); l!=common::sncab[j][b].end(); l++){
	  int ind = l->first;
	  double prf = l->second/common::Ns[b];
	  if (ind>=0 && ind<Na)	sum += prf*Cmp(j,ind);
	}
      A00(i,b) = sum/(M_PI*M_PI*common::Q);
    }
  }
}

void Physical::AddCore(const mesh1D& omega, const function2D<double>& Gp, const function2D<double>& Gm,
		       const function1D<double>& Energy, ofstream& log)
{
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
	A00c(i,b) = sum1+sum2;
      }
    }
    // Checking doping!
    for (int b=0; b<baths; b++){
      double suma = 0, sumc = 0;
      for (int i=0; i<omd.size(); i++) {
	suma += A00(i,b)*fed[i]*omd.Dh(i);
	sumc += A00c(i,b)*fed[i]*omd.Dh(i);
      }
      suma *= common::Ns[b];
      sumc *= common::Ns[b];
      double miss_nd = common::nalpha[b]-(suma+sumc);
      double core_fact = 1.;
      if (sumc!=0 && common::renorm_core){
	core_fact = (common::nalpha[b]-suma)/sumc;
	if (core_fact<0) core_fact=0;
	if (core_fact>10) core_fact = 10;
	log<<b<<" : "<<miss_nd<<" renormaliziang core part by "<<core_fact<<endl;
      }
      for (int i=0; i<omd.size(); i++) A00(i,b) += A00c(i,b)*core_fact;

      if (common::renorm){
	double suml=0, sumr=0;
	for (int i=0; i<omd.size(); i++){
	  suml += A00(i,b)*fed[i]*omd.Dh(i);
	  sumr += A00(i,b)*(1-fed[i])*omd.Dh(i);
	}
	int izero = omd.find_(0.0);
	double ml1=0, mr1=0;
	for (int i=0; i<izero; i++) {
	  ml1 += omd[i]*A00(i,b)*fed[i]*omd.Dh(i);
	  mr1 += omd[i]*A00(i,b)*(1-fed[i])*omd.Dh(i);
	}
	double ml2=0, mr2=0;
	for (int i=izero+1; i<omd.size(); i++) {
	  ml2 += omd[i]*A00(i,b)*fed[i]*omd.Dh(i);
	  mr2 += omd[i]*A00(i,b)*(1-fed[i])*omd.Dh(i);
	}
	double n0 = common::nalpha[b]/common::Ns[b];
	double C = (-ml2 + ml2*n0 + mr2*n0 - mr2*suml + ml2*sumr)/(ml1*mr2-ml2*mr1);
	double D = (ml1 - ml1*n0 - mr1*n0 + mr1*suml - ml1*sumr)/(ml1*mr2-ml2*mr1);

	if (1+C*omd[0]<0) C = -1/omd[0];
	if (1+D*omd.last()<0)  D = -1/omd.last();
	
	for (int i=0; i<izero; i++) A00(i,b) *= (1+C*omd[i]);
	for (int i=izero+1; i<omd.size(); i++) A00(i,b) *= (1+D*omd[i]);
      
	log<<"Renormalizing A["<<b<<"] by "<<C<<", "<<D<<"at negative and positive frequency"<<endl;
      }
    }
  }
  
  ofstream out("Aloc.imp"); out.precision(16);
  for (int i=0; i<omd.size(); i++){
    out<<setw(25)<<omd[i]<<" ";
    for (int b=0; b<baths; b++) out<<setw(25)<<A00(i,b)<<" ";
    out<<endl;
  }
}

void Physical::ReadAlocImp(const string& filename)
{
  ifstream inp(filename.c_str());
  inp.ignore(10000,'\n'); // comment line
  for (int i=0; i<omd.size(); i++){
    double tomd;
    inp>>tomd;
    if (fabs(tomd-omd[i])>1e-5 || !inp){
      cerr<<"Impurity spectral function "<<filename<<" is bad or not compatible with the bath mesh"<<endl; exit(1);
    }
    for (int b=0; b<baths; b++) inp>>A00[i][b];
    inp.ignore(10000,'\n'); // the rest of the line
  }
}

void RememberParams (int argc, char *argv[]){
  ofstream param ((common::outdir+"/history.nca").c_str(), ios::app);
  if (!param) cerr<<" Didn't suceeded to open params file!"<<(common::outdir+"/history.nca")<<endl;
  for (int i=0; i<argc; i++) param << argv[i] << " ";
  param << endl;
}

#endif
