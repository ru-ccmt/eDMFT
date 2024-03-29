// @Copyright 2007 Kristjan Haule
// 
#ifdef _MPI
#include <mpi.h>
#endif
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>
#include "parser.h"
#include "timer.h"
#include "util.h"
#include "complex.h"
#include "mesh.h"
#include "blas.h"
#include "function.h"
#include "classes.h"
#include "integrate.h"
#include "Common.h"
#include "cmpi.h"

using namespace std;


namespace Parameters{
  StringPar   Sig        ("Sig",     "Sigma.000", "\t\t# The name of the input Self-energies");
  StringPar   Ac         ("Ac",    "None",  "\t\t# The name of the input bath function if only the spectral function is given");
  StringPar   Delta      ("Delta", "None", "\t\t# The name of the input bath function if a retarded Delta is given");
  StringPar   out        ("out",           ".", "\t\t# The name of the output directory");
  StringPar   gloc       ("gloc",   "gloc.out", "\t# The name of the output Green's function");
  StringPar   sig        ("sig",     "sig.out", "\t\t# The name of the output Self-energy");
  StringPar   SigOut     ("SigOut","Sigma.000", "\t# The name of the output pseudo-particle Self-energy");
  StringPar   AlocOut    ("AlocOut", "Aloc.imp", "\t# The name of the output Physical Spectral function");
  StringPar   cix        ("cix",     "cix.dat", "\t\t# The name of the input file containing information about bands and their degeneracy");
  StringPar   fin        ("fin",     ".fin.0",  "\t\t# The name of the file marking finish of the calculation");
  Par<double> Ed         ("Ed",         -2., "\t\t# Energy level");
  Par<double> U          ("U",           4., "\t\t\t# Coulomb repulsion");
  Par<double> T          ("T",          0.2, "\t\t# Temperature ");
  Par<double> Q0         ("Q",            1, "\t\t\t# Default average Q in grand-canonical ansamble");
  Par<int>    prt        ("prt",          1, "\t\t# Weather to print intermediate results");
  Par<int>    prtSF      ("prtSF",        1, "\t\t# Weather to print out Spectral functions of auxiliary particles");
  Par<double> alpha      ("alpha",      0.5, "\t\t# The fraction of the new self energy to be used in the next iteration");
  Par<double> max_diff   ("max_diff",  1e-6, "\t# Difference in auxiliary functions");
  Par<int>    max_steps  ("max_steps",    1, "\t\t# Number of steps in auxiliary loop");
  Par<double> StartLambda("StartLambda", -1, "\t# Where to start searching for the lambda0");
  Par<double> EndLambda  ("EndLambda",    1, "\t\t# Where to stop searching for the lambda0");
  Par<double> dLambda    ("dLambda",    0.1, "\t\t# Step in searching for the lambda");
  Par<int>    followPeak ("followPeak",  -3, "\t# Wheather to determin zero frequency from the diverging pseudo-particle (-1=none, 0=b, 1=f, 2=a)");
  Par<double> MissDopSt  ("MissDopSt", -100, "\t# Missing doping due to projection and finite mesh starting at MissDopSt");
  Par<int>    N_ac       ("N_ac",        50, "\t\t# Number of points to add for calculating Delta to be used in core calculation");
  Par<double> dom_ac     ("dom_ac",      0.5,"\t\t# Distance between points that are added to the mesh for core calculation");
  Par<int>    acore      ("acore",         0,"\t\t# Weather to add core contribution to auxiliary self-energies");
  Par<int>    pcore      ("pcore",         1,"\t\t# Weather to add core contribution to physical spectral function");
  Par<int>    renorm_core("renorm_core",  0, "\t# Renormalization of core");
  Par<int>    renorm     ("renorm",       1, "\t\t# Renormalization of spectral function");
  Par<int>    OnlyCore   ("Ocore",        0, "\t\t\t# Only core is added to physical spectral function. It must already exist.");
}

int common::Na, common::Nc;
double common::U;
double common::T;
int common::baths;
function1D<double> common::Ed;
function1D<int> common::Ns;
function2D<double> common::Ms;
function1D<int> common::Mtot;
function1D<int> common::deg;
function1D<string> common::Eds;
function1D<double> common::nalpha;
function1D<double> common::miss_nd;
double common::beta;
double common::lambda0;
double common::Q;
double common::Q0;
double common::nd;
string common::outdir;
int common::totDeg;
function1D<double> common::sJc;
vector<vector<map<int,double> > > common::sncab;	   // index for hole diagrams
vector<vector<map<int,double> > > common::sncaf;	   // index for particle diagrams
deque<OCAd> common::ocaPhi;
int common::N_ac;
double common::dom_ac;
bool common::renorm_core, common::renorm;
int common::acore, common::pcore;
function1D<double> common::Sinfty;
function2D<double> common::moment;
double common::Fimp;
double common::Epot;
double common::TrLogGimp;

/********* Only for testing ****************/
inline void AddProduct(const double f0, const double f1, const double f2, double& sum, double dh)
{
 sum += f0*f1*f2*dh;
}
inline void AddProduct(const double f0, const double f1, double& sum, double dh)
{
  sum += f0*f1*dh;
}
inline void Setv(double& fi, const functionb<double>& f, int i)
{
  fi = f[i];
}
inline void Seti(double& fi, const functionb<double>& f, const intpar& ip)
{
  if (ip.p<0 || ip.p>1) { fi = 0; return;}
  fi = f[ip.i] + ip.p*(f[ip.i+1] - f[ip.i]);
}
inline void Setv(double& f12, double f1, double f2)
{
  f12 = f1*f2;
}
/******* Only for testing ******************/

ostream& operator<<(ostream& stream, const Amp& a){
  stream<<a.m<<" "<<a.p;
  return stream;
};


void DivideFile(double pa, int nom, int my_rank, int mpi_size, function1D<int>& strFreq, function1D<int>& endFreq)
{
  double S0 = 2/(1+pa)*nom/mpi_size;
  double dS = 4*(1-pa)/(1+pa)*nom/(mpi_size*mpi_size);
  
  int st=0;
  for (int j=0; j<mpi_size; j++){
    int ssize = static_cast<int>(S0 - dS*fabs((mpi_size-1.)/2.-j) + 0.5);
    strFreq[j] = st;
    st += ssize;
    if (st>nom) st=nom;
    endFreq[j] = st;
  }
}

int main(int argc, char *argv[])
{
  setvbuf (stdout, NULL, _IONBF, BUFSIZ);
  /*
  MPI::Init(argc, argv);
  int my_rank = MPI::COMM_WORLD.Get_rank();
  int mpi_size = MPI::COMM_WORLD.Get_size();
  */
  int my_rank, mpi_size;
  MPI_Init(argc, argv, my_rank, mpi_size);
  int Master = 0;

  using namespace Parameters;
  DynArray arr(100, &Sig, &Ac, &Delta, &cix, &out, &gloc, &sig, &SigOut, &AlocOut, &fin, &Ed, &U, &T, &Q0, &alpha, &max_diff, &max_steps,
	       &StartLambda, &EndLambda, &dLambda, &prt, &prtSF, &followPeak,
	       &MissDopSt, &N_ac, &dom_ac, &acore, &pcore,
	       &renorm_core, &renorm, &OnlyCore, NULL);
  
  int targc=1;
  char* targv[2] = {"", "PARAMS.oca"};
  ifstream inpfile(targv[1]);
  if (!inpfile) {cerr<<"Could not find input file "<<targv[1]<<"! Exiting...."<<endl;  MPI_Finish(); exit(1);}
  arr.ParsComLine(targc, targv);
  
  ofstream log(NameOfFile("oca_log", my_rank).c_str());
  
  arr.print (log);
  
  common::ParsInputFile(cix, log);
  common::SetParameters(Ed,U,T,Q0,out,N_ac,dom_ac,acore,pcore,renorm_core,renorm);

  Physical ph(common::Na, common::Nc, common::baths);
  Auxiliary au(common::Na, common::Nc, common::baths);

  if (static_cast<string>(Ac) != "None"){
    if (!ph.ReadBathFunction(Ac, log, true)) {MPI_Finish(); exit(1);}
  }else if (static_cast<string>(Delta) != "None"){
    if (!ph.ReadBathFunction(Delta, log, false)) {MPI_Finish(); exit(1);}
  }else{
    cerr<<"The bath input is not set! Please give input filename for Ac or Delta!"<<endl;
    exit(1);
  }
  if (!au.ReadSelfEnergy(Sig,Ed,T,U,ph.omega(),ph.Ac0(),log)) { MPI_Finish(); exit(1);}

//   au.CreateSigma000(ph.omega(), ph.Ac0());
//   au.Print(999);
//   return 0;


  
  StartLambda = StartLambda+au.minEnergy;
  EndLambda = EndLambda+au.minEnergy;
  
  if (followPeak==-3){
    StartLambda = StartLambda+au.minEnergy;
    EndLambda = EndLambda+au.minEnergy;
    common::Q0 = common::totDeg;
    followPeak=-1;
  }

  double pa=0.5;
  double pp=0.9;
  
  function1D<int> StrFreq(mpi_size), EndFreq(mpi_size);
  DivideFile(pa, au.omega().size(), my_rank, mpi_size, StrFreq, EndFreq);
  log<<"Frequencies divided in the following way: "<<endl;
  for (int i=0; i<mpi_size; i++) log<<i<<" "<<StrFreq[i]<<" - "<<EndFreq[i]<<endl;
  
  function1D<int> send_size(mpi_size);
  function1D<int> start(mpi_size);
  for (int i=0; i<mpi_size; i++){
    send_size[i] = (EndFreq[i]-StrFreq[i])*common::Na;
    start[i] = StrFreq[i]*common::Na;
  }


  Timer t; t.start();
  if (max_steps!=0) au.SetUpAverageAc(ph.omega(), ph.momega(), ph.Ac0(), ph.fe());
  int n=0; double diff=1;
  while (diff>max_diff && n++<max_steps){
    
    log<<"  "<<n<<".) KramarsKronig"<<endl;
    au.KramarsKronig();
    log<<"  "<<n<<".) DeterminSpectralFunctions ";    
    common::nd = au.DeterminSpectralFunctions(StartLambda,EndLambda,dLambda,followPeak, log);
    log<<"  "<<n<<".) nd:  "<<common::nd<<endl; 
    au.PrintNorm(log);
    if (prt && my_rank==Master) au.Print(n);

    log<<"  "<<n<<".) SelfConsistentStep"<<endl;
    au.SetSignToZero();

    au.Calc_NCA_Sigmab(ph.omega());
    au.Calc_NCA_Sigmaf(ph.omega());
    
    for (int i=StrFreq[my_rank]; i<EndFreq[my_rank]; i++){
      log<<i<<" ";
      au.Calc_OCA_SelfEnergy(i, ph.omega(),ph.pA(),ph.momega(),ph.mA());
      log.flush();
    }
    log<<endl;

    // The new self-energy in au.Sign() is communicated to all machines
    //MPI::COMM_WORLD.Allgatherv(MPI_IN_PLACE, send_size[my_rank], MPI_DOUBLE, au.Sign(), send_size.MemPt(), start.MemPt(), MPI_DOUBLE);
    MPI_Allgather_Self_energy(my_rank, send_size, start, au.Sign());

    
    if (prt && my_rank==Master) au.PrintSign();

    log<<"  "<<n<<".) Difference between steps: \t\t"<<COLOR(GREEN,(diff=au.Difference()))<<endl;
    au.DeterminSelfEnergies(Parameters::alpha);
  }

  t.stop();
  log<<"Time needed: "<<t.toSeconds()<<" average time per frequency: "<<t.toSeconds()/(EndFreq[my_rank]-StrFreq[my_rank])<<endl;
  
  t.start();

  // This happens only if max_steps==0 -> Only Physical spectral function is calculated and
  // no self-consistency step. In this case one needs to calculate (for the first time) auxiliary Spectral functions.
  if (max_steps==0){
    log<<"  "<<n<<".) KramarsKronig"<<endl;
    au.KramarsKronig();
    log<<"  "<<n<<".) DeterminSpectralFunctions ";
    
    common::nd = au.DeterminSpectralFunctions(StartLambda,EndLambda,dLambda,followPeak,log);
    log<<"  "<<n<<".) nd:  "<<common::nd<<endl; 
    au.PrintNorm(log);
    if (prt && my_rank==Master) au.Print(n);
  }

  function1D<int> dStrFreq(mpi_size), dEndFreq(mpi_size);
  DivideFile(pp, ph.omega().size(), my_rank, mpi_size, dStrFreq, dEndFreq);
  log<<"Frequencies divided in the following way: "<<endl;
  for (int i=0; i<mpi_size; i++) log<<i<<" "<<dStrFreq[i]<<" - "<<dEndFreq[i]<<endl;
  for (int i=0; i<mpi_size; i++){
    send_size[i] = (dEndFreq[i]-dStrFreq[i])*common::baths;
    start[i] = dStrFreq[i]*common::baths;
  }
  
  if (OnlyCore){
    ph.ReadAlocImp(AlocOut);
  } else{
    
    ph.Calculate_NCA_A00(au.omega(), au._Gp(), au._Gm(), au.Energ());
    if (my_rank==Master) ph.PrintA00();
    
    for (int i=dStrFreq[my_rank]; i<dEndFreq[my_rank]; i++){
      log<<i<<" ";
      ph.SetA00(i,au.Calc_OCA_A00(ph.omega()[i], ph.momega(), ph.mA()));
      log.flush();
    }
    log<<endl;

    // The new A00 is collected on root
    //MPI::COMM_WORLD.Gatherv((my_rank==Master) ? MPI_IN_PLACE : ph.P_A00()+start[my_rank], send_size[my_rank], MPI_DOUBLE, ph.P_A00(), send_size.MemPt(), start.MemPt(), MPI_DOUBLE, Master);
    MPI_Gather_A00(my_rank, Master, ph.P_A00(), send_size, start);
  }
  t.stop();
  
  log<<"Time needed: "<<t.toSeconds()<<" average time per frequency: "<<t.toSeconds()/(dEndFreq[my_rank]-dStrFreq[my_rank])<<endl;

  if (my_rank==Master){
    ph.AddCore(au.omega(), au._Gp(), au._Gm(), au.Energ(), log);
    
    ph.DeterminG00(log);
    ph.CalcSelfEnergy();

    ph.MissingDoping(MissDopSt, log);
    ofstream outg((common::outdir+"/"+static_cast<string>(gloc)).c_str());
    outg.precision(12);
    print(outg, ph.omd, ph.G00, 25);
    ofstream outs((common::outdir+"/"+static_cast<string>(sig)).c_str());
    outs.precision(12);
    print(outs, ph.omd, ph.Sig, 25);
    
    au.Print(common::outdir+"/"+static_cast<string>(SigOut).c_str(),prtSF);
    ph.Print(common::outdir+"/"+static_cast<string>(AlocOut).c_str());
  }
  
  //MPI::Finalize();
  MPI_Finish();
  return 0;
}
