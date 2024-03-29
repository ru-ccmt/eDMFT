// @Copyright 2007 Kristjan Haule
// 
#include <sys/stat.h>
#include <sys/types.h>
#include "parser.h"
#include "timer.h"
#include "util.h"
#include "complex.h"
#include "mesh.h"
#include "blas.h"
#include "function.h"
#include "Common.h"

using namespace std;

namespace Parameters{
  StringPar   Sig        ("Sig",  "sig.inp", "\t\t# The name of the input Self-energies");
  StringPar   Ac         ("Ac",    "None",  "\t\t# The name of the input bath function if only the spectral function is given");
  StringPar   Delta      ("Delta", "None", "\t\t# The name of the input bath function if a retarded Delta is given");
  StringPar   out        ("out",        ".", "\t\t# The name of the output directory");
  StringPar   gloc       ("gloc","gloc.out", "\t# The name of the output Green's function");
  StringPar   sig        ("sig",  "sig.out", "\t\t# The name of the output Self-energy");
  StringPar   susc       ("susc","Susc.out", "\t\t# The name of the output susceptibility");
  StringPar   cix        ("cix",  "cix.dat", "\t\t# The name of the input file containing information about bands and their degeneracy");
  Par<double> Ed         ("Ed",         -2., "\t\t# Energy level");
  Par<double> U          ("U",           4., "\t\t\t# Coulomb repulsion");
  Par<double> T          ("T",          0.2, "\t\t# Temperature ");
  Par<double> Q0         ("Q",            1, "\t\t\t# Default average Q in grand-canonical ansamble");
  Par<int>    prt        ("prt",          1, "\t\t# Weather to print intermediate results");
  Par<double> alpha      ("alpha",      0.5, "\t\t# The fraction of the new self energy to be used in the next iteration");
  Par<double> max_diff   ("max_diff",  1e-6, "\t# Criterium to finish the procedure");
  Par<int>    max_steps  ("max_steps",  300, "\t# Criterium to finish the procedure");
  Par<double> StartLambda("StartLambda", -1, "\t# Where to start searching for the lambda0");
  Par<double> EndLambda  ("EndLambda",    1, "\t\t# Where to stop searching for the lambda0");
  Par<double> dLambda    ("dLambda",    0.1, "\t\t# Step in searching for the lambda");
  Par<int>    followPeak ("followPeak",  -3, "\t# Wheather to determin zero frequency from the diverging pseudo-particle (-2: lamdba0==StartLambda, -1: Q==Q0, 0: follow b, 1: follow f, 2: foolow a)");
  Par<double> MissDopSt  ("MissDopSt", -100, "\t# Missing doping due to projection and finite mesh starting at MissDopSt");
  Par<int>    N_ac       ("N_ac",        50, "\t\t# Number of points to add for calculatin Delta to be used in core calculation");
  Par<double> dom_ac     ("dom_ac",      0.5,"\t\t# Distance between points that are added to the mesh for core calculation");
  Par<int>    acore      ("acore",         1,"\t\t# Weather to add core contribution to auxiliary self-energies");
  Par<int>    pcore      ("pcore",         1,"\t\t# Weather to add core contribution to physical spectral function");
  Par<int>    lorentz    ("lorentz",       1,"\t\t# Weather to subtract lorentz from diverging spectral functions and treat it analytically");
  Par<double> SrchLorentz("SearchLorentz", 2.5,"\t# How far from zero to search fro Lorentz");
  Par<double> LorentzMaxRatio("LorentzMaxRatio", 1,"\t# How far from zero to search for Lorentz");
  Par<int>    FirstLorentz("FirstLorentz", 0,"\t# First pseudoparticle which could be augmented with lorentz");
  Par<int>    LastLorentz("LastLorentz", 10000,"\t# Last pseudoparticle which could be augmented with lorentz");
  Par<int>    CmpDiff     ("CmpDiff",-1, "\t\t# When calculating Difference, only first CmpDiff particles should be taken into account (-1,all)");
  Par<int>    renorm_core("renorm_core",  0, "\t# Renormalization of core");
  Par<int>    renorm     ("renorm",       1, "\t\t# Renormalization of spectral function");
}

int common::Na;
int common::Nc;
double common::U;
//double common::J;
double common::T;
int common::baths;
function1D<double> common::Ed;
function1D<int> common::Ns;
function2D<double> common::Ms;
function1D<int> common::Mtot;
function1D<int> common::deg;
function2D<int> common::ncab;	  // index for hole diagrams
function2D<int> common::ncaf;	  // index for particle diagrams
function2D<double> common::prefactb; // prefactor for hole digrams
function2D<double> common::prefactf; // prefactor for particle diagrams
function2D<double> common::prefactG; // prefactor to calculate local Green's function
function1D<string> common::Eds;
function1D<double> common::nalpha;
function1D<double> common::miss_nd;
function1D<double> common::sJc;
function1D<double> common::Sinfty;
vector<vector<map<int,double> > > common::sncab;
vector<vector<map<int,double> > > common::sncaf;
vector<map<int,double> > common::suscb;
double common::beta;
double common::lambda0;
double common::Q;
double common::Q0;
double common::nd;
string common::outdir;
int common::totDeg;
int common::N_ac;
double common::dom_ac;
int common::acore, common::pcore;
bool common::SubtractLorentz;
double common::dlmin;
double common::LorentzMaxRatio;
double common::SearchLorentz;
bool common::renorm_core, common::renorm;
int common::FirstLorentz;
int common::LastLorentz;
bool common::cmp_susc;
function2D<double> common::moment;
double common::Fimp;
double common::Epot;
double common::TrLogGimp;

int main (int argc, char *argv[], char *env[]){
  using namespace Parameters;
  setvbuf (stdout, NULL, _IONBF, BUFSIZ);
  DynArray arr(90, &Sig, &Ac, &Delta, &cix, &out, &gloc, &sig, &susc, &Ed, &U, &T, &Q0, &alpha, &max_diff, &max_steps,
	       &StartLambda, &EndLambda, &dLambda,  &followPeak, &prt, &MissDopSt, &N_ac, &dom_ac,
	       &acore, &pcore, &lorentz, &SrchLorentz, &LorentzMaxRatio,&FirstLorentz,&LastLorentz,
	       &CmpDiff,&renorm_core,&renorm,NULL);
  
  if (argc<2) {
    arr.printl(clog);
    return 0;
  }
  arr.ParsComLine(argc, argv);
  arr.print (clog);

  common::ParsInputFile(cix);
  common::cmp_susc=false;

  common::SetParameters(Ed,U,T,Q0,out,N_ac,dom_ac,acore,pcore,lorentz,SrchLorentz,LorentzMaxRatio,
			FirstLorentz,LastLorentz,renorm_core,renorm);
  RememberParams (argc, argv);
  Physical ph(common::Na, common::Nc, common::baths);
  Auxiliary au(common::Na, common::Nc, common::baths);

  if (static_cast<string>(Ac) != "None"){
    if (!ph.ReadBathFunction(Ac, true)) exit(1);
  }else if (static_cast<string>(Delta) != "None"){
    if (!ph.ReadBathFunction(Delta, false)) exit(1);
  }else{
    cerr<<"The bath input is not set! Please give input filename for Ac or Delta!"<<endl;
    exit(1);
  }
  
  if (!au.ReadSelfEnergy(Sig,Ed,T,U,ph.omega(), ph.Ac0())) exit(1);

  //  au.CreateSigma000(ph.omega(), ph.Ac0());
  //  au.Print(999);
  //  return 0;

  
  au.SetUpAverageAc(ph.omega(), ph.momega(), ph.Ac0(), ph.fe());

  StartLambda = StartLambda+au.minEnergy;
  EndLambda = EndLambda+au.minEnergy;
  if (followPeak==-3){
    common::Q0 = common::totDeg;
    followPeak=-1;
  }

  int n=0; double diff=1;
  while (diff>max_diff && n++<max_steps){
    clog<<"  "<<n<<".) DeterminSpectralFunctions ";
    common::nd = au.DeterminSpectralFunctions(StartLambda,EndLambda,dLambda,followPeak);
    clog<<"  "<<n<<".) nd:  "<<common::nd<<endl; 
    au.PrintNorm(clog);
    if (prt) au.Print(n);

    clog<<"  "<<n<<".) SelfConsistentStep"<<endl;
    au.SetSignToZero();
    au.CalcSigmab(ph.omega());
    au.CalcSigmaf(ph.omega());
    
    diff = au.DeterminSelfEnergies(Parameters::alpha,CmpDiff);
    clog<<"  "<<n<<".) Difference between steps: \t\t"<<COLOR(GREEN,diff)<<endl;

    au.PrintCore("cores.dat");
  }

  ph.CalculateA00(au.omega(), au._Gp(), au._Gm(), au.Energ(), au.Lorentzm(), au.Lorentzp());
  
  ph.DeterminG00(1.0,clog);
  //ph.KramarsKronig();
  ph.CalcSelfEnergy();

  ph.MissingDoping(MissDopSt);
  
  au.Print(0);
  ofstream outg((common::outdir+"/"+static_cast<string>(gloc)).c_str());
  print(outg, ph.omd, ph.G00, 25);
  ofstream outs((common::outdir+"/"+static_cast<string>(sig)).c_str());
  print(outs, ph.omd, ph.Sig, 25);
  if (common::cmp_susc){
    ofstream outs((common::outdir+"/"+static_cast<string>(susc)).c_str());
    print(outs, ph.omd, ph.Chi, 25);
  }

  ph.Print0(common::outdir+"/Aloc.imp");

  return 0;
}
