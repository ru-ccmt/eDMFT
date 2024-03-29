// @Copyright 2007 Kristjan Haule
// 
#include <iostream>
#include <iomanip>
#include <fstream>
#include <list>
#include <vector>
#include <sstream>
#include "complex.h"
#include "mesh.h"
#include "function.h"
using namespace std;

class Lorentz
{
public:
  double x0, Gamma, Prefact;
  double Gamma2;
  Lorentz() : Prefact(0.0),Gamma(1),x0(0){};
  Lorentz(double x0_, double Gamma_, double Prefact_)
  { x0 = x0_; Gamma = Gamma_; Prefact = Prefact_; Gamma2 = Gamma*Gamma; }
  void Set(double x0_, double Gamma_, double Prefact_)
  { x0 = x0_; Gamma = Gamma_; Prefact = Prefact_; Gamma2 = Gamma*Gamma; }
  double operator() (double x) const {return Prefact/(sqr(x-x0)+Gamma2);}
  double realPart(double x) const {return Prefact*(x0-x)/(Gamma2+sqr(x-x0))/Gamma;}
  dcomplex analytic(const dcomplex& z){ return Prefact/Gamma/(x0-(z+dcomplex(0,Gamma)));}
  double integral() const { return M_PI*Prefact/Gamma;}
  //  l1(x-x1)**l2(x-x2-u)
  friend double convolution(const Lorentz& l1, const Lorentz& l2, double u);
  double g0(double a, double b) const
  {
    return 0.5*Prefact*log((a*a+Gamma2)/(b*b+Gamma2));
  }
  double g1(double a, double b) const
  {
    return Prefact*(atan(a/Gamma)-atan(b/Gamma))/Gamma;
  }
  string print()
  {
    stringstream str;
    str<<x0<<" "<<Gamma<<" "<<Prefact;
    return str.str();
  }
};

bool CheckStream(istream& inputf, int& n, int& m)
{
  istream input(inputf.rdbuf());
  input.seekg(0,ios::beg);
  
  string str; bool begincomm=false; n=0;
  getline(input,str); n++;
  if (!input.good()) {
    cerr << "ERROR: Wrong file format for hilbert or no data!" << endl;
    return false;
  }
  if (str.find('#')<string::npos) {
    begincomm=true;
    cout<<str<<endl;
    getline(input,str);
  };

  stringstream oneline(str);
  //  oneline << str << ends;
  m=0; double t;
  while (oneline){oneline>>t; m++;}
  m--;
  while (input){ getline(input,str); n++;}
  n--;
 
  clog << " Number of entries: "<< n <<endl;
  clog << " Number of columns: "<< m <<endl;
 
  inputf.seekg(0,ios::beg);
  if (begincomm) getline(inputf, str);
  if (!inputf){ cerr<<"Reopening didn't suceeded!"<<endl; return false;}
  return true;
}

bool ReadData(istream& input, mesh1D& om, function2D<double>& fi, int mdata, const vector<int>& cn)
{
  clog<<"Reading data with "<<om.size()<<" entries"<<endl;
  vector<double> data(mdata);
  int i=-1;
  while (input && ++i<om.size()){
    for (int j=0; j<mdata; j++) input>>data[j];
    input.ignore(500,'\n');
    for (int j=0; j<cn.size(); j++) fi[j][i] = data[cn[j]-1];
    om[i] = data[0];
  }
  return true;
}

string print_cn(const vector<int>& cn){
  stringstream str;
  for(int j=0; j<cn.size(); j++) str<<cn[j]<<" ";
  str<<ends;
  return str.str();
}

int main (int argc, char *argv[], char *env[])
{
  double mu=0.0;
  list<string> files;
  bool imagAx = false;
  bool rstdin = false;
  bool bosonic = false;
  double diom=0.1, maxS=100;
  double scale=1.0;
  vector<int> cn;
  cn.push_back(5);
  int i=0, Ni=0;
  while (++i<argc){
    string str(argv[i]);
    if (str=="-mu" && i<argc-1) mu = atof(argv[++i]);
    else if (str=="-Ni" && i<argc-1)  Ni = atoi(argv[++i]);
    else if (str=="-cn" && i<argc-1)  {
      stringstream cnum(argv[++i]);
      cn.clear();
      while(cnum){
	int it;
	cnum>>it;
	if (!cnum) break;
	cn.push_back(it);
      }
      clog<<"parsed columns: "<<print_cn(cn)<<endl;
    }
    else if (str=="-ms" && i<argc-1)  maxS = atof(argv[++i]);
    else if (str=="-d" && i<argc-1){
      string command;
      command = "perl -e'use Math::Trig; use Math::Complex; print " + string(argv[++i]) + "'>convert.temp";
      system(command.c_str());
      ifstream inp("convert.temp");
      inp>>diom;
      command = "rm -f convert.temp";
      system(command.c_str());
    } else if (str=="-s"){
      string str(argv[++i]);
      string command;
      command = "perl -e'use Math::Trig; use Math::Complex; print " + str + "'>convert.temp";
      system(command.c_str());
      ifstream inp("convert.temp");
      inp>>scale;
      command = "rm -f convert.temp";
      system(command.c_str());
    }
    else if (str=="-im")  imagAx = true;
    else if (str=="-si")  rstdin=true;
    else if (str=="-b")  bosonic=true;
    else {
      ifstream file(argv[i]);
      if (file) files.push_back(argv[i]);
    }
  }

  int n1, m1;
  fstream input;
  if (files.size()==1){
    input.open(files.begin()->c_str(), ios::in);
    if (!CheckStream(input,n1,m1)) exit(2);
  } else if(rstdin) {
    input.open("skrams.temp", ios::in | ios::out);
    input<<cin.rdbuf();
    input.seekg(0,ios::beg);
    if (!CheckStream(input,n1,m1)) exit(2);
  } else{
    clog<<"************* GENERAL KRAMARS-KRONIG ***************\n";
    clog<<"**                                                **\n";
    clog<<"**       Copyright Kristjan Haule,  9.2.2003      **\n";
    clog<<"****************************************************\n";
    clog<<"\n";
    clog<<"[... |] skrams file [-cn int,int,int...] [-Ni int -d double] [-si] [-ms double]\n" ;
    clog<<"Options:   file     Filename of the imaginary part of the function\n";
    clog<<"           -cn      Column number of data ("<<print_cn(cn)<<")\n";
    clog<<"           -Ni      Number of points on imaginary axes (if Ni=0, skrams works on real axes) ("<<Ni<<")\n";
    clog<<"           -d       Distance between points on imaginary axes (can contain any perl expresion) ("<<diom<<")\n";
    clog<<"           -si      Data will be read from standard input\n";
    clog<<"           -ms      Condition to substract Lorentz ("<<maxS<<")\n";
    clog<<"           -b       Bosonic (fermionic) matsubara frequencies! "<<bosonic<<"\n";
    clog<<"           -s       A constant to scale output (for example -pi for spectral function) "<<scale<<"\n";
    clog<<"*****************************************************\n"; 
    return 0;
  }
  function2D<double> fi(cn.size(),n1);
  mesh1D om(n1);
  for (int j=0; j<cn.size(); j++){
    if (cn[j]>m1){
      cerr<<"The column number "<<cn[j]<<" does not exist"<<endl;
      exit(1);
    }
  }
  ReadData(input, om, fi, m1, cn);
  om.SetUp(0.0);

//    Finds correct Lorentzian to substract from the imaginary part of Sigma
  vector<Lorentz> Slorentz(cn.size());
  for (int l=0; l<cn.size(); l++){
    int j;
    for (j=0; j<om.size(); j++){
      if (fabs(fi[l][j]*om.Dh(j))>maxS) break;
    }
    if (j<om.size()){
      int s = j-10>=0 ? j-10 : 0;
      int imax=s;
      for (j=s+1; j<om.size(); j++)
	if (fabs(fi[l][j])>fabs(fi[l][imax])) imax=j;
      
      double x0 = om[imax-1], x1 = om[imax], x2 = om[imax+1];
      double y0 = 1/fi[l][imax-1], y1 = 1/fi[l][imax], y2 = 1/fi[l][imax+1];
      double x02 = x0*x0, x12 = x1*x1, x22 = x2*x2;
      double a = ((x02-x12)*(y1 - y2) - (x12-x22)*(y0-y1))/((x2-x1)*(y0-y1)-(x1-x0)*(y1-y2))/2.;
      double P = (sqr(x0-a)-sqr(x1-a))/(y0-y1);
      double g2 = P*y0-sqr(x0-a);
      Slorentz[l].Set(a,sqrt(fabs(g2)),P);
      clog<<" Lorentz added for cn="<<cn[l]<<": "<<Slorentz[l].print()<<endl;
    }
    for (int i=0; i<om.size(); i++) fi[l][i] -= Slorentz[l](om[i]);
  }

  cout.precision(16);
  if (Ni<=0){
    for (int i=0; i<om.size(); i++){
      cout<<setw(25)<<om[i];
      for (int l=0; l<cn.size(); l++){
	double r = KramarsKronig(fi[l], om, om[i], i, fi[l][i]);
	cout<<setw(25)<<(r+Slorentz[l].realPart(om[i]))*scale<<setw(25)<<(fi[l][i]+Slorentz[l](om[i]))*scale;
      }
      cout<<endl;
    }
  } else{
    intpar p = om.Interp(0.0);
    cout<<0.0;
    dcomplex z=(0,1e-10);
    vector<double> S0(cn.size());
    for (int l=0; l<cn.size(); l++){
      S0[l] = fi[l](p);
      dcomplex r = KramarsKronig(fi[l], om, z, p.i, S0[l]);
      cout<<setw(25)<<(r + Slorentz[l].analytic(z))*scale;
    }
    cout<<endl;
    if (bosonic) z.imag()+=diom;
    else z.imag()+=0.5*diom;
    for (int i=0; i<Ni; i++){
      cout<<z.imag();
      for (int l=0; l<cn.size(); l++){
	dcomplex r = KramarsKronig(fi[l], om, z, p.i, S0[l]);
	cout<<setw(25)<<(r + Slorentz[l].analytic(z))*scale;
	z.imag()+=diom;
      }
      cout<<endl;
    }
  }

  if (rstdin) system("rm -f skrams.temp");

  return 0;
}
