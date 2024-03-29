#include <iostream>
#include <sstream>
#include <fstream>
#include <list>
#include <string>
#include <algorithm>
#include <set>
#include <deque>
#include <functional>
#include "interpolate.h"

using namespace std;
using namespace blitz;

bool CheckStream(istream& inputf, int& n, int& m)
{
  istream input(inputf.rdbuf());
  input.seekg(0,ios::beg);
  
  string str; int begincomm=0; n=0;
  getline(input,str); n++;
  if (!input.good()) {
    cerr << "ERROR: Wrong file format for hilbert or no data!" << endl;
    return false;
  }
  while (str.find('#')<string::npos and input.good()) {
    begincomm += 1;
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
  for (int i=0; i<begincomm; ++i) getline(inputf, str);
  if (!inputf){ cerr<<"Reopening didn't suceeded!"<<endl; return false;}
  return true;
}

bool ReadData(istream& input, Array<double,1>& om, Array<double,2>& fi, int mdata)
{
  clog<<"Reading data with "<<om.size()<<" entries"<<endl;
  Array<double,1> data(mdata);
  int i=-1;
  while (input && ++i<om.size()){
    for (int j=0; j<mdata; j++) input>>data(j);
    input.ignore(500,'\n');
    om(i) = data(0);
    fi(Range(0,mdata-2),i) = data(Range(1,mdata-1));
    //for (int j=1; j<mdata; j++) fi(j-1,i) = data(j);
  }
  return true;
}

class Fwidth{
  double w0, k;
public:
  Fwidth(double w0_, double k_) : w0(w0_), k(k_){}
  double operator()(double x)
  { return w0+k*abs(x); }
};

inline Array<double,1> GiveDh(const Array<double,1>& om){
  int N = om.size();
  Array<double,1> dh(N); 
  dh(0) = 0.5*(om(1)-om(0));
  dh(N-1)= 0.5*(om(N-1)-om(N-2));
  for (int i=1; i<N-1; ++i)
    dh(i) = 0.5*(om(i+1)-om(i-1));
  return dh;
}

Array<double,2> BroadenFunctions(const Array<double,1>& om, const Array<double,2>& fi, double width0, double kwidth, bool zeropad=true, int Ng_=200)
{
  int Ng = 2*(Ng_/2);
  Array<double,1> xg0, xg, fg;
  Fwidth fwidth(width0, kwidth);
  
  Array<double,1> sum(fi.extent(0));
  double om_a = om(0), om_b=om(om.size()-1);
  Array<double,2> fout(fi.extent(0), fi.extent(1));
  for (int i=0; i<om.size(); i++){
    double width = width0 + kwidth*fabs(om(i));
    if (not zeropad and (om(i)<width+om_a or om_b-width> om(i))){ // we need data outside known interval to broaden and we are not supposed to zeropad
      fout(Range::all(),i) = fi(Range::all(),i);
    }else{
      // Actual broadening
      double x0 = width/50;
      GiveTanMesh0(xg0, x0, width*20, Ng/2);
      int N = xg0.size();
      xg.resize(2*N);
      fg.resize(2*N);
      for (int i=0; i<N; ++i){
	xg(N-1-i) = -xg0(i);
	xg(N+i) = xg0(i);
      }
      for (int j=0; j<xg.size(); j++) // make normalized gaussian of width
      	fg(j) = exp(-xg(j)*xg(j)/(2*width*width))/sqrt(2*M_PI)/width;
      
      deque<double> epsx;
      for (int j=0; j<om.size(); j++) // any om_j is inside interval [om_i-xg[-1], om_i-xg[0] ] is added to epsx
	if (om(j) >= om(i)-xg(2*N-1) and om(j) <= om(i)-xg(0) ) epsx.push_back( om(j) );
      for (int j=0; j<xg.size(); j++)
	// xg[:]+om_i mesh is centered at om[i] and resolves gaussian.
	// if any point from mesh xg[:]+om_i is inside interval [om[0],om[-1]], we add it to mesh epsx
	if (xg(j)+om(i)>=om(0) and xg(j)+om(i)<=om(om.size()-1)) epsx.push_back( xg(j)+om(i) );
      sort(epsx.begin(),epsx.end());
      
      Array<double,1> eps(epsx.size());
      for (int j=0; j<epsx.size(); j++) eps(j) = epsx[j];
      Array<double,1> dh = GiveDh(eps);
      
      int posom = 0;
      int posxg = xg.size()-1;
      sum=0;
      double norm=0;
      for (int j=0; j<eps.size(); j++){
	intpar ip = InterpRight(om(i)-eps(j),posxg,xg);
	double wgh = linInterp(fg, ip)*dh(j);
	norm += wgh;
	intpar p = InterpLeft(eps(j),posom,om);
	for (int l=0; l<fi.extent(0); l++) sum(l) += linInterp(fi,l,p)*wgh;
      }
      for (int l=0; l<sum.size(); l++) fout(l,i) = sum(l)/norm;
    }
  }
  return fout;
}

double FindMaxDeriv(const Array<double,1>& om, const Array<double,2>& fi, double exclude=1.0)
{
  set<double,std::greater<double>> deriv;
  
  for (int i=0; i<om.size()-1; i++){
    if (fabs(om(i))>exclude){
      for (int l=0; l<fi.extent(0); l++){
	double df = (fi(l,i+1)-fi(l,i))/(om(i+1)-om(i));
	deriv.insert(fabs(df));
      }
    }
  }
  int k=0;
  double sum=0;
  for (double s : deriv){
    if (k>1 and k<=11) sum += s;
    ++k;
  }
  return (sum/10);
}

int main(int argc, char *argv[], char *env[])
{
  list<string> files;
  int Ng=200;
  double width0=0.01;
  double kwidth=0.0;
  double exclude=1.0;
  bool zeropad=true;
  bool width_set=false;
  int cn=2;
  int i=0;
  while (++i<argc){
    string str(argv[i]);
    if (str=="-w" and i<argc-1){  width0 = atof(argv[++i]); width_set=true;}
    if (str=="-k" and i<argc-1){  kwidth = atof(argv[++i]); width_set=true;}
    else {
      ifstream file(argv[i]);
      if (file) files.push_back(argv[i]);
    }
  }
  
  if (files.size()<1){
    clog<<"************** GAUSSIAN BROADENING *****************\n";
    clog<<"**                                                **\n";
    clog<<"**      Copyright Kristjan Haule, 12.10.2003      **\n";
    clog<<"****************************************************\n";
    clog<<"\n";
    clog<<"broad file [-w double] [-k double]\n" ;
    clog<<"Options:   file     Filename of the function\n";
    clog<<"           -w       Starting width of broadening w="<<width0<<"\n";
    clog<<"           -k       Linear increasing broadening with coefficient k="<<kwidth<<"\n";
    clog<<"           -z       pretend that function is zero outside z="<<zeropad<<"\n";
    clog<<"*****************************************************\n"; 
    return 0;
  }

  fstream input(files.begin()->c_str());
  int n1, m1;
  if (!CheckStream(input,n1,m1)) exit(2);
  Array<double,2> fi(m1-1,n1);
  Array<double,1> om(n1);
  ReadData(input, om, fi, m1);

  if (not width_set){
    kwidth=0.01*FindMaxDeriv(om,fi,exclude);
    clog<<"kwidth="<<kwidth<<endl;
  }

  Array<double,2> fout = BroadenFunctions(om,fi,width0,kwidth,zeropad,Ng);
  
  cout.precision(16);
  for (int i=0; i<om.size(); i++){
    cout<<setw(25)<<om(i)<<" ";
    for (int l=0; l<fi.extent(0); l++) cout<<setw(25)<<fout(l,i)<<" ";
    cout<<endl;
  }
  return 0;
}
