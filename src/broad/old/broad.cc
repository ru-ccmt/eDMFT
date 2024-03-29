// @Copyright 2007 Kristjan Haule
#include <iostream>
#include <sstream>
#include <fstream>
#include <list>
#include <algorithm>
#include <deque>
#include "smesh.h"
#include "sfunction.h"

using namespace std;

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

bool ReadData(istream& input, mesh1D& om, function2D<double>& fi, int mdata)
{
  clog<<"Reading data with "<<om.size()<<" entries"<<endl;
  vector<double> data(mdata);
  int i=-1;
  while (input && ++i<om.size()){
    for (int j=0; j<mdata; j++) input>>data[j];
    input.ignore(500,'\n');
    om[i] = data[0];
    for (int j=1; j<mdata; j++) fi[j-1][i] = data[j];
  }
  return true;
}

class Fwidth{
  double w0, k;
public:
  Fwidth(double w0_, double k_) : w0(w0_), k(k_){}
  double operator()(double x)
  {
    return w0+k*abs(x);
  }
};

int main(int argc, char *argv[], char *env[])
{
  list<string> files;
  int Ng=200;
  double width0=0.4;
  double kwidth=0.;
  int cn=2;
  int i=0;
  while (++i<argc){
    string str(argv[i]);
    if (str=="-w" && i<argc-1)  width0 = atof(argv[++i]);
    if (str=="-k" && i<argc-1)  kwidth = atof(argv[++i]);
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
    clog<<"broad file [-cn int] [-w double]\n" ;
    clog<<"Options:   file     Filename of the function\n";
    clog<<"           -w       Starting width of broadening "<<width0<<"\n";
    clog<<"           -k       Coefficients of broadeining "<<kwidth<<"\n";
    clog<<"*****************************************************\n"; 
    return 0;
  }

  fstream input(files.begin()->c_str());
  int n1, m1;
  if (!CheckStream(input,n1,m1)) exit(2);
  function2D<double> fi(m1-1,n1);
  mesh1D om(n1);
  ReadData(input, om, fi, m1);
  om.SetUp(0.0);

  Ng = 2*(Ng/2);
  mesh1D xg(Ng);
  function1D<double> fg(xg.size());
  
  Fwidth fwidth(width0, kwidth);

  function1D<double> sum(fi.size_N());
  cout.precision(16);
  for (int i=0; i<om.size(); i++){
    double width = fwidth(om[i]);
    //if (om[i]-om[0]>width && om[om.size()-1]-om[i]>width){
    {
      xg.MakeTanMesh(Ng/2, 0.0, width, width/50., width*20);
      xg.resize(Ng);
      for (int j=0; j<Ng/2; j++) xg[Ng/2+j] = xg[j];
      for (int j=0; j<Ng/2; j++) xg[Ng/2-j-1] = -xg[Ng/2+j];
      sort(xg.begin(),xg.end());
      xg.SetUp(0);
      for (int j=0; j<xg.size(); j++)
	fg[j] = exp(-xg[j]*xg[j]/(2*width*width))/sqrt(2*M_PI)/width;
      
      deque<double> epsx;
      for (int j=0; j<om.size(); j++){
	// any om_j is inside interval [om_i-xg[-1], om_i-xg[0] ] is added to epsx
	if (om[j] >= om[i]-xg[xg.size()-1] and om[j] <= om[i]-xg[0] ) epsx.push_back( om[j] );
      }
      for (int j=0; j<xg.size(); j++){
	// xg[:]+om_i mesh is centered at om[i] and resolves gaussian.
	// if any point from mesh xg[:]+om_i is inside interval [om[0],om[-1]], we add it to mesh epsx
	if (xg[j]+om[i] >= om[0] && xg[j]+om[i] <= om[om.size()-1]) epsx.push_back( xg[j]+om[i] );
      }
      sort(epsx.begin(),epsx.end());
      
      mesh1D eps(epsx.size());
      for (int j=0; j<epsx.size(); j++) eps[j] = epsx[j];
      eps.SetUp(0);
      
      tint posom = om.InitInterpLeft();
      tint posxg = xg.InitInterpRight();
      sum=0;
      double norm=0;
      for (int j=0; j<eps.size(); j++){
	double wgh = fg(xg.InterpRight(om[i]-eps[j],posxg))*eps.Dh(j);
	norm += wgh;
	intpar p = om.InterpLeft(eps[j],posom);
	for (int l=0; l<fi.size_N(); l++)	sum[l] += fi[l](p)*wgh;
      }
      cout<<setw(25)<<om[i]<<" ";
      for (int l=0; l<sum.size(); l++) cout<<setw(25)<<sum[l]/norm<<" ";
      cout<<endl;
    }
    //else{
    //  cout<<setw(25)<<om[i]<<" ";
    //  for (int l=0; l<sum.size(); l++) cout<<setw(25)<<fi[l][i]<<" ";
    //  cout<<endl;
    //}
  }
  
  return 0;
}
