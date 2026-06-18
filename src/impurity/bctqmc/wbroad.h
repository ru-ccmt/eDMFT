#include <iostream>
#include <sstream>
#include <fstream>
#include <list>
#include <string>
#include <algorithm>
#include <set>
#include <deque>
#include <functional>
#include <numeric>
#include <gsl/gsl_sf_bessel.h>
#include "interpolate.h"

using namespace std;
using namespace blitz;

class Gauss_inc_w{
  // This is Gaussian of width k*sqrt(x0), i.e.,
  //   1/N(x0,k) * exp(-(x-x0)^2/(k*|x0|))
  // Here N(x0,k) is normalization that should preserve the integral, i.e.,
  //  Integrate[ 1/N(x0,k) * exp(-(x-x0)^2/(k*|x0|)) f(x), x] = fn(x0)
  //  and Integrate[fn(x0), x0] = Integrate[f(x), x]
  // In order to preserve the integral, we should have
  //   N(x0,k) ~ Integrate[exp(-(x-x0)^2/(k*|x|)), x] = 4*k*NN(|x0|/k), where NN(y)=y*cosh(y)*K1(y)
  //   with K1(y) being modified Bessel (Neuman) function.
public:
  double k, x0, norm, width;
  Gauss_inc_w(double _k_, double _x0_) : k(fabs(_k_)), x0(_x0_)
  {
    double x0_k = fabs(x0)/k;
    double dn;
    if (x0_k > 100.){
      dn = 0.5*sqrt(0.5*M_PI*x0_k)*(1 + 1.0/x0_k*(3./8.- 1.0/x0_k*(15./128.-1.0/x0_k*105./1024.)));
    }else{
      double K1 = gsl_sf_bessel_K1(x0_k); // K_1(x)
      dn = x0_k*cosh(x0_k)*K1;
    }
    norm = 4*k*dn;
    width = 2*k*fabs(x0);
    if (width<0.05*k) width=0.1*k;
  }
  double operator()(double x){
    return exp(-(x-x0)*(x-x0)/width)/norm;
  }
};

template<class T>
bl::Array<T,2> BroadenFunctions(const bl::Array<double,1>& om, const bl::Array<T,2>& fi, double kwidth, ostream& log, bool zeropad=true, int Ng_=200)
{
  int Ng = 2*(Ng_/2);
  bl::Array<double,1> xg0, xg, fg;
  if (fi.extent(1) != om.extent(0)) throw std::runtime_error("BroadenFunctions: fi.shape(1) must match om.size()");
  if (om.extent(0) < 2) throw std::runtime_error("om needs at least 2 points");
  if (fi.extent(0) <= 0 || fi.extent(1) <= 0) throw std::runtime_error("fi must be non-empty");
  if (!(kwidth >= 0.0)) throw std::runtime_error("kwidth must be non-negative");
  
  bl::Array<T,1> _sum_(fi.extent(0));
  double om_a = om(0), om_b=om(om.size()-1);
  bl::Array<T,2> fout(fi.extent(0), fi.extent(1));
  
  // we first set output==input, and later modify only if necessary and successful
  fout(bl::Range::all(),bl::Range::all()) = fi(bl::Range::all(),bl::Range::all());
  
  //log<<"Inside BroadenFunctions with kwidth="<<kwidth<<" om.size="<<om.size()<<" fi.shape="<<fi.shape()<<endl;
  //log<<"# BroadenFunction sum(fi)="<<sum(fi)<<endl;
  for (int i=0; i<om.size(); i++){
    double width = kwidth*sqrt(fabs(om(i)));
    //
    double dom = (i == 0) ? (om(1)-om(0)) : (i==om.size()-1) ? (om(om.size()-1)-om(om.size()-2)) : (0.5*(om(i+1)-om(i-1)));
    if (width<dom/3.) // the width to broaden is much smaller than distances between points, hence no point of computing convolution. just copy data.
      continue;
    if (not zeropad and (om(i)<width+om_a or om_b-width> om(i)))
      // we need data outside known interval to broaden and we are not supposed to zeropad. Should not happen
      continue;
    
    { // Actual broadening
      double x0 = width/50;
      GiveTanMesh0(xg0, x0, width*20, Ng/2, false, log);
      int N = xg0.size();
      xg.resize(2*N);
      fg.resize(2*N);
      for (int i=0; i<N; ++i){
	xg(N-1-i) = -xg0(i);
	xg(N+i) = xg0(i);
      }
      Gauss_inc_w giw(kwidth, om(i) );
      for (int j=0; j<xg.size(); j++){ // make normalized gaussian of width
	//double arg = xg(j)*xg(j)/(2*width*width);
      	//fg(j) = arg<200 ? exp(-arg)/sqrt(2*M_PI)/width : 0.0;
	fg(j) = giw(xg(j)+om(i));
      }
      
      deque<double> ep1, ep2;
      for (int j=0; j<om.size(); j++) // any om[j] that satisfies: om[i]-xg[-1]<= om[j] <= om[i]-xg[0]
	if (om(j) >= om(i)-xg(2*N-1) and om(j) <= om(i)-xg(0) ) ep1.push_back( om(j) );
      for (int j=0; j<xg.size(); j++) // any xg[j] that satisfies: om[0] <= xg[j]+om[i] <= om[-1]
	if (xg(j)+om(i)>=om(0) and xg(j)+om(i)<=om(om.size()-1)) ep2.push_back( xg(j)+om(i) );
      // For the algorithm below, eps[j] needs to satisfy:
      //    om[i]-xg[-1] <= eps[j] <= om[i]-xg[0]   and  om[0] <= eps[j] <= om[-1]
      // The first set of mesh points above eps[j]~om[j] above requires:
      // om[i]-xg[-1] <= om[j] <= om[i]-xg[0]   and  om[0] <= om[j] <= om[-1]
      // The second set of mesh points above eps[j]~ xg[j]+om[i] above requires:
      // om[i]-xg[-1] <= xg[j]+om[i] <= om[i]-xg[0]   and  om[0] <= xg[j]+om[i] <= om[-1]
      // which is equivalen to
      // -xg[-1]<= xg[j] <= -xg[0] and om[0]<= xg[j]+om[i] <= om[-1]

      //log<<"x0="<<x0<<" width="<<width<<" xg.size="<<xg.size()<<" fg.size="<<fg.size()<<" ep1.size="<<ep1.size()<<" ep2.size="<<ep2.size()<<endl;
      //log.flush();
      bl::Array<double,1> eps(ep1.size()+ep2.size());
      {// combine two sorted lists
	int i1=0;
	int i2=0;
	for (int i=0; i<eps.size(); ++i){
	  if (i1>=ep1.size()){ eps(i) = ep2[i2]; i2+=1; continue;}
	  if (i2>=ep2.size()){ eps(i) = ep1[i1]; i1+=1; continue;}
	  if (ep1[i1]<ep2[i2]){
	    eps(i) = ep1[i1];
	    i1+=1;
	  }else{
	    eps(i) = ep2[i2];
	    i2+=1;
	  }
	}
      }
      //log<<" sorted(epsx)="<<eps.size()<<endl;
      //log<<"epsx=[";
      //for (int l=0; l<eps.size(); ++l) log<<eps(l)<<",";
      //log<<"]"<<endl;
      //log.flush();
      //bl::Array<double,1> eps(epsx.size());
      //for (int j=0; j<epsx.size(); j++) eps(j) = epsx[j];
      bl::Array<double,1> dh(eps.size()); dh = GiveDh(eps);
      //log<<"dhx=[";
      //for (int l=0; l<eps.size(); ++l) log<<dh(l)<<",";
      //log<<"]"<<endl;
      //log.flush();
      int posom = 0;
      int posxg = xg.size()-1;
      _sum_(bl::Range::all())=0;    
      double norm=0;
      // fout[i] = Intg(fi(eps)*gauss(om[i]-eps)*deps)/Intg(gauss(om[i]-eps)*deps)
      // norm = Intg(gauss(om[i]-eps)*deps)
      // sum(l) = Intg(fi(l,eps)*gauss(om[i]-eps)*deps)
      //log<<"   . xg.size="<<xg.size()<<" om.size="<<om.size()<<" eps.size="<<eps.size()<<endl;
      for (int j=0; j<eps.size(); j++){
	intpar ip_om_eps = InterpRight(om(i)-eps(j),posxg,xg);    // om[i]-xg[-1] <= eps[j] <= om[i]-xg[0]
	intpar ip_eps = InterpLeft(eps(j),posom,om);              // om[0] <= eps[j] <= om[-1]
	//log<<j<<" ip_om_eps="<<ip_om_eps<<" ip_eps="<<ip_eps<<" om-eps="<<om(i)-eps(j)<<" eps(j)="<<eps(j)<<" xg=["<<xg(0)<<","<<xg(xg.size()-1)<<"], om=["<<om(0)<<","<<om(om.size()-1)<<"]"<<endl;
	// Conditions:
	//xg[0] <= om[i]-eps[j] <= xg[-1]
	//om[0] <= eps[j] <= om[-1]
	if (ip_om_eps.p>=0 and ip_om_eps.p<=1.0 and ip_eps.p>=0 and ip_eps.p<=1.0){ // no extrapolation allowed
	  double wgh = linInterp(fg, ip_om_eps)*dh(j); // gauss(om[i]-eps)*deps
	  norm += wgh;
	  for (int l=0; l<fi.extent(0); l++) _sum_(l) += linInterp(fi,l,ip_eps)*wgh; // _sum_[l] += \sum_eps fi[l](eps)*gauss(om[i]-eps)*deps
	}
	//log<<j<<"   wgh="<<norm<<" ff="<<sum<<endl;
      }
      if (norm!=0)
	for (int l=0; l<fi.extent(0); l++) fout(l,i) = _sum_(l)/norm;
      //log<<" om("<<i<<")="<<om(i)<<" fout="<<fout(bl::Range::all(),i)<<endl;
      //log<<"# om["<<i<<"]="<<om(i)<<" sum(fout)="<<sum(fout(bl::Range::all(),i))<<"sum(fg)="<<sum(fg)<<" sum(xg)="<<sum(xg)<<" norm="<<norm<<endl;
    }
  }
  //log<<"Finished BroadenFunctions with kwidth="<<kwidth<<endl;
  return fout;
}
double FindMaxDeriv(const bl::Array<double,1>& om, const bl::Array<double,2>& fi, double exclude=1.0)
{//The function estimates a typical large derivative magnitude of the dataset fi with respect to om(i) outside exclude window.
  // It is not looking for the absolute maximum (which might be noisy), but rather an average of the 10 largest derivative magnitudes,
  // skipping the two very largest one to avoid outliers.
  int how_many_average=20;
  std::vector<double> deriv;
  deriv.reserve((om.size()-1) * fi.extent(0));
  for (int i=0; i<om.size()-1; i++){
    if (fabs(om(i))>exclude){
      for (int l=0; l<fi.extent(0); l++){
        double df = (fi(l,i+1)-fi(l,i))/(om(i+1)-om(i));
	deriv.push_back(std::abs(df));
      }
    }
  }
  std::partial_sort(deriv.begin(), deriv.begin()+how_many_average+2, deriv.end(), std::greater<double>());
  double sum = std::accumulate(deriv.begin()+2, deriv.begin()+2+how_many_average, 0.0)/how_many_average;
  //cout<<"deriv=[";
  //for (int i=0; i<how_many_average+1; ++i) cout<<deriv[i]<<",";
  //cout<<"] : "<< sum << endl;
  return sum;
}

/*
template<class T>
Array<T,2> BroadenFunctions(const Array<double,1>& om, const Array<T,2>& fi, double width0, double kwidth, ostream& log, bool zeropad=true, int Ng_=200)
{
  int Ng = 2*(Ng_/2);
  Array<double,1> xg0, xg, fg;
  
  Array<T,1> _sum_(fi.extent(0));
  double om_a = om(0), om_b=om(om.size()-1);
  Array<T,2> fout(fi.extent(0), fi.extent(1));
  
  // we first set output==input, and later modify only if necessary and successful
  fout(Range::all(),Range::all()) = fi(Range::all(),Range::all());
  
  //log<<"Inside BroadenFunctions with kwidth="<<kwidth<<" om.size="<<om.size()<<" fi.shape="<<fi.shape()<<endl;
  log<<"# BroadenFunction sum(fi)="<<sum(fi)<<endl;
  for (int i=0; i<om.size(); i++){
    double width = width0 + kwidth*fabs(om(i));
    //
    double dom = (i == 0) ? (om(1)-om(0)) : (i==om.size()-1) ? (om(om.size()-1)-om(om.size()-2)) : (0.5*(om(i+1)-om(i-1)));
    if (width<dom/3.) // the width to broaden is much smaller than distances between points, hence no point of computing convolution. just copy data.
      continue;
    if (not zeropad and (om(i)<width+om_a or om_b-width> om(i)))
      // we need data outside known interval to broaden and we are not supposed to zeropad. Should not happen
      continue;
    
    { // Actual broadening
      double x0 = width/50;
      GiveTanMesh0(xg0, x0, width*20, Ng/2, false, log);
      int N = xg0.size();
      xg.resize(2*N);
      fg.resize(2*N);
      for (int i=0; i<N; ++i){
	xg(N-1-i) = -xg0(i);
	xg(N+i) = xg0(i);
      }
      for (int j=0; j<xg.size(); j++){ // make normalized gaussian of width
	double arg = xg(j)*xg(j)/(2*width*width);
      	fg(j) = arg<200 ? exp(-arg)/sqrt(2*M_PI)/width : 0.0;
      }
      
      deque<double> ep1, ep2;
      for (int j=0; j<om.size(); j++) // any om[j] that satisfies: om[i]-xg[-1]<= om[j] <= om[i]-xg[0]
	if (om(j) >= om(i)-xg(2*N-1) and om(j) <= om(i)-xg(0) ) ep1.push_back( om(j) );
      for (int j=0; j<xg.size(); j++) // any xg[j] that satisfies: om[0] <= xg[j]+om[i] <= om[-1]
	if (xg(j)+om(i)>=om(0) and xg(j)+om(i)<=om(om.size()-1)) ep2.push_back( xg(j)+om(i) );
      // For the algorithm below, eps[j] needs to satisfy:
      //    om[i]-xg[-1] <= eps[j] <= om[i]-xg[0]   and  om[0] <= eps[j] <= om[-1]
      // The first set of mesh points above eps[j]~om[j] above requires:
      // om[i]-xg[-1] <= om[j] <= om[i]-xg[0]   and  om[0] <= om[j] <= om[-1]
      // The second set of mesh points above eps[j]~ xg[j]+om[i] above requires:
      // om[i]-xg[-1] <= xg[j]+om[i] <= om[i]-xg[0]   and  om[0] <= xg[j]+om[i] <= om[-1]
      // which is equivalen to
      // -xg[-1]<= xg[j] <= -xg[0] and om[0]<= xg[j]+om[i] <= om[-1]

      //log<<"x0="<<x0<<" width="<<width<<" xg.size="<<xg.size()<<" fg.size="<<fg.size()<<" ep1.size="<<ep1.size()<<" ep2.size="<<ep2.size()<<endl;
      //log.flush();
      Array<double,1> eps(ep1.size()+ep2.size());
      {// combine two sorted lists
	int i1=0;
	int i2=0;
	for (int i=0; i<eps.size(); ++i){
	  if (i1>=ep1.size()){ eps(i) = ep2[i2]; i2+=1; continue;}
	  if (i2>=ep2.size()){ eps(i) = ep1[i1]; i1+=1; continue;}
	  if (ep1[i1]<ep2[i2]){
	    eps(i) = ep1[i1];
	    i1+=1;
	  }else{
	    eps(i) = ep2[i2];
	    i2+=1;
	  }
	}
      }

      //log<<" sorted(epsx)="<<eps.size()<<endl;
      //log<<"epsx=[";
      //for (int l=0; l<eps.size(); ++l) log<<eps(l)<<",";
      //log<<"]"<<endl;
      //log.flush();
      
      //Array<double,1> eps(epsx.size());
      //for (int j=0; j<epsx.size(); j++) eps(j) = epsx[j];
      Array<double,1> dh(eps.size()); dh = GiveDh(eps);
      //log<<"dhx=[";
      //for (int l=0; l<eps.size(); ++l) log<<dh(l)<<",";
      //log<<"]"<<endl;
      //log.flush();
 

      int posom = 0;
      int posxg = xg.size()-1;
      _sum_(Range::all())=0;    
      double norm=0;
      // fout[i] = Intg(fi(eps)*gauss(om[i]-eps)*deps)/Intg(gauss(om[i]-eps)*deps)
      // norm = Intg(gauss(om[i]-eps)*deps)
      // sum(l) = Intg(fi(l,eps)*gauss(om[i]-eps)*deps)
      //log<<"   . xg.size="<<xg.size()<<" om.size="<<om.size()<<" eps.size="<<eps.size()<<endl;
      for (int j=0; j<eps.size(); j++){
	intpar ip_om_eps = InterpRight(om(i)-eps(j),posxg,xg);    // om[i]-xg[-1] <= eps[j] <= om[i]-xg[0]
	intpar ip_eps = InterpLeft(eps(j),posom,om);              // om[0] <= eps[j] <= om[-1]
	//log<<j<<" ip_om_eps="<<ip_om_eps<<" ip_eps="<<ip_eps<<" om-eps="<<om(i)-eps(j)<<" eps(j)="<<eps(j)<<" xg=["<<xg(0)<<","<<xg(xg.size()-1)<<"], om=["<<om(0)<<","<<om(om.size()-1)<<"]"<<endl;
	// Conditions:
	//xg[0] <= om[i]-eps[j] <= xg[-1]
	//om[0] <= eps[j] <= om[-1]
	if (ip_om_eps.p>=0 and ip_om_eps.p<=1.0 and ip_eps.p>=0 and ip_eps.p<=1.0){ // no extrapolation allowed
	  double wgh = linInterp(fg, ip_om_eps)*dh(j); // gauss(om[i]-eps)*deps
	  norm += wgh;
	  for (int l=0; l<fi.extent(0); l++) _sum_(l) += linInterp(fi,l,ip_eps)*wgh; // _sum_[l] += \sum_eps fi[l](eps)*gauss(om[i]-eps)*deps
	}
	//log<<j<<"   wgh="<<norm<<" ff="<<sum<<endl;
      }
      if (norm!=0)
	for (int l=0; l<fi.extent(0); l++) fout(l,i) = _sum_(l)/norm;
      //log<<" om("<<i<<")="<<om(i)<<" fout="<<fout(Range::all(),i)<<endl;
      //log<<"# om["<<i<<"]="<<om(i)<<" sum(fout)="<<sum(fout(Range::all(),i))<<"sum(fg)="<<sum(fg)<<" sum(xg)="<<sum(xg)<<" norm="<<norm<<endl;
    }
  }
  //log<<"Finished BroadenFunctions with kwidth="<<kwidth<<endl;
  return fout;
}
*/

/*
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
*/

bl::Array<double,2> FindRangeToRemove(const bl::Array<double,1>& om, const bl::Array<double,2>& fi, double exclude=3.0, double max_der=0.5)
{
  //analyzes a 2D dataset fi(l,i) defined on a 1D mesh om(i) and identifies, for each row l,
  // the range of om values where the function changes too rapidly (i.e., has large derivatives).
  // It returns those ranges as a 2-column matrix x_remove(l, 0:1) — roughly the left and right boundaries of “bad” regions to be excluded.
  //cout<<"om2="<<om<<endl;
  bl::Array<double,2> x_remove(fi.extent(0),2);
  for (int j=0; j<fi.extent(0); ++j){
    x_remove(j,0) = -1e100;
    x_remove(j,1) =  1e100;
    for (int i=1; i<om.size()-1; i++){
      if (om(i)>exclude){
        double df = 0.5*abs((fi(j,i+1)-fi(j,i))/(om(i+1)-om(i))) + 0.5*abs((fi(j,i)-fi(j,i-1))/(om(i)-om(i-1)));
	//cout<<"om["<<i<<"]="<<om(i)<<" f="<<fi(j,i)<<" df="<<df<<" df>max_der="<< (df>max_der) <<endl;
        if (df>max_der){
          x_remove(j,1) = om(i);
          break;
        }
      }
    }
    for (int i=om.size()-1; i>0; i--){
      if (om(i)<-exclude){
        double df = 0.5*abs((fi(j,i+1)-fi(j,i))/(om(i+1)-om(i))) + 0.5*abs((fi(j,i)-fi(j,i-1))/(om(i)-om(i-1)));
	//cout<<"om["<<i<<"]="<<om(i)<<" f="<<fi(j,i)<<" df="<<df<<" df>max_der="<< (df>max_der) <<endl;
        if (df>max_der){
          x_remove(j,0) = om(i);
          break;
        }
      }
    }
  }
  return x_remove;
}

/*
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

inline Array<double,1> GiveDh(const Array<double,1>& om){
  int N = om.size();
  Array<double,1> dh(N); 
  dh(0) = 0.5*(om(1)-om(0));
  dh(N-1)= 0.5*(om(N-1)-om(N-2));
  for (int i=1; i<N-1; ++i)
    dh(i) = 0.5*(om(i+1)-om(i-1));
  return dh;
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

  Array<double,2> fout = BroadenFunctions(om,fi,width0,kwidth,std::cout,zeropad,Ng);

  cout.precision(16);
  for (int i=0; i<om.size(); i++){
    cout<<setw(25)<<om(i)<<" ";
    for (int l=0; l<fi.extent(0); l++) cout<<setw(25)<<fout(l,i)<<" ";
    cout<<endl;
  }
  return 0;
}
*/
