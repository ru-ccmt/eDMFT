#include <blitz/array.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <limits>
#include <cassert>
#include <complex>
using namespace std;
using namespace blitz;

inline Array<double,1> GiveDh(const Array<double,1>& om){
  int N = om.size();
  Array<double,1> dh(N); 
  dh(0) = 0.5*(om(1)-om(0));
  dh(N-1)= 0.5*(om(N-1)-om(N-2));
  for (int i=1; i<N-1; ++i)
    dh(i) = 0.5*(om(i+1)-om(i-1));
  return dh;
}

#include "wbroad.h"


static inline vector<double> split_doubles(const string& line) {
    std::istringstream iss(line);
    vector<double> v;
    double x;
    while (iss >> x) v.push_back(x);
    return v;
}

void Load_g00_file_stackedG(const string& path,
                            Array<double,1>& xt,
                            Array<double,2>& G00_stacked) // (2*Nbath x Nw)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open file: " + path);

    vector<vector<double>> rows;
    string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto vals = split_doubles(line);
        if (!vals.empty()) rows.push_back(std::move(vals));
    }
    if (rows.empty()) throw std::runtime_error("Empty file: " + path);

    const int Nw = static_cast<int>(rows.size());
    const int ncols = static_cast<int>(rows[0].size());
    if (ncols < 1 || ((ncols - 1) % 4) != 0)
        throw std::runtime_error("Unexpected column count in " + path);
    const int Nbath = (ncols - 1) / 4;

    xt.resize(Nw);
    G00_stacked.resize(4*Nbath, Nw);
    
    for (int i = 0; i < Nw; ++i) {
        const auto& r = rows[i];
        xt(i) = r[0];
        int c = 1;
        // Re/Im G00 first
        for (int ip = 0; ip < 4*Nbath; ++ip) G00_stacked(ip,        i) = r[c++]; // Re G00
        if (c != ncols) throw std::runtime_error("Parse error at line " + std::to_string(i));
    }
}

Array<double,2> Broaden_Quantities(int n, const Array<double,1>& xt, const Array<double,2>& G00,
				   double Wtypical, double broad, double exclude, double max_der, double artificial_beta)
{
  int Nbath_unique = G00.extent(0)/4;
  Array<double,2> zg00 = G00.copy();
  double md = FindMaxDeriv(xt, zg00, Wtypical);
  double ky = 0.015*md*broad*n*n/16.;
  double kwidth = ky/(1.+ky);
  double _exclude_ = exclude;
  if (n==1) _exclude_ = 1000; // the lowest orders should not be cut off at all.

  cout<<"max_der="<<max_der<<" exclude="<<_exclude_<<endl;
  
  Array<double,2> x_remove = FindRangeToRemove(xt, zg00, _exclude_, max_der);
  Array<double,2> y_remove(Nbath_unique,2);
  for (int ib=0; ib<Nbath_unique; ++ib){
    y_remove(ib,0) = -1e100;
    y_remove(ib,1) = 1e100;
    for (int itype=0; itype<4; ++itype){// itype=Re(G00),Im(G00),Re(F00),Im(F00)
      int row = ib + Nbath_unique*itype;
      y_remove(ib,0) = std::max(y_remove(ib,0), x_remove(row,0));
      y_remove(ib,1) = std::min(y_remove(ib,1), x_remove(row,1));
    }
    //y_remove(ib,0) = max(x_remove(Range(4*ib,4*ib+3),0));
    //y_remove(ib,1) = min(x_remove(Range(4*ib,4*ib+3),1));
  }
  cout<<" broadening kwidth["<<n<<"]="<<kwidth<<" remove=";
  for (int j=0; j<Nbath_unique; ++j){
    cout<<"(";
    if (y_remove(j,0)>xt(0)) cout<<y_remove(j,0);
    cout<<",";
    if (y_remove(j,1)<xt(xt.extent(0)-1)) cout<<y_remove(j,1);
    cout<<") ";
  }
  cout<<endl;
  //cout<<" broadening kwidth["<<iorder<<"]="<<kwidth<<endl;
  zg00 = BroadenFunctions(xt,zg00,kwidth,cout);
  //
  for (int j=0; j<Nbath_unique; ++j){
    if (y_remove(j,1)<xt(xt.extent(0)-1)){
      double x_start = y_remove(j,1)-(1.0+0.15*y_remove(j,1))/artificial_beta; // because of broadening we should start cutting off a bit earlier than at y_remove, and when y_remove is larger spectra is broader, hence should cutoff earlier
      for (int i=0; i<xt.extent(0); ++i)
	zg00(Range(4*j,4*j+3),i) *= fermi((xt(i)-x_start)*artificial_beta);
    }
    if (y_remove(j,0)>xt(0)){
      double x_start = y_remove(j,0)+(1.0+0.15*fabs(y_remove(j,0)))/artificial_beta; // because of broadening we should start cutting off a bit earlier than at y_remove, and when y_remove is larger spectra is broader, hence should cutoff earlier
      for (int i=0; i<xt.extent(0); ++i)
	zg00(Range(4*j,4*j+3),i) *= 1-fermi((xt(i)-x_start)*artificial_beta);
    }
  }
  return zg00;
}
int main()
{
  Array<double,1> xt;
  Array<double,2> G00_st;
  int n=3;
  Load_g00_file_stackedG("g00_."+to_string(n), xt, G00_st);
  double broad=2.0, Wtypical=1.0, exclude=1.0, max_der=10., artificial_beta=2.;
  auto G00_br = Broaden_Quantities(n, xt, G00_st,
				   Wtypical, broad, exclude, max_der, artificial_beta);

  ofstream gout("_g00."s+to_string(n)); gout.precision(12);
  for (int i=0; i<xt.size(); ++i){
    gout<<xt(i)<<" ";
    for (int ip=0; ip<G00_br.extent(0); ++ip) gout<<G00_br(ip,i)<<" ";
    gout<<endl;
  }
}
