#ifndef _BASICS
#define _BASICS
#include <chrono>
#include <deque>
#include <vector>
#include <set>
#include <map>
#include <regex>
#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>
#define FMT_HEADER_ONLY
#include "fmt/core.h"
#include "fmt/format.h"
#include "fmt/ranges.h"
#include <blitz/array.h>
#include "interpolate.h"
using namespace blitz;
using namespace std;
using namespace std::chrono;

typedef std::chrono::duration<float> fsec;
typedef std::chrono::time_point<std::chrono::high_resolution_clock> fclock;
//typedef std::chrono::high_resolution_clock fclock;
fsec zero_t = std::chrono::duration<float,std::exa>::zero();

double get_cpu_time() {
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    return r.ru_utime.tv_sec + r.ru_utime.tv_usec / 1e6
         + r.ru_stime.tv_sec + r.ru_stime.tv_usec / 1e6;
}

inline int isign(int a){
  if (a>0) return 1;
  if (a<0) return -1;
  return 0;
}
int count_zeros(const Array<int,1>& dpseudo )
{
  int ni=0;
  for (int i=0; i<dpseudo.extent(0); ++i)
    if (dpseudo(i)==0) ni+=1;
  return ni;
}
deque<int> which_have_same_equals(const Array<int,1>& how_many_equal){
  int nmax=0;
  for (int i=0; i<how_many_equal.size(); ++i)
    if (how_many_equal(i)>nmax) nmax = how_many_equal(i);
  deque<int> candidates;
  for (int i=0; i<how_many_equal.size(); ++i)
    if (how_many_equal(i) == nmax) candidates.push_back(i);
  return candidates;
}
deque<int> which_have_same_equals_but_not_same(const Array<int,1>& how_many_equal, int two_n){
  int nmax=0;
  for (int i=0; i<how_many_equal.size(); ++i)
    if (how_many_equal(i)!=two_n and how_many_equal(i)>nmax) nmax = how_many_equal(i);
  deque<int> candidates;
  for (int i=0; i<how_many_equal.size(); ++i)
    if (how_many_equal(i) == nmax) candidates.push_back(i);
  return candidates;
}

template<int N>
bool operator!=(TinyVector<int,N> v1, TinyVector<int,N> v2){
  for (int i=0; i<N; ++i)
    if (v1[i]!=v2[i]) return true;
  return false;
}
template<int N>
ostream& operator<<(ostream& out, const TinyVector<int,N>& v)
{
  out << "[";
  for (int i=0; i<N-1; ++i) out<< v[i] << ",";
  out << v[N-1]<<"]";
  return out;
}
template<int N>
ostream& operator<<(ostream& out, const TinyVector<double,N>& v)
{
  out << "[";
  for (int i=0; i<N-1; ++i) out<< v[i] << ",";
  out << v[N-1]<<"]";
  return out;
}
ostream& operator<<(ostream& out, const set<int>& v){
  out << "(";
  out << fmt::format("{}", fmt::join(v, ","));
  out << ")";
  return out;
}
ostream& operator<<(ostream& out, const Array<int,1>& v){
  int N = v.extent(0);
  out << "[";
  for (int i=0; i<N-1; ++i) out<< v(i) << ",";
  out << v(N-1) << "]";
  return out;
}
ostream& operator<<(ostream& out, const Array<double,1>& v){
  int N = v.extent(0);
  out << "[";
  for (int i=0; i<N-1; ++i) out<< v(i) << ",";
  out << v(N-1) << "]";
  return out;
}
template<class T>
string to_string(const Array<T,1>& v)
{
  int N = v.extent(0);
  string out="[";
  for (int i=0; i<N-1; ++i) out += to_string(v(i))+",";
  out += to_string(v(N-1))+"]";
  return out;
}
ostream& operator<<(ostream& out, const Array<unsigned char,1>& v){
  int N = v.extent(0);
  out << "[";
  for (int i=0; i<N-1; ++i) out<< static_cast<int>(v(i)) << ",";
  out << static_cast<int>(v(N-1)) << "]";
  return out;
}
ostream& operator<<(ostream& out, const deque<int>& v){
  int N = v.size();
  out << "[";
  for (int i=0; i<N-1; ++i) out<< v[i] << ",";
  out << v[N-1] << "]";
  return out;
}
ostream& operator<<(ostream& out, const vector<int>& v){
  int N = v.size();
  out << "[";
  for (int i=0; i<N-1; ++i) out<< v[i] << ",";
  out << v[N-1] << "]";
  return out;
}
/*
ostream& operator<<(ostream& out, const Array<double,2>& v){
  int N = v.extent(0);
  int Nd=v.extent(1);
  out << "[";
  for (int i=0; i<N; ++i){
    out <<"(";
    for (int j=0; j<Nd-1; ++j) out<< v(i,j) << ",";
    if (i!=N-1) out <<v(i,Nd-1)<<"),";
    else out<<v(i,Nd-1)<<")";
  }
  out << "]";
  return out;
}
*/
ostream& operator<<(ostream& out, const Array<TinyVector<int,2>,1>& v){
  out << "[";
  for (int i=0; i<v.extent(0); ++i){
    out <<"("<<v(i)[0]<<","<<v(i)[1]<<")";
    if (i<v.extent(0)-1) out<<",";
  }
  out << "]";
  return out;
}

static std::string PReal(double x) {
  if (x == 0.0) return "0";
  static const std::regex re(R"(([-+]?\d*\.?\d+)[eE]([-+]?\d+))");
  auto s = fmt::format("{:.8g}", x);
  return std::regex_replace(s, re, "$1*10^$2");
}
static std::string PComplex(const std::complex<double>& v) {
  auto re = PReal(v.real());
  auto im = PReal(std::abs(v.imag()));
  if (v.imag() >= 0) 
    return re + "+" + im + "*I";
  else
    return re + "-" + im + "*I";
}
template<typename T>
std::string to_str(const T& x, int precision = 8)
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision);

    if constexpr (std::is_same_v<T, std::complex<double>>) {
        oss << x.real() << " " << x.imag();
    } else if constexpr (std::is_same_v<T, double>) {
        oss << x;
    } else {
        static_assert(std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>,
                      "Unsupported type for to_str()");
    }

    return oss.str();
}
ostream& operator<<(ostream& out, const Array<double,2>& v){
  int N = v.extent(0);
  int Nd=v.extent(1);
  out << "{";
  for (int i=0; i<N; ++i){
    out <<"{";
    for (int j=0; j<Nd-1; ++j) out<< v(i,j)<< ", ";
    if (i!=N-1) out << v(i,Nd-1)<<"},";
    else out<< v(i,Nd-1)<<"}";
  }
  out << "}";
  return out;
}
ostream& operator<<(ostream& out, const Array<complex<double>,2>& v){
  int N = v.extent(0);
  int Nd=v.extent(1);
  out << "{";
  for (int i=0; i<N; ++i){
    out <<"{";
    for (int j=0; j<Nd-1; ++j) out<< PComplex(v(i,j))<< ", ";
    if (i!=N-1) out << PComplex(v(i,Nd-1))<<"},";
    else out<< PComplex(v(i,Nd-1))<<"}";
  }
  out << "}";
  return out;
}
ostream& operator<<(ostream& out, const Array<TinyVector<int,2>,2>& v){
  int N = v.extent(0);
  int Nd=v.extent(1);
  out << "[";
  for (int i=0; i<N; ++i){
    out <<"{";
    for (int j=0; j<Nd-1; ++j)
                out << "("<<v(i,j)[0]<<","<<v(i,j)[1]<<"),";
    if (i!=N-1) out << "("<<v(i,Nd-1)[0]<<","<<v(i,Nd-1)[1]<<")},";
    else        out << "("<<v(i,Nd-1)[0]<<","<<v(i,Nd-1)[1]<<")}";
  }
  out << "]";
  return out;
}

ostream& operator<<(ostream& out, const Array<complex<double>,1>& v){
  int N = v.extent(0);
  out << "[";
  for (int i=0; i<N-1; ++i) out<< v(i) << ",";
  out << v(N-1) << "]";
  return out;
}

template<typename T, int N>
bool allclose(const Array<T,N>& a, const Array<T,N>& b, double rtol=1e-5, double atol=1e-8)
{
  if (a.shape() != b.shape()) {
    return false;
  }
  for (typename Array<T,N>::const_iterator it_a = a.begin(), it_b = b.begin(); it_a != a.end(); ++it_a, ++it_b){
    if (!(std::abs(*it_a - *it_b) <= atol + rtol * std::abs(*it_b))){
      return false;
    }
  }
  return true;
}
template<int N>
bool allclose(const blitz::Array<int,N>& a, const blitz::Array<int,N>& b)
{
  if (a.shape() != b.shape()){
    return false;
  }
  for (typename blitz::Array<int,N>::const_iterator it_a = a.begin(), it_b = b.begin(); it_a != a.end(); ++it_a, ++it_b){
    if (*it_a != *it_b) {
      return false;
    }
  }
  return true;
}
/*
ostream& operator<<(ostream& out, const set<int>& v){
  int N = v.size();
  out << "[";
  for (auto const &i: v) cout<< i << ",";
  out<<"]";
  return out;
}
*/

ostream& operator<<(ostream& out, const intpar& ip){
  out << "("<<ip.i<<","<<ip.p<<")";
  return out;
}
ostream& operator<<(ostream& out, const Array<int,2>& diag)
{
  int nn = diag.extent(0);
  int nm = diag.extent(1);
  out<<"[";
  for (int i=0; i<nn-1; ++i){
    out<<"(";
    for (int j=0; j<nm-1; ++j) out<<diag(i,j)<<",";
    out<<diag(i,nm-1)<<"),";
  }
  out<<"(";
  for (int j=0; j<nm-1; ++j) out<<diag(nn-1,j)<<",";
  out<<diag(nn-1,nm-1) << ")";
  out << "]";
  return out;
}
ostream& operator<<(ostream& out, const Array<bool,2>& diag)
{
  int nn = diag.extent(0);
  int nm = diag.extent(1);
  out<<"[";
  for (int i=0; i<nn-1; ++i){
    out<<"(";
    for (int j=0; j<nm-1; ++j) out<<diag(i,j)<<",";
    out<<diag(i,nm-1)<<"),";
  }
  out<<"(";
  for (int j=0; j<nm-1; ++j) out<<diag(nn-1,j)<<",";
  out<<diag(nn-1,nm-1) << ")";
  out << "]";
  return out;
}
ostream& operator<<(ostream& out, const Array<unsigned char,2>& diag)
{
  int nn = diag.extent(0);
  int nm = diag.extent(1);
  for (int i=0; i<nn; ++i){
    for (int j=0; j<nm; ++j) out<< static_cast<int>(diag(i,j));
    out<<endl;
  }
  return out;
}

ostream& operator<<(ostream& out, const vector<deque<int>>& d)
{
  int nn = d.size();
  out<<"[";
  for (int i=0; i<nn-1; ++i){
    int nm=d[i].size();
    out<<"(";
    for (int j=0; j<nm-1; ++j) out<<d[i][j]<<",";
    out<<d[i][nm-1]<<"),";
  }
  int nm=d[nn-1].size();
  out<<"(";
  for (int j=0; j<nm-1; ++j) out<<d[nn-1][j]<<",";
  out<<d[nn-1][nm-1] << ")";
  out << "]";
  return out;
}
ostream& operator<<(ostream& out, const map<int,int>& v)
{
  out<<"[";
  for (auto itt=v.begin(); itt!=v.end(); ++itt)
    out<<"("<<itt->first<<":"<<itt->second<<"), ";
  out << "]";
  return out;
}
ostream& operator<<(ostream& out, const vector<string>& v)
{
  out<<"[";
  for (int it=0; it<v.size()-1; ++it) out<<v[it]<<",";
  out << v[v.size()-1]<<"]";
  return out;
}
ostream& operator<<(ostream& out, const Array<deque<int>,2>& ftype_real)
{
  out<<"[";
  for (int i=0; i<ftype_real.extent(0); ++i){
    for (int j=0; j<ftype_real.extent(1); ++j)
      if (ftype_real(i,j).size()>0)
	out<<ftype_real(i,j)<<";";
    out<<"|";
  }
  out<<"]";
  return out;
}

ostream& operator<<(ostream& out, const Array<long,2>& A)
{
  int nn = A.extent(0);
  int nm = A.extent(1);
  out<<"[";
  for (int i=0; i<nn-1; ++i){
    out<<"(";
    for (int j=0; j<nm-1; ++j) out<<A(i,j)<<",";
    out<<A(i,nm-1)<<"),";
  }
  out<<"(";
  for (int j=0; j<nm-1; ++j) out<<A(nn-1,j)<<",";
  out<<A(nn-1,nm-1) << ")";
  out << "]";
  return out;
}


bool operator==(const Array<int,2>& diag1, const Array<int,2>& diag2){
  if (diag1.extent(0) != diag2.extent(0)) return false;
  if (diag1.extent(1) != diag2.extent(1)) return false;
  for (int i=0; i<diag1.extent(0); ++i)
    for (int j=0; j<diag1.extent(1); ++j)
      if ( diag1(i,j) != diag2(i,j) ) return false;
  return true;
}


bool operator==(const Array<unsigned char,1>& d1, const Array<unsigned char,1>& d2)
{
  if (d1.extent(0)!=d2.extent(0)) return false;
  for (int i=0; i<d1.extent(0); ++i)
    if ( d1(i) != d2(i) ) return false;
  return true;
}

bool operator==(const Array<int,1>& pseudo1, const Array<int,1>& pseudo2){
  if (pseudo1.extent(0) != pseudo2.extent(0)) return false;
  for (int i=0; i<pseudo1.extent(0); ++i)
    if (pseudo1(i) != pseudo2(i)) return false;
  return true;
}

void swap(int& a, int& b){
  //cout << "Swaping " << a << " and " << b << endl;
  int ca = a;
  a = b;
  b = ca;
}

void Print(ostream& out, const Array<int,2>& wtb){
  for (int i=0; i<wtb.extent(0); ++i){
    for (int j=0; j<wtb.extent(1); ++j)
      out<<setw(2)<<wtb(i,j)<<" ";
    out<<endl;
  }
}


complex<double> Product(const Array<complex<double>,1>& A)
{
  complex<double> res=1.0;
  for (int i=0; i<A.size(); ++i) res *= A(i);
  return res;
}

///////////////////////////////////////////////////
// This should allow structured binding
// for TinyVector<T,N>, for example
// auto [a1,a2,a3] = TinyVector<int,3>(b1,b2,b3)
////////////////////////////////////////////////////
/*
namespace blitz {
  template<std::size_t I, typename T, int N>
  constexpr T& get(TinyVector<T,N>& v)  noexcept
  {
    static_assert(I < N);
    return v[I];
  }
  
  template<std::size_t I, typename T, int N>
  constexpr const T& get(const TinyVector<T,N>& v)  noexcept
  {
    static_assert(I < N);
    return v[I];
  }
}
namespace std {
  template<typename T, int N>
  struct tuple_size<blitz::TinyVector<T,N>> : std::integral_constant<std::size_t,N> {};
  //
  template<std::size_t I, typename T, int N>
  struct tuple_element<I, blitz::TinyVector<T,N>> { using type = T; };
}
*/


inline double int_newt_cotes(const Array<double,1>& f, double a, double b, int N) noexcept
{
  // Implements Newton–Cotes integration rules when f is given on a discrete uniform mesh.
  //
  //                 Integrate[ f(x)dx , {x,a,b}]
  // 
  // User give number of points to be used, and interval [a,b], and we use high level 
  // routine to get as precise answer as possible for small price. Should be better than simpson
  double dh = (N>1) ? (b-a)/(N-1) : (b-a);
  if (N>=8){// O(dh^4)
    const double coef[]={17/48.,59/48.,43/48.,49./48.,1.0};
    double dsum=0;
    dsum += coef[0]*(f(0)+f(N-1)) + coef[1]*(f(1)+f(N-2)) + coef[2]*(f(2)+f(N-3)) + coef[3]*(f(3)+f(N-4));
    for (int i=4; i<N-4; ++i) dsum += f(i);
    return dsum*dh;
  }else{
    switch (N){
    case 1:
      return dh*f(0);
    case 2:
      return 0.5*dh*(f(0)+f(1));
    case 3:
      return (f(0)+4*f(1)+f(2))*dh/3.; // O(h^5)
    case 4:
      return (f(0)+3*(f(1)+f(2))+f(3))*3./8*dh; //O(h^5)
    case 5:
      return (7*f(0)+32*f(1)+12*f(2)+32*f(3)+7*f(4))*2./45.*dh; //O(h^7)
    case 6:
      return (19*f(0)+75*f(1)+50*(f(2)+f(3))+75*f(4)+19*f(5))*5./288.*dh; // O(h^7)
    case 7:
      return (41*f(0)+216*f(1)+27*f(2)+272*f(3)+27*f(4)+216*f(5)+41*f(6))/140*dh; //O(h^9)
    }
  }
  return 0;
}

struct BlockJob{
    int id;        // global index or whatever you need later
    int n;         // matrix dimension  (1,2,4,6,8,12,…)
    double cost;   // weight ∝ n²
};
ostream& operator<<(ostream& out, const BlockJob& v)
{
  out << "(" << "id="<<v.id<<",n="<<v.n<<",cost="<<v.cost<<")";
  return out;
}

/*
template<class T>
inline bool operator==(const blitz::TinyVector<T,3>& a, const blitz::TinyVector<T,3>& b){
    return (a[0] == b[0]) && (a[1] == b[1]) && (a[2] == b[2]);
}
*/

//namespace detail {
//  // equality (needed by std::unordered_map)
//  inline bool operator==(const TinyVector<int,3>& a, const TinyVector<int,3>& b) noexcept
//  { return a(0)==b(0) && a(1)==b(1) && a(2)==b(2);}
//
//// hash specialisation
//} // namespace detail
//
namespace std {
  struct TV4Hash
  {
    std::size_t operator()(const blitz::TinyVector<int,4>& v) const noexcept
    {
      // 64-bit FNV-1a magic constants
      uint64_t h = 14695981039346656037ull;          // offset basis
      constexpr uint64_t prime = 1099511628211ull;   // FNV prime
      h ^= static_cast<uint32_t>(v(0));  h *= prime;
      h ^= static_cast<uint32_t>(v(1));  h *= prime;
      h ^= static_cast<uint32_t>(v(2));  h *= prime;
      h ^= static_cast<uint32_t>(v(3));  h *= prime;
      return static_cast<std::size_t>(h);
    }
  };
  struct TV4Equal {
    bool operator()(const blitz::TinyVector<int,4>& a, const blitz::TinyVector<int,4>& b) const noexcept
    { return a(0)==b(0) && a(1)==b(1) && a(2)==b(2) && a(3)==b(3); }
  };
}

bool is_close(double phi1, double phi2, double rel_tol = 1e-10, double abs_tol = 1e-12) {
  //return std::abs(phi1 - phi2) <= std::max(rel_tol * std::max(std::abs(phi1), std::abs(phi2)), abs_tol);
  if (phi1==0) return phi2==0;
  return fabs(phi2/phi1-1)<1e-5;
}

string to_binary(const Array<int,1>& xfsgn)
{
  std::string s = "[";
  for (int i = 0; i < xfsgn.extent(0); ++i){
    s += fmt::format("{:b}", xfsgn(i));          // binary text
    if (i+1 < xfsgn.extent(0)) s += ", ";
  }
  s += "]";
  return s;
}

inline blitz::Range SRange(int lo, int hi, const std::string& msg) {
    if (hi < lo) {
        // Blitz doesn’t support empty slices; handle upstream.
        throw std::runtime_error("empty range "+msg);
    }
    return blitz::Range(lo, hi);
}


#endif// _BASICS
