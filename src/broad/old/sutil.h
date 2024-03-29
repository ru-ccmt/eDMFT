// @Copyright 2007 Kristjan Haule
// 
#ifndef MY_UTILS
#define MY_UTILS
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>

#define PlacementNewOffset 8

#ifdef _DEBUG
#define Assert(condition, message)\
{\
  if(!(condition)) std::cerr << (message) << std::endl;\
}
#define Assert2(condition, message, message1)\
{\
  if(!(condition)) std::cerr << (message) << " "<< (message1)<< std::endl;\
}
#define Assert3(condition, message, message1, message2)\
{\
  if(!(condition)) std::cerr << (message) << " " << (message1) << " "<< message2 << std::endl;\
}
#define Assert5(condition, message, message1, message2, message3, message4)\
{\
  if(!(condition)) std::cerr <<(message)<<" "<<(message1)<<" "<<(message2)<<" "<<(message3)<<" "<<(message4)<< std::endl;\
}
#define _LOG(x) x
#define CHECK(x) x
#else /* NO_ARG_CHECK */
#define Assert(condition, message)
#define Assert2(condition, message, message1)
#define Assert3(condition, message, message1, message2)
#define Assert5(condition, message, message1, message2, message3, message4)
#define _LOG(x)
#define CHECK(x)
#endif /* NO_ARG_CHECK */
#define HPoffset 8

//************************************************//
// Simple class used to exchange the data between //
// functions and meshes when doing linear	  //
// interpolation  				  //
//************************************************//
class intpar{
public:
  int i;
  double p;
  intpar(int i_, double p_) : i(i_), p(p_){}
  intpar(){};
};

//template <class T>
//inline T sqr(const T& x) {return x*x;}

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

inline double ferm_f(double x)
{
  if (x>100.) return 0.0;
  if (x<-36.) return 1.0;
  return 1.0/(1+exp(x));
}
inline double nbose_f(double x)
{
  if (x>100.) return 0.0;
  if (x<-36.) return -1.0;
  return 1.0/(exp(x)-1);
}
inline double exp_f(double x){
  if (x<-100) return 0;
  return exp(x);
}

template <class T>
double power(T p, int k)
{
  T s = 1;
  for (int i=0; i<k; i++) s*=p;
  return s;
}
inline std::string NameOfFile(const std::string& name, int n, int m=3){
  std::stringstream str;
  str<<name<<".";
  for (int i=1; i<m; i++){
    if (n<power(10,i)) str<<"0";
  }
  str<<n<<std::ends;
  return std::string(str.str());
}

inline std::string NameOfFile_(const std::string& name, int n1, int n2, int m1=2, int m2=3){
  std::stringstream str;
  str<<name<<".";
  for (int i=1; i<m1; i++){
    if (n1<power(10,i)) str<<"0";
  }
  str<<n1;
  str<<".";
  for (int i=1; i<m2; i++){
    if (n2<power(10,i)) str<<"0";
  }
  str<<n2<<std::ends;
  return std::string(str.str());
}
#define RED           "\033[31;1m"
#define GREEN         "\033[32;1m"
#define YELLOW        "\033[33;1m"
#define BLUE          "\033[34;1m"
#define PURPLE        "\033[35;1m"
#define CYAN          "\033[36;1m"
#define WHITE         "\033[37;1m"
#define URED          "\033[31;1;4m"
#define UGREEN        "\033[32;1;4m"
#define UYELLOW       "\033[33;1;4m"
#define UBLUE         "\033[34;1;4m"
#define UPURPLE       "\033[35;1;4m"
#define UCYAN         "\033[36;1;4m"
#define UWHITE        "\033[37;1;4m"
#define DEFAULTCOLOR "\033[0m"
 
#ifndef _PARALEL_
#define COLOR(color, arg) color << arg << DEFAULTCOLOR
#else
#define COLOR(color, arg) arg
#endif

#endif // MY_UTILS
