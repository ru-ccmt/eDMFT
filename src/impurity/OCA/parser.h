// @Copyright 2007 Kristjan Haule
// 
#include <cstdio>
#include <cstdarg>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <cmath>

#ifndef _PARSER_
#define _PARSER_

class Entity{
protected:
  std::string Name;
  std::string Desc;
  bool thisParSet;
public:
  Entity(const std::string& name, const std::string& desc="") : thisParSet(false)
  {
    Name = name;
    Desc = desc;
  }
  virtual void print(std::ostream& stream, int m = 1) = 0;
  virtual void printl(std::ostream& stream) = 0;
  std::string name(){return Name;}
  virtual void pars(std::istringstream& input) = 0;
  bool IsSet() const { return thisParSet;}
  virtual ~Entity(){};
};

template <class T>
class Par : public Entity{
  // type=0 : b=10.0
  // type=1 : b={0,1,3}  : ar_n=3 ar_x={0,1,3}
  // type=2 : b={0-10,2} : dx=2 x0=0 x1=10
  // type=3 : b={0-10:4} : nx=4 x0=0 x1=10
  T x, x0, x1, dx;
  int type;
  int nx, ix;
  std::list<T> ar_x;
  typename std::list<T>::iterator ar_p;
public:
  Par(const std::string& name, const std::string& desc="") : Entity(name,desc), x(0), x0(-1), type(0), ix(0), ar_p(ar_x.begin()) {}
  Par(const std::string& name, const T& x_, const std::string& desc="") : Entity(name,desc), x(x_), x0(x_), type(0), ix(0), ar_p(ar_x.begin()) {}
  Par& operator =  (const T& x_) {x = x_; return *this;}
  operator T() const {return x;} // implicit conversion to type T
  virtual void print(std::ostream& stream, int m = 1);
  virtual void printl(std::ostream& stream);
  virtual void pars(std::istringstream& input);
  bool next();
  bool next(T& v) { v=x; return next();}
  T LowerBound() {
    switch (type){
    case 0 : return x;
    case 1 : return ar_x.front();
    case 2 : return x0;
    case 3 : return x0;
    default : return 0;
    };
  }
  T UpperBound(){
    switch (type){
    case 0 : return x;
    case 1 : return ar_x.back();
    case 2 : return x1;
    case 3 : return x1;
    default : return 0;
    };
  }
  void reset (){ ix=0;}
  void Set (const std::string& line);
  void Set0 (const T& x_)
  {type=0; x=x0=x_; ix=0;}
  void Set2 (const T& x0_, const T& x1_, const T& dx_)
  {type=2; x0=x0_; x1=x1_; dx=dx_;ix=0;}
  void Set3 (const T& x0_, const T& x1_, int nx_)
  {type=3; x0=x0_; x1=x1_; nx=nx_; ix=0;}
  //  int ar_num () {return ar_n;}
  bool IsDefault() {return x==x0;}
};

class StringPar : public Entity{
  std::string str;
public:
  StringPar(const std::string& name, const std::string& DefVal): Entity(name), str(DefVal){};
  StringPar(const std::string& name, const std::string& DefVal, const std::string& descr):Entity(name,descr), str(DefVal){};
  operator std::string() const {return str;} // Implicit conversion to char*
  virtual void print(std::ostream& stream, int m=1);
  virtual void printl(std::ostream& stream);
  virtual void pars(std::istringstream& input);
  bool operator== (const std::string& str_) {return str == str_;}
  StringPar& operator=(const std::string& str_){ str = str_; return *this;}
};

class DynArray{
  Entity** ar;
  int N0;
public:
  DynArray(int N_);
  DynArray(int N_, Entity* p0, ...);
  void Add(Entity* p0 ...);
  Entity*& operator[] (int i){return ar[i];}
  void ParsComLine(int argc, char *argv[]);
  void ParsEntity(const std::string& line);
  void ParsFile(const std::string& filename);
  void print(std::ostream& stream);
  void print();
  void printl(std::ostream& stream);
  void printl();
  void prints(std::ostream& stream);
  void prints();
  ~DynArray();
};


template <class T>
inline void Par<T>::Set (const std::string& line){
  std::istringstream input(line);
  pars(input);
}

template <class T>
inline bool Par<T>::next() {
  switch (type) {
  case 0 :
    if (ix<1) {ix++; return true;}
    else return false;
  case 1 :
    if (ar_p!=ar_x.end()) {
      x = *ar_p;
      ar_p++;
      return true;
    } else return false;
  case 2 :
    if (fabs(x0+ix*dx-x1)>0.5*fabs(dx)){
      x = x0+(ix++)*dx;
      return true;
    } else return false;
  case 3 :
    if (ix<=nx){
      x = x0 + (ix++)*(T)(x1-x0)/nx;
      return true;
    } else return false;
  default : return false;
  }
}

template <class T>
inline void Par<T>::pars(std::istringstream& input){
  thisParSet = true;
  char tch;
  input>>tch;
  if (tch != '['){
    input.putback(tch);
    type = 0;
    input >> x;
  } else {
    input >> x0;
    x = x0;
    char tch1;
    input>>tch1;
    
    if (tch1==','){
      input.putback(tch1);
      type = 1;
      ar_x.push_back(x0);
      while (input.get() != ']' && input.good()){
	T t; input>>t;
	ar_x.push_back(t);
      }
      ar_p = ar_x.begin();
    } else if (tch1=='-'){
      input >> x1;
      char tch2;
      input>>tch2;
      if (tch2==','){
	type = 2;
	input >> dx;
      } else if (tch2==':'){
	type = 3;
	input >> nx;
      }
    }
  }
}

template <class T>
inline void Par<T>::print(std::ostream& stream, int m){
  if (m) stream << "\t";
  stream << Name;
  if (m) stream << " = " ;
  else stream << "=";
  switch (type)  {
  case 0 : stream << x; break;
  case 1 :
    stream << '[';
    for (typename std::list<T>::iterator i=ar_x.begin(); i!=ar_x.end(); i++){
      stream << *i;
      if ((++i)-- != ar_x.end()) stream << ',';
      else stream << ']';
    }
    break;
  case 2 : stream << '[' << x0 << "-" << x1 << "," << dx <<']'; break;
  case 3 : stream << '[' << x0 << "-" << x1 << ":" << nx <<']'; break;
  }
}

template <class T>
inline void Par<T>::printl(std::ostream& stream){
  print(stream);
  stream << " " << Desc;
}

inline void StringPar::print(std::ostream& stream, int m){
  if (m) stream << "\t";
  stream << Name;
  if (m) stream << " = ";
  else stream << "=";
  stream << str;
}

inline void StringPar::printl(std::ostream& stream){
  print(stream);
  stream << " " << Desc;
}

inline void StringPar::pars(std::istringstream& input){
  using namespace std;
  thisParSet = true;
  input>>str;

  if (str.find("\"")!=string::npos)
    str = str.substr(1,str.length()-2);
  if (str.find("\'")!=string::npos)
    str = str.substr(1,str.length()-2);
  //cout<<"str="<<str<<endl;
}

inline void DynArray::ParsEntity(const std::string& line){
  std::string::size_type en;
  if ((en=line.find('=')) != std::string::npos){
    std::string::size_type start  = line.find_first_not_of(" \t");
    std::string::size_type end    = line.find_last_not_of(" \t", en-1);
    std::string var = line.substr(start,end-start+1);
    std::istringstream input(line.substr(en+1,std::string::npos));
    int i = -1;
    //    std::clog<<" In line :"<<line<<": found en="<<en<<" start="<<start<<" end="<<end<<" var="<<var<<";"<<std::endl;
    while (ar[++i] != NULL && i<N0){
      if(ar[i]->name() == var){
	ar[i]->pars(input);
	break;
      }
    }
  }
}

inline void DynArray::ParsFile(const std::string& FileName){
  std::ifstream input(FileName.c_str());
  if (input.fail()) {std::cerr<<"Invalid input filename "<<FileName<<" in parser.h\n"; return;}
  std::string line;
  while (getline(input,line,'\n')){
    // removing comments
    std::string::size_type comp = line.find_first_of('#');
    if (comp<std::string::npos) line.erase(comp);
    if (line.length()!=0) ParsEntity (line);
  }
  input.close();
}

inline DynArray::DynArray(int N_){
  N0 = N_;
  ar = new Entity* [N0];
  for (int i = 0; i < N0; i++) ar[i]=NULL;
}

inline DynArray::DynArray(int N_, Entity* p0 ...){
    N0 = N_;
    ar = new Entity* [N0];
    int i;
    for (i = 0; i < N0; i++) ar[i]=NULL;
    
    va_list ap;
    va_start (ap, p0);
    
    if (i!=0) ar[0] = p0;
    for (i=1; i < N0; i++){
	Entity* p = va_arg (ap, Entity*);
	if (p==NULL) break;
	ar[i] = p;
    }    
    va_end(ap);
}

inline void DynArray::Add(Entity* p0 ...){
    va_list ap;
    va_start (ap, p0);
    int i = 0;
    
    while (ar[i]!=NULL) i++;

    if (i<N0) ar[i++] = p0;
    for (int j=i; j < N0; j++){
	Entity* p = va_arg (ap, Entity*);
	if (p==NULL) break;
	ar[i] = p;
    }        
    va_end(ap);    
}

inline DynArray::~DynArray(){
    delete[] ar;
}

inline void DynArray::ParsComLine(int argc, char *argv[]){
    int first_arg = 1;
    std::ifstream input(argv[1]);
    if (!input.fail()){
	ParsFile (argv[1]);
	first_arg = 2;
    }
    for (int i = first_arg; i < argc; i++) ParsEntity(argv[i]);
}

inline void DynArray::print(std::ostream& stream){
    int i = 0;
    while (ar[i] != NULL && i < N0){
	ar[i++]->print(stream);
	stream << std::endl;
    }
}
inline void DynArray::print(){
  print(std::cout);
}

inline void DynArray::printl(std::ostream& stream){
    int i = 0;
    while (ar[i] != NULL && i < N0){
	ar[i++]->printl(stream);
	stream << std::endl;
    }
}
inline void DynArray::printl(){
    printl (std::cout);
}
inline void DynArray::prints(std::ostream& stream){
    int i = 0;
    stream << "# ";
    while (ar[i] != NULL && i < N0){
	if (i%6==0 && i!=0) stream << std::endl << "# ";
	ar[i++]->print(stream,0);
	stream << "\t";
    }
    stream << std::endl;
    stream.flush();
}

inline void DynArray::prints(){
    prints(std::cout);
}

/*
int main(int argc, char *argv[]){
  double u_, t_;
  int V_, S_;
  std::string InF;
  
  StringPar   InFi("InFile", "InputFile", "\t# Input Filename");
  Par<double> t("t", 2.0, "\t\t# Hopping amplitude");
  Par<int>    S("s", 3);
  Par<double> dob("dob", 5.0, "\t\t# dobavitev");
  Par<double> u("u", 2.0);
  Par<int>    V("V", 2, "\t\t# Ha Ha Ha");
  
  DynArray arr(6, &u, &V, &t, &S,  NULL);
  arr.Add(&dob, &InFi);
  
  arr.ParsComLine(argc, argv);
  arr.print();
  
  t.Set("{0, 1, 2, 3}");
  u.Set("{0-15,4}");
  V.Set("{0-15:5}");
  InFi = "Ne parkiraj!";
  
  arr.printl();
  
  InF = InFi;
  t_ = t;
  S_ = S;
  u_ = u;
  V_ = V;
  
  std::cout << u_ << "," << V_ << "," << t_ << "," << S_ << "," << InF << "," << dob << std::endl;
  while (dob.next()){
    std::cout << "dob*2 : " << dob*2 << std::endl;
    dob = 4;
    std::cout << "dob*2 : " << dob*2 << std::endl;
  }
  return 0;
}
*/
#endif // _PARSER
