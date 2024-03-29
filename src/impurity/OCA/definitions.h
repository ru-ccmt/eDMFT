// @Copyright 2007 Kristjan Haule
// 
#ifndef _DEFINITIONS_
#define _DEFINITIONS_

class diag{
  int v[3];
public:
  explicit diag(int v0=-1, int v1=-1, int v2=-1)
  {
    v[0]=v0; v[1]=v1; v[2]=v2;
  }
  int& operator[](int i) {return v[i];}
  int operator[](int i) const {return v[i];}
  friend ostream& operator<<(ostream& stream, const diag& d);
  bool exist() const {return v[0]!=-1 || v[1]!= -1 || v[2]!=-1;}
};
ostream& operator<<(ostream& stream, const diag& d)
{
  stream<<setw(2)<<d[0]<<" "<<setw(2)<<d[1]<<" "<<setw(2)<<d[2];
  return stream;
}

class OCAd{
public:
  vector<int> states;
  int ib1, ib2;
  double f;
  OCAd(int i0=-1, int i1=-1, int i2=-1, int i3=-1, int ib1_=-1, int ib2_=-1, double f_=0) : states(4), ib1(ib1_), ib2(ib2_), f(f_)
  { states[0] = i0;    states[1] = i1;    states[2] = i2;    states[3] = i3;  }
  bool operator==(const OCAd& oc)
  { return (states[0]==oc.states[0] && states[1]==oc.states[1] && states[2]==oc.states[2] && states[3]==oc.states[3] && ib1==oc.ib1 && ib2==oc.ib2);}
};
std::ostream& operator<< (std::ostream& stream, const OCAd& r){
  stream<<setw(3)<<r.states[0]<<" "<<setw(3)<<r.states[1]<<" "<<setw(3)<<r.states[2]<<" "<<setw(3)<<r.states[3]<<"    ";
  stream<<setw(3)<<r.ib1<<" "<<setw(3)<<r.ib2<<"    "<<setw(20)<<r.f<<" ";
  return stream;
}
std::istream& operator>> (std::istream& stream, OCAd& r){
  for (int i=0; i<r.states.size(); i++) stream >> r.states[i];
  stream >> r.ib1 >> r.ib2 >> r.f;
  return stream;
}
inline double product(const double* A, const double* G, int size)
{
  double sum = 0;
  for (int i=0; i<size; i++) sum += A[i]*G[i];
  return sum;
}

#endif
