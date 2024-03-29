// @Copyright 2007 Kristjan Haule
// 
#ifndef MESH_
#define MESH_
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include "util.h"

/////////////////////// Hierarchy of mashes: ////////////////////////
//                             
//                      base clas  = mesh
//              |                                |
//           mesh1D                          meshProxy
//
//                            mesh2D

//*****************************************//
// Classes implemented in this header file //
//*****************************************//
class mesh;
class mesh1D;
class meshProxy;
class mesh2D;

//************************************************//
// Simple class used to exchange the data between //
// functions and meshes when doing linear	  //
// interpolation  				  //
//************************************************//
typedef int tint;
class intpar;
//class intpari;

class parMeshFor2Peaks{
public:
  double cf1, cf2, peakf;
  int mejaf1_a, mejaf1_b, mejaf2_a, mejaf2_b;
};

//****************************************//
// Clas compare is used in AddPoints for  //
// sorting alhorithm from stl.	    	  //
//****************************************//
class compare{
  std::vector<double> const* u;
public:
  compare(const std::vector<double>& u_);
  bool operator()(int a, int b) const;
};

//******************************************************************************//
//  Basic class for meshes: mesh1D, mesh2D. It is also used as proxy class for	//
//  mesh2D to prevent calling constructor and destructor. Memory for 2D mesh	//
//  is allocated in the fortran style for efficiency. Constructors for proxy	//
//  classes (class mesh) are called with placement new operator.		//
//  Class mesh can find specified real value in the mesh (find_) and can	//
//  calculate relative distance between specified point and the nearest point	//
//  in the mesh (class intpar). 						//
//******************************************************************************//
class mesh{
public:
  int N, N0;
protected:
  double *om;
  double *delta;
  double *dh;
public:
  double& operator[](int i) {Assert(i<N,"Out of range in mesh[]"); return om[i];}
  const double& operator[](int i) const {Assert(i<N,"Out of range in mesh[]"); return om[i];}
  int size() const {return N;}
  double last() const {return om[N-1];}
  int find_(double x, int& ai0) const;
  int find_(double x) const;
  int _find(double x, int ai0, int& bi) const;
  tint InitInterpLeft() const {return 0;}
  tint InitInterpRight() const {return N-2;}
  intpar InterpLeft(const double x, tint& a) const;
  intpar InterpRight(const double x, tint& a) const;
  intpar Interp(const double x) const;
//    void InterpLeft(const double x, tint& a, intpari& pi) const;
//    void InterpLeft(const double x, tint& a, intpari& pi, intpari& pj) const;
//    void InterpLeft(const double x, tint& a, intpar& p, intpari& pi, intpari& pj) const;
//    void InterpLeft(const double x, tint& a, intpar& p, bintpar& pi) const;
  template <int left> tint InitInterp() const;
  template <int left> intpar Interp(const double x, tint& a) const;
  
  void InterpBoth(const double x, int& ix, intpar& ip) const;
  int findBoth(const double x, int ix) const;
  void SetUp();
  double* begin() const {return om;}
  double* end() const {return om+N;}
  friend std::ostream& operator<<(std::ostream& stream, const mesh& m);
  double Delta(int i) const {Assert(i<N,"Out of range in mesh.delta[]"); return delta[i];}
  double Dh(int i) const {Assert(i<N,"Out of range in mesh.dh[]"); return dh[i];}
  const double* Dh() const {return dh;}
  int gFinda(double x) const;
  int gFindb(double x, int a) const;
  int gFindb(double x, int a, int b) const;
  double* MemoryPointer() {return om;}
protected:
  mesh(): N(0),N0(0),om(NULL),delta(NULL),dh(NULL){};
  ~mesh(){};
  void mcopy(double* om_, double* delta_, double* dh_);
  void mcopy(const mesh& m);
private:
  static const int dN = 10;
  int bisection(double x, int& klo, int& khi, int bi) const;
  friend class mesh2D;
  friend class meshProxy;
};

//************************************************//
// One dimentional mesh derived from class mesh.  //
// In addition to the base class, it has a	  //
// constructor that allocates the needed memory	  //
// and also destructor to prevent resorce leakes. //
// It has also an integer (center) that holds	  //
// position of the peak in the mesh		  //
//************************************************//
class mesh1D : public mesh{
  int center;
public:
  mesh1D(const mesh1D& m);
  mesh1D(int N_);
  mesh1D(){};
  ~mesh1D();
  mesh1D& operator=(const mesh1D& m);
  void resize(int n);
  void SetUp(int center_);
  void SetUp(double center_);
  void MakeCentralMesh(const mesh1D& m);
  double dcenter() const;
  int icenter() const;
  void MakeMesh(int N, int N1, int N2, double dcenter, double tanc,
		double tanw, double a0, double b0, double b1);
  void MakeMesh(int N_, double tanc, double tanw, double b0, double b1);
  void MakeMesh(int N_, double start, double end);
  //  void CopyAndAddPoints(double centerOm, const mesh1D& ome, const mesh1D& omc, int M, int* ind_1);
  void CopyAndAddPoints(double centerOm, const mesh1D& ome, const mesh1D& omc, int M, int iremove);
  template <class meshy, class Tind>
  void CopyAndAddPoints(double centerOm, const meshy& ome, const mesh1D& omc, int M, Tind& ind_1);
  void MeshFor2Peaks_Intervals(double peakOm, parMeshFor2Peaks& par) const;
  
  friend void MeshFor2Peaks_SmallMesh(mesh1D& m, int MaxPoints, const mesh1D& ome, const parMeshFor2Peaks& par);
};

class meshProxy : public mesh{
public:
  typedef mesh1D Cousin;
  void Initialize(int N_, double* om_, double* delta_, double* dh_);
  void Initialize(const mesh& m);
  void resize(int n);
  meshProxy& operator=(const meshProxy& m);
  ~meshProxy(){};
  void SetUp(double);//{mesh::SetUp();}
  double dcenter() const { return 0.0;}
};

//****************************************************************//
// Two dimentional mesh stored like and array of class mesh.	  //
// It is allocated in a fortran-like manner. It also holds	  //
// a pointer to the one dimentional array (Om) that represents	  //
// the value of the mesh in a distinct row. So the size of the	  //
// 1D mesh must be equal to the number of rows, while the number  //
// of columns can wory from one row to the another.		  //
// There is also a two dimentional array (ind_1) that holds the	  //
// information about the original position of points (in 2D mesh) //
// according to the parent 1D mesh.				  //
//****************************************************************//
class mesh2D{
private:
  int N, N0, Nd, Nd0;
  meshProxy* r;
  mesh1D const* Om;
  int **ind_1;
public:
  mesh2D() : N(0), N0(0), Nd(0), Nd0(0), ind_1(NULL), memory(NULL) {};
  mesh2D(int N_, int Nd_);
  ~mesh2D();
  void resize(int N_, int Nd_);
  meshProxy& operator[](int i) {Assert(i<N,"Out of range in mesh2D[]"); return r[i];}
  const meshProxy& operator[](int i) const {Assert(i<N,"Out of range in mesh2D[]"); return r[i];}
  void SetUp(const mesh1D& om, const mesh1D& Om, const mesh1D& omc, int M, double peakpos);
  const int * getInd_1(int k) {Assert(k<N,"Out of range in mesh2D.ind_1[]"); return ind_1[k];}
  int size_N() { return N;}
  int size_Nd() {return Nd;}
private:  
  void *memory;
  mesh2D(const mesh2D&){};
};

template <class meshx> void MeshFor2Peaks_FullMesh(meshx& m, double peakOm, const mesh1D& ome, int M);
template <class meshx> void MeshFor2Peaks_FullMesh(meshx& m, const mesh1D& ome, const parMeshFor2Peaks& par);
template <class meshx> void MeshFor2Peaks_SmallMesh(meshx& m, double peakOm, const mesh1D& ome, int MaxPoints);
template <class meshx> void MeshFor2Peaks_SmallMesh(meshx& m, int MaxPoints, const mesh1D& ome, const parMeshFor2Peaks& par);

// mesh ////////////////////////////////////////////////////////////////////
inline void mesh::mcopy(double* om_, double* delta_, double* dh_)
{
  std::copy(om_,om_+N,om);
  std::copy(delta_,delta_+N,delta);
  std::copy(dh_,dh_+N,dh);
}

inline void mesh::mcopy(const mesh& m)
{
  std::copy(m.om,   m.om+N,   om);
  std::copy(m.delta,m.delta+N,delta);
  std::copy(m.dh,   m.dh+N,   dh);
}

inline int mesh::bisection(double x, int& klo, int& khi, int bi) const
{
  int k;
  khi=bi-1;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (om[k] > x) khi=k;
    else klo=k;
  }
  return klo;
}

inline int mesh::find_(double x, int& ai0) const
{
  if (x<om[ai0+1]) return ai0;
  int ai1 = ai0+1;
  if (ai1+1<N) {
    if (x<om[ai1+1]) {
      ai0 = ai1;
      ai1 = ai1+1;
      return ai0;
    }
  } else return ai0;
  if (ai1+1<N-1) ai0 = ai1+1;
  else return ai1;
  if (ai1+dN<N){
    if (x<om[ai1+dN]) return bisection (x, ai0, ai1, ai1+dN+1);
  } else return bisection (x, ai0, ai1, N);
  if (ai1+dN<N-1) ai0 = ai1+dN;
  else ai0 = ai1+dN-1;
  return bisection (x, ai0, ai1, N);
}

inline int mesh::find_(double x) const
{
  int ai0=0;
  return find_(x,ai0);
}

inline int mesh::_find(double x, int ai0, int& bi) const
{
  if (x>om[bi]) return bi;
  if (x>om[bi-1]) return bi-1;
  if (bi-2>ai0){
    if (x>om[bi-2]) return bi-2;
  } else return ai0;
  if (bi-dN>ai0){
    if (x>om[bi-dN]) {
      int ai1, ai00;
      ai00 = bi-dN;
      return bisection(x,ai00,ai1,bi-1); 
    }
    bi=bi-dN+1;
  } else bi-=1; 
  {
    int ai1;
    return bisection (x,ai0,ai1,bi);
  }
}

inline int mesh::gFinda(double x) const{
  if (x<=om[0]) return 0;
  int k, klo=0;
  int khi=N-1;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (om[k] > x) khi=k;
    else klo=k;
  }
  return khi;
}

inline int mesh::gFindb(double x, int a = 0) const{
  if (x>=om[N-1]) return N;
  int k, klo=(a!=0) ? a-1 : a;
  int khi=N-1;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (om[k] > x) khi=k;
    else klo=k;
  }
  return khi;
}

inline int mesh::gFindb(double x, int a, int b) const{
  if (x>=om[b-1]) return b;
  int k, klo=(a!=0) ? a-1 : a;
  int khi=b-1;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (om[k] > x) khi=k;
    else klo=k;
  }
  return khi;
}

inline int mesh::findBoth(const double x, int ix) const
{
  if (x==om[ix]) return ix;
  if (x>om[ix]) return find_(x,ix);
  else return _find(x,0,ix);
}

inline void mesh::InterpBoth(const double x, int& ix, intpar& ip) const
{
  int i;
  if (ix>N-1) {ix = N-1; i = _find(x,0,ix);}
  else if (x>om[ix]) i = find_(x,ix);
  else if (x<om[ix]) i = _find(x,0,ix);
  else {
    ip.i = ix;
    ip.p = 0;
    return;
  }
  ip.i = ix = i;
  ip.p = (x-om[i])*delta[i];
  SUNCA_LOG(if ((ip.p>1.0 || ip.p<0) && i!=0 && i!=(N-2)) std::cerr<<"Variable p="<<ip.p<<" at i="<<i<<" flag1="<<debug_flag1<<" flag2="<<debug_flag2<<std::endl);
}

inline intpar mesh::InterpLeft(const double x, tint& a) const
{
  int i = find_(x,a);
#ifdef SUNCA_DEBUG
  double p = (x-om[i])*delta[i];
  if ((p>1.0 || p<0) && i!=0 && i!=(N-2)) {
    std::cerr<<"Variable p="<<p<<" at i="<<i<<" flag1="<<debug_flag1<<" flag2="<<debug_flag2<<std::endl;
  }
  return intpar(i,p);
#else
  return intpar(i,(x-om[i])*delta[i]);
#endif
}

inline intpar mesh::InterpRight(const double x, tint& a) const
{
  int i = _find(x,0,a);
#ifdef SUNCA_DEBUG
  double p = (x-om[i])*delta[i];
  if ((p>1.0 || p<0) && i!=0 && i!=(N-2)) {
    std::cerr<<"Variable p="<<p<<" at i="<<i<<" flag1="<<debug_flag1<<" flag2="<<debug_flag2<<std::endl;
  }
  return intpar(i,p);
#else
  return intpar(i,(x-om[i])*delta[i]);
#endif
}

inline intpar mesh::Interp(const double x) const
{
  tint a = InitInterpLeft();
  return InterpLeft(x,a);
}

template <int left>
tint mesh::InitInterp() const
{ std::cerr << "Bug!! InterpLeft<> should not be called"<<std::endl;}
template <int left>
intpar mesh::Interp(const double x, tint& a) const
{ std::cerr << "Bug!! InterpLeft<> should not be called"<<std::endl;}

template <>
inline tint mesh::InitInterp<1>() const
{ return InitInterpLeft();}
template <>
inline intpar mesh::Interp<1>(const double x, tint& a) const
{ return InterpLeft(x, a);}
template <>
inline tint mesh::InitInterp<-1>() const
{ return InitInterpRight();}
template <>
inline intpar mesh::Interp<-1>(const double x, tint& a) const
{ return InterpRight(x, a);}

inline void mesh::SetUp()
{
  delta[0] = 1/(om[1]-om[0]);
  dh[0] = 0.5*(om[1]-om[0]);
  for (int i=1; i<N-1; i++){
    delta[i] = 1/(om[i+1]-om[i]);
    dh[i] = 0.5*(om[i+1]-om[i-1]);
  }
 delta[N-1] = 0.0;
 dh[N-1] = 0.5*(om[N-1]-om[N-2]);
}

inline std::ostream& operator<<(std::ostream& stream, const mesh& m)
{
  int width = stream.width();
  for (int i=0; i<m.N; i++)
    stream << std::setw(3) << i << std::setw(width) << m[i] << std::endl;
  return stream;
}

// mesh1D //////////////////////////////////////////////////////////////////
inline mesh1D::mesh1D(int N_)
{
  N0=N=N_;
  om = new double[N0];
  delta = new double[N0];
  dh = new double[N0];
}

inline mesh1D::~mesh1D()
{
  delete[] om;
  delete[] delta;
  delete[] dh;
  om = NULL;
  delta = NULL;
  dh = NULL;
  N0=N=0;
}

inline void mesh1D::resize(int n)
{
  if (n>N0){
    if (om) delete[] om;
    if (delta) delete[] delta;
    if (dh) delete[] dh;
    om = new double[n];
    delta = new double[n];
    dh = new double[n];
    N0 = n;
  }
  N = n;
}

inline mesh1D::mesh1D(const mesh1D& m)
{
  resize(m.N);
  mcopy(m);
  center = m.center;
}

inline mesh1D& mesh1D::operator=(const mesh1D& m)
{
  resize(m.N);
  mcopy(m);
  center = m.center;
  return *this;
}

inline void mesh1D::SetUp(int center_)
{
  mesh::SetUp();
  center = center_;
}

inline void mesh1D::SetUp(double center_)
{
  mesh::SetUp();
  center = find_(center_);
}

inline void mesh1D::MakeCentralMesh(const mesh1D& m)
{
  resize(m.N-m.center-1);
  mcopy(m.om+m.center+1,m.delta+m.center+1,m.dh+m.center+1);
  double center = 0.5*(m.om[m.center]+m.om[m.center+1]);
  for (int i=0; i<N; i++) om[i] -= center;
  center = 0;
}

inline double mesh1D::dcenter() const
{
  return 0.5*(om[center]+om[center+1]);
}

inline int mesh1D::icenter() const
{
  return center;
}

inline void mesh1D::MakeMesh(int N_, int N1, int N2, double dcenter, double tanc,
			     double tanw, double a0, double b0, double b1){
  Assert(N1>=0,"N1 must be positive!");
  Assert(N2>0,"N2 must be positive!");
  Assert(a0<b0 && b0<b1, "Relation must hold: a0<b0<b1!");
  Assert(b0<tanw && tanw<b1,"Relation mesu hold: b0<tanw<b1!");
  Assert(a0>0, "a0 must be positive!");
  Assert(b0>0, "b0 must be positive!");
  Assert(-b1<dcenter || dcenter<b1, "Relation must hold -b1<dcenter<b1!");
  Assert(N_>=N1+N2, "Relation N>N1+N2 must hold!");
  resize(N_);
  const int Np = (N+N1)/2, N3_2=(N-N1-N2)/2, N1_2=N1/2, N2_2=N2/2;
  const double du = atan(((tanc-b0)/tanw)), b1n = atan((b1-tanc)/tanw)+du;
  double a0l_ = log (a0), b0l = log(b0);
  for (int i=0; i<N3_2; i++){
    om[Np+i] = exp(a0l_+(b0l-a0l_)*i/(N3_2-1));
  }
  for (int i=0; i<N1_2; i++){
    om[N/2+i] = a0*(-1+2.0*(i+1+N1_2)/static_cast<double>(N1+1));
  }
  for (int i=0; i<N2_2; i++) om[i+N-N2_2] = tanc+tanw*tan(b1n*(i+1)/N2_2-du);
  for (int i=0; i<N/2; i++) om[i] = -om[N-i-1];
  for (int i=0; i<N; i++) om[i] += dcenter;
  SetUp(N/2-1);
}

inline void mesh1D::MakeMesh(int N_, double tanc, double tanw, double b0, double b1){
  Assert(N_>=0,"N must be positive!");
  Assert(b0<b1, "Relation must hold: b0<b1!");
  Assert(b0<tanw && tanw<b1,"Relation mesu hold: b0<tanw<b1!");
  Assert(b0>0, "b0 must be positive!");
  resize(N_);
  const double du = atan(((tanc-b0)/tanw)), b1n = atan((b1-tanc)/tanw)+du;
  om[0]=0.0;
  for (int i=1; i<N; i++) om[i] = tanc+tanw*tan(b1n*(i-1)/(N-2)-du);
  SetUp(0);
}

inline void mesh1D::MakeMesh(int N_, double start, double end){
  resize(N_);
  double x=start, dh=(end-start)/(N-1.);
  for (int i=0; i<N; i++,x+=dh)
    om[i] = x;
  SetUp(N/2);
}

// compare /////////////////////////////////////////////////////////
inline compare::compare(const std::vector<double>& u_)
{
  u = &u_;
}

inline bool compare::operator()(int a, int b) const
{
  if ((*u)[a] < (*u)[b]) return true;
  else return false;
};

template <class meshx, class meshy, class Tind>
void CopyAndAddPoints(meshx& m, double centerOm, const meshy& ome, const mesh1D& omc, int M, Tind& ind_1){
/*************************************************************************
* Points are added in the neighbourhut of Om. They start for
* offset apart from Om. The number of added points is M
* (M/2 points are added to the left and M/2 to the right side).
* If Om is close to zero than every point from parent mesh is added
* otherwise each second point is added. Coefficients for integration
* are computed in the end.
***************************************************************************/
  const int M2=M/2; M=2*M2;
  const int Nd = ome.size()+M;

  m.resize(Nd);
  // Adding points in array omadd
  std::vector<double> omadd(Nd);
  int ixup, ixdo;
  if (fabs(centerOm-ome.dcenter())<2*omc[M2-1]){
    // every points from omc is added
    ixup = ome.find_(centerOm+omc[M2-1])+2; // after last point added
    ixdo = ome.find_(centerOm-omc[M2-1])-1; // before first point added
    for (int i=0; i<M2; i++){
      omadd[i+M2] = centerOm+omc[i];
      omadd[M2-i-1] = centerOm-omc[i];
    }
  } else {
    // only each second point from parent mesh is added
    ixup = ome.find_(centerOm+omc[2*(M2-1)])+2;
    ixdo = ome.find_(centerOm-omc[2*(M2-1)])-1;
    for (int i=0; i<M2; i++) {
      omadd[i+M2] = centerOm+omc[2*i];
      omadd[M2-i-1] = centerOm-omc[2*i];
    }
  }
  if (ixdo<0) ixdo=0;
  if (ixup>ome.size()-1) ixup = ome.size();
  
  // Setting up new mesh consisting of new and old points.
  // New points are just added to the end to preserve information
  // about parent greed. Array ind is calculated so that points can
  // be reached sequentially.
  std::vector<double> epst(Nd);
  std::copy(ome.begin(),ome.end(),epst.begin());
  std::copy(omadd.begin(),omadd.begin()+M,epst.begin()+ome.size());

  std::vector<int> ind(Nd);
  for (int i=0; i<ixup; i++) ind[i]=i;
  for (int i=ixup; i<ome.size(); i++) ind[M+i]=i;
  for (int i=0; i<M; i++) ind[ixup+i]=ome.size()+i;
  // sorting points
  compare comp(epst);
  std::sort(ind.begin()+ixdo,ind.begin()+ixup+M,comp);
  
  //  N = Nd;
  for (int i=0; i<Nd; i++){
    m[i] = epst[ind[i]];
    ind_1[ind[i]] = i;
  }
  int center = m.find_(ome.dcenter());
  // Calculating new integration coefficients and deltas
  m.SetUp(center);
}

template <class meshx>
void CopyAndAddPoints(meshx& m, double centerOm, const mesh1D& ome, const mesh1D& omc, int M, int iremove)
{/*************************************************************************
* Points are added in the neighbourhut of Om. They start for
* offset apart from Om. The number of added points is M
* (M/2 points are added to the left and M/2 to the right side).
* If Om is close to zero than every point from parent mesh is added
* otherwise each second point is added. Coefficients for integration
* are computed in the end.
***************************************************************************/
  const int M2=M/2; M=2*M2;
  const int Nd = (iremove>=0) ? ome.N+M-1 : ome.N+M;

  // Adding points in array omadd
  std::vector<double> omadd(M);
  int ixup, ixdo;
  if (fabs(centerOm-ome.dcenter())<2*omc[M2-1]){
    // every points from omc is added
    ixup = ome.find_(centerOm+omc[M2-1])+2; // after last point added
    ixdo = ome.find_(centerOm-omc[M2-1])-1; // before first point added
    for (int i=0; i<M2; i++){
      omadd[i+M2] = centerOm+omc[i];
      omadd[M2-i-1] = centerOm-omc[i];
    }
  } else {
    // only each second point from parent mesh is added
    ixup = ome.find_(centerOm+omc[2*(M2-1)])+2;
    ixdo = ome.find_(centerOm-omc[2*(M2-1)])-1;
    for (int i=0; i<M2; i++) {
      omadd[i+M2] = centerOm+omc[2*i];
      omadd[M2-i-1] = centerOm-omc[2*i];
    }
  }
  if (ixdo<0) ixdo=0;
  if (ixup>ome.N-1) ixup = ome.N;
  
  // Setting up new mesh consisting of new and old points.
  // New points are just added to the end to preserve information
  // about parent greed. Array ind is calculated so that points can
  // be reached sequentially.

  if (m.begin()==ome.begin()){
    mesh1D omt(ome);
    m.resize(Nd);
    if (iremove>=0){
      std::copy(omt.begin(), omt.begin()+iremove, m.begin());
      std::copy(omt.begin()+iremove+1, omt.end(), m.begin()+iremove);
      std::copy(omadd.begin(), omadd.begin()+M, m.begin()+omt.N-1);
    } else{
      std::copy(omt.begin(), omt.end(), m.begin());
      std::copy(omadd.begin(), omadd.begin()+M, m.begin()+omt.N);
    }
  } else {
    m.resize(Nd);
    if (iremove>=0){
      std::copy(ome.begin(), ome.begin()+iremove, m.begin());
      std::copy(ome.begin()+iremove+1, ome.end(), m.begin()+iremove);
      std::copy(omadd.begin(), omadd.begin()+M, m.begin()+ome.N-1);
    } else{
      std::copy(ome.begin(), ome.end(), m.begin());
      std::copy(omadd.begin(), omadd.begin()+M, m.begin()+ome.N);
    }
  }
  
  sort(m.begin(), m.end());
  int center = m.find_(ome.dcenter());
  // Calculating new integration coefficients and deltas
  m.SetUp(center);
}

inline void mesh1D::MeshFor2Peaks_Intervals(double peakOm, parMeshFor2Peaks& par) const
{
  par.peakf = 0.5*(om[center]+om[center+1]);
  if (par.peakf<peakOm){
    par.cf1 = par.peakf;
    par.cf2 = peakOm;
  } else{
    par.cf1 = peakOm;
    par.cf2 = par.peakf;
  }
  par.mejaf1_a = 0;
  par.mejaf1_b = gFindb(0.5*(par.cf2-par.cf1)+par.peakf);
  par.mejaf2_a = gFinda(-0.5*(par.cf2-par.cf1)+par.peakf);
  par.mejaf2_b = N;
  if (par.mejaf1_b<(size()-1) && om[par.mejaf2_a]+par.peakf>om[par.mejaf1_b]) par.mejaf1_b++;
}

template <class meshx>
inline void MeshFor2Peaks_FullMesh(meshx& m, const mesh1D& ome, const parMeshFor2Peaks& par){
  int N_ = par.mejaf1_b-par.mejaf1_a+par.mejaf2_b-par.mejaf2_a;
  m.resize(N_);
  
  int j=0;
  for (int i=par.mejaf1_a; i<par.mejaf1_b; i++) m[j++] = ome[i]-par.peakf+par.cf1;
  for (int i=par.mejaf2_a; i<par.mejaf2_b; i++) m[j++] = ome[i]-par.peakf+par.cf2;
  m.SetUp();
}

template <class meshx>
inline void MeshFor2Peaks_SmallMesh(meshx& m, int MaxPoints, const mesh1D& ome, const parMeshFor2Peaks& par){
  double p = (par.mejaf1_b-par.mejaf1_a+par.mejaf2_b-par.mejaf2_a)/static_cast<double>(MaxPoints);
  int mejaf1_ap, mejaf1_bp, mejaf2_ap, mejaf2_bp;
  double p1, p2, p3, p4;
  if (p>1.0){
    mejaf1_ap = static_cast<int>(ome.center-ome.center/p);
    p1 = (ome[par.mejaf1_a]-par.peakf)/(ome[mejaf1_ap]-par.peakf);
    mejaf1_bp = static_cast<int>(ome.center+(par.mejaf1_b-ome.center)/p);
    p2 = (0.5*(par.cf2-par.cf1)+par.peakf)/(ome[mejaf1_bp-1]-par.peakf);
    mejaf2_ap = static_cast<int>(ome.center-(ome.center-par.mejaf2_a)/p);
    p3 = (-0.5*(par.cf2-par.cf1)+par.peakf)/(ome[mejaf2_ap-1]-par.peakf);
    mejaf2_bp = static_cast<int>((par.mejaf2_b-ome.center)/p+ome.center);
    p4 = (ome[par.mejaf2_b-1]-par.peakf)/(ome[mejaf2_bp-1]-par.peakf);
  } else {
    p1 = p2 = p3 = p4 = 1.0;
    mejaf1_ap = par.mejaf1_a;
    mejaf1_bp = par.mejaf1_b; 
    mejaf2_ap = par.mejaf2_a;
    mejaf2_bp = par.mejaf2_b;
  }
  
  int N_ = mejaf1_bp-mejaf1_ap+mejaf2_bp-mejaf2_ap;
  m.resize(N_);
  
  int j=0;
  for (int i=mejaf1_ap; i<ome.center; i++) m[j++] = (ome[i]-par.peakf)*p1+par.cf1;
  for (int i=ome.center; i<mejaf1_bp; i++) m[j++] = (ome[i]-par.peakf)*p2+par.cf1;
  for (int i=mejaf2_ap; i<ome.center; i++) m[j++] = (ome[i]-par.peakf)*p3+par.cf2;
  for (int i=ome.center; i<mejaf2_bp; i++) m[j++] = (ome[i]-par.peakf)*p4+par.cf2;
  m.SetUp();
}

template <class meshx>
inline void MeshFor2Peaks_FullMesh(meshx& m, double peakOm, const mesh1D& ome, int M){
  parMeshFor2Peaks par;
  ome.MeshFor2Peaks_Intervals(peakOm, par);
  MeshFor2Peaks_FullMesh(m, ome, par);
}

template <class meshx>
inline void MeshFor2Peaks_SmallMesh(meshx& m, double peakOm, const mesh1D& ome, int M){
  parMeshFor2Peaks par;
  ome.MeshFor2Peaks_Intervals(peakOm, par);
  MeshFor2Peaks_SmallMesh(m, M, ome, par);
}

template <class meshy, class Tind>
inline void mesh1D::CopyAndAddPoints(double centerOm, const meshy& ome, const mesh1D& omc, int M, Tind& ind_1)
{
  resize(ome.size()+M);
  ind_1.resize(ome.size()+M);
  ::CopyAndAddPoints(*this, centerOm, ome, omc, M, ind_1);
}

// meshProxy //////////////////////////////////////////////////////////////
inline void meshProxy::Initialize(int N_, double* om_, double* delta_, double* dh_)
{
  N = N0 = N_; om = om_; delta = delta_; dh = dh_;
}

inline void meshProxy::Initialize(const mesh& m)
{
  N = N0 = m.N;
  om = m.om;
  delta = m.delta;
  dh = m.dh;
}

inline void meshProxy::resize(int n)
{
  if (n>N0) std::cerr<<"Can't resize meshProxy, to small MeshProxy!"<<std::endl;
  else N = n;
}

meshProxy& meshProxy::operator=(const meshProxy& m)
{
  resize(m.N);
  mcopy(m);
  return *this;
}

void meshProxy::SetUp(double)
{
  mesh::SetUp();
}

// mesh2D //////////////////////////////////////////////////////////////////
mesh2D::mesh2D(int N_, int Nd_) : N(N_), N0(N_), Nd(Nd_), Nd0(Nd_)
{
  //  N0=N=N_; Nd0=Nd=Nd_;
  // aloocating memory for the whole matrix
  memory = operator new (sizeof(meshProxy)*N0+3*N0*Nd0*sizeof(double)+HPoffset);

  Assert(memory!=NULL,"Out of memory");
  // calling constructors on that memory. First few bytes
  // are used to store an array of meshes (just pointers to the
  // real mesh)
  r = new (memory) meshProxy[N0];
  
  // real data comes afterwards (offset)
  int offset = sizeof(meshProxy)*N0+HPoffset;
  double *data = reinterpret_cast<double*>(static_cast<char*>(memory)+offset);
  for (int i=0; i<N0; i++){
    r[i].Initialize(Nd0, data+2*i*Nd0, data+(2*i+1)*Nd0, data+2*N0*Nd0+i*Nd0);
  }
  // ind_1 can be anywhere because it is not urgent
  // to be so fast
  ind_1 = new int*[N0];
  for (int i=0; i<N0; i++)
    ind_1[i] = new int[Nd0];
}

mesh2D::~mesh2D()
{
  for (int i=0; i<N0; i++)
    r[i].~meshProxy();
  operator delete (memory);
  memory = NULL;
  for (int i=0; i<N0; i++)
    delete[] ind_1[i];
  delete[] ind_1;
}

void mesh2D::SetUp(const mesh1D& ome, const mesh1D& Om_, const mesh1D& omc, int M, double peakpos)
{
  Om = &Om_;
  for (int i=0; i<Om_.size(); i++)
    CopyAndAddPoints(r[i],Om_[i]+peakpos,ome,omc,M,ind_1[i]);
}

void mesh2D::resize(int N_, int Nd_) 
{
  if (N_>N0 || Nd_>Nd0){   
    // freeing previously allocated memory
    for (int i=0; i<N0; i++) delete[] ind_1[i];
    delete[] ind_1;
    operator delete(memory);
    
    N0 = N = N_; Nd0 = Nd = Nd_;
    // aloocating memory for the whole matrix
    memory = operator new (sizeof(meshProxy)*N0+3*N0*Nd0*sizeof(double)+HPoffset);

    Assert(memory!=NULL,"Out of memory");
    // calling constructors on that memory. First few bytes
    // are used to store an array of meshes (just pointers to the
    // real mesh)
    r = new (memory) meshProxy[N0];
  
    // real data comes afterwards (offset)
    int offset = sizeof(meshProxy)*N0+HPoffset;
    double *data = reinterpret_cast<double*>(static_cast<char*>(memory)+offset);
    for (int i=0; i<N0; i++){
      r[i].Initialize(Nd0, data+2*i*Nd0, data+(2*i+1)*Nd0, data+2*N0*Nd0+i*Nd0);
    }
    
    // ind_1 can be anywhere because it is not necessary to be that fast
    ind_1 = new int*[N0];
    for (int i=0; i<N0; i++)
      ind_1[i] = new int[Nd0];
  } else{
    N = N_; Nd = Nd_;
  }
}

#endif
