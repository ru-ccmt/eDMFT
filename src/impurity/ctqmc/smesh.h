// @Copyright 2007 Kristjan Haule
// 
#ifndef MESH_
#define MESH_
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include "sutil.h"
#include "zeroin.h"

/////////////////////// Simplified verison of mesh1D is defined in this header  ////////////////////////

//*****************************************//
// Classes implemented in this header file //
//*****************************************//
class mesh;  // Basic class
class mesh1D;// 1D grid 
class MeshJoin; // Small class for root-finding routine
class mesh2D;

typedef int tint; // For position in mesh we use tint clss, which is at the moment just int. 
                  // We use class name for flexibility in implementation.
//*********************************************************//
// Additional nonstandard classes used in this header file //
//*********************************************************//
class intpar;// Simple class for exchanging data between functions and meshes when doing linear interpolation

//************************************************//
//  Declaration of basic class for holding grids  //
//************************************************//
class mesh{
protected:
  int N, N0;       // size, size of the allocated memory (might be larger than size)
  double *om;      // grid points
  double *delta;   // contains weights for interpolation
  double *dh;      // contains weights for integration
  static const int dN = 10; // when searching ordered table, searching is done first between a0 and a0+dN...
protected:
  mesh(): N(0),N0(0),om(NULL),delta(NULL),dh(NULL){}; // constructor is made protected such that mesh can not be instantiated
  ~mesh(){};                                          // This class is used only as base class
public:
  // OPERATORS
  double& operator[](int i) {Assert(i<N,"Out of range in mesh[]"); return om[i];}
  const double& operator[](int i) const {Assert(i<N,"Out of range in mesh[]"); return om[i];}
  // ROUTINES FOR SEARCHING ORDERED TABLE
  int find(double x) const;           // searching table if previous position is not known
  int find_(double x, int& ia) const; // searching table forward from previous position
  int find_(double x) const; // searching table forward from previous position
  int _find(double x, int& ia) const; // searching table backward from previous position
  int findBoth(double x, int& ia) const;  // searching table in both direction - point is usually close
  // LINEAR INTERPOLATION ROUTINES - INITIALIZATION
  tint InitInterpLeft() const {return 0;}    // Initialization of position for forward search
  tint InitInterpRight() const {return N-2;} // Initialization of position for backward search 
  // LINEAR INTERPOLATION ROUTINES 
  intpar Interp(const double x) const;               // Finds parameters for linear interpolation at point x
  intpar InterpLeft(const double x, tint& a) const;  // Finds parameters for linear interpolation when freqeuncy is increasing
  intpar InterpRight(const double x, tint& a) const; // Finds parameters for linear interpolation when freqeuncy is decreasing
  intpar InterpBoth(const double x, tint& a) const;  // If frequency is believed to be close to previous frequnecy
  // OTHER SHORT MEMBER FUNCTIONS
  int size() const {return N;}
  double last() const {return om[N-1];}
  double* begin() const {return om;}
  double* end() const {return om+N;}
  double& Delta(int i) const {Assert(i<N,"Out of range in mesh.delta[]"); return delta[i];}
  double& Dh(int i) const {Assert(i<N,"Out of range in mesh.dh[]"); return dh[i];}
  const double* MemPt() const {return om;}
  double* MemPt() {return om;}
protected:
  void mcopy(const mesh& m);
  void SetUp();
private:
  int bisection(double x, int& klo, int& khi, int bi) const;
  int find(double x, int& dummy) const;  // searching table if previous position is not known
  friend std::ostream& operator<<(std::ostream& stream, const mesh1D& m);
  friend class mesh2D;
  friend class meshProxy;
};

//************************************************//
// One dimentional grid derived from class mesh.  //
// In addition to the base class, it has a	  //
// constructor that allocates necessary memory	  //
// and also destructor to prevent resorce leakes. //
// It also has an integer (center) that holds	  //
// position of the peak in the mesh		  //
//************************************************//
class mesh1D : public mesh{
  int center;
public:
  // COSTRUCTORS AND DESTRUCTORS
  mesh1D() {};              // default constructor
  mesh1D(const mesh1D& m);  // copy constructor
  explicit mesh1D(int N_);  // another constructor
  ~mesh1D();                // destructor
  // OPERATORS
  mesh1D& operator=(const mesh1D& m);                   // needs to be changed from default
  // INITIALIZATION ROUTINES
  void MakeEquidistantMesh(int N_, double start, double end);
  void MakePositiveLogTanMesh(int N, double x0, double x1, double x2, double alpha=0);
  void MakePositiveLogTanMesh0(int N, double x0, double x1, double x2, double alpha=0);
  void MakeLogTanMesh(int N_, double x0, double x1, double x2, double alpha=0);
  void resize(int n);
  void SetUp(int center_);    // knowing om, initializes delta, dh
  void SetUp(double center_); // knowing om, initializes delta, dh
  // OTHER SHORT MEMBER FUNCTIONS
  double dcenter() const { return 0.5*(om[center]+om[center+1]);}     
  int icenter() const {return center;}
  
  friend void CopyAndAddPoints(mesh1D& m, double centerOm, const mesh1D& ome, const mesh1D& omc, int M, int iremove);
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


//////////////////// Implementation of mesh ///////////////////////////////////////
/////////////////// INITIALIZATION ROUTINES ///////////////////////////////////////
inline void mesh::mcopy(const mesh& m) 
{ // To copy arrays of the class
  std::copy(m.om,   m.om+N,   om);
  std::copy(m.delta,m.delta+N,delta);
  std::copy(m.dh,   m.dh+N,   dh);
}

inline void mesh::SetUp()
{ // Initialization of delta and dh arrays
  delta[0] = 1/(om[1]-om[0]);
  dh[0] = 0.5*(om[1]-om[0]);
  for (int i=1; i<N-1; i++){
    delta[i] = 1/(om[i+1]-om[i]);
    dh[i] = 0.5*(om[i+1]-om[i-1]);
  }
 delta[N-1] = 0.0;
 dh[N-1] = 0.5*(om[N-1]-om[N-2]);
}

///////////////////// SEARCHING ROUTINES ///////////////////////////////////////////
inline int mesh::bisection(double x, int& klo, int& khi, int bi) const
{// Basic routine for searching ordered table is bisection
  int k;
  khi=bi-1;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (om[k] > x) khi=k;
    else klo=k;
  }
  return klo;
}

inline int mesh::find(double x) const
{ // If nothing is known about the point x
  if (x<=om[0]) return 0;
  if (x>=om[N-2]) return N-2;
  int ai0=0, ai1=N-1;
  return bisection(x,ai0,ai1,N);
}
inline int mesh::find(double x, int& dummy) const
{ return find(x);}

inline int mesh::find_(double x, int& ai0) const
{ // This is most often called searching routine
  // It is used for searching table in increasing order
  if (x<om[ai0+1]) return ai0; // Checks weather x is stil in [ai0:ai0+1]
  int ai1 = ai0+1;             // Makes a step
  if (ai0>=N-2) return ai0;    // Needs to check for the end of the table
  if (x<om[ai1+1]){            // Checks weather x is in [ai0+1:ai0+2]
    ai0 = ai1;
    ai1 = ai1+1;
    return ai0;
  }
  if (ai1>=N-2) return ai1; // Again checks for the end of the table
  ai0 = ai1+1;              // makes another step
  if (ai1+dN<N){            // First uses bisection is small interval between [ai1:ai1+dN]
    if (x<om[ai1+dN]) return bisection (x, ai0, ai1, ai1+dN+1);
  } else return bisection (x, ai0, ai1, N);
  if (ai1+dN<N-1) ai0 = ai1+dN;
  else ai0 = ai1+dN-1;
  return bisection (x, ai0, ai1, N); // If still not found, use bisection on the rest of the grid
}

inline int mesh::find_(double x) const
{
  int ai0=0;
  return find_(x,ai0);
}

inline int mesh::_find(double x, int& bi) const
{ // This routine is used to search ordered table in decreasing order
  int ai0=0;
  if (x>om[bi]) return bi;     // Checks weather x is still in [bi:bi+1]
  if (bi<=0) return bi;        // Checks start of the table
  if (x>om[bi-1]) return bi-1; // Checks weather x is in [bi-1:bi]
  if (bi-2<=ai0) return ai0;   // If [bi-2:bi-1] is first interval (equal to [0:1]) we are done
  if (x>om[bi-2]) return bi-2; // Checks interbal [bi-2:bi-1]
  if (bi-dN>ai0){              // Bisection only between [bi-dN:ai0]
    if (x>om[bi-dN]){
      int ai1, ai00;
      ai00 = bi-dN;
      return bisection(x,ai00,ai1,bi-1); 
    }
    bi=bi-dN+1;
  } else bi-=1; 
  {                          // If everything else failed, search everywhere between [ai0:bi]
    int ai1;
    return bisection (x,ai0,ai1,bi);
  }
}

inline int mesh::findBoth(const double x, int& ix) const
{ // This routine is uded when point is usually close to the previous point but not
  // on left or right side of it
  if (x==om[ix]) return ix;
  if (x>om[ix]) return find_(x,ix);
  else return _find(x,ix);
}

inline intpar mesh::Interp(const double x) const
{
  int i = find(x);
#ifdef _DEBUG
  double p = (x-om[i])*delta[i];
  if ((p>1.0 || p<0) && i!=0 && i!=(N-2)) {std::cerr<<"Variable p="<<p<<" at i="<<i<<std::endl;}
  return intpar(i,p);
#else
  return intpar(i,(x-om[i])*delta[i]); // return value optimization
#endif
}

inline intpar mesh::InterpLeft(const double x, tint& a) const
{
  int i = find_(x,a);
#ifdef _DEBUG
  double p = (x-om[i])*delta[i];
  if ((p>1.0 || p<0) && i!=0 && i!=(N-2)) {std::cerr<<"Variable p="<<p<<" at i="<<i<<std::endl;}
  return intpar(i,p);
#else
  return intpar(i,(x-om[i])*delta[i]); // return value optimization
#endif
}

inline intpar mesh::InterpRight(const double x, tint& a) const
{
  int i = _find(x,a);
#ifdef _DEBUG
  double p = (x-om[i])*delta[i];
  if ((p>1.0 || p<0) && i!=0 && i!=(N-2)) {std::cerr<<"Variable p="<<p<<" at i="<<i<<std::endl;}
  return intpar(i,p);
#else
  return intpar(i,(x-om[i])*delta[i]); // return value optimization
#endif
}

inline intpar mesh::InterpBoth(const double x, tint& a) const
{
  int i = findBoth(x,a);
#ifdef _DEBUG
  double p = (x-om[i])*delta[i];
  if ((p>1.0 || p<0) && i!=0 && i!=(N-2)) {std::cerr<<"Variable p="<<p<<" at i="<<i<<std::endl;}
  return intpar(i,p);
#else
  return intpar(i,(x-om[i])*delta[i]); // return value optimization
#endif
}


inline std::ostream& operator<<(std::ostream& stream, const mesh& m)
{
  int width = stream.width();
  for (int i=0; i<m.size(); i++)
    stream << std::setw(3) << i << std::setw(width) << m[i] << std::endl;
  return stream;
}

////////////////// Implementation of mesh1D ///////////////////////////////////////
/////////////////// INITIALIZATION ROUTINES ///////////////////////////////////////

inline mesh1D::mesh1D(int N_) 
{ // Prefered constructor
  N0=N=N_;
  om = new double[N0];
  delta = new double[N0];
  dh = new double[N0];
}

inline mesh1D::~mesh1D() 
{ // Destructor cleans up
  if (om != NULL){
    delete[] om;
    delete[] delta;
    delete[] dh;
    om = NULL;
    delta = NULL;
    dh = NULL;
  }
  N0=N=0;
}

inline void mesh1D::resize(int n) 
{ // Resizing necessary if default constructor was called or if mesh should be resize
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

inline mesh1D::mesh1D(const mesh1D& m) : mesh()
{ // Copy constructor
  resize(m.N);
  mcopy(m);
  center = m.center;
}

inline mesh1D& mesh1D::operator=(const mesh1D& m) 
{ // operator=
  resize(m.N);
  mcopy(m);
  center = m.center;
  return *this;
}

inline void mesh1D::SetUp(int center_)
{ mesh::SetUp();  center = center_;}
inline void mesh1D::SetUp(double dcenter)
{ mesh::SetUp();  center = find(dcenter);}

inline void mesh1D::MakeEquidistantMesh(int N_, double start, double end)
{ // For building equidistant mesh
  resize(N_);
  double x=start, dh=(end-start)/(N-1.);
  for (int i=0; i<N; i++,x+=dh) om[i] = x;
  SetUp(N/2);
}

class MeshJoin{
  double dwt, xt;
public:
  MeshJoin(double dwt_, double xt_) : dwt(dwt_), xt(xt_) {};
  double operator()(double u){
    double tg=tan(u);
    return u-atan(tg/xt)-dwt*xt*tg/(xt*xt+tg*tg);
  }
};

inline void mesh1D::MakePositiveLogTanMesh(int N_, double x0, double x1, double x2, double alpha)
{ // For building mesh which is logarithmic at small frequency and tan at large frequency
  // Only positive frequnecy used
  resize(N_);
  double eta = log(x1/x0)/(x2/x1-1);
  int N1_min = static_cast<int>((1+eta*N)/(1+eta)+0.5);
  int N1 = static_cast<int>((1+alpha)*N1_min);
  if (N1>N-2) N1=N-2;
  
  int N2 = N-N1;
  double xt = x2/x1;
  double dwt = N2*(log(x1)-log(x0))/(N1-1);

  MeshJoin mj(dwt,xt);
  double ut = zeroin(1e-5,M_PI/2-1e-5,mj,1e-10);
  
  double a = atan(tan(ut)/xt);
  double b = dwt*sin(a)*cos(a);
  double w = x1/tan(a);

  resize(N);
  for (int i=0; i<N1; i++) om[i] = exp(log(x0)+i*(log(x1)-log(x0))/(N1-1));
  for (int i=0; i<N2; i++) om[N1 + i] = w*tan(a+(i+1)*b/N2);
  SetUp(0.0);
}

inline void mesh1D::MakePositiveLogTanMesh0(int N_, double x0, double x1, double x2, double alpha)
{
  resize(N_);
  MakePositiveLogTanMesh(N_-1,x0,x1,x2,alpha);
  resize(N_);
  for (int i=N_-1; i>0; i--) om[i]=om[i-1];
  om[0]=0;
}

inline void mesh1D::MakeLogTanMesh(int N_, double x0, double x1, double x2, double alpha)
{ // For building mesh which is logarithmic at small frequency and tan at large frequency
  // Mesh is symmetric around zero frequency
  int N2 = N_/2;
  MakePositiveLogTanMesh(N2,x0,x1,x2,alpha);
  double* tom = new double[N2];
  for (int i=0; i<N2; i++) tom[i] = om[i];
  resize(2*N2);
  for (int i=0; i<N2; i++) om[i] = -tom[N2-i-1];
  for (int i=0; i<N2; i++) om[N2+i] = tom[i];
  delete[] tom;
  SetUp(0.0);
}

// compare /////////////////////////////////////////////////////////
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
  if (memory != NULL){
    for (int i=0; i<N0; i++)
      r[i].~meshProxy();
    operator delete (memory);
    memory = NULL;
    for (int i=0; i<N0; i++)
      delete[] ind_1[i];
    delete[] ind_1;
  }
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
