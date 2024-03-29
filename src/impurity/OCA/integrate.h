// @Copyright 2007 Kristjan Haule
// 
#include <vector>

using namespace std;

#ifdef HEAVY_DEBUG_
ofstream dstream("integrate.debug.log");
ofstream* debugstream = &dstream;
#endif

class cmp{
  const vector<double>& f;
public:
  cmp(const vector<double>& f_) : f(f_){};
  bool operator()(int i, int j) {
    if (i<(int)f.size() && j<(int)f.size()) return f[i] < f[j];
    else {
      cerr<<" Comparing wrong objects!"<<endl;
      return false;
    }
  }
};

template <class T>
T Integrate3(const mesh1D& om0, const function1D<T>& f0, double x0,
	     const mesh1D& om1, const function1D<T>& f1, double x1,
	     const mesh1D& om2, const function1D<T>& f2, double x2)
{
  vector<double> peak(3);
  peak[0] = om0.dcenter() + x0;
  peak[1] = om1.dcenter() + x1;
  peak[2] = om2.dcenter() + x2;

  // Functions and meshes are sorted such that their peaks are in acsending order
  // Here we only calculate index array
  vector<int> index(3);
  int index_1[3];
  for (int i=0; i<3; i++) index[i]=i;
  cmp comp(peak);
  sort(index.begin(),index.end(),comp);
  for (int i=0; i<3; i++) index_1[index[i]]=i;
  // And here we make array of pointers to mimic sorted array of functions and meshes
  mesh1D const* omp[3];
  function1D<T> const* fp[3];
  omp[index_1[0]] = &om0;
  omp[index_1[1]] = &om1;
  omp[index_1[2]] = &om2;
  fp [index_1[0]] = &f0;
  fp [index_1[1]] = &f1;
  fp [index_1[2]] = &f2;
  // Sorted arrays for xi == xs and for peak == peaks
  double xs[3], peaks[3];
  xs [index_1[0]] = x0;
  xs [index_1[1]] = x1;
  xs [index_1[2]] = x2;
  for (int i=0; i<3; i++) peaks[index_1[i]] = peak[i];
  // Three separate intervals are defined such that
  // in each interval one of the meshes is most dense and
  // is being used to calculate convolution on that paricular mesh
  double start[3], end[3];
  start[0] = (*omp[0])[0]-xs[0];
  end  [0] = 0.5*(peaks[0]+peaks[1])-xs[0];
  start[1] = 0.5*(peaks[0]+peaks[1])-xs[1];
  end  [1] = 0.5*(peaks[1]+peaks[2])-xs[1];
  start[2] = 0.5*(peaks[1]+peaks[2])-xs[2];
  end  [2] = (*omp[2]).last()-xs[2];
  // The indexes for the start end end of each interval is calculated
  int is[3], ie[3];
  for (int i=0; i<3; i++) {
    is[i] = omp[i]->find_(start[i]);
    ie[i] = omp[i]->find_(end[i]);
  }
  is[0]=0;
  ie[2]=omp[2]->size()-1;
    
  T fx[3];
  T sum=0;
  tint pos1, pos2;
  intpar p1,p2;
  pos1 = omp[1]->InitInterpLeft();
  pos2 = omp[2]->InitInterpLeft();
  // Integral on the first interval
  for (int i=is[0]; i<=ie[0]; i++){
    double x = (*omp[0])[i];
    p1 = omp[1]->InterpLeft(x+xs[0]-xs[1],pos1);
    p2 = omp[2]->InterpLeft(x+xs[0]-xs[2],pos2);
    fx[0] = (*fp[0])[i];
    fx[1] = (*fp[1])(p1);
    fx[2] = (*fp[2])(p2);
    sum += fx[0]*fx[1]*fx[2]*omp[0]->Dh(i);
  }
  // Interstitial region between the first and the second interval
  if (ie[0]+1<omp[0]->size() && is[1]+1<omp[1]->size()){
    double mm0=(*omp[0])[ie[0]];
    double mm1=(*omp[0])[ie[0]+1];
    double mm2=(*omp[1])[is[1]];
    double mm3=(*omp[1])[is[1]+1];
    double dhm = 0.5*(mm3-mm1+xs[1]-xs[0]);
    double dhp = 0.5*(mm2-mm0+xs[1]-xs[0]);
    sum += fx[0]*fx[1]*fx[2]*dhm;
    double x = (*omp[1])[is[1]+1];
    pos1 = omp[0]->InitInterpLeft();
    pos2 = omp[2]->InitInterpLeft();
    p1 = omp[0]->InterpLeft(x+xs[1]-xs[0],pos1);
    p2 = omp[2]->InterpLeft(x+xs[1]-xs[2],pos2);
    fx[0] = (*fp[0])(p1);
    fx[1] = (*fp[1])[is[1]+1];
    fx[2] = (*fp[2])(p2);
    sum += fx[0]*fx[1]*fx[2]*dhp;
  }
  // Integral over the second interval
  for (int i=is[1]+1; i<=ie[1]; i++){
    double x = (*omp[1])[i];
    p1 = omp[0]->InterpLeft(x+xs[1]-xs[0],pos1);
    p2 = omp[2]->InterpLeft(x+xs[1]-xs[2],pos2);
    fx[0] = (*fp[0])(p1);
    fx[1] = (*fp[1])[i];
    fx[2] = (*fp[2])(p2);
    sum += fx[0]*fx[1]*fx[2]*omp[1]->Dh(i);
  }
  // And the second interstitial region
  if (ie[1]+1<omp[1]->size() && is[2]+1<omp[2]->size()){
    double mm0=(*omp[1])[ie[1]];
    double mm1=(*omp[1])[ie[1]+1];
    double mm2=(*omp[2])[is[2]];
    double mm3=(*omp[2])[is[2]+1];
    double dhm = 0.5*(mm3-mm1+xs[2]-xs[1]);
    double dhp = 0.5*(mm2-mm0+xs[2]-xs[1]);
    sum += fx[0]*fx[1]*fx[2]*dhm;
    double x = (*omp[2])[is[2]+1];
    pos1 = omp[0]->InitInterpLeft();
    pos2 = omp[1]->InitInterpLeft();
    p1 = omp[0]->InterpLeft(x+xs[2]-xs[0],pos1);
    p2 = omp[1]->InterpLeft(x+xs[2]-xs[1],pos2);
    fx[0] = (*fp[0])(p1);
    fx[1] = (*fp[1])(p2);
    fx[2] = (*fp[2])[is[2]+1];
    sum += fx[0]*fx[1]*fx[2]*dhp;
  }
  // Finally, integral over the third interval
  for (int i=is[2]+1; i<=ie[2]; i++){
    double x = (*omp[2])[i];
    p1 = omp[0]->InterpLeft(x+xs[2]-xs[0],pos1);
    p2 = omp[1]->InterpLeft(x+xs[2]-xs[1],pos2);
    fx[0] = (*fp[0])(p1);
    fx[1] = (*fp[1])(p2);
    fx[2] = (*fp[2])[i];
    sum += fx[0]*fx[1]*fx[2]*omp[2]->Dh(i);
  }
  return sum;
}

//typedef Hxy rT;

template <class T0, class T1, class T2>
T0 Integrate3_(const mesh1D& om0, const function1D<T0>& f0, double x0,
	       const mesh1D& om1, const function1D<T1>& f1, double x1,
	       const mesh1D& om2, const function1D<T2>& f2, double x2)
{
  vector<double> peak(3);
  peak[0] = om0.dcenter() + x0;
  peak[1] = om1.dcenter() + x1;
  peak[2] = om2.dcenter() + x2;

  // Functions and meshes are sorted such that their peaks are in acsending order
  // Here we only calculate index array
  vector<int> index(3);
  for (int i=0; i<3; i++) index[i]=i;
  cmp comp(peak);
  sort(index.begin(),index.end(),comp);
  int lind = index[0]*9+index[1]*3+index[2];
  switch (lind){
  case 5  : return IntegrateInside(om0, f0, x0, peak[0], om1, f1, x1, peak[1], om2, f2, x2, peak[2]);// 0,1,2
  case 7  : return IntegrateInside(om0, f0, x0, peak[0], om2, f2, x2, peak[2], om1, f1, x1, peak[1]);// 0,2,1
  case 11 : return IntegrateInside(om1, f1, x1, peak[1], om0, f0, x0, peak[0], om2, f2, x2, peak[2]);// 1,0,2
  case 15 : return IntegrateInside(om1, f1, x1, peak[1], om2, f2, x2, peak[2], om0, f0, x0, peak[0]);// 1,2,0
  case 19 : return IntegrateInside(om2, f2, x2, peak[2], om0, f0, x0, peak[0], om1, f1, x1, peak[1]);// 2,0,1
  case 21 : return IntegrateInside(om2, f2, x2, peak[2], om1, f1, x1, peak[1], om0, f0, x0, peak[0]);// 2,1,0
  default : cerr<<"Something very wrong!!\n"; return 0;// wrong!
  }
}

template <class T0, class T1, class T2>
T0 IntegrateInside(const mesh1D& om0, const function1D<T0>& f0, double x0, double peak0,
		   const mesh1D& om1, const function1D<T1>& f1, double x1, double peak1,
		   const mesh1D& om2, const function1D<T2>& f2, double x2, double peak2)
{
  // Three separate intervals are defined such that
  // in each interval one of the meshes is most dense and
  // is being used to calculate convolution on that paricular mesh
  int is[3], ie[3];
  is[0] = 0;
  ie[0] = om0.find_(0.5*(peak0+peak1)-x0);
  is[1] = om1.find_(0.5*(peak0+peak1)-x1);
  ie[1] = om1.find_(0.5*(peak1+peak2)-x1);
  is[2] = om2.find_(0.5*(peak1+peak2)-x2);
  ie[2] = om2.size()-1;
  
  T0 fx0;
  T1 fx1;
  T2 fx2;
  T0 sum=0;
  tint pos1, pos2;
  intpar p1,p2;
  pos1 = om1.InitInterpLeft();
  pos2 = om2.InitInterpLeft();
  // Integral on the first interval
  for (int i=is[0]; i<=ie[0]; i++){
    double x = om0[i];
    p1 = om1.InterpLeft(x+x0-x1,pos1);
    p2 = om2.InterpLeft(x+x0-x2,pos2);
    Setv(fx0, f0, i);
    Seti(fx1, f1, p1);
    Seti(fx2, f2, p2);
    AddProduct(fx0,fx1,fx2, sum, om0.Dh(i));
  }
  // Interstitial region between the first and the second interval
  if (ie[0]+1<om0.size() && is[1]+1<om1.size()){
    double mm0= om0[ie[0]];
    double mm1= om0[ie[0]+1];
    double mm2= om1[is[1]];
    double mm3= om1[is[1]+1];
    double dhm = 0.5*(mm3-mm1+x1-x0);
    double dhp = 0.5*(mm2-mm0+x1-x0);
    AddProduct(fx0,fx1,fx2,sum, dhm);
    double x = om1[is[1]+1];
    pos1 = om0.InitInterpLeft();
    pos2 = om2.InitInterpLeft();
    p1 = om0.InterpLeft(x+x1-x0,pos1);
    p2 = om2.InterpLeft(x+x1-x2,pos2);
    Seti(fx0, f0, p1);
    Setv(fx1, f1, is[1]+1);
    Seti(fx2, f2, p2);
    AddProduct(fx0,fx1,fx2, sum, dhp);
    AddProduct(fx0,fx1,fx2, sum, om1.Dh(is[1]+1));
  }
  // Integral over the second interval
  for (int i=is[1]+2; i<=ie[1]; i++){
    double x = om1[i];
    p1 = om0.InterpLeft(x+x1-x0,pos1);
    p2 = om2.InterpLeft(x+x1-x2,pos2);
    Seti(fx0, f0, p1);
    Setv(fx1, f1, i);
    Seti(fx2, f2, p2);
    AddProduct(fx0,fx1,fx2, sum, om1.Dh(i));
  }
  // And the second interstitial region
  if (ie[1]+1<om1.size() && is[2]+1<om2.size()){
    double mm0= om1[ie[1]];
    double mm1= om1[ie[1]+1];
    double mm2= om2[is[2]];
    double mm3= om2[is[2]+1];
    double dhm = 0.5*(mm3-mm1+x2-x1);
    double dhp = 0.5*(mm2-mm0+x2-x1);
    AddProduct(fx0,fx1,fx2, sum, dhm);
    double x = om2[is[2]+1];
    pos1 = om0.InitInterpLeft();
    pos2 = om1.InitInterpLeft();
    p1 = om0.InterpLeft(x+x2-x0,pos1);
    p2 = om1.InterpLeft(x+x2-x1,pos2);
    Seti(fx0, f0, p1);
    Seti(fx1, f1, p2);
    Setv(fx2, f2, is[2]+1);
    AddProduct(fx0,fx1,fx2, sum, dhp);
    AddProduct(fx0,fx1,fx2, sum, om2.Dh(is[2]+1));
  }
  // Finally, integral over the third interval
  for (int i=is[2]+2; i<=ie[2]; i++){
    double x = om2[i];
    p1 = om0.InterpLeft(x+x2-x0,pos1);
    p2 = om1.InterpLeft(x+x2-x1,pos2);
    Seti(fx0, f0, p1);
    Seti(fx1, f1, p2);
    Setv(fx2, f2, i);
    AddProduct(fx0, fx1, fx2, sum, om2.Dh(i));
  }
  return sum;
}

/*******************************************************************/
/* Calculates convolution of two functions where one function is   */
/* simple function with one peak (mesh is densed only at one point)*/
/* and combined function of two simple functions (has two peaks    */
/* with dense mesh at two different points)                        */
/* This algorith is faster than calculating convolution of three   */
/* simple functions because in this case only one interpolation is */
/* needed in the inner loop while othervise two interpolations     */
/* are necesary                                                    */
/*******************************************************************/
template <class T1, class T2>
T1 Integrate2c1(const mesh1D& omeps, const function1D<T1>& f0f1,
		const mesh1D& om2, const functionb<T2>& f2,
		double x0, double x1, double x2, double peak0, double peak1, double peak2)
{
  vector<double> peak(3);
  peak[0] = peak0;
  peak[1] = peak1;
  peak[2] = peak2;

  // Functions and meshes are sorted such that their peaks are in acsending order
  // Here we only calculate index array
  vector<int> index(3);
  for (int i=0; i<3; i++) index[i]=i;
  cmp comp(peak);
  sort(index.begin(),index.end(),comp);
  int lind = index[0]*9+index[1]*3+index[2];
  T1 dummy;// Here I didn't found yet any better solution! Please help if you can find more convenient solution!
  // A template can not differ only in return type, therefore IntegrateInside2 can not be template
  // of three types where one of them is return type only. Therefore I added another dummy argument
  // that has the same type as return type such that the compiler knows that return type might be different
  // from the other types involved. I think there must be better solution (see for example Stroustrup p.575)
  // where he returnes a template argument without using dummy argument for function get_temporary()!
  switch (lind){
    // f0 and f1 are both defined on common mesh omeps while f2 is on mesh om2
    // here we have two different cases:
    //             a) peak2 is inside peak0 and peak1
    //                In this case we need to distinguis three intervals. Function IntegrateInside3
    //                is called
    //             b) peak2 is not inside peak0 and peak1
    //                In this case we have can divide the whole region into two intervals only.
    //                IntegrateInside2 is called.
    //										            order of peaks
  case 5  : return IntegrateInside2(omeps, f0f1, x0, peak[1], om2, f2, x2, peak[2],dummy);         // 0,1,2
  case 11 : return IntegrateInside2(omeps, f0f1, x0, peak[0], om2, f2, x2, peak[2],dummy);         // 1,0,2
  case 19 : return IntegrateInside2(om2, f2, x2, peak[2], omeps, f0f1, x0, peak[0],dummy);         // 2,0,1
  case 21 : return IntegrateInside2(om2, f2, x2, peak[2], omeps, f0f1, x0, peak[1],dummy);         // 2,1,0
  case 7  : return IntegrateInside3(omeps, f0f1, x0, peak[0], peak[1], om2, f2, x2, peak[2],dummy);// 0,2,1
  case 15 : return IntegrateInside3(omeps, f0f1, x0, peak[1], peak[0], om2, f2, x2, peak[2],dummy);// 1,2,0
  default : cerr<<"Something very wrong!!\n"; return T1(0.0);// wrong!
  }
}

template <class T0, class T1, class rT>
rT IntegrateInside2(const mesh1D& om0, const functionb<T0>& f0, double x0, double peak0,
		    const mesh1D& om1, const functionb<T1>& f1, double x1, double peak1,
		    const rT&)
{
  // Two separate intervals are defined such that
  // in each interval one of the meshes is most dense and
  // is being used to calculate convolution on that paricular mesh
  int is[2], ie[2];
  is[0] = 0;
  ie[0] = om0.find_(0.5*(peak0+peak1)-x0);
  is[1] = om1.find_(0.5*(peak0+peak1)-x1);
  ie[1] = om1.size()-1;
  
  T0 fx0;
  T1 fx1;
  rT sum=0;
  tint pos0, pos1;
  intpar p0, p1;
  pos1 = om1.InitInterpLeft();
  // Integral on the first interval
  for (int i=is[0]; i<=ie[0]; i++){
    double x = om0[i];
    p1 = om1.InterpLeft(x+x0-x1,pos1);
    Setv(fx0, f0, i);
    Seti(fx1, f1, p1);
    AddProduct(fx0,fx1, sum, om0.Dh(i));
    //    cout<<x+x0<<" "<<fx0<<" "<<fx1<<" "<<fx0*fx1<<endl;
#ifdef HEAVY_DEBUG_
    *debugstream<<x+x0<<" ";
    PrintProduct(*debugstream,fx0,fx1);
#endif   
  }
  // Interstitial region between the first and the second interval
  if (ie[0]+1<om0.size() && is[1]+1<om1.size()){
    double mm0= om0[ie[0]];
    double mm1= om0[ie[0]+1];
    double mm2= om1[is[1]];
    double mm3= om1[is[1]+1];
    double dhm = 0.5*(mm3-mm1+x1-x0);
    double dhp = 0.5*(mm2-mm0+x1-x0);
    AddProduct(fx0,fx1, sum, dhm);
    double x = om1[is[1]+1];
    pos0 = om0.InitInterpLeft();
    p0 = om0.InterpLeft(x+x1-x0,pos0);
    Seti(fx0, f0, p0);
    Setv(fx1, f1, is[1]+1);
    AddProduct(fx0,fx1, sum, dhp);
    AddProduct(fx0,fx1, sum, om1.Dh(is[1]+1));
    //    cout<<x+x1<<" "<<fx0<<" "<<fx1<<" "<<fx0*fx1<<endl;
#ifdef HEAVY_DEBUG_
    *debugstream<<x+x1<<" ";
    PrintProduct(*debugstream,fx0,fx1);
#endif   
  }
  // Integral over the second interval
  for (int i=is[1]+2; i<=ie[1]; i++){
    double x = om1[i];
    p0 = om0.InterpLeft(x+x1-x0,pos0);
    Seti(fx0, f0, p0);
    Setv(fx1, f1, i);
    AddProduct(fx0,fx1, sum, om1.Dh(i));
    //    cout<<x+x1<<" "<<fx0<<" "<<fx1<<" "<<fx0*fx1<<endl;
#ifdef HEAVY_DEBUG_
    *debugstream<<x+x1<<" ";
    PrintProduct(*debugstream,fx0,fx1);
#endif   
  }
  return sum;
}
template <class T0, class T2, class rT>
rT IntegrateInside3(const mesh1D& omeps, const function1D<T0>& f0f1, double x0, double peak0, double peak1,
		    const mesh1D& om2, const functionb<T2>& f2, double x2, double peak2, const rT&)
{
  // Three separate intervals are defined such that
  // in each interval one of the meshes is most dense and
  // is being used to calculate convolution on that paricular mesh
  int is[3], ie[3];
  is[0] = 0;
  ie[0] = omeps.find_(0.5*(peak0+peak2)-x0);
  is[2] = om2.find_  (0.5*(peak0+peak2)-x2);
  ie[2] = om2.find_  (0.5*(peak1+peak2)-x2);
  is[1] = omeps.find_(0.5*(peak1+peak2)-x0);
  ie[1] = omeps.size()-1;
  
  T0 fx0;
  T2 fx2;
  rT sum=0;
  tint pos0, pos2;
  intpar p0,p2;
  pos2 = om2.InitInterpLeft();
  // Integral on the first interval over combined mesh epsom that has two peaks
  for (int i=is[0]; i<=ie[0]; i++){
    double x = omeps[i];
    p2 = om2.InterpLeft(x+x0-x2,pos2);
    Setv(fx0, f0f1, i);
    Seti(fx2, f2, p2);
    AddProduct(fx0,fx2, sum, omeps.Dh(i));
    //    cout<<x+x0<<" "<<fx0<<" "<<fx2<<" "<<fx0*fx2<<endl;
#ifdef HEAVY_DEBUG_
    *debugstream<<x+x0<<" ";
    PrintProduct(*debugstream,fx0,fx2);
#endif   
  }
  // Interstitial region between the first and the second interval
  if (ie[0]+1<omeps.size() && is[2]+1<om2.size()){
    double mm0= omeps[ie[0]];
    double mm1= omeps[ie[0]+1];
    double mm2= om2[is[2]];
    double mm3= om2[is[2]+1];
    double dhm = 0.5*(mm3-mm1+x2-x0);
    double dhp = 0.5*(mm2-mm0+x2-x0);
    AddProduct(fx0,fx2, sum, dhm);
    double x = om2[is[2]+1];
    pos0 = omeps.InitInterpLeft();
    p0 = omeps.InterpLeft(x+x2-x0,pos0);
    Seti(fx0, f0f1, p0);
    Setv(fx2, f2, is[2]+1);
    AddProduct(fx0,fx2, sum, dhp);
    AddProduct(fx0,fx2, sum, om2.Dh(is[2]+1));
    //    cout<<x+x2<<" "<<fx0<<" "<<fx2<<" "<<fx0*fx2<<endl;
#ifdef HEAVY_DEBUG_
    *debugstream<<x+x2<<" ";
    PrintProduct(*debugstream,fx0,fx2);
#endif   
  }
  // Integral over the second interval - noncombined mesh with single peak
  for (int i=is[2]+2; i<=ie[2]; i++){
    double x = om2[i];
    p0 = omeps.InterpLeft(x+x2-x0,pos0);
    Seti(fx0, f0f1, p0);
    Setv(fx2, f2, i);
    AddProduct(fx0,fx2, sum, om2.Dh(i));
    //    cout<<x+x2<<" "<<fx0<<" "<<fx2<<" "<<fx0*fx2<<endl;
#ifdef HEAVY_DEBUG_
    *debugstream<<x+x2<<" ";
    PrintProduct(*debugstream,fx0,fx2);
#endif   
  }
  // And the second interstitial region
  if (ie[2]+1<om2.size() && is[1]+1<omeps.size()){
    double mm0= om2[ie[2]];
    double mm1= om2[ie[2]+1];
    double mm2= omeps[is[1]];
    double mm3= omeps[is[1]+1];
    double dhm = 0.5*(mm3-mm1+x0-x2);
    double dhp = 0.5*(mm2-mm0+x0-x2);
    AddProduct(fx0,fx2, sum, dhm);
    double x = omeps[is[1]+1];
    pos2 = om2.InitInterpLeft();
    p2 = om2.InterpLeft(x+x0-x2,pos2);
    Setv(fx0, f0f1, is[1]+1);
    Seti(fx2, f2, p2);
    AddProduct(fx0,fx2, sum, dhp);
    AddProduct(fx0,fx2, sum, omeps.Dh(is[1]+1));
    //    cout<<x+x0<<" "<<fx0<<" "<<fx2<<" "<<fx0*fx2<<endl;
#ifdef HEAVY_DEBUG_
    *debugstream<<x+x0<<" ";
    PrintProduct(*debugstream,fx0,fx2);
#endif   
  }
  // Finally, integral over the third interval over combined mesh epsom
  for (int i=is[1]+2; i<=ie[1]; i++){
    double x = omeps[i];
    p2 = om2.InterpLeft(x+x0-x2,pos2);
    Setv(fx0, f0f1, i);
    Seti(fx2, f2, p2);
    AddProduct(fx0, fx2, sum, omeps.Dh(i));
    //    cout<<x+x0<<" "<<fx0<<" "<<fx2<<" "<<fx0*fx2<<endl;
#ifdef HEAVY_DEBUG_
    *debugstream<<x+x0<<" ";
    PrintProduct(*debugstream,fx0,fx2);
#endif   
  }
  return sum;
}

/*************************************************************************/
/* This function makes a combined function out of two simple functions   */
/* Each simple function can have one peak only (mesh is dense only at one*/
/* point while the combined mesh is densed at two points. The combined   */
/* function is obtained by routine Set and can be a product of the two   */
/* simple function or more complicated construction.                     */
/*************************************************************************/
template <class Th, class Ta, class Tg>
void MakeCombinedFunctions(mesh1D& omeps, function1D<Th>& f0f1,
			   const mesh1D& om0, const functionb<Tg>& f0, double x0,
			   const mesh1D& om1, const functionb<Ta>& f1, double x1)
{
  double peak0 = om0.dcenter() + x0;
  double peak1 = om1.dcenter() + x1;
  double meja = 0.5*(peak0+peak1);
  if (peak0<peak1)
    MakeCombinedFunctionsInside(omeps, f0f1, om0, f0, x0, om1, f1, x1, meja, 0.0);
  else
    MakeCombinedFunctionsInside(omeps, f0f1, om1, f1, x1, om0, f0, x0, meja, x1-x0);

  omeps.SetUp(om0.dcenter());
}
template <class Th, class T1, class T2>
void MakeCombinedFunctionsInside(mesh1D& omeps, function1D<Th>& f0f1,
				 const mesh1D& om0, const functionb<T1>& f0, double x0,
				 const mesh1D& om1, const functionb<T2>& f1, double x1,
				 double meja, double shift)
{
  int is[2], ie[2];
  is[0] = 0;
  ie[0] = om0.find_(meja-x0);
  is[1] = om1.find_(meja-x1);
  ie[1] = om1.size()-1;
  int totsize = ie[0]+ie[1]-is[0]-is[1]+2;
  omeps.resize(totsize);
  f0f1.resize(totsize);
  int j=0; intpar p0, p1;
  tint pos1 = om1.InitInterpLeft();
  T1 fx0;
  T2 fx1;
  for (int i=is[0]; i<=ie[0]; i++, j++){
    omeps[j] = om0[i]+shift;
    p1 = om1.InterpLeft(om0[i]+x0-x1,pos1);
    Setv(fx0, f0, i);
    Seti(fx1, f1, p1);
    Setv(f0f1[j], fx0, fx1);
    //    f0f1[j] = f0[i]*f1(p1);
    //    cout<<omeps[j]<<" "<<f0[i]<<" "<<f1(p1)<<" "<<f0f1[j]<<endl;
  }
  tint pos0 = om0.InitInterpLeft();
  p0 = om0.InterpLeft(meja-x0,pos0);
  p1 = om1.InterpLeft(meja-x1,pos1);
  omeps[j] = meja-x0+shift;
  Seti(fx0, f0, p0);
  Seti(fx1, f1, p1);
  Setv(f0f1[j], fx0, fx1);
  //  f0f1[j] = f0(p0)*f1(p1);
  j++;
  for (int i=is[1]+1; i<=ie[1]; i++, j++){
    omeps[j] = om1[i]+x1-x0+shift;
    p0 = om0.InterpLeft(om1[i]+x1-x0,pos0);
    Seti(fx0, f0, p0);
    Setv(fx1, f1, i);
    Setv(f0f1[j], fx0, fx1);
    //    f0f1[j] = f0(p0)*f1[i];
    //    cout<<omeps[j]<<" "<<f0(p0)<<" "<<f1[i]<<" "<<f0f1[j]<<endl;
  }
}

void OldMakeCombinedFunctions_(mesh1D& omeps, function1D<double>& f0f1,
			   const mesh1D& om0, const function1D<double>& f0, double x0,
			   const mesh1D& om1, const function1D<double>& f1, double x1)
{
  double peak[2];
  peak[0] = om0.dcenter() + x0;
  peak[1] = om1.dcenter() + x1;
  int is[2], ie[2];
  double meja = 0.5*(peak[0]+peak[1]);
  if (peak[0]<peak[1]){
    is[0] = 0;
    ie[0] = om0.find_(meja-x0);
    is[1] = om1.find_(meja-x1);
    ie[1] = om1.size()-1;
    int totsize = ie[0]+ie[1]-is[0]-is[1]+2;
    omeps.resize(totsize);
    f0f1.resize(totsize);
    int j=0; intpar p0, p1;
    tint pos1 = om1.InitInterpLeft();
    for (int i=is[0]; i<=ie[0]; i++, j++){
      omeps[j] = om0[i];
      p1 = om1.InterpLeft(om0[i]+x0-x1,pos1);
      f0f1[j] = f0[i]*f1(p1);
    }
    tint pos0 = om0.InitInterpLeft();
    p0 = om0.InterpLeft(meja-x0,pos0);
    p1 = om1.InterpLeft(meja-x1,pos1);
    omeps[j] = meja-x0;
    f0f1[j] = f0(p0)*f1(p1);
    j++;
    for (int i=is[1]+1; i<=ie[1]; i++, j++){
      omeps[j] = om1[i]+x1-x0;
      p0 = om0.InterpLeft(om1[i]+x1-x0,pos0);
      f0f1[j] = f0(p0)*f1[i];
    }
  }else{
    is[1] = 0;
    ie[1] = om1.find_(meja-x1);
    is[0] = om0.find_(meja-x0);
    ie[0] = om0.size()-1;
    int totsize = ie[0]+ie[1]-is[0]-is[1]+2;
    omeps.resize(totsize);
    f0f1.resize(totsize);
    int j=0; intpar p0, p1;
    tint pos0 = om0.InitInterpLeft();
    for (int i=is[1]; i<=ie[1]; i++, j++){
      omeps[j] = om1[i]+x1-x0;
      p0 = om0.InterpLeft(om1[i]+x1-x0,pos0);
      f0f1[j] = f0(p0)*f1[i];
    }
    p0 = om0.InterpLeft(meja-x0,pos0);
    tint pos1 = om1.InitInterpLeft();
    p1 = om1.InterpLeft(meja-x1,pos1);
    omeps[j] = meja-x0;
    f0f1[j] = f0(p0)*f1(p1);
    j++;
    for (int i=is[0]+1; i<=ie[0]; i++, j++){
      omeps[j] = om0[i];
      p1 = om1.InterpLeft(om0[i]+x0-x1,pos1);
      f0f1[j] = f0[i]*f1(p1);
    }
  }
  omeps.SetUp(om0.dcenter());
}
