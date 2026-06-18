#include <iostream>
#include "basics.h"

void GiveR(Array<double,2>& Rx, const Array<double,2>& Kx, const Array<double,1>& gx, double dh)
{//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This is used to calculate definite integral of the form
  // Int_0^beta dt_n Int_0^t_n dt_{n-1} Int_0^{t_{n-1}} dt_{n-2} .... Int_0^{t_2} g(t1)*g(t2-t1)*g(t3-t2)*g(t4-t3)....*g(t_n-t_{n-1})
  // where g functions are picewise constant functions.
  //
  // But this code only does one integral, which can then be used iteratively.
  // Here we  Integrate:
  //                     Rx(t2) = Intg[ K(t1)*g(t2-t1)dt1 , {dt1,0,t2}]
  // Here gx(x) is picewise constant function with Nbin points, which are not so many that continuous integration is acceptable.
  // Further Rx(t2) and K(t1) are stored as Rx[j,:]=Polynomial(n+1) and K[l,:]=Polynomial(n).
  // The first integral of constant function is a linear function, and the second is quadratic, etc., hence at step n, we have polynomial of order n in each interval.
  // 
  // More precisely, K(t1) = \sum_{m=0}^n K_m[l=int(t1)] * (t1-l)^m , and if we select t1=l+t, where l is int and t=[0,1], then K(l+t)=\sum_{m=0}^n K_m[l] * t^m
  //  hence
  //  R(j+u) = Int_0^{j+u} dt_1 (\sum_{m=0}^n K_m[int(t1)]*(t1-l)^m ) * g(j+u-t1) =
  //           =\sum_{l=0}^{j-1} Int_l^{l+1} dt_1 (\sum_{m=0}^n K_m[l]*(t1-l)^m ) * g(j+u-t_1)  + Int_j^{j+u} dt_1 (\sum_{m=0}^n K_m[j]*(t1-j)^m ) * g(j+u-t_1) =
  //           =\sum_{l=0}^{j-1} Int_0^1 dt (\sum_{m=0}^n K_m[l]*t^m ) * g(j+u-l-t)  + Int_0^{u} dt (\sum_{m=0}^n K_m[j]*t^m ) * g(u-t) =
  //           =\sum_{l=0}^{j-1} \sum_{m=0}^n K_m[l]/(m+1) * [ u^{m+1}*g[j-l] + (1-u^{m+1})*g[j-l-1] ]  +
  //                             \sum_{m=0}^n K_m[j]/(m+1) * u^{m+1} * g[0]
  //           =\sum_{l=0}^{j-1} \sum_{m=0}^n K_m[l]/(m+1) * [ g[j-l-1] + u^{m+1}*(g[j-l]-g[j-l-1]) ]  +
  //                             \sum_{m=0}^n K_m[j]/(m+1) * u^{m+1} * g[0]
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int Nbin = gx.extent(0);
  int n = Kx.extent(1);
  Array<double,1> Kxl(n);
  Rx.resize(Nbin,n+1);
  Rx=0;
  for (int j=0; j<Nbin; ++j){
    int lstart = (j==0)? 1 : 0;
    for (int l=lstart; l<j; ++l){
      for (int i=0; i<n; ++i) Kxl(i) = Kx(l,i)/(i+1);
      Rx(j,0) += gx(j-l-1)*sum(Kxl);
      Rx(j,Range(1,n)) += (gx(j-l)-gx(j-l-1))*Kxl(Range::all());
    }
    for (int i=0; i<n; ++i) Kxl(i) = Kx(j,i)/(i+1);
    Rx(j,Range(1,n)) += gx(0)*Kxl(Range::all());
    Rx(j,Range::all()) *= dh;
  }
}

/*
double g(double x){
  double alpha = 0.0000001;
  return exp(-alpha*x);
}
int main()
{
  int Nbin=10;
  double beta=100;
  
  double dh=beta/Nbin;

  Array<double,1> gx(Nbin);
  for (int i=0; i<Nbin; ++i){
    double xms=(i+0.5)*dh;
    gx(i) = g(xms);
  }
  Array<double,2> Kx(Nbin,1);
  for (int i=0; i<Nbin; ++i) Kx(i,0) = gx(i);

  Array<double,2> Rx;
  GiveR(Rx, Kx, gx, dh); // Rx = g*g  -> order 1

  int norder=7;
  for (int iord=1; iord<norder; ++iord){
    GiveR(Kx, Rx, gx, dh); // Kx = g*g*g
    GiveR(Rx, Kx, gx, dh); // Rx = g*g*g*g -> order 2

    double at_beta=0;
    for (int l=0; l<Rx.extent(1); ++l) // polynomial \sum_{m=0}^n Rx(Nbin-1,m)*t^m
      at_beta += Rx(Nbin-1,l);
    cout<<"at_beta["<<iord+1<<"]="<<at_beta<<endl;
    
  }

  
  for (int i=0; i<Nbin; ++i){
    double xm = (0.5+i)*dh;
    cout<<xm<<" "<< setprecision(7)<<Rx(i,0)<<" "<<Rx(i,1)<<" "<<Rx(i,2)<<" "<<Rx(i,3)<<endl;
  }
  
  return 0;
}
*/
