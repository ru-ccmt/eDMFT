#include <cstdint>
#include <iostream>
#include <deque>
#include <complex>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <blitz/array.h>

using namespace std;
namespace py = pybind11;
namespace bl = blitz;
void FromSlaterToMatrixU(py::array_t<complex<double>>& _UC_, py::array_t<double>& _gck_,
			 int l, py::array_t<complex<double>>& _T2C_)
{
  //Uc(l+1,nw,nw,nw,nw,nw)
  py::buffer_info info_U = _UC_.request();
  bl::Array<complex<double>,5> UC((complex<double>*)info_U.ptr,
				  bl::shape(info_U.shape[0],info_U.shape[1],info_U.shape[2],info_U.shape[3],info_U.shape[4]),
				  bl::neverDeleteData);
  //gckgck(0:3,0:6,0:6,0:3)
  py::buffer_info info_g = _gck_.request();
  bl::Array<double,4> gck((double*)info_g.ptr, bl::shape(info_g.shape[0],info_g.shape[1],info_g.shape[2],info_g.shape[3]),
			  bl::neverDeleteData);
  // T2C(_mw_,_mw_)
  py::buffer_info info_t = _T2C_.request();
  bl::Array<complex<double>,2> T2C((complex<double>*)info_t.ptr,
				   bl::shape(info_t.shape[0],info_t.shape[1]),
				   bl::neverDeleteData);
  
  //auto gck = _gck_.mutable_unchecked<4>();  
  if (gck.extent(0)!=4){
    std::cout<<"Dimension of gaunt coeff gck.shape()[0]="<<_gck_.request().shape[0]<<" which is wrong!"<<endl;
    return;
  }
  if (gck.extent(1)!=7){
    std::cout<<"Dimension of gaunt coeff gck.shape()[1]="<<_gck_.request().shape[1]<<" which is wrong!"<<endl;
    return;
  }
  if (gck.extent(2)!=7){
    std::cout<<"Dimension of gaunt coeff gck.shape()[2]="<<_gck_.request().shape[2]<<" which is wrong!"<<endl;
    return;
  }
  if (gck.extent(3)!=4){
    std::cout<<"Dimension of gaunt coeff gck.shape()[3]="<<_gck_.request().shape[3]<<" which is wrong!"<<endl;
    return;
  }
  int mw = 2*l+1;
  int nw = 0, ns=0;
  if (T2C.extent(0) == mw){
    nw = mw;
    ns = 1;
  }else if (T2C.extent(0) == 2*mw){
    nw = 2*(2*l+1);
    ns = 2;
  }else{
    std::cout<<"ERROR in atom_d.py: T2C has wrong shape"<<endl;
    return;
  }
  //auto UC = _UC_.mutable_unchecked<5>();
  if (UC.extent(0)!=l+1 or UC.extent(1)!=nw or UC.extent(2)!=nw or UC.extent(3)!=nw or UC.extent(4)!=nw){
    std::cout<<"ERROR Uc(l+1,nw,nw,nw,nw) has wrong shape "<<_UC_.request().shape[0]<<","<<_UC_.request().shape[1]<<","<<_UC_.request().shape[2];
    std::cout<<","<<_UC_.request().shape[3]<<","<<_UC_.request().shape[4]<<endl;
  }
  
  int shft = 3-l;
  int shft2 = 2*l+1;
  
  bl::Array<complex<double>,3> Sum1(nw,nw,shft2*2), Sum2(nw,nw,shft2*2);
  bl::Array<complex<double>,2> T2Cp(T2C.extent(1),T2C.extent(0));
  T2Cp = bl::conj(T2C.transpose(bl::secondDim,bl::firstDim));

  //cout<<"nw="<<nw<<" mw="<<mw<<" gck="<<endl;
  //for (int i1=0; i1<4; ++i1)
  //  for (int i2=0; i2<7; ++i2)
  //    for (int i3=0; i3<7; ++i3)
  // 	for (int i4=0; i4<4; ++i4)
  // 	  cout<<i1<<" "<<i2<<" "<<i3<<" "<<i4<<" "<<gck(i1,i2,i3,i4)<<endl;
  
  for (int k=0; k<l+1; k++){
    Sum1=0.0;
    for (int i4=0; i4<nw; i4++){
      for (int i1=0; i1<nw; i1++){
	for (int m4=0; m4<mw; m4++){
	  for (int m1=0; m1<mw; m1++){
	    for (int s=0; s<ns; s++) Sum1(i4,i1,m1-m4+shft2) += T2Cp(i4,m4+s*mw)*gck(l,shft+m4,shft+m1,k)*T2C(m1+s*mw,i1);
	  }
	}
      }
    }
    Sum2=0.0;
    for (int i3=0; i3<nw; i3++){
      for (int i2=0; i2<nw; i2++){
	for (int m3=0; m3<mw; m3++){
	  for (int m2=0; m2<mw; m2++){
	    for (int s=0; s<ns; s++) Sum2(i3,i2,m3-m2+shft2) += T2Cp(i3,m3+s*mw)*gck(l,shft+m2,shft+m3,k)*T2C(m2+s*mw,i2);
	  }
	}
      }
    }
    for (int i4=0; i4<nw; i4++){
      for (int i3=0; i3<nw; i3++){
	for (int i2=0; i2<nw; i2++){
	  for (int i1=0; i1<nw; i1++){
	    complex<double> csum=0.0;
	    for (int dm=0; dm<shft2*2; dm++) csum += Sum1(i4,i1,dm)*Sum2(i3,i2,dm);
	    UC(k,i4,i3,i2,i1) = csum;
	  }
	}
      }
    }
  }
}
void FromSlaterToMatrixUd(py::array_t<complex<double>>& _UC_, py::array_t<double>& _gck_,
			 double l, py::array_t<complex<double>>& _T2C_)
{
  //Uc(l+1,nw,nw,nw,nw,nw)
  py::buffer_info info_U = _UC_.request();
  bl::Array<complex<double>,5> UC((complex<double>*)info_U.ptr,
				  bl::shape(info_U.shape[0],info_U.shape[1],info_U.shape[2],info_U.shape[3],info_U.shape[4]),
				  bl::neverDeleteData);
  //gckgck(0:3,0:6,0:6,0:3)
  py::buffer_info info_g = _gck_.request();
  bl::Array<double,4> gck((double*)info_g.ptr, bl::shape(info_g.shape[0],info_g.shape[1],info_g.shape[2],info_g.shape[3]),
			  bl::neverDeleteData);
  // T2C(_mw_,_mw_)
  py::buffer_info info_t = _T2C_.request();
  bl::Array<complex<double>,2> T2C((complex<double>*)info_t.ptr,
				   bl::shape(info_t.shape[0],info_t.shape[1]),
				   bl::neverDeleteData);
  
  //auto gck = _gck_.mutable_unchecked<4>();  
  if (gck.extent(0)!=4){
    std::cout<<"Dimension of gaunt coeff gck.shape()[0]="<<_gck_.request().shape[0]<<" which is wrong!"<<endl;
    return;
  }
  if (gck.extent(1)!=7){
    std::cout<<"Dimension of gaunt coeff gck.shape()[1]="<<_gck_.request().shape[1]<<" which is wrong!"<<endl;
    return;
  }
  if (gck.extent(2)!=7){
    std::cout<<"Dimension of gaunt coeff gck.shape()[2]="<<_gck_.request().shape[2]<<" which is wrong!"<<endl;
    return;
  }
  if (gck.extent(3)!=4){
    std::cout<<"Dimension of gaunt coeff gck.shape()[3]="<<_gck_.request().shape[3]<<" which is wrong!"<<endl;
    return;
  }
  int mw = 2*l+1;
  int nw = 0, ns=0;
  if (T2C.extent(0) == mw){
    nw = mw;
    ns = 1;
  }else if (T2C.extent(0) == 2*mw){
    nw = 2*(2*l+1);
    ns = 2;
  }else{
    std::cout<<"ERROR in atom_d.py: T2C has wrong shape"<<endl;
    return;
  }
  //auto UC = _UC_.mutable_unchecked<5>();
  if (UC.extent(0)!=int(l)+1 or UC.extent(1)!=nw or UC.extent(2)!=nw or UC.extent(3)!=nw or UC.extent(4)!=nw){
    std::cout<<"ERROR Uc(l+1,nw,nw,nw,nw) has wrong shape "<<_UC_.request().shape[0]<<","<<_UC_.request().shape[1]<<","<<_UC_.request().shape[2];
    std::cout<<","<<_UC_.request().shape[3]<<","<<_UC_.request().shape[4]<<endl;
  }
  
  int shft = 3-l;
  int shft2 = 2*l+1;
  int _l_ = int(l);
  if (abs(l - round(l))>0.3){// l is half-integer
    shft = 2-_l_;
  }
	  
  bl::Array<complex<double>,3> Sum1(nw,nw,shft2*2), Sum2(nw,nw,shft2*2);
  bl::Array<complex<double>,2> T2Cp(T2C.extent(1),T2C.extent(0));
  T2Cp = bl::conj(T2C.transpose(bl::secondDim,bl::firstDim));

  //cout<<"nw="<<nw<<" mw="<<mw<<" gck="<<endl;
  //for (int i1=0; i1<4; ++i1)
  //  for (int i2=0; i2<7; ++i2)
  //    for (int i3=0; i3<7; ++i3)
  // 	for (int i4=0; i4<4; ++i4)
  // 	  cout<<i1<<" "<<i2<<" "<<i3<<" "<<i4<<" "<<gck(i1,i2,i3,i4)<<endl;
  
  for (int k=0; k<_l_+1; k++){
    Sum1=0.0;
    for (int i4=0; i4<nw; i4++){
      for (int i1=0; i1<nw; i1++){
	for (int m4=0; m4<mw; m4++){
	  for (int m1=0; m1<mw; m1++){
	    for (int s=0; s<ns; s++) Sum1(i4,i1,m1-m4+shft2) += T2Cp(i4,m4+s*mw)*gck(_l_,shft+m4,shft+m1,k)*T2C(m1+s*mw,i1);
	  }
	}
      }
    }
    Sum2=0.0;
    for (int i3=0; i3<nw; i3++){
      for (int i2=0; i2<nw; i2++){
	for (int m3=0; m3<mw; m3++){
	  for (int m2=0; m2<mw; m2++){
	    for (int s=0; s<ns; s++) Sum2(i3,i2,m3-m2+shft2) += T2Cp(i3,m3+s*mw)*gck(_l_,shft+m2,shft+m3,k)*T2C(m2+s*mw,i2);
	  }
	}
      }
    }
    for (int i4=0; i4<nw; i4++){
      for (int i3=0; i3<nw; i3++){
	for (int i2=0; i2<nw; i2++){
	  for (int i1=0; i1<nw; i1++){
	    complex<double> csum=0.0;
	    for (int dm=0; dm<shft2*2; dm++) csum += Sum1(i4,i1,dm)*Sum2(i3,i2,dm);
	    UC(k,i4,i3,i2,i1) = csum;
	  }
	}
      }
    }
  }
}

void FromSlaterToMatrixU_diagonal(py::array_t<complex<double>>& _UC_, py::array_t<double>& _gck_,
				  int l, py::array_t<complex<double>>& _T2C_)
{
  //Uc(l+1,nw,nw,nw,nw,nw)
  py::buffer_info info_U = _UC_.request();
  bl::Array<complex<double>,3> UC((complex<double>*)info_U.ptr,
				  bl::shape(info_U.shape[0],info_U.shape[1],info_U.shape[2]),
				  bl::neverDeleteData);
  //gckgck(0:3,0:6,0:6,0:3)
  py::buffer_info info_g = _gck_.request();
  bl::Array<double,4> gck((double*)info_g.ptr, bl::shape(info_g.shape[0],info_g.shape[1],info_g.shape[2],info_g.shape[3]),
			  bl::neverDeleteData);
  // T2C(_mw_,_mw_)
  py::buffer_info info_t = _T2C_.request();
  bl::Array<complex<double>,2> T2C((complex<double>*)info_t.ptr,
				   bl::shape(info_t.shape[0],info_t.shape[1]),
				   bl::neverDeleteData);
  
  //auto gck = _gck_.mutable_unchecked<4>();  
  if (gck.extent(0)!=4){
    std::cout<<"Dimension of gaunt coeff gck.shape()[0]="<<_gck_.request().shape[0]<<" which is wrong!"<<endl;
    return;
  }
  if (gck.extent(1)!=7){
    std::cout<<"Dimension of gaunt coeff gck.shape()[1]="<<_gck_.request().shape[1]<<" which is wrong!"<<endl;
    return;
  }
  if (gck.extent(2)!=7){
    std::cout<<"Dimension of gaunt coeff gck.shape()[2]="<<_gck_.request().shape[2]<<" which is wrong!"<<endl;
    return;
  }
  if (gck.extent(3)!=4){
    std::cout<<"Dimension of gaunt coeff gck.shape()[3]="<<_gck_.request().shape[3]<<" which is wrong!"<<endl;
    return;
  }
  int mw = 2*l+1;
  int nw = 0, ns=0;
  if (T2C.extent(0) == mw){
    nw = mw;
    ns = 1;
  }else if (T2C.extent(0) == 2*mw){
    nw = 2*(2*l+1);
    ns = 2;
  }else{
    std::cout<<"ERROR in atom_d.py: T2C has wrong shape"<<endl;
    return;
  }
  //auto UC = _UC_.mutable_unchecked<5>();
  if (UC.extent(0)!=l+1 or UC.extent(1)!=nw or UC.extent(2)!=nw){
    std::cout<<"ERROR Uc(l+1,nw,nw) has wrong shape "<<_UC_.request().shape[0]<<","<<_UC_.request().shape[1]<<","<<_UC_.request().shape[2]<<endl;
  }
  
  int shft = 3-l;
  int shft2 = 2*l+1;
  
  bl::Array<complex<double>,2> Sum1(nw,shft2*2), Sum2(nw,shft2*2);
  bl::Array<complex<double>,2> T2Cp(T2C.extent(1),T2C.extent(0));
  T2Cp = bl::conj(T2C.transpose(bl::secondDim,bl::firstDim));
  
  for (int k=0; k<l+1; k++){
    Sum1=0.0;
    for (int i1=0; i1<nw; i1++){
      for (int m4=0; m4<mw; m4++){
	for (int m1=0; m1<mw; m1++){
	  for (int s=0; s<ns; s++)
	    Sum1(i1,m1-m4+shft2) += T2Cp(i1,m4+s*mw)*gck(l,shft+m4,shft+m1,k)*T2C(m1+s*mw,i1);
	}
      }
    }
    Sum2=0.0;
    for (int i2=0; i2<nw; i2++){
      for (int m3=0; m3<mw; m3++){
	for (int m2=0; m2<mw; m2++){
	  for (int s=0; s<ns; s++)
	    Sum2(i2,m3-m2+shft2) += T2Cp(i2,m3+s*mw)*gck(l,shft+m2,shft+m3,k)*T2C(m2+s*mw,i2);
	}
      }
    }
    for (int i1=0; i1<nw; i1++){
      for (int i2=0; i2<nw; i2++){
	complex<double> csum=0.0;
	for (int dm=0; dm<shft2*2; dm++)
	  csum += Sum1(i1,dm)*Sum2(i2,dm);
	UC(k,i1,i2) = csum;
      }
    }
  }
}

py::array_t<double> FastIsing5d(py::array_t<double>& _Eimpc_, py::array_t<double>& _occ_, 
				py::array_t<complex<double>>& _UC_, py::array_t<double>& _FkoJ_)
{
  //Uc(l+1,nw,nw,nw,nw,nw)
  py::buffer_info info_U = _UC_.request();
  bl::Array<complex<double>,5> UC((complex<double>*)info_U.ptr,
				  bl::shape(info_U.shape[0],info_U.shape[1],info_U.shape[2],info_U.shape[3],info_U.shape[4]),
				  bl::neverDeleteData);
  py::buffer_info info_e = _Eimpc_.request();
  bl::Array<double,1> Eimpc((double*)info_e.ptr,info_e.shape[0],bl::neverDeleteData);
  py::buffer_info info_o = _occ_.request();
  bl::Array<double,1> occ((double*)info_o.ptr,info_o.shape[0],bl::neverDeleteData);
  py::buffer_info info_f = _FkoJ_.request();
  bl::Array<double,1> FkoJ((double*)info_f.ptr,info_f.shape[0],bl::neverDeleteData);

  py::array_t<double> _EUterms_({2});
  auto EUterms = _EUterms_.mutable_unchecked<1>();
  EUterms(0) = EUterms(1) = 0;
  
  for (int i=0; i<Eimpc.size(); i++) EUterms(0) += Eimpc(i)*occ(i);
  for (int i=0; i<Eimpc.size(); i++){
    for (int j=0; j<Eimpc.size(); j++){
      double dsum=0.0;
      dsum += (UC(0,i,j,j,i)-UC(0,i,j,i,j)).real()*FkoJ(0);
      for (int k=1; k<FkoJ.size(); k++)
	dsum += (UC(k,i,j,j,i)-UC(k,i,j,i,j)).real()*FkoJ(k);
      EUterms(1) += 0.5*dsum * occ(i)*occ(j);
    }
  }
  return _EUterms_;
}
py::array_t<double> FastIsing3d(py::array_t<double>& _Eimpc_, py::array_t<double>& _occ_, 
				py::array_t<complex<double>>& _UC_, py::array_t<double>& _FkoJ_,
				py::array_t<int>& _bath_)
{
  py::buffer_info info_U = _UC_.request();
  bl::Array<complex<double>,5> UC((complex<double>*)info_U.ptr,
				  bl::shape(info_U.shape[0],info_U.shape[1],info_U.shape[2],info_U.shape[3],info_U.shape[4]),
				  bl::neverDeleteData);
  py::buffer_info info_e = _Eimpc_.request();
  bl::Array<double,1> Eimpc((double*)info_e.ptr,info_e.shape[0],bl::neverDeleteData);
  py::buffer_info info_o = _occ_.request();
  bl::Array<double,1> occ((double*)info_o.ptr,info_o.shape[0],bl::neverDeleteData);
  py::buffer_info info_f = _FkoJ_.request();
  bl::Array<double,1> FkoJ((double*)info_f.ptr,info_f.shape[0],bl::neverDeleteData);
  py::buffer_info info_i = _bath_.request();
  bl::Array<int,1> bath((int*)info_i.ptr,info_i.shape[0],bl::neverDeleteData);
  //cout<<"Eimpc="<<Eimpc<<" occ="<<occ<<" FkoJ="<<FkoJ<<" bath="<<bath<<endl;
  py::array_t<double> _EUterms_({2});
  auto EUterms = _EUterms_.mutable_unchecked<1>();
  EUterms(0) = EUterms(1) = 0;
  
  for (int i=0; i<Eimpc.size(); i++) EUterms(0) += Eimpc(i)*occ(i);
  for (int i=0; i<bath.size(); i++){
    int ii = bath(i);
    int si = (i*2)/bath.size();
    for (int j=0; j<bath.size(); j++){
      int jj = bath(j);
      int sj = (j*2)/bath.size();
      double dsum=0.0;
      for (int k=0; k<FkoJ.size(); k++){
	dsum += UC(k,ii,jj,jj,ii).real()*FkoJ(k);
	if (si==sj) dsum -= UC(k,ii,jj,ii,jj).real()*FkoJ(k);
      }
      EUterms(1) += 0.5*dsum * occ(i)*occ(j);
      //cout<<"i="<<i<<" si="<<si<<" j="<<j<<" sj="<<sj<<" "<<dsum<<" "<<EUterms(1)<<" "<<UC(1,ii,jj,jj,ii).real()<<" "<<UC(1,ii,jj,ii,jj).real()<<endl;
    }
  }
  return _EUterms_;
}

PYBIND11_MODULE(dpybind,m) {
  m.doc() = "pybind11 wrap for small routines accompaning atom_d.py script for exact diagonalization";
  m.def("FromSlaterToMatrixU", &FromSlaterToMatrixU);
  m.def("FromSlaterToMatrixUd", &FromSlaterToMatrixUd);
  m.def("FromSlaterToMatrixU_diagonal", &FromSlaterToMatrixU_diagonal);
  m.def("FastIsing5d", &FastIsing5d);
  m.def("FastIsing3d", &FastIsing3d);
}
