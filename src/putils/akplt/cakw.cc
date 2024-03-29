#include <complex>
#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include "pybind11/stl.h"
#include <cstdint>

namespace py = pybind11;
using namespace std;

double Akw(int nbands, double omega, double mu, py::array_t<std::complex<double>>& ekom, py::array_t<std::complex<double> >& cohd, double small)
{
  auto _ekom_ = ekom.mutable_unchecked<1>();
  auto _cohd_ = cohd.mutable_unchecked<1>();

  double Ax = 0;
  //#pragma omp parallel for
  for (int ib=0; ib<nbands; ib++){
    complex<double> ekw=_ekom_(ib);
    if (ekw.imag() > -small) ekw=complex<double>(ekw.real(),-small);
    complex<double> gc = abs(_cohd_(ib))/(omega+mu-ekw);
    Ax += -gc.imag()/M_PI;
  }
  return Ax;
}

double Akw0(int nbands, double omega, py::array_t<std::complex<double>>& ekom, double small)
{
  auto _ekom_ = ekom.mutable_unchecked<1>();
  
  double Ax = 0;
  //#pragma omp parallel for
  for (int ib=0; ib<nbands; ib++){
    complex<double> ekw=_ekom_(ib);
    if (ekw.imag() > -small) ekw=complex<double>(ekw.real(),-small);
    complex<double> gc = 1.0/(omega-ekw);
    Ax += -gc.imag()/M_PI;
  }
  return Ax;
}

void Akw_orb(py::array_t<double>& Am, int nbands, double omega, py::array_t<std::complex<double>>& ekom, py::array_t<std::complex<double> >& cohd, double small)
{
  auto _ekw_ = ekom.mutable_unchecked<1>();
  auto _cohd_ = cohd.mutable_unchecked<2>();
  auto _Am_ = Am.mutable_unchecked<1>();
  
  double w=omega;
  int norb=_Am_.shape(0);
  //#pragma omp parallel for
  for (int ib=0; ib<nbands; ib++){
    complex<double> ew=_ekw_(ib);
    if (ew.imag() > -small) ew=complex<double>(ew.real(),-small);
    for (int iorb=0; iorb<norb; iorb++){
      double Ac = -abs(_cohd_(iorb,ib))*(1./(w-ew)).imag()/M_PI;
      _Am_(iorb) += Ac;
    }
  }
}


PYBIND11_MODULE(cakw,m) {
  m.doc() = "pybind11 wrap for small code that computes Akw";
  m.def("Akw", &Akw);
  m.def("Akw0", &Akw0);
  m.def("Akw_orb", &Akw_orb);
}
