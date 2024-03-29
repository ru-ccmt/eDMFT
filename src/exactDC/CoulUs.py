# @Copyright 2007 Kristjan Haule
from numpy import zeros,array
import gaunt
import dpybind
#import scipy

def CoulUsC2(l, T2C):
    # Precomputes gaunt coefficients for speed
    gck = gaunt.cmp_all_gaunt()
    gck = array(gck, order='C')   # This is essential for pybind11 code to work with fortran gck
    T2C = array(T2C, order='C')
    nw = len(T2C)
    UC = zeros((l+1, nw, nw, nw, nw), dtype=complex)
    dpybind.FromSlaterToMatrixU(UC, gck, l, T2C)
    return UC

def CoulUsC2_diagonal(l, T2C):
    # Precomputes gaunt coefficients for speed
    # shape(gck)=(4,7,7,4). It contains gck(l, m4, m1, k)
    gck = gaunt.cmp_all_gaunt()
    gck = array(gck, order='C')   # This is essential for pybind11 code to work with fortran gck
    T2C = array(T2C, order='C')
    nw = len(T2C)
    UC = zeros((l+1, nw, nw), dtype=complex)
    dpybind.FromSlaterToMatrixU_diagonal(UC, gck, l, T2C)
    return UC


