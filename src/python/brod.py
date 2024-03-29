#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
from scipy import *
from scipy import integrate, interpolate
import brd
import time
from numba import njit

def MakeTanMesh(N, tanc, tanw, b0, b1):
    if not(b0<b1): print("Relation must hold: b0<b1!")
    if not(b0<tanw and tanw<b1): print("Relation mesu hold: b0<tanw<b1!")
    if not(b0>0): print("b0 must be positive!")
    du = arctan(((tanc-b0)/tanw))
    b1n = arctan((b1-tanc)/tanw)+du
    m0 = [tanc + tanw * tan(b1n*(i-1)/(N-2)-du) for i in range(1,N)]
    return hstack( (-array(m0[::-1]), array([0]+m0) ) )

def Broad(width, kwidth, om, fw):
    " Broadens the data with gaussian of width=width"
    
    fwi = interpolate.interp1d(om, fw)
    fwn=zeros(len(om), fwi.dtype)

    for im in range(len(om)):
        w=width + kwidth*abs(om[im])
        if (om[im]-om[0]>w*4 and om[-1]-om[im]>w*4):  # Gaussian is fully within existing mesh
            x = brd.maketanmesh(200,0.0,w,w/50,w*20)
            x2,ni = brd.combinemesh(om[im],om,x)
            eps = x2[:ni]
            x3 = om[im]-eps
            tw = (2*w**2)
            gs = exp(-x3**2/tw)/(sqrt(2*pi)*w)
            norm = integrate.trapz(gs,x=x3)
            yn = integrate.trapz(fwi(eps) * gs, x=eps)/abs(norm)
        else:
            yn = fw[im]
        fwn[im] = yn
    return fwn

if __name__ == '__main__':
    import sys
    from numpy import *
    import optparse
    
    usage = """usage: %broad.py [ options ] filename

    Broadens existing datafile with gaussian broadening. 
    First column is assumed to be x, and all other columns are ys.
    """
    parser = optparse.OptionParser(usage)
    parser.add_option('-w',   dest='w', type=float, default=0.05, help="constant broadening of gaussian width w")
    parser.add_option('-k',   dest='k', type=float, default=0.1,  help="linearly increasing broadening width~ w+|x|*k")
    # Next, parse the arguments
    (options, args) = parser.parse_args()
    if len(args)!=1:
        print('Need exactly one argument: the name of input file to broaden')
        sys.exit(1)
    fname = args[0]
    
    data = loadtxt(fname).T
    om = data[0]

    fg = data[1::2,:] + data[2::2,:]*1j
    y=zeros(shape(fg), dtype=complex)
    for i in range(len(fg)):
        y[i,:] = Broad(options.w, options.k, om, fg[i])

    for iw in range(shape(y)[1]):
        print(om[iw], end=' ')
        for i in range(len(y)):
            print(y[i,iw].real, y[i,iw].imag, end=' ')
        print()
        
    
