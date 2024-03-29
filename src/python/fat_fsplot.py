#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
""" Plots 2D cut of the Fermi surface with three coherence factors (red/green/blue).
  The documentation was written several years after the code was used, hence 
  some details are possibly wrong. Please verify that it works for your case before using it.

  Using dmftgk you must first prepare eigvals.dat with only one frequency om=0 
  (I think there is first data-point for infinite frequencies followed by omega=0).
  Using fat_coh you also need to prepare coherence factors for red/green/blue color, which 
  are in different files, given in the command line as input.

  The mesh of k-points should be regular mesh in 2D.

  In this particular case only an irreducible piece of the 2D Fermi surface was computed
  and the k-mesh dimension in one direction was assigned to be nk1 = int((-1+sqrt(1+8*nk))/2), 
  where nk is the number of momentum points in the data file. 
  This needs to be adjusted for particular case, as it is not working in general.
"""
from numpy import *
import numpy as np
from pylab import *
import os, sys
import scipy


#code="""
#     using namespace std;
#
#     double Ry2eV = 13.6056920311654;
#
#     double Ax = 0;
#     for (int ib=0; ib<nbands; ib++){
#        //complex<double> ekom( data(1+2*ib),  min(data(2+2*ib), -small) );
#        complex<double> ekom( data(1+2*ib),  -small );
#        complex<double> cohf( dach(1+2*ib),  dach(2+2*ib) );
#        complex<double> gc = cohf/(omega+mu-ekom);
#        Ax += -gc.imag()/M_PI;
#     }
#     return_val = Ax;
#     
#     """
#code0="""
#     using namespace std;
#
#     double Ry2eV = 13.6056920311654;
#
#     double Ax = 0;
#     for (int ib=0; ib<nbands; ib++){
#        //complex<double> ekom( data(1+2*ib),  min(data(2+2*ib), -small) );
#        complex<double> ekom( data(1+2*ib),  -small );
#        complex<double> gc = 1./(omega+mu-ekom);
#        Ax += -gc.imag()/M_PI;
#     }
#     return_val = Ax;
#     
#     """
#code1="""
#     using namespace std;
#     bool found = false;
#     int ib=0;
#     for (ib=0; ib<nbands-1; ib++){
#        double ekom_i = data(1+2*ib);
#        double ekom_ip = data(1+2*ib+2);
#        if ( (omega+mu-ekom_i)*(omega+mu-ekom_ip)<0.0 ){
#           found = true;
#           break;
#        }
#     }
#     double cohf=0;
#     if (found){
#        double ekom_i = data(1+2*ib);
#        double ekom_ip = data(1+2*ib+2);
#        double x_i = omega+mu-ekom_i;
#        double x_ip = omega+mu-ekom_ip;
#        double coh_i  = dach(1+2*ib)*dach(1+2*ib)+dach(2+2*ib)*dach(2+2*ib);
#        double coh_ip = dach(1+2*ib+2)*dach(1+2*ib+2)+dach(2+2*ib+2)*dach(2+2*ib+2);
#        cohf = coh_i + (coh_ip-coh_i)*(0-x_i)/(x_ip-x_i);
#        //cout<<x_i<<" "<<x_ip<<" "<<cohf<<endl;
#     }
#     return_val = cohf;
#     
#     """

def Cols(x0,x1,x2,bright):
    col = array([x0, x1, x2])    # for red,green,blue
    col *= 1./sqrt(dot(col,col)) # normalized to lengh 1.
    col *= bright                # and brightness is determined by Akw/Akmax
    col = array([1-sqrt(col[1]**2+col[2]**2),1-sqrt(col[0]**2+col[2]**2),1-sqrt(col[0]**2+col[1]**2)])
    return col

if __name__ == '__main__':

    small = 1.e-2
    mu = float(open('EF.dat', 'r').readline().split()[0])
    print('mu=', mu)
    
    fdat = open('eigvals.dat', 'r')
    #if os.path.isfile('cohfactorsd.dat'):
    #    fcoh = open('cohfactorsd.dat', 'r')
    #else:
    #    fcoh = None

    print(sys.argv)

    if len(sys.argv)<4:
        print('Please give 3 files of coherence factors')
        sys.exit(0)
        
    print(sys.argv[1])
    print(sys.argv[2])
    print(sys.argv[3])
    
    fcoh1 = open(sys.argv[1], 'r')
    fcoh2 = open(sys.argv[2], 'r')
    fcoh3 = open(sys.argv[3], 'r')
    
    wdata = fdat.readlines()

    ikp=0
    Akom=[]
    Ch1kw=[]
    Ch2kw=[]
    Ch3kw=[]
    ii=0
    while (ii<len(wdata)): # over k-points
        
        data0 = wdata[ii].split()
        dach1 = next(fcoh1).split()
        dach2 = next(fcoh2).split()
        dach3 = next(fcoh3).split()
        ii += 1
        #print data
        (ikp, isym, nbands, nemin, nomega) = list(map(int, data0[1:6]))

        ekom = zeros(nbands, dtype=complex)
        
        Aom=zeros(nomega)
        Ch1w=zeros(nomega)
        Ch2w=zeros(nomega)
        Ch3w=zeros(nomega)
        om=zeros(nomega)
        for iom in range(nomega):
            data = array([float(x) for x in wdata[ii].split()])
            dach1 = array([float(x) for x in next(fcoh1).split()])
            dach2 = array([float(x) for x in next(fcoh2).split()])
            dach3 = array([float(x) for x in next(fcoh3).split()])
            ii += 1
            
            omega = float(data[0])
            om[iom] = omega
            
            #Ax = weave.inline(code, ['nbands', 'omega', 'mu', 'data', 'small', 'ikp', 'dach'],
            #                  type_converters=weave.converters.blitz, compiler = 'gcc')
            #Ax = weave.inline(code0, ['nbands', 'omega', 'mu', 'data', 'small', 'ikp'],
            #                  type_converters=weave.converters.blitz, compiler = 'gcc')
            #dach = dach1
            #Ch1 = weave.inline(code1, ['nbands', 'omega', 'mu', 'data', 'small', 'ikp', 'dach'],
            #                  type_converters=weave.converters.blitz, compiler = 'gcc')
            #dach = dach2
            #Ch2 = weave.inline(code1, ['nbands', 'omega', 'mu', 'data', 'small', 'ikp', 'dach'],
            #                  type_converters=weave.converters.blitz, compiler = 'gcc')
            #dach = dach3
            #Ch3 = weave.inline(code1, ['nbands', 'omega', 'mu', 'data', 'small', 'ikp', 'dach'],
            #                  type_converters=weave.converters.blitz, compiler = 'gcc')

            ekom = data[1::2]
            xi = omega+mu-ekom
            # code0
            Aom[iom] = -sum( (1.0/(xi+small*1j)).imag )/M_PI
            whr = where( np.diff(np.sign(xi)) != 0 )[0]  # finds the Fermi surface points
            if whr:
                cohf1 = dach1[1::2] + dach1[2::2]*1j
                cohf2 = dach2[1::2] + dach2[2::2]*1j
                cohf3 = dach3[1::2] + dach3[2::2]*1j
                # compute coherence factor with linear approximation for the point on the fermi surface
                Ch1w[iom] = sum([abs(cohf1[j])**2 + (abs(cohf1[j+1])**2 - abs(cohf1[j])**2)*(0.0-xi[j])/(xi[j+1]-xi[j]) for j in whr])
                Ch2w[iom] = sum([abs(cohf2[j])**2 + (abs(cohf2[j+1])**2 - abs(cohf2[j])**2)*(0.0-xi[j])/(xi[j+1]-xi[j]) for j in whr])
                Ch3w[iom] = sum([abs(cohf3[j])**2 + (abs(cohf3[j+1])**2 - abs(cohf3[j])**2)*(0.0-xi[j])/(xi[j+1]-xi[j]) for j in whr])
        
        Akom.append( Aom )
        Ch1kw.append( Ch1w )
        Ch2kw.append( Ch2w )
        Ch3kw.append( Ch3w )
        ikp += 1

    Akom = array(Akom).T   # Akom[w,k]
    Ch1kw = array(Ch1kw).T # Ch1kw[w,k]
    Ch2kw = array(Ch2kw).T # Ch2kw[w,k]
    Ch3kw = array(Ch3kw).T # Ch3kw[w,k]
    Akomax=max(Akom[0,:])  # over all k

    nk = shape(Akom)[1]
    nk1 = int((-1+sqrt(1+8*nk))/2)
    
    # From now on
    FS = zeros((nk1,nk1,3), dtype=float)
    ikk=0
    for i in range(nk1):
        for j in range(i,nk1):
            col = Cols(Ch1kw[0,ikk], Ch2kw[0,ikk], Ch3kw[0,ikk], Akom[0,ikk]/Akomax)
            
            FS[i,j,:] = col
            FS[j,i,:] = col
            ikk+=1
    
    
    #vmm = [0,20]
    #xmm = [0, shape(Akom)[1]]
    #ymm = [om[0], om[-1]]
    
    #imshow(Akom, interpolation='bilinear', cmap=cm.hot, origin='lower', vmin=vmm[0], vmax=vmm[1], extent=[xmm[0],xmm[1],ymm[0],ymm[1]], aspect=(xmm[1]-xmm[0])*0.8/(ymm[1]-ymm[0]) )
    #imshow(FS, interpolation='bilinear', cmap=cm.hot, origin='lower' , vmin=0.0, vmax=80.)
    print(shape(FS))
    
    imshow(FS, origin='lower',  interpolation='bilinear', extent=[0,1,0,1], aspect=1. )
    
    show()
    
    
    
