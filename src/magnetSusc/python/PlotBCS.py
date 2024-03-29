#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
from scipy import *
from pylab import *
from scipy import interpolate

fname = sys.argv[1]

for i in range(1,2):
    #gs = loadtxt('gs_symmetryr.'+str(i) ).transpose()
    gs = loadtxt( fname ).transpose()
    Max = max(map(max,gs))
    Min = min(map(min,gs))
    MM = max(abs(Min),Max)
    (norb_norb,nk2) = shape(gs)
    nk1=int(round(sqrt(nk2)))
    norb=int(round(sqrt(norb_norb)))
    gs_=gs.reshape((norb_norb,nk1,nk1))
    gs2 = zeros((norb_norb,nk1+1,nk1+1))
    gs2[:,:nk1,:nk1] = gs_[:,:,:]
    gs2[:,nk1,:nk1] = gs_[:,0,:]
    gs2[:,:nk1,nk1] = gs_[:,:,0]
    gs2[:,nk1,nk1] = gs_[:,0,0]
    print 'Min=', Min, 'Max=', Max, 'MM=', MM
    for j in range(norb):
        subplot(2,norb/2+1,j+2)
        
        imshow( gs2[j*norb+j], origin='lower',  interpolation='bilinear', extent=[0,2,0,2], aspect=1., vmin=-MM,vmax=MM )
        #imshow( gs2[j*norb+j], origin='lower',  interpolation='bilinear', extent=[0,2,0,2], aspect=1.)
        
        #subplot(2,norb,norb+j)
        #kxy = arange(nk1+1.)/float(nk1)
        
        #fgs = interpolate.RectBivariateSpline(kxy, kxy, gs2[j*norb+j,:,:], kx=2,ky=2,s=0 )
        #x=linspace(0,1,100)
        #print 'shape(fgs(x,0.0))=', shape(fgs(x,0.0))
        #print 'shape(fgs(0.0,x))=', shape(fgs(0.0,x))
        #plot(x, fgs(0.0,x).ravel())
        #plot(x, fgs(0.5,x).ravel())
        
        #plot(kxy, gs2[j*norb+j,0,:], 'o')
        #plot(kxy, gs2[j*norb+j,nk1/2,:], 'o')
        
        colorbar()
    show()
