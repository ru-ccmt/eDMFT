#!/usr/bin/env python

# @Copyright 2007 Kristjan Haule
from numpy import *
from pylab import *
import glob, os, sys
import scipy
import cakw

from matplotlib.colors import colorConverter
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib import colormaps

if __name__ == '__main__':
    small = 1e-4
    DY = 0
    
    if len(sys.argv)<2:
        intensity = 0.97
    else:
        intensity = float(sys.argv[1])

    fname1 = 'eigvals.dat'
    fname2 = 'eigvalsdn.dat'
    if len(sys.argv)>2:
        fname1 = sys.argv[2]
    if len(sys.argv)>3:
        fname2 = sys.argv[3]
    print('using', fname1, fname2)

    if os.path.isfile(fname1):
        fdat1 = open(fname1, 'r')
    elif os.path.isfile(fname1+'.gz'):
        import gzip
        fdat1 = gzip.open(fname1+'.gz', 'rt')
    else:
        print('ERROR: file', fname1, 'does not exist!')
        sys.exit(1)
    if os.path.isfile(fname2):
        fdat2 = open(fname2, 'r')
    elif os.path.isfile(fname2+'.gz'):
        import gzip
        fdat2 = gzip.open(fname2+'.gz', 'rt')
    else:
        print('ERROR: file', fname2, 'does not exist!')
        sys.exit(1)


    _cmap_ = cm.hot # color map from matplotlib
    _col_ = 'k'
    
    fEF = open('EF.dat', 'r')
    mu = float(fEF.read())

    print('mu=', mu)

    rxy = array(loadtxt('kindex.dat'), dtype=int)
    Nx, Ny = rxy[-1][:2]
    
    with open('kindex.dat','r') as fkin:
        fkin.readline()
        line = fkin.readline()[2:]
        xi = eval(line.split(' ',1)[1])
        line = fkin.readline()[2:]
        yi = eval(line.split(' ',1)[1])
        print('xi=', xi, 'yi=', yi)
    
    
    ikp=0
    Akm1=[]
    Akm2=[]
    while True:
        line1 = fdat1.readline()
        line2 = fdat2.readline()
        if not line1 or not line2:
            break
        data = line1.split()
        (ikp, isym, nbands, nemin, nomega) = list(map(int, data[1:6]))
        
        ekom = zeros(nbands, dtype=complex)
        dach=ones((nomega,nbands), dtype=complex)
        index=list(range(nomega))
        
        for iom in range(nomega):
            data1 = array(list(map(float, fdat1.readline().split())))
            data2 = array(list(map(float, fdat2.readline().split())))
            if iom==0:
                ekom1 = data1[1::2]+data1[2::2]*1j
                ekom2 = data2[1::2]+data2[2::2]*1j
                omega = float(data1[0])
                cohd=ones(nbands, dtype=complex)
                Aom1 = cakw.Akw(nbands,omega,mu,ekom1,cohd,small)
                Aom2 = cakw.Akw(nbands,omega,mu,ekom2,cohd,small)
        Akm1.append( Aom1 )
        Akm2.append( Aom2 )

    Akm1 = array(Akm1)
    Akm2 = array(Akm2)
    
    ht, bin_edges = histogram(vstack((Akm1.ravel(),Akm2.ravel())),bins=5000)
    xh = 0.5*(bin_edges[1:]+bin_edges[:-1])
    cums = cumsum(ht)/sum(ht)
    i = searchsorted(cums, intensity)
    print('with intensity=', intensity, 'we determine cutoff at Ak=', xh[i])
    vmm = [0,xh[i]]
    
    print('vm=', vmm, shape(Akm1), shape(Akm2) )
    Akom1 = zeros((Nx+1,Ny+1))
    Akom2 = zeros((Nx+1,Ny+1))
    for i,j,ii in rxy:
        Akom1[i,j] = Akm1[ii]
        Akom2[i,j] = Akm2[ii]
    Akom1 = Akom1.T
    Akom2 = Akom2.T
    
    #print(shape(Akom1), shape(Akom2))
    #plot(Akom1[:,49])
    #plot(Akom2[:,49])
    #plot(Akom1[:,49]-Akom2[:,49])
    #show()

    
    subplot(221)
    imshow(Akom1, cmap=colormaps['Reds'], vmin=vmm[0], vmax=vmm[1], extent=[xi[0],xi[1],yi[0],yi[1]],aspect=1.0)
    xticks([]);yticks([])
    subplot(223)
    imshow(Akom2, cmap=colormaps['Blues'], vmin=vmm[0], vmax=vmm[1], extent=[xi[0],xi[1],yi[0],yi[1]],aspect=1.0)
    xticks([]);yticks([])
    subplot(222)
    alphas = 1-exp(-3*(abs(Akom1)+abs(Akom2))/vmm[1])
    #alphas = ones(shape(Akom1))
    imshow(Akom1-Akom2, alpha=alphas, interpolation='bilinear', cmap=colormaps['coolwarm'],
               vmin=-vmm[1], vmax=vmm[1], extent=[xi[0],xi[1],yi[0],yi[1]],aspect=1.0)
    xticks([]);yticks([])
    show()
    

