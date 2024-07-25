#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
from scipy import *
from pylab import *
import glob, os, sys
import scipy
import cakw

from matplotlib.colors import colorConverter
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl


mingle_names={'W':'$W$','L':'$L$','LAMBDA':'$\Lambda$','\Gamma':'$\Gamma$','GAMMA':'$\Gamma$','DELTA':'$\Delta$','X':'$X$','Z':'$Z$','W':'$W$','K':'$K$',
              'R': '$R$', 'S':'$S$', 'T':'$T$', 'U':'$U$', 'Y':'$Y$'}

if __name__ == '__main__':

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

        
    small = 1e-4 # 0.01 # 1e-5
    #itensity = 0.2
    DY = 0 # 0.01318

    # colors
    if True:
        _cmap_ = cm.hot # color map from matplotlib
        _col_ = 'w'     # lines are of this color
    else:
        _cmap_ = cm.Purples
        _col_ = 'k'
    _col_='k'
    
    fEF = open('EF.dat', 'r')
    mu = float(fEF.read())

    print('mu=', mu)


    wg = glob.glob('*.klist_band')
    if len(wg)>0:
        fg = open(wg[0], 'r')
        wkpointi=[]
        wkpoints=[]
        for il,line in enumerate(fg):
            if line[:3]=='END': break
            com = line[:10].split()
            if com:
                legnd=line.split()[0]
                if legnd in mingle_names:
                    legnd = mingle_names[legnd]
                else:
                    legnd = '$'+legnd+'$'
                wkpoints.append(legnd)
                wkpointi.append(il)
        print(wkpointi)
        print(wkpoints)

    nkp = wkpointi[-1]+1
    print('nkp=', nkp)
    
    if os.path.isfile('cohfactorsd.dat'):
        fcoh = open('cohfactorsd.dat', 'r')
    else:
        fcoh = None
    
    ikp=0
    Akom1=[]
    Akom2=[]
    while True:
        line1 = fdat1.readline()
        line2 = fdat2.readline()
        if not line1 or not line2:
            break
        data = line1.split()
        if fcoh is not None:
            #dach = fcoh.next().split()
            line = fcoh.readline()
            if not line:
                break
            dach = line.split()
        
        (ikp, isym, nbands, nemin, nomega) = list(map(int, data[1:6]))
        
        ekom = zeros(nbands, dtype=complex)
        dach=ones((nomega,nbands), dtype=complex)
        index=list(range(nomega))
        omw=zeros(nomega,dtype=float)
        if fcoh is not None:
            for iom in range(nomega):
                datc = array(list(map(float,fcoh.readline().split())))
                omw[iom] = datc[0]
                dach[iom,:] = datc[1::2]+datc[2::2]*1j
                #print 'shape=', shape(dach), 'nbands=', nbands
            # need to sort frequency because open-mp mixes them up
            index=sorted(index, key=lambda i: omw[i])
            #for i in range(len(index)):
            #    print omw[index[i]],
            #print

        Aom1=zeros(nomega,dtype=float)
        Aom2=zeros(nomega,dtype=float)
        om=zeros(nomega,dtype=float)
        for iom in range(nomega):
            data1 = array(list(map(float, fdat1.readline().split())))
            ekom1 = data1[1::2]+data1[2::2]*1j
            data2 = array(list(map(float, fdat2.readline().split())))
            ekom2 = data2[1::2]+data2[2::2]*1j
            omega = float(data1[0])
            om[iom] = omega
            cohd = dach[index[iom]]
            Aom1[iom] = cakw.Akw(nbands,omega,mu,ekom1,cohd,small)
            Aom2[iom] = cakw.Akw(nbands,omega,mu,ekom2,cohd,small)

        Akom1.append( Aom1 )
        Akom2.append( Aom2 )
            
    
    Akom1 = array(Akom1).T
    Akom2 = array(Akom2).T
    print('shape(Akom1),shape(Akom2)=', shape(Akom1), shape(Akom2))

    #vmm = [0,max(list(map(max,Akom1)))*intensity]    

    ht, bin_edges = histogram(vstack((Akom1.ravel(),Akom2.ravel())),bins=5000)
    xh = 0.5*(bin_edges[1:]+bin_edges[:-1])
    cums = cumsum(ht)/sum(ht)
    i = searchsorted(cums, intensity)
    print('with intensity=', intensity, 'we determine cutoff at Ak=', xh[i])
    vmm = [0,xh[i]]
    
    #hist(Akom1.ravel(), bins='auto')
    #show()
    
    (ymin,ymax) = (om[0]+DY,om[-1]+DY)
    (xmin,xmax) = (0, shape(Akom1)[1]-1)
    #(xmin,xmax) = (0, nkp-1)
    
    print('xmin,xmax,ymin,ymax=', xmin, xmax, ymin, ymax)
    
    # make custom the colormaps
    cmap1 = mpl.colors.LinearSegmentedColormap.from_list('my_cmap1',['white','blue','navy'],N=256,gamma=1.0)
    cmap2 = mpl.colors.LinearSegmentedColormap.from_list('my_cmap2',['white','red','crimson'],N=256,gamma=1.0)
    # create the _lut array, with rgba values
    cmap2._init() 
    # create alpha array and fill the colormap with them.
    # here it is progressive, but you can create whathever you want
    alphas = np.linspace(0, 0.9, cmap2.N+3)
    cmap2._lut[:,-1] = alphas

    aspect=(xmax-xmin)*0.5/(ymax-ymin)
    imshow(Akom1, interpolation='bilinear', cmap=cmap1, origin='lower', vmin=vmm[0], vmax=vmm[1], extent=[xmin,xmax,ymin,ymax],aspect=aspect)
    imshow(Akom2, interpolation='bilinear', cmap=cmap2, origin='lower', vmin=vmm[0], vmax=vmm[1], extent=[xmin,xmax,ymin,ymax],aspect=aspect)

    for i in range(len(wkpointi)):
        print('wp=', wkpointi[i])
        plot([wkpointi[i],wkpointi[i]], [ymin,ymax], _col_+'-')
        
    plot([xmin,xmax],[0,0], _col_+':')

    dytck=0.005
    Ntck=5
    for j in range(len(wkpointi)-1):
        for ix in range(1,Ntck):
            x = wkpointi[j]+(wkpointi[j+1]-wkpointi[j])*ix/float(Ntck)
            plot([x,x],[-dytck,dytck], _col_+'-')
        
    axis([xmin,xmax,ymin,ymax])
    xticks( wkpointi, wkpoints, fontsize='x-large' )
    #colorbar()
    show()
    
