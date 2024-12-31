#!/usr/bin/env python
# @Copyright 2024 Kristjan Haule
import numpy as np
from pylab import *
import glob, os, sys
import argparse
from matplotlib.colors import colorConverter
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import plt_auxiliary as au

mingle_names={'W':'$W$','L':'$L$','LAMBDA':'$\Lambda$','GAMMA':'$\Gamma$','DELTA':'$\Delta$','X':'$X$','Z':'$Z$','W':'$W$','K':'$K$',
              'R': '$R$', 'S':'$S$', 'T':'$T$', 'U':'$U$', 'Y':'$Y$'}

if __name__ == '__main__':
    usage = 'Plots the spectral furnction after dmftp step has been executed'
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('fname', nargs='?', default='eigvals.dat', type=str, help='filename eigvals.dat or eigenvalues.dat. Default=eigvals.dat')
    parser.add_argument('-i', type=float, default=0.97, help='color intensity, a number slightly smaller than 1.(default 0.97). Percentage of points being used to find maximum value.')
    parser.add_argument('-b', type=float, default=1e-5, help='small broadening in calculating A(k,w), default 1e-5')
    parser.add_argument('-d', type=float, default=0, help='shift of zero in y axis when there is a gap and fermi level can be moved away from zero')
    parser.add_argument('-c', type=str, default='cm.hot', help='color map, default is cm.hot but could be changed to cm.Purples or any other matplotlib color map')
    parser.add_argument('-l', type=str, default='k', help='color of the lines.default k for white')
    parser.add_argument('-g', default=False, action='store_true', help='add color bar')
    args = parser.parse_args()
    
    fname1 = args.fname
    fbase,fend = os.path.splitext(fname1)
    fname2 = fbase+'dn'+fend
    intensity = args.i
    small = args.b
    DY = args.d
    _cmap_ = eval(args.c)
    _col_ = args.l
    
    mu = au.FindEF()
    
    fnames = [fname1,fname2]
    print('using fnames=',fnames, 'and chemical potential', mu)
    #print('cmap=', args.c, 'color=', args.l, 'intensity=', args.i, 'small=', args.b, 'DY=', args.d)
    
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
        print('high-symmetry points found', [(wkpointi[i],wkpoints[i]) for i in range(len(wkpoints))] )
    
    nkp = wkpointi[-1]+1
    print('nkp=', nkp)

    Akw=[]
    for spin in [0,1]:
        ekw = loadtxt(fnames[spin])
        nom = int(len(ekw)/nkp)
        if nkp*nom != len(ekw):
            print('ERROR data length in', fname1, 'seems to be incompatible with case.klist_band. We have nkp=', nkp, 'and num lines=', len(ekw), 'which is incomensurate with number of k-points')
            sys.exit(1)
        om = ekw[:nom,0]
        zekw = ekw[:,1::2]+ekw[:,2::2]*1j
        nbnd = shape(zekw)[1]
        zekw = reshape(zekw, (nkp,nom,nbnd))
        # ensure causality
        zekw = np.where(zekw.imag < -small, zekw, zekw.real - small * 1j)
        # Akom[ikp,iom]
        Akom = np.sum( (1.0/(om[np.newaxis,:,np.newaxis] + mu - zekw)).imag, axis=2)*(-1/pi)
        Akw.append(Akom)
        
    # finds how to scale the image, and what is the maximum value for density plot
    ht, bin_edges = histogram(vstack((Akw[0].ravel(),Akw[1].ravel())),bins=5000)
    xh = 0.5*(bin_edges[1:]+bin_edges[:-1])
    cums = cumsum(ht)/sum(ht)
    i = searchsorted(cums, intensity)
    print('with intensity=', intensity, 'we determine cutoff at max(Ak[:,:])=', xh[i])
    vmm = [0,xh[i]]
    
    (ymin,ymax) = (om[0]+DY,om[-1]+DY)
    (xmin,xmax) = (0, nkp-1)
    
    # make custom the colormaps
    cmap1 = mpl.colors.LinearSegmentedColormap.from_list('my_cmap1',['white','blue','navy'],N=256,gamma=1.0)
    cmap2 = mpl.colors.LinearSegmentedColormap.from_list('my_cmap2',['white','red','crimson'],N=256,gamma=1.0)
    # create the _lut array, with rgba values
    cmap2._init() 
    # create alpha array and fill the colormap with them.
    # here it is progressive, but you can create whathever you want
    alphas = np.linspace(0, 0.9, cmap2.N+3)
    cmap2._lut[:,-1] = alphas
    
    #print('xmin,xmax,ymin,ymax=', xmin, xmax, ymin, ymax)
    aspect=(xmax-xmin)*0.5/(ymax-ymin)
    imshow(Akw[0].T, interpolation='bilinear', cmap=cmap1, origin='lower', vmin=vmm[0], vmax=vmm[1], extent=[xmin,xmax,ymin,ymax],aspect=aspect)
    imshow(Akw[1].T, interpolation='bilinear', cmap=cmap2, origin='lower', vmin=vmm[0], vmax=vmm[1], extent=[xmin,xmax,ymin,ymax],aspect=aspect)
    
    for i in range(len(wkpointi)):
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
    if args.g:
        colorbar()
    show()
    
