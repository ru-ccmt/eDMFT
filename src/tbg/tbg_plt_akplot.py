#!/usr/bin/env python
# @Copyright 2024 Kristjan Haule
import numpy as np
from pylab import *
import glob, os, sys
import scipy
import argparse

mingle_names={'W':'$W$','L':'$L$','LAMBDA':'$\Lambda$','GAMMA':'$\Gamma$','DELTA':'$\Delta$','X':'$X$','Z':'$Z$','W':'$W$','K':'$K$',
              'R': '$R$', 'S':'$S$', 'T':'$T$', 'U':'$U$', 'Y':'$Y$'}

def ReadDataFile(fname):
    try:
        ekw = np.loadtxt(fname)
    except ValueError:
        # Fallback: parse line by line
        data_list = []
        with open(fname, "r") as f:
            for line in f:
                # Strip whitespace and skip empty or comment lines if needed
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                # Split and convert each field to float
                fields = line.split()
                row = [float(x) for x in fields]
                data_list.append(row)
        min_columns = np.min([len(line) for line in data_list])
        data_list = [line[:min_columns] for line in data_list]
        eks = np.array(data_list)
    return eks
    
if __name__ == '__main__':

    orb_plot=[0,0,0,0]
    
    usage = 'Plots the spectral furnction after dmftp step has been executed'
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('fname', nargs='?', default='eigvals.dat', type=str, help='filename eigvals.dat or eigenvalues.dat. Default=eigvals.dat')
    parser.add_argument('-i', type=float, default=0.97, help='color intensity, a number slightly smaller than 1.(default 0.97). Percentage of points being used to find maximum value.')
    parser.add_argument('-b', type=float, default=1e-4, help='small broadening in calculating A(k,w), default 1e-5')
    parser.add_argument('-d', type=float, default=0, help='shift of zero in y axis when there is a gap and fermi level can be moved away from zero')
    parser.add_argument('-c', type=str, default='viridis', help='color map, default is cm.hot but could be changed to cm.Purples or any other matplotlib color map')
    parser.add_argument('-l', type=str, default='k', help='color of the lines.default w for white')
    parser.add_argument('-g', default=False, action='store_true', help='add color bar')
    parser.add_argument('-a', type=float, default=0.4, help='aspect ratio of the plot, default 0.8')
    args = parser.parse_args()
    
    #fname = args.fname
    intensity = args.i
    small = args.b
    DY = args.d
    _cmap_ = args.c #eval(args.c)
    _col_ = args.l



    exec(open("params.py").read())
    unit = par_bs['unit']
    band = '_band'
    path = par_bs['path']

    
    npzfile = load('zek.npz')
    zek, om, mu = npzfile['zek'], npzfile['w'], npzfile['mu']
    
    npz3file = load('kpath'+band+'.npz')
    kpath, Npth = npz3file['kpath'], npz3file['Npth']
    Nband, nw, Nkp = shape(zek)

    print('Nkp=', Nkp)
    print('Npth=', Npth)

    wkpointi = [sum(Npth[:i]) for i in range(len(Npth))]+[len(kpath)]
    wkpoints = ['$'+x+'$' for x in par_bs['path_legend']]
    
    print('high-symmetry points found', [(wkpointi[i],wkpoints[i]) for i in range(len(wkpoints))] )

    npzfile = load('Uproject_band.npz')
    UAA, UAB = npzfile['UAA'], npzfile['UAB']

    # coh[ik,iband]=UAA[ifc,ik,iband]
    coh0 = transpose(abs(UAA[:,:,:])**2)  # [:nband,:nkp,:ifc]

    cohi = zeros((Nband,Nkp,4))
    cohi[:,:,3]=1.0  # alpha
    
    unique_vals = np.unique(np.array(orb_plot))
    unique_channels = unique_vals[unique_vals <= 2]
    
    for c in unique_channels:
        # Find all i where orb_plot[i] == c
        idx = np.where(orb_plot == c)[0]
        # Sum across the third axis for all those i
        cohi[:, :, c] = np.sum(coh0[:, :, idx], axis=2)

    Nfunc, Nik, Nband = shape(UAA)
    print('Nik=', Nik)
    
    
    nkp = wkpointi[-1]+1
    print('nkp=', nkp)
    
    nbnd,nom2,nkp2 = shape(zek)
    
    # ensure causality
    zek = np.where(zek.imag < -small, zek, zek.real - small * 1j)

    # Akom[:nom,:nkp,4]
    # cohi[:nbnd,:nkp,4]
    # zek[:nbnd,:nom,:nkp]
    #                inside [:nbnd,:nom,:nkp,4]
    Akom = np.sum( (cohi[:,None,:,:]/(om[None,:,None,None] + mu - zek[:,:,:,None])).imag, axis=0)*(-1/pi)
    
    # finds how to scale the image, and what is the maximum value for density plot
    ht, bin_edges = histogram(Akom[:,:,3].ravel(),bins=5000)
    xh = 0.5*(bin_edges[1:]+bin_edges[:-1])
    cums = cumsum(ht)/sum(ht)
    i = searchsorted(cums, intensity)
    print('with intensity=', intensity, 'we determine cutoff at max(Ak[:,:])=', xh[i])
    vmm = [0,xh[i]]

    Akom[:,:,0] /= Akom[:,:,3]
    Akom[:,:,1] /= Akom[:,:,3]
    Akom[:,:,2] /= Akom[:,:,3]
    Akom[:,:,3] /= vmm[1]
    
    #plot(xh, ht,'.-')
    #show()
    #sys.exit(0)
    
    (ymin,ymax) = (om[0]+DY,om[-1]+DY)
    (xmin,xmax) = (0, nkp-1)
    
    Akom/=vmm[1]
    Akom[Akom > 1.0] = 1.0
    imshow(Akom, interpolation='bilinear', origin='lower', extent=[xmin,xmax,ymin*unit,ymax*unit], aspect=(xmax-xmin)*args.a/(unit*(ymax-ymin)) )
    
    for i in range(len(wkpointi)):
        plot([wkpointi[i],wkpointi[i]], [ymin*unit,ymax*unit], _col_+'-')
        
    plot([xmin,xmax],[0,0], _col_+':')

    dytck=0.005
    Ntck=5
    for j in range(len(wkpointi)-1):
        for ix in range(1,Ntck):
            x = wkpointi[j]+(wkpointi[j+1]-wkpointi[j])*ix/float(Ntck)
            plot([x,x],[-dytck*unit,dytck*unit], _col_+'-')
        
    axis([xmin,xmax,ymin*unit,ymax*unit])
    xticks( wkpointi, wkpoints, fontsize='x-large' )
    if args.g:
        colorbar()
    ylabel('Energy[meV]')
    show()
    
