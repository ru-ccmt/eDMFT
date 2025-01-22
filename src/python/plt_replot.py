#!/usr/bin/env python
import numpy as np
from pylab import *
import glob, os, sys
import scipy
import argparse
mingle_names={'W':'$W$','L':'$L$','LAMBDA':'$\Lambda$','GAMMA':'$\Gamma$','DELTA':'$\Delta$','X':'$X$','Z':'$Z$','W':'$W$','K':'$K$',
              'R': '$R$', 'S':'$S$', 'T':'$T$', 'U':'$U$', 'Y':'$Y$'}

if __name__=='__main__':
    
    aspect=0.8
    _col_ = 'k'
    
    Akom = load('Akom.npy')
    rest = load('rest.npy')
    (vmin,vmax,xmin,xmax,ymin,ymax)=rest

    #plot(0.5*40*Akom[:,:,1].ravel()/Akom[:,:,3].ravel())
    #plot(12*Akom[:,:,0].ravel()/Akom[:,:,3].ravel())
    #plot(0.5*16*Akom[:,:,2].ravel()/Akom[:,:,3].ravel())
    #show()
    #sys.exit(0)
    
    Ak2 = zeros(shape(Akom))
    #Ak2[:,:,3] = Akom[:,:,3]/15.
    #Ak2[:,:,0] = Akom[:,:,0]*12*4/22.
    #Ak2[:,:,2] = Akom[:,:,2]*16*1.5/22.
    #Ak2[:,:,1] = Akom[:,:,1]*40/22.

    Ak2[:,:,3] = Akom[:,:,3]/8.
    Ak2[:,:,0] = Akom[:,:,0]*1.5
    Ak2[:,:,2] = Akom[:,:,2]*1.5
    Ak2[:,:,1] = Akom[:,:,1]*1.5

    

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
    for i in range(len(wkpointi)):
        plot([wkpointi[i],wkpointi[i]], [ymin,ymax], _col_+'-', lw=1, zorder=1)
        
    plot([xmin,xmax],[0,0], _col_+'-',zorder=1)

    dytck=0.005
    Ntck=5
    for j in range(len(wkpointi)-1):
        for ix in range(1,Ntck):
            x = wkpointi[j]+(wkpointi[j+1]-wkpointi[j])*ix/float(Ntck)
            plot([x,x],[-dytck,dytck], _col_+'-')
        
    axis([xmin,xmax,ymin,ymax])
    xticks( wkpointi, wkpoints, fontsize='x-large' )

    imshow(Ak2, interpolation='bilinear', origin='lower', vmin=vmin, vmax=vmax, extent=[xmin,xmax,ymin,ymax], aspect=(xmax-xmin)*aspect/(ymax-ymin), zorder=2 )
    
    show()
    
