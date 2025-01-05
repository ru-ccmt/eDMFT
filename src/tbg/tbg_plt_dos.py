#!/usr/bin/env python
from pylab import *
import os
import numpy as np
import argparse
from matplotlib import colors as mcolors

colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
cname = ['red','black','silver']

if __name__=='__main__':
    if os.path.exists('pDOS.out'):
        dat = np.loadtxt('pDOS.out').T
    elif os.path.exists('pDOS.dat'):
        dat = np.loadtxt('pDOS.dat').T
    else:
        print('No pDOS.out|dat exists')
        
    unit=50
    usage = 'Plots partial DOS for TBG'
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('-x', type=str, help='x range of the plot in the form -x0:10')
    parser.add_argument('-y', type=str, help='y range of the plot in the form -y0:10')
    parser.add_argument('-g', default=False, action='store_true', help='grid')
    args = parser.parse_args()

    xrng = [-unit,unit]
    yrng = [0,max(dat[1])*1.1/unit]
    
    if args.x!=None:
        w = args.x.split(':')
        xrng = [float(w[i]) if w[i]!='' else None for i in range(2)]
    if args.y!=None:
        w = args.y.split(':')
        yrng = [float(w[i]) if w[i]!='' else None for i in range(2)]
    
    ax = subplot(111)
    
    plot(dat[0]*unit, dat[1]/unit, color=colors[cname[1]], lw=2, label='total')
    plot(dat[0]*unit, dat[3]/unit, color=colors[cname[0]], lw=2, label='AA')
    plot(dat[0]*unit, dat[4]/unit, color=colors[cname[2]], lw=2, label='BA')
    xlim(xrng)
    ylim(yrng)
    #plot([0,0],[0,ymax], ':', color=colors['silver'], lw=1)
    #ax.yaxis.tick_right()
    #ax.xaxis.tick_top()
    legend(loc='best')
    ax = gca()
    ax.set_xlim(ax.get_xlim()[::-1])
    ylabel('partial DOS')
    xlabel('Energy[meV]')
    if args.g: ax.grid()
    show()
