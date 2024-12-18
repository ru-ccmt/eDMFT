#!/usr/bin/env python

from scipy import *
from pylab import *
import re, os, sys
import glob
import argparse

Ry2eV = 13.60569193

if __name__ == '__main__':
    usage = 'Plots bands by reading  case.output1 or case.outdmftp files'
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('fname', nargs='?', default=None, type=str, help='filename either case.output1 or case.outdmftp')
    parser.add_argument('-x', type=str, help='x range of the plot')
    parser.add_argument('-y', type=str, help='y range of the plot')
    parser.add_argument('-g', default=False, action='store_true', help='grid')
    parser.add_argument('-ef', type=float, default=None, help='Fermi energy in eV')
    parser.add_argument('-Ry', default=False, action='store_true', help='display in Rydberg units')
    args = parser.parse_args()

    if (args.fname==None):
        print('ERROR : Give file case.output1 or case.outdmftp file')
        sys.exit(1)
    if (args.ef==None):
        if os.path.exists('EF.dat'):
            EF = loadtxt('EF.dat')
        else:
            scf = glob.glob('*.scf')[0]
            with open(scf, 'r') as fscf:
                lines = fscf.readlines()
            for i in range(len(lines)-1,0,-1):
                if lines[i][:4]==':FER':
                    EF = float(lines[i].split()[-1])
                    EF *= Ry2eV
                    break
            print('EF=', EF)
    else:
        EF = args.ef

    fi = open(args.fname,'r')
    k_str = '     K='
    line = next(fi)
    first_band = []
    kpt = []
    kene = []
    while(fi):
        if (line[:7]==k_str):
            spl = line.split()
            vk = list(map(float,spl[1:4]))
            kname = ''
            if len(spl)>4: kname = spl[4]
            #print vk,kname
            kpt.append( [vk,kname] )
            next(fi)
            next(fi)
            ene=[]
            for i in range(100):
                line = next(fi)
                m = re.search('EIGENVALUES', line)
                if m is not None:
                    sp = line.split()
                    i_first = int(sp[0])
                    break
                spl = line.split()
                if len(spl)>0:
                    ene += list(map(float,spl))
            kene.append( ene )
            first_band.append(i_first)
            
        try:    
            line = next(fi)
        except StopIteration:
            break
    
    nkp = len(kene)
    i_first = max(first_band)
    last = [first_band[ik]+len(kene[ik]) for ik in range(nkp)]
    i_last = min(last)
    nbands = i_last-i_first

    enek = zeros((nkp,nbands))
    for ik in range(nkp):
        i_start = i_first-first_band[ik]
        i_end = i_start + nbands
        if args.Ry:
            enek[ik,:] = array(kene[ik][i_start:i_end])
        else:
            enek[ik,:] = array(kene[ik][i_start:i_end])*Ry2eV-EF
    xtcks=[]
    xlabels=[]
    for ik in range(nkp):
        if (kpt[ik][1]):
            xtcks.append(ik)
            xlabels.append('$'+kpt[ik][1]+'$')

    ax = subplot()
    for ib in range(nbands):
        plot(enek[:,ib], '-')
    
    xlim([0,nkp])
    
    x_ticks = True
    if args.x!=None:
        w = args.x.split(':')
        if (w[0]=='' and w[1]==''): # occurs if given -x:
            x_ticks = False
        else:
            xlim( list(map(float,w)) )
        
    if args.y!=None:
        y_lim = list(map(float,args.y.split(':')))
        if not args.Ry:
            ylim( y_lim )
        else:
            ylim( (y_lim[0]+EF)/Ry2eV, (y_lim[1]+EF)/Ry2eV)
    
    if x_ticks:
        xticks(xtcks, xlabels)
    else:
        bottom, top = ylim()
        for i in range(len(xtcks)):
            plot( [xtcks[i],xtcks[i]], [bottom,top], 'k:' )
            text(xtcks[i], bottom-(top-bottom)*0.1, xlabels[i], va='bottom', ha='center')
        ax.tick_params(bottom=False, top=True)
        ax.tick_params(labelbottom=False, labeltop=True)

    if args.Ry:
        left,right = xlim()
        plot([left,right], [EF/Ry2eV,EF/Ry2eV], 'r:')
    
    ylabel('Energy[eV]')
    if args.g!=None:
        grid()
    show()
    
        
