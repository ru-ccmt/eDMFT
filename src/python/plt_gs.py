#!/usr/bin/env python
from pylab import *
import numpy as np
import glob
import utils
from indmffile import Indmfl, ParsIndmfi
from functools import reduce
import re, os
import argparse
from imp2lattc import ImpurityLatticeConnection, SimplifySiginds

if __name__=='__main__':
    usage = 'Plots Green\'s function or self-energy by reading imp.?/Gf.out.* and case.gc? files'
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('fname', nargs='?', default='gc', type=str, help='filename either gc or sig.inp. Default=gc')
    parser.add_argument('-x', type=str, help='x range of the plot in the form -x0:10')
    parser.add_argument('-y', type=str, help='y range of the plot in the form -y0:10')
    parser.add_argument('-i', type=str, default=None, help='impurity icols to show in the form of dic like {0:[1,2,3]}')
    parser.add_argument('-l', type=str, default=None, help='lattice lcols to show in the form of dic like {0:[1,2,3], 1:[1,2,3]}')
    parser.add_argument('-n', type=str, default=None, help='lattice lcols_dn to show in the form of dic like {0:[1,2,3], 1:[1,2,3]}')
    parser.add_argument('-g', default=False, action='store_true', help='grid')
    parser.add_argument('-a', default='1', type=str, help='how many impurity steps displayed. Default one.')
    args = parser.parse_args()
    
    xrng = [None,None]
    yrng = [None,None]
    if args.x!=None:
        w = args.x.split(':')
        xrng = [float(w[i]) if w[i]!='' else None for i in range(2)]
    if args.y!=None:
        w = args.y.split(':')
        yrng = [float(w[i]) if w[i]!='' else None for i in range(2)]
    lastn = int(args.a)
        
    #print('fname=', args.fname)
    #print('xrng=', xrng)
    #print('yrng=', yrng)
    #print('lastn=', lastn)

    env = utils.W2kEnvironment()
    if args.fname=='gc':
        fnames = [env.case+'.gc', 'Gf.out', 'G']
    else:
        fnames = ['sig.inp', 'Sig.out','$\Sigma$']
        
    inl = Indmfl(env.case)
    inl.read()
    iSiginds = ParsIndmfi(env.case)
    _icols = SimplifySiginds(iSiginds)
    _lcols = SimplifySiginds(inl.siginds)
    _lcols_dn=None
    if os.path.isfile(env.case+'.indmfldn'):
        inldn = Indmfl(env.case, 'indmfldn')
        inldn.read()
        _lcols_dn = SimplifySiginds(inl.siginds)
    imp2latt = ImpurityLatticeConnection(_lcols, _icols, sys.stdout)
    # these are columns we can fill in with our impurity problems
    print('_icols=', _icols)
    print('_lcols=',_lcols)
    
    # First create simple columns for communicating with the user.
    icols={}
    for ii in _icols:
        icols[ii+1] = list(range(len(_icols[ii])))
    lcols={}
    for icix in _lcols:
        lcols[icix] = (array(_lcols[icix])-min(_lcols[icix])).tolist()
    lcols_dn=None
    if _lcols_dn is not None:
        lcols_dn={}
        for icix in _lcols_dn:
            lcols_dn[icix] = (array(_lcols_dn[icix])-min(_lcols_dn[icix])).tolist()
    
    print('imp2lattice connection=', imp2latt)
    for ii in imp2latt:
        print(' imp.'+str(ii)+'/Gf.out icols['+str(ii+1)+']=', icols[ii+1])
        print(' connected to ')
        for icix in imp2latt[ii]:
            print('  gc'+str(icix)+'   lcols['+str(icix)+']=',array(lcols[icix]))
            if lcols_dn is not None:
                print('  gc'+str(icix)+'dn lcols['+str(icix)+']=',array(lcols_dn[icix]))
    print('All available columns:')
    print(' icols=', icols)
    print(' lcols=', lcols)
    if lcols_dn is not None:
        print(' lcols_dn=', lcols_dn)
    changed=False
    if args.i is not None:
        icols = eval(args.i)
        changed = True
    if args.l is not None:
        lcols = eval(args.l)
        changed = True
    if lcols_dn is not None and args.n is not None:
        lcols_dn = eval(args.n)
        changed = True
    if changed:
        print('user changed to')
        print(' icols=', icols)
        print(' lcols=', lcols)
        if lcols_dn is not None:
            print(' lcols_dn=', lcols_dn)
        
    fig, axs = plt.subplots(2, len(imp2latt), sharex=True)
    if len(imp2latt) == 1:
        axs = axs.reshape(2, 1)
    # Below plotting impurity quantity
    for ii in imp2latt:
        il=0
        if ii+1 in icols and os.path.exists('imp.'+str(ii)):
            #for icix in icols:
            icix = ii+1
            dr = 'imp.'+str(ii)
            # finds all imp.?/Gf.out.?.? for all iterations
            fgimp = glob.glob(dr+'/'+fnames[1]+'.*')
            # finds available iterations, i.e, last two numbers
            itrs = [list(map(int,fg.split('.')[-2::])) for fg in fgimp]
            # sorts them in descending order
            itrs = sorted(itrs, key=lambda x:-x[0]-x[1]/1000.)
            its = itrs[:lastn]
            # repeating for the last few iterations, if needed
            for it in its:
                fname = dr+'/'+fnames[1]+'.'+str(it[0])+'.'+str(it[1])
                print('reading', fname, 'cols=', icols[icix])
                dat = loadtxt(fname).T
                om = dat[0]
                Gimp = dat[1::2]+dat[2::2]*1j
                for i,j in enumerate(icols[icix]):
                    axs[0][ii].plot(om,Gimp[j].imag, 'C'+str(il%10), label='imp['+str(it[0])+','+str(j)+']')
                    axs[1][ii].plot(om,Gimp[j].real, 'C'+str(il%10), label='imp['+str(it[0])+','+str(j)+']')
                    il+=1
        
        # Below plotting lattice quantity
        il=0
        foundOneYet,foundOneYetdn=False,False
        for icix in imp2latt[ii]: #lcols:
            if icix in lcols and not foundOneYet:
                fname = fnames[0]+str(icix)
                if not os.path.exists(fname):
                    continue
                else:
                    foundOneYet=True
                cls = array(lcols[icix])
                print('reading', fname, 'cols=', cls)
                dat = loadtxt(fname).T
                w = dat[0]
                Glat = array([dat[1+2*j]+dat[2+2*j]*1j for j in cls])
                for i,j in enumerate(cls):
                    axs[0][ii].plot(w,Glat[i].imag, 'C'+str(il%10)+'.', label='lat['+str(icix)+',$'+inl.legends[icix][j]+'$]')
                    axs[1][ii].plot(w,Glat[i].real, 'C'+str(il%10)+'.', label='lat['+str(icix)+',$'+inl.legends[icix][j]+'$]')
                    il+=1
                
            if lcols_dn is not None and icix in lcols_dn and foundOneYetdn:
                    fname = fnames[0]+str(icix)+'dn'
                    if not os.path.exists(fname):
                        continue
                    else:
                        foundOneYetdn=True
                    cls = array(lcols_dn[icix])
                    print('reading', fname, 'cols=', cls)
                    dat = loadtxt(fname).T
                    w = dat[0]
                    Glat = array([dat[1+2*j]+dat[2+2*j]*1j for j in cls])
                    for i,j in enumerate(cls):
                        axs[0][ii].plot(w,Glat[i].imag, 'C'+str(il%10)+'.', label='lat['+str(icix)+'dn,$'+inl.legends[icix][j]+'$]')
                        axs[1][ii].plot(w,Glat[i].real, 'C'+str(il%10)+'.', label='lat['+str(icix)+'dn,$'+inl.legends[icix][j]+'$]')
                        il+=1
                
        if xrng[0]!=None and xrng[1]!=None:
            axs[0][ii].set_xlim(xrng)
        elif xrng[0]!=None and xrng[1]==None:
            axs[0][ii].set_xlim([xrng[0],om[-1]])
        elif xrng[1]!=None and xrng[0]==None:
            axs[0][ii].set_xlim([om[0],xrng[1]])
        
        if yrng[0]!=None and yrng[1]!=None:
            axs[0][ii].set_ylim(yrng)
            axs[1][ii].set_ylim(yrng)
        elif yrng[0]!=None and yrng[1]==None:
            axs[0][ii].set_ylim(bottom=yrng[0])
            axs[1][ii].set_ylim(bottom=yrng[0])
        elif yrng[0]==None and yrng[1]!=None:
            axs[0][ii].set_ylim(top=yrng[1])
            axs[1][ii].set_ylim(top=yrng[1])
        
        axs[0][ii].legend(loc='best', fontsize='small')
        #axs[1].legend(loc='best', fontsize='small')
        axs[1][ii].set_xlabel('Energy[eV]')
        axs[0][ii].set_title('imp.'+str(ii))
    axs[1][0].set_ylabel('Re '+fnames[2])
    axs[0][0].set_ylabel('Im '+fnames[2])
    
    if args.g:
        for ii in imp2latt:
            axs[0][ii].grid()
            axs[1][ii].grid()
    show()
