#!/usr/bin/env python
import numpy as np
import os, glob, shutil, sys, re
import subprocess
import argparse
import utils
from saverage import saverage_sig

def ChangeMatsubaraFlag(fname, matsubara, x_range=None, nfrq=None):
    with open(fname, 'r') as fi:
        lines = fi.readlines()
    if x_range is None and nfrq is None:
        lines[1] = re.sub(r"^\s*\d+", str(matsubara), lines[1])
    else:
        dat = lines[1].split()
        dat[0] = str(matsubara)
        if x_range is not None:
            dat[4],dat[5] = str(x_range[0]),str(x_range[1])
        if nfrq is not None:
            dat[3] = str(nfrq)
        lines[1] = ' '.join(dat)
    # Write the updated content back to the file
    with open(fname, 'w') as fi:
        fi.writelines(lines)
    
if __name__ == '__main__':
    usage = """creates onreal subdirectory, copies imaginary axis calculation to the new directory
    runs maxentropy to produce self-energy on the real axis, and 
    runs lapw0,lapw1,dmft1 to get density of states
    run dmftp in case.klist_band exists.
    """
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('fname', nargs='?', default='gc', type=str, help='filename either gc or sig.inp. Default=gc')
    parser.add_argument('-y', type=str, default=None, help='frequency range (in eV) of the spectral plot in A(k,omega) -x0:10')
    parser.add_argument('-a', default='3', type=str, help='how many self-energy\'s of the last few steps shoud be averaged over for maxent.')
    parser.add_argument('-o', type=str, default='onreal', help='directory name for real axis calculation. Default: onreal')
    parser.add_argument('-m', type=str, default='maxent', help='directory name for maximum entropy calculation. Default: maxent')
    parser.add_argument('-n', default=None, type=int, help='how many point for A(k,omega) plot in frequency omega.')
    args = parser.parse_args()
    maxent_dir = args.m
    onreal_dir = args.o
    
    xrng = None
    if args.y is not None:
        w = args.y.split(':')
        xrng = [float(w[i]) if w[i]!='' else None for i in range(2)]
        if xrng[0]==None: xrng[0]=-3.0
        if xrng[1]==None: xrng[1]= 1.0
    lastn = int(args.a)
    
    w2k = utils.W2kEnvironment() 
    ChangeMatsubaraFlag(onreal_dir+'/'+w2k.case+'.indmfl', 0, xrng, args.n)
    
    fsigs = glob.glob('sig.inp.*.*')
    itrs = [list(map(int,fg.split('.')[-2::])) for fg in fsigs]
    # sorts them in descending order
    itrs = sorted(itrs, key=lambda x:-x[0]-x[1]/1000.)
    lastn = min(lastn, int(len(itrs)*0.1)+1) # should never be more than 10%
    
    ss = ['sig.inp.'+str(it[0])+'.'+str(it[1]) for it in itrs[:lastn]]
    print('averaging over last few sigmas: ', ss)
    
    saverage_sig(ss, 'sig.inpx')
    if not os.path.exists(maxent_dir):
        os.mkdir(maxent_dir)
    
    shutil.move('sig.inpx', maxent_dir+'/sig.inpx')
    if os.path.exists('maxent_params.dat'):
        shutil.move('maxent_params.dat', maxent_dir+'/maxent_params.dat')
    else:
        mparams="""params={'statistics': 'fermi', # fermi/bose
        'Ntau'      : 400,     # Number of time points
        'L'         : 30.0,    # cutoff frequency on real axis
        'x0'        : 0.01,   # low energy cut-off
        'bwdth'     : 0.004,    # smoothing width
        'Nw'        : 450,     # number of frequency points on real axis
        'gwidth'    : 2*15.0,  # width of gaussian
        'idg'       : 1,       # error scheme: idg=1 -> sigma=deltag ; idg=0 -> sigma=deltag*G(tau)
        'deltag'    : 0.01,   # error
        'Asteps'    : 4000,    # anealing steps
        'alpha0'    : 1000,    # starting alpha
        'min_ratio' : 0.001,    # condition to finish, what should be the ratio
        'iflat'     : 1,       # iflat=0 : constant model, iflat=1 : gaussian of width gwidth, iflat=2 : input using file model.dat
        'Nitt'      : 500,     # maximum number of outside iterations, 1000 can take too much time
        'Nr'        : 0,       # number of smoothing runs
        'Nf'        : 40,      # to perform inverse Fourier, high frequency limit is computed from the last Nf points
        }"""
        with open(maxent_dir+'/maxent_params.dat', 'w') as fo:
            fo.write(mparams)



    dmfe = utils.DmftEnvironment()
    cmd = 'cd '+maxent_dir+'; '+dmfe.ROOT+'/maxent_run.py sig.inpx'
    print('executing...', cmd)
    info=subprocess.call(cmd,shell=True,stdout=sys.stdout,stderr=sys.stderr)
    
    if not os.path.exists(onreal_dir):
        os.mkdir(onreal_dir)
    
    cmd = 'cd '+onreal_dir+'; '+dmfe.ROOT+'/dmft_copy.py ../'
    print('executing...', cmd)
    info=subprocess.call(cmd,shell=True,stdout=sys.stdout,stderr=sys.stderr)
    shutil.copyfile(maxent_dir+'/Sig.out', onreal_dir+'/sig.inp')
    shutil.copyfile('mpi_prefix.dat', onreal_dir+'/mpi_prefix.dat')
    fname = w2k.case+'.klist_band'
    if os.path.exists(fname):
        shutil.copyfile(fname, onreal_dir+'/'+fname)
    
    w2k = utils.W2kEnvironment() 
    ChangeMatsubaraFlag(onreal_dir+'/'+w2k.case+'.indmfl', 0, xrng, args.n)
    cmd = 'cd '+onreal_dir+'; '+w2k.WIENROOT+'/x_lapw -f '+w2k.case+' lapw0'
    print('executing...', cmd)
    info=subprocess.call(cmd,shell=True,stdout=sys.stdout,stderr=sys.stderr)
    cmd = 'cd '+onreal_dir+'; '+dmfe.ROOT+'/x_dmft.py lapw1'
    print('executing...', cmd)
    info=subprocess.call(cmd,shell=True,stdout=sys.stdout,stderr=sys.stderr)
    cmd = 'cd '+onreal_dir+'; '+dmfe.ROOT+'/x_dmft.py dmft1'
    print('executing...', cmd)
    info=subprocess.call(cmd,shell=True,stdout=sys.stdout,stderr=sys.stderr)

    if os.path.exists(w2k.case+'.klist_band'):
        shutil.copyfile(w2k.case+'.klist_band', onreal_dir+'/'+w2k.case+'.klist_band')
        cmd = 'cd '+onreal_dir+'; '+dmfe.ROOT+'/x_dmft.py lapw1 --band'
        print('executing...', cmd)
        info=subprocess.call(cmd,shell=True,stdout=sys.stdout,stderr=sys.stderr)
        cmd = 'cd '+onreal_dir+'; '+dmfe.ROOT+'/x_dmft.py dmftp'
        print('executing...', cmd)
        info=subprocess.call(cmd,shell=True,stdout=sys.stdout,stderr=sys.stderr)
    
