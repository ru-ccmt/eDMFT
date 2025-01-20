#!/usr/bin/env python
""" This script takes self-energy of DFT+DMFT code,
and performs max-entropy on auxiliary Green's function of the form
Gc(iom) = 1/(iom-Sigma+s_oo)
Then it performs Kramars-Kronig on resulting Dos, to obtain auxiliary G(omega).
From G(omega), one can obtain Sigma(omega) by
Sigma(omega) = omega+s_oo-1/G(omega)
which is written to 'Sig.out'
"""
import sys, shutil, os
os.environ["OMP_NUM_THREADS"] = "1"
import math, cmath
import numpy as np
from maxentropy import *


if len(sys.argv)<2:
    print('give input file Sigma(iom)')
    sys.exit(0)

Parallel, mpi_rank = False, 0
if len(sys.argv)>=2:
    Parallel = True
    from mpi4py import MPI
    # Get the MPI communicator
    comm = MPI.COMM_WORLD
    # Determine the rank and size
    mpi_rank = comm.Get_rank()
    mpi_size = comm.Get_size()
    # Each process prints a message
    print(f"Running in parallel mode rank={mpi_rank} size={mpi_size}")
    
Sfile = sys.argv[1]

with open(Sfile, 'r') as fi:
    firstlines = [next(fi),next(fi)]

Sdata = np.loadtxt(Sfile).transpose()
# number of baths
# number of baths
nb2=0; nz2=0;
for i in range(1,len(Sdata)):
    if sum(abs(Sdata[i]))>0:
        nb2 +=1
    else:
        nz2 +=1
nb = int(nb2/2)
nz = int(nz2/2)

s_oo_present=True
if not firstlines[0].lstrip().startswith('#') or not firstlines[1].lstrip().startswith('#'):
    s_oo_present=False
    s_oo = [Sdata[1+2*b,-1] for b in range(nb)]
    for b in range(nb):
        Sdata[1+2*b,:] -= s_oo[b]
    

# path to skrams executable
#exepath = os.environ.get('WIEN_DMFT_ROOT')
# need 'maxent_params.dat'
if not os.path.exists('maxent_params.dat'):
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
    with open('maxent_params.dat', 'w') as fo:
        fo.write(mparams)

    
exec(compile(open('maxent_params.dat', "rb").read(), 'maxent_params.dat', 'exec'))

iom = Sdata[0]
beta = math.pi/iom[0]
tau = np.linspace(0,beta,params['Ntau']+1)

print('nb=',nb, 'nz=',nz, 'beta=', beta, 's_oo_present=', s_oo_present)

if Parallel:
    per_proc = int(nb/mpi_size+0.99999)
    task = list(range(per_proc*mpi_rank,min(nb,per_proc*(mpi_rank+1))))
else:
    task = range(nb)
Sigt=[]
for b in task:
    Gm = 1/(iom*1j-Sdata[1+2*b]-Sdata[2+2*b]*1j)
    
    Gt = InverseFourier(Gm, iom, tau, beta, params['Nf'])
    np.savetxt('gt0.'+str(b), np.vstack((tau,Gt)).transpose())
    
    (Aw, omega) = MaximumEntropy(params, tau, Gt, sb='.'+str(b))
    
    # removes zero from mesh, because Kramars-Kronig cannot be done with zero
    izero = np.argmin(np.abs(omega))
    if np.abs(omega[izero])<1e-6:
        #omega_n = np.hstack([omega[:izero],omega[izero+1:]])
        #Aw_n = np.hstack([Aw[:izero],Aw[izero+1:]])
        omega_n = np.delete(omega, izero)
        Aw_n = np.delete(Aw, izero)
    else:
        omega_n = omega
        Aw_n = Aw

    Aw_r = KramarsKronig(omega_n, Aw_n)
    Gd = -pi*(Aw_r + Aw_n*1j)
    Sc = omega_n-1/Gd
    
    np.savetxt('sig.'+str(b), np.array([omega_n,np.real(Sc),np.imag(Sc)]).transpose())
    if b==0: Sigt.append(omega_n)
    Sigt.append(np.real(Sc))
    Sigt.append(np.imag(Sc))

if Parallel:
    Sigt_all = comm.gather(Sigt, root=0)
    if mpi_rank==0:
        Sigt = [item for sublist in Sigt_all for item in sublist]
    
if mpi_rank==0:
    for z in range(nz):
        Sigt.append( np.zeros(len(Sigt[0])) )
        Sigt.append( np.zeros(len(Sigt[0])) )
    
    Sigt = np.array(Sigt).T
    with open('Sig.out', 'w') as fo:
        if s_oo_present:
            print(firstlines[0].strip(), file=fo)
            print(firstlines[1].strip(), file=fo)
        else:
            for b in range(nb):
                Sigt[:,1+2*b] += s_oo[b]
        
        for sg in Sigt:
            for b in sg:
                print(b, end=' ', file=fo)
            print(file=fo)
