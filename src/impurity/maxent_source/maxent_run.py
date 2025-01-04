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
from scipy import *
import shutil
from maxentropy import *

if len(sys.argv)<2:
    print('give input file Sigma(iom)')
    sys.exit(0)

Sfile = sys.argv[1]

with open(Sfile, 'r') as fi:
    firstlines = [next(fi),next(fi)]

Sdata = loadtxt(Sfile).transpose()
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
exepath = os.environ.get('WIEN_DMFT_ROOT')
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
beta = pi/iom[0]
tau = linspace(0,beta,params['Ntau']+1)

print('nb=',nb, 'nz=',nz, 'beta=', beta, 's_oo_present=', s_oo_present)
Sigt=[]
for b in range(nb):
    Gm = 1/(iom*1j-Sdata[1+2*b]-Sdata[2+2*b]*1j)
    
    Gt = InverseFourier(Gm, iom, tau, beta, params['Nf'])
    savetxt('gt0.'+str(b), vstack((tau,Gt)).transpose())
    
    (Aw, omega) = MaximumEntropy(params, tau, Gt)
    savetxt('dos.out.'+str(b), vstack((omega,Aw)).transpose())

    shutil.copy2('gtn', 'gtn.'+str(b))
    
    # removes zero from mesh, because Kramars-Kronig can be be done with zero
    #izero = omega.tolist().index(0)
    izero = argmin(abs(omega))
    if abs(omega[izero])<1e-6:
        omega_n = hstack([omega[:izero],omega[izero+1:]])
        Aw_n = hstack([Aw[:izero],Aw[izero+1:]])
    else:
        omega_b = omega
        Aw_n = Aw
    savetxt('dosn', vstack((omega_n,Aw_n)).transpose())
    
    # Performs Kramars-Kronig
    cmd = exepath + '/skrams -cn 2 -s -pi dosn > Gc'
    print(cmd)
    print(os.popen(cmd).read())
    Gdata = loadtxt('Gc').transpose()
    om = Gdata[0]
    Sc = om-1/(Gdata[1]+Gdata[2]*1j)
    savetxt('sig.'+str(b), array([om,real(Sc),imag(Sc)]).transpose())
    if b==0: Sigt.append(om)
    Sigt.append(real(Sc))
    Sigt.append(imag(Sc))

for z in range(nz):
    Sigt.append( zeros(len(Sigt[0])) )
    Sigt.append( zeros(len(Sigt[0])) )

Sigt = array(Sigt).transpose()
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
