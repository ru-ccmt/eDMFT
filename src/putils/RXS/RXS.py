#!/usr/bin/env python

from scipy import *
import matplotlib.pyplot as plt
import struct1
import utils
import sys, os
from scipy import interpolate
import subprocess
import re
from w2k_atpar import readpotential, readlinearizatione, atpar, rint13, rint13g
import trafoso
from thompson import *
from scattering import Scattering
import shutil

#def ScatteringVector(theta, lmbda, abc, hkl):
#    a,b,c = abc
#    h,k,l = hkl
#    q = 2*pi*array([h/a, k/b, l/c])
#
#    # solving the equation k*q = q**2/2
#    # with k= 2*pi/lmbda * [sin(t)*cos(fi),sin(t)*sin(fi),cos(t)]
#    xi=0.5*lmbda*((h/a)**2+(k/b)**2+(l/c)**2)
#
#    if h==0 and k==0:
#        if (cos(theta) - lmbda/(2*abc[2]))>1e-4:
#            if Sprint: print 'If you scatter in z direction, there is only one solution, namely cos(theta)=lmbda/(2*c), and you did not satisfy it'
#            return [[],[],[],[]]
#        else:
#            if Sprint: print 'phi is still undertermined'
#            return [[],[],[],[]]
#        
#    if (theta==0):
#        if ( abs(xi*c/l - 1)>1e-6):
#            if Sprint: print 'No solution 1'
#            return [[],[],[],[]]
#        else:
#            ek  = array([0,0,1.])
#            ekp = array([-h/a,-k/b,1/lmbda-l/c])
#            ekp *= 1/sqrt(sum(kout**2))
#            sig = cross(ek,ekp)
#            sig *= 1.0/sqrt(sum(sig**2))
#            p_in = cross(sig,ek)
#            p_in *= 1.0/sqrt(sum(p_in**2))
#            p_out = cross(sig,ekp)
#            p_out *= 1.0/sqrt(sum(p_out**2))
#            return [sig,p_in,p_out]
#        
#    eta=(xi-l/c*cos(theta))/sin(theta)
#    ax = (h/a)**2 + (k/b)**2
#    
#    #  h/a*sqrt(1-x^2) + k/b*x = eta
#    #
#    #  eta^2*(k/b)^2 > ( (h/a)^2 + (k/b)^2 ) * (eta^2-(h/a)^2)
#    #  ( (h/a)^2 + (k/b)^2 + (l/c)^2 )*sin(theta)^2  >  xi^2-2*xi*l/c*cos(theta) + (l/c)^2
#    # four possible solutions, but not all are good
#    if h!=0:
#        bx = eta*k/b
#        cx = eta**2-(h/a)**2
#        sinf1 = (bx + sqrt(bx**2-ax*cx))/ax
#        sinf2 = (bx - sqrt(bx**2-ax*cx))/ax
#        sinfi = [ sinf1, sinf1, sinf2, sinf2]
#        cosfi = [ sqrt(1-sinf1**2),-sqrt(1-sinf1**2), sqrt(1-sinf2**2),-sqrt(1-sinf2**2)]
#    else:
#        bx = eta*h/a
#        cx = eta**2-(k/b)**2
#        cosf1 = (bx + sqrt(bx**2-ax*cx))/ax
#        cosf2 = (bx - sqrt(bx**2-ax*cx))/ax
#        cosfi = [ cosf1, cosf1, cosf2, cosf2]
#        sinfi = [ sqrt(1-cosf1**2),-sqrt(1-cosf1**2), sqrt(1-cosf2**2),-sqrt(1-cosf2**2)]
#        
#    if eta > max(abs(h/a),abs(k/b)) or eta < -max(abs(h/a),abs(k/b)):
#        if Sprint : print 'No solution 2', eta, h/a, k/b
#        return [[],[],[],[]]
#
#    #sinf1 = (bx + sqrt(bx**2-ax*cx))/ax
#    #sinf2 = (bx - sqrt(bx**2-ax*cx))/ax
#    #sinfi = [ sinf1, sinf1, sinf2, sinf2]
#    #cosfi = [ sqrt(1-sinf1**2),-sqrt(1-sinf1**2), sqrt(1-sinf2**2),-sqrt(1-sinf2**2)]
#    # is this solution good?
#    good = []
#    for i in range(4):
#        if ( abs(h/a*cosfi[i] + k/b*sinfi[i]-eta)<1e-6): good.append(i)
#    
#    pi_in=[]
#    pi_out=[]
#    sigm=[]
#    eks=[]
#    for i in good:
#        ek=array([sin(theta)*cosfi[i],sin(theta)*sinfi[i],cos(theta)])
#        ekp = (2*pi/lmbda*ek - q)/(2*pi/lmbda)
#
#        sig = cross(ek,ekp)
#        sig *= 1.0/sqrt(sum(sig**2))
#
#        p_in = cross(sig,ek)
#        p_in *= 1.0/sqrt(sum(p_in**2))
#        
#        p_out = cross(sig,ekp)
#        p_out *= 1.0/sqrt(sum(p_out**2))
#        pi_in.append(p_in)
#        pi_out.append(p_out)
#        sigm.append(sig)
#        eks.append(ek)
#        
#        if Sprint: 
#            print 'k_in = ', ek*2*pi/lmbda
#            print 'k_out= ', ekp*2*pi/lmbda
#            print 'k_in-k_out=', (ek-ekp)*2*pi/lmbda
#            print 'q=', 2*pi*array([h/a,k/b,l/c])
#            print 'norm(ek)=', sum(ek**2)
#            print 'norm(ekp)=', sum(ekp**2)
#            print 'sig=', sig, sum(sig**2)
#            print 'pi_in=', p_in, sum(p_in**2)
#            print 'pi_out=', p_out, sum(p_out**2)
#    return (sigm,pi_in,pi_out,eks)
#
#
#
#def Azimuthal_single_valued(theta, lmbda, hkl, abc):
#    """
#    To compute the azimuthal angle, we compute
#         cz = dot(ek,eq) , which is like projection to the new z-axis  (cos(t))
#         sz_cp = dot(ek,psi_is_zero) , which is like projection to the x axis (sin(t)cos(psi))
#         cpsi = sz_cp/sqrt(1-cz**2) , which is like  cos(psi)
#         psi = arccos(cpsi)   is the final result
#         
#      But we need to figure out where is psi vanishing, to make it single valued. We do that by requiring
#        ek * (eq x psi_iz_zero ) == 0
#      If h=0 or k=0, we choose psi_is_zero in x or y direction, respectively.
#      In both cases the equation
#           ek * (eq x psi_iz_zero ) == 0
#      is satisfied when
#           cos(theta) = lmbda*l/(2*c)
#    If h!=0 and k!=0 then the equations get more involved and I did not yet worked them out
#    """
#    # psi will be counted from the axis which has no projection of q=(h,k,l)
#    psi_is_zero = zeros(3)
#    which_axis=1
#    if hkl[0]==0:
#        psi_is_zero[0]=1
#    elif hkl[1]==0:
#        psi_is_zero[1]=1
#    else:
#        print 'Not yet worked out!'
#        sys.exit(0)
#    #should_vanish = dot(ek,cross(eq,psi_is_zero))
#    theta0 = arccos(lmbda*hkl[2]/(2*abc[2]))
#    if theta<theta0:
#        return -1
#    else:
#        return 1
#

def ReadIndmfl(case):
    def divmodulo(x,n):
        "We want to take modulo and divide in fortran way, so that it is compatible with fortran code"
        return ( sign(x)* (abs(x)/n) , sign(x)*mod(abs(x),n))

    findmfl = open(case+'.indmfl', 'r')
    lines = [line.split('#')[0].strip() for line in findmfl.readlines()] # strip comments
    findmfl.close()
    lines = (line for line in lines if line)  # strip blank lines & create generator expression
    
    dat = lines.next().split()[:4]
    hybr_emin, hybr_emax, Qrenorm, projector = (float(dat[0]),float(dat[1]),int(dat[2]),int(dat[3]))
    dat = lines.next().split()[:6]
    matsubara, broadc, broadnc, om_npts, om_emin, om_emax = (int(dat[0]),float(dat[1]),float(dat[2]),int(dat[3]),float(dat[4]),float(dat[5]))
    natom = int(lines.next())

    Rmt2 = zeros(natom,dtype=float)
    (atms, Lsa, qsplita, icpsa) = ([],[],[],[])
    for ia in range(natom):
        dat = lines.next().split()
        iatom, nL, locrot_shift = map(int, dat[:3])
        atms.append(iatom)
        if len(dat)>3: Rmt2[ia] = float(dat[3])
            
        (shift,locrot) = divmodulo(locrot_shift,3)
        if locrot==-1:
            locrot=3
        elif locrot==-2:  # with locrot -2, we have different rotation matrix for each atom
            locrot=3*nL
                    
        (Ls, qsplit, icps) = (zeros(nL,dtype=int), zeros(nL,dtype=int), zeros(nL,dtype=int))
        for il in range(nL):
            (Ls[il], qsplit[il], icps[il]) = map(int, lines.next().split()[:3])
        
        Lsa.append( Ls )
        qsplita.append( qsplit )
        icpsa.append( icps )
        
        new_zx = [[float(x) for x in lines.next().split()[:3]] for loro in range(abs(locrot))]
        vec_shift = [float(x) for x in lines.next().split()] if shift else None

    
    Ncix, max_dim, max_components = map(int,lines.next().split()[:3])
    Siginds=[]
    T2Cs=[]
    for ici in range(Ncix):
        (icix,dim,ncomponents) = map(int,lines.next().split()[:3])
        legend = lines.next()
        #print 'legend=', legend
        Sigind = zeros((dim,dim),dtype=int)
        for j in range(dim):
            Sigind[j,:] = map(int,lines.next().split())
        T2C = zeros((dim,dim),dtype=complex)
        for j in range(dim):
            w = array(map(float,lines.next().split()[:2*dim]))
            T2C[j,:] = w[::2]+w[1::2]*1j
        T2Cs.append(T2C)
        Siginds.append(Sigind)
        #print 'Sigind=', Sigind
        #print 'T2c=', T2c
        
            
    print 'atms=', atms
    print 'Rmt2=', Rmt2
    print 'Lsa=', Lsa
    print 'qsplita=', qsplita
    print 'icpsa=', icpsa
    return (atms, Lsa, icpsa, Siginds, T2Cs)

def read_projector():
    fi = open('projectorw.dat', 'r')
    line = fi.next()
    (Nv, Nrmax) = map(int, line[1:].split()[:2])

    wave=[]
    Type=[]
    for ic in range(Nv):
        line = fi.next()
        (Nr,Nr0,jatom_1,lc) = map(int,line[1:].split()[:4])
        Ar = zeros(Nr)
        Br = zeros(Nr)
        Rx = zeros(Nr)
        for ir in range(Nr):
            (Rx[ir], Ar[ir], Br[ir]) = map(float, fi.next().split()[:3])
        wave.append( (Rx,Ar,Br) )
        Type.append( (jatom_1-1, lc) )
        
    return (wave, Type)

def SelectValence(atoms, valence_l, wave,Type, plot_waves=False):
    valence_waves=[]
    valence_jatom=[]
    for i,(jatom,l) in enumerate(Type):
        Rx, Ag, Bg = wave[i]
        
        if jatom in atoms and l==valence_l:
            print 'Found desired valence wave function', jatom, l
            valence_waves.append( (Ag,Bg/speed_of_light) )
            valence_jatom.append( jatom )
        
        Nr = struct.npt[jatom]
        r0 = struct.r0[jatom]
        Rmt = struct.rmt[jatom]
        dx = log(Rmt/r0)/(Nr-1)
        #Rx = r0 * exp( arange(Nr)*dx ) # radial points
        
        norm = rint13g(1.0,1.0/speed_of_light**2,Ag,Bg,Ag,Bg,dx,r0)
        print 'iat=', jatom, chemical_symbol[l], 'norm=', norm
        
        if plot_waves:
            plt.plot(Rx, Ag, label=str(jatom)+':'+str(l))
    
    if plot_waves:
        plt.title('Valence wave functions from projector')
        plt.legend(loc='best')
        plt.show()
        
    return (valence_waves, valence_jatom)

def read_core(case, struct):
    # reading case.inc
    fc = open(case+'.inc', 'r')
    n_kappa_occup =[]
    iprint=[]
    for jatom in range(struct.nat):
        line = fc.next()
        data = line.split()
        nc = int(data[0])
        iprint.append( int(data[2]) )
        #nc = int(line[:3])
        nko = zeros((nc,3),dtype=int)
        for il in range(nc):
            line = fc.next()
            n, kappa, occup = (int(line[:1]), int(line[2:4]), int(line[5:6]))
            nko[il,:] = array([n,kappa,occup])
        n_kappa_occup.append(nko)
        
    print 'iprint=', iprint
    
    # reading case.corewf
    fi = open(case+'.corewf','r')
    wavec=[]
    cType=[]
    cEne=[]
    for jatom in range(struct.nat):
        line = fi.next().split()
        nc0 = int(line[0])
        nc = len(n_kappa_occup[jatom])
        noccup0 = n_kappa_occup[jatom][0,2]
        #print noccup0, nc, nc0
        waveci=[]
        cTypei=[]
        cEnei=[]
        if (noccup0 and iprint[jatom]):
            for ic in range(nc):
                line = fi.next()
                m1 = re.search('CORE STATES = (\d\S)', line)
                m2 = re.search('CORE ENERGY=(.*) Ry', line)
                if m1 is not None and m2 is not None:
                    ctype = m1.group(1)
                    cenergy = float(m2.group(1))

                ucore = zeros((2,struct.npt[jatom]))
                for i in range(2):
                    n=0
                    while (n<struct.npt[jatom]):
                        line = fi.next()
                        dat = [line[3:3+19], line[3+19:3+2*19], line[3+2*19:3+3*19], line[3+3*19:3+4*19]]
                        dat = filter(lambda s: len(s)>1, dat) # throw away some empty strings
                        dat = map(float,dat)
                        ucore[i,n:n+len(dat)]= dat[:]
                        n += len(dat)
                waveci.append(ucore)
                print ctype, cenergy, shape(ucore)
                cTypei.append( ctype )
                cEnei.append(cenergy)
                
        wavec.append(waveci)
        cType.append(cTypei)
        cEne.append(cEnei)

    return (wavec, cType, cEne, iprint, n_kappa_occup)

def SelectCore(atoms, core, wavec, cType, cEne, iprint, n_kappa_occup, struct, plot_waves=False):
    core_waves=[]
    core_jatom=[]
    core_Ene=[]
    for jatom in range(struct.nat):
        if iprint[jatom]:
            for ic in range(len(n_kappa_occup[jatom])):
                Nr = struct.npt[jatom]  # Number of radial points from struct file
                r0 = struct.r0[jatom]
                Rmt = struct.rmt[jatom]
                dx = log(Rmt/r0)/(Nr-1)
                Rx = r0 * exp( arange(Nr)*dx ) # radial points
                Ag = wavec[jatom][ic][0]
                Bg = wavec[jatom][ic][1]
                norm = rint13g(1.0,1.0,Ag,Bg,Ag,Bg,dx,r0)

                if jatom in atoms and cType[jatom][ic]==core:
                    print 'Found desired core wave function'
                    core_waves.append( (Ag, Bg) )
                    core_jatom.append( jatom )
                    core_Ene.append(cEne[jatom][ic])
                    if plot_waves:
                        plt.plot(Rx, Ag, label=cType[jatom][ic])
                    
                print cType[jatom][ic], 'norm=', norm
            if plot_waves:
                plt.title('Core wave fnunctions')
                plt.legend(loc='best')
                plt.show()
    return (core_waves, core_jatom, core_Ene)

def SphericalIntegrals(ii, l1,m1,l2,m2):
    # ii=1 : <Y_{l1,m1}|e_r|Y_{l2,m2}>
    # ii=2 : <Y_{l1,m1}|r\nabla|Y_{l2,m2}>
    # ii=3 : <(r*\nabla)Y_{l1,m1}|e_r|(r*\nabla)Y_{l2,m2}>
    #
    def alm(l,m):
        return sqrt((l+m+1.)*(l+m+2.)/((2.*l+1.)*(2.*l+3.)))
    def flm(l,m):
        return sqrt((l+m+1.)*(l-m+1.)/((2.*l+1.)*(2.*l+3.)))
    cl=[0.5, -l2/2., l2*(l2+2.)/2.]
    dl=[0.5, (l2+1.)/2., (l2-1.)*(l2+1.)/2.]

    if abs(m1-m2)>1 or abs(l1-l2)>1:
        return zeros(3)
    
    if m1==m2:
        if l1==l2+1:
            return cl[ii]*2*flm(l2,m2)*array([0,0,1])
        elif l1==l2-1:
            return dl[ii]*2*flm(l1,m1)*array([0,0,1])
        else:
            return zeros(3)
    elif m1==m2-1:
        if l1==l2+1:
            return cl[ii]*alm(l2,-m2)*array([1,1j,0])
        elif l1==l2-1:
            return -dl[ii]*alm(l1,m1)*array([1,1j,0])
        else:
            return zeros(3)
    elif m1==m2+1:
        if l1==l2+1:
            return -cl[ii]*alm(l2,m2)*array([1,-1j,0])
        elif l1==l2-1:
            return dl[ii]*alm(l1,-m1)*array([1,-1j,0])

def Print(A):
    for i in range(shape(A)[0]):
        for j in range(shape(A)[1]):
            print "%7.4f %7.4f   " % (A[i,j].real, A[i,j].imag),
        print

def Get_momentum_P_quadrupole(lv, lc, valence_jatom, struct, valence_waves, core_waves, plot_waves=False):
    cPrv = zeros((3,3,2*lc+1,2*lv+1),dtype=complex)
    vPrc = zeros((3,3,2*lv+1,2*lc+1),dtype=complex)
    for ia,jatom in enumerate(valence_jatom):
        Nr = struct.npt[jatom]
        r0 = struct.r0[jatom]
        Rmt = struct.rmt[jatom]
        dx = log(Rmt/r0)/(Nr-1)
        Rx = r0 * exp( arange(Nr)*dx ) # radial points

        Av = valence_waves[ia][0]
        Bv = valence_waves[ia][1]
        Ac = core_waves[ia][0]
        Bc = core_waves[ia][1]
                            
        normc = rint13g(1.0,1.0, Ac, Bc, Ac, Bc, dx, r0)
        normv = rint13g(1.0,1.0, Av, Bv, Av, Bv, dx, r0)
        print 'norms=', normc, normv
        
        fAv = interpolate.UnivariateSpline(Rx, Av, s=0)
        fBv = interpolate.UnivariateSpline(Rx, Bv, s=0)
        dAv = fAv.derivative(n=1)(Rx)
        dBv = fBv.derivative(n=1)(Rx)

        fAc = interpolate.UnivariateSpline(Rx, Ac, s=0)
        fBc = interpolate.UnivariateSpline(Rx, Bc, s=0)
        dAc = fAc.derivative(n=1)(Rx)
        dBc = fBc.derivative(n=1)(Rx)

        PdP = rint13g(1.0,1.0, Ac, Bc, Rx*dAv, Rx*dBv, dx, r0)
        P_P = rint13g(1.0,1.0, Ac, Bc, Av, Bv, dx, r0)
        dPP = rint13g(1.0,1.0, Rx*dAc, Rx*dBc, Av, Bv, dx, r0)

        cRv = PdP-P_P
        vRc = dPP+2*P_P
        print 'PdP=', PdP, 'P_P=', P_P, 'dPP=', dPP
        print '<core|r d/dr|valence>', cRv, '<valence|d/dr r|core>', vRc

        for mc in range(-lc,lc+1):
            for mv in range(-lv,lv+1):
                Yc_er_er_Yv = zeros((3,3),dtype=complex)
                Yc_er_nabla_Yv = zeros((3,3),dtype=complex)
                #Yv_er_nabla_Yc = zeros((3,3),dtype=complex)
                Yv_nabla_er_Yc = zeros((3,3),dtype=complex)
                for m in range(-(lc+1),lc+2):
                    Yc_er_Ycp1 = SphericalIntegrals(0, lc,mc, lc+1,m)
                    Ycp1_er_Yv = SphericalIntegrals(0, lc+1,m, lv,mv)
                    Ycp1_nabla_Yv = SphericalIntegrals(1, lc+1,m, lv,mv)
                    Yc_er_er_Yv[:,:] += outer(Yc_er_Ycp1, Ycp1_er_Yv)
                    Yc_er_nabla_Yv[:,:] += outer(Yc_er_Ycp1, Ycp1_nabla_Yv)
                Yv_er_er_Yc = conjugate(Yc_er_er_Yv)
                cPrv[:,:,mc+lc,mv+lv] = Yc_er_er_Yv * cRv + Yc_er_nabla_Yv * P_P
        cPrv[:,:,:,:] *= -1j # because P = -i*nabla
        for i in range(3):
            for j in range(3):
                vPrc[i,j,:,:] = transpose(conjugate(cPrv[i,j,:,:]))
        
        if plot_waves:
            plt.plot(Rx, Ac*dAv, label='<core|d/dr|valence>')
            plt.plot(Rx, dAc*Av, label='<valence|d/dr|core>')
            plt.plot(Rx, Ac*Av/Rx, '-.', label='<core|1/r|valence>')
    if plot_waves:
        plt.legend(loc='best')
        plt.show()
    return (cPrv, vPrc)
    
def Get_momentum_P_dipole(lv,lc, valence_waves, core_waves, jatom, struct, plot_waves=False):
    cPv = zeros((3,2*lc+1,2*lv+1),dtype=complex)
    vPc = zeros((3,2*lv+1,2*lc+1),dtype=complex)

    Nr = struct.npt[jatom]
    r0 = struct.r0[jatom]
    Rmt = struct.rmt[jatom]
    dx = log(Rmt/r0)/(Nr-1)
    Rx = r0 * exp( arange(Nr)*dx ) # radial points

    Av = valence_waves[0]
    Bv = valence_waves[1]
    Ac = core_waves[0]
    Bc = core_waves[1]
                        
    normc = rint13g(1.0,1.0, Ac, Bc, Ac, Bc, dx, r0)
    normv = rint13g(1.0,1.0, Av, Bv, Av, Bv, dx, r0)
    print 'norms=', normc, normv
    
    fAv = interpolate.UnivariateSpline(Rx, Av, s=0)
    fBv = interpolate.UnivariateSpline(Rx, Bv, s=0)
    dAv = fAv.derivative(n=1)(Rx)
    dBv = fBv.derivative(n=1)(Rx)

    fAc = interpolate.UnivariateSpline(Rx, Ac, s=0)
    fBc = interpolate.UnivariateSpline(Rx, Bc, s=0)
    dAc = fAc.derivative(n=1)(Rx)
    dBc = fBc.derivative(n=1)(Rx)

    PdP = rint13g(1.0,1.0, Ac, Bc, dAv, dBv, dx, r0)
    PrP = rint13g(1.0,1.0, Ac, Bc, Av/Rx, Bv/Rx, dx, r0)
    dPP = rint13g(1.0,1.0, dAc, dBc, Av, Bv, dx, r0)


    for mc in range(-lc,lc+1):
        for mv in range(-lv,lv+1):
            cPv[:,mc+lc,mv+lv] = SphericalIntegrals(0, lc,mc, lv,mv)*(PdP-PrP)+SphericalIntegrals(1, lc,mc, lv,mv)*PrP
            vPc[:,mv+lv,mc+lc] = SphericalIntegrals(0, lv,mv, lc,mc)*(dPP-PrP)+SphericalIntegrals(1, lv,mv, lc,mc)*PrP
    
    cPv[:,:,:] *= -1j # because P = -i*nabla
    vPc[:,:,:] *= -1j # because P = -i*nabla
    
    print '<core|d/dr|valence>=', PdP, '<valence|d/dr|core>=', dPP, '<core|1/r|valence>=', PrP
    
    if plot_waves:
        plt.plot(Rx, Ac*dAv, label='<core|d/dr|valence>')
        plt.plot(Rx, dAc*Av, label='<valence|d/dr|core>')
        plt.plot(Rx, Ac*Av/Rx, '-.', label='<core|1/r|valence>')
        plt.legend(loc='best')
        plt.show()
    return (cPv, vPc)

def Print_Component_of_P(ii, cPv, lc, lv):
    ch=['x','y','z']
    print ch[ii]+'-component <core|P|valence>:'
    for mc in range(-lc,lc+1):
        for mv in range(-lv,lv+1):
            print "%6.3f %6.3f   " % (cPv[ii,mc+lc,mv+lv].real, cPv[ii,mc+lc,mv+lv].imag),
        print
    print ch[ii]+'-component <valence|P|core>:'
    for mc in range(-lc,lc+1):
        for mv in range(-lv,lv+1):
            print "%6.3f %6.3f   " % (vPc[ii,mv+lv,mc+lc].real, vPc[ii,mv+lv,mc+lc].imag),
        print

def Double_Core_For_spin(cPv, lc, lv):
    # core part needs to have spin part
    cPv_spin = zeros((3,2*(2*lc+1),val_size),dtype=complex)
    for ii in range(3):
        cPv_spin[ii,0:2*lc+1,0:2*lv+1] = cPv[ii,0:2*lc+1,0:2*lv+1]
        cPv_spin[ii,2*lc+1:2*(2*lc+1),2*lv+1:2*(2*lv+1)] = cPv[ii,0:2*lc+1,0:2*lv+1]
    return cPv_spin

#def Thompson_f(atom, q):
#    # Found at:
#    # http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php
#    Zthms={}
#    Zthms['Ni']=[12.8376,3.8785,7.292,  0.2565, 4.4438, 12.1763, 2.38,    66.3421, 1.0341]
#    Zthms['Ni2+']=[11.4166,3.6766,7.4005, 0.2449, 5.3442, 8.873,   0.9773,  22.1626, 0.8614]
#    Zthms['Ni3+']=[10.7806,3.5477,7.75868,0.22314,5.22746,7.64468, 0.847114,16.9673, 0.386044]
#    cs = Zthms[atom]
#    a = cs[0::2][:-1]
#    b = cs[1::2]
#    c = cs[-1]
#    
#    if type(q) is ndarray:
#        res = zeros(len(q))
#    else:
#        res = 0.0
#    for i in range(len(a)):
#        res += a[i]*exp(-b[i]*(q/(4*pi))**2)
#    return res + c

def ferm(x):
    y = 1/(exp(x)+1.)
    for i in range(len(x)):
        if x[i]>100: y[i]=0.0
        if x[i]<-100: y[i]=1.0
    return y

llc=0
def Prepare_G_widetilde(gfilename, beta, wbroad, kbroad=0.0):    
    data = loadtxt(gfilename).transpose()
    om = data[0]
    gim = data[2::2]
    gre = data[1::2]
    
    fcout = open('cout.dat', 'w')
    gcf = zeros((len(gim),len(om)),dtype=complex) # final Green's function
    # Now multiply with the proper fermi function
    fermi = ferm(-om*beta)
    for ib in range(len(data)/2):
        gim[ib,:] = gim[ib,:] * fermi
    savetxt('gim.dat', vstack((om,gim)).transpose())
    # Now we perform Kramars-Kronig in the terminal, and we also broaden the result
    for ib in range(len(data)/2):
        end_name = '.dat'
        cmd = dmfe.ROOT+'/skrams -cn '+str(ib+2)+' gim.dat > gre'+str(ib+2)+end_name
        #subprocess.call(cmd,shell=True,stdout=sys.stdout,stderr=sys.stderr)
        subprocess.call(cmd,shell=True,stdout=fcout,stderr=fcout)
        if wbroad>0:
            cmd = dmfe.ROOT+'/broad -w '+str(wbroad)+' -k '+str(kbroad)+' gre'+str(ib+2)+end_name+' > gre'+str(ib+2)+'.broad'
            end_name = '.broad'
            #subprocess.call(cmd,shell=True,stdout=sys.stdout,stderr=sys.stderr)
            subprocess.call(cmd,shell=True,stdout=fcout,stderr=fcout)
        dt = loadtxt('gre'+str(ib+2)+end_name).transpose()
        gcf[ib,:] = dt[1]+dt[2]*1j
        
        os.remove('gre'+str(ib+2)+'.dat')
        if wbroad>0:
            os.remove('gre'+str(ib+2)+end_name)
    os.remove('gim.dat')
    #global llc
    #shutil.copyfile('gim.dat','gim.dat'+str(llc))
    #llc+=1
    return (om,gcf)


if __name__ == '__main__':
    chemical_symbol=['S','P','D','F','G','H']
    l_from_chemical_symbol = {'S':0, 'P':1, 'D':2, 'F':3, 'G':4, 'H':5}
    au2ang = 0.52917721067
    speed_of_light = 137.0359895
    Ry2eV = 13.6058
    
    Sprint=False
    plot_waves=False
    lmbda=1.4855
    #hkl=[1,0,1]
    hkl=[1,0,5]
    
    #hkl=[0,1,1]
    #hkl=[1,1,1]
    #hkl=[0,1,5]
    
    atoms = [0,1]
    wbroad=1.6
    kbroad=0.0
    Ecc=8329.5
    beta=100
    #EF=8.559000606
    Z_for_4p = 0.525 #*0.8
    core = '1S'
    valence_l = 1
    om = linspace(-5,40,400)
    core_hole=['imp.0/core_hole.dat', 'imp.1/core_hole.dat']
    #### Finished with input
    #######################################

    
    Uch=[]
    Pch=[]
    U_aver=0
    for ia in range(len(core_hole)):
        dat = loadtxt(core_hole[ia]).transpose()
        #for i in range(len(dat[0])):
        #    dat[0][i] += 0.5*i#????
        #dat[0] *= 0.0
        Uch.append(dat[0])
        Pch.append(dat[1])
        U_aver += sum(dat[0])/len(dat[0])
    U_aver *= 1./len(core_hole)
    print 'Uch=', Uch
    print 'Pch=', Pch
    Uch = array(Uch)
    Pch=array(Pch)
    Uch -= U_aver

    w2k = utils.W2kEnvironment()
    dmfe = utils.DmftEnvironment()
    struct = struct1.Struct(w2k.case)
    
    a = struct.a * au2ang
    b = struct.b * au2ang
    c = struct.c * au2ang

    print 'a,b,c [A]=', a, b, c
    print 'wavelength of light [A]=', lmbda
    print 'h,k,l=', hkl

    qs=array([2*pi*hkl[0]/a,2*pi*hkl[1]/b,2*pi*hkl[2]/c])
    q_norm = sqrt(sum(qs**2))

    f_thompson = 0j
    f_dynamic = 0j
    for iat in range(len(struct.mult)):
        phs = 0
        for mu in range(struct.mult[iat]):
            qr = dot(hkl,struct.pos[iat][mu])
            phs += exp(1j*2*pi*qr)
            print struct.aname[iat], struct.pos[iat][mu], qr
        print 'phase sum for above atoms =', phs
        fth = Thompson_f(struct.aname[iat], q_norm)
        print 'Thompson part for ', struct.aname[iat], 'is', fth, 'with phase is', fth*phs
        f_dyn = fc_dynamic(struct.aname[iat])
        if iat in atoms:
            Fc_dynamic = f_dyn
        print 'Dynamic-atomic part for', struct.aname[iat],'is', f_dyn, 'with phase is', f_dyn*phs
        #if struct.aname[iat]=='Nd':
        #    fth -= 2.1  # correction due to f' of Nd.
        #    # see : http://physics.nist.gov/PhysRefData/FFast/html/form.html
        #    print "correcting for energy dependent f': fth=", fth
        f_thompson += fth*phs
        f_dynamic += f_dyn*phs
    
    
    #f_thompson = Thompson_f('Ni2+', q_norm)
    print 'total f_Thompson=', f_thompson, 'for q=', q_norm
    print 'dynamic part withouth correlated atom ', f_dynamic
    print 'dynamic part for correlated atom', Fc_dynamic

    # How to define azimuthal angle. Where does it vanish?
    eq = qs/q_norm
    psi_is_zero = zeros(3)  # the new x-axis for azimuthal angle. It must be orthogonal to eq.
    which_axis=1
    if 0 in hkl:
        which_axis = hkl.index(0)
    psi_is_zero[which_axis]=1
    psi_is_zero = cross(cross(eq,psi_is_zero),eq)
    psi_is_zero *= 1/sum(psi_is_zero**2)
    print 'psi_is_zero=', psi_is_zero



    if True:
        scatter = Scattering(hkl,(a,b,c),lmbda,Sprint)
        # computing sigma and pi polarizations for many thetas, and then sorting with respect to azimuthal angle
        psi=[]
        theta=[]
        sig=[]
        pis=[]
        for t in linspace(0,pi,30):
            if (Sprint):
                print '**********************************************'
                print 'theta=', t/pi, '*pi'
            (sigm, pi_in, pi_out, eks) = scatter.ScatteringVector(t)
            if len(sigm)>0:
                if (Sprint):
                    print '********* Results for theta=', t/pi, '*pi:'
                    print 's=', sigm[0], 'p=', pi_out[0], 'ek=', eks[0]
                for isol in range(len(sigm)):
                    s, p, ek = sigm[isol], pi_out[isol], eks[isol]
                    if (Sprint): print 'ss=', sigm, 'pi,pi=', pi_in,pi_out
                    sig.append(s)
                    pis.append(p)
                    cz = dot(ek,eq)             # cos(z)
                    sz_cp = dot(ek,psi_is_zero) # sin(z)*cos(psi)
                    cpsi = sz_cp/sqrt(1-cz**2)  # cos(psi)
                    ps = arccos(cpsi)/pi        # psi
                    ps = ps * scatter.Azimuthal_single_valued(t)
                    psi.append(ps)
                    theta.append(t)
        # now sorting to have azimuthal angle sorted
        indx = range(len(psi))
        indx = sorted(indx, key=lambda i: psi[i])
        psi = [psi[indx[i]] for i in range(len(indx))]
        theta = [theta[indx[i]] for i in range(len(indx))]
        sig = [sig[indx[i]] for i in range(len(indx))]
        pis = [pis[indx[i]] for i in range(len(indx))]
        # Now done with polarizations.

        
    else:
        # computing sigma and pi polarizations for many thetas, and then sorting with respect to azimuthal angle
        psi=[]
        theta=[]
        sig=[]
        pis=[]
        for t in linspace(0,pi,30):
            if (Sprint):
                print '**********************************************'
                print 'theta=', t/pi, '*pi'
            (sigm, pi_in, pi_out, eks) = ScatteringVector(t, lmbda, [a,b,c], hkl)
            if len(sigm)>0:
                if (Sprint):
                    print '********* Results for theta=', t/pi, '*pi:'
                    print 's=', sigm[0], 'p=', pi_out[0], 'ek=', eks[0]
                for isol in range(len(sigm)):
                    s, p, ek = sigm[isol], pi_out[isol], eks[isol]
                    if (Sprint): print 'ss=', sigm, 'pi,pi=', pi_in,pi_out
                    sig.append(s)
                    pis.append(p)
                    cz = dot(ek,eq)             # cos(z)
                    sz_cp = dot(ek,psi_is_zero) # sin(z)*cos(psi)
                    cpsi = sz_cp/sqrt(1-cz**2)  # cos(psi)
                    ps = arccos(cpsi)/pi        # psi
                    ps = ps * Azimuthal_single_valued(t, lmbda, hkl, (a,b,c))
                    psi.append(ps)
                    theta.append(t)
        # now sorting to have azimuthal angle
        indx = range(len(psi))
        indx = sorted(indx, key=lambda i: psi[i])
        psi = [psi[indx[i]] for i in range(len(indx))]
        theta = [theta[indx[i]] for i in range(len(indx))]
        sig = [sig[indx[i]] for i in range(len(indx))]
        pis = [pis[indx[i]] for i in range(len(indx))]
        # Now done with polarizations.
    
    # Core is geting ready
    (wavec, cType, cEne, iprint, n_kappa_occup) = read_core(w2k.case, struct)
    (core_waves, core_jatom, Ec) = SelectCore(atoms, core, wavec, cType, cEne, iprint, n_kappa_occup, struct, plot_waves=plot_waves)
    # now projector
    (wave, Type) =  read_projector()
    (valence_waves, valence_jatom) = SelectValence(atoms, valence_l, wave,Type, plot_waves=plot_waves)

     
    (atms, Lsa, icpsa, Siginds, T2Cs) = ReadIndmfl(w2k.case)
    # atoms in indmfl, in structure and those that should be used here, should be connected.
    # which atoms (equivalent and inequivalent) will be considered?
    latom=0
    cixs=[]
    siginds=[]
    T2cs=[]
    ifirst=[]
    for iat in range(len(struct.mult)):
        if iat in atoms: ifirst.append(latom)
        for mu in range(struct.mult[iat]):
            latom += 1
            if iat in atoms:
                iatom = atms.index(latom)
                print 'Lsa[iatom]=', Lsa[iatom]
                for j in range(len(Lsa[iatom])):
                    if Lsa[iatom][j]==valence_l:
                        icix = icpsa[iatom][j]
                        cixs.append(icix)
                        siginds.append(Siginds[icix-1])
                        T2cs.append( T2Cs[icix-1])
    
    # In these files, we expect the green's function of the valence electrons
    gcfiles = [w2k.case+'.gc'+str(icix) for icix in cixs]
    siginds = [ array(Sigind) - min(ravel(Sigind)) for Sigind in siginds]


    print 'cixs=', cixs
    print 'ifirst=', ifirst
    #print siginds
    #print T2cs
    print gcfiles
    print 'Ec=', Ec
    Ec = array(Ec)*Ry2eV #+ EF
    Ec0 = sum(Ec)/len(Ec)
    Ec = array(Ec)-Ec0
    print 'Ec0=', Ec0

    # core and valence l
    lv = valence_l
    lc = l_from_chemical_symbol[core[1]]
    print 'lv=', lv, 'lc=', lc

    # core is always relativistic
    cT2C = matrix(trafoso.trafoso(l=lc))
    print 'T2c in the core:', cT2C
    
    vdim = len(T2cs[0])
    
    Fc = zeros((len(atoms),vdim,vdim,len(om)),dtype=complex)
    iat=0
    for ia,jatom in enumerate(atoms):
        first_atom = ifirst[ia]
        
        (cPv, vPc) = Get_momentum_P_dipole(lv,lc, valence_waves[ia], core_waves[ia], jatom, struct, plot_waves)
        
        for ii in range(3):
            Is_Hermitian = sum(abs(cPv[ii,:,:]-conjugate(transpose(vPc[ii,:,:]))))
            if (Is_Hermitian>1e-3):
                print 'ERROR: It seems the momentum matrix elements appears non-hermitian!', Is_Hermitian

        cPv_spin = zeros(shape(cPv),dtype=complex)
        
        T2v = T2cs[first_atom]

        valence_spin_deg=1
        #if len(T2v)==(2*valence_l+1): valence_spin_deg=2
        
        core_size=2*(2*lc+1)
        core_deg=1
        if lc>0:
            # core part needs to have spin part, because core splits into lc-1/2 and lc+1/2
            cPv_spin = Double_Core_For_spin(cPv, lc, lv, double_valence)
            if len(T2v)==(2*valence_l+1):
                dim = 2*valence_l+1
                T2C_spin = zeros((2*dim,2*dim),dtype=complex)
                T2C_spin[0:dim,0:dim] = T2v[:,:]
                T2C_spin[dim:2*dim,dim:2*dim] = T2v[:,:]
                T2v = T2C_spin
            elif len(T2v)==2*(2*valence_l+1):
                pass # It already contains spin
            else:
                print 'ERROR: There is a problem with dimensions of T2cs matrix', len(T2v), (2*valence_l+1), 2*(2*valence_l+1)
                sys.exit(0)
            
            for ii in range(3):
                cPv_spin[ii,:,:] = cT2C * matrix(cPv_spin[ii,:,:]) * matrix(T2v).H        
        else:
            core_size = 1
            core_deg = 2
            for ii in range(3):
                cPv_spin[ii,:,:] = matrix(cPv[ii,:,:]) * matrix(T2v).H
        
        dir=['x','y','z']
        for ii in range(3):
            print 'Left part of the Scattering tenzor for '+dir[ii]+' direction is:'
            Print(cPv_spin[ii,:,:])

        # now over all equivalent atoms in the unit cell
        for im in range(struct.mult[jatom]):
            iat = first_atom + im
            
            gfilename = gcfiles[iat]
            gfilename = gcfiles[iat]
            w,gc = Prepare_G_widetilde(gcfiles[iat], beta, wbroad, kbroad)


            #plt.plot(w, real(gc[0]))
            #plt.plot(w, imag(gc[0]))
            #plt.show()
            
            pos=struct.pos[jatom][im]
            print 'Scattering computed for atom', struct.aname[jatom], 'at position', pos
            exp_ph=exp(dot(hkl,pos)*2*pi*1j)
            print 'Phase=', exp_ph
            sigind = siginds[iat]
            
            Gv = zeros((vdim,vdim,len(om)),dtype=complex)
            for i in range(vdim):
                for j in range(vdim):
                    for ni in range(len(Uch[ia])):
                        #w_shift = w-Ec[ia]+nch[ni]*Uch
                        w_shift = w-Ec[ia]+Uch[ia,ni]
                        
                        gr = interpolate.UnivariateSpline(w_shift, real(gc[sigind[i,j]]),s=0)
                        gi = interpolate.UnivariateSpline(w_shift, imag(gc[sigind[i,j]]),s=0)
                        
                        ind = range(len(om))
                        ind_small = filter(lambda ii: om[ii]>w_shift[0] and om[ii]<w_shift[-1], ind)
                        om_small = array([om[ii] for ii in ind_small])
                        gc_small = (gr(om_small) + gi(om_small)*1j)
                        gc_full = zeros(len(om),dtype=complex)
                        for ii in range(len(ind_small)):
                            gc_full[ind_small[ii]] = gc_small[ii]
                            
                        Gv[i,j,:] += conjugate(gc_full)*Pch[ia,ni]  # because we require advanced, not retarder G, and the sign convention seems to be non-standard

                        #if i==j:
                        #    plt.plot(w_shift, real(gc[sigind[i,j]])*Pch[ia,ni], label='n='+str(Uch[ia,ni]/7.)+' '+dir[i])
                        #    plt.plot(w_shift, imag(gc[sigind[i,j]])*Pch[ia,ni], label='n='+str(Uch[ia,ni]/7.)+' '+dir[i])
                        #    plt.plot(om, real(gc_full)*Pch[ia,ni], label='n='+str(Uch[ia,ni]/7.)+' '+dir[i])
                        #    plt.plot(om, imag(gc_full)*Pch[ia,ni], label='n='+str(Uch[ia,ni]/7.)+' '+dir[i])
                        #    plt.legend(loc='best')
                        #    plt.show()

            #for i in range(vdim):
            #    plt.plot(om,real(Gv[i,i,:]))
            #    plt.plot(om,imag(Gv[i,i,:]))
            #plt.show()
            
            ii=1
            for i in range(vdim):
                axs=[]
                for j in range(vdim):                          #  cPv_spin[i,lmc,lmv]
                    p_G = dot(cPv_spin[i,:,:],Gv[:,:,:])       #  p_G[lmc,lmv,iom]
                    p_c = transpose(conjugate(cPv_spin[j,:,:]))#  p_c[lmv,lmc]
                    p_G = p_G.transpose((2,0,1))               #  p_G[iom,lmc,lmv]
                    p_G_p = dot(p_G,p_c)                       #  p_G_p[iom,lmc,lmcp] = p_G[iom,lmc,:]*[:,lmcp]
                    for iom in range(len(om)):
                        Fc[ia,i,j,iom] += p_G_p[iom].trace() * (core_deg * valence_spin_deg * Z_for_4p * Ry2eV) * exp_ph
                    
                    if i==j:
                        Fc[ia,i,i,:] += exp_ph*Fc_dynamic  # value at infinity from NIST
                        
                    #if j==0:
                    #    ax=plt.subplot(vdim,vdim,ii)
                    #else:
                    #    ax=plt.subplot(vdim,vdim,ii,sharey=ax)
                    #ii+=1
                    #plt.plot(om, real(Fc[ia,i,j,:]/exp_ph), label='Re dFc '+dir[i]+dir[j])
                    #plt.plot(om, imag(Fc[ia,i,j,:]/exp_ph), label='Im dFc '+dir[i]+dir[j])
            #plt.subplots_adjust(hspace=0,wspace=0)
            #plt.show()

    A_XRAS = zeros(len(om))
    fig, axs = plt.subplots(vdim, vdim, sharex='col', sharey='row')
    #fig.suptitle('$F^{resonant}(E)$', fontsize=20)
    fig.subplots_adjust(hspace=0,wspace=0)
    col=['g','b','c','r','m']
    for i in range(vdim):
        for j in range(vdim):
            for ia in range(len(atoms)-1,-1,-1):
                deg = struct.mult[ia]
                phs = (-1)**(ia+1)
                axs[i,j].plot(om, real(phs*Fc[ia,i,j,:])/deg, col[ia]+'-', lw=1,label='Re Ni$_'+str(2-ia)+'$')
                axs[i,j].plot(om, imag(phs*Fc[ia,i,j,:])/deg, col[ia]+'-', lw=2,label='Im Ni$_'+str(2-ia)+'$')

                A_XRAS += imag(phs*Fc[ia,i,j,:])/deg
                
            #axs[i,j].text(1,2,'$ F_{'+dir[i]+dir[j]+'}$',fontsize='xx-large')
            if (i!=vdim-1): plt.setp(axs[i,j].get_xticklabels(), visible=False)
            if (j!=0): plt.setp(axs[i,j].get_yticklabels(), visible=False)
            axs[i,j].set_xlim([0.01,39])
            axs[i,j].set_ylim([-8.8,5.5])
            axs[i,j].set_xticks([0,10,20,30])
            if (i==0 and j==2):
                axs[i,j].legend(loc=(0.7,0.3))
            if (i==2 and j==1):
                axs[i,j].set_xlabel('$Energy(eV)$',fontsize='xx-large')
    savetxt('XRAS.dat',vstack((om,A_XRAS)).transpose())
    plt.show()
    
    #plt.plot((om-15.4+8350)/1000., A_XRAS, 'g-', lw=2)
    #plt.xlim([8.3305,8.3695])
    #plt.xlabel('Energy (keV)')
    #plt.ylabel('$Im(F_{xx}+F_{yy}+F_{zz})$')
    #plt.show()

    
    Fct = zeros((vdim,vdim,len(om)),dtype=complex)
    for ia in atoms:
        Fct[:,:,:] += Fc[ia,:,:,:]

    if False and wbroad>0.0:
        Fc2 = zeros(shape(Fct),dtype=complex)
        data = om
        for i in range(vdim):
            for j in range(vdim):
                data = vstack( (data,real(Fct[i,j,:]),imag(Fct[i,j,:])) )
        savetxt('in_broad.dat',data.transpose())
        cmd = '~/dbin/broad -w 0.5 in_broad.dat > out_broad'
        subprocess.call(cmd,shell=True,stdout=sys.stdout,stderr=sys.stderr)
        data = loadtxt('out_broad').transpose()
        print sum(abs(data[0]-om))
        cdat = data[1::2]+data[2::2]*1j
        for i in range(vdim):
            for j in range(vdim):
                Fc2[i,j,:] = cdat[i*vdim+j,:]
        Fct = Fc2

    if hkl==[1,0,5]:
        expdat = loadtxt('exp_1_0_5.dat').transpose()
    elif hkl==[1,0,1] or hkl==[0,1,1]:
        expdat = loadtxt('exp_1_0_1.dat').transpose()
        Ecc += 0.8
    elif hkl==[0,1,5]:
        expdat = loadtxt('exp_0_1_5.dat').transpose()

    for i in range(len(psi)):
        ps = psi[i]
        sg = sig[i]
        po = pis[i]
        fc = zeros(len(om),dtype=complex)
        for i in range(3):
            for j in range(3):
                fc += sg[i]*Fct[i,j,:]*sg[j]
        plt.plot(om,real(fc),'-', label='Re')
        plt.plot(om,imag(fc),'-',label='Im')
        plt.plot(om,ones(len(om))*f_thompson,':', lw=2, label='thompson')
        plt.plot(om,ones(len(om))*real(f_dynamic), ':', lw=2, label='dynamic Re')
        plt.plot(om,ones(len(om))*imag(f_dynamic), ':', lw=2, label='dynamic Im')
        plt.legend(loc='best')
        plt.show()
        fc += f_thompson + f_dynamic
        
        plt.title('$ s*fc*s \;\; for\;\; \\psi='+str(ps)+'*\\pi $')
        plt.plot(om+Ecc, abs(fc)**2/21, label='$|F|^2$', lw=2)
        #plt.plot(expdat[0],(expdat[1]-5.0)*9.5, 'o', label='exp') # [1,0,5]

        if hkl==[1,0,5]:
            plt.plot(expdat[0],(expdat[1])*4.08, 'o', alpha=0.5, label='Exp') # [1,0,5]
            plt.xlim([8330,8370])
        elif hkl==[1,0,1] or hkl==[0,1,1]:
            #plt.plot(expdat[0], (expdat[1]-3.5)*0.5, 'o', label='exp')  # [0,1,1]
            plt.plot(expdat[0], (expdat[1])*1.72, 'o', alpha=0.5, label='Exp')  # [0,1,1]
        elif hkl==[0,1,5]:
            plt.plot(expdat[0],(expdat[1])*4.13, 'o', alpha=0.5, label='Exp') # [1,0,5]
            
        #plt.grid()
        plt.legend(loc='best')
        plt.xlabel('$Energy (eV)$', fontsize='xx-large')
        plt.ylabel('$I_{RXS}$', fontsize='xx-large')
        plt.show()
        
