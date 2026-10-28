#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
from numpy import *
import sys,os,re
import json
import ast
import numpy as np
from scipy.special import zeta
from scipy import interpolate
import builtins

# All variables set through exec call need to be global for Python3 to set them. Hence we are forced to work with global variables.
cix='actqmc.cix'
Delta='Delta.inp'
Gf='Gf.out'
Sig='Sig.out'
#beta,mu=0,0
#nf,mu,Ekin,Epot,TrSigmaG,mom=0,0,0,0,0,0

def ReadPARAMS2(fparams, fout):
    with open(fparams, 'r') as fp:
        lines = fp.readlines()

    stripped = []
    for line in lines:
        n = line.find('#')
        if n>=0:
            line = line[:n]
        stripped.append(line)
    content = ''.join(stripped).strip()

    if content.startswith('{'):
        params = json.loads(content)
        mapping = {'fSig': 'Sig', 'fGf': 'Gf', 'Delta': 'Delta', 'cix': 'cix', 'beta': 'beta', 'mu': 'mu'}
        for key_in, key_out in mapping.items():
            if key_in in params:
                globals()[key_out] = params[key_in]
                print('executed:', key_out+'='+repr(params[key_in]), file=fout)
    else:
        for line in stripped:
            dat = line.split()
            if len(dat)>=2:
                m = re.search('cix|Delta|Gf|Sig', dat[0])
                if m is not None :
                    cmn = dat[0]+'="'+dat[1]+'"'
                    exec(cmn, globals())
                    print('executed:', cmn, file=fout)
                m = re.search('beta|mu',dat[0])
                if m is not None :
                    cmn = dat[0]+'='+dat[1]
                    exec(cmn, globals())
                    print('executed:', cmn, file=fout)
    print('beta=', beta, 'mu=', mu, 'cix=', cix, 'Delta=', Delta, file=fout)
    return (cix,Delta,Gf,Sig,beta,mu)

def ReadCix(cix, fout):
    fcx = open(cix,'r')
    next(fcx);next(fcx)
    header = next(fcx).split()
    Deg={}

    if len(header) == 4:
        (Nclust,Nm,Ns,maxm) = list(map(int, header))
        next(fcx)
        for i in range(Ns):
            dat = list(map(int,next(fcx).split()))
            dim = dat[1]
            if (dim>1): 
                print('WARNING: The impurity functional cmpEimp2.py works for now only when baths are diagonal. You should not trust total free energy. Needs some development...!')
            for j in range(dim):
                for k in range(dim):
                    if k==j:
                        ii = dat[2+k+j*dim]
                        if ii in Deg:
                            Deg[ii] += 1
                        else:
                            Deg[ii]=1
        next(fcx)
        ss = next(fcx)
    else:
        # New Bold-CTQMC cix format:
        # cluster_size, number of states, number of psi^+ operators,
        # number of baths, number of blocks, maximum_matrix_size
        Ns = int(header[3])
        next(fcx)
        for i in range(Ns):
            dat = next(fcx).split()
            ii = int(dat[3])
            if ii in Deg:
                Deg[ii] += 1
            else:
                Deg[ii] = 1
        next(fcx)
        ss = next(fcx)

    print('Deg=', Deg, file=fout)
    if ss[:11]=='FL_FROM_IFL':
        for i in range(Ns):
            next(fcx)
        next(fcx)
        ss = next(fcx)
    Eks = list(map(float,ss.split()))
    print('Eks=', Eks, file=fout)
    if len(Deg)>0:
        deg = zeros(builtins.max(list(Deg.keys())) + 1, dtype=int)
        for ii,val in Deg.items():
            deg[ii] = val
    else:
        deg = zeros(0, dtype=int)
    return (deg, array(Eks))


def ReadGandSigandDelta(Gf, Sig, Delta, fout):
    Delta_data = loadtxt(Delta).transpose()
    omD = Delta_data[0]
    Dlta = Delta_data[1::2]+Delta_data[2::2]*1j
    
    Gf_data = loadtxt(Gf).transpose()
    om = Gf_data[0]
    Gfc = Gf_data[1::2]+Gf_data[2::2]*1j

    Sig_data = loadtxt(Sig).transpose()
    om3 = Sig_data[0]
    Sigc = Sig_data[1::2]+Sig_data[2::2]*1j
    
    print('shapes=', shape(om), shape(omD), shape(om3), file=fout)
    #print 'Diff=', sum(abs(om-omD))
    #if sum(abs(om-omD))>1. or sum(abs(om-om3))>1.:
    #    print >> fout, 'frequency mesh of ',Gf,' and ',Sig, 'and', Delta,' are not equal', om-om2, om-om3
    
    first_line = open(Gf, 'r').readline()
    dat = first_line.split()
    for d in dat:
        m = re.search('nf|mu|T|TrSigmaG|Ekin|Epot|mom', d)
        if m is not None:
            print('executed:', d, file=fout)
            exec(d, globals())
            
    print('nf=', nf, 'mu=', mu, 'Ekin=', Ekin, 'Epot=', Epot, file=fout)
    return (om, omD, Gfc, Sigc, Dlta, nf, mu, Ekin, Epot, TrSigmaG, mom) 


def ferm(x):
    if x>300:
        return 0
    elif x<-300:
        return 1
    else:
        return 1/(exp(x)+1)

def GetHighFrequency(CC,om):
    " Approximates CC ~  A/(i*om-C) "
    A = 1./imag(1/(CC[-1]*om[-1]))
    C = -A*real(1./CC[-1])
    return (A, C)

def ReturnHighFrequency(A,C,beta,E):
    " Returns the value of Tr(A/(iom-C) 1/(iom-E) )"
    return A*(ferm(E*beta)-ferm(C*beta))/(E-C)

def LogGimp(omega,Gfc,EimpS,beta,Deg,fout):
    " Tr(log(-G_imp))"
    def NonIntF0(beta,Ene):
        if beta*Ene>200: return 0.0
        if beta*Ene<-200: return Ene
        return -log(1+exp(-Ene*beta))/beta
    
    lnGimp=0.
    for i in range(len(Deg)): 
        print('EimpS=', EimpS[i], ', real(iom-1/G)', real(-1/Gfc[i][-1]), file=fout)
        if EimpS[i] < 1000 :
            Gd = Gfc[i,:]
            eimps = EimpS[i]
            
            CC = omega*1j-eimps-1/Gd
            A,C = GetHighFrequency(CC,omega)

            ff = log(-Gd)-log(-1./(omega*1j-eimps)) - 1./(omega*1j-eimps)*A/(omega*1j-C)
            lnGimp_i = 2*sum(ff.real)/beta+NonIntF0(beta,eimps)+ReturnHighFrequency(A,C,beta,eimps)
            lnGimp_i *= Deg[i]
            
        else:
            sum_imp = 0.0
            Fnint_imp = 0.0
            lnGimp_i = 0.0
        lnGimp += lnGimp_i
        print('lnGimp: Deg[i]=', Deg[i], 'Tr(log(-Gimp))=', lnGimp_i/Deg[i], 'Tr(log(-Gimp))*deg=', lnGimp_i, file=fout)
    return lnGimp


def LogGimpSig(omega, Gfc, Sigc, Vdc, EimpS, beta, Deg, fout):
    " Tr(log(1+(Sig-Vdc)G))"
    def NonIntF0(beta,Ene):
        if beta*Ene>200: return 0.0
        if beta*Ene<-200: return Ene
        return -log(1+exp(-Ene*beta))/beta
    
    lnGSimp=0.
    for i in range(len(Deg)):#len(Gfc)):
        Gd = Gfc[i,:]
        Sg = Sigc[i,:]
        eimps = EimpS[i]
        s_oo = Sigc[i,-1].real
        vdc = Vdc[i]
        
        #print 'shape(Gd)=', shape(Gd)
        ff0 = log((Sg-vdc)*Gd+1.)
        ff1 = (s_oo-vdc)/(omega*1j-eimps)
        
        lnGimp_i = 2*sum((ff0-ff1).real)/beta + (s_oo-vdc)*ferm(eimps*beta)
        lnGimp_i *= Deg[i]

        lnGSimp += lnGimp_i
        print('lnGSimp: Deg[i]=', Deg[i], 'Tr(log(1+Sig*Gimp))=', lnGimp_i/Deg[i], 'Tr(log(1+Sig*Gimp))*deg=', lnGimp_i, file=fout)

        #plot(omega, real(ff0-ff1), label='diff-real')
        #plot(omega, imag(ff0-ff1), label='diff-imag')
        #legend(loc='best')
        #show()
    return lnGSimp
    
def fTrDeltaG(om, omD, Dlta, Gfc, Deg, EimpS, beta, fout):
    def GetInterpolation(omD, Delta):
        "Returns derivative  d Delta/d om_i"
        Dr = interpolate.UnivariateSpline(omD, Delta.real, s=0,k=3)
        Di = interpolate.UnivariateSpline(omD, Delta.imag, s=0,k=3)
        return (Dr, Di)
    lngh = builtins.min(len(Dlta),len(Gfc),len(Deg))
    GDf=0.0
    for i in range(lngh):
        if len(omD)!=len(om):
            # interpolating Delta
            Dr, Di = GetInterpolation(omD, Dlta[i])
            Delta = Dr(om)+Di(om)*1j
        else:
            Delta = Dlta[i]
        (A,C) = GetHighFrequency(Delta,om)
        print('A_d=', A, 'C_d=', C, file=fout)
        eimps = EimpS[i]
        GD1 = Gfc[i]*Delta - 1./(om*1j-eimps)*A/(om*1j-C)
        dGD = Deg[i]*(2*sum(GD1.real)/beta + ReturnHighFrequency(A,C,beta,eimps))
        GDf += dGD
        #print >> fout, 'dTr(Delta*G)=', dGD/Deg[i]
    return GDf

def fTrSigmaG(om, Gfc, Sigc, Deg, EimpS, beta, fout):
    SDf=0.0
    for i in range(len(Deg)):
        if EimpS[i] < 1000 :
            Gd = Gfc[i,:]
            eimps = EimpS[i]
            Sg = Sigc[i,:]
            s_oo = Sigc[i,-1].real
            Sg0 = Sg[:] - s_oo   # Sigma-s_oo
            (A,C) = GetHighFrequency(Sg0,om)
            #print '::: A,C, for Sigma=', A, C
            C=1.
            ff = Gd[:]*Sg0[:] - 1./(om*1j-eimps)*A/(om*1j-C)
            SDf_i = 2*sum(ff.real)/beta + ReturnHighFrequency(A,C,beta,eimps)
            SDf_i *= Deg[i]
        else:
            SDf_i = 0.0
        SDf += SDf_i
        #print >> fout, 'Tr((Sigma-s_oo)*G): Deg[i]=', Deg[i], 'Tr((Sigma-s_oo)*G)=', SDf_i/Deg[i], 'Tr((Sigma-s_oo)*G)*deg=', SDf_i
    return SDf
            
def CmpImpurityCorrection2(om, omD, Dlt, Gfc, Deg, EimpS, beta, fout):
    def GetInterpolation(omD, Delta):
        "Returns derivative  d Delta/d om_i"
        Dr = interpolate.UnivariateSpline(omD, Delta.real, s=0,k=3)
        Di = interpolate.UnivariateSpline(omD, Delta.imag, s=0,k=3)
        return (Dr, Di)

    Delta_inf=[]
    for b in range(len(Dlt)):
        Delta_inf.append( [Dlt[b,-1].real*omD[-1]**2, Dlt[b,-1].imag*omD[-1] ] )
    
    cxx = 0.0
    lngh = builtins.min(len(Dlt),len(Gfc),len(Deg))
    
    oms = om
    Gfs = Gfc
    if oms[0]<omD[0] and abs(oms[0]-omD[0])<1e-6: oms[0]=omD[0]
    
    print('oms[-1]=', oms[-1], 'omD[-1]=', omD[-1], file=fout)
    
    for i in range(lngh):
        #print('shape(omD)=', shape(omD), 'shape(Dlt[i])=', shape(Dlt[i]), 'len(oms)=', len(oms))
        Dr, Di = GetInterpolation(omD, Dlt[i])
        dDlt=zeros(len(oms),dtype=complex)
        for ix,x in enumerate(oms):
            if omD[0]<=x<=omD[-1]:
                dDlt[ix] = Dr.derivatives(x)[1]+Di.derivatives(x)[1]*1j
            else:
                dDlt[ix] = -2*Delta_inf[i][0]/x**3 - Delta_inf[i][1]/x**2 * 1j
        
        CC = oms*dDlt
        #savetxt( 'brisx', vstack((oms, real(CC), imag(CC))).transpose() )
        #print 'shape(dDlt)=', shape(dDlt), 'shape(CC)=', shape(CC), 'shape(Gfs)=', shape(Gfs)
        (A,C) = GetHighFrequency(CC,oms)
        print('A=', A, 'C=', C, file=fout)
        eimps = EimpS[i]
        GD1 = Gfs[i]*CC - 1./(oms*1j-eimps)*A/(oms*1j-C)
        dcxx = Deg[i]*(2*sum(GD1.real)/beta + ReturnHighFrequency(A,C,beta,eimps))
        cxx += dcxx
        print('dcxx=', dcxx/Deg[i], file=fout)
        
    print('cxx=', cxx, file=fout)
    return cxx

# Normal-state matrix implementation for the direct-basis Bold CIX format.
# The legacy CTQMC functional below remains the path for its four-field CIX.
def ReadBoldCix(filename):
    with open(filename) as stream:
        next(stream); next(stream)
        header = next(stream).split()
        if len(header) == 4:
            return None
        if len(header) != 6:
            raise ValueError('Unrecognized CIX header in impurity functional')
        norbitals, nbaths = int(header[2]), int(header[3])
        next(stream)
        entries = []
        for ib in range(nbaths):
            row = next(stream).split('#', 1)[0].split()
            if int(row[0]) != ib:
                raise ValueError('CIX bath entries are out of order')
            first, second, channel = row[1:4]
            entries.append((int(first.rstrip('*')), int(second.rstrip('*')),
                            abs(int(channel.rstrip('*'))),
                            -1 if channel.startswith('-') else 1,
                            channel.endswith('*'),
                            not first.endswith('*') and second.endswith('*')))
        next(stream)
        levels = np.array([float(v) for v in next(stream).split('#', 1)[0].split()])
    if not entries or builtins.max(e[2] for e in entries) >= len(levels):
        raise ValueError('CIX component labels exceed its impurity levels')
    if builtins.any(a < 0 or b < 0 or a >= norbitals or b >= norbitals
           for a, b, *_ in entries):
        raise ValueError('CIX orbital index out of range')
    return norbitals, entries, levels


def BoldAtomicEnergy(layout, moments, atomic_energy):
    """Header moments already include directed-channel multiplicities."""
    levels = layout[2]
    moments = np.asarray(moments, dtype=float)
    if moments.shape != levels.shape or not np.isfinite(moments).all():
        raise ValueError('Atomic density components do not match CIX levels')
    return float(atomic_energy - np.dot(levels, moments))


def _bold_matrix(values, layout):
    n, entries, levels = layout
    values = np.asarray(values)
    if values.shape[-1] != len(levels) or not np.isfinite(values).all():
        raise ValueError('Matrix components do not match CIX or are not finite')
    result = np.zeros(values.shape[:-1] + (n, n), dtype=complex)
    seen = set()
    for a, b, c, sign, conjugate, normal in entries:
        if not normal:
            raise ValueError('Bold impurity functional requires normal baths; '
                             'anomalous baths need a separate Nambu functional')
        if (a, b) in seen:
            raise ValueError('Duplicate physical matrix entry in CIX')
        seen.add((a, b))
        result[..., a, b] = sign * (values[..., c].conj() if conjugate else values[..., c])
    if not builtins.all((i, i) in seen for i in range(n)):
        raise ValueError('Bold impurity functional requires all physical diagonal entries')
    return result


def _bold_read_table(filename, ncomponents):
    data = np.loadtxt(filename, ndmin=2)
    if data.shape[1] != 1 + 2*ncomponents or not np.isfinite(data).all():
        raise ValueError('Invalid impurity table: ' + str(filename))
    if len(data) < 8 or not np.all(np.diff(data[:, 0]) > 0):
        raise ValueError('Impurity frequencies must be increasing: ' + str(filename))
    return data[:, 0], data[:, 1::2] + 1j*data[:, 2::2]


def _bold_header(filename):
    with open(filename) as stream:
        header = next(stream)
    result = {}
    for token in header.lstrip('#').split():
        if '=' in token:
            key, value = token.split('=', 1)
            result[key] = ast.literal_eval(value)
    return result


def _bold_sum(values, omega, beta, fout, name):
    """Sum a real O(omega^-2) trace with a fitted analytic infinite tail.

    Fit the trace, rather than individual entries, to preserve basis invariance.
    Different fitting windows supply a tail-sensitivity diagnostic, not a
    statistical error bar. The input must contain every positive Matsubara point.
    """
    values = np.asarray(values, dtype=float)
    expected = (2*np.arange(len(omega)) + 1)*np.pi/beta
    if values.shape != omega.shape or not np.isfinite(values).all():
        raise ValueError('Invalid Matsubara summand: ' + name)
    if not np.allclose(omega, expected, rtol=2e-10, atol=1e-11):
        raise ValueError('Bold functional requires a complete positive Matsubara mesh')
    n = len(omega)
    if n < 32:
        raise ValueError('At least 32 Matsubara points are required for tail checks')
    scale = omega[-1]
    def tail(first):
        w = omega[first:]
        fit = np.column_stack(((scale/w)**2, (scale/w)**4))
        coef = np.linalg.lstsq(fit, values[first:], rcond=None)[0]
        return 2/beta*builtins.sum(coef[j]*scale**(2*j+2)*(beta/(2*np.pi))**(2*j+2)
                          *zeta(2*j+2, n+0.5) for j in range(2))
    correction = tail(builtins.max(n//2, n-3000))
    alternate = tail(builtins.max(3*n//4, n-1500))
    sensitivity = abs(correction-alternate)
    print(name, 'tail=', correction, 'tail_window_sensitivity=', sensitivity, file=fout)
    if sensitivity > 1e-3:
        print('WARNING:', name, 'tail window sensitivity exceeds 1 meV; '
              'increase the frequency cutoff', file=fout)
    result = 2/beta*np.sum(values) + correction
    if not np.isfinite(result):
        raise ValueError('Non-finite Matsubara trace: ' + name)
    return float(result)


def _bold_delta_on_mesh(omega, source_omega, delta):
    """Interpolate by frequency values and use moments beyond the input mesh."""
    target = np.array(omega, copy=True)
    for edge in (source_omega[0], source_omega[-1]):
        target[np.isclose(target, edge, rtol=1e-10, atol=1e-11)] = edge
    spline = interpolate.CubicSpline(source_omega, delta, axis=0, extrapolate=False)
    inside = (target >= source_omega[0]) & (target <= source_omega[-1])
    if np.any(target < source_omega[0]):
        raise ValueError('Hybridization does not cover the lowest Matsubara frequency')
    values = np.empty((len(target), delta.shape[1]), complex)
    derivative = np.empty_like(values)
    values[inside] = spline(target[inside])
    derivative[inside] = spline(target[inside], 1)
    if not inside.all():
        w = source_omega[builtins.max(len(source_omega)//2, len(source_omega)-1000):]
        d = delta[-len(w):]
        scale = w[-1]
        re = np.linalg.lstsq(np.column_stack(((scale/w)**2, (scale/w)**4)), d.real, rcond=None)[0]
        im = np.linalg.lstsq(np.column_stack((scale/w, (scale/w)**3)), d.imag, rcond=None)[0]
        x = target[~inside, None]; u = scale/x
        values[~inside] = u**2*re[0]+u**4*re[1] + 1j*(u*im[0]+u**3*im[1])
        derivative[~inside] = (-2*u**2*re[0]-4*u**4*re[1]
                               +1j*(-u*im[0]-3*u**3*im[1]))/x
    return values, derivative


def _bold_log_g(omega, green, reference_hamiltonian, beta, fout):
    hermitian = (reference_hamiltonian + reference_hamiltonian.conj().T)/2
    if not np.allclose(hermitian, reference_hamiltonian, rtol=1e-8, atol=1e-9):
        raise ValueError('High-frequency reference Hamiltonian is not Hermitian')
    eigenvalues = np.linalg.eigvalsh(hermitian)
    sign, logabs = np.linalg.slogdet(green)
    if np.any(sign == 0):
        raise ValueError('Singular impurity Green function in matrix trace-log')
    # Real part of log(det G/det Gref), independent of logarithm branches.
    relative = logabs + np.log(np.abs(1j*omega[:, None]-eigenvalues)).sum(axis=1)
    free_reference = -np.logaddexp(0., -beta*eigenvalues).sum()/beta
    return _bold_sum(relative, omega, beta, fout, 'Tr(log(-Gimp))') + free_reference


def _give_bold_functional(dire, names, layout, Vdc, fout):
    cix, delta_file, green_file, sigma_file, beta, mu = names
    ncomp = len(layout[2])
    omega, g = _bold_read_table(os.path.join(dire, green_file), ncomp)
    omega_s, sigma = _bold_read_table(os.path.join(dire, sigma_file), ncomp)
    omega_d, delta = _bold_read_table(os.path.join(dire, delta_file), ncomp)
    if omega.shape != omega_s.shape or not np.allclose(omega, omega_s, rtol=1e-10, atol=1e-11):
        raise ValueError('Green function and self-energy frequency meshes differ')
    meta = _bold_header(os.path.join(dire, green_file))
    beta, mu = float(beta), float(mu)
    if not np.isclose(beta*float(meta['T']), 1., rtol=1e-9):
        raise ValueError('PARAMS and Green-function temperatures disagree')
    if not np.isclose(mu, float(meta['mu']), rtol=1e-9, atol=1e-10):
        raise ValueError('PARAMS and Green-function chemical potentials disagree')
    moments = np.asarray(meta['mom'], float)
    multiplicity = np.bincount([e[2] for e in layout[1]], minlength=ncomp)
    if np.any(multiplicity == 0):
        raise ValueError('Unused CIX component in impurity functional')
    rho = _bold_matrix(moments/multiplicity, layout)
    if not np.allclose(rho, rho.conj().T, atol=1e-9):
        raise ValueError('Atomic density matrix is not Hermitian')
    density_trace = float(np.trace(rho).real)
    nd = float(meta.get('nd', density_trace))
    if not np.isclose(density_trace, nd, rtol=1e-9, atol=1e-10):
        # Older Bold output can have an inconsistent reduced moment header.
        # Preserve the direct pseudo-particle number for the mu*N term;
        # report the mismatch rather than silently renormalizing rho.
        print('WARNING: atomic density trace=', density_trace,
              'differs from pseudo-particle nd=', nd, file=fout)
    green = _bold_matrix(g, layout)
    sigma_matrix = _bold_matrix(sigma, layout)
    sigma_inf = _bold_matrix(sigma[-1].real, layout)
    levels = _bold_matrix(layout[2], layout)
    reference = levels-mu*np.eye(layout[0])+sigma_inf
    log_g = _bold_log_g(omega, green, reference, beta, fout)
    delta, derivative = _bold_delta_on_mesh(omega, omega_d, delta)
    delta_matrix = _bold_matrix(delta, layout)
    derivative_matrix = _bold_matrix(derivative, layout)
    trace_product = lambda a, b: np.einsum('nab,nba->n', a, b).real
    delta_g = _bold_sum(trace_product(delta_matrix, green), omega, beta, fout, 'Tr(Delta G)')
    cxx = _bold_sum(omega*trace_product(derivative_matrix, green), omega, beta, fout, 'Tr(omega dDelta/domega G)')
    sigma_g_dynamic = _bold_sum(trace_product(sigma_matrix-sigma_inf, green), omega, beta, fout, 'Tr((Sigma-Sigma_inf) G)')
    sigma_g_static = float(np.trace(sigma_inf@rho).real)
    sigma_g = sigma_g_dynamic+sigma_g_static
    atomic_energy = float(meta['Epot'])
    interaction = BoldAtomicEnergy(layout, moments, atomic_energy)
    one_body = float(np.dot(layout[2], moments))
    energy = delta_g+atomic_energy-cxx-mu*nd
    dc_matrix = _bold_matrix(np.asarray(Vdc), layout)
    dc_density = float(np.trace(dc_matrix@rho).real)
    zeorb = delta_g+atomic_energy-cxx-one_body
    # Preserve the established energy-based estimator combination. This is
    # Phi_DMFT + T*Simp; entropy is supplied separately when available.
    phi_energy = -log_g+energy+sigma_g
    print('Bold matrix functional: density=atomic; estimator combination=legacy energy-based', file=fout)
    print('nd=', nd, 'nf=', meta['nf'], 'density_difference=', meta['nf']-nd, file=fout)
    print('TrDeltaG=', delta_g, 'solver Ekin=', meta['Ekin'], file=fout)
    print('cxx=', cxx, 'lnGimp=', log_g, 'Eimp=Fimp+TSimp=', energy, file=fout)
    print('Atomic interaction energy=', interaction, 'header TrSigmaG=', meta['TrSigmaG'], file=fout)
    print('1/2*Tr(Sigma*G)=', sigma_g/2, '(diagnostic; not substituted for atomic energy)', file=fout)
    print('Vdc*nd=', dc_density, 'zeorb[eV]=', zeorb, file=fout)
    print('Phi_DMFT+TSimp=', phi_energy, file=fout)
    if 'Fimp' in meta:
        print('solver Fimp=', meta['Fimp'], '(diagnostic; not substituted)', file=fout)
    result = (phi_energy, log_g, energy, dc_density, zeorb)
    if not np.isfinite(result).all():
        raise ValueError('Non-finite Bold impurity functional')
    return result


def GiveImpFunctional(dire, fparams, Vdc, fout=sys.stdout):

    if not (os.path.exists(dire+'/'+fparams) and os.path.getsize(dire+'/'+fparams)>0):
        print('Give input PARAMS file', file=fout)
        sys.exit(0)
    
    names = ReadPARAMS2(dire+'/'+fparams, fout)
    (cix,Delta,Gf,Sig,beta,mu_qmc) = names
    layout = ReadBoldCix(os.path.join(dire, cix))
    if layout is not None:
        return _give_bold_functional(dire, names, layout, Vdc, fout)
    (Deg, Eks) = ReadCix(dire+'/'+cix, fout)
    
    
    (om, omD, Gfc, Sigc, Dlta, nf, mu, Ekin, Epot, TrSigmaG, mom) = ReadGandSigandDelta(dire+'/'+Gf, dire+'/'+Sig, dire+'/'+Delta, fout)
    s_oo = Sigc[:,-1].real
    print('s_oo=', s_oo, file=fout)
    
    lengh = builtins.min(len(Eks),len(s_oo))
    EimpS = Eks[:lengh]-mu_qmc+s_oo[:lengh]
    # Tr(log(-Gimp))
    lnGimp = LogGimp(om,Gfc,EimpS,beta,Deg,fout) 

    #lnGS = LogGimpSig(om, Gfc, Sigc, Vdc, EimpS, beta, Deg, fout)

    trsigg =  fTrSigmaG(om, Gfc, Sigc, Deg, EimpS, beta, fout) # numericall Tr((Sigma-s_oo)*G)
    trsigg_oo = sum(mom[:]*s_oo[:])

    # Tr(Delta*G)
    TrDeltaG = fTrDeltaG(om, omD, Dlta, Gfc, Deg, EimpS, beta, fout)
    print('TrDeltaG=', TrDeltaG, '<k>/T=', Ekin, file=fout)
    # Tr(om*dDelta/dom)
    cxx = CmpImpurityCorrection2(om, omD, Dlta, Gfc, Deg, EimpS, beta, fout)

    # Eimp = Tr(Delta*G)+1/2*Tr(Sigma*G)+eimp*nimp-Tr(om*dDelta/dom)
    #Eimp = Ekin+Epot-cxx-mu*nf                # Impurity Energy
    Eimp = TrDeltaG+Epot-cxx-mu*nf                # Impurity Energy

    # Potthof2[G] == Phi[G]-1/2*Tr(Sigma*G)+T*S_imp = Eimp-Tr(log(-G_imp))+1/2*Tr(Sigma*G)
    #Potthof2 = Eimp-lnGimp+TrSigmaG            # Phi[G]-1/2*Tr(Sigma*G)+T*S_imp

    #nd_numeric=zeros(len(Deg))
    #for i in range(len(Deg)):
    #    Gd = Gfc[i,:]
    #    eimps = EimpS[i]
    #    ff = Gd[:]-1./(om*1j-eimps)
    #    nd_numeric[i] = (2*sum(ff.real)/beta + ferm(eimps*beta))*Deg[i]
    #print 'nd_numeric=', nd_numeric
    
    #Eimp1 = lnGS - lnGimp + Eimp  # Tr(log(1+Sigma*G))-Tr(log(-Gloc))+E_imp

    logZatom=0
    filename=dire+'/ctqmc.log'
    if os.path.isfile(filename):
        fw=open(filename,'r')
        lines=fw.readlines()
        for line in lines:
            m = re.search('log\(Zatom\)=(\d+\.\d+)',line)
            if m is not None:
                logZatom = float(m.groups()[0])
        fw.close()
        #print 'logZatom=', logZatom
    P0=0
    filename=dire+'/histogram.dat'
    if os.path.isfile(filename):
        fw=open(filename,'r')
        next(fw)# comment
        P0 = float(next(fw).split()[1])
        #print 'P0=', P0

    Eimp2 = TrSigmaG   # 1/2*Tr(Sigma*G)
    print('lnGimp=', lnGimp, file=fout)
    print("%-12s "*6 % ('Eimp', 'Ekin', 'Epot', 'cxx', 'mu*nf', 'lnGimp'), file=fout)
    print("%-12.6f "*6 % (Eimp, TrDeltaG, Epot, cxx, mu*nf, lnGimp), file=fout)
    print('#',"%-12.6f "*6 % (Eimp+Ekin-TrDeltaG, Ekin, Epot, cxx, mu*nf, lnGimp), file=fout)
    print('Eimp=', Eimp, file=fout)
    if logZatom and P0>1e-7:
        logZ = logZatom-log(P0)
        Fimp = -logZ/beta
        Entropy = (Eimp-Fimp)*beta
        WARN=''
        if (P0<1e-5): WARN='   ** WARNING P0 is very small, hence Entropy is not accurate'
        print('Entropy=', Entropy, 'Fimp=', Fimp, 'logZatom=', logZatom, 'P0=', P0, WARN, file=fout)
    #print >> fout, 'P[Sigma]=', Potthof2

    print('1/2*Tr(Sigma*G)=', 0.5*(trsigg+trsigg_oo), 'deimp*nf=', sum(Eks*mom), file=fout)

    Vdc_nd = sum(array(mom[:])*array(Vdc[:]))
    zeorb = TrDeltaG+Epot-cxx - sum(Eks*mom) # note that Epot ~ 1/2*Tr(Sigma*G)+Eks*mom
    print('Vdc*nf=', Vdc_nd, file=fout)
    Phi_DMFT = -lnGimp+Eimp+(trsigg+trsigg_oo)
    print('Phi_DMFT=', Phi_DMFT, file=fout)
    print('Phi_DMFT - Vdc*nf=', Phi_DMFT-Vdc_nd, file=fout)
    Ry2eV = 13.60569193
    print('Phi_DMFT - Vdc*nf=', (Phi_DMFT-Vdc_nd)/Ry2eV, 'Ry', file=fout)
    print('Epot - Vdc*nf=', (TrSigmaG-Vdc_nd)/Ry2eV, 'Ry', file=fout)
    print('zeorb[eV]=', zeorb, file=fout)
    return (Phi_DMFT, lnGimp, Eimp, Vdc_nd, zeorb)


if __name__ == '__main__':
    import re
    
    if len(sys.argv)<2:
        dire = '.'
    else:
        dire = sys.argv[1]
    Ry2eV = 13.60569193

    fi = open(dire+'/Eimp.inp','r')
    line = next(fi)  # Eimp
    line = next(fi)  # Overlap
    line = next(fi)  # Vdc
    Vdc = list(map(float,re.sub(r'#.*',' ',line).split()))

    (Phi_DMFT, lnGimp, Fimp, Vdc_nd, zeorb) = GiveImpFunctional(dire, 'PARAMS', Vdc)
    print('ev: Phi_DMFT-Vdc*nd=', Phi_DMFT-Vdc_nd, 'Fimp+TS=', Fimp, 'logGimp=', lnGimp, 'zeorb=', zeorb)
    print('Ry: Phi_DMFT-Vdc*nd=', (Phi_DMFT-Vdc_nd)/Ry2eV, 'FImp+TS=', Fimp/Ry2eV, 'logGimp=', lnGimp/Ry2eV, 'zeorb=', zeorb/Ry2eV)
    
