#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
# 
import sys
#from scipy import *
#from numpy import *
import numpy as np
import math
from math import pi, tan, sqrt
from numpy import random
from scipy import interpolate
from scipy import integrate
from scipy import optimize
import maxent_routines as me

def KramarsKronig_inside(fi, om, Dh, x, i0, S0):
    """
    A simplified version of the Kramers–Kronig calculation for real-axis spectral function.
    We can get the real part, or the complex Matsubara Green's function if x is purely imaginary.
    Parameters
    ----------
    fi : 1D NumPy array
        Imaginary part of the spectral function at each mesh point.
    om : 1D NumPy array
        The mesh points (frequencies or energies).
    Dh : 1D NumPy array
        Trapezoid integration weights for each om[i].
    x  : float or complex
        The point at which we compute the real part or Matsubara function.
    i0 : int
        Index where x == om[i0] or x is very close to om[i0]
    S0 : float
        Value S0=fi[i0] of more generally S0=fi(om[i0]) (the imaginary part at x) used for subtraction.
    Returns
    -------
    float or complex
        The real part of the function at x, computed via a Kramers–Kronig relation.
    """
    N  = om.size
    # sum over all points except om[i0]
    sum_val  = np.sum( (fi[:i0]-S0)*Dh[:i0]/(om[:i0]-x) )
    sum_val += np.sum( (fi[i0+1:]-S0)*Dh[i0+1:]/(om[i0+1:]-x) )
    # correction for point om[i0], for which x==om[i0] and it would give NaN
    if i0 == 0:
        sum_val += Dh[i0]*(fi[i0+1]-S0)/(om[i0+1]-x)
    elif i0==N-1:
        sum_val += Dh[i0]*(fi[i0-1]-S0)/(om[i0-1]-x)
    else:
        sum_val += Dh[i0]*0.5*( (fi[i0-1]-S0)/(om[i0-1]-x) + (fi[i0+1]-S0)/(om[i0+1]-x) )
    # Log term for the subtracted portion, valid if x is strictly inside the domain
    if x!=om[-1] and x!=om[0]:
        sum_val += S0 * math.log(abs((om[-1]-x)/(x-om[0])))
    return sum_val/math.pi

def KramarsKronig(om, fi):
    """Kramers–Kronig calculation for real-axis spectral function, it return the real part.
    Parameters
    ----------
    fi : 1D NumPy array
        Imaginary part of the spectral function at each mesh point.
    om : 1D NumPy array
        The mesh points (frequencies or energies). It should not contain frequency om=0.
    Returns:
    real_part : 1D Numpy array.
    """
    def GiveDh(om):
        Dh = np.empty_like(om)  # Preallocate output array of the same shape and type as om
        # Compute edge values separately
        Dh[0] = 0.5 * (om[1]-om[0])
        Dh[-1] = 0.5 * (om[-1]-om[-2])
        # Vectorized computation for inner elements
        Dh[1:-1] = 0.5 * (om[2:]-om[:-2])
        return Dh
    Dh = GiveDh(om)
    real_part = np.empty_like(om)
    for i0,x in enumerate(om):
        S0 = fi[i0] if not(i0==0 or i0==len(om)-1) else 0.0 # What we subtract (log) would diverge. Hence avoid it here.
        # at the two end-points the imaginary part should vanish, and this should be irrelevant. But when it does not, we want to avoid numeric problem.
        real_part[i0] = KramarsKronig_inside(fi, om, Dh, om[i0], i0, S0)
    return real_part


def Broad(width, om, fw):
    " Broadens the data with gaussian of width=width"
    def MakeTanMesh(N, tanc, tanw, b0, b1):
        if not(b0<b1): print("Relation must hold: b0<b1!")
        if not(b0<tanw and tanw<b1): print("Relation mesu hold: b0<tanw<b1!")
        if not(b0>0): print("b0 must be positive!")
        du = np.arctan(((tanc-b0)/tanw))
        b1n = np.arctan((b1-tanc)/tanw)+du
        m0 = [tanc + tanw * tan(b1n*(i-1)/(N-2)-du) for i in range(1,N)]
        return np.hstack( (-np.array(m0[::-1]), np.array([0]+m0) ) )

    if width<1e-5: return fw
    
    w=width
    x = MakeTanMesh(200,0.0,w,w/50,w*20)
    fwi = interpolate.interp1d(om, fw)
    fwn=[]
    for im in range(len(om)):
        x1 = [t for t in om if t>=om[im]-x[-1] and t<=om[im]-x[0]]
        x2 = [t for t in x+om[im] if t>=om[0] and t<=om[-1]]
        eps = sorted(np.hstack((x1, x2)))
        x3 = om[im]-eps
        gs = np.exp(-x3**2/(2*w**2))/(sqrt(2*pi)*w)
        yn = integrate.trapezoid(fwi(eps) * gs, x=eps)
        fwn.append(yn)
    return np.array(fwn)

def GiveTanMesh(x0,L,Nw):
    def fun(x,x0,L,Nw):
        "x[0]=d, x[1]=w"
        d=x[0]
        w=x[1]
        #print 'd=', d, 'w=', w
        return np.array([L-w/tan(d), x0-w*tan(pi/(2*Nw)-d/Nw) ])
    
    xi=x0/L
    d0 = Nw/2.*(tan(pi/(2*Nw))-sqrt(tan(pi/(2*Nw))**2 - 4*xi/Nw))
    w0 = L*d0

    sol=optimize.root(fun, [d0,w0], args=(x0,L,Nw) )
    (d,w) = sol.x
    om = w*np.tan(np.linspace(0,1,2*Nw+1)*(pi-2*d) -pi/2+d)
    return om

def MaximumEntropy(p, tau, Gt, sb='', log=sys.stdout):
    beta = tau[-1]
    random.seed( 1 ) # seed for random numbers

    if 'x0' in p:
        omega = GiveTanMesh(p['x0'],p['L'],p['Nw'])
    else:
        omega = np.linspace(-p['L'],p['L'],2*p['Nw']+1)
    dom = np.array([0.5*(omega[1]-omega[0])]+[0.5*(omega[i+1]-omega[i-1]) for i in range(1,len(omega)-1)]+[0.5*(omega[-1]-omega[-2])])

    fsg=1
    if p['statistics']=='fermi':
        Gt = -Gt
        fsg=-1
        normalization = Gt[0]+Gt[-1]
        Ker = me.initker_fermion(omega,dom,beta,tau)        
    elif p['statistics']=='bose':
        normalization = integrate.trapezoid(Gt,x=tau)
        Ker = me.initker_boson(omega,dom,beta,tau)
    
    print('beta=', beta, sb,file=log)
    print('normalization'+sb+'=', normalization, sb, file=log)

    # Set error
    if p['idg']:
        sxt = np.ones(len(tau))/(p['deltag']**2)
    else:
        sxt = Gt*p['deltag']
        for i in range(len(sxt)):
            if sxt[i]<1e-5: sxt[i]=1e-5
        sxt = 1./sxt**2
    
    # Set model
    if p['iflat']==0:
        model = normalization*np.ones(len(omega))/sum(dom)
    elif p['iflat']==1:
        model = np.exp(-omega**2/p['gwidth'])
        model *= normalization/(model @ dom)
    else:
        dat=np.loadtxt('model.dat').transpose()
        fm=interpolate.interp1d(dat[0],dat[1])
        model = fm(omega)
        model *= normalization/(model @ dom)
        #np.savetxt('brisi_test', np.vstack((tau, fsg*dot(model,Ker))).transpose())

        
    print('Model normalization'+sb+'=', (model@dom), sb, file=log)

    # Set starting Aw(omega)
    Aw = random.rand(len(omega))
    Aw = Aw * (normalization/(Aw@dom))
    print('Aw normalization'+sb+'=', (Aw@dom), sb, file=log)

    dlda = me.initdlda(omega,dom,Ker,sxt)
    
    temp=10.
    rfac=1.
    alpha=p['alpha0']
    
    for itt in range(p['Nitt']):
        print(str(itt)+sb, 'Restarting maxent with rfac=', rfac, 'alpha=', alpha, file=log)
        iseed = random.randint(0,2**16-1)
        me.maxent(Aw,rfac,alpha,temp,Ker,sxt,Gt,model,dom,p['Asteps'],iseed)
        S = me.entropy(Aw,model,dom)
        Trc = me.lambdac(alpha,Aw,omega,dom,dlda)
        chi2 = sum(sxt*(Gt-Aw@Ker)**2)
        
        ratio = -2*S*alpha/Trc
        #print('Finished maxent with chi2=', chi2, 'alpha=',alpha,'-2*alpha*S=',-2*alpha*S,'Trace=',Trc, 'S=', S)
        print(str(itt)+sb,'Finished maxent with chi2={:.3f} alpha={:.5f} -2*alpha*S={:.5f} Trace={:.5f} S={:.5f}'.format(chi2,alpha,-2*alpha*S,Trc,S), file=log)
        print(str(itt)+sb,'   ratio=', ratio, file=log)
        
        np.savetxt('dos_'+str(itt)+sb, np.vstack((omega,Aw)).transpose())
        temp=0.001
        rfac=0.05
    
        if abs(ratio-1)<p['min_ratio'] and chi2<sqrt(len(Gt)): break
    
        if (abs(ratio)<0.05):
            alpha *= 0.5
        else:
            alpha *= (1.+0.001*(random.rand()-0.5))/ratio
        
    for itt in range(p['Nr']):
        print(sb,'Smoothing itt ', itt, file=log)
        Aw = Broad(p['bwdth'],omega,Aw)
        Aw *= (normalization/(Aw@dom)) # Normalizing Aw
        
        np.savetxt('dos_'+str(p['Nitt'])+sb, np.vstack((omega,Aw)).transpose())
        
        temp=0.005
        rfac=0.005
        iseed = random.randint(0,2**31-1)
        me.maxent(Aw,rfac,alpha,temp,Ker,sxt,Gt,model,dom,p['Asteps'],iseed)
        
        S = me.entropy(Aw,model,dom)
        Trc = me.lambdac(alpha,Aw,omega,dom,dlda)
        ratio = -2*S*alpha/Trc
        print(sb,'Finished smoothing run with alpha=',alpha,'-2*alpha*S=',-2*alpha*S,'Trace=',Trc, file=log)
        print(sb,'   ratio=', ratio, file=log)
        
    np.savetxt('gtn'+sb, np.vstack((tau, fsg*(Aw@Ker))).transpose())
    Aw = Broad(p['bwdth'],omega,Aw)
    np.savetxt('dos.out'+sb, np.vstack((omega,Aw)).transpose())
    return (Aw, omega)

def Pade(om, Gm, x, gamma, Norder):
    zn = om[:Norder]*1j
    gn = Gm[:Norder]
    Pt = me.padecof(gn, zn)
    Gx = np.array([me.padeg(w+gamma*1j, zn, Pt) for w in x])
    return Gx

def InverseFourier(Gm, om, tau, beta, Nf=40, stat='fermi'):
    """Inverse Fourier transform which
       computes G(tau) from G(iom)
    """
    def FindHighFrequency(Gm,om,Nf):
        S=0.; Sx=0.; Sy=0.; Sxx=0.; Sxy=0;
        for j in range(len(om)-Nf,len(om)):
            x = om[j]
            y = Gm[j].imag * x
            x2= x**2
            Sy += y
            Sx += 1/x2
            Sxx += 1/x2**2
            Sxy += y/x2
            S += 1
    
        dd = S*Sxx-Sx**2
        a = (Sxx*Sy-Sx*Sxy)/dd
        bx = (S*Sxy-Sx*Sy)/dd
        ah = -a;
        if abs(ah-1.0)<1e-3: ah=1.0
        return ah
    
    Gtau = np.zeros(len(tau),dtype=float)
    # Correction 1/omega^2 (sometimes usefull)
    df = Gm[-1].real*om[-1]/pi
    print('df=', df)
    if stat=='fermi':
        ah = FindHighFrequency(Gm,om,Nf)
        for it,t in enumerate(tau):
            Gtau[it] = me.fourpart(t,Gm,om,ah,beta)
        Gtau[0] += df
        Gtau[-1] -= df
    else:
        #ah=0
        #om[0]=1e-12
        for it,t in enumerate(tau):
            Gtau[it] = me.fourpartb(t,Gm,om,beta)
            #Gtau[it] = me.fourpart(t,Gm,om,ah,beta)
        #om[0]=0.
        Gtau[0] += df
        Gtau[-1] += df
        
    return Gtau


####################################
# Limited use -- only for rear cases
####################################
def MaximumEntropyTest(p, tau, Gt):
    
    def MEStep(alpha,rfac,Aw,temp,Ker,sxt,Gt,model,f0,Asteps,itt,reset=True):
        if (reset):
            temp=0.001
            rfac=0.05
        print('Restarting maxent with rfac=', rfac, 'alpha=', alpha)
        iseed = random.randint(0,2**31-1)
        me.maxent(Aw,rfac,alpha,temp,Ker,sxt,Gt,model,f0,Asteps,iseed)
        S = me.entropy(Aw,model,f0)
        Trc = me.lambdac(alpha,Aw,omega,dlda)
        ratio = -2*S*alpha/Trc
        print('Finished maxent with alpha=',alpha,'-2*alpha*S=',-2*alpha*S,'Trace=',Trc)
        print('   ratio=', ratio)
        np.savetxt('dos_'+str(itt), np.vstack((omega,Aw)).transpose())
        return ratio
    
    beta = tau[-1]

    random.seed( 1 ) # seed for random numbers
    
    omega = np.linspace(-p['L'],p['L'],2*p['Nw']+1)
    f0,f1,f2 = me.initf0(omega)
    fsg=1
    if p['statistics']=='fermi':
        Gt = -Gt
        fsg=-1
        normalization = Gt[0]+Gt[-1]
        Ker = me.initker_fermion(omega,beta,tau)        
    elif p['statistics']=='bose':
        normalization = integrate.trapezoid(Gt,x=tau)
        dom = np.array([0.5*(omega[1]-omega[0])]+[0.5*(omega[i+1]-omega[i-1]) for i in range(1,len(omega)-1)]+[0.5*(omega[-1]-omega[-2])])
        Ker = me.initker_boson(omega,dom,beta,tau)
    
    print('beta=', beta)
    print('normalization=', normalization)

    # Set error
    if p['idg']:
        sxt = np.ones(len(tau))/(p['deltag']**2)
    else:
        sxt = Gt*p['deltag']
        for i in range(len(sxt)):
            if sxt[i]<1e-5: sxt[i]=1e-5
        sxt = 1./sxt**2
    
    # Set model
    if p['iflat']==0:
        model = normalization*np.ones(len(omega))/sum(f0)
    elif p['iflat']==1:
        model = exp(-omega**2/p['gwidth'])
        model *= normalization/(model@f0)
    else:
        dat=np.loadtxt('model.dat').transpose()
        fm=interpolate.interp1d(dat[0],dat[1])
        model = fm(omega)
        model *= normalization/(model@f0)
        #np.savetxt('brisi_test', np.vstack((tau, fsg*dot(model,Ker))).transpose())

        
    print('Model normalization=', (model@f0))

    # Set starting Aw(omega)
    Aw = random.rand(len(omega))
    Aw = Aw * (normalization/(Aw@f0))
    print('Aw normalization=', (Aw@f0))

    dlda = me.initdlda(omega,Ker,sxt)
    
    temp=10.
    rfac=1.
    alpha=p['alpha0']

    for itt in range(10):
        ratio = MEStep(alpha,rfac,Aw,temp,Ker,sxt,Gt,model,f0,p['Asteps'],itt,itt!=0)
        if abs(ratio-1)<p['min_ratio']: break
        if (ratio<0.05):
            if ratio>0:
                alpha *= 1.1
            else:
                alpha /= 10.
        else:
            alpha *= (1.+0.001*(random.rand()-0.5))/ratio
            
    if abs(ratio-1)>2*p['min_ratio']:
        alpha=1.
        p['Asteps'] *= 1.5
        for itt in range(3,6):
            ratio = MEStep(alpha,rfac,Aw,temp,Ker,sxt,Gt,model,f0,p['Asteps'],itt)
    
    for itt in range(p['Nr']):
        print('Smoothing itt ', itt)
        Aw = Broad(p['bwdth'],omega,Aw)
        Aw *= (normalization/(Aw@f0)) # Normalizing Aw
        ratio = MEStep(alpha,rfac,Aw,temp,Ker,sxt,Gt,model,f0,p['Asteps'],p['Nitt']+itt)
        
        
    np.savetxt('gtn', np.vstack((tau, fsg*(Aw@Ker))).transpose())
    Aw = Broad(p['bwdth'],omega,Aw)
    np.savetxt('dos.out', np.vstack((omega,Aw)).transpose())
    return (Aw, omega)

def PadeTest(om, Gm, x, gamma, Norder):
    from scipy import poly1d
    from pylab import plot, show

    zn = om[:Norder]*1j
    gn = Gm[:Norder]
    an = me.padecof(gn, zn)

    print('zn=', zn)
    print('an=', an)
    
    Aq_0 = poly1d([0])
    Aq_1 = poly1d([an[0]])
    Bq_0 = poly1d([1.])
    Bq_1 = poly1d([1.])
    for i in range(Norder-1):
        Aq_2 = Aq_1 +poly1d([1,-zn[i]])*an[i+1]*Aq_0
        Aq_0 = Aq_1
        Aq_1 = Aq_2
        Bq_2 = Bq_1 +poly1d([1,-zn[i]])*an[i+1]*Bq_0
        Bq_0 = Bq_1
        Bq_1 = Bq_2
        
    poles = sorted( roots(Bq_2), key=real )
    ezeros = sorted( roots(Aq_2), key=real )
    Bqp = poly1d(poles, r=True)
    Cbq = Bq_2(0.0)/Bqp(0.0)
    print('ratio=', Cbq)
    
    wgh=[]
    print('poles=')
    for i in range(len(poles)):
        Rq = poly1d(poles[:i] + poles[i+1:], r=True)
        wi = Aq_2(poles[i])/(Cbq*Rq(poles[i]))
        wgh.append(wi)
        if (poles[i].imag<0): used='Yes'
        else: used='No'
        print("%2d %12.4g %12.4g   %12.4g %12.4g  used=%s" % (i+1,poles[i].real, poles[i].imag, wi.real, wi.imag, used))
    wgh=np.array(wgh)
    #print 'np.zeros='
    #for i in range(len(ezeros)):
    #    print i+1, ezeros[i]
    #print 'weights=', wgh

    yt = np.zeros(len(om), dtype=complex)
    for i in range(len(poles)):
        if (poles[i].imag<0):
            yt += real(wgh[i])/(om*1j-poles[i])
    
    normy = sum(real(yt))
    normg = sum(real(Gm))
    if normy>1e-6:
        print('Warning: norm mismatch: ratio=', normg/normy)
        print('Renormalizing')
        wgh *= (normg/normy)
    else:
        print('Warning: Not enough poles. Bailing out')
    
    yt = np.zeros(len(om), dtype=complex)
    for i in range(len(poles)):
        if (poles[i].imag<0):
            yt += real(wgh[i])/(om*1j-poles[i])
    yr = np.zeros(len(x), dtype=complex)
    for i in range(len(poles)):
        if (poles[i].imag<0):
            yr += real(wgh[i])/(x-poles[i]+gamma*1j)
    
    np.savetxt('pade.corrected', np.vstack((x,real(yr),imag(yr))).transpose())
    G0 = Aq_2(zn)/Bq_2(zn)
    print('G0=', G0)
    plot(om, real(Gm), 'o')
    plot(om, real(yt), ':')
    show()
    plot(x, imag(yr))
    #xlim([0,2.4])
    #ylim([0,0.6])
    show()
    
    
    Gx = np.array([me.padeg(w+gamma*1j, zn, an) for w in x])
    return Gx
    
def test1a(filename,stat,deltag=0.006,L=4):
    params={'statistics': stat,    # fermi/bose                                                                                     
            'L'         : L,       # cutoff frequency on real axis                                                                  
            'gwidth'    : 2*L,     # width of gaussian                                                                              
            'Nw'        : 201,     # number of frequencies                                                                          
            'idg'       : 1,       # error scheme: idg=1 -> sigma=deltag ; idg=0 -> sigma=deltag*G(tau)                             
            'deltag'    : deltag,  # error                                                                                          
            'Asteps'    : 4000,    # anealing steps                                                                                 
            'alpha0'    : 1000,    # starting alpha                                                                                 
            'min_ratio' : 0.001,    # condition to finish, what should be the ratio                                                 
            'iflat'     : 1,       # iflat=0 : constant model, iflat=1 : gaussian of width gwidth, iflat=2 : input using file model$
            'Nitt'      : 1000,    # maximum number of outside iterations                                                           
            'Nr'        : 0,       # number of smoothing runs                                                                       
            'bwdth'     : 0.03,    # smoothing width                                                                                
            }

    Gtdat = np.loadtxt(filename).transpose()
    tau = Gtdat[0]
    Gt  = Gtdat[1]
    (Aw, omega) = MaximumEntropy(params, tau, Gt)

def test1(om,Gm,stat,deltag=0.006,Ntau=400,L=4):
    params={'statistics': stat,    # fermi/bose
            'L'         : L,       # cutoff frequency on real axis
            'gwidth'    : 2*L,     # width of gaussian
            'Nw'        : 301,     # number of frequencies
            'idg'       : 1,       # error scheme: idg=1 -> sigma=deltag ; idg=0 -> sigma=deltag*G(tau)
            'deltag'    : deltag,  # error
            'Asteps'    : 4000,    # anealing steps
            'alpha0'    : 1000,    # starting alpha
            'min_ratio' : 0.001,    # condition to finish, what should be the ratio
            'iflat'     : 1,       # iflat=0 : constant model, iflat=1 : gaussian of width gwidth, iflat=2 : input using file model.dat
            'Nitt'      : 1000,    # maximum number of outside iterations
            'Nr'        : 0,       # number of smoothing runs
            'bwdth'     : 0.03,    # smoothing width
            }
    

    if stat=='bose':
        beta = 2*pi/om[1]
    else:
        beta = pi/om[0]
    tau = np.linspace(0,beta,Ntau+1)
    Gt = InverseFourier(Gm, om, tau, beta,  40, stat)
    
    np.savetxt('Gtau.dat', np.vstack((tau,Gt)).transpose() )
    (Aw, omega) = MaximumEntropy(params, tau, Gt)
    return (Aw,omega)

def test3(filename='Giom.dat',stat='fermi',deltag=0.006,Ntau=200,L=20.):
    params={'statistics': stat,    # fermi/bose
            'L'         : L,       # cutoff frequency on real axis
            'gwidth'    : 2*L,     # width of gaussian
            'Nw'        : 200,     # number of frequencies
            'idg'       : 1,       # error scheme: idg=1 -> sigma=deltag ; idg=0 -> sigma=deltag*G(tau)
            'deltag'    : deltag,   # error
            'Asteps'    : 4000,    # anealing steps
            'alpha0'    : 1000,    # starting alpha
            'min_ratio' : 0.001,   # condition to finish, what should be the ratio
            'iflat'     : 1,       # iflat=0 : constant model, iflat=1 : gaussian of width gwidth, iflat=2 : input using file model.dat
            'Nitt'      : 1000,    # maximum number of outside iterations
            'Nr'        : 0,       # number of smoothing runs
            'bwdth'     : 0.03,    # smoothing width
            }

    Gdat = np.loadtxt(filename).transpose()
    om = Gdat[0]
    Gm = Gdat[1]+Gdat[2]*1j
    beta = pi/om[0]
    tau = np.linspace(0,beta,Ntau+1)
    
    Gt = InverseFourier(Gm, om, tau, beta, Nf=40)
    (Aw, omega) = MaximumEntropy(params, tau, Gt)

    
def test5(filename, Norder, L=6.):
    data = np.loadtxt(filename).transpose()
    tau = data[0]
    Gt = data[1]
    beta = tau[-1]
    x = np.linspace(-L,L,801)
    gamma = 0.005
    
    fGt = interpolate.UnivariateSpline(tau, Gt, s=0)
    tau2 = np.linspace(0,beta,len(tau)*10)
    fGt2 = fGt(tau2)
    np.savetxt('Gt_interpolated', np.vstack((tau2,fGt2)).transpose())
    Gm=[]
    om=[]
    for im in range(Norder):
        omi = 2*im*pi/beta
        om.append(omi)
        Gm.append( integrate.trapezoid( fGt2*cos(tau2*omi), x=tau2 ) + 0j )
    Gm = np.array(Gm)
    om = np.array(om)

    np.savetxt('Gm', np.vstack((om,real(Gm))).transpose())

    Gx = Pade(om, Gm, x, gamma, Norder)
    np.savetxt('dos.pade', np.vstack((x,real(Gx),imag(Gx))).transpose())

    Ax = np.zeros(len(x),dtype=float)
    for ix,xx in enumerate(x):
        if (abs(xx)>1e-12):
            Ax[ix] = imag(Gx[ix])/(pi*xx)
        else:
            Ax[ix] = 0.5*imag(Gx[ix-1])/(pi*x[ix-1]) + 0.5*imag(Gx[ix+1])/(pi*x[ix+1]) 
    np.savetxt('Ax', np.vstack((x,Ax)).transpose())
    omega = x
    dom = np.array([0.5*(omega[1]-omega[0])]+[0.5*(omega[i+1]-omega[i-1]) for i in range(1,len(omega)-1)]+[0.5*(omega[-1]-omega[-2])])
    Ker = me.initker_boson(omega,dom,beta,tau)
    np.savetxt('gtn', np.vstack((tau, (Ax@Ker))).transpose())

def test5b(filename, Norder, L=6.):
    data = np.loadtxt(filename).transpose()
    tau = data[0]
    Gt = data[1]
    beta = tau[-1]
    x = np.linspace(-L,L,801)
    gamma = 0.005
    
    fGt = interpolate.UnivariateSpline(tau, Gt, s=0)
    tau2 = np.linspace(0,beta,len(tau)*10)
    fGt2 = fGt(tau2)
    np.savetxt('Gt_interpolated', np.vstack((tau2,fGt2)).transpose())
    Gm=[]
    om=[]
    for im in range(Norder*3):
        omi = 2*im*pi/beta
        om.append(omi)
        Gm.append( integrate.trapezoid( fGt2*cos(tau2*omi), x=tau2 ) + 0j )
    Gm = np.array(Gm)
    om = np.array(om)

    np.savetxt('Gm', np.vstack((om,real(Gm))).transpose())

    Gx = PadeTest(om, Gm, x, gamma, Norder)
    np.savetxt('dos.pade', np.vstack((x,real(Gx),imag(Gx))).transpose())

    Ax = np.zeros(len(x),dtype=float)
    for ix,xx in enumerate(x):
        if (abs(xx)>1e-12):
            Ax[ix] = imag(Gx[ix])/(pi*xx)
        else:
            Ax[ix] = 0.5*imag(Gx[ix-1])/(pi*x[ix-1]) + 0.5*imag(Gx[ix+1])/(pi*x[ix+1]) 
    np.savetxt('Ax', np.vstack((x,Ax)).transpose())
    omega = x
    dom = np.array([0.5*(omega[1]-omega[0])]+[0.5*(omega[i+1]-omega[i-1]) for i in range(1,len(omega)-1)]+[0.5*(omega[-1]-omega[-2])])
    Ker = me.initker_boson(omega,dom,beta,tau) # (ker,w,dw,beta,t,nt,nw)
    np.savetxt('gtn', np.vstack((tau, (Ax@Ker))).transpose())

def test6(filename,Norder=100,L=10):
    
    Sdat = np.loadtxt(filename).transpose()
    om = Sdat[0]
    Sm = Sdat[1::2]+Sdat[2::2]*1j
    beta = pi/om[0]
    print('beta=', beta)

    for i in range(len(Sm)):
        #Gm = 1/(om*1j-Sm[i])
        Gm = Sm[i]
        x = np.linspace(-L,L,501)
        gamma = 0.001
    
        for norder in range(10,Norder):
            Gx = Pade(om, Gm, x, gamma, norder)
            #Sx = x-1/Gx
            Sx = Gx
            np.savetxt('dos.pade.'+str(i)+'.Norder'+str(norder), np.vstack((x,real(Sx),imag(Sx))).transpose())


def test6b(fin,fout,Norder=100,L=10,Nmax=150): #filename,Norder=100,L=10):
    
    Gdat = np.loadtxt(fin).transpose()
    om = Gdat[0,:Nmax]
    Gm = Gdat[1,:Nmax]+Gdat[2,:Nmax]*1j
    beta = pi/om[0]
    
    x = np.linspace(-L,L,501)
    gamma = 0.01
    
    Gx = PadeTest(om, Gm, x, gamma, Norder)
    np.savetxt(fout, np.vstack((x,real(Gx),imag(Gx))).transpose())


def test2(fin, fout, Norder):
    # for the one band model on Bethe lattice
    def Gz(z):
        return 2*(z - sign(z.real+1e-16)*sqrt(z**2-1.))
    
    x = np.linspace(-6,6,1001)
    gamma = 0.001
    data = np.loadtxt(fin).transpose()
    om=data[0]
    Sm = data[2]*1j
    beta = pi/om[0]
    print('beta=',beta)

    Sx = Pade(om, Sm, x, gamma, Norder)
    Ax = -imag(Gz(x-Sx))/pi

    np.savetxt('sig.pade', np.vstack((x,real(Sx),imag(Sx))).transpose())
    np.savetxt(fout, np.vstack((x,Ax)).transpose())

    
if __name__ == '__main__':
    import sys

    dat=np.loadtxt('Gf.out').transpose()
    om=dat[0]
    exec(compile(open('maxent_params.dat', "rb").read(), 'maxent_params.dat', 'exec'))
    beta = pi/om[0]
    tau = np.linspace(0,beta,params['Ntau']+1)

    for i in range(int(len(dat)/2)):
        Gm=dat[2*i+1]+dat[2*i+2]*1j
        Gt = InverseFourier(Gm, om, tau, beta,  params['Nf'], params['statistics'])
        np.savetxt('Gtau.dat.'+str(i), np.vstack((tau,Gt)).transpose() )
        (Aw, omega) = MaximumEntropy(params, tau, Gt)
        np.savetxt('Aw.dat.'+str(i), np.vstack((omega,Aw)).transpose() )
    sys.exit(0)


    Gw = np.loadtxt('Gimp').transpose()
    om = Gw[0]

    for i in range(2,6):
        Gm = Gw[2*i+1]+Gw[2*i+2]*1j
        (Aw,omega) = test1(om, Gm, 'fermi',0.001,400,15)
        np.savetxt('dos.out.'+str(i), np.vstack((omega,Aw)).transpose())
    
    
