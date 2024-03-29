# @Copyright 2007 Kristjan Haule
from numpy import *
from scipy import interpolate
from scipy import integrate
from scipy import special

def CmpSlaterInt(lmbda,Rx,Ag,Bg,l,Nk,prnt=False):
    Ry2eV = 13.60569193
    CIN = 1/137.0359895**2

    fA = interpolate.UnivariateSpline(Rx,Ag,s=0)
    fB = interpolate.UnivariateSpline(Rx,Bg,s=0)
    r = linspace(0,Rx[-1],Nk)
    An = fA(r)
    Bn = fB(r)
    
    ul2 = An**2 + CIN*Bn**2
    
    #if prnt: print >> fh_info, 'Integral of the wave function for l=',l,'is', integrate.romb(ul2,dx=r[1]-r[0])

    if lmbda==0.0:
        Fk=[]
        for k in range(0,2*l+2,2):
            U_inside=zeros(len(r))
            for ir in range(1,len(r)):
                U_inside[ir] = integrate.simps(ul2[:ir+1]*r[:ir+1]**k, x=r[:ir+1])
            
            U_outside=zeros(len(r))
            U_outside[1:] = 2*U_inside[1:]*ul2[1:]/r[1:]**(k+1)
            
            Fk.append( 2*Ry2eV*integrate.romb(U_outside, dx=r[1]-r[0]) )
    else:
        Fk=[]
        for k in range(0,2*l+2,2):
            U_inside=zeros(len(r))
            for ir in range(1,len(r)):
                y1 = special.iv(0.5+k,lmbda*r[1:ir+1])/sqrt(r[1:ir+1])
                if k==0:
                    y0 = lmbda*sqrt(2/pi)
                else:
                    y0 = 0
                y2 = hstack( ([y0],y1) )
                U_inside[ir] = integrate.simps(ul2[:ir+1]*y2, x=r[:ir+1])
            
            U_outside=zeros(len(r))
            y4 = special.kv(0.5+k,lmbda*r[1:])/sqrt(r[1:])
            U_outside[1:] = 2*(2*k+1)*U_inside[1:]*ul2[1:]*y4
            
            Fk.append( 2*Ry2eV*integrate.romb(U_outside, dx=r[1]-r[0]) )
    return Fk

    #UJs=SlaterF2J(Fk,l)
    #if prnt:
    #    print 'Fk=', Fk
    #    print 'U,J=', UJs
    #return UJs

def GetWaveF(finput):
    fp = open(finput, 'r')
    ncase, Nrmax = list(map(int, next(fp)[1:].split()[:2]))
    
    Rx_=[]
    Ag_=[]
    Bg_=[]
    ls_=[]
    for icase in range(ncase):
        (Nr,Nr0,jatom,l) = list(map(int,next(fp)[1:].split()[:4]))
        Rx=[]
        Ag=[]
        Bg=[]
        for ir in range(Nr):
            rx, ag, bg = list(map(float,next(fp).split()))
            Rx.append(rx)
            Ag.append(ag)
            Bg.append(bg)
        Rx_.append(Rx)
        Ag_.append(Ag)
        Bg_.append(Bg)
        ls_.append(l)
    return (Rx_, Ag_, Bg_, ls_)

def SlaterF2J(Fk,l):
    if l==0:
        return [Fk[0],0]
    elif l==1:
        # F2 for p-electrons
        J2 = Fk[1]/5.
        return [Fk[0],J2]
    elif l==2:
        # F2 and F4 for d-electrons
        J2 = Fk[1]*1.625/14.
        J4 = Fk[2]*1.625/(14.*0.625)
        return [Fk[0], J2, J4]
    elif l==3:
        J2 = Fk[1]*539.76/6435
        J4 = Fk[2]*539.76/(0.668*6435.)
        J6 = Fk[3]*539.76/(0.494*6435.)
        return [Fk[0], J2, J4, J6]

if __name__ == '__main__':
    Nk=2**8+1
    (Rx_, Ag_, Bg_, ls_) = GetWaveF('projectorw.dat')
    for icase in range(len(Rx_)):
        Rx = Rx_[icase]
        Ag = Ag_[icase]
        Bg = Bg_[icase]
        l = ls_[icase]

        lmbda = 0.0
        epsilon = 1.0
            
        print('lambda=', lmbda, 'epsilon=', epsilon)
        Fk = CmpSlaterInt(lmbda,Rx,Ag,Bg,l,Nk,True)
        UJs=SlaterF2J(Fk,l)
        print('Fk=', (array(Fk)/epsilon).tolist())
        print('U,J=', (array(UJs)/epsilon).tolist())
        
 
