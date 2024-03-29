#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
from numpy import *
from scipy import interpolate
from scipy import integrate
from scipy import special
#from pylab import *
import optparse

def SlaterF2J(Fk,l):
    if l==0:
        return Fk
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
    usage = """usage: %prog [ options ]
    Given an input screening length lambda and real space projector (projectorw.dat) it calculates the screened 
    Coulomb interaction for the orbital. Note that it returns U,J1,J2 and this should be further screened
    by dielectric constant, which is here set to 1. This screening is only Yukawa type screening.
    """
    parser = optparse.OptionParser(usage)
    parser.add_option("-n", "--npts",  dest="Nk",    type="int", default=8, help="Number of points in the radial integration will be 2**NK+1.")
    parser.add_option("-i", "--inw", dest="finput",  default='projectorw.dat', help="filename of the input file containing projector", metavar="FILE")
    parser.add_option("-l", "--lambda", dest="lmbda",  type="float", default=0.0,  help="The screening parameter lambda (in units of inverse bohr radius)")
    
    #parser.add_option("-o", "--sout", dest="outsig", default='sig.inp', help="filename-part of the output file (sig.inp)")
    #parser.add_option("-d", "--outDC", dest="outDC", default='Edc.dat', help="filename of the output DC-file (Edc.dat)")
    #parser.add_option("-l", "--lext", dest="m_extn", default='', help="For magnetic calculation, it can be 'dn'.")
    # Next, parse the arguments
    (options, args) = parser.parse_args()
    
    fp = open(options.finput, 'r')
    ncase, Nrmax = list(map(int, next(fp)[1:].split()[:2]))
    
    #k=8
    Nk=2**options.Nk+1
    
    CIN = 1/137.0359895**2
    Ry2eV = 13.60569193
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
        fA = interpolate.UnivariateSpline(Rx,Ag,s=0)
        fB = interpolate.UnivariateSpline(Rx,Bg,s=0)
        r = linspace(0,Rx[-1],Nk)
        An = fA(r)
        Bn = fB(r)
    
        ul2 = An**2 + CIN*Bn**2
    
        print('Integral of the wave function for l=',l,'is', integrate.romb(ul2,dx=r[1]-r[0]))

        if options.lmbda==0.0:
            Fk=[]
            for k in range(0,2*l+2,2):
                U_inside=zeros(len(r))
                for ir in range(1,len(r)):
                    U_inside[ir] = integrate.simps(ul2[:ir+1]*r[:ir+1]**k, x=r[:ir+1])
                
                U_outside=zeros(len(r))
                U_outside[1:] = 2*U_inside[1:]*ul2[1:]/r[1:]**(k+1)
                
                Fk.append( 2*Ry2eV*integrate.romb(U_outside, dx=r[1]-r[0]) )
            print('Fk=', Fk)
            print('U,J=', SlaterF2J(Fk,l))
        else:
            lmbda = options.lmbda
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
            print('Fk=', Fk)
            print('U,J=', SlaterF2J(Fk,l))
            
