#!/usr/bin/env python
from scipy import *
from pylab import *
import sys
from scipy import integrate

def ReadCix(filename):
    fi = open(filename)
    [fi.next() for i in range(2)]
    (Nc,Ns,Nb,mmax) = map(int,fi.next().split()[:4])
    fi.next()
    [fi.next() for ib in range(Nb)]
    [fi.next() for i in range(3)]
    Nv=[]
    Sz=[]
    vsize=[]
    for i in range(Ns):
        dat = fi.next().split()
        Nv.append(int(dat[1]))
        Sz.append(float(dat[3]))
        vsize.append( int(dat[4]))
    return (Nv, vsize, Sz)

def Print(A):
    for i in range(shape(A)[0]):
        for j in range(shape(A)[1]):
            if (A.dtype==complex): print "%6.3f+%6.3f*i " % (real(A[i,j]),imag(A[i,j])),
            else: print "%6.3f " % (real(A[i,j]),),
        print

if __name__ == '__main__':

    if len(sys.argv)<=1:
        dirc = './'
    else:
        dirc = sys.argv[1]+'/'
    
    U=7.
    broad = U/20.
    
    dat = loadtxt(dirc+'Probability.dat')
    datT = loadtxt(dirc+'TProbability.dat')
    (Nv, vsize, Sz) = ReadCix(dirc+'actqmc.cix')

    Pc = zeros(len(Nv))
    ii=0
    for i in range(len(Nv)):
        p=0
        for j in range(vsize[i]):
            p += float(dat[ii][2])
            ii += 1
        Pc[i]=p

    xm = linspace(-max(20*broad,U/2), max(20*broad,U/2)+U*max(Nv), 2**10+1)
    
    
    Ns = set(Nv)
    Ps={}
    ii=0
    for i in range(len(Nv)):
        N = Nv[i]
        if Ps.has_key(N):
            Ps[N] += Pc[i]
        else:
            Ps[N] = Pc[i]
            
        
    G0 = zeros(len(xm), dtype=complex)
    for i in range(len(Nv)):
        G0 += Pc[i]/(xm-U*Nv[i]+broad*1j)
            
    if True:
        # Each atomic state is considered, and transition probabilities between them
        TP = zeros((len(Nv),len(Nv)))
        for ii in range(len(datT)):
            TP[datT[ii,0]-1,datT[ii,1]-1] = datT[ii,2]
        # due to numerics, the transition probability is not really symmeric, but it should be
        TP = 0.5*(transpose(TP)+TP)
        
        print 'Transition probability for each atomic state:'
        Print(TP)
        
        V = zeros((len(Nv),len(Nv)),dtype=complex)
        Gm = zeros(len(xm),dtype=complex)
        for ix,x in enumerate(xm):
            for i in range(len(Nv)):
                V[i,i] = x-U*Nv[i]+broad*1j
            G = linalg.inv(V-TP)
            dsum=0
            for i in range(len(Nv)):
                dsum += Pc[i]*G[i,i]
            Gm[ix] = dsum
        
            
        #plot(xm, -imag(Gm)/pi, label='G_solid')
        #plot(xm, -imag(G0)/pi, label='G_atom')
        #legend(loc='best')
        #show()
    
    else:
    # Only each valence is considered
    
        Nn = max(Nv)+1
        TT = zeros((Nn,Nn))
        for ii in range(len(datT)):
            i0 = int(datT[ii,0]-1)
            i1 = int(datT[ii,1]-1)
            N0 = Nv[i0]
            N1 = Nv[i1]
            tp = datT[ii,2]
            TT[N0,N1] += tp
            #print N0, N1, tp
            
        print 'Transition probability for each valence:'
        print TT
        
        V = zeros((Nn,Nn),dtype=complex)
        Gm = zeros(len(xm),dtype=complex)
        for ix,x in enumerate(xm):
            for n in range(Nn):
                V[n,n] = x-U*n+broad*1j
            G = linalg.inv(V-TT)
            dsum=0
            for n in range(Nn):
                dsum += Ps[n]*G[n,n]
            Gm[ix] = dsum

    Ac = -imag(Gm)/pi
    norm = integrate.simps(Ac, xm)
    print 'norm=', norm

    maxim=[]
    imax=[]
    dAc = gradient(Ac)
    for ix in range(len(xm)-1):
        if dAc[ix]*dAc[ix+1]<0:
            if Ac[ix]>Ac[ix+2] and Ac[ix]>Ac[ix-2]: # local maximum
                # f(x) = f0 + (f1-f0)*(x-x0)/(x1-x0) = 0
                #  x = x0 - f0*(x1-x0)/(f1-f0)
                x0 = xm[ix]  - dAc[ix]*(xm[ix+1]-xm[ix])/(dAc[ix+1]-dAc[ix])
                maxim.append(x0)
                if Ac[ix]>Ac[ix+1]:
                    imax.append(ix)
                else:
                    imax.append(ix+1)
    print 'maximums=', maxim
    #print 'imaxs=', imax
    weight=[]
    ia = 0
    for i in range(len(imax)):
        if i<len(imax)-1:
            ib = (imax[i]+imax[i+1])/2
        else:
            ib = len(xm)
        tweight = integrate.simps(Ac[ia:ib], xm[ia:ib])/norm
        weight.append(tweight)
        print "peak=%7.4f  weight=%7.5f  P=%7.5f" % (maxim[i], tweight, Ps[i])
        ia = ib
        
    
    fo = open(dirc+'core_hole.dat','w')
    for n in range(len(maxim)):
        print >> fo, maxim[n], weight[n], Ps[n]
    
    #plot(xm, dAc, 'o-')
    
    
    plot(xm, Ac, '-', label='Ac_solid')
    plot(xm, -imag(G0)/pi, '-.', label='Ac_atom')
    legend(loc='best')
    grid()
    show()
    

        
    sys.exit(0)
    y = zeros(len(xm))
    for n in Ps.keys():
        y += exp(-(xm-n)**2/0.05)*Ps[n]
    plot(w,y)
    show()
    
    sys.exit(0)
     

       
