#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule

from scipy import *
from pylab import *
from scipy import interpolate
import os, sys
#from mayavi.mlab import *
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import scipy

from distutils.version import StrictVersion
if StrictVersion(scipy.__version__) > StrictVersion('0.19.0'):
    import weave
else:
    import scipy.weave as weave

code0="""
     using namespace std;
     double Ax = 0;
     for (int ib=0; ib<ekw.extent(0); ib++){
        complex<double> gammac(0,small);
        complex<double> gc = 1./(omega+mu-ekw(ib)+gammac);
        Ax += -gc.imag()/M_PI;
     }
     return_val = Ax;
     """

def ReadEigenvectors(ifso, first):
    fL = open('UL.dat_', 'r')
    fR = open('UR.dat_', 'r')
    dat = fL.next().split()
    fR.next()
    
    nkp,nsymop,nom,norbitals = map(int,dat[1:5])  # dimensions
    dims = map(int,dat[5:5+norbitals])
    fL.next()
    fR.next()
    
    maxbnd = max(ifso)-first
    
    coh=zeros((sum(dims),maxbnd+1,nkp),dtype=float)
    for ik in range(nkp):
        for isym in range(nsymop):
            for iom in range(nom):
                line1 = fL.next()
                fR.next()
                dat = line1.split()
                omega = float(dat[0])
                nbands = int(dat[1])
                for ibnd in range(nbands):
                    datL = fL.next().split()
                    datR = fR.next().split()
                    ii = int(datL[0])
                    if ii in ifso and iom==0:
                        gl = array(map(float,datL[1:]))
                        gr = array(map(float,datR[1:]))
                        UL = gl[::2]+gl[1::2]*1j
                        UR = gr[::2]+gr[1::2]*1j
                        coh[:,ii-first,ik] = abs(UL*UR)
    print 'shape(coh)=(norbs,nband,nkp)=', shape(coh)
    return coh


def PrintDataExplorer():
    ff = open('FermiData.dat','w')
    print >> ff, '-100',nk1,nk1,1,len(ifs), ' : sorting number, grid size, number of bands'
    for ik1 in range(nk1):
        for ik2 in range(nk1):
            print >> ff, -0.5+2*ik1/float(nk1) , -0.5+2*ik2/float(nk1), 
            for j,i in enumerate(ifs):
                print >> ff, Ebnd[i,ik1*nk1+ik2].real,
            print >> ff
    ff.close()
    ff = open('fermi.general','w')
    print>> ff, 'file = FermiData.dat'
    print>> ff, 'grid =  '+str(nk1)+' x '+str(nk1)
    print>> ff, 'format = ascii'
    print>> ff, 'interleaving = field'
    print>> ff, 'majority = column'
    print>> ff, 'header = lines 1'
    print>> ff, 'field = locations,',' '.join(["field"+str(i)+',' for i in range(len(ifs))])
    print>> ff, 'structure = 2-vector,', 'scalar, '*len(ifs)
    print>> ff, 'type = float,', 'float, '*len(ifs)
    print>> ff
    print>> ff, 'end'
    ff.close()
    
    for iorb in range(len(coh)):
        ff = open('Coherence_'+str(iorb)+'.dat','w')
        print >> ff, '-100',nk1,nk1,1,len(ifs), ' : sorting number, grid size, number of bands'
        for ik1 in range(nk1):
            for ik2 in range(nk1):
                print >> ff, -0.5+2*ik1/float(nk1) , -0.5+2*ik2/float(nk1), 
                for j,i in enumerate(ifs):
                    print >> ff, coh[iorb,i,ik1*nk1+ik2],
                print >>ff
        ff.close()
        ff = open('coherence.general_'+str(iorb),'w')
        print>> ff, 'file = '+'Coherence_'+str(iorb)+'.dat'
        print>> ff, 'grid =  '+str(nk1)+' x '+str(nk1)
        print>> ff, 'format = ascii'
        print>> ff, 'interleaving = field'
        print>> ff, 'majority = column'
        print>> ff, 'header = lines 1'
        print>> ff, 'field = locations,',' '.join(["field"+str(i)+',' for i in range(len(ifs))])
        print>> ff, 'structure = 2-vector,', 'scalar, '*len(ifs)
        print>> ff, 'type = float,', 'float, '*len(ifs)
        print>> ff
        print>> ff, 'end'
        ff.close()
        

if __name__ == '__main__':

    COHERENCE=True
    onlyHalfPoints=True
    orb_plot=range(5)
    orb_plot=[0,0,1,1,0,0,2,2,2,2,2,2,2,2,0,0,1,1,0,0,2,2,2,2,2,2,2,2,]
    norb_plot=max(orb_plot)+1
    
    #fname = 'eigvals.dat'
    #fname = 'eigenvalues.dat'
    if len(sys.argv)<2:
        print 'Give filename for eigenvalues'
        sys.exit(0)

    fname = sys.argv[1]

    small = 0.005
    intensity=0.99
    
    fEF = open('EF.dat', 'r')
    mu = float(fEF.next().split()[0])
    print 'mu=', mu
    
    fdat = open(fname, 'r')
    
    formt=None

    ikp=0
    Ek=[]
    Nmin=[]
    Nmax=[]
    try:
        while True:
            #data = fdat.next().split()
            strn = fdat.next()
            data = strn[1:].split()
            if formt==None:
                comm = strn[0]
                if comm=='#':                                                                                                                                                                
                    formt='dmft1'                                                                                                                                                            
                elif comm=='!':                                                                                                                                                              
                    formt='dmftgk'                                                                                                                                                           
                else:                                                                                                                                                                        
                    print 'Can not recognize format of', fname
                    sys.exit(0)
            if formt=='dmft1':
                (ikp, isym, nbands, nemin, nomega) = map(int, data[:5])
            else:
                (ikp, isym, tot_nbands, tot_nemin, nomega) = map(int, data[:5])
            wgh = float(data[5])
            
            Ekw=[]
            Nemin=[]
            Nemax=[]
            for iom in range(nomega):
                dat = fdat.next().split()
                omega = float(dat[0])
                if formt=='dmft1':
                    data = array(map(float, dat[1:]))
                else:
                    nbands = int(dat[1])
                    nemin = int(dat[2])
                    data = array(map(float, dat[3:]))
                ekw = data[::2]+data[1::2]*1j - mu
            
                if len(ekw)!=nbands: print 'ERROR nbands=', nbands, 'len(ekw)=', len(ekw), 'ekw=', ekw.tolist()

                Ekw.append( ekw )
                Nemin.append( nemin )
                Nemax.append( nemin+nbands )
            Ek.append( Ekw[0] )
            Nmin.append(Nemin[0])
            Nmax.append(Nemax[0])
            ikp += 1
    except StopIteration:
        pass

    nkp=ikp-1
    
    first=min(Nmin)
    last=max(Nmax)+1
    print 'first=', first
    print 'last=', last
    
    if (onlyHalfPoints):
        nk1 = int((-1+sqrt(1+8*nkp))/2)
        # Store energy at EF for all bands
        Ebnd = zeros((last-first,nk1*nk1),dtype=complex)    
        ik=0
        for ik1 in range(nk1):
            for ik2 in range(ik1,nk1):
                ika = ik1*nk1+ik2
                ikb = ik2*nk1+ik1
                for i in range(len(Ek[ik])):
                    Ebnd[Nmin[ik]+i-first,ika] = Ek[ik][i]
                    Ebnd[Nmin[ik]+i-first,ikb] = Ek[ik][i]
                ik+=1
    else:
        nk1 = int((sqrt(nkp)))
        print nkp, len(Ek), len(Nmin)
        
        # Store energy at EF for all bands
        Ebnd = zeros((last-first,nkp),dtype=complex)    
        for ik in range(nkp):
            for i in range(len(Ek[ik])):
                #if Nmin[ik]+i-first<0 or Nmin[ik]+i-first>len(Ebnd): print 'ERROR'
                Ebnd[Nmin[ik]+i-first,ik] = Ek[ik][i]
            
    # Which bands cross EF?
    ifs=[]
    ifso=[]
    for i in range(first,last):
        Eband = Ebnd[i-first,:].real
        if min(Eband)<0.0 and max(Eband)>0.0:
            ifs.append(i-first)
            ifso.append(i)
            
    print 'ifs=', ifs
    print 'ifso=', ifso
    
    # Some bands are not traced through entire BZ. They might hence 
    # remind zero in some region, because of initial condition for Ebnd=zeros
    # Will set Ebnd to average energy for this band to avoid artificial EF crossings
    for i in ifs:
        Ecenter = sum(Ebnd[i,:].real)/nkp
        for j in range(nkp):
            if Ebnd[i,j]==0.0: Ebnd[i,j]=sign(Ecenter)
            
    if COHERENCE:
        coh = ReadEigenvectors(ifso,first)
        norb=len(coh)
        fcohr=[]
        for i in ifs:
            
            cohi = zeros((norb_plot,nk1,nk1))
            for iorb in range(norb):
                if onlyHalfPoints:
                    cohi_ = zeros((nk1,nk1))
                    ik=0
                    for ik1 in range(nk1):
                        for ik2 in range(ik1,nk1):
                            cohi_[ik1,ik2] = coh[iorb,i,ik]
                            cohi_[ik2,ik1] = coh[iorb,i,ik]
                            ik+=1
                else:
                    cohi_ = coh[iorb,i,:].reshape(nk1,nk1)
                cohi[ orb_plot[iorb], :,:] += cohi_[:,:]
                
            fcoh=[]
            for iorb in range(norb_plot):
                fcoh.append( interpolate.RectBivariateSpline(range(nk1), range(nk1), cohi[iorb,:,:],s=0 ) )
            fcohr.append( fcoh )
            
        print 'shape(fcohr)=', norb, len(ifs)

    print 'shape(Ebnd)=', shape(Ebnd)


    segw,zw=[],[]
    for j,i in enumerate(ifs):
        Ebc = zeros((nk1+1,nk1+1))
        Ebc0 = real(Ebnd[i,:]).reshape((nk1,nk1))
        Ebc[:nk1,:nk1] = Ebc0
        Ebc[nk1,:nk1] = Ebc0[0,:]
        Ebc[:nk1,nk1] = Ebc0[:,0]
        Ebc[nk1,nk1] = Ebc0[0,0]
        
        cs = contour(Ebc, (0,), colors='k', origin='lower', linewidths = 1, hold='on', extent=(0,nk1,0,nk1))
        p = cs.collections[0].get_paths()
        segu,zu=[],[]
        for path in p:
            print 'new path:'
            segs,zs=[],[]
            x0,y0=None,None
            for (vertex,code) in path.iter_segments():
                xt,yt = vertex[0],vertex[1]
                if x0!=None:
                    segs.append( [(x0,y0),(xt,yt)])
                x0,y0=xt,yt
                print xt, yt
                if COHERENCE: zs.append( [ fcohr[j][iorb](xt,yt)[0,0]  for iorb in range(norb_plot)] )

            segu.append( segs )
            zu.append( zs )
        segw.append( segu )
        zw.append( zu )
    plt.axes().set_aspect('equal')
    show()

    if COHERENCE:
        zwn=[[[[zw[j][p][t][iorb]  for t in range(len(zw[j][p]))] for p in range(len(zw[j]))] for j in range(len(zw))] for iorb in range(norb_plot)]
        
        cmaps=['Reds', 'Greens', 'Purples', 'Blues', 'Purples', 'BuGn']
        for iiorb in range(norb_plot):
            print iiorb, 'Using cmap=', cmaps[iiorb]
            for j,i in enumerate(ifs):
                for p in range(len(segw[j])):
                    print j,p, len(zwn[iiorb][j][p]), len(segw[j][p])
                    lc = LineCollection(segw[j][p], cmap=plt.get_cmap(cmaps[iiorb]),norm=plt.Normalize(0, 1))
                    lc.set_array(array(zwn[iiorb][j][p]))
                    lc.set_linewidth( abs(array(zwn[iiorb][j][p]))*10 )
                    plt.gca().add_collection(lc)
            
    plt.xlim(0, nk1)
    plt.ylim(0, nk1)
    plt.axes().set_aspect('equal')
    show()
    
    
