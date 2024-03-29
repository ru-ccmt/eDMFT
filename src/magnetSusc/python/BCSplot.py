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

def bilinear_interpolation(x, y, points):
    '''Interpolate (x,y) from values associated with four points.

    The four points are a list of four triplets:  (x, y, value).
    The four points can be in any order.  They should form a rectangle.

        >>> bilinear_interpolation(12, 5.5,
        ...                        [(10, 4, 100),
        ...                         (20, 4, 200),
        ...                         (10, 6, 150),
        ...                         (20, 6, 300)])
        165.0

    '''
    # See formula at:  http://en.wikipedia.org/wiki/Bilinear_interpolation

    points = sorted(points)               # order points by x, then by y
    (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points

    if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
        raise ValueError('points do not form a rectangle')
    if not x1 <= x <= x2 or not y1 <= y <= y2:
        print 'x,y=', x,y, 'x1,x2=', x1, x2, 'y1,y2=', y1,y2
        raise ValueError('(x, y) not within the rectangle')

    return (q11 * (x2 - x) * (y2 - y) +
            q21 * (x - x1) * (y2 - y) +
            q12 * (x2 - x) * (y - y1) +
            q22 * (x - x1) * (y - y1)
            ) / ((x2 - x1) * (y2 - y1) + 0.0)

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


def ReshapeExtend(data1D,nk1):
    c_ = data1D.reshape(nk1,nk1)
    res = zeros((nk1+1,nk1+1),dtype=data1D.dtype)
    res[:nk1,:nk1] = c_[:,:]
    res[nk1,:nk1]  = c_[0,:]
    res[:nk1,nk1]  = c_[:,0]
    res[nk1,nk1]   = c_[0,0]
    return res

def ReadEigenvectors(ifso, first, fUL, fUR):
    fL = open(fUL, 'r')
    fR = open(fUR, 'r')
    dat = fL.next().split()
    fR.next()
    
    nkp,nsymop,nom,norbitals = map(int,dat[1:5])  # dimensions
    dims = map(int,dat[5:5+norbitals])
    fL.next()
    fR.next()
    
    maxbnd = max(ifso)-first
    
    coh=zeros((sum(dims),maxbnd+1,nkp),dtype=float)
    UL=zeros( (sum(dims),maxbnd+1,nkp),dtype=complex)
    UR=zeros( (sum(dims),maxbnd+1,nkp),dtype=complex)
    
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
                        tUL = gl[::2]+gl[1::2]*1j
                        tUR = gr[::2]+gr[1::2]*1j
                        coh[:,ii-first,ik] = abs(tUL*tUR)
                        UL[:,ii-first,ik] = tUL
                        UR[:,ii-first,ik] = tUR
    return (coh, UL, UR)

def Print(A):
    for i in range(shape(A)[0]):
        for j in range(shape(A)[1]):
            print "%7.3f " % A[i,j].real,
        print

def FindGap(En):
    for i in range(len(En)):
        if En[i]>0: break
    #print En[i-1], En[i]
    return real(En[i]-En[i-1])


def Interpolate((x,y),(ix,iy),Ex):
    "Given point (x,y) and bracketed index of the array Ex[ix+[0,1],iy+[0,1]], interpolates Ex"
    points=[(ix,iy,Ex[ix,iy]),(ix+1,iy,Ex[ix+1,iy]),(ix,iy+1,Ex[ix,iy+1]),(ix+1,iy+1,Ex[ix+1,iy+1])]
    return bilinear_interpolation(x, y, points)


def GetBogGap((xt,yt), Ebc, Dlt, ifBog, Delta_magnitude):
    ix,iy = floor(x),floor(y)
    NB = len(ifBog)
    epsk=zeros(NB)
    delt=zeros((NB,NB),dtype=complex)
    for ib in range(NB):
        epsk[ib] = Interpolate((xt,yt),(ix,iy),Ebc[ib,:,:])
        for jb in range(len(ifBog)):
            delt[ib,jb]=Interpolate((xt,yt),(ix,iy),Dlt[ib,jb,:,:])

    HBog=zeros((2*NB,2*NB),dtype=complex)
    for ib in range(NB):
        HBog[ib,ib] = epsk[ib]
        HBog[NB+ib,NB+ib] = -epsk[ib]
        for jb in range(NB):
            HBog[ib,NB+jb]=delt[ib,jb]*Delta_magnitude
            HBog[ib+NB,jb]=delt[jb,ib].conj()*Delta_magnitude
    lambd = linalg.eigvalsh(HBog)
    gap = FindGap(lambd)/(2*Delta_magnitude)
    return gap


if __name__ == '__main__':

    Delta_magnitude=1.
    PlotAbs=False
    
    if len(sys.argv)<2:
        print 'Give input BCS eigenvector.' 

    fBCS = sys.argv[1]
    fname = 'eigenvalues.dat'
    fUL = 'UL.dat_'
    fUR = 'UR.dat_'
    
    gs0 = loadtxt(fBCS).transpose()
    (norb_norb,nk2) = shape(gs0)
    nk1=int(round(sqrt(nk2)))
    norb=int(round(sqrt(norb_norb)))
    
    gs_=gs0.reshape((norb_norb,nk1,nk1))
    gs2 = zeros((norb_norb,nk1+1,nk1+1))
    gs2[:,:nk1,:nk1] = gs_[:,:,:]
    gs2[:,nk1,:nk1]  = gs_[:,0,:]
    gs2[:,:nk1,nk1]  = gs_[:,:,0]
    gs2[:,nk1,nk1]   = gs_[:,0,0]
    gs=[]
    for i in range(norb_norb):
        kxy = arange(nk1+1.)/float(nk1)
        gs.append( interpolate.RectBivariateSpline(kxy, kxy, gs2[i,:,:], kx=2,ky=2,s=0) )
    
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
    nk1 = int((sqrt(nkp)))
    print nkp, len(Ek), len(Nmin)
    

    first=min(Nmin)
    last=max(Nmax)+1
    print 'first=', first
    print 'last=', last
    
    ibnd_first = max([Nmin[ik]-first for ik in range(nkp)])
    ibnd_last = min([Nmin[ik]-first+len(Ek[ik]) for ik in range(nkp)])
    ifBog = range(ibnd_first,ibnd_last)
    print 'ibnd_first=', ibnd_first, 'ibnd_last=', ibnd_last, 'ifBog=', ifBog
    
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
        


    (coh,UL,UR) = ReadEigenvectors(ifso,first,fUL,fUR)

    gsk_ = zeros((len(ifBog),len(ifBog),nkp),dtype=complex)
    
    kx,ky = mgrid[0:nk1,0:nk1]
    kx = kx.ravel()/float(nk1)
    ky = ky.ravel()/float(nk1)
    gsk=[]
    cohk=[]
    for ik in range(nkp):
        gst = array([gs[i](kx[ik],ky[ik])[0,0] for i in range(norb_norb)]).reshape(norb,norb)
        gs_ij = dot( dot(UL[:,:,ik].T, gst), UR[:,:,ik] )
        gs_ii = [ gs_ij[ibn,ibn] for ibn in ifs ]
        coh_ii = [ sum([UL[iorb,ibn,ik]*UR[iorb,ibn,ik] for iorb in range(norb)]) for ibn in ifs ]                  # coherence factors for the Fermi surface
        gsk_[:,:,ik]=gs_ij[ifBog[0]:ifBog[-1]+1,ifBog[0]:ifBog[-1]+1]
        gsk.append( gs_ii )
        cohk.append( coh_ii )
        
    Ebc=zeros( (len(ifBog),nk1+1,nk1+1), dtype=float)
    Dlt=zeros( (len(ifBog),len(ifBog),nk1+1,nk1+1),dtype=complex)
    for i,ib in enumerate(ifBog):
        Ebc[i,:,:] = real(ReshapeExtend(Ebnd[ib,:],nk1))
        for j,jb in enumerate(ifBog):
            Dlt[i,j,:,:]=ReshapeExtend(gsk_[i,j,:],nk1)

    
    gsk=array(gsk)
    cohk=array(cohk)
    
    print 'nk1=', nk1

    print 'shape(gsk)=', shape(gsk)
    fDelta=[]
    fcoh=[]
    for i,ib in enumerate(ifs):
        gs1 = ReshapeExtend(gsk[:,i],nk1)
        fDelta.append( interpolate.RectBivariateSpline(range(nk1+1), range(nk1+1), real(gs1), s=0 )  )
        coh1 = ReshapeExtend(cohk[:,i],nk1)
        fcoh.append( interpolate.RectBivariateSpline(range(nk1+1), range(nk1+1), abs(coh1), s=0 )  )
        #imshow(real(gs1),origin='lower',  interpolation='bilinear', extent=[0,1,0,1], aspect=1. )
        #colorbar()
        #ebn = ReshapeExtend(Ebnd[ib,:],nk1)
        #contour(ebn, (0,), colors='k', origin='lower', linewidths = 1, hold='on', extent=(0,1,0,1))
        #show()
        #print 'shape(gs1)=', shape(gs1), gs1[0,0], gs1[nk1/4,nk1/4], gs1[nk1/2,nk1/2]
        
    print 'shape(Ebnd)=', shape(Ebnd)
    
    segw=[]
    for j,i in enumerate(ifs):
        #Ebc = real(Ebnd[i,:]).reshape((nk1,nk1))
        Ebc_ = ReshapeExtend(Ebnd[i,:],nk1)
        cs = contour(Ebc_, (0,), colors='k', origin='lower', linewidths = 1, hold='on', extent=(0,nk1,0,nk1))
        p = cs.collections[0].get_paths()
        segu=[]
        for path in p:
            segs=[]
            x0,y0=None,None
            for (vertex,code) in path.iter_segments():
                xt,yt = vertex[0],vertex[1]
                if x0!=None:
                    segs.append( [(x0,y0),(xt,yt)])
                x0,y0=xt,yt
                #print xt, yt
            segu.append( segs )
        segw.append( segu )
    plt.axes().set_aspect('equal')
    show()


    zdl,zgp,zch=[],[],[]
    for j,i in enumerate(ifs):
        udl,uch,ugp=[],[],[]
        for p in range(len(segw[j])):
            wdl=[]
            wch=[]
            wgp=[]
            wx=[]
            wy=[]
            for si,((x0,y0),(xt,yt)) in enumerate(segw[j][p]):
                x=0.5*(x0+xt)
                y=0.5*(y0+yt)
                ch=fcoh[j](x,y)[0,0]   # How much of this band after unfolding (total character)
                dl=fDelta[j](x,y)[0,0] # projection of gap Delta to bands / only diagonal component
                gap = GetBogGap((x,y), Ebc, Dlt, ifBog, Delta_magnitude)
                wdl.append(dl)
                wch.append(ch)
                wgp.append(gap)
                wx.append(x)
                wy.append(y)
            if sum(wch)/len(wch)>0.3:
                fsp=open('fspath.'+str(j)+'.'+str(p),'w')
                print >> fsp, '#  x   y   |gap|   gap   character'
                for si in range(len(wx)):
                    print >> fsp, wx[si], wy[si], wgp[si], wdl[si], wch[si]
            udl.append(wdl)
            ugp.append(wgp)
            uch.append(wch)
        zdl.append(udl)
        zgp.append(ugp)
        zch.append(uch)



    cmaps=['RdBu','Purples','Blues','Greens','Reds']
    if PlotAbs:
        Zmax = max([max([max(zgp[j][p]) for p in range(len(zgp[j]))]) for j in range(len(zgp))])
        Zmin=0
        ZM = [Zmin,Zmax]
        imap=4
    else:
        Zmax = max([max([max(zdl[j][p]) for p in range(len(zdl[j]))]) for j in range(len(zdl))])
        Zmin = min([min([min(zdl[j][p]) for p in range(len(zdl[j]))]) for j in range(len(zdl))])
        ZM = [-max(abs(Zmin),Zmax),max(abs(Zmin),Zmax)]
        imap=0
    
    print 'ZM=', ZM
    for j,i in enumerate(ifs):
        for p in range(len(segw[j])):
            print j,p, len(zdl[j][p]), len(segw[j][p])
            lc = LineCollection(segw[j][p], cmap=plt.get_cmap(cmaps[imap]),norm=plt.Normalize(ZM[0],ZM[1]))
            if PlotAbs:
                lc.set_array(array(zgp[j][p]))               # sets color
            else:
                lc.set_array(array(zdl[j][p]))               # sets color
            lc.set_linewidth( abs(array(zch[j][p]))*10 ) # sets linewidth
            plt.gca().add_collection(lc)
            
    plt.xlim(0, nk1)
    plt.ylim(0, nk1)
    
    fig = gcf()
    axcb = fig.colorbar(lc)
    axcb.set_label('SC gap')
    show()
