#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule

from scipy import *
from scipy import linalg, integrate
import sys
import LinLogMesh
import krams
from pylab import *
import scipy

from distutils.version import StrictVersion
if StrictVersion(scipy.__version__) > StrictVersion('0.19.0'):
    import weave
else:
    import scipy.weave as weave

class FermGF:
    def __init__(self, Gf, nomv):
        self.g = Gf
        self.nomv = nomv  # we enumerate Gf for only first nomv Matsubara points. Gf might contain more points.
    def __getitem__(self,m):
        """Gives G(iom) for positive and negative frequencies. We have array for positive Matsubara points. 
           Negative Matsubara points are obtained from analytic property G(-iom)=G(iom)^*
           We enumerate matsubara points {0:-2*nomv+1, 1:-2*nomv+3,.... nomv-1:1, nomv:1,... , 2*nomv-1: 2*nomv-1}
        """
        if m>=self.nomv:
            return self.g[m-self.nomv]
        else:
            return self.g[self.nomv-1-m].conjugate()
    def size(self):
        return 2*self.nomv
    
def ReadKlist(fklist, ReadBS=False):
    fk = open(fklist,'r')
    data = fk.readlines()
    nkp = [line[:3]=='END' for line in data].index(True)
    if data[nkp][:3]!='END': 
        print 'wrong klist ', fklist
    kp=[]
    for i in range(nkp):
        kp.append( map(int, [data[i][10:15], data[i][15:20], data[i][20:25], data[i][25:30]]) )

    if (ReadBS):
        BS = [map(float,line.split()) for line in data[nkp+1:nkp+4]]
        BSI = matrix(array(BS).T).I
        return (array(kp), array(BS), array(BSI))
    else:
        return array(kp)


def CheckKpoints():
    BSI = matrix(BS.T).I
    for ik in range(len(kps)):
        k = array(kps[ik][:3])/float(kps[ik][3])
        print ik, k, array(dot(BSI,k))[0]

def PrintM(A):
    if A.dtype=='float':
        for i in range(len(A)):
            for j in range(len(A)):
                print "%8.3f  " % A[i,j],
            print
    else:
        for i in range(len(A)):
            for j in range(len(A)):
                print "%8.3f %8.3f   " % (A[i,j].real,A[i,j].imag),
            print

        
class K_index:
    def __init__(self, BSI, kps):
        self.BSI = BSI
        self.SCALE = kps[0][3]
        self.ind1={}
        for ik,k in enumerate(kps):
            wik = tuple(map(int, dot(BSI,k[:3])))
            self.ind1[wik] = ik
    def __call__(self, ik):
        wik = tuple(map(int, dot(self.BSI,ik[:3])%self.SCALE))
        return self.ind1[wik]
        
    
def Get_k_m_q_index_Python(kps, k_index):
    k_m_q = zeros((len(kps),len(kps)),dtype=int)
    for ik,k in enumerate(kps):
        print 'ik', ik, 'done'
        for iq,q in enumerate(kps):
            k_m_q[ik,iq] = k_index(k[:3]-q[:3])
    return k_m_q


def Get_k_m_q_index_CPP(kps):
    SCALE = int(kps[0][3])
    ind1 = zeros((SCALE,SCALE,SCALE),dtype=int)
    k_m_q = zeros((len(kps),len(kps)),dtype=int)
    
    support_code="""
    #line 93 "Suscept.py"
    using namespace blitz;
    inline void k_index(TinyVector<int,3>& kind, const Array<long int,1>& kq, const Array<double,2>& BSI, int SCALE)
    {
        for (int i=0; i<3; i++)
           kind(i) = int(BSI(i,0)*kq(0)+BSI(i,1)*kq(1)+BSI(i,2)*kq(2)+SCALE) % SCALE;
    }
    """
    code="""
    #line 102 "Suscept.py"
    using namespace blitz;

    TinyVector<int,3> kind(3);
    Array<long int,1> kq(3);
    for (int ik=0; ik<kps.extent(0); ik++){
        kq = kps(ik,Range::all());
        k_index(kind,kq,BSI,SCALE);
        ind1(kind) = ik;
    }
    for (int ik=0; ik<kps.extent(0); ik++){
        for (int iq=0; iq<kps.extent(0); iq++){
            kq = kps(ik,Range::all())-kps(iq,Range::all());
            k_index(kind,kq,BSI,SCALE);
            k_m_q(ik,iq) = ind1(kind);
        }
    }
    """
    weave.inline(code, ['k_m_q','kps','BSI','ind1', 'SCALE'],support_code=support_code,type_converters=weave.converters.blitz, compiler='gcc')
    return k_m_q



def ReadGk(filegk, nkp, nsymop, nom, cixdm):
    
    gkdat = loadtxt(filegk)
    gkdat2 = gkdat.reshape(nkp,nsymop,nom,2*cixdm*cixdm+1)
    om = zeros(nom,dtype=float)
    gk = zeros((cixdm,cixdm,nom,nkp),dtype=complex)
    
    for ik in range(nkp):
        for isym in range(nsymop):
            for iom in range(nom):
                gg = gkdat2[ik,isym,iom]
                om[iom] = gg[0]
                gs = (gg[1::2]+gg[2::2]*1j).reshape(cixdm,cixdm)
                gk[:,:,iom,ik] = gs
    gk *= 1./nsymop
                    
    return (gk,om)

def ReadGlc(fileglc, nom, cixdm):
    om = zeros(nom,dtype=float)
    Gf = zeros((cixdm,nom),dtype=complex)
    
    gdat = loadtxt(fileglc)
    for iom in range(nom):
        om[iom] = gdat[iom][0]
        Gf[:,iom] = gdat[iom][1::2]+gdat[iom][2::2]*1j

    return (om,Gf)

def ReadVertex(fvertex):
    fi=open(fvertex)
    fi.next()  # comment # beta, Nvfl, nomv, nOm nom
    data=fi.next().split()
    (beta, Nvfl, nomv, nOm, nom) = [float(data[0])] + map(int,data[1:])
    
    bfl_index=zeros(Nvfl,dtype=int) # bath index 
    for ib in range(Nvfl):
       t, bfl_index[ib] = map(int,fi.next().split())
       if t!=ib:
           print "Wrong format of ",fvertex
           sys.exit(1)
    
    fi.next()  # comment # b0 b1 Om om1
    VertexH=zeros((Nvfl,Nvfl,2*nOm-1,2*nomv,2*nomv),dtype=complex)
    VertexF=zeros((Nvfl,Nvfl,2*nOm-1,2*nomv,2*nomv),dtype=complex)
    for i0 in range(Nvfl):
       for i1 in range(Nvfl):
          for iOm in range(2*nOm-1):
             dOm=nOm-1-iOm
             sm0 = max(0, -dOm)
             em0 = min(2*nomv, 2*nomv-dOm)
             for im0 in range(sm0,em0):
                 data = fi.next().split()
                 ti0,ti1,Om,om1 = map(int,data[:2]) + map(float,data[2:4])
                 if ti0!=i0 or ti1!=i1:
                     print "Wrong format of ", fvertex
                     sys.exit(1)
                 sm1 = max(0,-dOm)
                 em1 = min(2*nomv,2*nomv-dOm)
                 for im1 in range(sm1,em1):
                     data = map(float,fi.next().split())
                     VertexH[i0,i1,iOm,im0,im1] = data[1]+data[2]*1j
                     VertexF[i0,i1,iOm,im0,im1] = data[3]+data[4]*1j
    
    return (beta,Nvfl,nomv,nOm,bfl_index,VertexH,VertexF)

def Cmp_chi_0(gf,nOm):
    norb = len(gf)
    nomv = gf[0].nomv
    chi0=zeros((2*nOm-1,norb*2*nomv),dtype=complex)
    for iOm in range(2*nOm-1):
        dOm=nOm-1-iOm
        for i in range(norb):
            for im in range(2*nomv):
                chi0[iOm,2*nomv*i+im]= -gf[i][im]*gf[i][im+dOm]
    return chi0


def Cmp_chi_loc(VertexH, VertexF, orb_ud, iOm, nomv):
    norb = len(orb_ud)
    
    ChiS = zeros((norb,norb,2*nomv,2*nomv),dtype=complex)
    ChiC = zeros((norb,norb,2*nomv,2*nomv),dtype=complex)
    for iorb1,(up1,dn1) in enumerate(orb_ud):
        VF_uu = 0.5*(VertexF[up1,up1,iOm]+VertexF[dn1,dn1,iOm]) 
        for iorb2,(up2,dn2) in enumerate(orb_ud):
            VH_uu = 0.5*(VertexH[up1,up2,iOm]+VertexH[dn1,dn2,iOm]) 
            VH_ud = 0.5*(VertexH[up1,dn2,iOm]+VertexH[dn1,up2,iOm]) 
            V_ud = VH_ud
            if iorb1!=iorb2:
                V_uu = VH_uu
            else:
                V_uu = VH_uu - VF_uu
            ChiS[iorb1,iorb2] = V_uu-V_ud
            ChiC[iorb1,iorb2] = V_uu+V_ud
            
    return (ChiS,ChiC)


def Cmp_chi0_Q(gk,nomv,dOm,norb,Qi,k_m_q):
    
    support_code="""
    using namespace std;
    complex<double> gkm(int iorb, int jorb, int im, int ik, blitz::Array<complex<double>,4>& gk, int nomv){
        if (im>=nomv) return gk(iorb,jorb,im-nomv,ik);
        else return conj(gk(jorb,iorb,nomv-1-im,ik));
    }
    """
    codeBubQ="""
    #line 269 "Suscept.py"
    using namespace std;
    for (int iorb=0; iorb<norb; iorb++){
       for (int jorb=0; jorb<norb; jorb++){
          for (int im=0; im<2*nomv; im++){
             complex<double> csum=0;
             for (int ik=0; ik<nkp; ik++){
                int ikq=k_m_q(ik,Qi);
                csum += -gkm(iorb,jorb,im,ik,  gk,nomv) * gkm(jorb,iorb,im+dOm,ikq, gk,nomv);
             }
             BubQ(iorb,jorb,im) = csum/static_cast<double>(nkp);
          }
       }
    }
    """
    nkp = len(k_m_q)
    if shape(gk)[3]!=nkp : print 'gk does not contain enough k-points', nkp, shape(gk)[3]
    BubQ=zeros((norb,norb,2*nomv),dtype=complex)
    weave.inline(codeBubQ, ['gk','nomv','dOm','norb','nkp','Qi','k_m_q','BubQ'],support_code=support_code,type_converters=weave.converters.blitz, compiler='gcc')
    return BubQ



    
def CmpRealAxisBubble(filemesh, filegkr, fbubbleReal, Qlist, k_index):

    def Cmp_ChiQ_real(iOm,Q,gkr,k_m_q,nkp,norb,dx_,idxl,Nd,zero_ind,k_index):
        level = (iOm-1)/Nd                      # which linear mesh should we use?
        om_idx=idxl[level]                      # index for the linear mesh on this level
        izero= (len(om_idx)-1)/2                # zero_ind on this level
        dOm = om_idx.index(zero_ind+iOm)-izero  # om-Om in integer notation is i-dOm
        om_idx=array(om_idx)                    
        
        Qi = k_index(Q)
        codeBub="""
            #line 239 "Suscept.py"
            using namespace std;
            for (int iorb=0; iorb<norb; iorb++){
               for (int jorb=0; jorb<norb; jorb++){
                  complex<double> csum=0.0;
                  for (int iom=izero; iom<izero+dOm+1; iom++){
                     double dom =  (iom==izero || iom==izero+dOm) ? dx/2.0 : dx;  // trapezoid rule has 1/2 at the last interval
                     complex<double> csum2=0.0;
                     for (int ik=0; ik<nkp; ik++){
                        int ikq  = k_m_q(ik,Qi);
                        int iom1 = om_idx(iom);
                        int iom2 = om_idx(iom-dOm);
                        complex<double> irhok=(gkr(iorb,jorb,iom1,ik)-conj(gkr(jorb,iorb,iom1,ik)))/2.0;             // (G_k-G_k^+)/2
                        complex<double> irhoq=(gkr(jorb,iorb,iom2,ikq)-conj(gkr(iorb,jorb,iom2,ikq)))/2.0;
                        csum2 -= irhok*irhoq; //  rho_k * rho_{k+q}
                     }
                     csum += csum2*dom;
                  }
                  ImBub(iorb,jorb) = csum/(nkp*M_PI);
               }
            }
            """
        ImBub=zeros((norb,norb),dtype=complex)
        dx=float(dx_)
        weave.inline(codeBub, ['ImBub','gkr','norb','dx','nkp','Qi','k_m_q','dOm','izero','om_idx'],type_converters=weave.converters.blitz, compiler='gcc')
        return ImBub


    #print 'Qlist=', Qlist
    ##########################
    # Reading real axis mesh #
    ##########################
    fm = open(filemesh, 'r')
    lined = fm.next().split()
    (mdelta, mmax) = map(float,lined[:2])
    Nd = int(lined[2])
    (oml,idxl) = LinLogMesh.LinLogMeshGen(mdelta,mmax,Nd)

    #print 'Qlist=', Qlist
    #######################################
    # Reading real axis Green's function  #
    #######################################
    # Reading some basic information from G_k
    fg = open(filegkr,'r')
    first_line = fg.next()
    nkp,nsymop,nom,cixdm,norbitals = map(int,first_line.split()[1:6])
    fg.close()
    (gkr,omr) = ReadGk(filegkr, nkp, nsymop, nom, cixdm)
    
    print 'shape(gkr)=', shape(gkr)
    if sum(abs(omr-oml))>1e-5: print 'Mesh in '+filegkr+' and in '+filemesh+' are not compatible.'
    
    zero_ind = oml.tolist().index(0.0)
    Oml=oml[zero_ind:]
    
    norb=cixdm

    ###############################
    # Computing real axis Bubble  #
    ###############################
    # This is the zero frequency value
    chi0r0 = zeros((len(Qlist),norb,norb), dtype=float)
    
    #print 'Qlist=', Qlist

    for iq,Q in enumerate(Qlist):
        print 'Q=', Q
        ImBub = zeros((norb,norb,len(Oml)),dtype=complex)
        for iOm in range(1,len(Oml)):
            ImBub[:,:,iOm] = Cmp_ChiQ_real(iOm,Q,gkr,k_m_q,nkp,norb,Oml[iOm]-Oml[iOm-1],idxl,Nd,zero_ind,k_index)
            
        ImBubr = real(ImBub)
        Bub = zeros((norb,norb,len(Oml)), dtype=complex)
        for iorb in range(norb):
            for jorb in range(norb):
                tOm = hstack( (-Oml[::-1][:-1], Oml[1:]) )
                ImB = hstack( (-ImBubr[iorb,jorb,::-1][:-1], ImBubr[iorb,jorb,1:]) )
                
                Bub[iorb,jorb,0] = integrate.trapz(ImB/tOm, x=tOm)/pi
                izero=len(tOm)/2
                for i in range(izero,len(tOm)):
                    Bub[iorb,jorb,i-izero+1] = krams.kramarskronig(ImB, tOm, i) + ImB[i]*1j
        
        chi0r0[iq,:,:] = Bub[:,:,0].real

        fq = open(fbubbleReal+str(iq), 'w')
        for iOm,Omx in enumerate(Oml):
            print >> fq, Omx,
            for iorb in range(norb):
                for jorb in range(norb):
                    print >> fq, Bub[iorb,jorb,iOm].real, Bub[iorb,jorb,iOm].imag,
            print >> fq
        fq.close()
        print 'Chi_q-real axis: Q point', iq, 'finished'
        
    return chi0r0

def ReadRealAxisBubble0(fbubbleReal, NQlist):
    chi0r0 = []
    for iq in range(NQlist):
        fq = open(fbubbleReal+str(iq), 'r')
        data = array(map(float, fq.next().split())[1::2])
        norb = sqrt(len(data))
        chi0r0.append( data.reshape(norb,norb) )
    return array(chi0r0)

def ReadRealAxisBubble(fbubbleReal, iq):
    data = loadtxt(fbubbleReal+str(iq))
    (nOm, nn) = shape(data)
    norb = int(sqrt((nn-1)/2))
    Oml = data[:,0]
    chi0r=zeros( (nOm, norb, norb), dtype=complex)
    for iOm in range(len(data)):
        chi0r[iOm,:,:] = (data[iOm,1::2]+data[iOm,2::2]*1j).reshape(norb,norb)
    return (chi0r,Oml)


def Print_2D(file_name, Obj_2D, mesh, dOm=0):
    if len(shape(Obj_2D))!=2 or shape(Obj_2D)[0]!=shape(Obj_2D)[1]:
        print "Wrong size of 2D_Object at "+file_name
        sys.exit(1)
    elif len(Obj_2D)!=len(mesh):
        print "Wrong size of 2D_mesh at "+file_name
        sys.exit(1)
    nom=len(mesh)
    fi=open(file_name,'w')
    for i in range(nom+dOm):
        print >>fi, "%.16f %.16f %.16f %.16f %.16f %.16f %.16f" %(mesh[i], Obj_2D[i,i].real, Obj_2D[i,i].imag, Obj_2D[i,nom-1+dOm-i].real, Obj_2D[i,nom-1+dOm-i].imag, Obj_2D[nom/2,i].real, Obj_2D[nom/2,i].imag)

def findex2(i1,i2,n1,n2):
    return i1*n2+i2

if __name__ == '__main__':
    
    smallShift=0
    if len(sys.argv)<2:
        print 'ERROR : need input filename'
        print 'The input file should contain: '
        print  'case.klist     # filename with k-list'
        print  'Qlist.dat      # filename with Qlist'
        print  'rmesh.dat      # real axis mesh'
        print  'G_k1r_         # file with real axis k-dependent Grens function'
        print  'G_local1r_     # file with real axis local Grens function'
        print  'chi0_real.     # name of the Bubble on real axis'
        print  'G_k1i_         # imaginary axis k-dependent Greens function'
        print  'G_local1i_     # imaginary axis local Greens function'
        print  'tvertex.dat    # ctqmc local vertex function'
        print  '100            # inverse temperature for bose function in Sq(omega)'
        sys.exit(1)
    fin = open(sys.argv[1], 'r')
    try:
        fklist      = fin.next().split()[0] # case.klist
        fQlist      = fin.next().split()[0] # case.qlist
        filemesh    = fin.next().split()[0] # rmesh.dat
        filegkr     = fin.next().split()[0] # G_k1r_
        fileglcr    = fin.next().split()[0] # G_local1r_
        fbubbleReal = fin.next().split()[0] # chi0_real.
        filegk      = fin.next().split()[0] # G_k1i_
        fileglc     = fin.next().split()[0] # G_local1i_
        fvertex     = fin.next().split()[0] # tvertex.dat
        beta_Sq     = float(fin.next().split()[0]) # 100
        line=fin.next()
    except StopIteration:
        line=''
    if line: smallShift = float(line.split()[0])
    
    print 'smallShift=', smallShift
    
    (kps, BS, BSI) = ReadKlist(fklist,True)
    Qlist = ReadKlist(fQlist,False)
    
    k_m_q = Get_k_m_q_index_CPP(kps)
    
    QCmpRealAxis = False
    
    if QCmpRealAxis:
        k_index = K_index(BSI,kps)
        chi0r0 = CmpRealAxisBubble(filemesh, filegkr, fbubbleReal, Qlist, k_index)
    else:
        chi0r0 = ReadRealAxisBubble0(fbubbleReal, len(Qlist))
        #for iq in range(len(Qlist)):
        #    print iq
        #    PrintM(chi0r0[iq])
    
    ##############################################
    # Reading green's function on imaginary axis #
    ##############################################
    # Reading some basic information from G_k
    fg = open(filegk,'r')
    first_line = fg.next()
    nkp,nsymop,nom,cixdm,norbitals = map(int,first_line.split()[1:6])
    second_line = fg.next()
    R_a = array(map(float,second_line.split()[1:1+3*norbitals]))
    R_a = R_a.reshape(norbitals,3)
    
    print 'R_a=', R_a
    print 'nkp,nom,cixdm=', nkp, nom, cixdm
    
    print 'ReadVertex'
    (beta,Nvfl,nomv,nOm,bfl_index,VertexH,VertexF)= ReadVertex(fvertex)
    
    print 'Reading Gloc on imaginary axis'
    (om,Gf) = ReadGlc(fileglc, 2*nomv, cixdm)
    
    print 'Reading Gk on imaginary axis'
    (gk,omi) = ReadGk(filegk, nkp, nsymop, nom, cixdm)
    
    gf = [ FermGF(Gf[i],nomv) for i in range(len(Gf))]
    
    ######################################
    # Computing Buble on imaginary axis  #
    ######################################
    Chi0_loc = Cmp_chi_0(gf,nOm)
    
    print 'bfl_index=', bfl_index
    
    orb_ud=[[] for _ in range(max(bfl_index)+1)]
    for i in range(len(bfl_index)): orb_ud[bfl_index[i]].append(i)
    print 'orb_ud=', orb_ud
    norb=len(orb_ud)
    dim = norb*2*nomv
    
    k_index = K_index(BSI,kps)
    
    fEU = open('EffU.dat','w')
    #####################################
    # Computing chi on imaginary axis   #
    #   and approximating it with U_eff #
    #####################################
    print 'Computing Chi on imaginary axis'
    for iOm in range(nOm-1,-1,-1):
        dOm=nOm-1-iOm
        (chi_S_loc, chi_C_loc) = Cmp_chi_loc(VertexH, VertexF, orb_ud, iOm, nomv)
        
        Chi_S_loc=zeros((dim,dim),dtype=complex)
        Chi_C_loc=zeros((dim,dim),dtype=complex)
        for i in range(norb):
            for j in range(norb):
                Chi_S_loc[2*nomv*i:2*nomv*(i+1),2*nomv*j:2*nomv*(j+1)]=chi_S_loc[i,j,:,:]
                Chi_C_loc[2*nomv*i:2*nomv*(i+1),2*nomv*j:2*nomv*(j+1)]=chi_C_loc[i,j,:,:]
        
        for im in range(dim):
            Chi_S_loc[im,im] +=Chi0_loc[iOm,im]
            Chi_C_loc[im,im] +=Chi0_loc[iOm,im]
        
        # chi = (chi0^{-1}+Gamma)^{-1}  hence
        # Gamma = chi^{-1} - chi0^{-1}
        GammaS=linalg.inv(Chi_S_loc)
        GammaC=linalg.inv(Chi_C_loc)
        for im in range(dim):
            GammaS[im,im] -= 1/Chi0_loc[iOm,im]
            GammaC[im,im] -= 1/Chi0_loc[iOm,im]

        
        omv=(2*(arange(2*nomv)-nomv)+1)*pi/beta
        for i1 in range(norb):
            for i2 in range(norb):
                fx = open('GammaS.'+str(i1)+'.'+str(i2),'w')
                fy = open('GammaC.'+str(i1)+'.'+str(i2),'w')
                for iom in range(2*nomv):
                    gs=GammaS[findex2(i1,iom,norb,2*nomv),findex2(i2,iom,norb,2*nomv)]
                    print >> fx, omv[iom], gs.real, gs.imag
                    gc=GammaC[findex2(i1,iom,norb,2*nomv),findex2(i2,iom,norb,2*nomv)]
                    print >> fy, omv[iom], gc.real, gc.imag
                fx.close()

        for im in range(dim):
            GammaS[im,im] += smallShift

        ###########???
        GT = zeros((dim,dim),dtype=complex)
        for i in range(dim):
            GT[:,i] = Chi0_loc[iOm,:]*GammaS[:,i]
        print 'chi0*G_local=', sorted(linalg.eigvals(GT).real)[:5]
        
        cs=zeros((norb,norb),dtype=float)
        ds=zeros((norb,norb),dtype=float)
        for i1 in range(norb):
            for i2 in range(norb):
                cs[i1,i2]=-GammaS[findex2(i1,nomv,norb,2*nomv),findex2(i2,nomv,norb,2*nomv)].real
                ds[i1,i2]= GammaC[findex2(i1,nomv,norb,2*nomv),findex2(i2,nomv,norb,2*nomv)].real
        print 'Gamma_Spin='
        PrintM(cs)
        print 'Gamma_Charge='
        PrintM(ds)
        ##########???
        
        Chi0_Q = zeros((dim,dim),dtype=complex)
        for iQ,Q in enumerate(Qlist):
            Qi = k_index(Q)
            chi0_Q = Cmp_chi0_Q(gk,nomv,dOm,norb,Qi,k_m_q)  # -G(iom) * G(iOm-iom)

            for i in range(norb):
                for j in range(norb):
                    for im in range(2*nomv):
                        Chi0_Q[2*nomv*i+im, 2*nomv*j+im]=chi0_Q[i,j,im]
            
            #  chi = (chi0^{-1}+Gamma)^{-1}
            Chi0_Q_inv = linalg.inv(Chi0_Q)
            Chi_S_Q = linalg.inv( GammaS + Chi0_Q_inv )
            Chi_C_Q = linalg.inv( GammaC + Chi0_Q_inv )



            GT = zeros((dim,dim),dtype=complex)
            GT = dot(Chi0_Q,GammaS)
            print iQ,Q,'chi0*G=', sorted(linalg.eigvals(GT).real)[:3]

            chiS_sum=zeros((norb,norb),dtype=complex)
            chi0_sum=zeros((norb,norb),dtype=complex)
            for i in range(norb):
                for j in range(norb):
                    chiS_sum[i,j]=sum(Chi_S_Q[i*2*nomv:(i+1)*2*nomv,j*2*nomv:(j+1)*2*nomv],axis=None)/beta
                    chi0_sum[i,j]=sum(Chi0_Q[i*2*nomv:(i+1)*2*nomv,j*2*nomv:(j+1)*2*nomv],axis=None)/beta
            #print Q
            #PrintM(chiS_sum*0.5) # spin-susceptibility with closing vertex
            #print
            #PrintM(chi0_sum*0.5)

            # The sum on imaginary axis is smaller than integral on the real axis, because we
            # use verly limited number of Matsubara points. We need to correct for that.
            # The correction is:
            #  correct == chi0(omega=0) * chi0(iOmega=0)^{-1}
            # We then correct chi on imaginary axis. Because there was double sum, we need to correct by
            #  chiS_sum -> correct * chiS_sum * correct
            # Finally, we define effective U, which is approximation for vertex Gamma at low energy by
            #  chi = (chi0^{-1}+ U_eff), hence   U_eff = chi^{-1}-chi0^{-1}
            #################################################################################
            BubQ_r  =  matrix(chi0r0[iQ,:,:])                           # chi0(Omega=0) on real axis
            correct =  matrix(chi0r0[iQ,:,:])*matrix(chi0_sum.real).I   # correction for finite Matsubara sum
            chiS_correct   = correct * matrix(chiS_sum.real) * correct  # corrected spin susceptibility
            Eff_U = chiS_correct.I - BubQ_r.I                           # Effective U
            # U = chi^{-1}-chi0r0^{-1}
            # ( chi0r^{-1} + U ) ^{-1} = chi

            
            print 'ChiS=', [chiS_sum[i,i].real for i in range(norb)]
            print 'ChiS_correct=', [chiS_correct[i,i] for i in range(norb)]
            #print 'Bubbl-real=', [BubQ_r[i,i] for i in range(norb)]
                
            # Reading real axis chi
            (chi0r, Oml) = ReadRealAxisBubble(fbubbleReal, iQ)
            chiQ = zeros(shape(chi0r), dtype=complex)
            SQ = zeros(len(Oml), dtype=float)
            for iOm in range(len(Oml)):
                chi = linalg.inv( linalg.inv(chi0r[iOm,:,:])+Eff_U )
                if iOm>0: SQ[iOm] = sum(chi, axis=None).imag/(1-exp(-beta_Sq*Oml[iOm]))
                chiQ[iOm,:,:] = chi


            print 'chiQ-final=', [chiQ[0,i,i].real for i in range(norb)]
                
            savetxt('Sqw.'+str(iQ) , vstack( (Oml[1:],SQ[1:]) ).transpose() )
            

            print >> fEU, iQ,
            for i in range(norb):
                for j in range(norb):
                    print >> fEU, Eff_U[i,j],
            print >> fEU
            
    fEU.close()
