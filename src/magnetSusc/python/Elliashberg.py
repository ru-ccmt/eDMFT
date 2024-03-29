#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule

from scipy import *
from scipy import linalg, integrate
import sys
import krams
import time
import scipy

from distutils.version import StrictVersion
if StrictVersion(scipy.__version__) > StrictVersion('0.19.0'):
    import weave
else:
    import scipy.weave as weave

def findex4(i1,i2,i3,i4,n1,n2,n3,n4):
    return ((i1*n2+i2)*n3+i3)*n4+i4
def findex3(i1,i2,i3,n1,n2,n3):
    return (i1*n2+i2)*n3+i3
def findex2(i1,i2,n1,n2):
    return i1*n2+i2

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

def Cmp_chi_0_F(gf,nOm):
    norb = len(gf)
    nomv = gf[0].nomv
    chi0=zeros((2*nOm-1,norb,norb,2*nomv),dtype=complex)
    for iOm in range(2*nOm-1):
        dOm=nOm-1-iOm
        for i in range(norb):
            for j in range(norb):
                for im in range(2*nomv):
                    chi0[iOm,i,j,im]= -gf[i][im]*gf[j][im+dOm]
    return chi0


def Cmp_chi_loc_H(VertexH, VertexF, orb_ud, iOm, nomv):
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

def Cmp_chi_loc_F(VertexH, VertexF, orb_ud, iOm, nomv):
    norb = len(orb_ud)
    ChiSC = zeros((norb,norb,2*nomv,2*nomv),dtype=complex)
    for iorb1,(up1,dn1) in enumerate(orb_ud):
        for iorb2,(up2,dn2) in enumerate(orb_ud):
            if iorb1!=iorb2:
                ChiSC[iorb1,iorb2] = 0.5*(VertexF[up1,up2,iOm]+VertexF[dn1,dn2,iOm]) 
    return ChiSC


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

def Cmp_chi0_Q2(gk,nomv,dOm,norb,Qi,k_m_q):
    
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
    for (int iorb1=0; iorb1<norb; iorb1++){
       for (int iorb2=0; iorb2<norb; iorb2++){
          for (int iorb3=0; iorb3<norb; iorb3++){
             for (int iorb4=0; iorb4<norb; iorb4++){
                for (int im=0; im<2*nomv; im++){
                   complex<double> csum=0;
                   for (int ik=0; ik<nkp; ik++){
                      int ikq=k_m_q(ik,Qi);
                      csum += -gkm(iorb3,iorb1,im,ik,  gk,nomv) * gkm(iorb2,iorb4,im+dOm,ikq, gk,nomv);
                   }
                   BubQ(iorb1,iorb2,iorb3,iorb4,im) = csum/static_cast<double>(nkp);
                }
             }
          }
       }
    }
    """
    nkp = len(k_m_q)
    if shape(gk)[3]!=nkp : print 'gk does not contain enough k-points', nkp, shape(gk)[3]
    BubQ=zeros((norb,norb,norb,norb,2*nomv),dtype=complex)
    weave.inline(codeBubQ, ['gk','nomv','dOm','norb','nkp','Qi','k_m_q','BubQ'],support_code=support_code,type_converters=weave.converters.blitz, compiler='gcc')
    return BubQ

def Cmp_chi0_Qpp(gk,nomv,norb,Qi,mQi,beta):
    " G_{a1,a3}(k,iom) * G_{a2,a4}(-k,-iom)"
    codeBubQ="""
    #line 324 "Elliashberg.py"
    using namespace std;
    for (int iorb1=0; iorb1<norb; iorb1++){
       for (int iorb2=0; iorb2<norb; iorb2++){
          for (int iorb3=0; iorb3<norb; iorb3++){
             for (int iorb4=0; iorb4<norb; iorb4++){
                complex<double> csum=0.;
                for (int im=0; im<nomv; im++){
                   csum += gk(iorb1,iorb3,im,Qi) * conj(gk(iorb4,iorb2,im,mQi));
                }
                BubQ(iorb1,iorb2,iorb3,iorb4) = csum.real()/beta;
             }
          }
       }
    }
    """
    nkp = shape(gk)[3]
    BubQ=zeros((norb,norb,norb,norb),dtype=float)
    weave.inline(codeBubQ, ['gk','nomv','norb','Qi','mQi','BubQ', 'beta'],type_converters=weave.converters.blitz, compiler='gcc')
    return BubQ

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

def fEnforceC4Symmetry(GammaS,symm,rest,nomv,norb):
    dim=shape(GammaS)[0]
    for ioma in range(2*nomv):
        for iomb in range(2*nomv):
            in0a=findex2(symm[0],ioma,norb,2*nomv)
            in0b=findex2(symm[0],iomb,norb,2*nomv)
            in1a=findex2(symm[1],ioma,norb,2*nomv)
            in1b=findex2(symm[1],iomb,norb,2*nomv)
            # Gamma[(2,oma),(2,omb)]==Gamma[(3,oma),(3,omb)]
            G0 = 0.5*(GammaS[in0a,in0b]+GammaS[in1a,in1b])
            GammaS[in0a,in0b]=G0
            GammaS[in1a,in1b]=G0
            # Gamma[(2,oma),(3,omb)]==Gamma[(3,oma),(2,omb)]
            G1 = 0.5*(GammaS[in0a,in1b]+GammaS[in1a,in0b])
            GammaS[in0a,in1b]=G1
            GammaS[in1a,in0b]=G1
            for j,i in enumerate(rest):
                ira=findex2(i,ioma,norb,2*nomv)
                irb=findex2(i,iomb,norb,2*nomv)
                # Gamma[(2,oma),(4,omb)]==Gamma[(3,oma),(4,omb)]
                G0 = 0.5*(GammaS[in0a,irb]+GammaS[in1a,irb])
                GammaS[in0a,irb]=G0
                GammaS[in1a,irb]=G0
                # Gamma[(4,omb),(2,oma)]==Gamma[(4,omb),(3,oma)]
                G1 = 0.5*(GammaS[irb,in0a]+GammaS[irb,in1a])
                GammaS[irb,in0a]=G1
                GammaS[irb,in1a]=G1
    return GammaS
    

def gEnforceC4Symmetry(GammaF,symm,rest,nomv,norb):
    cave = 0.5*(GammaF[symm[0],symm[1],:,:]+GammaF[symm[1],symm[0],:,:])
    GammaF[symm[0],symm[1],:,:] = cave
    GammaF[symm[1],symm[0],:,:] = cave
    for j,i in enumerate(rest):
        cave = 0.5*(GammaF[symm[0],i,:,:]+GammaF[symm[1],i,:,:])
        GammaF[symm[0],i,:,:] = cave
        GammaF[symm[1],i,:,:] = cave
        cave = 0.5*(GammaF[i,symm[0],:,:]+GammaF[i,symm[1],:,:])
        GammaF[i,symm[0],:,:] = cave
        GammaF[i,symm[1],:,:] = cave
    return GammaF

def CheckC4Symmetry_Gkf(Gf,gk, symm,rest,imk,kps):
    print 'Locally Diff=', sum(abs(Gf[symm[0],:]-Gf[symm[1],:]))
    for ik in range(len(kps)):
        k0 = kps[ik]
        k_rot = [ [k0[1],-k0[0],k0[2],k0[3]],  # rotating clockwise
                  [-k0[1],k0[0],k0[2],k0[3]] ] # rotating counter-clockwise
        ik_rot = [k_index(k_rot[0]),k_index(k_rot[1])]
        
        gxz = gk[symm[0],symm[0],0,ik]
        gyz = gk[symm[1],symm[1],0,ik_rot[0]]
        diff1 = abs((gxz-gyz)/(gxz+gyz))
        if abs(diff1)>5e-2:
            print 'diff-00-11:', ik, gxz, gyz, abs((gxz-gyz)/(gxz+gyz))
        #g01 = gk[symm[0],symm[1],0,ik]
        #g10 = gk[symm[1],symm[0],0,ik_rot[0]]
        #diff1 = abs((g01-g10)/(g01+g10))
        #if abs(diff1)>5e-2:
        #    print 'diff-01-10:', ik, g01, g10, abs((g01-g10)/(g01+g10))
        
        imkc =[ik_rot[j] for j in imk]
        for j,i in enumerate(rest):
            g1 = gk[symm[0],i,0,ik]
            g2 = gk[symm[1],i,0,imkc[j]]
            if abs(g1)>1e-5 or abs(g2)>1e-5 :
                diff2 = abs((g1-g2)/(g1+g2))
                if diff2>5e-2: print 'diff-00-i:', ik, i, g1, g2, abs((g1-g2)/(g1+g2)), diff2
            g1 = gk[i,symm[0],0,ik]
            g2 = gk[i,symm[1],0,imkc[j]]
            if abs(g1)>1e-5 or abs(g2)>1e-5:
                diff2 = abs((g1-g2)/(g1+g2))
                if diff2>5e-2:
                    print 'diff-i-00:', ik, i, g1, g2, abs((g1-g2)/(g1+g2)), diff2


def EnforceC4Symmetry_Gkf(Gf,gk, symm,rest,imk,kps):
    G_ave= 0.5*(Gf[symm[0],:]+Gf[symm[1],:])
    Gf[symm[0],:]=G_ave
    Gf[symm[1],:]=G_ave

    # Check C4 symmetry
    # gk(cixdm,cixdm,nom,nkp)
    for ik in range(len(kps)):
        k0 = kps[ik]
        k_rot = [ [k0[1],-k0[0],k0[2],k0[3]],  # rotating clockwise
                  [-k0[1],k0[0],k0[2],k0[3]] ] # rotating counter-clockwise
        ik_rot = [k_index(k_rot[0]), k_index(k_rot[1])]
        
        g_ave = 0.5*(gk[symm[0],symm[0],:,ik]+gk[symm[1],symm[1],:,ik_rot[0]])
        gk[symm[0],symm[0],:,ik] = g_ave
        gk[symm[1],symm[1],:,ik_rot[0]] = g_ave
        g_ave = 0.5*(gk[symm[0],symm[1],:,ik]+gk[symm[1],symm[0],:,ik_rot[0]])
        gk[symm[0],symm[1],:,ik] = g_ave
        gk[symm[1],symm[0],:,ik_rot[0]] = g_ave
        
        imkc =[ik_rot[j] for j in imk]
        for j,i in enumerate(rest):
            g1 = gk[symm[0],i,:,ik]
            g2 = gk[symm[1],i,:,imkc[j]]
            g_ave = 0.5*(g1+g2)
            gk[symm[0],i,:,ik] = g_ave
            gk[symm[1],i,:,imkc[j]] = g_ave
            g1 = gk[i,symm[0],:,ik]
            g2 = gk[i,symm[1],:,imkc[j]]
            g_ave = 0.5*(g1+g2)
            gk[i,symm[0],:,ik] = g_ave
            gk[i,symm[1],:,imkc[j]] = g_ave
    return (Gf, gk)


if __name__ == '__main__':
    
    smallShift=0.0
    #Symmetrize=True
    CheckTimeReversalSymmetry=False
    #CheckC4Symmetry=True
    #EnforceC4Symmetry=True
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
    
    fklist      = fin.next().split()[0] # case.klist
    fQlist      = fin.next().split()[0] # case.qlist
    #filemesh    = fin.next().split()[0] # rmesh.dat
    #filegkr     = fin.next().split()[0] # G_k1r_
    #fileglcr    = fin.next().split()[0] # G_local1r_
    #fbubbleReal = fin.next().split()[0] # chi0_real.
    filegk      = fin.next().split()[0] # G_k1i_
    fileglc     = fin.next().split()[0] # G_local1i_
    fvertex     = fin.next().split()[0] # tvertex.dat
    smallShift  = float(fin.next().split()[0]) # smallShift
    print 'smallShift=', smallShift
    exec(fin.next())  # CheckC4Symmetry
    exec(fin.next())  # EnforceC4Symmetry
    print 'CheckC4Symmetry=', CheckC4Symmetry
    print 'EnforceC4Symmetry=', EnforceC4Symmetry
    if CheckC4Symmetry or EnforceC4Symmetry:
        exec(fin.next()) # symm=[2,3]
        exec(fin.next()) # rest=[0,1,4]
        exec(fin.next()) # imk=[0,1,1]
        print 'symm=', symm
        print 'rest=', rest
        print 'imk=', imk
    
    (kps, BS, BSI) = ReadKlist(fklist,True)
    Qlist = ReadKlist(fQlist,False)
    
    k_m_q = Get_k_m_q_index_CPP(kps)
    
    print 'filegk=', filegk
    print 'fileglc=', fileglc
    
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
    
    k_index = K_index(BSI,kps)


    if CheckTimeReversalSymmetry:
        for ik in range(nkp):
            imk = k_index(-kps[ik])
            for i in range(cixdm):
                for j in range(cixdm):
                    if abs(gk[i,j,0,ik]-gk[j,i,0,imk])>1e-3:
                        print ik, i, j, gk[i,j,0,ik],gk[j,i,0,imk]
                    
    
    if CheckC4Symmetry:
        CheckC4Symmetry_Gkf(Gf,gk, symm,rest,imk,kps)
        
    
        # Enforces C4 symmetry for Green's function. In particularly, makes xz components equal to yz components
    if EnforceC4Symmetry:
        (Gf, gk) = EnforceC4Symmetry_Gkf(Gf,gk, symm,rest,imk,kps)
    
    
    ii=(omi[:]/om[0]-1)/2.
    for iw in range(len(ii)):
        if abs(ii[iw]-iw)>1e-4: break
    nw=iw
    
    gf = [ FermGF(Gf[i],nomv) for i in range(len(Gf))]
    
    ######################################
    # Computing Buble on imaginary axis  #
    ######################################
    Chi0_loc = Cmp_chi_0(gf,nOm)
    chi0_loc_F = Cmp_chi_0_F(gf,nOm)

    print 'bfl_index=', bfl_index
    
    orb_ud=[[] for _ in range(max(bfl_index)+1)]
    for i in range(len(bfl_index)): orb_ud[bfl_index[i]].append(i)
    print 'orb_ud=', orb_ud
    norb=len(orb_ud)
    dim = norb*2*nomv

    
    #for iom in range(2*nomv):
    #    in1=findex2(symm[0],iom,norb,2*nomv)
    #    in2=findex2(symm[1],iom,norb,2*nomv)
    #    c1=Chi0_loc[0,in1]
    #    c2=Chi0_loc[0,in2]
    #    if abs(c1)>1e-5 or abs(c2)>1e-5:
    #        if abs((c1-c2)/(c1+c2))>1e-5:
    #            print 'ERR-1:', iom,in1,in2,c1,c2
    #
    #    c1=chi0_loc_F[0,symm[0],symm[0],iom]
    #    c2=chi0_loc_F[0,symm[1],symm[1],iom]
    #    if abs(c1)>1e-5 or abs(c2)>1e-5:
    #        if abs((c1-c2)/(c1+c2))>1e-5:
    #            print 'ERR-2:', iom,in1,in2,c1,c2
    #
    #    for j,i in enumerate(rest):
    #        c1=chi0_loc_F[0,symm[0],i,iom]
    #        c2=chi0_loc_F[0,symm[1],i,iom]
    #        c3=chi0_loc_F[0,i,symm[0],iom]
    #        c4=chi0_loc_F[0,i,symm[1],iom]
    #
    #        if abs(c1)>1e-5 or abs(c2)>1e-5:
    #            if abs((c1-c2)/(c1+c2))>1e-5:
    #                print 'ERR-3:', iom,in1,in2,c1,c2
    #        if abs(c3)>1e-5 or abs(c4)>1e-5:
    #            if abs((c3-c4)/(c3+c4))>1e-5:
    #                print 'ERR-4:', iom,in1,in2,c3,c4
    #
    #sys.exit(0)
    
    
    omv=(2*(arange(2*nomv)-nomv)+1)*pi/beta
    #########################################

    fog = open('Chi0pp.dat','w')
    fos = open('Chi0pp_diagonal.dat','w')
    for iQ,Q in enumerate(Qlist):
        t1 = time.clock()
        Qi = k_index(Q)
        mQi = k_index(-Q)
        
        chi0pp = Cmp_chi0_Qpp(gk,nw,norb,Qi,mQi,beta)  # G_{a3,a1}(k,iom) * G*_{a2,a4}(-k,iom)
        
        t2 = time.clock()
        print 'iQ=', iQ, 'index(Q)=', Qi, 'index(-Q)=', mQi, 't[s]=',t2-t1
        
        print >> fos, iQ,
        for i1 in range(norb):
            print >> fos, chi0pp[i1,i1,i1,i1],
            for i2 in range(norb):
                for i3 in range(norb):
                    for i4 in range(norb):
                        print >> fog, chi0pp[i1,i2,i3,i4],
        print >> fog
        print >> fos
    fog.close()
    fos.close()
    
    #####################################
    # Computing chi on imaginary axis   #
    #####################################
    print 'Computing Chi on imaginary axis'
    iOm=0
    
    (chi_S_loc, chi_C_loc) = Cmp_chi_loc_H(VertexH, VertexF, orb_ud, iOm, nomv)
    
    Chi_S_loc=zeros((dim,dim),dtype=complex)
    Chi_C_loc=zeros((dim,dim),dtype=complex)
    for i in range(norb):
        for j in range(norb):
            Chi_S_loc[2*nomv*i:2*nomv*(i+1),2*nomv*j:2*nomv*(j+1)]=chi_S_loc[i,j,:,:]
            Chi_C_loc[2*nomv*i:2*nomv*(i+1),2*nomv*j:2*nomv*(j+1)]=chi_C_loc[i,j,:,:]
    
    # Adding the bubble, which is subtracted out in ctqmc for numerical reasons.
    # Need to add the bubble back.
    for im in range(dim):
        Chi_S_loc[im,im] +=Chi0_loc[iOm,im]
        Chi_C_loc[im,im] +=Chi0_loc[iOm,im]
    
    # Here we compute vertex for the case where the left-hand-side orbitals (incoming and outgoing)
    #  are both equal and the right hand side orbitals are equal (see picture):
    #   i1,om1         i2,om2
    #   -->------|---|--->--
    #   i1,om1-Om|   | i2,om2-Om
    #   --<------|---|---<--
    # chi = (chi0^{-1}+Gamma)^{-1}  hence
    # Gamma = chi^{-1} - chi0^{-1}
    GammaS=linalg.inv(Chi_S_loc)
    GammaC=linalg.inv(Chi_C_loc)
    #print 'GammaS_: ', GammaS[findex2(2,nomv,norb,2*nomv),findex2(2,nomv,norb,2*nomv)], GammaS[findex2(3,nomv,norb,2*nomv),findex2(3,nomv,norb,2*nomv)]
    for im in range(dim):
        GammaS[im,im] -= 1/Chi0_loc[iOm,im]
        GammaC[im,im] -= 1/Chi0_loc[iOm,im]

        GammaS[im,im] += smallShift
        
    #!!#!!
    if EnforceC4Symmetry:
        GammaS = fEnforceC4Symmetry(GammaS,symm,rest,nomv,norb)
        GammaC = fEnforceC4Symmetry(GammaC,symm,rest,nomv,norb)

    ##############??
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
    print 'Gamma_Spin(Ueff)='
    PrintM(cs)
    print 'Gamma_Charge(Ueff)='
    PrintM(ds)

    # Here we add the term for the case of unequal orbitals on the left and right-hand side
    #   i1,om1        i1,om2
    #   --->-----|---|--->---
    #   i2,om1-Om|   |i2,om2-Om
    #   ---<-----|---|---<---
    # We compute here for the case where we strictly have i1!=i2
    GammaF = zeros((norb,norb,2*nomv,2*nomv),dtype=complex)
    chi_loc_F = Cmp_chi_loc_F(VertexH, VertexF, orb_ud, iOm, nomv)
    
    #if EnforceC4Symmetry:
    #    chi_loc_F = gEnforceC4Symmetry(chi_loc_F,symm,rest,nomv,norb)
    
    for i in range(norb):
        for j in range(norb):
            if i!=j:
                # Add bubble to ctqmc-chi
                for im in range(2*nomv): chi_loc_F[i,j,im,im] += chi0_loc_F[iOm,i,j,im] 
                # Gamma=chi^-1-chi0^-1
                GammaF[i,j,:,:] = linalg.inv( chi_loc_F[i,j,:,:] )
                for im in range(2*nomv): GammaF[i,j,im,im] -= 1/chi0_loc_F[iOm,i,j,im]

    if EnforceC4Symmetry:
        GammaF = gEnforceC4Symmetry(GammaF,symm,rest,nomv,norb)
        
    cs=zeros((norb,norb),dtype=float)
    for i1 in range(norb):
        for i2 in range(norb):
            cs[i1,i2]=-GammaF[i1,i2,nomv,nomv].real
    print 'Gamma_Fock(Ueff)='
    PrintM(cs)
    

    # Here we combine the above two terms which contribute to the vertex. They are both functions of two
    # orbitals, but once combined, we can only express the sum in terms of four dimensional tensor. 
    # Our definition is:
    # 
    #  Gamma(i1,i2,om1, i3,i4, om2)==
    #
    #  --->-----|-----|--->----
    #  i1,om1   |Gamma| i3,om2
    #  ---<-----|-----|--<-----
    #  i2,om1-Om        i4,om2-Om
    #  
    # There are two kinds of Gamma: the magnetic and the charge Gamma.
    # We break the spin-orbital index i1 into orbital b1 and spin s1 : i1=(b1,s1)
    # For magnetic,  we have GammaM = Gamma((b1,up),(b2,up),(b3,up),(b4,up)) - Gamma((b1,up),(b2,up),(b3,dn),(b4,dn))
    # and for charge we have GammaC = Gamma((b1,up),(b2,up),(b3,up),(b4,up)) + Gamma((b1,up),(b2,up),(b3,dn),(b4,dn))
    dim2 = norb*norb*2*nomv
    Gamma_QS = zeros((dim2,dim2),dtype=complex) 
    Gamma_QC = zeros((dim2,dim2),dtype=complex) 
    for i1 in range(norb):
        for i3 in range(norb):
            for im in range(2*nomv):
                for jm in range(2*nomv):
                    Gamma_QS[findex3(i1,i1,im,norb,norb,2*nomv),findex3(i3,i3,jm,norb,norb,2*nomv)]=GammaS[findex2(i1,im,norb,2*nomv),findex2(i3,jm,norb,2*nomv)]
                    Gamma_QC[findex3(i1,i1,im,norb,norb,2*nomv),findex3(i3,i3,jm,norb,norb,2*nomv)]=GammaC[findex2(i1,im,norb,2*nomv),findex2(i3,jm,norb,2*nomv)]

    # Adding the second term where uniqual orbitals are found on the left-hand side.
    for i1 in range(norb):
        for i2 in range(norb):
            if i1!=i2:
                for im in range(2*nomv):
                    for jm in range(2*nomv):
                        Gamma_QS[findex3(i1,i2,im,norb,norb,2*nomv),findex3(i1,i2,jm,norb,norb,2*nomv)]+=GammaF[i1,i2,im,jm]
                        Gamma_QC[findex3(i1,i2,im,norb,norb,2*nomv),findex3(i1,i2,jm,norb,norb,2*nomv)]+=GammaF[i1,i2,im,jm]
    
    Gamma_QS_inv = linalg.inv( Gamma_QS )
    Gamma_QC_inv = linalg.inv( Gamma_QC )

    cs=zeros((norb,norb),dtype=float)
    cd=zeros((norb,norb),dtype=float)
    for i1 in range(norb):
        for i2 in range(norb):
            cs[i1,i2]=-Gamma_QS_inv[findex3(i1,i1,nomv,norb,norb,2*nomv),findex3(i2,i2,nomv,norb,norb,2*nomv)].real
            cd[i1,i2]=Gamma_QC_inv[findex3(i1,i1,nomv,norb,norb,2*nomv),findex3(i2,i2,nomv,norb,norb,2*nomv)].real
    print 'Gamma_PH_Spin='
    PrintM(cs)
    print 'Gamma_PH_Charge='
    PrintM(cd)
    
    Chi0_Q = zeros((dim2,dim2),dtype=complex)
    print 'dim2=', dim2
    t_chi0=0
    t_inv=0

    print 'Gamma_QS:', Gamma_QS[findex3(2,2,nomv,norb,norb,2*nomv),findex3(2,2,nomv,norb,norb,2*nomv)], Gamma_QS[findex3(3,3,nomv,norb,norb,2*nomv),findex3(3,3,nomv,norb,norb,2*nomv)]
    print 'Gamma_QS_inv:', Gamma_QS_inv[findex3(2,2,nomv,norb,norb,2*nomv),findex3(2,2,nomv,norb,norb,2*nomv)], Gamma_QS_inv[findex3(3,3,nomv,norb,norb,2*nomv),findex3(3,3,nomv,norb,norb,2*nomv)]

    fopp = open('Gpp.dat','w')
    fopm = open('Gpm.dat','w')
    fomm = open('Gmm.dat','w')
    for iOm in range(nOm-1,-1,-1):
        dOm=nOm-1-iOm
        for iQ,Q in enumerate(Qlist):
            t1 = time.clock()
            Qi = k_index(Q)
            print 'iOm=', iOm, 'iQ=', iQ
            
            Chi0_Q[:,:] = 0
            
            chi0_Q = Cmp_chi0_Q2(gk,nomv,dOm,norb,Qi,k_m_q)  # -G_{a3,a1}(iom) * G_{a2,a4}(iom-iOm)
            
            for i1 in range(norb):
                for i2 in range(norb):
                    for i3 in range(norb):
                        for i4 in range(norb):
                            for im in range(2*nomv):
                                Chi0_Q[findex3(i1,i2,im,norb,norb,2*nomv),findex3(i3,i4,im,norb,norb,2*nomv)]=chi0_Q[i1,i2,i3,i4,im]

            t2 = time.clock()
            #G*(chi0^-1+G)^{-1}*G
            # #  chi = (chi0^{-1}+Gamma)^{-1}
            # Chi0_Q_inv = linalg.inv(Chi0_Q)
            # #print 'inverse2'
            # Chi_QS = linalg.inv( Gamma_QS + Chi0_Q_inv ) 
            # #print 'double products start'
            # Gpp_QS = dot(dot(Gamma_QS, Chi_QS), Gamma_QS)
            # print 'double products finished'
            # Chi_QC = linalg.inv( Gamma_QC + Chi0_Q_inv )
            # Gpp_QC = dot(dot(Gamma_QC, Chi_QC), Gamma_QC)
            # Gpp_QS = Gamma_QS - linalg.inv( Gamma_QS_inv+Chi0_Q )  # equivalent to -G*chi*G
            # Gpp_QC = Gamma_QC - linalg.inv( Gamma_QC_inv+Chi0_Q )  # equivalent to -G*chi*G
            # Gpp = Gpp_QS*1.5 - Gpp_QC*0.5 + 0.5*Gamma_QS + 0.5*Gamma_QC

            Gpp_QS = linalg.inv( Gamma_QS_inv+Chi0_Q )  # equivalent to G-G*chi*G
            Gpp_QC = linalg.inv( Gamma_QC_inv+Chi0_Q )  # equivalent to G-G*chi*G
            Gpp = -Gpp_QS*1.5 + Gpp_QC*0.5 + 2*Gamma_QS # equivalent to 0.5*Gs+0.5*Gc+1.5*Gs*Chis*Gs-0.5*Gc*Chic*Gc
            
            
            t3 = time.clock()
            
            print 'Q=', Q, 'Omega=0:'
            cd=zeros((norb,norb),dtype=float)
            for i1 in range(norb):
                for i2 in range(norb):
                    cd[i1,i2]=Gpp[findex3(i1,i1,nomv,norb,norb,2*nomv),findex3(i2,i2,nomv,norb,norb,2*nomv)].real
            PrintM(cd)
            
            for i1 in range(norb):
                for i2 in range(norb):
                    for i3 in range(norb):
                        for i4 in range(norb):
                            #im=nomv
                            #print >> fopp, Gpp[(i1*norb+i2)*2*nomv+im,(i3*norb+i4)*2*nomv+im].real,
                            i1i2p = findex3(i1,i2,nomv,  norb,norb,2*nomv)
                            i1i2m = findex3(i1,i2,nomv-1,norb,norb,2*nomv)
                            i3i4p = findex3(i3,i4,nomv,  norb,norb,2*nomv)
                            i3i4m = findex3(i3,i4,nomv-1,norb,norb,2*nomv)
                            print >> fopp, Gpp[i1i2p, i3i4p].real,
                            print >> fopm, Gpp[i1i2p, i3i4m].real,
                            print >> fomm, Gpp[i1i2m, i3i4m].real,
            print >> fopp
            print >> fopm
            print >> fomm
            t_chi0 += t2-t1
            t_inv += t3-t2
            
            print 't(chi0)=', t_chi0, ' t(inv)=', t_inv, ' t(total)=', t_chi0+t_inv
    fopp.close()
    fopm.close()
    fomm.close()

