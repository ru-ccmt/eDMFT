#!/usr/bin/env python
from pylab import *
from numpy import *
import numpy as npy
from scipy import linalg
from itertools import chain, product
from timeit import default_timer as timer
import sys
import builtins

def PrintM(A):
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            print("%12.8f %12.8f   " % (A[i,j].real, A[i,j].imag), end=' ')
        print()

def OrthogonalizeDiag(UmA, Olap):
    olp,vlp = linalg.eigh(-Olap)  # Note: vlp.H * Olp * vlp = olp
    olp *= -1
    print('Corresponding eigensystem of the overlap=', olp)
    #check = dot(dot(vlp.T.conj(), Olap), vlp)
    #print 'check=', check
    # V = vlp * 1/sqrt(o) * vlp.H
    V = zeros(shape(Olap),dtype=complex)
    V[:,:] = vlp[:,:]
    for i in range(len(olp)): V[:,i] *= 1./sqrt(olp[i])
    V = dot(V,vlp.T.conj())
    print('Transformation=')
    PrintO(V.T)
    Uma_new =  dot(V, UmA)
    return (Uma_new, V)

def OrthogonalizeDiag2(UmA, which, Print=True):
    U, s, Vh = linalg.svd(UmA)
    if Print:
        print('singular values '+which+'=', s)
    if min(abs(s))<1e-7:
        print('ERROR : One the functions is linear dependent, hence you have a bad basis. Please change LocFunc_...', end=' ')
    #print 'shape(U)=', shape(U), 'shape(Vh)=', shape(Vh)
    Uma_new = dot(U, Vh[:len(s),:])
    #print 'Vh=', [abs(dot(Vh.T.conj()[:,:len(s)], Vh[:len(s),:]))[l,l] for l in range(len(Vh))]
    return Uma_new


def gram_schmidt(Olap):
    N = len(Olap)
    alpha = zeros((N,N), dtype=complex)
    alpha[0,0] = 1/sqrt(abs(Olap[0,0]))
    for i in range(1,N):
        m = zeros(i+1, dtype=complex)
        b = dot(dot(alpha[:i,:i],alpha[:i,:i].T.conj()),Olap[:i,i])
        m[i] = 1.
        m[:i] = -b[:]
        d2 = dot(dot(conj(m),Olap[:(i+1),:(i+1)]),m)
        alpha[:(i+1),i] = m[:]/sqrt(abs(d2))
    return alpha.T

def OrthogonalizeGramSch(Uma, Olap):
    alpha = gram_schmidt(Olap)
    print('The transformation matrix according to GramSchmidt=')
    print(alpha.tolist())
    Uma_new = dot(alpha, Uma)
    return Uma_new
    
def PrintO(Olp):
      for i in range(len(Olp)):
        print("%12.8f "*len(Olp[i]) % tuple(abs(Olp[i,:])))
      
if __name__ == '__main__':
      #Nr = 200      # how many real space points. They will be in the interval [-Nr/2..Nr/2]*Rm for one sheet of graphene, and [-Nr/2..Nr/2]*Rp for the other sheet.
      #Ndivide = 4   # how many points inside single unit cell of original graphene
      #how_many = 1  # how many shells of emergent lattice to keep in the plot
      #optimal_r = 5.  # width of gaussian is exp(- (5/|Rs| * |r|)^2 )
      BA_region_needed = True
      exec(open("tbg_params.py").read())
      for p in list(par_pr.keys()):
        ss = p + '=' + str(par_pr[p])
        print(ss)
        exec(ss)
      unit = par_bs['unit']
      band=''
      if par_bs['PlotBands']:
            band = '_band'
      Eks     = load('ene'+band+'.npy')
      Psis    = load('vector'+band+'.npy')
      npzfile = load('basis'+band+'.npz')
      
      print('shape(Psis)=', shape(Psis))
      
      k_basis, Rm, Rp, Rs = npzfile['k_basis'], npzfile['Rm'], npzfile['Rp'], npzfile['Rs']
      rBAm, rBAp, bm, bp, bs, b0 = npzfile['rBAm'], npzfile['rBAp'], npzfile['bm'], npzfile['bp'], npzfile['bs'], npzfile['b0']
      Nm, Np = npzfile['Nm'], npzfile['Np']

      Rs_BA = (Rs[0]+Rs[1])/3.
      print('real space vectors of upper layer graphene, Rm=', Rm.tolist())
      print('real space vectors of lower layer graphene, Rp=', Rp.tolist())
      print('real space vectors of the emergent unit-c., Rs=', Rs.tolist())
      print('reciprocal vectors of upper layer graphene, bm=', bm.tolist())
      print('reciprocal vectors of lower layer graphene, bp=', bp.tolist())
      print('reciprocal vectors of the emergent unit-c., bs=', bs.tolist())
      print('Nm=', Nm.tolist(), 'allclose(Nm*bs,bm)=', allclose(Nm@bs,bm), 'allclose(Nm.T*Rm,Rs)=', allclose(Nm.T @ Rm,Rs))
      print('Np=', Np.tolist(), 'allclose(Np*bs,bp)=', allclose(Np@bs,bp), 'allclose(Np.T*Rp,Rs)=', allclose(Np.T @ Rp,Rs))
      print('vector between A and B site of the emergent unit cell Rs_BA=',Rs_BA.tolist())
      print('R_{BA}-R_{AA}=Rs_BA=', Rs_BA)
      print('Rs_BA in terms of Rm unit cell vectors, rr=', Rs_BA@bm.T/(2*pi) )
      print('Rs_BA in terms of Rp unit cell vectors, qq=', Rs_BA@bp.T/(2*pi) )
      print('vector from A to B sublattice on upper layer, rBAm=', rBAm.tolist())
      print('vector from A to B sublattice on lower layer, rBAp=', rBAp.tolist())
      #
      # These are all Dirac cones in the original graphene for upper layer
      Km = array([(2*bm[0]+bm[1])/3.,(-bm[0]+bm[1])/3.,-(2*bm[1]+bm[0])/3.])
      # These are all Dirac cones in the original graphene on lower layer
      Kp = array([(2*bp[0]+bp[1])/3.,(-bp[0]+bp[1])/3.,-(2*bp[1]+bp[0])/3.])
      print('Location of Dirac cones in the upper layer=', Km.tolist())
      print('Location of Dirac cones in the lower layer=', Kp.tolist())
      # Construsting local functions as a good guess
      LcF_Km = [Km,-Km] # all 6 Dirac cones in the upper layer, not just the first three unique
      LcF_Kp = [Kp,-Kp] # all 6 Dirac cones in the lower layer, not just the firts three unique
      #
      # Explanation of the algorithm:
      #             top-layer               bottom-layer
      # phi1 = exp( i*Km_j*(r-R_A)) - i * exp( i*Kp_j*(r-R_A))
      # phi3 = exp(-i*Km_j*(r-R_A)) + i * exp(-i*Kp_j*(r-R_A))
      # phi2 = exp( i*Km_j*(r-R_B)) + i * exp( i*Kp_j*(r-R_B))
      # phi4 = exp(-i*Km_j*(r-R_B)) - i * exp(-i*Kp_j*(r-R_B))
      ### here Km_j (Kp_j) is the Dirac vector on top (bottom) layer, and j runs from 0..2 for three vectors.
      #
      # We first compute LocFunc_phasem, which contains only the second part of the exponent:
      # LocFunc_phasem[0] = exp( i*Km_j*(-R_A)); LocFunc_phasep[0] = -i*exp( i*Kp_j*(-R_A))
      # LocFunc_phasem[1] = exp(-i*Km_j*(-R_A)); LocFunc_phasep[1] = +i*exp(-i*Kp_j*(-R_A))
      # LocFunc_phasem[2] = exp( i*Km_j*(-R_B)); LocFunc_phasep[2] = +i*exp( i*Kp_j*(-R_B))
      # LocFunc_phasem[3] = exp(-i*Km_j*(-R_B)); LocFunc_phasep[3] = -i*exp(-i*Kp_j*(-R_B))
      #
      # We later compute the more expensive exponents for all real space vectors within first unit cell
      # LF_rmA[0,r_A] = exp( i*Km_j*(r_A-R_A)); LF_rpA[0,r_A] = -i*exp( i*Kp_j*(r_A-R_A))
      #
      # phi1[r_A] = LF_rmA[0,r_A] +  LF_rpA[0,r_A]
      # phi1[r_B] = LF_rmB[0,r_B] +  LF_rpB[0,r_B]
      # phi3[r_A] = LF_rmA[1,r_A] +  LF_rpA[1,r_A]
      # phi3[r_B] = LF_rmB[1,r_B] +  LF_rpB[1,r_B]
      # phi2[r_A] = LF_rmA[2,r_A] +  LF_rpA[2,r_A]
      # phi2[r_B] = LF_rmB[2,r_B] +  LF_rpB[2,r_B]
      # phi4[r_A] = LF_rmA[3,r_A] +  LF_rpA[3,r_A]
      # phi4[r_B] = LF_rmB[3,r_B] +  LF_rpB[3,r_B]

      N_loc_func = 4
      LocFunc_Km = zeros((4,3,2))
      LocFunc_phasem = zeros((4,3),dtype=complex)
      R_AB = [-rBAm/2,rBAm/2] # R_AB[0] (R_AB[1]) is vector to A (B) site on the top layer. Note that rBAm is vector between A&B subplatices, and origin is set to be in-between
      ssgn = [1,1,1,1,1,1,1,1]
      for i in range( N_loc_func ):
          R = R_AB[(int(i/2))%2]
          KS = LcF_Km[i % 2]        # the three K vectors for m-layer that host Dirac cone. when i=0 or 1, we have the other set of three vectors
          for j,k in enumerate(KS): # over the three vectors, either + or -
              LocFunc_Km[i,j,:] = k   # remember Km
              LocFunc_phasem[i,j] = ssgn[i] * exp(-1j*dot(R,k))
              # LocFunc_phasem[i=0,j=0..2]= exp(-i*R_A*Km_j)
              # LocFunc_phasem[i=1,j=0..2]= exp(+i*R_A*Km_j)
              # LocFunc_phasem[i=2,j=0..2]= exp(-i*R_B*Km_j)
              # LocFunc_phasem[i=3,j=0..2]= exp(+i*R_B*Km_j)
      
      LocFunc_Kp = zeros((4,3,2))
      LocFunc_phasep = zeros((4,3),dtype=complex)
      R_AB = [-rBAp/2,rBAp/2]
      ssgn = [-1j,1j,1j,-1j,1j,-1j,-1j,1j]
      for i in range( N_loc_func ):
          R = R_AB[(int(i/2))%2]
          KS = LcF_Kp[i % 2]
          for j,k in enumerate(KS):
              LocFunc_Kp[i,j,:] = k
              LocFunc_phasep[i,j] = ssgn[i] * exp(-1j*dot(R,k))
              # LocFunc_phasep[i=0,j=0..2]= -i*exp(-i*R_A*Kp_j)
              # LocFunc_phasep[i=1,j=0..2]=  i*exp(+i*R_A*Kp_j)
              # LocFunc_phasep[i=2,j=0..2]=  i*exp(-i*R_B*Kp_j)
              # LocFunc_phasep[i=3,j=0..2]= -i*exp(+i*R_B*Kp_j)

      print('LocFunc_phasem=')
      for i in range(N_loc_func):
        print(', '.join([f"{c.real:11.8f} + {c.imag:11.8f}j" for c in LocFunc_phasem[i,:]]))
      print('LocFunc_phasep=')
      for i in range(N_loc_func):
        print(', '.join([f"{c.real:11.8f} + {c.imag:11.8f}j" for c in LocFunc_phasep[i,:]]))


      _R0_ = linalg.inv(array([[ 1.0,1.0/sqrt(3.)], [-1.0,1.0/sqrt(3.)]]).T)  # direct vectors of original unrotated graphene
      ds = 0.5*(linalg.norm(Rs[0])+linalg.norm(Rs[1]))  # how long is the unit cell vector of the emergent unit cell
      # distances between neighboring AA sites
      print('distances between AA sites in emergent lattice, ds=', ds)
      Ris = array([[sqrt(3.)/2,0.5],[0.,1.],[-sqrt(3.)/2,0.5]])*2./ds
      # r @ Ris.T measures how many emergent unit cells away is point r.
      print('(R0*ds/sqrt(3))*Ris.T=', ((_R0_ @ Ris.T)*ds/sqrt(3.)).tolist())
      
      Rm_AB = [-rBAm/2,rBAm/2]
      Rp_AB = [-rBAp/2,rBAp/2]
      # R_ij == long list of integer points (i,j), where (i,j) \in [-Nr/2,Nr/2]
      # for example: R_ij  = array(list(product(range(-int(Nr/2),int(Nr/2)+1),range(-int(Nr/2),int(Nr/2)+1))))
      R_ij = array(meshgrid(arange(-Nr//2, Nr//2+1),arange(-Nr//2,Nr//2+1))).T.reshape(-1, 2)
      #
      rmA = R_ij @ Rm + Rm_AB[0]  # points on the A sublattice on the upper sheet of graphene
      rmB = R_ij @ Rm + Rm_AB[1]  # points on the B sublattice on the upper sheet of graphene
      rpA = R_ij @ Rp + Rp_AB[0]  # points on the A sublattice on the lower sheet of graphene
      rpB = R_ij @ Rp + Rp_AB[1]  # points on the B sublattice on the lower sheet of graphene
      
      _lengths_ = abs(rmA @ Ris.T) # check how far away are rmA points
      rmA = rmA[ (_lengths_[:,0]<how_many) & (_lengths_[:,1]<how_many) & (_lengths_[:,2]<how_many) ] # only keep those not too far
      _lengths_ = abs(rmB @ Ris.T) # check how far away are rmB points
      rmB = rmB[ (_lengths_[:,0]<how_many) & (_lengths_[:,1]<how_many) & (_lengths_[:,2]<how_many) ] # only keep those not too far
      _lengths_ = abs(rpA @ Ris.T) # check how far away are rpA points
      rpA = rpA[ (_lengths_[:,0]<how_many) & (_lengths_[:,1]<how_many) & (_lengths_[:,2]<how_many) ] # only keep those not too far
      _lengths_ = abs(rpB @ Ris.T) # check how far away are rpB points
      rpB = rpB[ (_lengths_[:,0]<how_many) & (_lengths_[:,1]<how_many) & (_lengths_[:,2]<how_many) ] # only keep those not too far
      
      print('Number of real space vectors on upper sheet to A sublattice r_mA=', len(rmA), 'to B sublattice r_mB=', len(rmB))
      print('Number of real space vectors on lower sheet to A sublattice r_pA=', len(rpA), 'to B sublattice r_pB=', len(rpB))
      
      # now find minimum and maximum coordinate of all those arrays
      all_arrays = concatenate((rmA, rmB, rpA, rpB), axis=0)
      rmin = npy.min(all_arrays, axis=0)
      rmax = npy.max(all_arrays, axis=0)
      L = rmax-rmin
      Nxy = array(npy.round(L)*Ndivide,dtype=int)
      print('min(r)=', rmin, 'max(r)=', rmax, 'L=', L, 'Nxy=', Nxy)
      
      r_range = optimal_r/norm(Rs)
      Nfunc = len(LocFunc_Km)
      print('Number of DMFT projectors we construct, Nfunc=', Nfunc)

      LF_rmA = zeros( (Nfunc, len(rmA) ), dtype=complex ) # AA - local function
      LF_rmB = zeros( (Nfunc, len(rmB) ), dtype=complex ) # AA - local function
      LF_rpA = zeros( (Nfunc, len(rpA) ), dtype=complex ) # AA - local function
      LF_rpB = zeros( (Nfunc, len(rpB) ), dtype=complex ) # AA - local function
      if BA_region_needed:
          LF_qmA = zeros( (Nfunc, len(rmA) ), dtype=complex ) # BA - local function
          LF_qmB = zeros( (Nfunc, len(rmB) ), dtype=complex ) # BA - local function
          LF_qpA = zeros( (Nfunc, len(rpA) ), dtype=complex ) # BA - local function
          LF_qpB = zeros( (Nfunc, len(rpB) ), dtype=complex ) # BA - local function
      for ifc in range(Nfunc):
          if BA_region_needed:
              # note that we slightly move the origin on AB region from Rs_BA => Rs_BA-rBAm and Rs_BA => Rs_BA+rBAp for the two sheets.
              # This is so that the local functions for both regions have the same form, otherwise we would need to adjust the functions for AB region,
              # which would be different from AA region. Note that below, where we project to band-structure, we keep correct definition of Rs_BA. So this is
              # just a convenient way of redefiniting the local functions.
              kRm_BA = LocFunc_Km[ifc,:,:] @ (Rs_BA - rBAm)   # kR_BA[ikc] -- for functions centered on BA instead of AA regions
              LocFunc_phasem_BA = LocFunc_phasem[ifc,:]*exp(-kRm_BA*1j) # LocFunc_phase_BA[ikc]
              kRp_BA = LocFunc_Kp[ifc,:,:] @ (Rs_BA + rBAp)   # kR_BA[ikc] -- for functions centered on BA instead of AA regions
              LocFunc_phasep_BA = LocFunc_phasep[ifc,:]*exp(-kRp_BA*1j) # LocFunc_phase_BA[ikc]

          nr = (linalg.norm(rmA,axis=1)*r_range)**2
          gauss = exp(-nr)                               # gauss[ir] = exp(-|rmA[ir]|^2 * r_range^2)
          exp_kr = exp( LocFunc_Km[ifc,:,:]@rmA.T * 1j ) # exp_kr[Km_j,ir] = exp(Km_j*r_mA *i) = exp( LocFunc_K[ifc,Km_j,:]*r.T[:,ir] *i)
          LF_rmA[ifc,:] = LocFunc_phasem[ifc,:] @ exp_kr # LF_rmA[ifc,ir] = \sum_j exp(Km_j*(r_mA-R_[AorB]) *i) = \sum_j LocFunc_phasem[ifc,Km_j]*exp( kr[Km_j,ir] *i )
          LF_rmA[ifc,:] *= gauss # adding gaussian
          if BA_region_needed:
              # BA region
              LF_qmA[ifc,:] = LocFunc_phasem_BA @ exp_kr # LF[ifc,ir] = LocFunc_phase[ifc,ikc] * exp(kR_BA*i) * exp( kr[ikc,ir]*i)
              LF_qmA[ifc,:] *= gauss
          
          nr = (linalg.norm(rmB,axis=1)*r_range)**2
          gauss = exp(-nr)                               # gauss[ir] = exp(-|rmB[ir]|^2 * r_range^2)
          exp_kr = exp( LocFunc_Km[ifc,:,:]@rmB.T * 1j ) # exp_kr[Km_j,ir] = exp(Km_j*r_mB *i) = exp( LocFunc_Km[ifc,ikc,:]*r.T[:,ir] *i )
          LF_rmB[ifc,:] = LocFunc_phasem[ifc,:] @ exp_kr # LF_rmB[ifc,ir] = \sum_j exp(Km_j*(r_mB-R_[AorB])) =\sum_j LocFunc_phase[ifc,Km_j]*exp( kr[Km_j,ir] *i)
          LF_rmB[ifc,:] *= gauss  # adding gaussian
          if BA_region_needed:
              # BA region
              LF_qmB[ifc,:] = LocFunc_phasem_BA @ exp_kr
              LF_qmB[ifc,:] *= gauss
          
          nr = (linalg.norm(rpA,axis=1)*r_range)**2
          gauss = exp(-nr)                               # gauss[ir] = exp(-|rpA[ir]|^2 * r_range^2)
          exp_kr = exp( LocFunc_Kp[ifc,:,:]@rpA.T * 1j ) # exp_kr[Kp_j,ir] = exp(Kp_j*r_pA *i) = exp( LocFunc_Kp[ifc,Kp_j,:]*r.T[:,ir] *i )
          LF_rpA[ifc,:] = LocFunc_phasep[ifc,:] @ exp_kr # LF_rpA[ifc,ir] = \sum_j (+-i)*exp(Kp_j*(r_pA-R_[AorB]) *i) = \sum_j LocFunc_phase[ifc,Kp_j]*exp( kr[Kp_j,ir] *i)
          LF_rpA[ifc,:] *= gauss  # adding gaussian
          if BA_region_needed:
              # BA region
              LF_qpA[ifc,:] = LocFunc_phasep_BA @ exp_kr
              LF_qpA[ifc,:] *= gauss
          
          nr = (linalg.norm(rpB,axis=1)*r_range)**2
          gauss = exp(-nr)                               # gauss[ir] = exp(-|rpB[ir]|^2 * r_range^2)
          exp_kr = exp( LocFunc_Kp[ifc,:,:]@rpB.T * 1j ) # exp_kr[Kp_j,ir] = exp(Kp_j*r_pB *i) = exp( LocFunc_K[ifc,Kp_j,:]*r.T[:,ir] *i )
          LF_rpB[ifc,:] = LocFunc_phasep[ifc,:] @ exp_kr # LF_rpB[ifc,ir] = \sum_j (+-i)*exp(Kp_k*(r_pB-R_[AorB]) *i) = \sum_j LocFunc_phase[ifc,Kp_j]*exp( kr[Kp_j,ir] *i)
          LF_rpB[ifc,:] *= gauss
          if BA_region_needed:
              # BA region
              LF_qpB[ifc,:] = LocFunc_phasep_BA @ exp_kr
              LF_qpB[ifc,:] *= gauss
          
          nmA = LF_rmA[ifc,:] @ LF_rmA[ifc,:].conj()
          nmB = LF_rmB[ifc,:] @ LF_rmB[ifc,:].conj()
          npA = LF_rpA[ifc,:] @ LF_rpA[ifc,:].conj()
          npB = LF_rpB[ifc,:] @ LF_rpB[ifc,:].conj()
          normc = sqrt( nmA.real + nmB.real + npA.real + npB.real )
          if (abs(normc)<1e-6):
              print('ERROR : normc[ifc='+str(ifc)+']=0, hence the choosen function is singular, normc=', normc)
              sys.exit(1)
          LF_rmA[ifc,:] *= 1./normc
          LF_rmB[ifc,:] *= 1./normc
          LF_rpA[ifc,:] *= 1./normc
          LF_rpB[ifc,:] *= 1./normc

          if BA_region_needed:
              nmA = LF_qmA[ifc,:] @ LF_qmA[ifc,:].conj()
              nmB = LF_qmB[ifc,:] @ LF_qmB[ifc,:].conj()
              npA = LF_qpA[ifc,:] @ LF_qpA[ifc,:].conj()
              npB = LF_qpB[ifc,:] @ LF_qpB[ifc,:].conj()
              normd = sqrt( nmA.real + nmB.real + npA.real + npB.real )
              if (abs(normd)<1e-6):
                  print('ERROR : normd[ifc='+str(ifc)+']=0, hence the choosen function is singular, normd=', normd)
                  sys.exit(1)
              LF_qmA[ifc,:] *= 1./normd
              LF_qmB[ifc,:] *= 1./normd
              LF_qpA[ifc,:] *= 1./normd
              LF_qpB[ifc,:] *= 1./normd
              print('norm[ifc=',ifc,']=', normc, normd)
          else:
              print('norm[ifc=',ifc,']=', normc)
          print('LF_rmA[%d]=%10.5f  LF_rmB[%d]=%10.5f' % (ifc,linalg.norm(LF_rmA[ifc,:]),ifc,linalg.norm(LF_rmB[ifc,:])) )
          print('LF_rpA[%d]=%10.5f  LF_rpB[%d]=%10.5f' % (ifc,linalg.norm(LF_rpA[ifc,:]),ifc,linalg.norm(LF_rpB[ifc,:])) )
          if BA_region_needed:
              print('LF_qmA[%d]=%10.5f  LF_qmB[%d]=%10.5f' % (ifc,linalg.norm(LF_qmA[ifc,:]),ifc,linalg.norm(LF_qmB[ifc,:])) )
              print('LF_qpA[%d]=%10.5f  LF_qpB[%d]=%10.5f' % (ifc,linalg.norm(LF_qpA[ifc,:]),ifc,linalg.norm(LF_qpB[ifc,:])) )

      Nband = len(Eks[0])

      Nrr = len(rmA)
      Niik = k_basis.shape[0]
      Nk = int(len(k_basis[0])/2)

      Check_normalization = False
      _times_ = zeros(4)
      print('Nband=', Nband, 'Nk=', Nk, 'Niik=', Niik, 'Nrr=', Nrr)
      #
      QmA = zeros((Nband,Nrr),dtype=complex) # Kohn-Sham wave function in real space projected on A sites and upper plane
      QmB = zeros((Nband,Nrr),dtype=complex) # Kohn-Sham wave function in real space projected on B sites and upper plane
      QpA = zeros((Nband,Nrr),dtype=complex) # Kohn-Sham wave function in real space projected on A sites and lower plane
      QpB = zeros((Nband,Nrr),dtype=complex) # Kohn-Sham wave function in real space projected on B sites and lower plane
      #
      UAA = zeros( (Niik,Nfunc,Nband), dtype=complex)  # Current projector for AA region, which still needs to be orthogonalized
      if BA_region_needed:
          UAB = zeros( (Niik,Nfunc,Nband), dtype=complex)  # Current projector for AB region, which still needs to be orthogonalized
      #
      OtAA = zeros( (Nband,Nband), dtype=complex)
      OtAB = zeros( (Nband,Nband), dtype=complex)
      for iik in range(Niik):
          _t1_ = timer()
          ks = k_basis[iik,:Nk,:] # ks[k_basis,2]
          ps = k_basis[iik,Nk:,:] # ps[p_basis,2]
          
          exp_k_BA = exp( ks@Rs_BA * 1j ) # exp_k_BA[k_basis] = exp( ks[k_basis,2]*Rs_BA[2] * i )
          exp_p_BA = exp( ps@Rs_BA * 1j ) # exp_p_BA[p_basis] = exp( ps[p_basis,2]*Rs_BA[2] * i )
      
          cn = 1./sqrt(Nrr)
          exp_kA = exp( ks@rmA.T *1j) * cn   # exp_kA[k_basis,ir] = exp( ks[k_basis,:]*rmA.T[:,ir] * 1j)
          exp_kB = exp( ks@rmB.T *1j) * cn   # exp_kB[k_basis,ir] = exp( ks[k_basis,:]*rmB.T[:,ir] * 1j)
          exp_pA = exp( ps@rpA.T *1j) * cn   # exp_pA[k_basis,ir] = exp( ps[p_basis,:]*rpA.T[:,ir] * 1j)
          exp_pB = exp( ps@rpB.T *1j) * cn   # exp_pB[k_basis,ir] = exp( ps[p_basis,:]*rpB.T[:,ir] * 1j)
          _t2_ = timer()
          # changing Kohn-Sham wave function to real space by Psi(k,n)*e^{i*k*r}
          Psi_mA = Psis[iik,:, 0:2*Nk:2] # Psi_ks projected on A site and upper layer
          QmA = Psi_mA @ exp_kA          # QmA(n_band,R^-_A) = \sum_k_A A_{iik, n_band, k_A} * exp^{i k_A * R^-_A}
          Psi_mB = Psis[iik,:, 1:2*Nk:2] # Psi_ks projected on B site and upper layer
          QmB = Psi_mB @ exp_kB          # QmB(n_band,R^-_B) = \sum_k_B A_{iik, n_band, k_B} * exp^{i k_B * R^-_B}
          Psi_pA = Psis[iik,:,2*Nk  ::2] # Psi_ks projected on A site and lower layer
          QpA = Psi_pA @ exp_pA          # QpA(n_band,R^+_A) = \sum_p_A A_{iik, n_band, p_A} * exp^{i p_A * R^+_A}
          Psi_pB = Psis[iik,:,2*Nk+1::2] # Psi_ks projected on B site and lower layer
          QpB = Psi_pB @ exp_pB          # QpB(n_band,R^+_B) = \sum_p_B A_{iik, n_band, p_B} * exp^{i p_B * R^+_B}
          if Check_normalization:
              # overlap. Just checking normalization, but not essential
              OtAA += QmA @ QmA.T.conj()
              OtAA += QmB @ QmB.T.conj()
              OtAA += QpA @ QpA.T.conj()
              OtAA += QpB @ QpB.T.conj()
          _t3_ = timer()
          # Computing Projectors: <phi_{LF}| psi_{n k}>, which is done by real space sum,
          # like: <phi_{LF}(r)|psi_{n}(r)>
          UmA = LF_rmA @ QmA.T  # UmA[ifc,nband] = LF_rmA[ifc,R^-_A] * QmA.T[R^-_A,n_band]
          UmB = LF_rmB @ QmB.T  # UmB[ifc,nband] = LF_rmB[ifc,R^-_B] * QmB.T[R^-_B,n_band]
          UpA = LF_rpA @ QpA.T  # UpA[ifc,nband] = LF_rpA[ifc,R^+_A] * QpA.T[R^+_A,n_band]
          UpB = LF_rpB @ QpB.T  # UpB[ifc,nband] = LF_rpB[ifc,R^+_B] * QpB.T[R^+_B,n_band]
          UAA[iik,:,:] = UmA + UmB + UpA + UpB # This is <local_function|Kohn-Sham function>
          _t4_ = timer()
          if BA_region_needed:
              # projection to BA region: exp_kA[k_basis,ir]*exp_k_BA[k_basis]
              exp_kA *= exp_k_BA[:, None]  # exp_kA[k_basis,ir] = exp_kA[k_basis,ir]*exp_k_BA[k_basis]
              exp_kB *= exp_k_BA[:,None]   # exp_kB[k_basis,ir] = exp_kB[k_basis,ir]*exp_k_BA[k_basis]
              exp_pA *= exp_p_BA[:,None]   # exp_pA[p_basis,ir] = exp_pA[p_basis,ir]*exp_p_BA[p_basis]
              exp_pB *= exp_p_BA[:,None]   # exp_pB[p_basis,ir] = exp_pB[p_basis,ir]*exp_p_BA[p_basis]
              QmA = Psi_mA @ exp_kA
              QmB = Psi_mB @ exp_kB 
              QpA = Psi_pA @ exp_pA
              QpB = Psi_pB @ exp_pB
              # Projectors: <phi_{LF}| psi_{n k}>
              VmA = LF_qmA @ QmA.T  # VmA[ifc,nband] = LF_qmA[ifc,R^-_A] * QmA.T[R^-_A,n_band]
              VmB = LF_qmB @ QmB.T  # VmB[ifc,nband] = LF_qmB[ifc,R^-_B] * QmB.T[R^-_B,n_band]
              VpA = LF_qpA @ QpA.T  # VpA[ifc,nband] = LF_qpA[ifc,R^+_A] * QpA.T[R^+_A,n_band]
              VpB = LF_qpB @ QpB.T  # VpB[ifc,nband] = LF_qpB[ifc,R^+_B] * QpB.T[R^+_B,n_band]
              UAB[iik,:,:] = VmA + VmB + VpA + VpB
              if Check_normalization:
                  # overlap. Just checking normalization, but not essential
                  OtAB += QmA @ QmA.T.conj()
                  OtAB += QmB @ QmB.T.conj()
                  OtAB += QpA @ QpA.T.conj()
                  OtAB += QpB @ QpB.T.conj()
          _t5_ = timer()
          _times_[0] += _t2_-_t1_
          _times_[1] += _t3_-_t2_
          _times_[2] += _t4_-_t3_
          _times_[3] += _t5_-_t4_
          Nband_2 = int(Nband/2)
          #print(iik, Eks[iik][Nband_2-2], Eks[iik][Nband_2-1], Eks[iik][Nband_2], Eks[iik][Nband_2+1])
          if (iik%100==0):
              print('%4d'%(iik), 't_exp=%8.5f t_P*e=%8.5f t_LF*Q=%8.5f t_BA=%8.5f' % tuple(_times_))
      print('%4d'%(iik), 't_exp=%8.5f t_P*e=%8.5f t_LF*Q=%8.5f t_BA=%8.5f' % tuple(_times_))
      
      if Check_normalization:
          print('Checking normalization of eigenvectors: OtAA=', ('%10.6f '*Nband)%tuple([real(OtAA[i,i]/Niik) for i in range(Nband)]))
          if BA_region_needed:
              print('Checking normalization of eigenvectors: OtAB=', ('%10.6f '*Nband)%tuple([real(OtAB[i,i]/Niik) for i in range(Nband)]))

      # Renormalizing the projector
      OlpAA = zeros((Nfunc,Nfunc),dtype=complex)
      Nfunc=builtins.min(8,N_loc_func)  # Now we rearange them and orhogonalize them
      UAA_new = zeros((Nfunc,Niik,Nband),dtype=complex)
      for iik in range(Niik):
          UAA_new[:,iik,:] = OrthogonalizeDiag2(UAA[iik,:Nfunc,:],'AA',False)
          
      if BA_region_needed:
          OlpAB = zeros((Nfunc,Nfunc),dtype=complex)
          UAB_new = zeros((Nfunc,Niik,Nband),dtype=complex)
          for iik in range(Niik):
              UAB_new[:,iik,:] = OrthogonalizeDiag2(UAB[iik,:Nfunc,:],'AB',False)

              
      savez('LocalFunc'+band, LocFunc_Km=LocFunc_Km, LocFunc_Kp=LocFunc_Kp, LocFunc_phasem=LocFunc_phasem, LocFunc_phasep=LocFunc_phasep) # SOlapAA=SOlapAA, SOlapAB=SOlapAB
      
      cc = zeros((Nfunc,Nfunc),dtype=complex)
      w = linspace(-Lpl,Lpl,Npl) # frequency on real axis
      pDOS = zeros((2,Nfunc,len(w)),dtype=float)
      DOS = zeros((3,len(w)),dtype=float)
      UU  = zeros((Nfunc,Nband),dtype=float)
      UU2 = zeros((Nfunc,Nfunc,Nband), dtype=complex)
      Gc  = zeros((2,Nfunc,Nfunc,len(w)), dtype=complex)
      for iik in range(Niik):
          cc += UAA_new[:,iik,:]@UAA_new[:,iik,:].T.conj()   # cc[ifc1,ifc2]=\sum_{iband} UAA_new[ifc1,iik,iband]*UAA_new.T[iband,iik,ifc2]^*
          # need reasonable broadening even when we have many k-points
          broad_values = npy.minimum(wbroad*npy.abs(w), 1.0/unit)
          om = w + broad_values*1j
          ek = Eks[iik,:]
          denom0 = 1.0/(om[npy.newaxis, :] - ek[:,npy.newaxis]) # denom0[iband,om] = 1/(om-ek)
          denom1 = -1./pi * denom0.imag  # denom1[iband,om] = -1/pi*im( 1/(om-ek) )
          DOS[0,:] += sum(denom1,axis=0)
          
          # UU2[ifb,ifc,ibnd] = UAA_new[ifb,ik,ibnd] * UAA_new[ifc,ik,ibnd]^*
          # Gc[ifb,ifc,iw] = sum_ibnd UAA_new[ifb,ik,ibnd] * denom[ibnd,iw] * UAA_new[ifc,ik,ibnd]^* = \sum_ibnd UU2[ifb,ifc,ibnd]*denom[ibnd,iw]
          #for ibnd in range(Nband):
          #    UU2[:,:,ibnd] = tensordot(UAA_new[:,iik,ibnd], UAA_new[:,iik,ibnd].conj(), axes=0)
          UU2 = npy.einsum('ib,jb->ijb', UAA_new[:,iik,:], UAA_new[:, iik, :].conj() )
          Gc[0,:,:,:] += UU2 @ denom0
          if BA_region_needed:
              #for ibnd in range(Nband):
              #    UU2[:,:,ibnd] = tensordot(UAB_new[:,iik,ibnd], UAB_new[:,iik,ibnd].conj(), axes=0)
              UU2 = npy.einsum('ib,jb->ijb', UAB_new[:,iik,:],UAB_new[:,iik,:].conj() )
              Gc[1,:,:,:] += dot(UU2, denom0)
          
          # UU[ifb,ifc,ibnd] = A[ifb,ibnd] * A[ifc,ibnd]
          # UAA_new[ifc,iik,iband]*UAA_new[ifc,iik,iband]^*
          # pDOS[ifc,iw] += UU[ifc,ibnd] * denom[ibnd,iw]
          UU[:,:] = (UAA_new[:,iik,:] * UAA_new[:,iik,:].conj()).real
          pDOS[0,:,:] += UU @ denom1
          #
          if BA_region_needed:
              UU[:,:] = (UAB_new[:,iik,:] * UAB_new[:,iik,:].conj()).real
              pDOS[1,:,:] += UU @ denom1

      Gc *= 1./Niik
      pDOS *= 1./Niik
      DOS  *= 1./Niik
      DOS[1,:]  = sum(pDOS[0,:,:],axis=0) # AA region
      DOS[2,:]  = sum(pDOS[1,:,:],axis=0) # AA region
      cc *= 1./Niik
      print('check normalization:')
      for ifc in range(Nfunc):
            print('  '+('%12.8f '*Nfunc) % tuple(cc[ifc,:].real)) 

      if BA_region_needed:
          savez('Uproject'+band, UAA=UAA_new, UAB=UAB_new)
      else:
          savez('Uproject'+band, UAA=UAA_new)
          
      fout = open('pDOS.dat','w')
      print('# total  AA+AB  AA   AB  ', end=' ', file=fout)
      for i in range(Nfunc) :
            print('funAA['+str(i)+']  ', end=' ', file=fout)
      for i in range(Nfunc) :
            print('funAB['+str(i)+']  ', end=' ', file=fout)
      print(file=fout)
      for iw in range(len(w)):
            print(w[iw], end=' ', file=fout)
            print(DOS[0,iw], DOS[1,iw]+DOS[2,iw], DOS[1,iw], DOS[2,iw], end=' ', file=fout)
            for i in range(len(pDOS)):
              for ifc in range(Nfunc):
                print(pDOS[i,ifc,iw], end=' ', file=fout)
            print(file=fout)
      fout.close()

      Nregions = 2 if BA_region_needed else 1
      #Delta+Eimp = w - Gc^{-1}

      with open('Gloc.dat', 'w') as fg, open('Delta.dat', 'w') as gd, open('Delta.inp', 'w') as fd:
          for iw in range(len(w)):
                print(w[iw], end=' ', file=fg)
                for iAA in range(Nregions):
                  for i in range(Nfunc):
                    for j in range(Nfunc):
                      print(Gc[iAA,i,j,iw].real, Gc[iAA,i,j,iw].imag, '  ', end=' ', file=fg)
                print(file=fg)
                print(w[iw], end=' ', file=fd)
                print(w[iw], end=' ', file=gd)
                for iAA in range(Nregions):
                  DltE = w[iw]*ones((Nfunc,Nfunc)) - linalg.inv(Gc[iAA,:,:,iw])
                  for i in range(Nfunc):
                    print(DltE[i,i].real, DltE[i,i].imag, '  ', end=' ', file=fd)
                    for j in range(Nfunc):
                      print(DltE[i,j].real, DltE[i,j].imag, '  ', end=' ', file=gd)
                print(file=fd)
                print(file=gd)
      
