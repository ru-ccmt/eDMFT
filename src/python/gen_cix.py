import sys, math
from functools import reduce
import operator
from scipy import *
import numpy as np
from numpy import *
from numba import njit
import time
from gclusters import *

# GENERATE CIX FILE FOR BHQMC


def ComputeHamiltonianMomentumSpace(ek,Ks,U,basis,Nm,masks,ind_basis,q1q2q3):
    """ Note: we just changed bit meaning so that order is [[up,dn],[up,dn],[up,dn],...]
    """
    #DEBUG=False
    Ham = zeros((len(basis),len(basis)),dtype=float)
    for ii,state in enumerate(basis):
        #if state in [3,12,48,192]: DEBUG=True
        #if state in [18,33,72,132]: DEBUG=True
        masks_state = state&masks # here we apply all masks on a state at once, and later use them
        for ik in range(Nm): # over all cluster momenta
            if masks_state[2*ik]: # we have down in this k
                Ham[ii,ii] += ek[ik]
            if masks_state[2*ik+1]: # we have up in this k
                Ham[ii,ii] += ek[ik]
        # We compute U*n_{up}n_{dn} = U/N * sum_{q1,q2,q3} psi^+_{q1+q2-q3,up} psi^+_{q3,dn} psi_{q2,dn} psi_{q1,up}
        #     U/N psi^+_{i4,up} psi_{i1,up} psi^+_{i3,dn} psi_{i2,dn}  -> 
        #     U/N psi^+_{2*i4} psi_{2*i1} psi^+_{2*i3+1} psi_{2*i2+1} 
        for i3 in range(Nm):  
            for i2 in range(Nm):
                if i2==i3:
                    isig = 1
                elif i2<i3:
                    isig = sign2(state,masks[2*i2+1],masks[2*i3+1])
                else:
                    isig = sign2(state,masks[2*i3+1],masks[2*i2+1])
                intermediate_state = state^masks[2*i2+1]^masks[2*i3+1]
                for i1 in range(Nm):
                    i4 = q1q2q3[i1,i2,i3] # Ks[i4] = q1+q2-q3
                    #if masks_state[i1+Nm] and masks_state[i2]:
                    if masks_state[2*i1] and masks_state[2*i2+1]:  # checking that i1_up, and i2_dn exist
                        #print('i1=', i1, 'i2=', i2, 'i3=', i3, 'i4=', i4)
                        # We require psi^+_{q1+q2-q3,up} psi^+_{q3,dn} psi_{q2,dn} psi_{q1,up}
                        if i1==i4 and i2==i3:
                            Ham[ii,ii] += U/Nm
                        #elif not(masks_state[i3]) and not(masks_state[i4+Nm]):
                        elif not(masks_state[2*i3+1]) and not(masks_state[2*i4]):
                            if i1==i4:
                                jsig=1
                            elif i1<i4:
                                jsig = sign2(intermediate_state,masks[2*i1],masks[2*i4])
                            else:
                                #if DEBUG:
                                #    print('intermediate_state=', intermediate_state, PrintSt(intermediate_state,masks),'masks=', masks[2*i4],masks[2*i1], 'state=',state)
                                #    print('because state=', PrintSt(state,masks),'i2=',i2,'i3=',i3,'mask[i2]=',masks[2*i2+1],'mask[i3]=',masks[2*i3+1])
                                #    print('here=', state^masks[2*i2+1]^masks[2*i3+1])
                                jsig = sign2(intermediate_state,masks[2*i4],masks[2*i1])
                            #
                            newstate = state^masks[2*i1]^masks[2*i2+1]^masks[2*i3+1]^masks[2*i4]
                            #newstate = intermediate_state^masks[2*i1]^masks[2*i4]
                            jj = ind_basis[newstate]
                            #if DEBUG:
                            #    print('checking:', PrintSt(state, masks)+'->'+PrintSt(newstate, masks), 'i1=', i1, 'i2=', i2, 'i3=', i3, 'i4=', i4,
                            #              'isig=', isig, 'jsig=', jsig)
                            Ham[ii,jj] += U/Nm*isig*jsig
    #if not np.allclose(Ham.T,Ham):
    #    print('ERROR Ham not Hermitian', Ham)
    return Ham



def PrintSt(state, masks):
    symbol = [['0','d'],['u','2']]
    Nm = int(len(masks)/2)
    state_masks = array((state&masks).astype(bool), dtype=int) # gives 0 or 1 for each bin that is occupied or empty
    # knowing which beans are occupied, we know that if bin j is occupied and j+N is empty, we have up spin at j
    # if bin j and bin j+N is occupied, we have doubly occupied state, which we call '2'
    st = ''.join([symbol[state_masks[2*j]][state_masks[2*j+1]] for j in range(Nm)])
    return st
    
def Precompute_N_Sz_momentum(Nm,Ni,Szi,momentum):
    """ Computes N and Sz for all states in direct basis. It computes it for all states in the basis at once.
    Should be very efficient.
    """
    Nstates = 1<<(2*Nm)
    states = arange(Nstates)
    Ns = zeros(Nstates,dtype=int)
    Sz = zeros(Nstates,dtype=int)
    Kstate = zeros((Nstates,2), dtype=float)
    for bit_index in range(2*Nm):
        w2 = ((states >> bit_index) & 1)    # we check if bit_index is set for all states in the basis at once
        Ns += w2*Ni[bit_index]              # w2 is either 0 or 1 and is array of the size of all basis. Ni is occupancy.
        Sz += w2*Szi[bit_index]             # Szi is spin value associated with bit bit_index
        Kstate += outer(w2,momentum[bit_index,:])  # Kstate[i,2] += w2[i]*momentum[bit_index,2]
    return (Ns,Sz,Kstate)

def FindEquivalentK(k, K_int, B_int):
    dK = k-K_int
    which = where(sum(abs(dK),axis=1)==0)[0]
    if len(which)>0:
        return which[0]
    for ix in [-1,0,1,2,-2]:
        for iy in [-1,0,1,2,-2]:
            if ix==0 and iy==0: continue
            _dK_ = dK + B_int[0]*ix + B_int[1]*iy
            which = where(sum(abs(_dK_),axis=1)==0)[0]
            if len(which)>0:
                return which[0]
    return -1

def Find_q1_p_q2_m_q3(Ks,Nm):
    """ Find index array for momentum q1+q2-q3
    If Ks[i1]=q1 and Ks[i2]=q2 and Ks[i3]=q3 than Ks[q1q2q3[i1,i2,i3]] is equivalent to q1+q2-q3
    """
    K_int = array(np.round(Ks*Nm),dtype=int)
    B_int = array([[1,0],[0,1]])*Nm
    #print('B_int=', B_int)
    #print('K_int=', K_int)
    q1q2q3 = zeros((len(Ks),len(Ks),len(Ks)),dtype=int)
    for i3,q3 in enumerate(K_int):
        for i2,q2 in enumerate(K_int):
            for i1,q1 in enumerate(K_int):
                k = q1+q2-q3
                dK = k-K_int
                i4 = FindEquivalentK(k, K_int, B_int)
                q1q2q3[i1,i2,i3] = i4
                #if sum(abs(k-K_int[i4]))!=0:
                #    print(i4, 'q1=', q1, 'q2=', q2, 'q3=', q3, 'k=', k, '=', K_int[i4])
    return q1q2q3

def Find_q1_m_q2(Ks,Nm):
    """ Find index array for momentum q1+q2
    If Ks[i1]=q1 and Ks[i2]=q2 than Ks[momentum_plus[i1,i2]] is equivalent to q1+q2
    """
    K_int = array(np.round(Ks*Nm),dtype=int)
    B_int = array([[1,0],[0,1]])*Nm
    momentum_plus = zeros((len(Ks),len(Ks)),dtype=int)
    for i2,q2 in enumerate(K_int):
        for i1,q1 in enumerate(K_int):
            k = q1+q2
            dK = k-K_int
            i3 = FindEquivalentK(k, K_int, B_int)
            momentum_plus[i1,i2] = i3
    return momentum_plus

@njit
def k_within_01(k):
    """ 
    Returns number between 0-1 from any floating point number.

    Momenta are stored in units of 2*pi. To bring them to equivalent canonical form (within 1BZ), we can 
    subtract or add any integer. Hence 1BZ is defined as number between 0-1.
    """
    for j in range(len(k)):
        if k[j]>=1 or k[j]<0:
            while (k[j]>=1): k[j]-=1
            while (k[j]< 0): k[j]+=1
    return k

def CombineEquiv(Kunique,Ks):
    """ Given a lits of k-points (Kunique) it finds groups of k-points which are equivalent.
    k-points are given in units of 2*pi, hence k+(1,0) or k+(0,1) is an equivalent k-point.
    Input:
       Kunique  -- list of 2D k-points
       Ks       -- list of cluster k-points, which are of course unique
    Output:
       grp      -- list or groups [[0,1,2],[3],[4,5,6]] means 0,1,2 are equivalent, 3 is its own group and 4,5,6 are equivalent
       kt       -- index to k-point in Ks. Namely, eack k-point in Kunique can be found also in cluster Ks momenta. 0,1,2 are all equivalent to an entry in Ks, and this gives number of that index.
    """
    Nm = len(Ks)
    K_int = array(np.round([k_within_01(k)*Nm for k in Ks]),dtype=int)
    Kequiv = copy(Kunique)
    for i,k in enumerate(Kequiv):
        Kequiv[i] = k_within_01(k)
    grp=[]
    kt=[]
    rest = range(0,len(Kequiv))
    while len(rest)>0:
        k = Kequiv[rest[0]]
        g = np.where(np.all(Kequiv[rest]==k, axis=1))[0]
        if len(g)>0:
            g2 = array(rest)[g]
            rest = [item for i,item in enumerate(rest) if i not in g]
            grp.append(g2.tolist())

            dK = array(np.round(k*Nm),dtype=int) - K_int
            which = where(sum(abs(dK),axis=1)==0)[0]
            if len(which)>0:
                ik = which[0]
            else:
                print('ERROR k=', k, 'which seems absent in cluster momenta Ks=', Ks)
                sys.exit(0)
            kt.append(ik)
    return (grp,kt)

def Get_momenta_Sigind(Ks):
    Nm = len(Ks)
    grp=[]
    rest = range(0,Nm)
    while len(rest)>0:
        k0 = Ks[rest[0]]
        k1 = k0[::-1]        # equivalent k with 90 degree rotation
        k2 = [-k0[0],k0[1]]  # x mirror
        dk0 = np.all(Ks[rest]==k0,axis=1)
        dk1 = np.all(Ks[rest]==k1,axis=1)
        dk2 = np.all(Ks[rest]==k2,axis=1)
        g = np.where( logical_or(logical_or( dk0, dk1), dk2) )[0]
        if len(g)>0:
            g2 = array(rest)[g]
            rest = [item for i,item in enumerate(rest) if i not in g]
            grp.append(g2.tolist())
    Sigind = zeros(Nm, dtype=int)
    for ig,g in enumerate(grp):
        for c in g:
            Sigind[c]=ig
    return Sigind.tolist()
            
def largest_linear_dim_from_Nm(Nm: int) -> int:
    if Nm <= 0:
        raise ValueError("Nm must be a positive integer")
    d = math.isqrt(Nm)                  # floor(sqrt(Nm))
    while Nm % d != 0:                  # find largest divisor ≤ sqrt(Nm)
        d -= 1
    return Nm // d                      # larger side of the closest-to-square pair

def k_to_ij(k, N1):
    """Map fractional k to integer grid indices (i,j) with wrapping."""
    k = np.asarray(k, dtype=float)
    kmod = np.mod(k, 1.0)           # [0,1)
    i = int(np.rint(kmod[0] * N1)) % N1
    j = int(np.rint(kmod[1] * N1)) % N1
    return (i, j)
def Get_minus_K(Ks,Nm):
    """Given list of momenta in 1BZ, returns index to -k for each entry
    """
    N1 = largest_linear_dim_from_Nm(Nm)
    ind_k={}
    for ik,k in enumerate(Ks):
        ck_xy = k_to_ij(k, N1)
        ind_k[ck_xy] = ik
    minus_k={}
    for ik,k in enumerate(Ks):
        ck_xy = k_to_ij(-k, N1)
        imk = ind_k[ck_xy]
        minus_k[ik] = imk
    return minus_k


def main2(U,a,Eimp=None,t0h=None,tp0=None,superc=0,ifDCA=0):
    tms=[[] for i in range(8)] # just for meassuring time it takes, contains string describing what is being timed
    tmi=zeros(len(tms))         # list of times 
    #tmx=zeros(len(tms)+1)       # list of times for one n, Sz_min

    if Eimp is None and t0h is None and tp0 is None:
        print('Error! at least one of Eimp, t0h, tp0 should be given!')
        print('Use t0h, tp0 for tight binding model.')
        print('Use Eimp to give Eimp directly, e.g. in eDMFT.')
        return 0
    tm0 = time.time()
    (Nm, Rxy, Ks, expK, shx, shy) = FindClusterQuantities(a)
    print('Number of sites: Nm=', Nm)
    print('Real space vectors in the cluster: Rxy=', Rxy.tolist())
    print('cluster momenta Ks=', Ks.tolist())
    for ik in range(1,Nm):
        dr = abs(sum(expK[ik,:]))
        if abs(dr)>1e-6:
            print('ERROR at K=', Ks[ik], 'sum_r(e^{ik*r})=', abs(sum(expK[ik,:])))
    for ir in range(1,Nm):
        dk = abs(sum(expK[:,ir]))
        if abs(dk)>1e-6:
            print('ERROR at R=', Rxy[ir], 'sum_k(e^{ik*r})=', abs(sum(expK[:,ir])))

    # which bit is used for what
    Ni = ones(2*Nm, dtype=int)                             # each bit stands for occupancy 1
    #Szi = hstack((-ones(Nm,dtype=int),ones(Nm,dtype=int))) # first Nm bits are for down, and next for up spin
    #print('Szi=', Szi)
    Szi = np.array([[1,-1] for i in range(Nm)]).flatten()
    #print('Szi=', Szi)
    #momentum = vstack((Ks,Ks))
    #print('K=', momentum)
    momentum = np.repeat(Ks, repeats=2, axis=0)
    #print('K=', momentum)
    masks = array([(1<<i) for i in range(2*Nm)],dtype=int) # masks that check bit i
    print('occupancy in bits: Ni=', Ni.tolist())
    print('Sz in bits: Szi=', Szi.tolist())
    print('momentum of bits=', momentum.tolist())
    print('masks=', masks.tolist())
    # creates 2^{2*Nm} direct states
    Nstates = 1<<(2*Nm)
    print('number of all direct states=', Nstates)

    q1q2q3 = Find_q1_p_q2_m_q3(Ks,Nm)
    momentum_plus = Find_q1_m_q2(Ks,Nm)
    tm1 = time.time()
    tmi[0] += tm1-tm0
    tms[0] = 't_cluster_prepare'
    Ns,Sz,Kstate = Precompute_N_Sz_momentum(Nm,Ni,Szi,momentum)
    tm2 = time.time()
    tmi[1] += tm2-tm1
    tms[1] = 't_precompute_N_Sz'
    if Eimp is None:
        if Nm==4:
            ek = [-t0h*(cos(k[0]*2*pi)+cos(k[1]*2*pi))-tp0*cos(k[0]*2*pi)*cos(k[1]*2*pi) for k in Ks]
            if ifDCA:# nearest neighbor only!
                ek[0]=-8.0/np.pi*t0h
                ek[3]=8.0/np.pi*t0h
        elif Nm==8:
            ek = [-2*t0h*(cos(k[0]*2*pi)+cos(k[1]*2*pi))-2*tp0*cos(k[0]*2*pi)*cos(k[1]*2*pi) for k in Ks]
        elif Nm==2:
            ek = [-2*t0h*cos(k[1]*2*pi) for k in Ks]
        else:
            ek = [-2*t0h*(cos(k[0]*2*pi)+cos(k[1]*2*pi))-4*tp0*cos(k[0]*2*pi)*cos(k[1]*2*pi) for k in Ks]
    else:
        ek=Eimp
    print('ek=', ek)
    
    tms[2] = 't_find_k_grps'
    tms[3] = 't_cmp_Ham'
    tms[4] = 't_diag_Ham'
    tms[5] = 't_cmp_S2'
    tms[6] = 't_diag_S2'
    tms[7] = 't_store'
    basises = {}
    Hks = {}
    Es = {}
    ind_basises = {}
    S2a = {}
    Sa = {}
    Nstates = len(Ns)
    for n in range(2*Nm+1): # over blocks of constant occupancy N
        Sz_min = n%2        # minimum Sz when we have N particles in units of 1/2.
        Sz_max = n
        #for iSz in range(Sz_min,Sz_max+1,2):
        for iSz in range(n,-n-1,-2):#range(-n,n+1,2):
            tm3s = time.time()
            # we find all states from precomputed N,Sz array, which have N and Sz_min
            basis0 = np.where(logical_and(Ns==n , Sz==iSz))[0]
            #print('n=', n, 'Sz=', 0.5*iSz)
            #print('basis0=', basis0, [PrintSt(st,masks) for st in basis0])
            Kall = Kstate[basis0]
            #print('Kall=', Kall.tolist())
            #Kunique = np.unique(Kall, axis=0)
            
            grp,kt = CombineEquiv(Kall,Ks)
            
            if len(grp)>0: print('grp(K)=', grp)
                
            tm3e = time.time()
            tmi[2] += tm3e-tm3s
            for ik,this_k in enumerate(grp):
                tm4s = time.time()
                iik = kt[ik]
                #this_k = np.where(np.all(Kall==k, axis=1))[0]
                kall = Kall[this_k]
                k = k_within_01(kall[0]).tolist()
                basis = basis0[this_k]
                ind_basis = get_index(basis) # if we have state, we can get index in the list of basis, i.e., ind_basis[state]
                print('---- n='+str(n), 'Sz=', 0.5*iSz, 'k=', k, 'Nbasis=', len(basis),'----')
                PrintStates(basis,masks,n,iSz,ks=kall)
                Ham = ComputeHamiltonianMomentumSpace(ek,Ks,U,basis,Nm, masks,ind_basis,q1q2q3) # Just Hamiltonian in direct basis
                tm4e = time.time()
                tmi[3] += tm4e-tm4s
                print('H=', Ham)
                
                #print('groupsh=', groupsh)
                Ei, v = linalg.eigh(Ham)
                TkE = v.T.conj()
                tme5 = time.time()
                tmi[4] += tme5-tm4e
                S2 = ComputeS2(basis, Nm, masks, ind_basis, Szi)
                S2k = TkE @ S2 @ TkE.T.conj()   # This is S2 in eigenbasis
                tme6 = time.time()
                tmi[5] += tme6-tme5
                
                groupsi = FindCoupledBlocs(S2k)
                (TkS, S2i) = DiagonalizeBlocks(S2k, TkE, groupsi)
                tme7 = time.time()
                tmi[6] += tme7-tme6
                #if all(abs(S2k.imag)<1e-6) and all(abs(np.round(S2k.real*4)-S2k*4)<1e-6):
                #    S2k = np.round(S2k.real*4)/4.
                Si = [0.5*(sqrt(4*S2i[i]+1)-1) for i in range(len(S2i))]
                #print('S2=', S2k)
                print('Si=', Si)
                print('E=', [round(x.real*4)/4 if (abs(x.imag)<1e-8 and abs(round(x.real*4)-(x*4))<1e-8) else x for x in Ei])
                Hks[(n,iSz,iik)] = Ham
                basises[(n,iSz,iik)] = basis
                Es[(n,iSz,iik)] = Ei
                S2a[(n,iSz,iik)] = S2
                Sa[(n,iSz,iik)] = Si
                ind_basises[(n,iSz,iik)] = ind_basis
                tme8 = time.time()
                tmi[7] += tme8-tme7

    print('# Timings:')
    tsort = sorted(range(len(tms)), key= lambda i: -tmi[i])
    for i in tsort:
        print(tms[i]+'='+str(tmi[i])+' ', end=' ')
    print()
    
    max_size=0
    n_all=0
    global_index={}
    for n in range(2*Nm+1): # over blocks of constant occupancy N
        #for iSz in range(-n,n+1,2):
        for iSz in range(n,-n-1,-2):
            for ik in range(Nm):
                if (n,iSz,ik) in basises:
                    basis = basises[(n,iSz,ik)]
                    global_index[(n,iSz,ik)]=list(range(n_all+1,n_all+len(basis)+1))
                    n_all += len(basis)
                    if len(basis)>max_size: max_size=len(basis)

    Sigind = Get_momenta_Sigind(Ks)
    msigind = max(Sigind)
    eks = zeros(msigind+1)
    for i,ei in enumerate(ek):
        eks[Sigind[i]] = ei
    Sigind = np.repeat(Sigind, repeats=2, axis=0)
    glb_flip = np.repeat(list(range(Nm)), repeats=2, axis=0)


    ############################################
    n_Sz_ik_list = []
    for n in range(2*Nm+1): # over blocks of constant occupancy N
        for iSz in range(n,-n-1,-2): #range(-n,n+1,2):
            for ik in range(Nm):
                if (n,iSz,ik) in basises:
                    n_Sz_ik_list.append((n,iSz,ik))
    
    
    for n_iSz_iik in n_Sz_ik_list:
        (n,iSz,ik) = n_iSz_iik
        basis = basises[n_iSz_iik]
        for i,state in enumerate(basis):
            ii = global_index[n_iSz_iik][i]
            print(ii, 'n_iSz_ik=', n_iSz_iik, 'state=', state)
    #############################################

    
    if superc==0:
        # fcix = open('bctqmc_U_'+str(U)+'.cix', 'w')
        fcix = open('bctqmc.cix', 'w')
        print('# Cix file for cluster DMFT with Bold-CTQMC', file=fcix)
        print('# cluster_size, number of states, number of psi^+ operators, number of baths, number of blocks, maximum_matrix_size ', file=fcix)
        print(Nm, n_all, 2*Nm, 2*Nm, len(basises), max_size, file=fcix)
        print('# c^_i c_i Sigind spin   moment', file=fcix)
        for i in range(2*Nm):
            print(f"{i}  {i}  {i}*    {Sigind[i]}    {Szi[i]/2:4}    {Ks[i//2].tolist()}", file=fcix)
        print('# cluster energies for non-equivalent baths, eps[k]', file=fcix)
        print(('{:f} '*len(eks)).format(*eks), file=fcix)
        print(f"# iblck n k  Sz  ", " ".join([f"F^+_{t:1}" for t in range(2*Nm)]), file=fcix)
        i_block=1

        for n_iSz_iik in n_Sz_ik_list:
            (n,iSz,ik) = n_iSz_iik
            basis = basises[n_iSz_iik]
            ind0 = global_index[(n,iSz,ik)]
            for i,ii in enumerate(ind0):
                state = basis[i]
                masks_state = state&masks
                psi_dagger=zeros(2*Nm, dtype=int)
                for diks in range(2*Nm):
                    #s = 2*int(diks/Nm)-1
                    s = 1-2*(diks % 2)
                    #dik = diks % Nm
                    dik = diks//2
                    if not(masks_state[diks]):
                        new_state = state^masks[diks]
                        sgn = sign3(state,masks[diks])
                        #sgn = sign4(state, masks[diks], Nm)
                        ind_basis = ind_basises[(n+1,iSz+s,momentum_plus[ik,dik])]
                        j = ind_basis[new_state]
                        jj = global_index[(n+1,iSz+s,momentum_plus[ik,dik])][j]
                        #if diks==1:
                        #    print('start=',ii,'end=',jj,'sign=',sgn, 'state=', state, '=',bin(state), 'mask=', masks[diks],'masks=', masks)
                        psi_dagger[diks]=sgn*jj
                        
                sst=PrintSt(state, masks)
                pd = " ".join([f"{t:5}" for t in psi_dagger])
                print(f"{ii:<3d} {i_block:3} {n} {ik} {iSz/2:4.1f}  {pd} # {sst}", file=fcix)
            i_block += 1
        
        print('# All blocks of local Hamiltonian & S^2 follow. First line: iblock, dim, indices_for_block [energy_list,S2_list]', file=fcix)
        i_block=1
        for n in range(2*Nm+1): # over blocks of constant occupancy N
            #for iSz in range(-n,n+1,2):
            for iSz in range(n,-n-1,-2):
                for ik in range(Nm):
                    if (n,iSz,ik) in basises:
                        basis = basises[(n,iSz,ik)]
                        ind0 = global_index[(n,iSz,ik)]
                        print(i_block, len(ind0), '    ', ' '.join([str(x) for x in ind0]), '  ',
                              Es[(n,iSz,ik)].tolist(), Sa[(n,iSz,ik)], '# iblock, dim, indices, Energies, Spin', file=fcix)
                        Ham = Hks[(n,iSz,ik)]
                        S2 = S2a[(n,iSz,ik)]
                        for i in range(len(Ham)):
                            print(' ', ('{:10.6f} '*len(Ham)).format(*Ham[i,:]), ' '*5, ('{:10.6f} '*len(S2)).format(*S2[i,:]), file=fcix)
                        i_block += 1
        
        print('# UCoulomb : (m1,s1) (m2,s2) (m3,s2) (m4,s1)  Uc[m1,m2,m3,m4]', file=fcix)
        # We compute \sum_i U*n_{i,up}n_{i,dn} = 1/2*U \sum_{i,s} n_{i,s} n_{i,sb}  , where sb = \bar{s} and means opposite
        # 1/2*U/N^2 \sum_{s,K1,K2,K3,K4} psi^+_{K1,s} psi^+_{K2,sb} \psi_{K3,sb} \psi_{K4,s} \sum_{i} e^{(-K1-K2+K3+K4)R_i} 
        # 1/2*U/N \sum_{K1,K2,K3} psi^+_{K1,s} psi^+_{K2,sb} \psi_{K3,sb} \psi_{K1+K2-K3,s}
        # The definition we use is 1/2\sum_{abcd} U_{abcd} \psi^+_a \psi^+_b \psi_c \psi_d
        #
        # hence U_{abcd} = U/N delta_{a=K1,s} \delta_{b=K2,sb} \delta_{c=K3,sb} \delta_{d=K1+K2-K3,s}
        #
        for i3 in range(Nm):
            for i2 in range(Nm):
                for i1 in range(Nm):
                    i4 = q1q2q3[i1,i2,i3] # Ks[i4] = q1+q2-q3
                    print('{:4d}{:4d}{:4d}{:4d}   {:f}'.format(2*i1,2*i2+1,2*i3+1,2*i4, U/Nm), file=fcix)
                    print('{:4d}{:4d}{:4d}{:4d}   {:f}'.format(2*i1+1,2*i2,2*i3,2*i4+1, U/Nm), file=fcix) 

    else:
        ###############################################################################
        print('********** Superconductivity on ********************')
        WantToSeeSC = False
        # just index for -k from list Ks
        minus_k = Get_minus_K(Ks, Nm)
        # SC amplitude for d-wave SC on all momenta
        Nambu_d = [ round(sqrt(Nm)/4.*(np.cos(k[0]*2*pi)-np.cos(k[1]*2*pi))) for k in Ks]
        # pair of indices for ([k,up],[-k,dn]) where d-wave is nonzero
        SC_amplitude={}
        for ik in range(len(Nambu_d)):
            if (Nambu_d[ik]!=0):
                k = Ks[ik]
                mk = Ks[minus_k[ik]]
                ks_indices = (2*ik, 2*minus_k[ik]+1)
                print('Nambu_d[ik='+str(ik)+']=', Nambu_d[ik], '(k,up)=', ks_indices[0], '(-k,dn)=', ks_indices[1], 'k=', k, '-k=', mk)
                SC_amplitude[ks_indices] = Nambu_d[ik]
        print('SC_amplitude=', SC_amplitude)
        
        # Now we deal with normal state groups, which we did not optimize before
        # current groups for normal state grps
        grps=[[] for ig in range(len(basises))]
        for ig,n_iSz_iik in enumerate(basises):
            grps[ig] = n_iSz_iik
        #print('grps=', grps)
        # index to get index of the group from state, i.e., ind_block[state]=ig, for which grps[ig]=n_iSz_iik
        ind_block = zeros(Nstates, dtype=int)
        for ig,n_iSz_iik in enumerate(basises):
            basis = basises[n_iSz_iik]
            for state in basis:
                ind_block[state] = ig

        # This is psi^+ operator specified between blocks
        # Fbc[diks, iig] tells you: starting from block iig = (n,iSz,ik), applying a creation operator with momentum index dik and spin ±1 takes you to block ii or −1 if that state doesn’t exist.
        index_of_qnums = { n_iSz_iik: ig for ig, n_iSz_iik in enumerate(grps) }
        Fbc = zeros((2*Nm,len(grps)),dtype=int)
        for iig in range(len(grps)):
            (n,iSz,ik) = grps[iig]
            for diks in range(2*Nm):
                ds = 1-2*(diks % 2) # 1 or -1
                dik = diks//2      # momentum
                k_p_dik = momentum_plus[ik,dik]
                Psi_dagger_n_iSz_iik = (n+1,iSz+ds,k_p_dik)
                #ii=-1
                #if Psi_dagger_n_iSz_iik in grps:
                #    ii = grps.index(Psi_dagger_n_iSz_iik)
                ii = index_of_qnums.get(Psi_dagger_n_iSz_iik, -1)
                Fbc[diks,iig] = ii
                
        connected_states={}
        group = list(range(len(grps))) # group[ig] will contain groups that need to be combined due to superconductivity
        #    group[ig] is effectively a label assigning each block ig to some SC “super-group” (blocks that need to be merged).
        #    in normal state, we would have groups=[0,1,2,3,4,5,....] and in superconductivity groups=[0,1,0,1,0,2...]
        for iig,n_iSz_iik in enumerate(basises):
            # This loop ensures that blocks whose states are connected by pair creation/annihilation (SC) share the same group label.
            to_group=[]
            basis = basises[n_iSz_iik]
            for istate,state in enumerate(basis):
                masks_state = masks&state # here we apply all masks on a state at once, and later use them
                masks_state2 = masks_state.astype(bool).astype(int) # masks_state2[i] is 0 or 1 if masks[i] is set in state.
                #print('(n,iSz,iik)', (n,iSz,iik), 'state=', state, 'masks_state=', masks_state2)
                con_states=[]
                for (ks,mkms) in SC_amplitude:
                    if (masks_state2[ks]==0 and masks_state2[mkms]==0):
                        new_state = state^masks[ks]^masks[mkms]
                        con_states.append(new_state)
                        ig_new = ind_block[new_state]
                        if ig_new not in to_group:
                            to_group.append(ig_new)
                if con_states:
                    connected_states[state] = con_states
            for ig in to_group:
                group[ig] = group[iig]
            if to_group:
                print(iig, '(n,iSz,iik)=', n_iSz_iik, 'groups with ', to_group, '=[(n,iSz,iik)]->', [grps[ig] for ig in to_group])

        for itt in range(2):
            # Now we do two refinement passes (itt=0,1) to enforce: For each final SC block, applying F⁺ with any bath/momentum/spin index must map one into a unique SC block.
            #
            # We have combined states from grps which couple due to superconductivity in group.
            # But the issue we face is that F^+ operator on each constructed block might not give unique block.
            # To ensure that F^+ operator acting on a block gives unique block, we make two passes through this loop itt=0,1
            # and combine those block that F^+ operator mixes.
            
            for iig in range(len(group)): # this is just printing of current groups
                (n,iSz,ik) = grps[iig]
                print('group['+str(iig)+']=', group[iig], 'n_iSz_iik=', grps[iig], 'Fbc=', Fbc[:,iig])

            # Now we create lists dictionary of lists, where values in groups contains the set of groups that need to be combined
            groups={}
            for iig,ig in enumerate(group):
                if ig in groups:
                    groups[ig].append(iig)
                else:
                    groups[ig] = [iig]
            max_size2=0
            print('Found the following groups to combine in SC state')
            for grp in groups:
                dims = [len(basises[grps[g]]) for g in groups[grp]]
                tot_dim = sum(dims)
                if tot_dim>max_size2: max_size2=tot_dim
                print('grp['+str(grp)+']=', groups[grp], '  dimensions=', dims, 'total_size=', tot_dim)
                
            # And now we pass through these groups and check whether F^+ acting on a group gives unique group.
            # If not, we need to combine more groups together, so that this property is satisfied.
            for grp in groups:
                for ib in range(2*Nm):
                    fg = [group[Fbc[ib,g]] if (Fbc[ib,g] >= 0) else -1  for g in groups[grp]]
                    fg_vals = sorted(list({x for x in fg if x >= 0}))
                    print('Fbc['+str(ib)+',:]=', Fbc[ib,groups[grp]], 'which corresponds to group=', fg)
                    if len(fg_vals) >= 2:
                        # F⁺ from states in this block send you to multiple different SC blocks → we merge those SC blocks.
                        smallest = fg_vals[0]
                        for iig in range(len(group)):
                            if group[iig] in fg_vals[1:]:
                                #print('Found group['+str(iig)+']=', group[iig], 'setting it to', smallest)
                                group[iig] = smallest
                        if itt==0:
                            print('WARNING itt=', itt, 'multiple values in the group', fg_vals, 'need to set them to', smallest)
                        else:
                            print('ERROR   itt=', itt, 'multiple values in the group', fg_vals, 'need to set them to', smallest)
                            
            # The double pass of this loop is because merging groups in one round might create new conflicts that only show up after updating group;
            # the second pass ensures closure. On the second pass we print errors instead of warnings if something still isn’t closed.
            
        
        print('connected_states=', connected_states)

        #Sigind = Get_momenta_Sigind(Ks)
        #msigind = max(Sigind)
        #eks = zeros(msigind+1)
        #for i,ei in enumerate(ek):
        #    eks[Sigind[i]] = ei
        #Sigind = np.repeat(Sigind, repeats=2, axis=0)
        
        iblock = zeros(Nstates+1, dtype=int)
        for i_block,grp in enumerate(groups):
            for g in groups[grp]:
                n_iSz_iik = grps[g]
                (n,iSz,ik) = n_iSz_iik
                basis = basises[n_iSz_iik]
                for i,state in enumerate(basis):
                    ii = global_index[n_iSz_iik][i]
                    iblock[ii] = i_block+1

        if superc==1:
            fcix = open('bctqmc_U_'+str(U)+'_SC1.cix', 'w')
            print('# Cix file for cluster DMFT with Bold-CTQMC', file=fcix)
            print('# cluster_size, number of states, number of psi^+ operators, number of baths, number of blocks, maximum_matrix_size ', file=fcix)
            print(Nm, n_all, 2*Nm, 2*Nm+2*len(SC_amplitude), len(groups), max_size2, file=fcix)
            print('# c^+_i,c_i Sigind  spin   moment', file=fcix)
            for i in range(2*Nm):
                print(f"{i:2}  {i}  {i}*     {Sigind[i]}    {Szi[i]/2:4}    {Ks[i//2].tolist()}", file=fcix)
                
            msigind = max(Sigind)
            for i,kk in enumerate(SC_amplitude):
                ic = 2*Nm+2*i
                isigind = (msigind+abs(SC_amplitude[kk]))*np.sign(SC_amplitude[kk])
                print(f"{ic:2}  {kk[0]}  {kk[1]}     {isigind:2}      {0}     [0,0]", file=fcix)
                print(f"{ic+1:2}  {kk[1]}* {kk[0]}*    {isigind:2}      {0}     [0,0]", file=fcix)
            
            sc_amp = set([abs(SC_amplitude[sc]) for sc in SC_amplitude])
            print('# cluster energies for non-equivalent baths, eps[k]', file=fcix)
            print(('{:f} '*len(eks)).format(*eks), '0 '*len(sc_amp), file=fcix)
            print(f"# iblck n k  Sz  ", " ".join([f"F^+_{t:1}" for t in range(2*Nm)]), file=fcix)
            
            for n_iSz_iik in n_Sz_ik_list:
                (n,iSz,ik) = n_iSz_iik
                basis = basises[n_iSz_iik]
                for i,state in enumerate(basis):
                    ii = global_index[n_iSz_iik][i]
                    i_block = iblock[ii]
                    masks_state = state&masks
                    psi_dagger=zeros(2*Nm, dtype=int)
                    for diks in range(2*Nm):
                        s = 1-2*(diks % 2) # 1 or -1
                        dik = diks//2      # momentum
                        if not(masks_state[diks]):
                            new_state = state^masks[diks]
                            sgn = sign3(state,masks[diks])
                            k_p_dik = momentum_plus[ik,dik]
                            ind_basis = ind_basises[(n+1,iSz+s,k_p_dik)]
                            j = ind_basis[new_state]
                            jj = global_index[(n+1,iSz+s,k_p_dik)][j]
                            psi_dagger[diks]=sgn*jj
                            
                    sst=PrintSt(state, masks)
                    pd = " ".join([f"{t:5}" for t in psi_dagger])
                    print(f"{ii:<3d} {i_block:3} {n} {ik} {iSz/2:4.1f}  {pd} # {sst}", file=fcix)
            
            print('# All blocks of local Hamiltonian & S^2 follow. First line: iblock, dim, indices_for_block [energy_list,S2_list]', file=fcix)
            for i_block,grp in enumerate(groups):
                dims = [len(basises[grps[g]]) for g in groups[grp]]
                tot_dim = sum(dims)
                #print('dims=', dims, 'tot_dim=', tot_dim)
                Ham = zeros((tot_dim,tot_dim))
                S2 = zeros((tot_dim,tot_dim))
                Es_tot = zeros(tot_dim)
                S2_tot = zeros(tot_dim)
                states_ind=[]
                for ig,g in enumerate(groups[grp]):
                    n_iSz_iik = grps[g]
                    istart = sum(dims[:ig],dtype=int)
                    iend = istart+dims[ig]
                    #print('istart=', istart, 'iend=', iend)
                    Ham[istart:iend,istart:iend] = Hks[n_iSz_iik]
                    S2[istart:iend,istart:iend] = S2a[n_iSz_iik]
                    Es_tot[istart:iend] = Es[n_iSz_iik]
                    S2_tot[istart:iend] = Sa[n_iSz_iik]
                    states_ind += [global_index[n_iSz_iik][i] for i in range(len(basises[n_iSz_iik]))]
                    
                if WantToSeeSC:
                    print('states_ind=', states_ind)
                    for ig,g in enumerate(groups[grp]):
                        n_iSz_iik = grps[g]
                        istart = sum(dims[:ig],dtype=int)
                        for ist,state in enumerate(basises[n_iSz_iik]):
                            ii = global_index[n_iSz_iik][ist]
                            if state in connected_states:
                                i_ = states_ind.index(ii)
                                con_sts = connected_states[state]
                                for i in range(len(con_sts)):
                                    con_state = con_sts[i]
                                    _n_Sz_k_ = grps[ind_block[con_state]]
                                    j = ind_basises[_n_Sz_k_][con_state]
                                    jj = global_index[_n_Sz_k_][j]
                                    j_ = states_ind.index(jj)
                                    Ham[i_,j_] += 100.
                                    Ham[j_,i_] += 100.
                                    #print('   state=', ii, '-('+str(i_)+') connected to', con_state, '_n_Sz_k_=', _n_Sz_k_, 'j=', j,'j_=', j_, 'jj=', jj)
                                    
                print(i_block+1, tot_dim, '  ', *states_ind, '  ',
                          Es_tot.tolist(), S2_tot.tolist(), '# iblock, dim, indices, Energies, Spin', file=fcix)
                for i in range(len(Ham)):
                    print(' ', ('{:8.4f} '*len(Ham)).format(*Ham[i,:]), ' '*5, ('{:8.4f} '*len(S2)).format(*S2[i,:]), file=fcix)
            
            print('# UCoulomb : (m1,s1) (m2,s2) (m3,s2) (m4,s1)  Uc[m1,m2,m3,m4]', file=fcix)
            # We compute \sum_i U*n_{i,up}n_{i,dn} = 1/2*U \sum_{i,s} n_{i,s} n_{i,sb}  , where sb = \bar{s} and means opposite
            # 1/2*U/N^2 \sum_{s,K1,K2,K3,K4} psi^+_{K1,s} psi^+_{K2,sb} \psi_{K3,sb} \psi_{K4,s} \sum_{i} e^{(-K1-K2+K3+K4)R_i} 
            # 1/2*U/N \sum_{K1,K2,K3} psi^+_{K1,s} psi^+_{K2,sb} \psi_{K3,sb} \psi_{K1+K2-K3,s}
            # In SC state we redefine psi^+_K == c^+_{K,1} and psi_{-K} = c^+_{K,2}
            # Hence
            # 1/2*U/N \sum_{K1,K2,K3} ( c^+_{K1,1} c_{K2,2} c^+_{K3,2} c_{K1+K2-K3,1} + c_{K1,2} c^+_{K2,1} c_{K3,1} c^+_{K1+K2-K3,2} )
            # 1/2*U/N \sum_{K1,K2,K3} (-c^+_{K1,1} c^+_{K3,2} c_{K2,2} c_{K1+K2-K3,1} - c^+_{K2,1} c^+_{K1+K2-K3,2} c_{K1,2} c_{K3,1} )+1/2*U/N \sum_{K1,K2}(c^+_{K1,1} c_{K1,1} + c^+_{K2,1}c_{K2,1})
            # -1/2*U/N \sum_{K1,K2,K3} (c^+_{K1,1} c^+_{K2,2} c_{K3,2} c_{K1+K2-K3,1} + c^+_{K1,1} c^+_{K2,2} c_{K3,2} c_{K1+K2-K3,1} )+ U n_{up,tot}
            # The definition we use is 1/2\sum_{abcd} U_{abcd} \psi^+_a \psi^+_b \psi_c \psi_d
            #
            # hence U_{abcd} = -U/N delta_{a=K1,s} \delta_{b=K2,sb} \delta_{c=K3,sb} \delta_{d=K1+K2-K3,s}
            # and we have extra potential shift U n_{up,tot}
            for i3 in range(Nm):
                for i2 in range(Nm):
                    for i1 in range(Nm):
                        i4 = q1q2q3[i1,i2,i3] # Ks[i4] = q1+q2-q3
                        print('{:4d}{:4d}{:4d}{:4d}   {:f}'.format(2*i1,2*i2+1,2*i3+1,2*i4, U/Nm), file=fcix)
                        print('{:4d}{:4d}{:4d}{:4d}   {:f}'.format(2*i1+1,2*i2,2*i3,2*i4+1, U/Nm), file=fcix)
                    
        elif superc==2:
            fcix = open('bctqmc_U_'+str(U)+'_SC2.cix', 'w')
            print('# Cix file for cluster DMFT with Bold-CTQMC', file=fcix)
            print('# cluster_size, number of states, number of psi^+ operators, number of baths, number of blocks, maximum_matrix_size ', file=fcix)
            print(Nm, n_all, 2*Nm, 2*Nm+2*len(SC_amplitude), len(groups), max_size2, file=fcix)
            print('# c^_i c_i Sigind spin   moment', file=fcix)
            for i in range(2*Nm):
                if i%2==0:
                    print(f"{i:2}  {i}  {i}*     {Sigind[i]}    {Szi[i]/2:4}    {Ks[i//2].tolist()}", file=fcix)
                else:
                    print(f"{i:2}  {i}  {i}*    -{Sigind[i]}*   {Szi[i]/2:4}    {Ks[i//2].tolist()}", file=fcix)
                
            msigind = max(Sigind)
            for i,kk in enumerate(SC_amplitude):
                ic = 2*Nm+2*i
                isigind = (msigind+abs(SC_amplitude[kk]))*np.sign(SC_amplitude[kk])
                print(f"{ic:2}  {kk[0]}  {kk[1]}*    {isigind:2}      {0}     [0,0]", file=fcix)
                print(f"{ic+1:2}  {kk[1]}  {kk[0]}*    {isigind:2}      {0}     [0,0]", file=fcix)
            
            sc_amp = set([abs(SC_amplitude[sc]) for sc in SC_amplitude])
            print('# cluster energies for non-equivalent baths, eps[k]', file=fcix)
            print(('{:f} '*len(eks)).format(*eks), '0 '*len(sc_amp), file=fcix)
            print(f"# iblck n k  Sz  ", " ".join([f"F^+_{t:1}" for t in range(2*Nm)]), file=fcix)
            
            for n_iSz_iik in n_Sz_ik_list:
                (n,iSz,ik) = n_iSz_iik
                basis = basises[n_iSz_iik]
                for i,state in enumerate(basis):
                    ii = global_index[n_iSz_iik][i]
                    i_block = iblock[ii]
                    masks_state = state&masks
                    psi_dagger=zeros(2*Nm, dtype=int)
                    for diks in range(2*Nm):
                        s = 1-2*(diks % 2) # 1 or -1
                        dik = diks//2      # momentum
                        if (s==1 and not(masks_state[diks])) or (s==-1 and masks_state[diks]):
                            new_state = state^masks[diks]
                            sgn = sign3(state,masks[diks])
                            if s==1:
                                k_pm_dik = momentum_plus[ik,dik]
                                n_pm_one = n+1
                                Sz_pm_s = iSz+s
                            elif s==-1:
                                k_pm_dik = momentum_plus[ik,minus_k[dik]]
                                n_pm_one = n-1
                                Sz_pm_s = iSz-s
                            ind_basis = ind_basises[(n_pm_one,Sz_pm_s,k_pm_dik)]
                            j = ind_basis[new_state]
                            jj = global_index[(n_pm_one,Sz_pm_s,k_pm_dik)][j]
                            psi_dagger[diks]=sgn*jj
                            
                    sst=PrintSt(state, masks)
                    pd = " ".join([f"{t:5}" for t in psi_dagger])
                    print(f"{ii:<3d} {i_block:3} {n} {ik} {iSz/2:4.1f}  {pd} # {sst}", file=fcix)
            
            print('# All blocks of local Hamiltonian & S^2 follow. First line: iblock, dim, indices_for_block [energy_list,S2_list]', file=fcix)
            for i_block,grp in enumerate(groups):
                dims = [len(basises[grps[g]]) for g in groups[grp]]
                tot_dim = sum(dims)
                #print('dims=', dims, 'tot_dim=', tot_dim)
                Ham = zeros((tot_dim,tot_dim))
                S2 = zeros((tot_dim,tot_dim))
                Es_tot = zeros(tot_dim)
                S2_tot = zeros(tot_dim)
                states_ind=[]
                for ig,g in enumerate(groups[grp]):
                    n_iSz_iik = grps[g]
                    istart = sum(dims[:ig],dtype=int)
                    iend = istart+dims[ig]
                    #print('istart=', istart, 'iend=', iend)
                    Ham[istart:iend,istart:iend] = Hks[n_iSz_iik]
                    S2[istart:iend,istart:iend] = S2a[n_iSz_iik]
                    Es_tot[istart:iend] = Es[n_iSz_iik]
                    S2_tot[istart:iend] = Sa[n_iSz_iik]
                    states_ind += [global_index[n_iSz_iik][i] for i in range(len(basises[n_iSz_iik]))]
                    
                if WantToSeeSC:
                    print('states_ind=', states_ind)
                    for ig,g in enumerate(groups[grp]):
                        n_iSz_iik = grps[g]
                        istart = sum(dims[:ig],dtype=int)
                        for ist,state in enumerate(basises[n_iSz_iik]):
                            ii = global_index[n_iSz_iik][ist]
                            if state in connected_states:
                                i_ = states_ind.index(ii)
                                con_sts = connected_states[state]
                                for i in range(len(con_sts)):
                                    con_state = con_sts[i]
                                    _n_Sz_k_ = grps[ind_block[con_state]]
                                    j = ind_basises[_n_Sz_k_][con_state]
                                    jj = global_index[_n_Sz_k_][j]
                                    j_ = states_ind.index(jj)
                                    Ham[i_,j_] += 100.
                                    Ham[j_,i_] += 100.
                                    #print('   state=', ii, '-('+str(i_)+') connected to', con_state, '_n_Sz_k_=', _n_Sz_k_, 'j=', j,'j_=', j_, 'jj=', jj)
                                    
                print(i_block+1, tot_dim, '  ', *states_ind, '  ',
                          Es_tot.tolist(), S2_tot.tolist(), '# iblock, dim, indices, Energies, Spin', file=fcix)
                for i in range(len(Ham)):
                    print(' ', ('{:8.4f} '*len(Ham)).format(*Ham[i,:]), ' '*5, ('{:8.4f} '*len(S2)).format(*S2[i,:]), file=fcix)
            
            print('# UCoulomb : (m1,s1) (m2,s2) (m3,s2) (m4,s1)  Uc[m1,m2,m3,m4]', file=fcix)
            # We compute \sum_i U*n_{i,up}n_{i,dn} = 1/2*U \sum_{i,s} n_{i,s} n_{i,sb}  , where sb = \bar{s} and means opposite
            # 1/2*U/N^2 \sum_{s,K1,K2,K3,K4} psi^+_{K1,s} psi^+_{K2,sb} \psi_{K3,sb} \psi_{K4,s} \sum_{i} e^{(-K1-K2+K3+K4)R_i} 
            # 1/2*U/N \sum_{K1,K2,K3} psi^+_{K1,s} psi^+_{K2,sb} \psi_{K3,sb} \psi_{K1+K2-K3,s}
            # In SC state we redefine psi^+_K == c^+_{K,1} and psi_{-K} = c^+_{K,2}
            # Hence
            # 1/2*U/N \sum_{K1,K2,K3} ( c^+_{K1,1} c_{K2,2} c^+_{K3,2} c_{K1+K2-K3,1} + c_{K1,2} c^+_{K2,1} c_{K3,1} c^+_{K1+K2-K3,2} )
            # 1/2*U/N \sum_{K1,K2,K3} (-c^+_{K1,1} c^+_{K3,2} c_{K2,2} c_{K1+K2-K3,1} - c^+_{K2,1} c^+_{K1+K2-K3,2} c_{K1,2} c_{K3,1} )+1/2*U/N \sum_{K1,K2}(c^+_{K1,1} c_{K1,1} + c^+_{K2,1}c_{K2,1})
            # -1/2*U/N \sum_{K1,K2,K3} (c^+_{K1,1} c^+_{K2,2} c_{K3,2} c_{K1+K2-K3,1} + c^+_{K1,1} c^+_{K2,2} c_{K3,2} c_{K1+K2-K3,1} )+ U n_{up,tot}
            # The definition we use is 1/2\sum_{abcd} U_{abcd} \psi^+_a \psi^+_b \psi_c \psi_d
            #
            # hence U_{abcd} = -U/N delta_{a=K1,s} \delta_{b=K2,sb} \delta_{c=K3,sb} \delta_{d=K1+K2-K3,s}
            # and we have extra potential shift U n_{up,tot}
            for i3 in range(Nm):
                for i2 in range(Nm):
                    for i1 in range(Nm):
                        i4 = q1q2q3[i1,i2,i3] # Ks[i4] = q1+q2-q3
                        print('{:4d}{:4d}{:4d}{:4d}   {:f}'.format(2*i1,2*i2+1,2*i3+1,2*i4, -U/Nm), file=fcix)
                        print('{:4d}{:4d}{:4d}{:4d}   {:f}'.format(2*i1+1,2*i2,2*i3,2*i4+1, -U/Nm), file=fcix)
                    
            
                    
    print('N_all=', n_all)
    #print('You should also print S^2, so that one can block diagonalize G for all frequencies')
if __name__=='__main__':
    tp=-0.
    U = 8
    superc=0
    ifDCA=0# if the cix file is for DCA. currently it only works for a = array([[2,0],[0,2]])with NN hopping without SuperC.
    t0h = 1
    tp0 = tp # it seems for cellular DMFT and 4-sites. For 8 sites CDMFT we don't have translational invariance!

    
    #a=array([[1,0],[0,2]])     # N=2 site cluster
    a = array([[2,0],[0,2]])   # N=4 site cluster
    #a = array([[2,2],[-2,2]])   # N=8 site cluster
    #a = array([[4,0],[0,4]])   # N=16 site cluster

    main2(U,a,None,t0h,tp0,superc,ifDCA)

