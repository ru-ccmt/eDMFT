import numpy as np
from numpy import *
from numba import njit

def FindEquivalent1(shR,Nm,a,Rxy):
    """ Given a vector shR (which is for example shifted vector in the cluster) we need to find 
        which vector in Rxy is equivalent to this vector.
        
        Nm  --  size of the cluster.
        a   --  unit vectors of the cluster.
        Rxy -- all unique points in the cluster
    """
    dR = Rxy-shR
    which = where(sum(abs(dR),axis=1)==0)[0] # check if dR[j]==shR for any j and return j for which this is true.
    if len(which)>0:
        return which[0]
    # Even if shR is not it Rxy, its equivalent vector might be. So, go over all equivalent vectors.
    for ix in [-1,0,1]:
        for iy in [-1,0,1]:
            if ix==0 and iy==0: continue # this was done above, so no need to look at it again
            dR = Rxy-(shR+ix*a[0]+iy*a[1])
            which = where(sum(abs(dR),axis=1)==0)[0]
            if len(which)>0:
                return which[0]
    #for j in range(Nm):
    #    finite = True
    #    for ix in range(-1,2):
    #        for iy in range(-1,2):
    #            ci = shR + ix*a[0]+iy*a[1] # we add to shR any direct vectors, and see if is equivalent to one of the existing vectors
    #            finite = finite and any(Rxy[j]-ci)
    #    if not finite:
    #        return j
    return -1

def FindEquivalent2(cand_i, a, Rxy):
    """ Given a vector cand_i and all Rxy vectors inside the cluster, it finds if cand_i is unique 
        vector, or is just equivalent to another vector in Rxy.
    """
    finite=True
    # if cand[i]+n0*a[0]+n1*a[1] is equal to existing Rxy, than they are equivalent
    for ix in range(-1,2):
        for iy in range(-1,2):
            ci = cand_i+ix*a[0]+iy*a[1] # we add to cand[i] any direct vectors, and see if is equivalent to one of the existing vectors
            for r in Rxy:
                finite = finite and any(r-ci) # if ci is in Rxy, any gives 0, and finite is false. Than we return false.
                if not finite:
                    return finite
    return finite

def FindClusterQuantities(a):
    """Input is the cluster direct vectors a[0] and a[1] in integer form.
    The output is:
                  Nm      -- number of sites in the cluster
                  Rxy     -- (integer) vectors to all the sites inside the cluster
                  Ks      -- cluster momenta
                  expK    -- exponents for Fourier transform
                  shx     -- if we translate the system in x direction for one lattice spacing, we get different
                             vector from Rxy. Which?
    """
    positiveRxy=True
    Nm = int(abs(cross(a[0],a[1])))  # number of points in the cluster comes from volume, i.e., cross product
    bs = linalg.inv(a).T             # This is the reciprocal vector
    Nxy = max(max(abs(a[0,:])),max(abs(a[1,:]))) # we will look for all vectors inside the cluster by going from -Nxy to Nxy
    # First find all vectors which are shorter than a[0] and a[1]. They are inside the cluster. But many are equivalent.
    L2 = max(linalg.norm(a[0])**2,linalg.norm(a[1])**2) # this is length^2 of unit vectors.
    cand=[]
    #print('Nxy=',Nxy)
    for i in range(Nxy,-Nxy-1,-1):
        for j in range(Nxy,-Nxy-1,-1):
            if (i**2+j**2 < L2-1e-10):
                cand.append([i,j]) # here we found a short enough vector
    cand.sort(key=lambda x: linalg.norm(x)) # now we sort them by length
    Rxy=[cand[0]] # start with shortest zero vector
    for i in range(1,len(cand)):
        # if cand[i]+n0*a[0]+n1*a[1] is equal to existing Rxy, than they are equivalent
        unique = FindEquivalent2(cand[i], a, Rxy)
        if unique: # if cand[i] is unique, than save it to Rxy
            Rxy.append(cand[i])
    if (len(Rxy)!=Nm):
        print('ERROR: could not find all R points of the cluster a=',a)
        print('Rxy=', Rxy)
    Rxy=array(Rxy,dtype=int)
    # Now search for unique momentum points
    #print('bs[0]=', bs[0], 'bs[1]=', bs[1])
    Ks=[]
    loop = list(range(0,Nxy+1))+list(range(-1,-Nxy-1,-1)) # all possible integers to search for momentum points
    print('unit vectors=', a.tolist())
    print('real space vectors in the cluster:', Rxy.tolist())
    for ix in loop:
        for iy in loop:
            K = bs[0]*ix+bs[1]*iy  # small reciprocal vectors, which are within (1,0) and (0,1) are all K-points
            # any vector of the for (kx,ky) with |kx|<1 and |ky|<1 is cluster vector. Of course this is in units of 2*pi.
            if (abs(K[0])<1-1e-10) and (abs(K[1])<1-1e-10):
                dkx = array([np.sign(K[0]),0]) # reciprocal vector in first direction
                dky = array([0,np.sign(K[1])]) # reciprocal vector in second direction
                # we can guess three equivalent vectors in 2D, which might be withing |kx|<1 and |ky|<1
                K1eq = K-dkx  # subtracting reciprocal gives equivalent vector
                K2eq = K-dky  # equivalent
                K3eq = K-dkx-dky # equivalent
                #
                unique = True
                for k in Ks: # is K1eq or K2eq or K3eq in accepted Ks? If yes, this is not unique k-point
                    unique = unique and any(K1eq-k) and any(K2eq-k) and any(K3eq-k)
                if unique:
                    Ks.append(K)
                #print('K=', K, 'dkx=', dkx, 'dky=', dky, 'K1eq=', K1eq, 'K2eq=', K2eq, 'K3eq=', K3eq, 'unique=', unique)
    print('cluster momentum points using reciprocal vectors bs=', bs.tolist())
    for i,k in enumerate(Ks):
        print(i, k.tolist())
    if (len(Ks)!=Nm):
        print('ERROR: could not find all momentum vectors')
        print('Have Ks=', array(Ks))
    Ks=array(Ks)
    
    phi = dot(Ks,Rxy.T)*2*pi # all possible phases for Fourier transform
    expK = exp(phi*1j) # [ik,ir]
    if sum(abs(imag(expK)))<1e-10:# correct for the case of small cluster, all phases are real
        expK = array(expK.real, dtype=float)
    #print('e^{k*r}=', expK)
    # Finds what happens with Rxy for elementary shift in x direction
    shx=zeros(Nm,dtype=int)
    dx = array([1,0],dtype=int)
    for i in range(Nm):
        shR = Rxy[i]+dx # shifted vector for Rxy[i] for one lattice spacing
        shx[i]=FindEquivalent1(shR,Nm,a,Rxy) # what is the index of equivalent lattice vector
        #print('R['+str(i)+']=', Rxy[i], 'and +dx=',shR, 'which is equivalent to', shx[i], '=', Rxy[shx[i]]) 
    # Finds what happens with Rxy for elementary shift in y direction
    shy=zeros(Nm,dtype=int)
    dy = array([0,1],dtype=int)
    for i in range(Nm):
        shR = Rxy[i]+dy
        shy[i]=FindEquivalent1(shR,Nm,a,Rxy)

    if (positiveRxy): # Here we want to have all Rxy positive, hence finding equivalent positions
        for i in range(len(Rxy)):
            r=Rxy[i]
            if r[0]<0 or r[1]<0:
                for ixy in range(1,4):
                    ix = ixy%2
                    iy = int(ixy/2)
                    rn = r+a[0]*ix+a[1]*iy
                    if rn[0]>=0 and rn[1]>=0:
                        r=rn
                        break
                #print('R_old=', Rxy[i], 'R_new=', r)
                Rxy[i]=r

    # just temporary to be compatible with C++
    #if (Nm==8):
    #    Rxy = array([[0,0],[1,0],[0,1],[1,1],[2,0],[2,1],[1,2],[1,3]],dtype=int)
    #    shx = array([1,4,3,5,6,7,0,2],dtype=int)
    #    shy = array([2,3,4,6,5,0,7,1],dtype=int)
    #    Ks = array([[0,0],[0.5,0],[0,0.5],[0.5,0.5],[0.25,0.25],[-0.25,0.25],[0.25,-0.25],[-0.25,-0.25]])
    #    phi = dot(Ks,Rxy.T)*2*pi
    #    expK = exp(phi*1j) # [ik,ir]

    return (Nm, Rxy, Ks, expK, shx, shy)


##################################################
## Common to get_cix.py and supp_get_cix.py
##################################################
@njit
def sign(a):
    n = len(a)
    s = 0
    for i in range(n-1):
        for j in range(n-1-i):
            if (a[j+1] < a[j]):  # compare the two neighbors
                a[j],a[j+1] = a[j+1],a[j] # swap a[j] and a[j+1]
                s+=1
    return 1-2*(s%2)

@njit
def sign2(state, mask_min, mask_max):
    """ Here we check how many electrons are in state between mask_min and mask_max.
        If even (odd) number - sign is positive (negative).
    """
    #print('sign2 start with state={:08b} mmin={:08b} mmax={:08b}'.format(state, mask_min, mask_max))
    mask = mask_min<<1       # mask starts with mask_min shifted to the left by one site.
    n=0
    while (mask < mask_max): # and mask is increased until mask_max
        if (mask & state):   # electron exists in state at mask
            n+=1             # number of electrons in-between
        #print('mask={:08b} n={:d}'.format(mask,n))
        mask = mask<<1       # shift mask further to the left by one site
    sgn=1-2*(n%2)
    #print('sign=', sgn)
    return sgn
@njit
def sign3(state, mask_max):
    """ Here we check how many electrons are in state from the begining to mask_max.
        If even (odd) number - sign is positive (negative).
    """
    #print('sign2 start with state={:08b} mmin={:08b} mmax={:08b}'.format(state, mask_min, mask_max))
    mask = 1       # mask starts with very first position
    n=0
    while (mask < mask_max): # and mask is increased until mask_max
        if (mask & state):   # electron exists in state at mask
            n+=1             # number of electrons in-between
            #print(f"s3: mask={mask:>8b} n={n}")
        mask = mask<<1       # shift mask further to the left by one site
    sgn=1-2*(n%2)
    #print('sign=', sgn)
    return sgn
@njit
def sign4(state, mask_min, Nm):
    """ Here we check how many electrons are in state from the end to max_min.
        If even (odd) number - sign is positive (negative).
    """
    #print('sign2 start with state={:08b} mmin={:08b} mmax={:08b}'.format(state, mask_min, mask_max))
    mask = 1<<(2*Nm)       # mask starts with very first position
    n=0
    while (mask > mask_min): # and mask is increased until mask_max
        if (mask & state):   # electron exists in state at mask
            n+=1             # number of electrons in-between
            #print(f"s4: mask={mask:>8b} n={n}")
        mask = mask>>1       # shift mask further to the left by one site
    sgn=1-2*(n%2)
    #print('sign=', sgn)
    return sgn


def PrintStates(basis,masks,n,Sz2,ks=None,k=None):
    symbol = [['0','d'],['u','2']]
    Nm = int(len(masks)/2)
    for i,state in enumerate(basis):
        state_masks = array((state&masks).astype(bool), dtype=int) # gives 0 or 1 for each bin that is occupied or empty
        # knowing which beans are occupied, we know that if bin j is occupied and j+N is empty, we have up spin at j
        # if bin j and bin j+N is occupied, we have doubly occupied state, which we call '2'
        st = ''.join([symbol[state_masks[2*j]][state_masks[2*j+1]] for j in range(Nm)])
        if ks is not None:
            print(f"{i:3d} {st}  n={n}  Sz={0.5*Sz2} k={ks[i].tolist()}")
        elif k is not None:
            print(f"{i:3d} {st}  n={n}  Sz={0.5*Sz2} k={k}")
        else:
            print(f"{i:3d} {st}  n={n}  Sz={0.5*Sz2}")

def get_index(basis):
    ind_basis={}
    for i,state in enumerate(basis):
        ind_basis[state] = i
    return ind_basis

def Shifts2(state, ir, ups, dns, masks, shift):
    """ Given state encoded with bits, and one of the translations of the cluster,
        returns translated state, and the sign compared to original state.
        ir    -- index of translation in cluster positions Rxy
        shift -- shift[ir,i]=j means that position i is translated to position j under translation ir.
        ups   -- which sites are occupied with spin up
        dns   -- which sites are occupied with spin dn
        masks -- masks for all bits, starting with downs, and following with ups.
    """
    Nm = int(len(masks)/2) # number of sites, because first half are for down, and second half for up
    newstate = 0
    sgn=1
    if len(ups)>0:  # over all up spins in this state
        orderu=zeros(len(ups),dtype=int)   # what will be the order in which up spins are created
        for i0,i in enumerate(ups):        # original up spins were properly ordered
            # under shift R[ir] site i will be translated to site j.
            j = shift[ir,i] # this is jump R[ir]: i->j
            newstate = newstate^masks[j+Nm] # we put spin up at site j. Note up corresponds to mask[j+Nm]
            orderu[i0]=j  # remember the order in which we add up electrons
            #print('{:016b} {:016b}'.format(state,newstate), 'Rxy=', Rxy, 'i=',(i%Nm), 'j=', j)
        sgn *= sign(orderu) # contribution to the sign is the order in which we add up spins.
    if len(dns)>0:
        orderd=zeros(len(dns),dtype=int)
        for i0,i in enumerate(dns):
            j = shift[ir,i] # this is jump R[ir]: i->j
            newstate = newstate^masks[j]
            orderd[i0]=j
            #print('{:016b} {:016b}'.format(state,newstate), 'Rxy=', Rxy, 'i=',(i%Nm), 'j=', j)
        sgn *= sign(orderd)
    #print('state=', state, 'newstate=', newstate, 'orderu=', orderu, 'orderd=', orderd, 'sgn=', sgn)
    # To understand the sign, we note that in original state we assumed that electrons are added in ordered way, i.e.,
    # ... psi_3^\dagger psi_1^\dagger psi_0^\dagger|vacuum>
    # and if order in up spins (orderu) is not correctly ordered, or order in down spins is not correctly ordered,
    # there is a minus sign.
    return newstate, sgn

def FindCoupledBlocs(A, debug=False):
    """For a generic matrix A, which is somewhat sparse, we find blocks which can be diagonalized separately.
    If the matrix is block diagonal, than we will have multiple groups, and diagonalization is simplified. 
    If the matrix is not block diagonal, all columns and rows will be in one group.
    """
    group=[]
    remaining=list(range(shape(A)[0]))
    while len(remaining):
        i = remaining[0]
        c = where(abs(A[i,:])>1e-8)[0]
        if debug: print('c=', c, shape(c))
        if len(c)>0:
            c = set(c)
            c.add(i)
        else:
            c={i}
        if debug: print('set c=', c)
        
        remaining = [e for e in remaining if (e not in c)]
        added=-1
        todelete=[]
        for j in range(len(group)):
            #is_in_group_j = isin(c, group[j])
            #ais_in_group_j = any(is_in_group_j)
            #print('comparing c=',c,' with gr['+str(j)+']=', group[j], ' got: ', is_in_group_j)
            #if ais_in_group_j:
            if c.intersection(group[j]):
                if added==-1:
                    group[j] = group[j].union(c)
                else:
                    group[j] = group[j].union(group[added])
                    todelete.append(added)
                added=j
        if added==-1:
            group.append(c)
        if todelete:
            for j in todelete[::-1]:
                del group[j]
        if debug: print('i=', i, 'coupled with', c, 'group=', group, 'remain=', remaining)
    return group

def FindBlocks(NumPsis, Si):
    """
    Finds the index to later resort the transformation matrix and eigenvalues of S^2.
    The S^2 operator was diagonalized, but eigenvectors are not sorted in blocks yet. 
    Input:
        NumPsis[k]   -- contains the size of the block which has good momentum k
        Si[:]        -- is array of 2*S values, which are given in the current eigenbasis of good Sz,k,S^2.
    Output:
        grps[:][:]   -- the blocks with good k and S^2, and their current place in eigenvector order
        S2grp[:]     -- Spin 2*S after reordering
        newNumPsis[:]-- dimension of each block with good Sz,k,S^2
    """
    grps=[]
    S2grp=zeros(sum(NumPsis),dtype=int)
    newNumPsis=[]
    il=0
    for ik in range(len(NumPsis)):
        if NumPsis[ik]==0: continue
        Si_block = Si[il:il+NumPsis[ik]] # the current block with good momentum, and possible different spin
        S = Si_block[0]  # first spin in the block
        igrp={S:0}       # dictionary that gives index, whe key is spin
        grp=[[]]         # groups of indices with the same spin
        for l in range(NumPsis[ik]): # over all eigenvalues of this momentum
            S = Si_block[l]          # what is spin of this one
            if S in igrp:            # if this spin already has a group
                grp[igrp[S]].append(l) # just add this eigevalue index to the group with the same spin
            else:
                igrp[S]=len(grp)     # otherwise start a new group. The dictionary now has new spin S
                grp.append([l])      # and we have new group with the single eigenvalue with index l
        
        #print('Sib=', Si_block, 'grp=', grp, 'igrp=', igrp)
        # Now that we have groups with the same spin, we just flatten the groups, because now we only care about the order.
        flat_grp=[]
        for g in grp: flat_grp += g
        # we could calculate spin of each reordered list from igrp[S]. But it is easier to use input Si_block
        # and use fancy indexing in Python to reorder Si_block
        S2grp[il:il+NumPsis[ik]] = Si_block[flat_grp]
        # How big are the blocks now
        newNumPsis += [len(g) for g in grp]
        # remember the reordered list for each momentum
        grps.append( [i+il for i in flat_grp] )
        il += NumPsis[ik]
    return (grps,S2grp,newNumPsis)

def DiagonalizeBlocks(Hk, Tkin, groups, name=''):
    def RoundIfCorrect(x):
        if abs(x.imag)<1e-8 and abs(round(x.real*4)-(x*4))<1e-8:
            return round(x.real*4)/4.
        else:
            return x
    Ei=zeros(len(Tkin),dtype=float)
    Tkout = zeros(shape(Tkin), Tkin.dtype)
    for ig,g in enumerate(groups):
        if len(g)>1:
            gi = list(g)
            Hk_small = Hk[gi,:][:,gi]
            # Hk_small @ v == v @ np.diag(w)
            W, v = linalg.eigh(Hk_small)
            # TkE[l,:] = \sum_lp v.H[l,lp]*Tkin[lp,:]
            dTk = v.T.conj() @ Tkin[gi,:]
            Tkout[gi,:] = dTk
            W = [RoundIfCorrect(w) for w in W]
            Ei[gi] = W
            #Wp = [RoundIfCorrect(w) for w in W]
            if name: print('eigvals('+name+')['+str(gi)+']=',W)
        else:
            l = g.pop()
            Tkout[l,:] = Tkin[l,:]
            Ei[l] = RoundIfCorrect(Hk[l,l])
            if name: print('eigvals('+name+')['+str(l)+']=',Ei[l])
    return (Tkout, Ei)

def ComputeS2(basis, Nm, masks, ind_basis, Szi):
    """
    This routine computes S^2 for all states in the basis. We know that S^2 operator can only
    mix states in the basis.
    Essentially we evaluate S^2 = S_z^2 + 1/2 \sum_{i,j\in sites} (S_i^+ S_j^- + S_i^- S_j^+), 
    which can be written as S^2 = S_z^2 + 1/2 \sum_{i==j} (S_i^+ S_j^- + S_i^- S_j^+) + 
                                        +   \sum_{i\ne j} S_i^- S_j^+
    Input:
       basis     -- states in this N,Sz basis
       Nm        -- number of sites
       masks     -- masks for up & down spins in position basis
       ind_basis -- index that maps state with its index in the basis
       Szi       -- Sz for each bit on the state.
    """
    S2 = zeros((len(basis),len(basis)),dtype=float)
    for istate,state in enumerate(basis):
        # S^2=S_z^2+1/2\sum_{i}(S_i^+ S_i^- + S_i^- S_i^+) + \sum_{i!=j}(S_i^- S_j^+)
        # If site i has double-occupancy, S_i^+ and S_i^- operators vanish.
        # If site i has single-occupancy than 1/2*(S_i^+ S_i^- + S_i^- S_i^+) gives the same state back with factor 1/2.
        dd2=0
        masks_state = masks&state # here we apply all masks on a state at once, and later use them
        masks_state2 = masks_state.astype(bool).astype(int) # masks_state2[i] is 0 or 1 if masks[i] is set in state.
        for i in range(Nm):
            # this means that only up spin or only down spin is present at position i.
            if (masks_state2[2*i]+masks_state2[2*i+1]==1): dd2 += 0.5
        sz = 0.5*sum(masks_state2 * Szi)
        #masks_state2 = array(array(masks_state,dtype=bool),dtype=int)
        #print(istate, 'masks_state=', masks_state, masks_state2)
        S2[istate,istate] += (sz**2 + dd2) # diagonal contributions from S_z^2 and \sum_{i==j} (S_i)^2

        for i in range(Nm):
            if (masks_state2[2*i]==1 and masks_state2[2*i+1]==0): # S^-_i gives nonzero, because we have up at position i.
                isig = 1 # sign2(state,masks[2*i],masks[2*i+1])
                nstate = state^masks[2*i]^masks[2*i+1]            # flipping spin at i
                for j in range(Nm):
                    if (i!=j and masks_state2[2*j+1]==1 and masks_state2[2*j]==0): # S^+_j gives nonzero, because at position j is down spin.
                        jsig = 1 # sign2(nstate,masks[2*j],masks[2*j+1])
                        mstate = nstate^masks[2*j+1]^masks[2*j] # flipping also spin at j
                        #sts.append([mstate, isig*jsig])
                        jstate = ind_basis[mstate]
                        S2[istate,jstate] += isig*jsig
    return S2

if __name__=='__main__':
    #a=array([[1,0],[0,2]])     # N=2 site cluster
    a = array([[2,0],[0,2]])   # N=4 site cluster
    #a = array([[2,2],[-2,2]])   # N=8 site cluster
    #a = array([[4,0],[0,4]])   # N=16 site cluster
    
    (Nm, Rxy, Ks, expK, shx, shy) = FindClusterQuantities(a)
    print('shx=', shx)
    print('shy=', shy)
    
