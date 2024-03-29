#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
# 

from scipy import *
from scipy import linalg
import sys
import copy
import optparse

def union(data1, data2):
    " Takes a union of two lists"
    res = data1
    for d in data2:
        if d not in res: res.append(d)
    return res

def funique(data):
    ndata=[]
    for i in data:
        if i not in ndata: ndata.append(i)
    return ndata

class compare:
    "Compares two-site states in the basis"
    def __init__(self,baseNKS):
        self.baseNKS = baseNKS
    def __call__(self,x,y): # special function
        return int(self.baseNKS[x][0]-self.baseNKS[y][0])
    
def wdot(sta,stb):
    dsum=0.0
    for i in range(len(sta)):
        for j in range(len(stb)):
            #print 'sta,stb=', sta[i][0], stb[j][0]
            if sta[i][0] == stb[j][0]:
                dsum += sta[i][1]*stb[j][1]
    #if abs(dsum)>0:
    #    print 'dsum=', dsum, 'sta0=', sta, 'sbt0=', stb
    return dsum


def wdot0(sta,base,pairs,wsmall=1e-10):
    candidates=[]
    for i in range(len(sta)):
        new_candidates = pairs[sta[i][0]] 
        for k in new_candidates:
            if k not in candidates: candidates.append(k)
    
    wprod={}
    for candidate in candidates:
        w = wdot(sta,base[candidate])
        if abs(w)>wsmall:
            if candidate in wprod: wprod[candidate] += w
            else: wprod[candidate] = w
            
    return wprod



def ReadOneSiteCix(fcix):
    "Reading cix file for one site DMFT"
    f = open(fcix,'r')
    s0 = next(f) # CIX file for ctqmc!
    s1 = next(f) # cluster_size, number of states, number of baths, maximum_matrix_size
    (Nc0,Ns0,Nb0,Nm0) = list(map(int,f.next().split()))
    if Nc0!= 1:
        print('Wrong cix file. Input cix should be for single site!')
    s2 = next(f) # baths, dimension, symmetry
    baths=zeros((Nb0,3),dtype=int)
    for ib in range(Nb0):
        (iib,dim,ibs,ibg) = list(map(int,f.next().split()))
        if iib!=ib:
            print('Something wrong reading cix file (1)!')
        baths[ib]=(dim,ibs,ibg)
    #print baths
    s3 = next(f) # cluster energies for non-equivalent baths, eps[k]
    Nunique = max(baths[:,2])+1
    eps = list(map(float,f.next().split()))
    #print eps
    s4 = next(f) # N   K   Sz size

    mN = zeros(Ns0,dtype=int)
    mSz = zeros(Ns0,dtype=float)
    msize = zeros(Ns0,dtype=int)
    Fi = zeros((Ns0,Nb0),dtype=int)
    mEne = zeros((Ns0,Nm0),dtype=float)
    mS2 = zeros((Ns0,Nm0),dtype=float)
    for i in range(Ns0):
        data = f.next().split()
        ii = int(data[0])
        if ii!=i+1:
            print('Something wrong reading cix file (2)!')
        mN[i] = int(data[1])
        mSz[i] = float(data[3])
        msize[i] = int(data[4])
        for ib in range(Nb0):
            Fi[i,ib] = int(data[5+ib])-1
        for j in range(msize[i]):
            mEne[i,j] = float(data[5+Nb0+j])
        for j in range(msize[i]):
            mS2[i,j] = float(data[5+Nb0+msize[i]+j])

    s5 = next(f) # matrix elements

    Fm=[]
    for i in range(Ns0):
        fm=[]
        for ib in range(Nb0):
            fm.append([])
        Fm.append(fm)
        
    for i in range(Ns0):
        for ib in range(Nb0):
            data = f.next().split()
            ii = int(data[0])
            if ii!=i+1:
                print('Something wrong reading cix file (3)!')
            ij = int(data[1])
            if ij!=Fi[i,ib]+1:
                print('Something wrong reading cix file (4)!')
            if ij>0:
                sizei = int(data[2])
                sizej = int(data[3])
                if sizei!=msize[i]:
                    print('Something wrong reading cix file (5)!')
                if sizej!=msize[ij-1]:
                    print('Something wrong reading cix file (6)!')
                fm = zeros((sizei,sizej),dtype=float)
                cii=0
                for im1 in range(sizei):
                    for im2 in range(sizej):
                        fm[im1,im2]=data[4+cii]
                        cii+=1
                Fm[i][ib] = fm
    f.close()

    return (Nc0, Ns0, Nb0, Nm0, baths, eps, mN, mSz, msize, mEne, mS2, Fi, Fm)


def RemoveBaths(Nb0,Ns0,baths,eps,Ek,Fi,Fm, PRINT):
    """###########################################
       # removing baths which have energy > 1000 #
       ###########################################
    """
    bkeep=[]
    for ib,b in enumerate(baths):
        if abs(eps[b[1]])<1000.: bkeep.append(ib)
        print('baths= ', ib, b, eps[b[1]])
    
    n_Nb0 = len(bkeep)
    n_baths = zeros((n_Nb0,3),dtype=int)
    for i,ib in enumerate(bkeep):
        n_baths[i,:] = baths[ib,:]
    n_baths[:,1]-=min(n_baths[:,1])
    n_baths[:,2]-=min(n_baths[:,2])

    n_Nunique = max(n_baths[:,2])+1
    n_eps = zeros(n_Nunique,dtype=float)
    n_Ek = zeros(n_Nb0*2,dtype=float)
    for i,ib in enumerate(bkeep):
        n_eps[n_baths[i][1]] = eps[baths[ib][1]]
        #n_Ek[i] = Ek[ib]
        #n_Ek[i+n_Nb0] = Ek[ib+Nb0]
        n_Ek[2*i] = Ek[2*ib]
        n_Ek[2*i+1] = Ek[2*ib+1]
        
    n_Fi = zeros((Ns0,n_Nb0),dtype=int)
    n_Fm = [[] for i in range(Ns0)]
    for i in range(Ns0):
        for ib,b in enumerate(bkeep):
            n_Fi[i,ib] = Fi[i,b]
            n_Fm[i].append( Fm[i][b] )

    #Nb0 = n_Nb0
    #baths = n_baths
    #Nunique = n_Nunique
    #eps = n_eps
    #Fi = n_Fi
    #Fm = n_Fm
    #Ek = n_Ek

    if PRINT:
        print()
        print('baths=', n_baths)
        print('Nunique=', n_Nunique)
        print('eps=', n_eps)
        print('Ek=', n_Ek)
    
    return (n_Nb0,n_baths,n_eps,n_Ek,n_Fi,n_Fm)


def CorrectSingleSiteEnergies(Ns0,Nb0,msize,mN,mEne,baths,Uc,eps,Ek,Fi,Fm,PRINT=False):
    """ Adds Coulomb-U to one site energies.
        Creates energies for the two site problem.
    """
    Ek_local = [0.5*(Ek[2*ib]+Ek[2*ib+1]) for ib in range(Nb0)]  # from input Ek
    Ek_local_orig = [eps[baths[ib][1]] for ib in range(Nb0)]  # from origonal SS-cix file
    mE = [mEne[i] + Uc*mN[i]*(mN[i]-1)/2. for i in range(len(mEne))] # Adding Coulomb U!

    for ib in range(Nb0):  # Correcting single-particle energies
        R = [array([0]) for i in range(Ns0)]
        for i in range(Ns0):
            j = Fi[i,ib]
            if (j>=0):
                r = zeros(msize[j],dtype=float)
                for im1 in range(msize[i]):
                    for im2 in range(msize[j]):
                        r[im2] += Fm[i][ib][im1,im2]**2
                R[j]=r
        for i in range(len(R)):
            mE[i] += R[i]*(Ek_local[ib]-Ek_local_orig[ib])
    
    return array(mE)


def Create2SiteBasis(Ns0, mN, mSz, mEne, PRINT=False):
    """ Creates a basis for the two site problem from the basis of a one-site problem.
           --  base[int]=[ [(i1,i2),coeff1], [(i3,i4),coeff2],... ]
                 is a basis constructed as a direct product of one-site states,
                 but in bonding--anti-bonding representation, such that momentum is a good
                 quantum number (only (0,..), and (pi,...)
                 
                 Quantum mechanically, we would write
                    base[int] = coeff1*|i1>x|i2> + coeff2|i2>x|i3>+...
                 
                 Here (i1,i2) is a tuple of pointers to one-site basis, i.e., i1 is index in the one-site basis.
                 Coeff1 is the prefactor for the direct product.
                 
           --  baseNKS[int] = [N, K, Sz, Energy]
                 Here Energy contains only the on-site part, but not the hopping part.
                 The hopping part is added later after exact diagonalization.
                 
        It also creates an index pointer
    """ 
    sq2 = 1./sqrt(2.)
    base=[]
    baseNKS=[]
    for i1 in range(Ns0):
        #for i2 in range(i1,Ns0):
        for i2 in range(0,i1+1):
            if i1!=i2:
                (n1,n2) = (mN[i1],mN[i2])
                (sz1,sz2) = (mSz[i1],mSz[i2])
                # K=0
                Ene = mEne[i1,:]+mEne[i2,:]
                base.append( [ [(i1,i2),sq2],[(i2,i1),(-1)**(n1*n2)*sq2] ] )
                baseNKS.append( [n1+n2,0,sz1+sz2,Ene] )
                # K=Pi
                Ene = mEne[i1,:]+mEne[i2,:]
                base.append( [ [(i1,i2),sq2],[(i2,i1),-(-1)**(n1*n2)*sq2] ] )
                baseNKS.append( [n1+n2,1,sz1+sz2,Ene] )
            else:
                base.append( [[(i1,i1),1.0]] )
                baseNKS.append( [mN[i1]*2,0,mSz[i1]*2,mEne[i1,:]*2] )
                
    # index for sorting the base
    indbase = list(range(len(base)))
    # creates cmp function by class compare
    comp = compare(baseNKS)
    # sorts base
    indbase.sort(cmp=comp)
    wbase=[]
    wbaseNKS=[]
    for a in range(len(base)):
        wbase.append( base[indbase[a]] )
        wbaseNKS.append( baseNKS[indbase[a]] )
    base=wbase
    baseNKS=wbaseNKS

    # Index pointer from two-site base to two one-site bases.
    pairs={}
    for i in range(Ns0):
        for j in range(Ns0):
            pairs[(i,j)]=[]
    for i in range(len(base)):
        for a in range(len(base[i])):
            pairs[ base[i][a][0] ].append(i)


    if PRINT:
        print('base=')
        for i in range(len(base)):
            print(i, baseNKS[i], base[i])

    return (base, baseNKS, pairs)



def ComputeFdag(base,pairs,Nb0,mN,Fi,Fm,PRINT=False):
    sq2 = 1./sqrt(2.)
    # Creating F^dagger in this new even-odd basis
    Fp_K=[]
    for i in range(len(base)):
        st = base[i]
        #print 'st=', st
        fmn1=[]
        for ib in range(Nb0):
            fmn0=[]
            for ik in range(2):
                #print 'st=', st # starting with this state
                Fst=[]
                for a in range(len(st)):
                    (st1,st2) = st[a][0]  # composed of st1 at site 1 and st2 at site 2
                    n1 = mN[st1]          # number of electrons on site 1
                    if Fi[st1,ib]>0:      # adding electron at site 1
                        # The starting state is a product state of |n_1> at site 1 and state |n_2> at site 2.
                        # psi^+_k |n_1> x |n_2> = 1/sqrt(2)[ (psi^+_1|n_1>) x |n_2> + e^{ik} |n_1> x (psi^+_2|n_2> x |n_2> )] 
                        Fst.append( [(Fi[st1,ib],st2),st[a][1]*sq2*Fm[st1][ib][0,0]] )
                    if Fi[st2,ib]>0:      # adding electron at site 2
                        Fst.append( [(st1, Fi[st2,ib]),st[a][1]*sq2*Fm[st2][ib][0,0]*(-1)**n1*(-1)**ik] )
                
                wprod = wdot0(Fst,base,pairs)
                
                fmn0.append( wprod )
                #print 'prodct=', i, ib, ik, wprod, Fst
                
            fmn1.append(fmn0)    
        Fp_K.append(fmn1)

    # Creating F from F^dagger
    Fm_K=[]
    for i in range(len(base)):
        f2=[]
        for ib in range(Nb0):
            f1=[]
            for ik in range(2):
                f1.append({})
            f2.append(f1)
        Fm_K.append(f2)
    for i in range(len(base)):
        for ib in range(Nb0):
            for ik in range(2):
                fp = Fp_K[i][ib][ik]
                for a in list(fp.keys()):
                    Fm_K[a][ib][ik][i]=fp[a]


    if PRINT:
        print('F^\dagger=')
        for i in range(len(base)):
            for ib in range(Nb0):
                for ik in range(2):
                    print(i, ib, ik, '+', Fp_K[i][ib][ik])
        print('F = ')
        for i in range(len(base)):
            for ib in range(Nb0):
                for ik in range(2):
                    print(i, ib, ik, '-', Fm_K[i][ib][ik])
    
    return (Fp_K,Fm_K)

def ComputeOccup(base,Nb0,Fp_K,Fm_K,TEST=False):
    """Here we compute all_NN = F^dagger*F
       We also compute F*F^dagger, and we can check if F^dagger*F+F*F^dagger=1
    """
    all_NN=[]
    for i in range(len(base)):
        NN2=[]
        for ib in range(Nb0):
            NN1=[]
            for ik in range(2):
                # Here we compute F^dagger * F
                NN={}
                fm=Fm_K[i][ib][ik]
                for a in list(fm.keys()):
                    fpfm = Fp_K[a][ib][ik]
                    for b in list(fpfm.keys()):
                        if b in NN:NN[b]+=fpfm[b]*fm[a]
                        else: NN[b] = fpfm[b]*fm[a]
                NN1.append(NN)
                # Here we compute F * F^dagger
                MM={}
                fp=Fp_K[i][ib][ik]
                for a in list(fp.keys()):
                    fmfp = Fm_K[a][ib][ik]
                    for b in list(fmfp.keys()):
                        if b in MM: MM[b]+=fmfp[b]*fp[a]
                        else: MM[b]=fmfp[b]*fp[a]

                # Here we merge F^dagger * F + F * F^dagger and check identity
                if (TEST):
                    akeys = union(list(NN.keys()),list(MM.keys()))
                    print('##i=', i, 'akeys=', akeys) #?????
                    canticom={}
                    for k in akeys:
                        cc = 0
                        if k in NN: cc+=NN[k]
                        if k in MM: cc+=MM[k]
                        if abs(cc)>1e-10: canticom[k]=cc
                    print('f^+f+ff^+', i, ib, ik, canticom)
            NN2.append(NN1)
        all_NN.append(NN2)
    return all_NN






def CreateDiagPseudos(Nb0, Fp_K, all_NN, Ek, base, baseNKS, wsmall, PRINT=False):
    """ Here we first apply F^+_K to empty state, to generate all atomic states with
        occupancy unity. We than apply F^+_K to those with occupancy 1 and create
        all with occupancy 2, etc.
        This creates blocks of superstates: pseudo[i] = [[i0],[i1,i2], [i3,i4],...]
        which contains groups of states, such as two dimensional superstate [i1,i2].
        Such state [i1,i2] or [i3,i4] are current superstates.
        
        We than create indexes to these superstates. We show example of a one band model,
        where a link contains 16 states.

        In this case, pseudo contains the following data:
        [[0],[1],[2],[3],[4],[5],[6,9],[7,10],[8],[11],[12],[13],[14],[15]]

        and index table contains:
        blci1[0]=[0]
        ....
        blci1[6]=[6,9]
        blci1[7]=[7,10]
        blci1[8]=[8]
        blci1[9]=[11]
        blci1[10]=[12]
        blci1[11]=[13]
        blci1[12]=[14]
        blci1[13]=[15]

        and its inverse:
        blci[0]=0
        blci[1]=1
        blci[2]=2
        blci[3]=3
        blci[4]=4
        blci[5]=5
        blci[6]=6
        blci[7]=7
        blci[8]=8
        blci[9]=6
        blci[10]=7
        blci[11]=9
        blci[12]=10
        blci[13]=11
        blci[14]=12
        blci[15]=13

        Next we create Hamiltonian, using the on-site terms baseNKS[:][3] and hopping.
        The latter is approximated by intra-orbital terms only (within the same orbital).
        In this case, hopping is block diagonal in the basis of superstates, and it takes the form:
        Hopp[ip,jp]=(F^+_{b,K} F_{b,K})_{ip,jp}*Ek[K,b]

        This is because a superstate contains all terms that can be reached with operator
                <m,i|F_{b,K}^+ F_{b,K}|m,j>

        Finally, we diagonalize the two-site Hamiltonian.
        For convenience, we make sure the largest component of the eigenvector is positive (add minus sign if needed).
        
        
        
    """
    Ek_hop = zeros(len(Ek),dtype=float)
    for ib in range(Nb0):
        Ek_local = 0.5*(Ek[2*ib]+Ek[2*ib+1])
        for ik in range(2):
            Ek_hop[2*ib+ik] = Ek[2*ib+ik]-Ek_local
    
    blc_ind=0
    blci = zeros(len(base),dtype=int)
    blci1={0:[0]}
    pseudo=[[0]]
    Tr={0:array([[1.]])}
    Eham={0:[0.0]}
    for ll in range(2*Nb0):
        pseudon=[]
        pseudo_old = pseudo[:]
        # Here we start with all states at valence N-1, and from them we generate
        # states with valence N. We generate by applying F_k^+ to state at N-1.
        for pse in pseudo:
            for ib in range(Nb0):
                for ik in range(2):
                    block=[]
                    for ii in pse:
                        nblck = list(Fp_K[ii][ib][ik].keys())
                        # Instead of block += nblck
                        for k in nblck:
                            if k not in block: block.append(k)
                    if block:
                        
                        #block.sort()
                        #if block not in pseudon: pseudon.append(block)
                        
                        block = set(block) # create a set such that we can take union and intersection.
                        saved=False
                        for p in pseudon: # over all current superstates at this occupancy
                            if p & block: # does this superstate contain parts of the block
                                p = p | block # than we need to take a union
                                saved=True
                                break
                        if not saved:     # if it is entirely different superstate
                            pseudon.append(block) # we just save it in entirety

                        
        #pseudo = sorted(pseudon)
        pseudo = sorted([sorted(list(p)) for p in pseudon])
        if not pseudo: break

        # diagonalizing hopping term
        for ps in pseudo:
            # first create index tp later update Fp
            blc_ind+=1
            blci1[blc_ind]=ps
            for ip in ps:
                blci[ip] = blc_ind
            
            ### Hopping part of the Hamiltonian
            hamilt = zeros((len(ps),len(ps)),dtype=float)
            for ib in range(Nb0):
                for ik in range(2):
                    for i,ip in enumerate(ps):
                        for j,jp in enumerate(ps):
                            if jp in all_NN[ip][ib][ik]:
                                hamilt[i,j] += all_NN[ip][ib][ik][jp]*Ek_hop[2*ib+ik]
            ### Potential part of the Hamiltonian
            for i,ip in enumerate(ps):
                hamilt[i,i] += baseNKS[ip][3][0] # For now only one dimensional states are treeted????
                
            ee = linalg.eigh(hamilt)
            
            if (PRINT): print('Hamilt=', hamilt, 'diag=', ee)
            
            # Here we add a phase factor to eigenvectors, such that the largest component
            # of the eigenvector is positive.
            for i in range(len(ps)):
                v = ee[1][:,i]
                max_val = sorted(v,key= lambda x: -abs(x))[0]
                if (max_val<0): v *= -1.
            
            Tr[blc_ind]=ee[1]
            Eham[blc_ind]=ee[0]

            
            # Here we check if there is any degeneracy left. In this case we might be able to choose
            # eigenvectors more efficiently!!!!
            degene=[]
            i=0
            while i<len(ee[0]):
                degen=[]
                j=i+1
                for j in range(i+1,len(ee[0])):
                    if abs(ee[0][i]-ee[0][j])<wsmall:
                        degen.append(j)
                        #print 'WARNING: degeneracy'
                        #print 'rs=', ps, ee[0] #, hop #, matrix(Tr[blc_ind]).T*matrix(hop)*matrix(Tr[blc_ind])
                    else: break
                if degen:
                    degene.append([i]+degen)
                i=j
            if degene:
                #print 'degen=', ps, degene, ee[0]
                for degen in degene:
                    if len(degen)==2:
                        v0 = ee[1][degen[0]]
                        v1 = ee[1][degen[1]]
                        vs0 = (v0+v1)/sqrt(2.)
                        vs1 = (v0-v1)/sqrt(2.)
                        e0 = ee[0][degen[0]]
                        e1 = ee[0][degen[1]]

                        # Maybe I should not do that!
                        #ee[1][degen[0]]=vs0
                        #ee[1][degen[1]]=vs1
                        
                        #print 'v0=', v0, 'E0=', e0
                        #print 'v1=', v1, 'E1=', e1
                        #print 'v0=', vs0, 'E0=', e0
                        #print 'v1=', vs1, 'E1=', e1
                        #Tr[blc_ind]=ee[1]
                        #Eham[blc_ind]=ee[0]
                
                if (PRINT):
                    for ii in ps:
                        for ib in range(Nb0):
                            for ik in range(2):
                                print('d', ii, ib, ik, Fp_K[ii][ib][ik])
                            
    
    
    if PRINT: print('blci1=', blci1)

    # Finally, writing Energy in more convenient form for later printing
    Energy=zeros(len(base),dtype=float)
    EVector=[[] for i in range(len(base))]
    for k in list(blci1.keys()):
        for i1,j1 in enumerate(blci1[k]):
            Energy[j1] = Eham[k][i1]
            vec={}
            for i2,j2 in enumerate(blci1[k]):
                vec[j2] = Tr[k][i2,i1]
            EVector[j1] = vec
    if PRINT:
        for i in range(len(base)):
            print('Energy[',i,']=', Energy[i], 'EVector[',i,']=', EVector[i])

    
    return (pseudo, blci, blci1, Tr, Eham, Energy)


def Transform_F_To_Eigensybase(Fp_K, Tr, Nb0, base, blci, blci1, wsmall, PRINT):
    """On the input, the matrix elements of creation operator are written in original two-site
       basis, but not yet in the eigenbasis.
       We computed the eigenbasis in "CreateDiagPseudos", but did not yet transform F.
       Here we do the transformation, and return F written in the eigenbasis.
    """
    Fp_Knew=[]
    for ib in range(Nb0):
        fpKn2=[]
        for ik in range(2):
            fpKn1=[]
            for c in range(len(base)): fpKn1.append({})
            for k in list(blci1.keys()): # over all current superstates
                # finds to which states are the states in blci1 connected through Fp_K.
                blc_inds=[]
                for i in blci1[k]: # over all components of the current superstate k
                    blc_inds += [blci[ip] for ip in list(Fp_K[i][ib][ik].keys())]
                blc_inds = funique(blc_inds)
                blc_inds.sort()
                #print 'Superstate ', k, 'is connected to superstates', blc_inds, 'for ib=', ib, 'and ik=', ik
                for j in blc_inds: # over all components of superstate in sector N+1
                    if PRINT: print('F^+[ib=',ib,'ik=',ik, '] connects superstates ', k, 'with superstate', j, 'which are composed of the following basis states', blci1[k],blci1[j], ', respectively')
                    Fpm=zeros((len(blci1[k]),len(blci1[j])),dtype=float)
                    for i1,j1 in enumerate(blci1[k]):
                        for i2,j2 in enumerate(blci1[j]):
                            if j2 in Fp_K[j1][ib][ik]:
                                Fpm[i1,i2] = Fp_K[j1][ib][ik][j2]
                    Fpm_new = matrix(Tr[k]).T*matrix(Fpm)*matrix(Tr[j])
                    if PRINT:
                        print('Fpm_old=', Fpm)
                        print('Fpm_new=', Fpm_new)
                    for i1,j1 in enumerate(blci1[k]):
                        for i2,j2 in enumerate(blci1[j]):
                            if (abs(Fpm_new[i1,i2])>wsmall):
                                fpKn1[j1][j2] = Fpm_new[i1,i2]
            fpKn2.append(fpKn1)
        Fp_Knew.append(fpKn2)
    return Fp_Knew


def CreateFinalSuperstatesAndIndeces(base, Nb0, Fp_Knew, wsmall, PRINT=False):
    """Creates superstates from the form of F^+ in eigenbasis.
       It also creates indeces to these superstates.
    """
    wblc_ind=0
    wblci = [[] for i in range(len(base))]
    wblci[0] = (0,0)
    wblci1={0:[0]}
    
    pseudo=[[0]]
    all_pseudo=[]
    all_pseudo += pseudo
    # We will do the same operation (of superstate creation) as above,
    # but this time we have F^+ in eigenbase.
    # We again start with all states at valence N-1, and from them we generate
    # states with valence N. We generate by applying F_k^+ to state at N-1.
    for ll in range(2*Nb0):
        pseudon=[]
        for pse in pseudo:
            for ib in range(Nb0):
                for ik in range(2):
                    block=[]
                    for ii in pse:
                        nblck = list(Fp_Knew[ib][ik][ii].keys())
                        # Instead of block += nblck 
                        for k in nblck:
                            if k not in block and abs(Fp_Knew[ib][ik][ii][k])>wsmall:
                                block.append(k)
                    if block:
                        block = set(block) # create a set such that we can take union and intersection.
                        saved=False
                        for p in pseudon: # over all current superstates at this occupancy
                            if p & block: # does this superstate contain parts of the block
                                p = p | block # than we need to take a union
                                saved=True
                                break
                        if not saved:     # if it is entirely different superstate
                            pseudon.append(block) # we just save it in entirety
                        
                        #print 'll=', ll, 'pse=', pse, 'ib,ik=', 2*ib+ik, 'block=', block, pseudon

        pseudo = sorted([sorted(list(p)) for p in pseudon])
        if not pseudo: break
        all_pseudo += pseudo

        for ps in pseudo:
            wblc_ind+=1
            wblci1[wblc_ind]=ps
            #print 'ind=', wblc_ind, 'ps=', ps
            wblc_inside=0
            for ip in ps:
                wblci[ip] = (wblc_ind,wblc_inside)
                wblc_inside+=1

    if PRINT:
        print('wblci=', wblci)
        print('wblci1=', wblci1)
    return (all_pseudo, wblci, wblci1)



def CreateIndexFi2(all_pseudo, Nb0, Fp_Knew, wblci,wsmall):
    """We create also an index table Fi2 from Fp_Knew,
       which is convenient for printing of cix file
    """
    Fi2 = zeros((len(all_pseudo),Nb0,2),dtype=int)
    for i,p in enumerate(all_pseudo):
        for ib in range(Nb0):
            for ik in range(2):
                ifin=[]
                for ip in p:
                    fis = list(Fp_Knew[ib][ik][ip].keys())
                    ifin=[]
                    for t in list(Fp_Knew[ib][ik][ip].keys()):
                        if abs(Fp_Knew[ib][ik][ip][t])>wsmall:
                            ifin.append(wblci[t][0])
                        
                ifin = funique(ifin)
                if len(ifin)>1: print('ERROR: Not unique superstate F^+', ifin, 'ib,ik=', 2*ib+ik, 'p=', p)
                if ifin:
                    Fi2[i,ib,ik]=ifin[0]
                else:
                    Fi2[i,ib,ik]=-1
    return Fi2

def FindMaxsize(all_pseudo, PRINT=False):
    maxsize=0
    for i,p in enumerate(all_pseudo):
        if len(p)>maxsize: maxsize=len(p)
    if PRINT: print('maxsize=', maxsize)
    return maxsize

def Print_Fp_old_Fp_New(Fp_K, Fp_Knew, Nb0, base):
    "Printing Fp in two-site basis and in eigenbasis"
    for ib in range(Nb0):
        for ik in range(2):
            for i in range(len(base)):
                print('old_f^+', ib, ik, i, Fp_K[i][ib][ik])
    
    for ib in range(Nb0):
        for ik in range(2):
            for i in range(len(base)):
                print('f^+', ib, ik, i, Fp_Knew[ib][ik][i])


def CheckNKSconsistency(all_pseudo, baseNKS):
    " Just checking NKS consistency "
    for i,p in enumerate(all_pseudo):
        pn = copy(p)
        NKS = baseNKS[pn.pop()][:3] # Since we always mixed states with the same N,K,Sz, we can just take this from superstate
        for j in pn:
            if baseNKS[j][:3] != NKS[:3]:
                print('ERROR: Combining states which do not have the same NKS!')
                print('ERROR:', [baseNKS[j] for j in p])

def PrintHeader(lcix, bathk_ind, all_pseudo, Nb0, Fi2, Fj2, baseNKS, Energy, SUPERC=False, WithIndex=False):
    "Prints header only"
    for i,p in enumerate(all_pseudo):
        NKS = baseNKS[p[0]][:3] # Since we always mixed states with the same N,K,Sz, we can just take this from superstate
        Energ = [ Energy[j] for j in p]
        print("%3d " % (i+1), end=' ', file=lcix)
        if WithIndex: print("%3d " % (i+1), end=' ', file=lcix)
        print("%2d %2d %4.1f %2d " % (NKS[0], NKS[1], NKS[2], len(p)), end=' ', file=lcix)

        if (SUPERC):
            for (ib,ik) in bathk_ind:
                if (ib<Nb0/2):
                    print("%3d" % (Fi2[i,ib,ik]+1), end=' ', file=lcix)
                else:
                    print("%3d" % (Fj2[i,ib,ik]+1), end=' ', file=lcix)
        else:
            for (ib,ik) in bathk_ind:
                print("%3d" % (Fi2[i,ib,ik]+1), end=' ', file=lcix)
                
        print("  ", end=' ', file=lcix)
        for ei in Energ:
            print(ei, end=' ', file=lcix)
        print("  ", end=' ', file=lcix)
        for iq in range(len(p)):
            print(0, end=' ', file=lcix)
        print(file=lcix)


def PrintOneMatrixElement(lcix, i, p, ib, ik, Fi2, Fp_Knew, all_pseudo, wblci):
    ifinal = Fi2[i,ib,ik]
    print("%3d %3d " % (i+1, ifinal+1), end=' ', file=lcix)
    Mw=zeros((len(p),len(p)),dtype=float)
    if ifinal>=0:
        q = all_pseudo[ifinal]
        print("%2d %2d" % (len(p), len(q)), end=' ', file=lcix)
        
        Fp = zeros((len(p),len(q)),dtype=float)
        for ii,ip in enumerate(p):
            for iq in list(Fp_Knew[ib][ik][ip].keys()):
                if wblci[iq][0]!=ifinal: print('ERROR in ifinal!')
                jj = wblci[iq][1]
                Fp[ii,jj] = Fp_Knew[ib][ik][ip][iq]

        Mw = dot(Fp, transpose(Fp))
        
        for ii,ip in enumerate(p):
            for jj,iq in enumerate(q):
                print(Fp[ii,jj], end=' ', file=lcix)
    else:
        print("%2d %2d" % (0, 0), end=' ', file=lcix)
    print(file=lcix)
    return Mw
        
def PrintMatrixElements(lcix, bathk_ind, all_pseudo, Nb0, Fi2, Fj2, Fp_Knew, Fm_Knew, wblci, SUPERC):
    "Prints matrix elements for F^dagger operator"
    Nocc=[]  # occupancy
    for ibk in range(len(bathk_ind)):
        Nocc.append([list([]) for _ in range(len(all_pseudo))])
        
    for i,p in enumerate(all_pseudo):
        for ibk,(ib,ik) in enumerate(bathk_ind):
            if (SUPERC):
                if (ib<Nb0/2):
                    Mw =PrintOneMatrixElement(lcix, i, p, ib, ik, Fi2, Fp_Knew, all_pseudo, wblci)
                else:
                    Mw =PrintOneMatrixElement(lcix, i, p, ib, ik, Fj2, Fm_Knew, all_pseudo, wblci)
            else:
                Mw =PrintOneMatrixElement(lcix, i, p, ib, ik, Fi2, Fp_Knew, all_pseudo, wblci)
            Nocc[ibk][i]=identity(len(p))-Mw
    return Nocc

def PrintOccupancy(lcix, bathk_ind, Nocc, all_pseudo,baths):

    cases=sorted([baths[ib][1]*2+ik for (ib,ik) in bathk_ind])
    Nineq = cases[-1]+1
    bequal = [list([]) for _ in range(Nineq)]
    for ibk,(ib,ik) in enumerate(bathk_ind):
        bequal[baths[ib][1]*2+ik].append(ibk)
        
    for i,p in enumerate(all_pseudo):

        for ib,b in enumerate(bequal):
            print("%2d %2d %2d" % (i+1, len(p), len(p)), end=' ', file=lcix)

            nocc = zeros((len(p),len(p)),dtype=float)
            for ibk in b: # sum over equivalent
                nocc += Nocc[ibk][i][:,:]
            
            for ii in range(len(p)):
                for jj in range(len(p)):
                    print("%10.6f " % nocc[ii,jj], end=' ', file=lcix)
            print(file=lcix)
            
def PrintCixFile(outfile, bathk_ind, Fi2, Fj2, Fp_Knew, Fm_Knew, all_pseudo, baths, epsk2q, baseNKS, Energy, wblci, Nb0, maxsize, SUPERC):
    "Printing Cix file"
    lcix = open(outfile, 'w')
    print('# CIX file for ctqmc! ', file=lcix)
    print('# cluster_size, number of states, number of baths, maximum_matrix_size', file=lcix)

    if SUPERC:
        print(2, len(all_pseudo), len(bathk_ind)/2, maxsize, file=lcix)
        print('# baths, dimension, symmetry', file=lcix)
        
        sc_bath=[] # combine spin up and down into one entry
        for i in range(len(bathk_ind)/2):
            sc_bath.append( [bathk_ind[2*i],bathk_ind[2*i+1]] ) 

        # off-diagonal elements will come at the very end
        (ib,ik)=bathk_ind[-1]
        lastb_ind = baths[ib][1]*2+ik+1
        
        for i in range(len(sc_bath)):
            ib1 = sc_bath[i][0][0] # spin-up
            ib2 = sc_bath[i][1][0] # spin-dn
            ik = sc_bath[i][0][1]  # k is either 0 or pi
            #print 'ib=', ib1, 'ib2=', ib2, 'ik=', ik, 'cc=', baths[ib1][1]*2+ik, baths[ib2][1]*2+ik, 
            print(("%-2d" % i), '  ', baths[ib][0]*2, '  ', end=' ', file=lcix) # index, dimension
            if ik==0:  # select symmetry for for s+-
                off_index = lastb_ind
            else:
                off_index = -lastb_ind
                lastb_ind+=1
            print(baths[ib1][1]*2+ik, ("%2d "%off_index), ("%2d "%off_index), '-'+str((baths[ib2][1]*2+ik))+'*', '  ', end=' ', file=lcix)
            #print >> lcix, baths[ib1][2]*2+ik, '  # ik=', ik  # symmetry
    else:
        print(2, len(all_pseudo), len(bathk_ind), maxsize, file=lcix)
        print('# baths, dimension, symmetry', file=lcix)
        
        for i,(ib,ik) in enumerate(bathk_ind):
            print(("%-2d" % i), '  ', baths[ib][0], baths[ib][1]*2+ik, '  ', baths[ib][2]*2+ik, file=lcix)   #, '# ib=', ib, 'ik=', ik
        
    print('# cluster energies for non-equivalent baths, eps[k]', file=lcix)
    for ib in range(len(epsk2q)):
        print(epsk2q[ib], end=' ', file=lcix)
    if (SUPERC):
        for i in range(len(sc_bath)/2):
            print(0, end=' ', file=lcix)
    print(file=lcix)
    print('#     N  K  Sz size', file=lcix)
    
    PrintHeader(lcix, bathk_ind, all_pseudo, Nb0, Fi2, Fj2, baseNKS, Energy, SUPERC)
    print('# matrix elements', file=lcix)
    Nocc = PrintMatrixElements(lcix, bathk_ind, all_pseudo, Nb0, Fi2, Fj2, Fp_Knew, Fm_Knew, wblci, SUPERC)

    print('HB1', file=lcix)
    print('# number of operators needed', file=lcix)
    print('1', file=lcix)
    print('# Occupancy ', file=lcix)
    PrintOccupancy(lcix, bathk_ind, Nocc, all_pseudo,baths)
    print('# Data for HB1', file=lcix)

    if SUPERC:
        print(2, len(all_pseudo), len(bathk_ind)/2, maxsize, file=lcix)
    else:
        #print >> lcix, 2, len(all_pseudo), Nb0*2, maxsize
        print(2, len(all_pseudo), len(bathk_ind), maxsize, file=lcix)
        
    print('#      ind N  K  Sz  size', file=lcix)

    PrintHeader(lcix, bathk_ind, all_pseudo, Nb0, Fi2, Fj2, baseNKS, Energy, SUPERC, True)
    print('# matrix elements', file=lcix)
    PrintMatrixElements(lcix, bathk_ind, all_pseudo, Nb0, Fi2, Fj2, Fp_Knew, Fm_Knew, wblci, SUPERC)



def Test_FpF_Expensive(Nb0, Fp_Knew, base):
    # f_a^+ f_a + f_a f_a^+ = 1
    for ib in range(Nb0):
        for ik in range(2):
            Fpt = Fp_Knew[ib][ik]
            Fd = zeros((len(base),len(base)),dtype=float)
            for i in range(len(base)):
                wp = Fpt[i]
                for p in list(wp.keys()):
                    Fd[i,p]=wp[p]
            
            Fd = matrix(Fd)
            ID = Fd*Fd.T + Fd.T*Fd
            print('abs(1-F^+F+FF^+)', sum(abs(ID-identity(len(base)))))


def CheckIndexFi(Fi2):
    for ib in range(shape(Fi2)[1]):
        for ik in range(2):
            fi = Fi2[:,ib,ik].tolist()
            
            while -1 in fi:
                fi.remove(-1)
            
            if len(set(fi))!=len(fi):
                for p in fi:
                    if fi.count(p)>1: break
                print('ERROR: Pseudo', p, 'appears multiple times for ib=', ib, 'and ik=', ik)
                
    
def Transpose_Fp(Fp_Knew, all_pseudo, Nb0, wblci, wsmall):
    "Transposes F^\dagger to obtain F!"
    Fm_Knew=[]
    for ib in range(Nb0):
        fm2=[]
        for ik in range(2):
            fm=[{} for i in range(len(Fp_Knew[ib][ik]))]
            for i,p in enumerate(all_pseudo):
                for ii,ip in enumerate(p):
                    for iq in list(Fp_Knew[ib][ik][ip].keys()):
                        #print 'ip=', ip, 'iq=', iq, 'fp=', Fp_Knew[ib][ik][ip][iq]
                        #print 'len(fm)=', len(fm)
                        fm[iq][ip] = Fp_Knew[ib][ik][ip][iq]
            
            fm2.append( fm )
        Fm_Knew.append( fm2 )
        
    Fj2 = zeros((len(all_pseudo),Nb0,2),dtype=int)
    for i,p in enumerate(all_pseudo):
        for ib in range(Nb0):
            for ik in range(2):
                ifin=[]
                for ip in p:
                    fis = list(Fm_Knew[ib][ik][ip].keys())
                    ifin=[]
                    for t in list(Fm_Knew[ib][ik][ip].keys()):
                        if abs(Fm_Knew[ib][ik][ip][t])>wsmall:
                            ifin.append(wblci[t][0])
                        
                ifin = funique(ifin)
                if len(ifin)>1: print('ERROR: Not unique superstate F^+', ifin, 'ib,ik=', 2*ib+ik, 'p=', p)
                if ifin:
                    Fj2[i,ib,ik]=ifin[0]
                else:
                    Fj2[i,ib,ik]=-1
    
    return (Fm_Knew, Fj2)

def PrintSingleSiteCixFile(outfile,Uc,Nb0,Nm0,baths,eps,mN,mSz,msize,mEne,mS2,Fi,Fm,Ek):
    "Writing single-site cix file"

    def PrintHeader(WithIndex=False):
        #print 'Fi=', Fi
        for i in range(Ns0):
            print("%3d " % (i+1), end=' ', file=f)
            if WithIndex: print("%3d " % (i+1), end=' ', file=f)
            print("%2d %2d %4.1f %2d " % (mN[i], 0, mSz[i], msize[i]), end=' ', file=f)
            for ib in range(Nb0):
                print("%3d" % (Fi[i,ib]+1), end=' ', file=f)
            for j in range(msize[i]):
                print(mEne[i][j], end=' ', file=f)
            for j in range(msize[i]):
                print(mS2[i,j], end=' ', file=f)
            print(file=f)

    def PrintMatrixElm():
        for i in range(Ns0):
            for ib in range(Nb0):
                j = Fi[i,ib]
                print("%3d %3d " % (i+1, j+1), end=' ', file=f)
                if (j>=0):
                    print("%2d %2d" % (msize[i], msize[j]), end=' ', file=f)
                    for im1 in range(msize[i]):
                        for im2 in range(msize[j]):
                            print(Fm[i][ib][im1,im2], end=' ', file=f)
                else:
                    print("%2d %2d" % (0, 0), end=' ', file=f)
                print(file=f)
    
    Ns0 = len(mN)
    # Lets correct onsite energies with Ek
    Ek_local = [0.5*(Ek[2*i]+Ek[2*i+1]) for i in range(Nb0)]  # from input Ek
    ###Ek_local_orig = [eps[baths[i][1]] for i in range(Nb0)]  # from origonal SS-cix file

    print('##Ek_local     =', Ek_local)
    ###print '##Ek_local_orig=', Ek_local_orig
    print('##eps=', eps)
    print('##baths=', baths)
    print('##degi=', baths)

    
    # correcting eps with Ek
    eps_local=zeros(len(eps),dtype=float)
    degi=zeros(len(eps),dtype=int)
    for i in range(Nb0):
        eps_local[baths[i][1]]+=Ek_local[i]
        degi[baths[i][1]]+=1
    eps_local = [eps_local[i]/degi[i] for i in range(len(eps))]

    print('##eps_local=', eps_local)

    # Adding Coulomb term U to on-site energies
    ###mE = [mEne[i] + Uc*mN[i]*(mN[i]-1)/2. for i in range(len(mEne))]

    #print 'Ek_local=', Ek_local
    
    ###for ib in range(Nb0):
    ###    R = [array([0]) for i in range(Ns0)]
    ###    for i in range(Ns0):
    ###        j = Fi[i,ib]
    ###        if (j>=0):
    ###            r = zeros(msize[j],dtype=float)
    ###            for im1 in range(msize[i]):
    ###                for im2 in range(msize[j]):
    ###                    r[im2] += Fm[i][ib][im1,im2]**2
    ###            R[j]=r
    ###    #for i in range(Ns0): Nt[i] += R[i]
    ###    for i in range(len(R)):
    ###        mE[i] += R[i]*(Ek_local[ib]-Ek_local_orig[ib])
            
    f = open(outfile,'w')
    print('# CIX file for ctqmc!', file=f)
    print('# cluster_size, number of states, number of baths, maximum_matrix_size', file=f)
    print(1,Ns0,Nb0,Nm0, file=f)
    print('# baths, dimension, symmetry', file=f)
    for ib in range(Nb0):
        print(ib, baths[ib][0], baths[ib][1], baths[ib][2], file=f)
    print('# cluster energies for non-equivalent baths, eps[k]', file=f)
    for i in range(len(eps)):
        print(Ek_local[i], end=' ', file=f) 
    print(file=f)
    print('# N   K   Sz size', file=f)

    PrintHeader()
    print('# matrix elements', file=f)
    PrintMatrixElm()
    print('HB1', file=f)
    print('# number of operators needed', file=f)
    print(0, file=f)
    print('# Data for HB1', file=f)
    print(1,Ns0,Nb0,Nm0, file=f)
    print('# ind   N   K   Jz size', file=f)
    PrintHeader(True)
    print('# matrix elements', file=f)
    PrintMatrixElm()
    
                
def main(Uc, outfile, infile, Ek0, verbose, wsmall, ssc, SUPERC):

    
    #################################
    # Reading the one-site cix file #
    #################################
    (Nc0, Ns0, Nb0, Nm0, baths, eps, mN, mSz, msize, mEne, mS2, Fi, Fm) = ReadOneSiteCix(infile)


    # building Ek for each k-orbital-and-spin
    Ek = zeros(2*Nb0,dtype=float)
    for ib in range(Nb0):
        for ik in range(2):
            Ek[2*ib+ik] = Ek0[2*baths[ib][1]+ik]
    print('##Ek0=', Ek0)
    print('##Ek=', Ek)
    
    ###########################################
    # removing baths which have energy > 1000 #
    ###########################################
    (Nb0,baths,eps,Ek,Fi,Fm) = RemoveBaths(Nb0,Ns0,baths,eps,Ek,Fi,Fm, verbose>100)

    print('##Ek=', Ek)
    
    mEne = CorrectSingleSiteEnergies(Ns0,Nb0,msize,mN,mEne,baths,Uc,eps,Ek,Fi,Fm,verbose>110)
    
    PrintSingleSiteCixFile(ssc,Uc,Nb0,Nm0,baths,eps,mN,mSz,msize,mEne,mS2,Fi,Fm,Ek)

    #############################
    # Adding the Coulomb U term #
    #############################
    epsk2q = Ek0
    
    #print '##mEne=', mEne
    
    ###########################
    # Creates two-site basis  # 
    ###########################
    (base, baseNKS, pairs) = Create2SiteBasis(Ns0, mN, mSz, mEne, verbose>100)
    
    #############################################
    # Creates matrix elements of creation and   #
    # anhilation operator in the two-site basis #
    #############################################
    (Fp_K,Fm_K) = ComputeFdag(base,pairs,Nb0,mN,Fi,Fm,verbose>120)

    #########################################
    # Here we compute all_NN = F^dagger*F   #
    # We also compute F*F^dagger, and we    #
    # can check if F^dagger*F+F*F^dagger=1  #
    #########################################
    all_NN = ComputeOccup(base,Nb0,Fp_K,Fm_K,verbose>1000)
    
    ##################################################################################
    # Here we generate all superstates (pseudo) and index tables from original basis #
    # to this new basis. We also diagonalize the two-site problem and save           #
    # eigenvalues (Eham) and eigenvectors(Tr)                                        #
    ##################################################################################
    (pseudo, blci, blci1, Tr, Eham, Energy) = CreateDiagPseudos(Nb0, Fp_K, all_NN, Ek, base, baseNKS, wsmall, verbose>130)

    ####################################################################
    # Here we transform creation operator F to the two-side eigenbasis #
    ####################################################################
    Fp_Knew = Transform_F_To_Eigensybase(Fp_K, Tr, Nb0, base, blci, blci1, wsmall, verbose>1000)
    
    #########################
    # Just printing of F^+  #
    #########################
    if verbose>140: Print_Fp_old_Fp_New(Fp_K, Fp_Knew, Nb0, base)

    ######################################################################
    # We generate new superstates based on the form of F^+ in eigenbasis #
    # We also generate indeces to these superstates                      #
    ######################################################################
    (all_pseudo, wblci, wblci1) = CreateFinalSuperstatesAndIndeces(base, Nb0, Fp_Knew, wsmall, verbose>150)

    #########################
    # Just checking baseNKS #
    #########################
    if (verbose>1000): CheckNKSconsistency(all_pseudo, baseNKS)

    ##################################################
    # Index table Fi for fast rejection of QMC steps #
    # is constructed from Fp_Knew                    #
    ##################################################
    Fi2 = CreateIndexFi2(all_pseudo, Nb0, Fp_Knew, wblci,wsmall)

    CheckIndexFi(Fi2)
    
    ###############################
    # We need maximum matrix size #
    #   for printing              #
    ###############################
    maxsize = FindMaxsize(all_pseudo, verbose>120)


    (Fm_Knew, Fj2) = Transpose_Fp(Fp_Knew, all_pseudo, Nb0, wblci, wsmall)


    bathk_ind=[]
    NEW=True
    if NEW:
        for i in range(Nb0/2):
            bathk_ind.append([i,0])          # spin-up,k=0
            bathk_ind.append([Nb0/2+i,0])    # spin-dn,k=0
            bathk_ind.append([i,1])          # spin-up,k=pi
            bathk_ind.append([Nb0/2+i,1])    # spin-dn,k=pi
    else:
        for i in range(Nb0):
            bathk_ind.append([i,0])
            bathk_ind.append([i,1])
            


    #######################################################################
    #               Below is printing for ctqmc  solver                   #
    #######################################################################
    PrintCixFile(outfile, bathk_ind, Fi2, Fj2, Fp_Knew, Fm_Knew, all_pseudo, baths, epsk2q, baseNKS, Energy, wblci, Nb0, maxsize, SUPERC)
    #PrintCixFile(outfile, Fj2, Fm_Knew, all_pseudo, baths, epsk2q, baseNKS, Energy, wblci, Nb0, maxsize)
    
    if verbose>3000: Test_FpF_Expensive(Nb0, Fp_Knew, base)
    
if __name__ == '__main__':
    """
    Given an input cix file for a single-site DMFT (limited to ising case for the moment)
    creates cix file for two site CDMFT.
    
    If SUPERC=True, it creates input file for s+- superconductivity

    It is important to specify the on-site coulomb repulsion Uc,
    and the non-local part of the impurity levels Ek.

    For Ek, the order in the 5 band model is the following:
     Ek[0]=E_{x^2-y^2,k=0}
     Ek[1]=E_{x^2-y^2,k=pi}
     Ek[2]=E_{z^2,k=0}
     Ek[3]=E_{z^2,k=pi}
     Ek[4]=E_{xz,k=0}      
     Ek[5]=E_{xz,k=pi}
     Ek[6]=E_{yz,k=0}
     Ek[7]=E_{yz,k=pi}
     Ek[8]=E_{xy,k=0}       
     Ek[9]=E_{xy,k=pi}
    """

    usage = """usage: %prog [ options ]

    Given an input cix file for a single-site DMFT (limited to ising case for the moment)
    creates cix file for two site CDMFT.
    
    If SUPERC=True, it creates input file for s+- superconductivity

    It is important to specify the on-site coulomb repulsion Uc,
    and the non-local part of the impurity levels Ek.

    Ek should be in the followin format: Ek=[E1_{k=0},E1_{k=pi},E2_{k=0},E2_{k=pi},....]
    
    """

    parser = optparse.OptionParser(usage)
    parser.add_option("-E", "--Ek",  dest="Ek",  type="string", default="[0,0,0,0,0]", help="Impurity levels of the two site cluster Ek")
    parser.add_option("-U", "--Uc",  dest="Uc",  type="float", default=0.0, help="On-site Hubbard Coulomb repulsion U")
    parser.add_option("-o", "--out", dest="out", type="string", default='link_actqmc.cix', help="filename of the output 2-site cix file")
    parser.add_option("-i", "--inp", dest="inp", type="string", default='actqmc.cix', help="filename of the input single-site cix file")
    parser.add_option("-s", "--ssc", dest="ssc", type="string", default='SS_actqmc.cix', help="filename of the intermediate single-site cix file")
    parser.add_option("-S", "--SUPER", action="store_true", default=False, dest="SUPER", help="cix file superconductivity")
    # Next, parse the arguments
    (options, args) = parser.parse_args()

    Ek0 = eval(options.Ek)
    print('Ek=', Ek0)
    print('U=', options.Uc)
    print('inp=', options.inp)
    print('out=', options.out)
    print('ssc=', options.ssc)
    print('super=', options.SUPER)
    
    wsmall = 1e-10
    verbose = 100       # debug information

    #### for three band model
    #  [ E_{xz,k=0}, E_{xz,k=pi}, E_{yz,k=0], E_{yz,k=pi}, E_{xy,k=0}, E_{xy,k=pi} ]
    #Ek0 = [0, 0, 0.1, 0.1, 0.00922, 0, 0, -0.1, -0.1, -0.00922]
    #Ek0=[-0.48224573,0.38333269,0.40048207,0.95556531,-0.00999893,-0.00469948]
    #Ek0=[-0.71890552,0.1466729, 0.16382228,0.71890552,-0.24665872,-0.24135927]
    #Ek0=[-0.71890552,0.71890552,0.1466729,-0.24665872,0.16382228,-0.24135927]
    #Ek0=[0.,1.4378107,0.86557738,0.47224477,0.88272563,0.47754649]

    
    #### for one band model
    #Ek0 = [0.5,-0.5]  # Ek=[E_{K=0},E_{K=pi}]
    #infile = 'actqmc_one_band.cix'
    
    main(options.Uc, options.out, options.inp, Ek0, verbose, wsmall, options.ssc, options.SUPER)
    
