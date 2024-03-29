# @Copyright 2007 Kristjan Haule
# 
from scipy import *
from scipy import linalg
import copy
import sys

def StringToMatrix(cfstr):
    mm=[]
    for line in cfstr.split('\n'):
        line = line.strip()
        if line:
            data = array(map(float,line.split()))
            mm.append( data[0::2]+data[1::2]*1j )
    mm=array(mm)
    return mm

def PrintM(H,n=0):
    if type(H[0,0])==float64:
        for i in range(len(H)):
            for j in range(len(H[i])):
                print " "*n, "%11.8f   " % H[i,j],
            print
    else:
        for i in range(len(H)):
            for j in range(len(H[i])):
                print " "*n, "%11.8f %11.8f   " % (H[i,j].real, H[i,j].imag),
            print



def FindMixing(Delta):
    w,V=linalg.eigh(Delta)
    # Delta_diagonal = V^+ Delta V
    V = transpose(V)
    #print 'eigenvalues of Delta=', w

    print 'beginning-V='
    PrintM(V)
    
    # We sort eigenvectors such that they are as close as possible to unity
    U=zeros(shape(V),dtype=float)
    for i in range(len(V)):
        ind=argmax(abs(V[:,i]))
        sn = sign(V[ind,i].real)
        U[i,:] = sn*V[ind,:].real 
        V = delete(V, ind, 0)
        #PrintM( V )
        #print
    # Checking that V was real, and hence U is unitary
    if sum(sum(abs(matrix(U).T * matrix(U)-identity(len(U)))))>1e-5:
        print 'It seems you have mixing which is not real! Not supported yet!'
        sys.exit(0)
        
    # Zeroing small matrix elements, which we will neglect
    for i in range(len(U)):
        for j in range(len(U[i])):
            if abs(U[i,j])<1e-3: U[i,j]=0
    # Renormalizing the remaining matrix, such that it is unitary again
    (u_,s_,v_) = linalg.svd(U)
    #print 'singular-values=', s_
    U = dot(u_ , v_)
    #print 'U='
    #PrintM(U)
    # Saving U in condensed way
    Vmix=[[] for i in range(len(U))]
    for i in range(len(U)):
        for j in range(len(U[i])):
            if abs(U[i,j])>=1e-3:
                Vmix[i].append( (j,U[i,j] ) )
    return Vmix

def FindMixing0(Olap):
    Vmix=[[] for i in range(len(Olap))]
    for i in range(len(Olap)):
        for j in range(len(Olap[i])):
            if abs(Olap[i,j])>=1e-3:
                if abs(Olap[i,j]-1.)<1e-4: Olap[i,j]=1.
                Vmix[i].append( (j,Olap[i,j].real ) )
    return Vmix


def ReadOneSiteCix(fcix):
    "Reading cix file for one site DMFT"
    f = open(fcix,'r')
    s0 = f.next() # CIX file for ctqmc!
    s1 = f.next() # cluster_size, number of states, number of baths, maximum_matrix_size
    (Nc0,Ns0,Nb0,Nm0) = map(int,f.next().split())
    if Nc0!= 1:
        print 'Wrong cix file. Input cix should be for single site!'
    s2 = f.next() # baths, dimension, symmetry
    baths=zeros((Nb0,3),dtype=int)
    for ib in range(Nb0):
        (iib,dim,ibs,ibg) = map(int,f.next().split())
        if iib!=ib:
            print 'Something wrong reading cix file (1)!'
        baths[ib]=(dim,ibs,ibg)
    #print baths
    s3 = f.next() # cluster energies for non-equivalent baths, eps[k]
    Nunique = max(baths[:,2])+1
    eps = map(float,f.next().split())
    #print eps
    s4 = f.next() # N   K   Sz size

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
            print 'Something wrong reading cix file (2)!'
        mN[i] = int(data[1])
        mSz[i] = float(data[3])
        msize[i] = int(data[4])
        for ib in range(Nb0):
            Fi[i,ib] = int(data[5+ib])-1
        for j in range(msize[i]):
            mEne[i,j] = float(data[5+Nb0+j])
        for j in range(msize[i]):
            mS2[i,j] = float(data[5+Nb0+msize[i]+j])

    s5 = f.next() # matrix elements

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
                print 'Something wrong reading cix file (3)!'
            ij = int(data[1])
            if ij!=Fi[i,ib]+1:
                print 'Something wrong reading cix file (4)!'
            if ij>0:
                sizei = int(data[2])
                sizej = int(data[3])
                if sizei!=msize[i]:
                    print 'Something wrong reading cix file (5)!'
                if sizej!=msize[ij-1]:
                    print 'Something wrong reading cix file (6)!'
                fm = zeros((sizei,sizej),dtype=float)
                cii=0
                for im1 in range(sizei):
                    for im2 in range(sizej):
                        fm[im1,im2]=data[4+cii]
                        cii+=1
                Fm[i][ib] = fm
    f.close()

    return (Nc0, Ns0, Nb0, Nm0, baths, eps, mN, mSz, msize, mEne, mS2, Fi, Fm)


def CorrectSingleParticleEnergies(Ns0,Nb0,msize,mEne,baths,eps,Fi,Fm):
    """ Adds Coulomb-U to one site energies.
        Creates energies for the two site problem.
    """
    E_local_orig = [eps[baths[ib][1]] for ib in range(Nb0)]  # from origonal SS-cix file
    for ib in range(Nb0):  # Correcting single-particle energies
        R = zeros( shape(mEne) )
        for i in range(Ns0):
            j = Fi[i,ib]
            if (j>=0):
                for im1 in range(msize[i]):
                    for im2 in range(msize[j]):
                        R[j,im2] += Fm[i][ib][im1,im2]**2
        
        for i in range(len(R)):
            for im in range(msize[i]):
                mEne[i,im] -= R[i,im]*E_local_orig[ib]
    
    return mEne

def Simplify(st):
    for i in range(len(st)):
        for j in range(i+1,len(st)):
            if st[i][:2]==st[j][:2]:
                st[i][2] += st[j][2]
                st[j][2]=0.0
    
    newst = filter(lambda x: abs(x[2])>1e-5, st)
    return newst

def dot_state(st1,st2):
    # state is stored as (index-of-superstate,index-within-superstate,coefficient)
    sm=0
    for i,im,ic in st1:
        for j,jm,jc in st2:
            if i==j and im==jm:
                sm += ic*jc
    return sm

def StateSum(st1,st2):
    return Simplify(st1+st2)
    

def IsDifferent(newst,stNp):
    for st in stNp:
        overlap=dot_state(newst,st)
        #print 'overlap newst,st=', overlap
        if abs(abs(overlap)-1)<1e-3: return False
        if abs(overlap)>1e-3:
            print 'This should not happen! We have two non-orthogonal states', st, newst
    return True

def ArrangePhase(st):
    # We want the largest entry to be positive
    ind=argmax([abs(a[2]) for a in st])
    if st[ind][2]<0:
        for i in range(len(st)): st[i][2]*=-1 # change the sign of all coefficients
    return st



def CreateFp(xgroup, xgroupold,Um,Umold,Fi,Fm,Nb0):
    small=1e-4
    
    #invgroup={}
    #for ig,group in enumerate(xgroup):
    #    for j in group:
    #        invgroup[j[0]]=ig

    Fmr={}
    for igold,gold in enumerate(xgroupold):
        Uold=matrix(Umold[igold])
        for ig,g in enumerate(xgroup):
            U=matrix(Um[ig])
            Fnew=[]
            for ib in range(Nb0):
                FF = zeros((len(gold),len(g)))
                #print 'gold=', gold, 'g=', g
                for iio,io in enumerate(gold):
                    j=Fi[io[0],ib]
                    if j>=0:
                        for ii,i in enumerate(g):
                            #print 'j=', j, invgroup[j], 'i=', i
                            if j==i[0]:
                                #print 'Success'
                                FF[iio,ii] += Fm[io[0]][ib][io[1],i[1]]
                #print 'Fnew[ib='+str(ib)+',igo='+str(igold)+',ig='+str(ig)+']=', FF.tolist()
                ### Uold[ieold,jmold]*FF[jmold,jmnew]*U[jnnew,ienew]
                Fnew_ = Uold*FF*U.T
                Fnew.append(Fnew_)
                #print 'Fnew[ib='+str(ib)+',igo='+str(igold)+',ig='+str(ig)+']=', Fnew_.tolist()

            for ib in range(Nb0):
                Fn_ = zeros(shape(Fnew[ib]))
                for iib,mix in Vmix[ib]:
                    Fn_ += Fnew[iib]*mix
                if sum(sum(abs(Fn_)))>1e-4:
                    #print 'Fnew[ib='+str(ib)+',igo='+str(igold)+',ig='+str(ig)+']=', Fn_.tolist()
                    if Fmr.has_key((ib,igold,ig)):
                        print 'ERROR: Superstates are not unique', (ib,igold,ig)
                    for ix in range(shape(Fn_)[0]):
                        for iy in range(shape(Fn_)[1]):
                            if abs(Fn_[ix,iy])<small: Fn_[ix,iy]=0.0
                    Fmr[(ib,igold,ig)] = Fn_
                    
    for gb in  sorted(Fmr.keys(),key=lambda x: x[2]):
        (ib,igold,ignew)=gb
        #ignew=Fir[gb]
        print 'Fnew[ib='+str(ib)+',igo='+str(igold)+',ign='+str(ignew)+']=', Fmr[gb].tolist()
    print
    return Fmr


def CreateEigenstates(sts,mEne,Enew):
    def overlap(ta,tb):
        for t in range(ta):
            if t in tb: return True
        return False
    
    groups=[]
    ingroup={}
    for i in range(len(sts)):
        found=False
        state=sts[i]
        includes = set([state[m][0] for m in range(len(state))])
        for l in range(len(groups)):
            if includes & groups[l]:
                ingroup[l].append(i)
                groups[l]= groups[l] | includes
                found=True
                break
        if not found:
            groups.append( includes )
            ingroup[len(groups)-1] = [i]
    #print 'groups=', groups
    #print 'ingroup=', ingroup
    
    AEnergies=[]
    AStates=[]
    xgroup=[]
    Um=[]
    #AU=[]
    #egroup= [ [] for ig in range(len(groups))]
    
    for ig in sorted(ingroup.keys()):
        
        #print 'group: ', ig, ' contains:'
        #for ist in ingroup[ig]:
        #    print '   ', sts[ist]
        
        jmBasis=set()
        for ist in ingroup[ig]:
            jmBasis.update( set([tuple(pst[:2]) for pst in sts[ist]]) )
        jmBasis=list(jmBasis)
        jmBasis = sorted(jmBasis, key=lambda x: x[0])
        jmBasis = sorted(jmBasis, key=lambda x: x[1])

        xgroup.append( jmBasis )
        size=len(jmBasis)
        
        
        if len(ingroup[ig]) != len(jmBasis):
            print 'It should not happen! The dimension of old and new basis should be the same'
            print ingroup[ig]
            print jmBasis
            sys.exit(0)
            
        Tr=zeros((size,size))
        for ii,ist in enumerate(ingroup[ig]):
            for pst in sts[ist]:
                jj = jmBasis.index( tuple(pst[:2]) )
                Tr[ii,jj] = pst[2]
        Tr=matrix(Tr)
        U, s, Vh = linalg.svd(Tr)
        Tr = matrix(U)*matrix(Vh)
        
        #print 'How good='
        #print sum(sum(abs(Tr*Tr.T)-identity(size)))
        
        if sum(sum(abs(Tr*Tr.T)-identity(size)))>1e-4:
            print 'ERROR: Transformation matrix should be unitary!'
            print Tr
            sys.exit(0)
            
        #print 'jmBasis=', jmBasis
        #print 'Tr=', Tr

        # Hamo is the energy from the cix file withouth crystal fields, i.e., Hunds only
        Hamo=zeros((size,size))
        for ii,(i,im) in enumerate(jmBasis):
            Hamo[ii,ii] = mEne[i,im]
        # Hamn is due to crystal fields, which will be presented in (j,mj) basis
        Hamn=zeros((size,size))
        for ii,ist in enumerate(ingroup[ig]):
            Hamn[ii,ii] = Enew[ist]
        Hamn = Tr.T * Hamn * Tr  # To (j,mj) basis. Note that Tr_{alpha,mj}
        Ham = Hamn+Hamo          # total Hamiltonian
        Es,Vs=linalg.eigh(Ham)
        Us=Vs.T
        
        # Here we arrange the phase such that the largest contribution is positive
        for i in range(len(Us)):
            ind=argmax(abs(Us[i,:]))
            if Us[i,ind]<0: Us[i,:] *= -1
        
        Ua=Us*Tr.T
        
        #print 'Energies=', Es
        #print 'States in (j,m) basis:'
        #PrintM(Us)
        #print 'States in psi^+_alpha basis:'
        #PrintM(array(Ua))
        #print
    
        #AU.append(Ua)
        Um.append(Us)
        
        for ie in range(size):
            AEnergies.append( Es[ie] )
            state=[]
            for j,(i,im) in enumerate(jmBasis):
                if abs(Us[ie,j])>1e-6:
                    state.append( [i,im,Us[ie,j]] )
            AStates.append(state)
            #egroup[ig].append( len(AStates)-1 )
    print 'xgroup=', xgroup
    
    return (AStates,AEnergies,xgroup,Um)


def CreateNewStates(Nc0, Ns0, Nb0, Nm0, baths, eps, mN, mSz, msize, mEne, mS2, Fi, Fm,Eimp):
    maxN = max(mN)
    
    FM={}
    STS=[[[0,0,1]]]
    ENE=[0]

    glst=0
    gpseudo=[[glst]]
    glst+=1
    stN=[ [[0,0,1]] ] # stored as (index,coefficient)
    Enew=[0]
    #ingroupold={0: [0]}
    #AStatesold=[[[0,0,1]]]
    #egroupold=[[0,],]
    xgroupold=[[(0,0)]]
    Umold=[identity(1)]
    gshiftold=0
    for n in range(1,maxN+1):
        stNp=[]
        Enewp=[]
        father=[]
        for ist,statea in enumerate(stN):
            for b in range(Nb0):
                newst=[]
                for i,im,cm in statea:
                    for ind,vb in Vmix[b]:
                        j=Fi[i,ind]
                        if j>=0:
                            for im2 in range(msize[j]):
                                newst.append( [j,im2,vb*cm*Fm[i][ind][im,im2]] )
                newst = Simplify(newst)
                #print 'n=', n, 'applying F['+str(b)+'] on |'+str(ist)+'> and got ', newst
                if newst and IsDifferent(newst,stNp):
                    newst = ArrangePhase(newst)
                    stNp.append(newst)
                    Enewp.append(Enew[ist]+Eimp[b])
                    father.append( (ist,b) )
        #print 'Finished level N'
        print 'States at particle number '+str(n)+' are:'
        print "%2s %13s %2s %2s" % ('#','Energy', 'fi', 'fb')
        for ii,si in enumerate(stNp):
            print "%2d %13.10f %2d %2d" % (ii, Enewp[ii], father[ii][0], father[ii][1]),
            print si
        
        stN = stNp
        Enew = Enewp
        (AStates, AEnergies,xgroup,Um) = CreateEigenstates(stNp,mEne,Enewp)

        #print 'egroupold=', egroupold
        #print 'egroup=', egroup
        #print 'group=', groups
        gshift = len(gpseudo)
        for grp in xgroup:
            gpseudo.append( (arange(len(grp))+glst).tolist() )
            glst += len(grp)
            
        print 'Eigenstates at particle number '+str(n)+' are:'
        for ii,si in enumerate(AStates):
            print "%2d %13.10f" % (ii, AEnergies[ii]),
            print si
        print

        
        Fmr = CreateFp(xgroup, xgroupold,Um,Umold,Fi,Fm,Nb0)
        
        print 'gshift=', gshift, 'gshiftold=', gshiftold
        print 'gpseudo=', gpseudo

        
        for bg in Fmr.keys():
            FM[(bg[0],bg[1]+gshiftold,bg[2]+gshift)] = Fmr[bg]

        STS += AStates
        ENE += AEnergies
        
        #ingroupold=ingroup
        #AStatesold=AStates
        #egroupold=egroup
        xgroupold=xgroup
        Umold=Um
        gshiftold=gshift
        
        
        
    #print 'stNp=', stNp

    return (STS,ENE,FM,gpseudo)










        

#def PrintOccupancy(lcix, bathk_ind, Nocc, all_pseudo,baths):
#
#    cases=sorted([baths[ib][1]*2+ik for (ib,ik) in bathk_ind])
#    Nineq = cases[-1]+1
#    bequal = [list([]) for _ in range(Nineq)]
#    for ibk,(ib,ik) in enumerate(bathk_ind):
#        bequal[baths[ib][1]*2+ik].append(ibk)
#        
#    for i,p in enumerate(all_pseudo):
#
#        for ib,b in enumerate(bequal):
#            print >> lcix, "%2d %2d %2d" % (i+1, len(p), len(p)),
#
#            nocc = zeros((len(p),len(p)),dtype=float)
#            for ibk in b: # sum over equivalent
#                nocc += Nocc[ibk][i][:,:]
#            
#            for ii in range(len(p)):
#                for jj in range(len(p)):
#                    print >> lcix, "%10.6f " % nocc[ii,jj],
#            print >> lcix
            
def PrintCixFile(outfile, STS,ENE,FM,gpseudo, Eimp, mN, mSz,mS2,Vmix,baths):

    def PrintHeader(lcix, WithIndex, gpseudo, FM, ENE, STS, Vmix, eps, nbaths, dim, mN, mSz, mS2):
        if WithIndex:
            print >> lcix, '#     ind  N  K  Sz size'
        else:
            print >> lcix, '#     N  K  Sz size'
        
        for i,gp in enumerate(gpseudo):
            Jz=[]
            J2=[]
            n=0
            for l,g in enumerate(gp):
                jz=0
                j2=0
                for st in STS[g]:
                    n=mN[st[0]]
                    jz += mSz[st[0]]*st[2]**2
                    j2 += mS2[st[0]][st[1]]*st[2]**2
                Jz.append(jz)
                J2.append(j2)
                
            print >> lcix, "%3d " % (i+1),
            if WithIndex: print >> lcix, "%3d " % (i+1),
            print >> lcix, "%2d %2d %5.2f %2d " % (n, 0, Jz[0], len(gp)),
            
            for ib in range(Nb0):
                j=-1
                for ign in range(len(gpseudo)):
                    bg = (ib,i,ign)
                    if FM.has_key( bg ):
                        j=ign
                        break
                print >> lcix, "%3d" % (j+1),
            
            print >> lcix, "  ",
            for l,g in enumerate(gp):
                print >> lcix, "%14.8f" % ENE[g],
            print >> lcix, "   ",
            for l,g in enumerate(gp):
                print >> lcix, "%6.4f" % J2[l],
            print >> lcix, ' # ',
            for l,g in enumerate(gp):
                print >> lcix, STS[g],
            print >> lcix

    def PrintMatrixElements(lcix,gpseudo,FM,Nb0):
        for ig,gp in enumerate(gpseudo):
            for ib in range(Nb0):
                ifinal=-1
                for ign in range(len(gpseudo)):
                    bg = (ib,ig,ign)
                    if FM.has_key( bg ):
                        ifinal=ign
                        #print 'xig=', ig, 'ib=', ib, 'ifinal=', ifinal, 'bg=', bg
                        break
                        
                print >> lcix, "%3d %3d " % (ig+1, ifinal+1),
                if ifinal>=0:
                    q = gpseudo[ifinal]
                    print >> lcix, "%2d %2d" % (len(gp), len(q)),
                    #print 'ig=', ig, 'ib=', ib, 'ifinal=', ifinal, 'bg=', bg
                    #if shape(FM[bg])!= (len(gp),len(q)):
                    #    print "ERROR SHAPES ARE NOT CORRECT!"
                    for ii in range(len(gp)):
                        for jj in range(len(q)):
                            print >> lcix, "%11.8f " % FM[bg][ii,jj],
                else:
                    print >> lcix, "%2d %2d" % (0, 0),
                print >> lcix
    ###################
    # Start routine   #
    ###################
    dim = [len(STS[i]) for i in range(len(STS))]
    #print 'dim=', dim
    nbaths=[]
    for ib in range(len(Vmix)):
        ibdominant = sorted(Vmix[ib],key=lambda x: -abs(x[1]))[0][0]
        nbaths.append( baths[ibdominant] )
    nbaths=array(nbaths)
    nsym=max(nbaths[:,1])+1
    
    eps = zeros(nsym)
    deg = zeros(nsym)
    for ib in range(len(Vmix)):
        sym=nbaths[ib][1]
        eps[sym] += Eimp[ib]
        deg[sym] += 1
    eps=eps/deg
    
    "Printing Cix file"
    lcix = open(outfile, 'w')
    print >> lcix, '# CIX file for ctqmc! '
    print >> lcix, '# cluster_size, number of states, number of baths, maximum_matrix_size'

    print >> lcix, 1, len(gpseudo), len(Vmix), max(dim)
    print >> lcix, '# baths, dimension, symmetry'
        
    for ib in range(len(Vmix)):
        print >> lcix, ("%-2d" % ib), '  ', nbaths[ib][0], nbaths[ib][1], '  ', nbaths[ib][2],   '# ib=', Vmix[ib]
        
    print >> lcix, '# cluster energies for non-equivalent baths, eps[k]'
    for ib in range(len(eps)):
        print >> lcix, eps[ib],
    print >> lcix

    PrintHeader(lcix, False, gpseudo, FM, ENE, STS, Vmix, eps, nbaths, dim, mN, mSz, mS2)
    print >> lcix, '# matrix elements'
    PrintMatrixElements(lcix,gpseudo,FM,Nb0)
        
    print >> lcix, 'HB1'
    print >> lcix, '# number of operators needed'
    print >> lcix, '1'
    print >> lcix, '# Occupancy '
    
    for ig,gp in enumerate(gpseudo):
        aocc=zeros((nsym,len(gp),len(gp)))
        ncheck=zeros((len(gp),len(gp)))
        for ib in range(Nb0):
            j=-1
            Mw=zeros((len(gp),len(gp)))
            for ign in range(len(gpseudo)):
                bg = (ib,ig,ign)
                if FM.has_key( bg ):
                    j=ign
                    break
            if j>=0:
                Mw = dot(FM[bg],FM[bg].T)
            occ = identity(len(gp))-Mw

            sym=nbaths[ib][1]
            aocc[sym] += occ
            ncheck += occ
            
        for sym in range(nsym):
            print >> lcix, "%2d %2d %2d" % (ig+1, len(gp), len(gp)),
            for ii in range(len(gp)):
                for jj in range(len(gp)):
                    print >> lcix, "%9.5f " % aocc[sym][ii,jj],
            print >> lcix
        
    
    print >> lcix, '# Data for HB1'
    print >> lcix, 1, len(gpseudo), len(Vmix), max(dim)
    PrintHeader(lcix, True, gpseudo, FM, ENE, STS, Vmix, eps, nbaths, dim, mN, mSz, mS2)
    print >> lcix, '# matrix elements'
    PrintMatrixElements(lcix,gpseudo,FM,Nb0)
    


if __name__ == '__main__':
    infile='actqmc.cix' 
    Eimp_new=[0,-0.0133061900015,-0.0664306100007]
    FromTransformation=True
    
    T0_str="""
     0.92421967  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.03141579  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.37731110  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.04967272  0.00000000    0.00000000  0.00000000  
     0.00000000  0.00000000    0.84369327  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.02221431  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.53359847  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.05441373  0.00000000  
    -0.00000000 -0.00000000   -0.00000000 -0.00000000    0.75592895 -0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000   -0.65465367 -0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000  
     0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.65465367  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.75592895  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000  
     0.05441374  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.00000000 -0.00000000    0.53359847 -0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000   -0.02221432 -0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000   -0.84369327 -0.00000000   -0.00000000 -0.00000000  
     0.00000000  0.00000000   -0.04967271  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.37731110  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.03141578  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.92421967  0.00000000  
     0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    1.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000  
     0.37796447  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.92582010  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000  
     0.00000000  0.00000000    0.53452248  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.84515425  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000  
     0.00000000  0.00000000    0.00000000  0.00000000    0.65465367  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.75592895  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000  
     0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.75592895  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.65465367  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000  
     0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.84515425  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.53452248  0.00000000    0.00000000  0.00000000  
     0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.92582010  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.37796447  0.00000000  
     0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    1.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000  
    """
    Tn_str="""
     0.88052456  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.16514575  0.00000000   -0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.35947264  0.00000000   -0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.26111836  0.00000000    0.00000000  0.00000000  
     0.00000000  0.00000000    0.80380526  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.11677569  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.00000000  0.00000000   -0.50837108  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.00000000  0.00000000   -0.28604086  0.00000000  
     0.00000000  0.00000000    0.00000000 -0.00000000    0.75592895 -0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.65465367  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000
     0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.65465367  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.75592895  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000  
     0.28604083  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.50837108  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.11677568  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.80380527  0.00000000    0.00000000  0.00000000  
     0.00000000  0.00000000   -0.26111838  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.35947264  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.16514576  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.88052455  0.00000000  
     0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    1.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000
     0.37796447  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.92582010  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000
     0.00000000  0.00000000    0.53452248  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.84515425  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000
     0.00000000  0.00000000    0.00000000  0.00000000    0.65465367  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.75592895  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000
     0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.75592895  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.65465367  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000
     0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.84515425  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.53452248  0.00000000    0.00000000  0.00000000
     0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.92582010  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.37796447  0.00000000
     0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    1.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000
    """
    Tn_str="""
     0.87982170  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.16638974  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.35918571  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.26308527  0.00000000    0.00000000  0.00000000
     0.00000000  0.00000000    0.80316365  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.11765531  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.50796529  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.28819550  0.00000000
    -0.00000000 -0.00000000   -0.00000000 -0.00000000    0.75592895 -0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000   -0.65465367 -0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000
     0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.65465367  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.75592895  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000
     0.28819547  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.50796529  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.11765531  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.80316366  0.00000000    0.00000000  0.00000000
     0.00000000  0.00000000   -0.26308529  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.35918570  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.16638975  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000   -0.87982170  0.00000000
     0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    1.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000
     0.37796447  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.92582010  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000
     0.00000000  0.00000000    0.53452248  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.84515425  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000
     0.00000000  0.00000000    0.00000000  0.00000000    0.65465367  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.75592895  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000
     0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.75592895  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.65465367  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000
     0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.84515425  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.53452248  0.00000000    0.00000000  0.00000000
     0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.92582010  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.37796447  0.00000000
     0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    1.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000
    """

    
    Delta_str="""
        0.20711067     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000      -0.10779232    -0.00000000      -0.00000000     0.00000000   
        0.00000000    -0.00000000      -0.17777680    -0.00000000      -0.00000000     0.00000000      -0.00000000    -0.00000000       0.00000000    -0.00000000      -0.10779232     0.00000000   
        0.00000000     0.00000000      -0.00000000    -0.00000000       0.29490355     0.00000000       0.00000000     0.00000000      -0.00000000     0.00000000       0.00000000    -0.00000000   
        0.00000000    -0.00000000      -0.00000000    -0.00000000      -0.00000000    -0.00000000       0.29490355     0.00000000      -0.00000000    -0.00000000       0.00000000     0.00000000   
       -0.10779232     0.00000000       0.00000000     0.00000000      -0.00000000    -0.00000000      -0.00000000     0.00000000      -0.17777683    -0.00000000      -0.00000000    -0.00000000   
       -0.00000000    -0.00000000      -0.10779232    -0.00000000       0.00000000     0.00000000       0.00000000    -0.00000000      -0.00000000     0.00000000       0.20711063    -0.00000000   
    """

    if FromTransformation:
        T0 = StringToMatrix(T0_str)
        Tn = StringToMatrix(Tn_str)
        Olap = dot( conjugate(T0), transpose(Tn) )
        Vmix0 = FindMixing0(Olap)
        Vmix = Vmix0[:6]
    else:
        Delta = StringToMatrix(Delta_str)
        Vmix = FindMixing(Delta)

    
    print 'Vmix=', Vmix
    
    #################################
    # Reading the one-site cix file #
    #################################
    (Nc0, Ns0, Nb0, Nm0, baths, eps, mN, mSz, msize, mEne, mS2, Fi, Fm) = ReadOneSiteCix(infile)
    
    for i in range(len(mEne)):
        print i+1,
        for im in range(msize[i]):
            print mEne[i,im],
        print

    mEne = CorrectSingleParticleEnergies(Ns0,Nb0,msize,mEne,baths,eps,Fi,Fm)

    for i in range(len(mEne)):
        print i+1,
        for im in range(msize[i]):
            print mEne[i,im],
        print

    Eimp = [Eimp_new[baths[ib][1]] for ib in range(Nb0)] 

    (STS,ENE,FM,gpseudo) = CreateNewStates(Nc0, Ns0, Nb0, Nm0, baths, eps, mN, mSz, msize, mEne, mS2, Fi, Fm, Eimp)
    
    print '######## Final print #######'

    print 'len(ENE)=', len(ENE)
    for ig,gp in enumerate(gpseudo):
        for l,g in enumerate(gp):
            print ig, l, ENE[g], STS[g]

    print 'F^+ operator'

    for ig in range(len(gpseudo)):
        for ib in range(Nb0):
            for ign in range(len(gpseudo)):
                bg = (ib,ig,ign)
                if FM.has_key( bg ):
                    print '(ib,igo,ign)=',bg
                    # dim (dim(gold),dim(gnew))
                    PrintM(FM[bg],1)
        

    PrintCixFile('test.cix', STS,ENE,FM,gpseudo, Eimp, mN, mSz,mS2,Vmix,baths)
