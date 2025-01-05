#!/usr/bin/env python
from numpy import linalg
from pylab import *
from tbg_diophantine import solve, lllhermite
from sympy import Matrix
from itertools import chain, product
from timeit import default_timer as timer
from numba import njit, jit

def tn_sin_cos(itheta):
    " angles for commensurate rotations"
    tn = 1./(2*sqrt(3.))/(itheta + 0.5) # tan(theta/2)
    sn = tn/sqrt(1.+tn*tn)              # sin(theta/2)
    cn = 1./sqrt(1.+tn*tn)              # cos(theta/2)
    return [sn,cn]

def RM(t):
    " 2D rotation "
    return array([[cos(t), -sin(t)], [sin(t), cos(t)]])


class TwistGraphene:
    def __init__(self, itheta, w0_w1=1.0):  # w1_w0 = w1/w0 is ratio between the AA and AB inter-layer-hopping (see Vishwanath)
        self.t  = 2730.                     # hooping in meV (tuned so that alpha=w/(v0*k_theta(theta=1.05))=0.606 (see Vishwanath)
        # hopping between the two layers is modeled by tq = t0 * exp(-(al*d_i)^gam * q^gam) with t0~2eV/A, al=0.13, gam=1.25 d_i=3.34 (see A. McDonnald PRB)
        # self.tq0 = t0/Vcell~1140meV, but to get w0 == t_q(q=kD)*Vcell = 110meV (see A. McDonalld), we slightly modify self.tq0=1066meV
        self.tq1 = 1066.
        self.tq0 = self.tq1*w0_w1
        self.al_di_gm = 0.379               # (al*d_i)^gam
        # 
        (sn,cn) = tn_sin_cos(itheta)        # [sin(theta/2), cos(theta/2)]
        self.theta = arcsin(sn)*2           # rotation angle
        print('theta=', self.theta*180/pi, 'sin(theta/2)=', sn)
        self.mR = RM(-self.theta/2.)        # rotation of the basis counter-clockwise for theta/2
        self.pR = RM( self.theta/2.)        # rotation of the basis clockwise for theta/2
        self.mM = self.pR                   # this is rotation for vectors. When basis is rotated counterclockwise, the vectors must be rotated clockwise
        self.pM = self.mR                   # i.e., active versus passive rotation.
        # Choice for the cartesian coordinates follows A. McDonnald, i.e., one graphene is rotated for +theta/2 and the other for -theta/2 compared to the cartesian coordinates
        # self.b0 is reciprocal unit vector for the fake unrotated unit cell (between the two graphene sheets)
        self.b0 =   [array([ 1.0,1.0/sqrt(3.)])*2*pi, array([-1.0,1.0/sqrt(3.)])*2*pi] # reciprocal unit vectors prior to rotation
        self.R0 = linalg.inv(transpose(self.b0))*(2*pi)  # direct vectors
        # reciprocal unit vector for small emergent unit cell at comensurate angle. It is rotated for 30 degrees compared to fake large BZ (convince yourself that this is indeed correct).
        _bs_ = [array([2.0/sqrt(3.),0])*2*pi, array([-1.0/sqrt(3.),1.0])*2*pi] 
        # reciprocal and direct unit vectors for the two sheets
        self.bm = [dot( self.mR, self.b0[0] ), dot( self.mR, self.b0[1] ) ] # left sheet reciprocal : rotated counter-clocwise for theta/2
        self.Rm = linalg.inv(transpose(self.bm))*(2*pi)                     # left sheet direct vectors
        self.bp = [dot( self.pR, self.b0[0] ), dot( self.pR, self.b0[1] ) ] # right sheet reciprocal : rotated clockwise for theta/2
        self.Rp = linalg.inv(transpose(self.bp))*(2*pi)                     # right sheet direct vectors
        # emergent unit cell with proper angle added
        self.bs = [_bs_[0]*2*sn, _bs_[1]*2*sn]                              # reciprocal vectors of the emergent unit cell
        self.Rs = linalg.inv(transpose(self.bs))*(2*pi)                     # real space vectors of the emergent unit cell
        # nearest neighbors on each graphene sheet, properly taking into account rotation of both sheets
        self.Rnnm = array([(self.Rm[0] + self.Rm[1])/3., (self.Rm[1]-2*self.Rm[0])/3., (self.Rm[0]-2*self.Rm[1])/3.]) # nearest neighbors on left sheet
        self.Rnnp = array([(self.Rp[0] + self.Rp[1])/3., (self.Rp[1]-2*self.Rp[0])/3., (self.Rp[0]-2*self.Rp[1])/3.]) # nearest neighbors on right sheet
        # positions of special points K==K1 and K'==K2 on the unrotated (K1, & K2) and the two rotated layers (K1m & K2m, K1p & K2p)
        self.K0 = zeros((3,2))
        self.K0[0,:] =  (2*self.b0[1]+self.b0[0])/3.
        self.K0[1,:] =  (2*self.b0[0]+self.b0[1])/3.
        self.K0[2,:] =  (self.K0[1]-self.b0[0])
        self.K1m = (2*self.bm[1]+self.bm[0])/3.
        self.K2m = (2*self.bm[0]+self.bm[1])/3.
        self.K1p = (2*self.bp[1]+self.bp[0])/3.
        self.K2p = (2*self.bp[0]+self.bp[1])/3.

        self.rBAm = (self.Rm[0]+self.Rm[1])/3. # BA vector on minus-rotated layer
        self.rBAp = (self.Rp[0]+self.Rp[1])/3. # BA vector on plus-rotated layer
        self.r_BAm = dot(self.rBAm, transpose(self.bs))/(2*pi) # BA vector on minus layer, but in terms of Rs
        self.r_BAp = dot(self.rBAp, transpose(self.bs))/(2*pi) # BA vector on plus layer, expressed in terms of Rs
        # note
        # rBAm = self.r_BAm[0]*self.Rs[0]+self.r_BAm[1]*self.Rs[1] == dot(self.r_BAm,self.Rs)
        # rBAp = self.r_BAp[0]*self.Rs[0]+self.r_BAp[1]*self.Rs[1] == dot(self.r_BAp,self.Rs)
        
    def tq_w0(self, q):
        "AA-AA hopping"
        return self.tq0 * exp(-self.al_di_gm * q**1.25)
    def tq_w1(self, q):
        "AA-AB hopping"
        return self.tq1 * exp(-self.al_di_gm * q**1.25)
    
    def zkm(self, k):
        return self.t * sum(exp(dot(self.Rnnm,k)*1j))
    def zkp(self, k):
        return self.t * sum(exp(dot(self.Rnnp,k)*1j))
        
    def H0m(self, k):
        "H0(-theta/2,k)"
        z_tk = self.t * sum(exp(dot(self.Rnnm,k)*1j))
        return array([[0.,z_tk],[conj(z_tk),0.]])
    def H0p(self, k):
        "H0(+theta/2,k)"
        z_tk = self.t * sum(exp(dot(self.Rnnp,k)*1j))
        return array([[0.,z_tk],[conj(z_tk),0.]])

    def Vkp(self, k, G1, p, G2):
        k_Gk = dot(k + G1, self.bs)
        tw0 = self.tq_w0(norm(k_Gk))
        tw1 = self.tq_w1(norm(k_Gk))
        phase_m = dot(G1, self.r_BAm)
        phase_p = dot(G2, self.r_BAp)
        exp_m = exp(2*pi*1j*phase_m)
        exp_p = exp(2*pi*1j*phase_p)
        return array([[tw0, tw1*conj(exp_p)],[tw1*exp_m, tw0*exp_m*conj(exp_p)]])
    
    def n0_max(self):
        return int(1/(2*sin(self.theta/2))+1)

    def get_Nm_Np(self):
        """ reciprocal vectors of both sheets (bm, bp) can be expanded in terms of reciprocal vectors of the emergent BZ (bz). 
        At the commensurate angle, the expansion has to have integer coefficients. This routine gives such integer coefficients
        """
        Nm = zeros((2,2),dtype=int)
        Np = zeros((2,2),dtype=int)
        for i in range(2):
            for j in range(2):
                nm = dot(self.bm[i],self.Rs[j])/(2*pi)
                rnm = round(nm)
                if sum(abs(rnm-nm))>1e-6:
                    print('ERROR : It seems the angle is not commensurate and reciprocal vectors of the rotated layer are not included in the emergent BZ')
                    print('nm=', nm)
                Nm[i,j] = round(nm)
                np = dot(self.bp[i],self.Rs[j])/(2*pi)
                rnp = round(np)
                if sum(abs(rnp-np))>1e-6:
                    print('ERROR : It seems the angle is not commensurate and reciprocal vectors of the rotated layer are not included in the emergent BZ')
                    print('np=', np)
                Np[i,j] = rnp
        return (Nm, Np)
    
    def Find_k_p_vcts(self, cutoffk = 10000.):
        """ Given k-point in the emergent BZ, there are many momentum points in the first BZ of the two graphene sheets that 
        couple with the given k-point. Here we produce sorted list of momenta for both sheets (n_k minus rotation, n_p plus rotation)
        which couple (are equivalent) to Gamma point in the emergent BZ.
        """
        n0 = self.n0_max()
        # First produce decomposition of the reciprocal vector of the emergent unit cell in terms of reciprocal vectors of the real graphene lattice
        bsRm = dot(self.bs, transpose(self.Rm))/(2*pi)  # bs_i = \sum_j bsRm_{ij} bm_j
        bsRp = dot(self.bs, transpose(self.Rp))/(2*pi)  # bs_i = \sum_j bsRp_{ij} bp_j

        n_k=[]
        n_p=[]
        kd = norm(self.K0[0])
        normb = norm(self.b0[0])
        for n1 in range(-2*n0,2*n0):
            for n2 in range(-2*n0,2*n0):
                # Use this decomposition to quickly find all momenta (in integer form), which are equivalent to a given momentum in the emergent BZ.
                # nkm is for sheet rotated in negative direction
                nkm = dot([n1,n2],bsRm)  # \sum_i n_i bs_i = \sum_j (n * bsRm)_j bm_j
                # nkp is for sheet rotated in positive direction
                nkp = dot([n1,n2],bsRp)  # \sum_i n_i bs_i = \sum_j (n * bsRp)_j bp_j
                # Now check if each of the two vectors is in their respective BZ. In integer form, the vector just needs to have all components in [0,1] interval.
                if IsInBZ(nkm):
                    n_k.append((n1,n2))
                if IsInBZ(nkp):
                    n_p.append((n1,n2))
        # Now we will sort all these k-points according to parameter t*zk, which is proportinal to energy
        Ek=[]
        for nk in n_k:
            kc = dot(nk,self.bs)
            Ek.append( abs(self.zkm(kc)) )
        # and also compute such energies for the other sheet
        Ep=[]
        for np in n_p:
            pc = dot(np,self.bs)
            Ep.append( abs(self.zkp(pc)) )
        # Now sort such k-points by energy
        indk = sorted(list(range(len(Ek))),key= lambda i: Ek[i])
        indp = sorted(list(range(len(Ep))),key= lambda i: Ep[i])
        # and produce new list, which is sorted, and includes only energies below the cutoff.
        n_k_2=[]
        n_p_2=[]
        for ii in range(len(indk)):
            ik = indk[ii]
            ip = indp[ii]
            if Ek[ik] < cutoffk:
                n_k_2.append( array(n_k[ik], dtype=int) )
            if Ep[ip] < cutoffk:
                n_p_2.append( array(n_p[ip], dtype=int) )
            #print ii, Ek[ik], Ep[ip], n_k[ik], n_p[ip], norm(dot(n_k[ik],self.bs))/kd, norm(dot(n_p[ip],self.bs))/kd
            
        if (len(n_k_2)<100): # for small number of k-points, we want to optimize further the mesh
            # We want to make sure that every point that appears in k-mesh also appears in p-mesh, if possible
            nks = [tuple(k) for k in n_k_2]
            nps = [tuple(k) for k in n_p_2]
            nk_extra = []
            np_extra = []
            for k in nks:
                if (k not in nps) and (k in n_p):
                    np_extra.append(k)
            for p in nps:
                if (p not in nks) and (k in n_k):
                    nk_extra.append(p)
                    
            n_k_2 += [array(k,dtype=int) for k in nk_extra]
            n_p_2 += [array(p,dtype=int) for p in np_extra]
        return (array(n_k_2, dtype=int), array(n_p_2, dtype=int))  # notice of change
    
    def Find_G1_G2_vectors(self, n_pk, Nm, Np, how_far=[-1,0,1]):
        """Finding possible reciprocal vectors G1,G2 which connect arbitrary k and p momentum points.
        input: n_pk is a dictionary of p-k points expressed in emergent BZ.
        
           We are solving the equation 

           Gm + Gp = p-k

           where Gm and Gp are arbitrary reciprocal vectors for the counter-clockwise rotated (Gm) and clockwise-rotated (Gp) sheet, 
           and p,k correspond to the same vector in the emergent BZ, but extent beyond first BZ. 
           Consequently, p-k must be a reciprocal vector of the emergent BZ. Likewise, Gm and Gp can also be expressed as reciprocal vectors of 
           the emergent BZ, hence we can express the equation in the basis of the emergent BZ: Gm+Gp-(p-k)=0 => (n_Gm+n_Gp-(n_p-n_k))*(bs)=0,
           where all entries have to be integers.
           
           The equation (in terms of the emergent BZ recoprocal vectors) takes the form:
           
           (n_Gm,n_Gp) * (Nm,Np) = n_p-n_k

           where (n_Gm,n_Gp) is a list of 4 unknown integers, (Nm,Np) is a 4x2 matrix of integers, and n_p-n_k is 2-component integer vector.
           
           This is so-called "system of linear DIOPHANTINE equations", i.e., see:
              http://sites.math.rutgers.edu/~sk1233/courses/ANT-F14/lec3.pdf
              http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.94.6106&rep=rep1&type=pdf

           and can be solved by transforming matrix (Nm,Np) into Hermite-Normal-Form. The linear combination of the rows of the pivotal matrix P 
           than form a solution of the problem.
        """
        # Construct a sympy matrix, which can be used to compute Hermite-Normal-Form of the same matrix
        Gn = zeros((4,2),dtype=int)
        Gn[:2,:] = Nm
        Gn[2:,:] = Np
        Gn = Matrix(Gn)
        # Now we call routine to compute Hermite-Normal-Form of matrix [Nm,Np], and get pivotal matrix P. The rows of this matrix are solutions.
        hnf, P, rank = lllhermite(Gn)
        Pm = array(P).astype(np.int64) # transforming it back to faster numpy form

        if hnf[0,1]!=0 or hnf[0,0]!=1 or hnf[1,1]!=1:
            print('ERROR for your angle, Hermite-Normal-Form is not diagonal. You will need to implement more complete Diophantine algorithm')
            sys.exit(0)
        
        for npk in list(n_pk.keys()):
            GpG = npk[0]*Pm[0,:] + npk[1]*Pm[1,:] # one solution
            # The last two rows can be added with an arbitrary integer coefficient. Trying to add [-1,0,1] of each row.
            other_forms = array([GpG + Pm[2,:]*i + Pm[3,:]*j for (i,j) in product(how_far,how_far)])
            # now checking that these values indeed solve the equation
            check = [ sum(abs(dot(gpg[:2],Nm) + dot(gpg[2:],Np) - npk)) for gpg in other_forms]
            if (sum(check)!=0):
                print('ERROR', 'npk=', npk, 'but the difference G1-G2=', [ dot(gpg[:2],Nm) + dot(gpg[2:],Np) for gpg in other_forms])
                sys.exit(1)
            # Finally, saving all these forms into dictionary.
            n_pk[npk] = other_forms 
            #nrm = [norm(gpg) for gpg in other_forms]
            #print 'npk=', npk, 'c0=', GpG, 'other_forms=', other_forms, 'nrm=', nrm
        return n_pk
    
    def Find_If_Out_of_BZ(self, ks_lattice, n_k, Nm, debug=True):
        """ when k+k_s drops out of the first BZ, we bring it back in with n_k_reciprocal
            note that everything is here expressed in terms of emergent BZ vectors self.bs
        """
        small=1e-6
        Nm_inv = linalg.inv(Nm)
        n_k_reciprocal = zeros((len(n_k),2),dtype=int) 
        for ik in range(len(n_k)):
            nk = ks_lattice + n_k[ik]
            km = dot(nk,Nm_inv) # k expressed in terms of reciprocal vectors bm.
            dk = zeros(2, dtype=int)
            for l in range(2):
                if km[l]>1+small:
                    dk[l]=-1
                if km[l]<-small:
                    dk[l]+=1
            if sum(abs(dk))>0:
                n_k_reciprocal[ik,:] = dot(dk,Nm) # now trasform from reciprocal bm expansion to reciprocal bs expansion.
                if (debug):
                    nk_new = nk + n_k_reciprocal[ik]
                    km_new = dot(nk_new,Nm_inv) # k expressed in terms of reciprocal vectors bm.
                    print(ik, nk, km, n_k_reciprocal[ik], km_new)
        return n_k_reciprocal

    def Numba_Optimized_npk(self, n_pk):
        """Since numba can not deal with dictionaries, we will create index array and list, which combined work like dictionary.
        so, instead on n_pk[key], we have self.n_pk_index and self.n_pk_value, so that
           n_pk[key] == self.n_pk_value[self.n_pk_index[key[0],key[1]]]
        """
        _keys_ = n_pk.keys()
        k0=sorted(set([k[0] for k in _keys_]))
        k1=sorted(set([k[1] for k in _keys_]))
        self.n_pk_index = zeros((k0[-1]-k0[0]+1,k1[-1]-k1[0]+1),dtype=int)
        self.n_pk_value = [[] for ii in range(len(n_pk))]
        for ii,_k_ in enumerate(n_pk):
            self.n_pk_index[_k_[0],_k_[1]]=ii
            self.n_pk_value[ii] = n_pk[_k_]
        #for _k_ in n_pk:
        #    aclose = allclose(n_pk[_k_], n_pk_value[n_pk_index[_k_[0],_k_[1]]])
        #    print(_k_, '  ', aclose, shape(n_pk[_k_]))
    
    def BuildIterlayerH(self, Hk, ks_lattice, n_pk, n_k, n_p, n_k_reciprocal, n_p_reciprocal, Nm, Np):
        FAST=True
        dim_k = len(n_k)
        if (not FAST): # readable code
            for ik,ip in product(range(len(n_k)),range(len(n_p))):
                npk = n_p[ip]-n_k[ik]
                GpG = n_pk[tuple(npk)]
                k = ks_lattice + n_k[ik] + n_k_reciprocal[ik] # k in interger representation
                p = ks_lattice + n_p[ip] + n_p_reciprocal[ip]
                for gpg in GpG:
                    n_Gk = dot(gpg[:2],Nm)       # Gk in integer representation
                    # we have:
                    # np - nk = Gp + Gk
                    #    but we want to have
                    # k + G1 == p + G2
                    #    hence
                    # G1 =  Gk - n_k_reciprocal[ik]
                    # G2 = -Gp - n_p_reciprocal[ip]
                    G1 = n_Gk - n_k_reciprocal[ik]
                    k_Gk = dot(k + G1, self.bs) # k+Gk in cartesian coordinates
                    tw0  = self.tq_w0(norm(k_Gk))  # AA or BB hopping in meV
                    if abs(tw0) > 1e-4:
                        n_Gp = dot(gpg[2:],Np)          # Gp in integer representation
                        G2 = -n_Gp - n_p_reciprocal[ip]
                        
                        phase_m = dot(G1, tw.r_BAm)     # e^(i*r^-_{BA}*Gk)
                        phase_p = dot(G2, tw.r_BAp)     # e^(i*r^+_{BA}*Gp)
                        
                        Vkp = tw.Vkp(k, G1, p, G2) # off-diagonal H
                        
                        Hk[2*ik:2*(ik+1),2*(dim_k+ip):2*(dim_k+ip+1)] += Vkp
                        Hk[2*(dim_k+ip):2*(dim_k+ip+1),2*ik:2*(ik+1)] += Vkp.conj().T
                        
                        if (False):
                            str_Gk = "Gk=[%4d,%4d]" % tuple(G1)
                            str_Gp = "Gp=[%4d,%4d]" % tuple(G2)
                            str_kGk = 'k+Gk=[%6.2f,%6.2f]' % tuple(k+G1)
                            str_pGp = 'p+Gp=[%6.2f,%6.2f]' % tuple(p+G2)
                            str_phase = "phase=%6.3f/3, %6.3f/3" % (phase_m*3,phase_p*3)
                            str_w0 = 'tw0=%-6.2f' % (tw0,)
                            print('nk=', n_k[ik], 'np=', n_p[ip], str_Gk, str_Gp, str_kGk, str_pGp, str_w0, str_phase, 'Vkp=', Vkp.tolist())
        else:  # maximally optimized code, but equivalent to above
            self_bs = array(self.bs, dtype=float64)
            self_r_BAm = array(self.r_BAm, dtype=float64)
            BuildIterlayerH_faster(Hk, ks_lattice, self.n_pk_index, self.n_pk_value,
                                       n_k, n_p, n_k_reciprocal, n_p_reciprocal, Nm, Np,
                                       self_bs, self.tq0, self.tq1, self.al_di_gm, self_r_BAm, self.r_BAp)
    
def IsInBZ(nk):
    small = 1e-7
    return (nk[0] > -small) and (nk[0] < 1-small) and (nk[1] > -small) and (nk[1] < 1-small)


def FindUnique_p_m_k(n_k, n_p):
    # Find how many unique p-k combinations we can have
    n_pk={}
    for ik,ip in product(range(len(n_k)),range(len(n_p))):
        npk = tuple( n_p[ip]-n_k[ik] )
        if npk not in n_pk:
            n_pk[npk]=1
    return n_pk

def PrintH(A):
    print('Real=')
    for i in range(A.shape[0]):
        print("%6.1f "*(A.shape[1]) % tuple(A[i,:].real))
    print('Imag=')
    for i in range(A.shape[0]):
        print("%6.1f "*(A.shape[1]) % tuple(A[i,:].imag))

def Give_k_Path(path, Nkp, bs):
    lng=[]
    for l in range(len(path)-1):
        k0 = array(path[l])
        k1 = array(path[l+1])
        lng.append( norm(dot(k1-k0,bs)) )
    lng = array(lng)*(Nkp-1)/sum(lng)
    Npth = array(round_(lng),dtype=int)
    kpath=[]
    for l in range(len(path)-1):
        k0 = array(path[l])
        k1 = array(path[l+1])
        nn = Npth[l]
        if (l==len(path)-2): nn += 1
        for i in range(nn):
            k = k0 + (k1-k0)*i/(Npth[l]-0.0)
            kpath.append(k)
    kpath = array(kpath)
    return (kpath, Npth)

@jit(nopython=True)
def BuildIterlayerH_faster(Hk, ks_lattice, n_pk_index, n_pk_value,
                            n_k, n_p, n_k_reciprocal, n_p_reciprocal, Nm, Np,
                            self_bs, self_tq0, self_tq1, self_al_di_gm, self_r_BAm, self_r_BAp):
    for ik in range(len(n_k)):
        for ip in range(len(n_p)):
            npk = n_p[ip]-n_k[ik]
            #GpG = n_pk[tuple(npk)]
            GpG = n_pk_value[n_pk_index[npk[0],npk[1]]]
            k = ks_lattice + n_k[ik] + n_k_reciprocal[ik] # k in interger representation
            p = ks_lattice + n_p[ip] + n_p_reciprocal[ip]
            for gpg in GpG:
                n_Gk = gpg[0]*Nm[0]+gpg[1]*Nm[1]   # Gk in integer representation
                # we have:
                # np - nk = Gp + Gk
                #    but we want to have
                # k + G1 == p + G2
                #    hence
                # G1 =  Gk - n_k_reciprocal[ik]
                # G2 = -Gp - n_p_reciprocal[ip]
                G1 = n_Gk - n_k_reciprocal[ik]
                k_Gk = (k + G1) @ self_bs    # k+Gk in cartesian coordinates
                t_exp = exp(-self_al_di_gm * norm(k_Gk)**1.25)
                tw0, tw1 = self_tq0 * t_exp,  self_tq1 * t_exp
                if abs(tw0) > 1e-4:
                    phase_m = (G1 + 0.0) @ self_r_BAm
                    n_Gp = gpg[2]*Np[0]+gpg[3]*Np[1]       # Gp in integer representation
                    G2 = -n_Gp - n_p_reciprocal[ip]
                    phase_p = (G2 + 0.0) @ self_r_BAp
                    exp_m = exp(2*pi*1j*phase_m)
                    exp_p = exp(2*pi*1j*phase_p)
                    Vkp = array([[tw0, tw1*conj(exp_p)],[tw1*exp_m, tw0*exp_m*conj(exp_p)]])
                    
                    Hk[2*ik:2*(ik+1),2*(dim_k+ip):2*(dim_k+ip+1)] += Vkp
                    Hk[2*(dim_k+ip):2*(dim_k+ip+1),2*ik:2*(ik+1)] += Vkp.conj().T
            
                    #if (False):
                    #    #phase_m = dot(G1, self_r_BAm)     # e^(i*r^-_{BA}*Gk)
                    #    #phase_p = dot(G2, self_r_BAp)     # e^(i*r^+_{BA}*Gp)
                    #    str_Gk = "Gk=[%4d,%4d]" % tuple(G1)
                    #    str_Gp = "Gp=[%4d,%4d]" % tuple(G2)
                    #    str_kGk = 'k+Gk=[%6.2f,%6.2f]' % tuple(k+G1)
                    #    str_pGp = 'p+Gp=[%6.2f,%6.2f]' % tuple(p+G2)
                    #    str_phase = "phase=%6.3f/3, %6.3f/3" % (phase_m*3,phase_p*3)
                    #    str_w0 = 'tw0=%-6.2f' % (tw0,)
                    #    print('nk=', n_k[ik], 'np=', n_p[ip], str_Gk, str_Gp, str_kGk, str_pGp, str_w0, str_phase, 'Vkp=', Vkp.tolist())



if __name__ == '__main__':
    # This allows user to change parameters by executing params.py
    exec(open("params.py").read())
    # this uses only par_bs
    for p in list(par_bs.keys()):
        ss = p + '=' + str(par_bs[p])
        print(ss)
        exec(ss)
    
    #shift_mesh = 1./3.
    tw = TwistGraphene(itheta, w0_w1)

    if PlotBands:
        kpath, Npth = Give_k_Path(path, Nkpt, tw.bs)
        #print 'kpath=', kpath
    else:
        N1 = int(sqrt(Nkpt))
        kpath=[]
        if shift_mesh != 0.0:
            for i,j in product(range(N1),range(N1)):
                kpath.append( [ (i+shift_mesh)/(N1+0.0), (j+shift_mesh)/(N1+0.0)] )
        else:
            for i,j in product(range(N1),range(N1)):
                kpath.append( [i/(N1+0.0), j/(N1+0.0)] )
        kpath = array(kpath)
        
        
    
    (Nm, Np) = tw.get_Nm_Np()
    #print('Nm=', Nm)
    #print('Np=', Np)
            
    _t1_ = timer()
    n_k, n_p = tw.Find_k_p_vcts(Ecutoff)
    _t2_ = timer()
    print('#nk=', len(n_k), '#np=', len(n_p), 't(generate-k)=', _t2_-_t1_)

    print('n_k=', shape(n_k))
    print('n_p=', shape(n_p))

    dim_k = len(n_k)
    dim_p = len(n_p)

    n_pk = FindUnique_p_m_k(n_k, n_p)
    _t3_ = timer()
    print('#n_pk=', len(n_pk), 't(find-unique)=', _t3_-_t2_)
    
    n_pk = tw.Find_G1_G2_vectors(n_pk, Nm, Np)

    tw.Numba_Optimized_npk(n_pk)
        
    nh = (dim_k+dim_p)*2 # size of Hamiltonian
    imin = max(0,int(nh/2)-npr)
    imax = min(  int(nh/2)+npr,nh)

    save('n_k', n_k)
    save('n_p', n_p)
    
    Eks=[]
    Psis=[]
    k_basis=[]
    for iik in range(len(kpath)):
        
        ks_lattice = kpath[iik]
    
        # when k+k_s drops out of the first BZ, we bring it back in with k_reciprocal
        # note that everything is here expressed in terms of emergent BZ vectors tw.bs
        n_k_reciprocal = tw.Find_If_Out_of_BZ(ks_lattice, n_k, Nm, False)
        n_p_reciprocal = tw.Find_If_Out_of_BZ(ks_lattice, n_p, Np, False)
        
        Hk = zeros(( (dim_k+dim_p)*2, (dim_k+dim_p)*2 ), dtype=complex)
        
        print(iik, 'H-size=', len(Hk))

        _t5_ = timer()
        
        k_bas=[]
        for ik in range(len(n_k)):
            nk = ks_lattice + n_k[ik] + n_k_reciprocal[ik]
            kc = dot(nk, tw.bs)  # k in cartesian
            Hk[2*ik:2*(ik+1),2*ik:2*(ik+1)] += tw.H0m(kc)
            k_bas.append(kc)
        for ip in range(len(n_p)):
            np = ks_lattice + n_p[ip] + n_p_reciprocal[ip]
            pc = dot(np, tw.bs) # p in cartesian
            Hk[2*(dim_k+ip):2*(dim_k+ip+1),2*(dim_k+ip):2*(dim_k+ip+1)] += tw.H0p(pc)
            k_bas.append(pc)
        k_basis.append(k_bas)
        _t6_ = timer()

        print(iik, 't(H-diagonal)=', _t6_-_t5_)
        
        #print 'Hk0='
        #PrintH(Hk)
        
        tw.BuildIterlayerH(Hk, ks_lattice, n_pk, n_k, n_p, n_k_reciprocal, n_p_reciprocal, Nm, Np)
        
        _t7_ = timer()
        print(iik, 't(H-offdiagonal)=', _t7_-_t6_)
        
        diff = Hk - Hk.conj().T
        print(iik, 'Is Hermitian=', (sum(abs(diff)) < 0.01))
        #print 'Hk='
        #PrintH(Hk)

        if False: #PlotBands:
            ek = linalg.eigvalsh(Hk)
        else:
            ek,psik = linalg.eigh(Hk)
            psi = psik.T
            Psis.append( psi[imin:imax] )

        _t8_ = timer()
        print(iik, 't(diagonalize)=', _t8_-_t7_)
        
        Eks.append(  ek[imin:imax] )
        
        print(iik, 'ek=', ek[imin : imax])
        
        #print
        #print 'Hk='
        #PrintH(Hk)

    ne = len(Eks[0])
    e0m = max([Eks[iik][int(ne/2)-1] for iik in range(len(Eks))])
    e0p = min([Eks[iik][int(ne/2)]   for iik in range(len(Eks))])
    e0 = (e0m+e0p)/2.
    print('e0=', e0)
    
    k_basis = array(k_basis)
    Psis = array(Psis)
    Eks = array(Eks)
    
    if PlotBands:
        fout = open('res1.dat', 'w')
        for iik in range(len(kpath)):
            evs = (Eks[iik]-e0)/unit
            print(iik, "%12.7f "*(imax-imin) % tuple(evs), file=fout)
        fout.close()

        savez('kpath_band', kpath=kpath, Npth=Npth)
        save('ene_band', (Eks-e0)/unit)
        save('vector_band', Psis)
        savez('basis_band', k_basis=k_basis, Rm=tw.Rm, Rp=tw.Rp, Rs=tw.Rs, rBAm=tw.rBAm, rBAp=tw.rBAp, bm=tw.bm, bp=tw.bp, bs=tw.bs, b0=tw.b0, Nm=Nm, Np=Np)
        
    else:
        fout = open('ek.dat', 'w')
        for iik in range(len(kpath)):
            evs = (Eks[iik]-e0)/unit
            kp = kpath[iik]*N1
            print(int(kp[0]), int(kp[1]), "%12.7f "*(imax-imin) % tuple(evs), file=fout)
        fout.close()

        save('ene', (Eks-e0)/unit)
        save('vector', Psis)
        savez('basis', k_basis=k_basis, Rm=tw.Rm, Rp=tw.Rp, Rs=tw.Rs, rBAm=tw.rBAm, rBAp=tw.rBAp, bm=tw.bm, bp=tw.bp, bs=tw.bs, b0=tw.b0, Nm=Nm, Np=Np)
