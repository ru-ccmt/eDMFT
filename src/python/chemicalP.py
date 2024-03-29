# @Copyright 2007 Kristjan Haule  
#from scipy import *
import numpy as np
from numpy import pi
from scipy import interpolate
import os,sys
import occupi
import occupr


def FindSignChange(fComputeCharge, old_mu, sdmu, fh_info, args):
    """ Routine brackets the chemical potential.
        Input:
          old_mu                 -- start of looking for mu
          sdmu                   -- the interval will be looked in the steps of sdmu in the following algorithm:
                                      1) try [mu, mu+sdmu]
                                      2) try [mu+sdmu, mu+2*sdmu]
                                      3) try [mu+2*sdmu, mu+2*2*sdmu]
                                      ....
    """
    print('Looking for the chemical potential. Starting from old_mu=', old_mu, file=fh_info)
    
    curr_mu = old_mu

    curr_dens = fComputeCharge(*(curr_mu,)+args) 
    
    #print >> fh_info, '(mu,dens)=', curr_mu, curr_dens

    if abs(curr_dens)<1e-6:
        dtmu = 1e-4
        return (curr_mu-dtmu, curr_mu+dtmu)
    elif (curr_dens<0):
        tdmu = sdmu
    else:
        tdmu = -sdmu
    
    while (True):
        new_dens = fComputeCharge(*(curr_mu+tdmu,)+args) #Egns, lmt, om, fh_info, noccb, max_metropolis_steps, use_tetra, LowerBound)
        print('(mu,dens)=', curr_mu+tdmu, new_dens, file=fh_info)
        if curr_dens*new_dens<0: break
        curr_dens = new_dens
        curr_mu += tdmu
        tdmu *= 2.
        
    print('(mu0,mu1), (dens0,dens1)', curr_mu, curr_mu+tdmu, curr_dens, new_dens, file=fh_info)
    return (curr_mu, curr_mu+tdmu)

def ReadEigenvals(filename, Qimag, gamma):
    """ Reads the input file, created by fortran code, containing the LDA+DMFT eigenvalues.
        Input:
              filename
              gamma                    -- minimum broadening
        Output:
           Ek[nkp][nbands,nomega]      -- eigenvalues
           omega[nomega]               -- frequency
           wkp[nkp]                    -- k-point weights
           nbands[nkp]                 -- number of bands kept
           nemin[nkp]                  -- number of bands not treated dynamically, being fully filled
    """

    #print 'Qimag:=', Qimag
    
    # Reading all lines of the file
    with open(filename, 'r') as fw:
        adata = fw.readlines()
    
    ii=0
    Ek=[]
    wkp=[]
    nbands=[]
    nemin=[]
    while ii<len(adata):
        data = adata[ii].split()
        ii+=1
        if len(data)<6: break
        (ikp, isym, nbandsk, nemink, nomega) = list(map(int, tuple(data[1:6])))
        nbands.append(nbandsk)
        nemin.append(nemink)
        wk = float(data[6])
        
        omega = np.zeros((nomega))
        tEk = np.zeros((nbandsk,nomega), dtype=complex, order='F')
        for iom in range(nomega):
            data = list(map(float, adata[ii].split()))
            ii+=1
            omega[iom] = data[0]
            for ip in range(nbandsk):
                if data[2+2*ip]>-gamma: data[2+2*ip]=-gamma
                tEk[ip,iom] = data[1+2*ip] + data[2+2*ip]*1j
        Ek.append(tEk)
        wkp.append(wk)

    if Qimag !=0 :

        print('... Found matsubara and processing imaginary axis ...', Qimag)
        maxbands= max(nbands)
        
        Ekz = np.zeros((maxbands,len(Ek),np.shape(Ek[0])[1]), dtype=complex, order='F')
        for ikp in range(len(Ek)):
            Ekz[:nbands[ikp],ikp,:] = Ek[ikp][:,:]

        Ek = Ekz
        
    return (Ek, omega, wkp, nbands, nemin)
    

def ChargeDiff(EF, valC, NOE):
    charge = 0
    for i in range(len(valC)):
        ch = valC[i](EF)
        print('EF=', EF, 'Charge ', i, '=', ch)
        charge += ch
    return charge - NOE
    
class ValCharge:
    def __init__(self, Qimag, wkp, omega, Ek, nemin, nbands, com, L=-1000):
        self.Qimag = Qimag
        self.wkp = np.array(wkp)
        self.L = L

        self.omega = omega
        self.Ek = Ek
        self.nemin = nemin
        self.nbands = nbands
        
        if (self.Qimag): # Imaginary axis
            # Creating big mesh with all matsubara points!
            # The mesh omega has only few points on logarithmic mesh
            n = int(0.5*(omega[-1]/omega[0]-1) + 0.4999)
            T = omega[-1]/(2*n+1)/pi
            self.lom=np.zeros(n+1, dtype=float)
            for i in range(n+1): self.lom[i] = (2*i+1)*pi*T
            print('Large matsubara mesh of size', n+1, 'was created')

            # Creating the array Ek0 which is used to converge the matsubara sum.
            # We will computer the following sum:    1/(om+mu-Ek) - 1/(om+mu-Ek0)
            self.Ek0 = np.zeros((np.shape(Ek)[0],np.shape(Ek)[1]), dtype=float, order='F')
            for ikp in range(np.shape(Ek)[1]):
                for i in range(np.shape(Ek)[0]):
                    if com==0:
                        self.Ek0[i,ikp] = Ek[i,ikp,-1].real
                    else:
                        self.Ek0[i,ikp] = Ek[i,ikp,0].real
            
        

    def __call__(self, mu):
        """ Computes the valence charge by call to fortran function occupx
            The latter perform integral over frequency analytically.
            The sum over k-points is done by special-points method.
            The return value if the difference between the valence_charge(mu)
            Input:
               mu                   --  current mu
               Ek[nkp][nbands,nom]  -- eigenvalues on real axis
               Ek[maxbands,nkp,nom] -- eigenvalues on imaginary axis
               omega[nomega]        -- frequency
               nemin                -- the number of bands which are not treated dynamically and are thus fully filled
               wkp                  -- weight of each k-point
               L                    -- the lower cut-off where the integrals starts; should be negative enough
            Output:
               charge-NOE           -- the missing charge
        """
        if (self.Qimag): # Imaginary axis
            
            # charge1 is sum of the fermi functions
            # charge2 is correction due to finite number of matsubara points used in the sum
            (charge, dcharge1, dcharge2) = occupi.occupi(mu, self.Ek, self.Ek0, self.wkp, self.omega, self.nbands)
            # interpolate on big mesh
            tckr = interpolate.splrep(self.omega, charge, s=0)
            lcharge = interpolate.splev(self.lom, tckr)
            # sum over the big mesh to get density
            nt1 = sum(lcharge)
            # bands which are not included in DMFT calculation and are fully filled
            nts = sum(self.wkp * np.array(self.nemin)) - 1.  # -1 because skip is one less than the first included
            return nts + dcharge1 + dcharge2 + nt1
        
        else: # Real axis
            
            charge=0
            for ikp in range(len(self.Ek)):
                charge += (occupr.occupr(mu, self.Ek[ikp], self.omega, self.L) + self.nemin[ikp] - 1.)*self.wkp[ikp]
            return charge

            
if __name__ == '__main__':
    gamma = 1e-10
    mu_current = 10.7957625653
    filename='eigvals.dat'
    NOE = 12
    sdmu=0.1
    mix_mu = 1.0
    com = 1 # when computing density on imaginary axis, we can subtract Sigma_oo (com=0) or Sigma(0) (com=1)
    Qimag = True
    mix_mu = 1.
    
    (Ek, omega, wkp, nbands, nemin) = ReadEigenvals(filename, Qimag, gamma)

    valC = ValCharge(Qimag, wkp, omega, Ek, nemin, nbands, com)

    print(np.shape(Ek))
    
    print(valC(mu_current))# , Ek, omega, nemin, NOE, nbands))
    
    (mu0, mu1) = FindSignChange(valC, mu_current, sdmu, sys.stdout, args=())#Ek, omega, nemin, NOE, nbands))

    print((mu0, mu1))
    
    mu_new = optimize.brentq(valC, mu0, mu1, args=(Ek, omega, nemin, NOE, nbands))

    dmu = mu_new - mu_current
    mu = mu_current*(1-mix_mu) + mix_mu*mu_new
    print(mu_new)
    print('New chemical potential found at ', dmu, ' -> Chemical potential becomes ', mu)

