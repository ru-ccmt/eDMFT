#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import os, sys
from numpy import *
import utils
from w2k_atpar import readpotential, readlinearizatione, atpar, rint13
from wstruct import Struct
import optparse
import re
from indmffile import Indmfl

def SolveForContinuousFunction(A,Ae,Aee,Rx,Nr0,Nr):
    """
    Routine extends solution beyond Rmt by making value and derivative continuous and make it vanish ar R2
    """
    dA = (A[Nr0]-A[Nr0-1])/(Rx[Nr0]-Rx[Nr0-1])    # Radial derivative of ul(r)
    dAe= (Ae[Nr0]-Ae[Nr0-1])/(Rx[Nr0]-Rx[Nr0-1])  # Radial derivative of dot{ul(r)}
    dAee=(Aee[Nr0]-Aee[Nr0-1])/(Rx[Nr0]-Rx[Nr0-1])# Radial derivative of dot{dot{ul(r)}}

    am=array([[A[Nr0-1], Ae[Nr0-1], Aee[Nr0-1]],  # 3x3 equation for continuous solution
              [A[Nr-1],  Ae[Nr-1],  Aee[Nr-1]],   # for r<Rmt, we take solution of Dirac equation
              [dA, dAe, dAee]])                   # outside we use cobination of ul(r),dot{ul(r)},dot{dot{ul(r)}}
    bm=array([A[Nr0-1], 0., dA])                  # 3 Eqs.: value at Rmt, value at R2, derivative at Rmt
    cm=linalg.solve(am,bm)

    return cm

def SolveForLocalizedFunction(A,Ae,Aee,Rx,Nr):
    dA = (A[Nr-1]-A[Nr-2])/(Rx[Nr-1]-Rx[Nr-2])    # Radial derivative of ul(r)
    dAe= (Ae[Nr-1]-Ae[Nr-2])/(Rx[Nr-1]-Rx[Nr-2])  # Radial derivative of dot{ul(r)}
    dAee=(Aee[Nr-1]-Aee[Nr-2])/(Rx[Nr-1]-Rx[Nr-2])# Radial derivative of dot{dot{ul(r)}}

    am=array([[ Ae[Nr-1], Aee[Nr-1]],
              [ dAe,     dAee     ]])
    bm=array([-A[Nr-1],-dA])
    cm=linalg.solve(am,bm)
    return cm

def SolveForLocalizedFunction3(A,Ae,Nr):
    return -A[Nr-1]/Ae[Nr-1]

def ReadIndmfl(case):
    inl = Indmfl(case)
    inl.read(case+'.indmfl')
    atms = list(inl.atoms.keys())
    Rmt2 = [atm_dat[3] for iatom,atm_dat in inl.atoms.items()] # inl.atoms[iatom] contains (locrot_shift, new_xyz, shift_vec, Rmt2)
    return (atms, inl.Lsa, inl.icpsa, Rmt2)
        
def FindChemicalPotential(case, updn, logf=sys.stdout):
    # Looking for the LDA chemical potential                                                                                                                                                                
    Ry2eV = 13.60569193
    EF_found = False

    fname = case+".scf2"
    if os.path.isfile(fname):
        fscf = open(fname, 'r')
        print("opening", fname, file=logf)
    elif os.path.isfile(fname+updn):
        fscf = open(fname+updn, 'r')
        print("opening", fname+updn, file=logf)
    else:
        fscf = None

    if fscf is not None:
        lines = fscf.readlines()
        for line in lines:
            if re.match(r':FER', line) is not None:
                EF = float(line[38:])
                print('EF=', EF*Ry2eV, 'eV = ', EF, 'Ry', file=logf)
                EF_found = True

    # The previous DMFT chemical potential                                                                                                                                                                  
    if  os.path.isfile('EF.dat'):
        fmu = open('EF.dat','r')
        mu = float(next(fmu))
        EF = mu/Ry2eV
        print('Found DMFT-EF=', EF*Ry2eV, 'eV', file=logf)
        EF_found = True

    if not EF_found:
        fname = case+".scf"
        if os.path.isfile(fname):
            print("opening", fname, file=logf)
            f = open(fname, 'r')
            lines = f.readlines()
            if not EF_found:
                for line in lines:
                    if re.match(r':FER', line) is not None:
                        EF = float(line[38:])
                        print('EF=', EF*Ry2eV, file=logf)
                        EF_found = True

    if not EF_found:
        raise Exception("Failed to determine chemical potential.")

    return EF

def main(case, atms, Lsa, icpsa, Rm2, localize=1, Emu='EF', logf=sys.stdout):
    so=''
    if os.path.isfile(case+".inso") and os.path.getsize(case+".inso")>0 :
        print('Found '+case+'.inso file, hence assuming so-coupling exists. Switching -so switch!', file=logf)
        so='so'
    potential_file = case+'.vsp'
    vector_file = case+'.vector'+so
    struct = Struct(case)
    struct.ReadStruct(case+'.struct')
    rel=True
    if struct.mode[:4]=='NREL': rel=False
    
    atm_l_case={}
    latom=0
    for jatom in range(struct.nat):
        ll_case=[]
        if (latom+1) in atms:           # This atom listed in indmfl file. Note here latom starts from 0, while in indmfl starts with 1
            ia = atms.index(latom+1)    # It appears as the ia consequitive atom in indmfl file
            for il in range(len(Lsa[ia])):
                if icpsa[ia][il]>0: # icix>0, hence it is correlated
                    ll_case.append(Lsa[ia][il])
        latom += struct.mult[jatom]
        if ll_case:
            atm_l_case[jatom]= (ll_case, Rm2[ia])
    print('atm_l_case=', atm_l_case, file=logf)

    Nrmax=0
    for jatom in list(atm_l_case.keys()):
        Rmt2 = atm_l_case[jatom][1]
        for lc in atm_l_case[jatom][0]:
            Nr0 = struct.jrj[jatom]  # Number of radial points from struct file
            Rmt = struct.rmt[jatom]
            if Rmt2>Rmt:
                # creating larger mesh to Rmt2, which can be larger then Rmt
                r0 = struct.r0[jatom]
                if Rmt2<Rmt: Rmt2 = Rmt 
                Nr = 1 + int( (Nr0-1)*log(Rmt2/r0)/log(Rmt/r0) + 0.99)  # Number of radial points after extending mesh to Rmt2
            else:
                Nr=Nr0
            if Nr>Nrmax: Nrmax=Nr

    lmax2 = 3

    Vr_all = readpotential(potential_file, max(struct.jrj), struct.jrj, struct.nat)               # Reading actual potential
    
    if Emu!='EF':
        E_all = readlinearizatione(vector_file, struct.nat, lmax2)      # Linearization energies from vector file
    else:
        updn=''
        EF = FindChemicalPotential(case, updn, logf=logf)
    
    fout = open('projectorw.dat', 'w')
    print('#', sum([len(atm_l_case[jatom][0]) for jatom in list(atm_l_case.keys())]), Nrmax, file=fout)
    
    #####################################
    projw=[]
    for jatom in list(atm_l_case.keys()):
        Rmt2 = atm_l_case[jatom][1]
        for lc in atm_l_case[jatom][0]:
            Nr0 = struct.jrj[jatom]  # Number of radial points from struct file
            # creating larger mesh to Rmt2, which can be larger then Rmt
            r0 = struct.r0[jatom]
            Rmt = struct.rmt[jatom]
            if Rmt2<Rmt: Rmt2 = Rmt
            Nr = 1 + int( (Nr0-1)*log(Rmt2/r0)/log(Rmt/r0) + 0.99)  # Number of radial points after extending mesh to Rmt2
            dx = log(Rmt/r0)/(Nr0-1)
            Rx = r0 * exp( arange(Nr)*dx )                        # First radial point
            
            Vr = hstack( (Vr_all[:,jatom], ones(Nr-Nr0)*Vr_all[-1,jatom] ) )             # Extending it in the interstitial with constant

            if Emu=='EF':
                Elc=EF
            else:
                El = E_all[:,jatom]
                print('El0=', El, file=logf)
                lapw = ones(len(El), dtype=bool)
                for l in range(len(El)):
                    if El[l]>150:
                        El[l]-=200
                        lapw[l]=False
                print('linearization E=', El, file=logf)
                Elc = El[lc]
                
            print('dx=', dx, file=logf)
            print('Rx[:]=', Rx[0], Rx[Nr0-1], Rx[Nr-1], file=logf)
            print('using E=', Elc, file=logf)
                
            A,B,Ae,Be,Aee,Bee,Pei = atpar(rel,lc,Vr,Rx,Elc,struct.r0[jatom],dx,struct.Znuc[jatom])
            if localize==1:
                if Nr>Nr0:
                    # Below we extend solution beyond Rmt by making value and derivative continuous and make it vanish ar Rmt2
                    cm = SolveForContinuousFunction(A,Ae,Aee,Rx,Nr0,Nr)
                    Ag = hstack(( A[:Nr0] ,  cm[0]*A[Nr0:]+cm[1]*Ae[Nr0:]+cm[2]*Aee[Nr0:] ))
                    cm = SolveForContinuousFunction(B,Be,Bee,Rx,Nr0,Nr)
                    Bg = hstack(( B[:Nr0] ,  cm[0]*B[Nr0:]+cm[1]*Be[Nr0:]+cm[2]*Bee[Nr0:] ))
                
                    # We then normalize the vawe function such that it is normalized within Rmt
                    #overlap = rint13(rel,Ag[:Nr0],Bg[:Nr0],Ag[:Nr0],Bg[:Nr0],dx,struct.r0[jatom])
                    overlap = rint13(rel,Ag[:Nr],Bg[:Nr],Ag[:Nr],Bg[:Nr],dx,struct.r0[jatom])
                    Ag *= 1/sqrt(overlap)  # normalization inside Rmt, and not Rmt2
                    Bg *= 1/sqrt(overlap)
                    A *= 1/sqrt(overlap)
                    B *= 1/sqrt(overlap)
                    print('overlap=', overlap, file=logf)
                else:
                    Ag = A
                    Bg = B
            elif localize==2:
                cm = SolveForLocalizedFunction(A,Ae,Aee,Rx,Nr)
                Ag = A + cm[0]*Ae + cm[1]*Aee
                cm = SolveForLocalizedFunction(B,Be,Bee,Rx,Nr)
                Bg = B + cm[0]*Be + cm[1]*Bee
                overlap = rint13(rel,Ag[:Nr],Bg[:Nr],Ag[:Nr],Bg[:Nr],dx,struct.r0[jatom])
                Ag *= 1/sqrt(overlap)  # normalization inside Rmt, and not Rmt2
                Bg *= 1/sqrt(overlap)
                print('cm=', cm, file=logf)
                print('overlap=', overlap, file=logf)
            elif localize==3:
                cm = SolveForLocalizedFunction3(A,Ae,Nr)
                Ag = A + cm*Ae
                cm = SolveForLocalizedFunction3(B,Be,Nr)
                Bg = B + cm*Be
                overlap = rint13(rel,Ag[:Nr],Bg[:Nr],Ag[:Nr],Bg[:Nr],dx,struct.r0[jatom])
                Ag *= 1/sqrt(overlap)  # normalization inside Rmt, and not Rmt2
                Bg *= 1/sqrt(overlap)
                print('cm=', cm, file=logf)
                print('overlap=', overlap, file=logf)
            else:
                print('Dot yet implemented!', file=logf)
                sys.exit(1)
            ############################################################################
            print('#', Nr, Nr0, jatom+1, lc, file=fout)
            for ir in range(Nr):
                print(Rx[ir], Ag[ir], Bg[ir], file=fout)
            projw.append( (Rx,Ag,Bg,A,B) )
    fout.close()
    return projw


if __name__ == '__main__':
    usage = """
    usage: %prog [ options ]

    Saved the projector wave function into projectorw.dat file.
    There are two Rmt used in this code. The first Rmt1 is muffin-thin sphere
    from structure file. The second is Rmt2 and is written in case.indmfl file.
    """
    parser = optparse.OptionParser(usage)
    parser.add_option("-l", "--localize",  dest="localize",    type="int", default=1, help="How localized should be the function")
    # Next, parse the arguments
    (options, args) = parser.parse_args()

    w2k = utils.W2kEnvironment()

    atms, Lsa, icpsa, Rm2 = ReadIndmfl(w2k.case)

    choice = input("""Should linearization energy for the projector be\n    [1]: EF  or\n    [2]: DFT linearization energy.\n  Enter [1|2]: """)
    if choice.strip()=='1':
        Emu='EF'
    else:
        Emu='Emu'
        
    projw = main(w2k.case, atms, Lsa, icpsa, Rm2, options.localize, Emu)
    from pylab import *

    for Rx,Ag,Bg,A,B in projw:
        subplot(2,1,1)
        plot(Rx, A, 'r-')
        plot(Rx, Ag, 'y-')
        subplot(2,1,2)
        plot(Rx, B, 'r-')
        plot(Rx, Bg, 'y-')
        show()
        
        ### self.atoms.keys()  -> atms
        ### self.Lsa  -> Lsa
        ### self.icpsa -> icpsa
        ### [self.atoms[iatom][3] for iatom in self.atoms.keys()]  -> Rmt2
        
