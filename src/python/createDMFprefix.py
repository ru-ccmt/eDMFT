#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
from numpy import *
import re, sys
import builtins

def FindNCpu(Nk,Ncpu_max):
    for Ncpu in range(Ncpu_max,1,-1):
        if (Nk-int(Nk/Ncpu)*Ncpu < Nk/Ncpu):
            Nrest = Nk-(Nk/Ncpu)*Ncpu
            break
    if Nrest>0 and Ncpu==Ncpu_max:
        for Ncpu in range(Ncpu_max-1,1,-1):
            if (Nk-int(Nk/Ncpu)*Ncpu < Nk/Ncpu):
                Nrest = Nk-(Nk/Ncpu)*Ncpu
                break
        return (Ncpu, Nrest)
    else:
        return (Ncpu, Nrest)

def GetMachines(fmachine):
    machns=[]
    with open(fmachine,'r') as fm:
        for line in fm:
            if line.strip():
                machns.append(line.strip())
    return machns

def GetNumberOfKpoints(fklist):
    Nk=0
    with open(fklist,'r') as fk:
        for line in fk:
            if line[:3]=='END': break
            Nk += 1
    return Nk

def FindBestDistributionCores(machns, Nk):
    Ncpu = len(machns)
    mold = machns[:]
    mnew = []
    OMP = int(round(Ncpu/float(Nk)))

    #print 'cpus=', len(machns), 'Nk=', Nk, 'OMP=', OMP
    
    if OMP>=16:  # 1 jobs per cpu (16 cores)
        OMP=16
    elif OMP>=8: # 2 jobs per cpu
        OMP=8
    elif OMP>5:  # 3 jobs per cpu
        OMP=5
    elif OMP>=4: # 4 jobs per cpu
        OMP=4
    elif OMP>3:  # 5 jobs per cpu
        OMP=3
    elif OMP>=2:  # 8 jobs per cpu
        OMP=2
    else:
        OMP=1    # 7-8 jobs per cpu
    
    while  len(mnew)<Nk :
        Of = len(mold)/float(Nk-len(mnew))
        Om = int(round(Of))
        #print 'Of=', Of, 'Om=', Om
        if Om>0:
            add = mold[::Om]
        else:
            add = mold[:]
        mnew = mnew + add
        for i in range(len(add)): mold.remove(add[i])
        #print 'mold=', mold
        #print 'mnew=', mnew
        if not mold: break  # mold empty

    #print 'OMP=', OMP
    #print 'mnew[:Nk]=', mnew
    return (mnew[:Nk], OMP)
    

if __name__ == '__main__':
    
    if (len(sys.argv)<2):
        print('Give two arguments: case.klist, mpi_prefix.dat!')
        sys.exit(0)

    machine_out='wmachines'
    
    fklist=sys.argv[1]
    fprefix=sys.argv[2]
    mpi_prefix = open(fprefix).readline().strip()
    
    Nkp = GetNumberOfKpoints(fklist)
    # Finds the name of the machinefile from mpi_run.dat directive
    m = re.search('-machinefile\s+([\w|.|/]*)',mpi_prefix)
    if m is not None:
        fmachine = m.group(1)
        machns = GetMachines(fmachine)
        Ncpu = len(machns)
        newmach, OMP = FindBestDistributionCores(machns, Nkp)

        if OMP<=1:
            OMP = 1
            Nn = builtins.min(Nkp,Ncpu)
            mpin = re.sub(r'(-n[p]?)\s*(\d+)', r'\1 '+str(Nn), mpi_prefix)
        else:
            # creates a new machinefile (wmachines) with fraction of cores in the list
            with open(machine_out,'w') as fm:
                for i in range(len(newmach)):
                    print(newmach[i], file=fm)
            # replaces -machinefile directive in mpi_prefix.dat with the name of the new machinefile
            mpin = re.sub('-machinefile\s+([\w|.|/]*)', '-machinefile '+machine_out, mpi_prefix)
            # reduced the number of cores for MPI run (-np directive) in mpi_prefix.dat
            m = re.search('-n[p]?\s*(\d+)',mpin)
            if m is not None:
                np = int(m.group(1))
                mpin = re.sub('(-n[p]?)\s*(\d+)', '-\1 '+str(len(newmach)), mpin)
                
    else:
        m = re.search('-n[p]?\s+(\d+)',mpi_prefix)
        OMP = 1
        if m is not None:
            Ncpu = int(m.group(1))
            # we can not fine tune
            Nn = builtins.min(Nkp,Ncpu)
            mpin = re.sub(r'(-n[p]?\s*)(\d+)', r'\1 '+str(Nn), mpi_prefix)
    

    m = re.search('-env\s*OMP_NUM_THREADS\s*(\d+)',mpin)
    if m is not None:
        mpin = re.sub('-env\s*OMP_NUM_THREADS\s*\d+', '-env OMP_NUM_THREADS '+str(OMP), mpin)
    else:
        mpin += ' -env OMP_NUM_THREADS '+str(OMP)
    print(mpin)
    
