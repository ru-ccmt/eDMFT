#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
from numpy import *
#from pylab import *
import sys

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
                #print machns[-1]
    return machns

def GetNumberOfKpoints(fklist):
    Nk=0
    with open(fklist,'r') as fk:
        for line in fk:
            if line[:3]=='END': break
            Nk += 1
    return Nk

def FindBestDistributionCores(machns, Nk):
    mold = machns[:]
    mnew = []

    OMP = int(round(len(machns)/float(Nk)))
    
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
        #print mold
        #print mnew
        if not mold: break
    return (mnew[:Nk], OMP)

if __name__ == '__main__':
    
    if (len(sys.argv)<2):
        print('Give two arguments: klist, machinefile!')
        sys.exit(0)

    fklist=sys.argv[1]
    fmachine=sys.argv[2]
    
    Nk = GetNumberOfKpoints(fklist)
    machns = GetMachines(fmachine)
    Ncpu_max = len(machns)
    print('Nk=', Nk, 'Ncpu_max=', Ncpu_max)
    print('machns=', machns)
    (Ncpu, Nrest) = FindNCpu(Nk,Ncpu_max)
    print('Ncpu=', Ncpu, 'Nrest=', Nrest)
    
    newmach, OMP = FindBestDistributionCores(machns, Nk)    


    fo = open('.machines','w')
    print('# machine file for Wien2K', file=fo)
    print('granularity:1', file=fo)
    
    shft=0
    #if Nrest>0: shft=1
    for i in range(shft,len(newmach)):
        print("1:"+newmach[i], file=fo)
    #if Nrest>0:
    #    print >> fo, 'residue:', newmach[0]

