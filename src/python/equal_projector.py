#!/usr/bin/env python
""" Small utility script to create an averaged projector, which might be useful when symmetry breaking 
on the lattice is not expected and projectors should be the same, but projectors are different due to 
numerical issues or not sufficiently converged charge. This is used only in connection with projector 5 and 6.

This utility reads all projectors and averages over them, and prints an averaged projector to projectorw.dat_new.
"""
# @Copyright 2007 Kristjan Haule
from numpy import *
import sys

if len(sys.argv)<2:
    print('Give the name of the projector file')
    sys.exit(0)

pfile = sys.argv[1]
fi=open(pfile,'r')
line1 = fi.readline()
Np, Nr = list(map(int,line1[1:].split()))

proj=[]
lines=[]
for p in range(Np):
    line2 = fi.readline()
    lines.append(line2)
    Nr = int(line2[1:].split()[0])
    dat=[]
    for i in range(Nr):
        line = fi.readline()
        dat.append( list(map(float,line.split())) )
    proj.append(dat)
proj = array(proj)
fi.close()

projf = zeros((shape(proj)[1],shape(proj)[2]))
for i in range(Np):
    projf[:,:] += proj[i,:,:]
projf *= 1.0/Np


fo = open(pfile+'_new', 'w')
print(line1, end=' ', file=fo)

for p in range(Np):
    print(lines[p], end=' ', file=fo)
    for i in range(Nr):
        print(projf[i,0], projf[i,1], projf[i,2], file=fo)
fo.close()

#from pylab import *
#plot(projf[:,0],projf[:,1])
#plot(projf[:,0],projf[:,2])
#show()





