#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
from scipy import *

path=[[0,0,0],[-0.5,-0.5,0],[-1,0,0],[0,0,0]]

BS=[[1,1,0], [-1,1,0], [0,0,2]]
Ni = [8,8,8]
dN = [2,2,2]

BS = array(BS,dtype=float)

toSCALE=1./min(filter(lambda x: x>1e-6, abs(BS.flatten())))
SCALE = int(round(max(Ni)*toSCALE))

path = array(path)
l=0
for j in range(len(path)-1):
    for i in range(0,Ni[0],dN[j]):
        l+=1
        dp = path[j+1]-path[j]
        k = path[j] + dp*i/float(Ni[0])
        kint = map(round, k*SCALE)
        NAME = '   '+str(l)
        print "%-10s%5d%5d%5d%5d%5.1f" % (NAME, kint[0], kint[1], kint[2], SCALE, 1.0)
print 'END'

#for i in range(3):
#    print "%10.5f "*3 % tuple(BS[i,:])

