#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
from scipy import *

BS=[[1,1,0], [-1,1,0], [0,0,2]]
Ni = [100,100]



BS = array(BS,dtype=float)

toSCALE=1./min(filter(lambda x: x>1e-6, abs(BS.flatten())))
SCALE = int(round(max(Ni)*toSCALE))

l=0
for i in range(Ni[0]):
    for j in range(Ni[1]):
        l+=1
        x = array( [i/float(Ni[0]),j/float(Ni[1]), 0.0] )
        k = dot(BS.T, x)
        kint = map(round, k*SCALE)
        NAME = '   '+str(l)
        print "%-10s%5d%5d%5d%5d%5.1f" % (NAME, kint[0], kint[1], kint[2], SCALE, 1.0)
print 'END'

for i in range(3):
    print "%10.5f "*3 % tuple(BS[i,:])

