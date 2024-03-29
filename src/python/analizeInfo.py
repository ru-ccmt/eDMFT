#!/usr/bin/env python
from scipy import *
import glob
import sys

fi = open('info.iterate','r')

data = fi.readlines()
iic=[]
for line in data:
    sp = line.split()
    if sp[0]=='#': n=-1
    else: n=int(sp[2])
    iic.append(n)
ind=[]
for j in range(len(iic)-1):
    if iic[j+1]<=iic[j]: ind.append(j)
ind.append(len(data)-1)

for j in ind:
    print(data[j], end=' ')


if len(sys.argv)>1 and sys.argv[1]=='-f':
    scf_ = glob.glob('*.scf')
    if len(scf_)>0:
        scf = scf_[0]
        fs = open(scf,'r')
        data = fs.readlines()
        i=0
        index=1
        fos = open('Forces.dat', 'w')
        while i<len(data):
            line=data[i]
            #print '# ', line,
            if line[7:7+56]=='TOTAL FORCE WITH RESPECT TO THE GLOBAL COORDINATE SYSTEM':
                if (index in ind):
                    Print=True
                else:
                    Print=False
                if (Print): print(index, line, end=' ', file=fos)
                for j in range(i+1,len(data)):
                    line = data[j]
                    if line[:4]==':FGL':
                        if Print: print(line, end=' ', file=fos)
                    else:
                        break
                index+=1
                i=j
            i+=1
        fos.close()
