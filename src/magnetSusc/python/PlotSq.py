#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule

from scipy import *
from pylab import *
from scipy import interpolate
import glob, re

head = 'Sqw.'
om_max = 0.4
itensity=0.015

files = glob.glob(head+'*')
ii={}
for fil in files:
    m=re.match(head+'(\d+)', fil)
    if m is not None: ii[int(m.group(1))]=m.group(1)

#print sorted(ii.keys())

Sq=[]
for i in sorted(ii.keys()):
    file = head+str(ii[i])
    data = loadtxt(file).transpose()
    om = data[0]
    for l in range(len(om)):
        if om[l]>om_max: break
    Om=om[:l]
    Sq.append( interpolate.UnivariateSpline(Om,data[1,:l],s=0))

Nx = round(Om[-1]/(Om[1]-Om[0]))
omx = linspace(Om[0],Om[-1], Nx)
pSq=[]
for i in range(len(Sq)):
    pSq.append( Sq[i](omx) )
pSq=array(pSq)


vmm = [0,max(map(max,pSq))*itensity]
(ymin,ymax) = (Om[0],Om[-1])
(xmin,xmax) = (0, len(ii.keys()))
print vmm

imshow(pSq.transpose(), cmap=cm.jet, interpolation='bilinear', extent=[xmin,xmax,ymin,ymax], origin='lower', vmin=vmm[0], vmax=vmm[1],aspect=100.)
show()
