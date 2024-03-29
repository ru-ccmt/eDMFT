#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import numpy as np
import sys

def findNbands(Emin,Emax,enefiles,strfile, log=sys.stdout):
    Ry2eV = 13.6056923
    # Find 'nat' in the structure file
    
    with open(strfile,'r') as fs:
        next(fs)
        line = next(fs)
        lattic = line[:4]
        nat = int(line[4+23:4+23+3])
    
    print('Number of all atoms found in struct file', nat, file=log)
    nemin=10000
    nemax=0
    for enefile in enefiles:
        # Find nemin,nemax in energy file
        fi = open(enefile,'r')
        for i in range(nat):
            next(fi) # linearization Energy
            next(fi) # linearization Energy
        try:
            for k in range(1,1000000):
                line = next(fi)
                S,T,Z = float(line[:19]),float(line[19:2*19]),float(line[2*19:3*19])
                KNAME = line[3*19:3*19+10]
                N, NEn = int(line[67:67+6]), int(line[67+6:67+6*2])
                nemin_=1
                nemax_=0
                for ii in range(NEn):
                    line = next(fi).split()
                    num, e1 = int(line[0]), float(line[1])
                    e1 *= Ry2eV
                    if (e1<Emin): nemin_ += 1
                    if (e1<Emax): nemax_ += 1
                nemin = min(nemin,nemin_)
                nemax = max(nemax,nemax_)
        except StopIteration:
            fi.close()
        print('file:', enefile, 'nemin=', nemin, 'nemax=', nemax, file=log)
    print('Finally set nemin=', nemin, 'nemax=', nemax, file=log)
    return (nemin,nemax)

if __name__ == '__main__':
    import os, sys, re, glob
    import utils
    Ry2eV = 13.6056923

    if len(sys.argv)<3:
        exmin=-10
        exmax= 10
    else:
        exmin=float(sys.argv[1])
        exmax=float(sys.argv[2])
    print('Energy window:', exmin, exmax)

    w2k = utils.W2kEnvironment()
    
    # looking for EF
    if os.path.isfile('EF.dat'):
        EF = float(open('EF.dat').read())
    else:
        fname = w2k.case+".scf2"
        if os.path.isfile(fname) or os.path.isfile(fname+'up'):
            if os.path.isfile(fname):
                fscf = open(fname, 'r')
            else:
                fscf = open(fname+'up', 'r')
            
            lines = fscf.readlines()
            for line in lines:
                if re.match(r':FER', line) is not None:
                    EF = float(line[38:])*Ry2eV
                    print('EF from scf file : ', EF)
                    break
        else:
            EF =float(open(w2k.case+'.indmf1').readlines()[1].split()[1])
            print('EF from indmf1 file : ', EF)
        
    print('EF=', EF)
    
    #Emin,Emax = -1.331295, 18.668705
    Emin, Emax = EF+exmin, EF+exmax
    print('Emin, Emax=', Emin, Emax)
    strfile = w2k.case+'.struct'
    
    enefiles = glob.glob(w2k.case+'.energy'+'*')
    enefiles = [fil for fil in enefiles if os.path.getsize(fil)>0] # Remove empty files
    for fil in enefiles:
        if re.match(w2k.case+'.energyso', fil): # Spin-orbit on, remove non-spin-orbit files
            enefiles = [fil for fil in enefiles if re.match(w2k.case+'.energyso', fil) is not None] # Remove empty files
            break

    print('enefiles=', enefiles)

    (nemin,nemax) = findNbands(Emin,Emax,enefiles,strfile)
    print('nemin,nemax=', nemin, nemax)
    print('Replace second line of '+w2k.case+'.indmfl with')
    print(nemin,nemax,1,5,'# hybridization nmin, nmax, renormalize for interstitials, projection type')

