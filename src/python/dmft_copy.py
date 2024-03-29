#!/usr/bin/env python
""" Used for fast copying all relevant files to new working directory.
They are copied to the current directory from a source directory, which is the argument.
"""
import sys
import re, os, glob, shutil
from os.path import getsize

# DFT type quantities (w2k files)
f0 = ['.struct','.in0','.clmsum','.inm','.inM','.clmup','.clmdn']
f1 = ['.in1','.in1c','.klist','.inso','.in2c']
f2 = ['.in2','.kgen']
fc = ['.inc','.scf2','.incup','.incdn']
# DMFT type quantities
d0 = ['params.dat', 'EF.dat', 'sig.inp', 'Sigma.000', 'projectorw.dat', 'Eorb_imp.dat']
d1 = ['.indmfl','.indmfldn','.indmfi']
dn = ['sfx.*']

usage = """
 usage: % dmft_copy.py dir_name [options]

   dir_name      -- the name of the source directory with a wien2k or dmft run
                    files will be copied to the current directory
   options are:  -l  -- the files for DFT execution will be copied
                 -d  -- the files for DMFT but not DFT will be compied
                 -a  -- all relevant files for DFT+eDMFT are copied (Default)

 DFT file are  : """ + str(f0+f1+f2+fc) + """
 DMFT files are: """ + str(d0+d1+dn)+"\n"   

if len(sys.argv)<2 or (len(sys.argv)==2 and sys.argv[1]=='-h' or sys.argv[1]=='--help'):
    print(usage)
    sys.exit(0)
    
cpdr = sys.argv[1]

if not os.path.exists(cpdr):
    print('The directory '+cpdr+' does not exists!')
    sys.exit(1)

if len(sys.argv)==2:
    opt = '-a'
else:
    opt = sys.argv[2]
    if opt not in ['-l','-d','-a']:
        print('Switch can be either  -l,  -d or  -a')

strcase = glob.glob(cpdr+'/*.struct')
if not strcase:
    print('No structure file present. Nothing to copy...!')
    sys.exit(0)
strcase = strcase[0]
if not strcase:
    print('Can not find structure file!')
    sys.exit(1)

(dirName, fileName) = os.path.split(strcase); (fileBaseName, fileExtension)=os.path.splitext(fileName)
case = fileBaseName

print('case=', case)


w0 = [case+x for x in f0]
w1 = [case+x for x in f1]
w2 = [case+x for x in f2]
wc = [case+x for x in fc]

if opt in ['-l', '-a']:
    for f in w0+w1+w2+wc:
        fle = glob.glob(cpdr+'/'+f)
        if fle:
            print('Copying ... ', fle[0])
            shutil.copy2(fle[0], '.')
            

d1_ = [case+x for x in d1]
if opt in ['-d', '-a']:
    for f in d0+d1_:
        fle = glob.glob(cpdr+'/'+f)
        if fle:
            print('Copying ... ', fle[0])
            shutil.copy2(fle[0], '.')

    for f in dn:
        fle = glob.glob(cpdr+'/'+f)
        for fl in fle:
            print('Copying ... ', fl)
            shutil.copy2(fl, '.')
