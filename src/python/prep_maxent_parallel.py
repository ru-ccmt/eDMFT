#!/usr/bin/env python
import os, shutil, sys
import glob

if len(sys.argv)<2:
    print('Give input file for maxent case. If puting files back, also give directory with computed data.')
    sys.exit(0)
    
sdir = 'to_maxent'

inpfile = sys.argv[1]
cases =  open(inpfile,'r').readlines()
cases = [line.strip() for line in cases if line.strip()]

if len(sys.argv)==2:
    if not os.path.exists(sdir):
        os.mkdir(sdir)
    
    Is_params_here = os.path.isfile('maxent_params.dat')
    for path in cases:
        path = path.strip()
        pth = os.path.normpath(path)
        cc = pth.split(os.sep)
        case = '_'.join(cc[:-1])
        ndir = sdir+'/'+case
        os.makedirs(ndir, exist_ok=True)
        if os.path.exists(path+'/sig.inpx'):
            shutil.copy2(path+'/sig.inpx', ndir)
        else:
            print('ERROR: File', path+'/sig.inpx', 'missing')
        if Is_params_here:
            shutil.copy2('maxent_params.dat', ndir)
        else:
            shutil.copy2(path+'/maxent_params.dat', ndir)

elif len(sys.argv)==3:
    jobdir = sys.argv[2]
    for path in cases:
        path = path.strip()
        pth = os.path.normpath(path)
        cc = pth.split(os.sep)
        case = '_'.join(cc[:-1])
        ndir = jobdir+'/'+sdir+'/'+case
        if not os.path.exists(ndir):
            print('WARNING: Can not fine dir', ndir)
        else:
            fls = glob.glob(ndir+'/*')
            for file in fls:
                cmp = file.split(os.sep)
                if cmp[-1]!='sig.inpx':
                    cmd = 'copy '+file+' -> '+path+'/'
                    print(cmd)
                    shutil.copy2(file, path)
else:
    print('Not yet assigned job')
