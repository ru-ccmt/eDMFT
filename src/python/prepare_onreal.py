#!/usr/bin/env python
import os, shutil, sys
import glob
import subprocess

if len(sys.argv)<2:
    print('Give input file for onreal case.')
    sys.exit(0)
    
onreal = 'onreal'
dbin = os.environ['WIEN_DMFT_ROOT']

inpfile = sys.argv[1]
cases =  open(inpfile,'r').readlines()
cases = [line.strip() for line in cases if line.strip()]

cdir = os.getcwd()
if len(sys.argv)==2:
    for path in cases:
        path = path.strip()
        pth = os.path.normpath(path)
        cc = pth.split(os.sep)
        case_dir = '/'.join(cc[:-1])
        case = cc[0]
        ndir = case_dir+'/'+onreal
        os.makedirs(ndir, exist_ok=True)
        print(ndir)
        os.chdir(ndir)
        with open('onohup.dat', 'w') as fh_info:    
            cmd = dbin+'/dmft_copy.py ../'
            subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)
            shutil.copy2('../../'+case+'/'+case+'.klist_band', '.')
            shutil.copy2('../maxent/Sig.out', 'sig.inp')

            with open(case+'.indmfl', 'r') as inl:
                dmfl = inl.readlines()
                dmfl[1] = '0 0.025 0.025 300 -3.000000 3.000000  #\n'
            with open(case+'.indmfl', 'w') as inl:
                inl.writelines(dmfl)
            with open(case+'.indmfldn', 'r') as inl:
                dmfl = inl.readlines()
                dmfl[1] = '0 0.025 0.025 300 -3.000000 3.000000  #\n'
            with open(case+'.indmfldn', 'w') as inl:
                inl.writelines(dmfl)
            
            sklist  = open(case+'.klist').readlines()
            nk = len([line for line in sklist if line.strip()])-1
            sklistb = open(case+'.klist_band').readlines()
            nkb =len([line for line in sklistb if line.strip()])-1
            nproc = int(round(max(nk,nkb)/2+0.49999))
            
            scr="""#!/bin/bash
set -x
########################################################################
# SUN Grid Engine job wrapper
# parallel job on opteron queue
########################################################################
#$ -N j"""+case[2:]+"""
#$ -pe orte """+str(nproc)+"""
#$ -q gki36m5
#$ -j y
#$ -M haule@physics.rutgers.edu
#$ -m e
#$ -v WIEN_DMFT_ROOT,WIENROOT,LD_LIBRARY_PATH,PATH
########################################################################
# DON'T remove the following line!
source $TMPDIR/sge_init.sh
########################################################################
export SMPD_OPTION_NO_DYNAMIC_HOSTS=1
export MODULEPATH=/opt/apps/modulefiles:/opt/intel/modulefiles:/opt/gnu/modulefiles:/opt/sw/modulefiles
export SCRATCH="."
export OMP_NUM_THREADS=1
export PATH=.:$PATH
module load iompi/wien/19
module load intel/2024
module load intel/ompi
export WIEN_DMFT_ROOT=/home/haule/dbin
export LD_LIBRARY_PATH=/opt/intel/24.0/ompi/arpack/lib:/opt/intel/24.0/ompi/fftw-3.3.10-mpi/lib:/opt/intel/24.0/ompi/lib:/opt/intel/oneapi/2024.0/lib
echo "mpirun -n $NSLOTS" > mpi_prefix.dat

$WIENROOT/x_lapw lapw0 -f """+case+""" >& nohup.dat 
$WIEN_DMFT_ROOT/x_dmft.py lapw1        >& nohup.dat
$WIEN_DMFT_ROOT/x_dmft.py dmft1        >& nohup.dat
$WIEN_DMFT_ROOT/x_dmft.py lapw1 --band >& nohup.dat
$WIEN_DMFT_ROOT/x_dmft.py dmftp        >& nohup.dat
$WIEN_DMFT_ROOT/x_dmft.py dmftp -l dn  >& nohup.dat

"""
            with open('submit_real.scr','w') as fo:
                fo.write(scr)
        os.chdir(cdir)    
