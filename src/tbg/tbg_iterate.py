#!/usr/bin/env python
import sys, subprocess
import re, os, glob, shutil, socket, time
from os.path import getsize
from numpy import *


if __name__ == '__main__':
    exec(open("params.py").read())
    dirs = par['dirs']
    log = open('info.iterate','w')

    ROOT = '.'
    #ROOT = os.environ['TBG_ROOT']
    ROOT = os.environ['WIEN_DMFT_ROOT']
    
    for itt in range(par['Nitt']):
        exec(open("params.py").read())  # updating parameters if changed
        Nc = par['Nc']
        #n0 = Nc-16
        #dch = par['dch']
        Sig_DC = iparams0["U"][0] * (Nf-par['dch'])
        print('Sig_DC=', Sig_DC, file=log)
        # DMFT Self-consistency condition
        sname = 'tbg_dmft_scc.py'
        print(time.strftime("%a %b %d %H:%M:%S %Z %Y")+'>     '+sname, file=log)
        cmd = ROOT+'/'+sname
        print(('#<'+sname+'>: '), cmd, file=log)
        log.flush()
        info=subprocess.call(cmd,shell=True,stdout=log)
        
        if os.path.isfile('EF.dat'):
            mu = float(open('EF.dat').readline())
            print(':FER ', mu, file=log)

        # Exact diagonalization of the atom
        Eds=[]
        sname = ROOT+'/atom_d.py'
        print(time.strftime("%a %b %d %H:%M:%S %Z %Y")+'>     '+sname, file=log)
        for idi in range(len(dirs)):
            line = open(dirs[idi]+'/Eimp','r').readlines()[0]
            Ed = array(eval(line))
            Eds.append(Ed)
            Edc = Ed-Ed[0] # what we put into cix is only the splitting
            cmd = 'cd '+dirs[idi]+'; '+sname+' l=1.5 "Eimp='+str(Edc.tolist())+'"'
            print(('#<'+sname+'>: '), cmd, file=log)
            log.flush()
            info=subprocess.call(cmd,shell=True,stdout=log)
        
        # Preparing CTQMC input file
        for idi in range(len(dirs)):
            mu_QMC = -Eds[idi][0] + Sig_DC
            fp = open(dirs[idi]+'/PARAMS','w')
            for p in list(iparams0.keys()):
                print(p, '  ', iparams0[p][0], iparams0[p][1], file=fp)
            print('mu    ', mu_QMC, file=fp)
            fp.close()
        
        sname = ROOT+'/ctqmc'
        print(time.strftime("%a %b %d %H:%M:%S %Z %Y")+'>     '+sname, file=log)
        for idi in range(0,1): #len(dirs)):
            cmd = 'cd '+dirs[idi]+';'+mpi_prefix+' '+sname+' >& nohup.dat '
            print(('#<'+sname+'>: '), cmd, file=log)
            log.flush()
            info=subprocess.call(cmd,shell=True,stdout=log)
        
            Sigdata = loadtxt(dirs[idi]+'/'+iparams0['Sig'][0])
            for i in range(1,Sigdata.shape[1],2):
                Sigdata[:,i] -= Sig_DC
            savetxt(dirs[idi]+'/'+par['Sigma'], Sigdata)
            
            shutil.copy2(dirs[idi]+'/'+iparams0['Sig'][0], dirs[idi]+'/Sig.out.'+str(itt))
            shutil.copy2(dirs[idi]+'/Gf.out', dirs[idi]+'/Gf.out.'+str(itt))
            shutil.copy2(dirs[idi]+'/'+par['Sigma'], dirs[idi]+'/'+par['Sigma']+'.'+str(itt))
            
