#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import sys, subprocess
import re, os, glob, shutil, socket, time
from os.path import getsize
from scipy import *
import numpy
import copy
import glob
from datetime import date,datetime

def Step(cmd, fh_info):
    dt = datetime.now().strftime("%B %d  %H:%M:%S %Y")
    logs = open(':log', 'a')
    logs.write(dt+' '+cmd+'\n')
    logs.close()
    print(cmd)
    fh_info.flush()
    subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)

def Diff(case):
    fs = open(case+'.scf', 'r')
    dat = fs.readlines()
    fs.close()

    DEne=[]
    DChr=[]
    for line in dat:
        if re.match(r':ENE', line):
            DEne.append(float(line.split()[-1]))
            #print line.split()[-1]
        if re.match(r':DIS', line):
            DChr.append(float(line.split()[-1]))
            #print line.split()[-1]
    if len(DChr)>1:
        drho = DChr[-1]
    else:
        drho = 1.0
    if len(DEne)>1:
        dene = abs(DEne[-1]-DEne[-2])
    else:
        dene = 1.0
    if dene==0: dene = 1.0    
    print(':ENERGY convergence: '+str(dene))
    print(':CHARGE convergence: '+str(drho))
    return (drho, dene)

def CheckError(step):
    if os.path.exists(step+'.error') and os.path.getsize(step+'.error')>0 :
        message = open(step+'.error','r').readlines()
        print(message)
        return True
    else:
        return False

def scf(case, wopt, outfile='.scf'):
    def SearchScf(ist,KEY,dat):
        for i in range(ist,len(dat)):
            line=dat[i].strip()
            if line[:len(KEY)]==KEY:
                break
        return (line,i)
    
    def GetSum(KEY, dat_up, dat_dn):
        line_up, iu_ = SearchScf(0,KEY,dat_up)
        line_dn, id_ = SearchScf(0,KEY,dat_dn)
        if iu_<len(dat_up) and id_<len(dat_dn):
            del dat_up[iu_]
            del dat_dn[id_]
            aSUM = 0.5*(float(line_up[30:50])+float(line_dn[30:50]))
            return line_up[:30]+("%20.9f"%aSUM)+line_up[50:]
        else:
            return ' '
    
    # First merger case.scf0, case.scf1 and case.scfso
    dat=[]
    scfs = ['0','1','so']
    if wopt['updn']: scfs = ['0','1up','1dn','so']
    for i in scfs:
        if os.path.exists(case+'.scf'+i):
            dat += '\n------- '+case+'.scf'+i+' --------\n'
            dat += open(case+'.scf'+i,'r').readlines()
    
    # case.scf2 is complicated when both up & dn are present
    # We need to average over case.scf2 and case.scf2dn
    if wopt['m_extn'] and os.path.exists(case+'.scf2'+wopt['m_extn']):
        lineo, io = SearchScf(0,':NATO',dat)
        nat = int(lineo.split(':')[2].split()[0])
        
        up__ = 'up' if wopt['updn'] else ''
        dat_up = open(case+'.scf2'+up__,'r').readlines()
        dat_dn = open(case+'.scf2'+wopt['m_extn'],'r').readlines()
        fu=case+'.scf2'+up__
        fd=case+'.scf2'+wopt['m_extn']
    
        line_s1 = GetSum(':SUM ', dat_up, dat_dn)
        line_s2 = GetSum(':XSUM ', dat_up, dat_dn)
        line_s3 = GetSum(':YSUM ', dat_up, dat_dn)
        line_s4 = GetSum(':ZSUM ', dat_up, dat_dn)
        line_s5 = GetSum(':WSUM ', dat_up, dat_dn)
        line_s6 = GetSum(':QSUM ', dat_up, dat_dn)
        line_s7 = GetSum(':VSUM ', dat_up, dat_dn)
        
        dat += '\n------- '+fu+' and '+fd+' --------\n'
    
        lineu, iu_ = SearchScf(0,':NOE ',dat_up)
        lined, id_ = SearchScf(0,':NOE ',dat_dn)
        dat += '\n------- '+fu+'  -------\n'
        dat += dat_up[:iu_]
        dat += lineu
        dat += '\n------- '+fd+'  -------\n'
        dat += dat_dn[:id_]
        dat += lined
    
        lineu, _iu = SearchScf(iu_,':FER ',dat_up)
        lined, _id = SearchScf(id_,':FER ',dat_dn)
        dat += '\n------- '+fu+'  -------\n'
        dat += dat_up[iu_+1:_iu]
        dat += lineu
        dat += '\n------- '+fd+'  -------\n'
        dat += dat_dn[id_+1:_id]
        dat += lined
        iu_,id_ = _iu,_id
        
        lineu, _iu = SearchScf(iu_,':CHA ',dat_up)
        lined, _id = SearchScf(id_,':CHA ',dat_dn)
        dat += '\n------- '+fu+'  -------\n'
        dat += dat_up[iu_+1:_iu]
        dat += '\n------- '+fd+'  -------\n'
        dat += dat_dn[id_+1:_id]
        iu_,id_ = _iu,_id
        chup = float(lineu.split()[-1])
        chdn = float(lined.split()[-1])
        chtot = 0.5*(chup+chdn)
        dat += lineu[:40]+("%14.6f"%chtot)
    
    
        lineu, _iu = SearchScf(iu_,':DRHO',dat_up)
        lined, _id = SearchScf(id_,':DRHO',dat_dn)
            
        dat += '\n------- '+fu+'  -------\n'
        dat += dat_up[iu_+1:_iu]
        dat += '\n------- '+fd+'  -------\n'
        dat += dat_dn[id_+1:_id]
        iu_,id_ = _iu,_id
        
        drho_up = float(lineu.split()[1])
        drho_dn = float(lined.split()[1])
        #drho_up = float(lineu[23:39])        
        #drat_up = float(lineu[56:])
        #drho_dn = float(lined[23:39])        
        #drat_dn = float(lined[56:])
        drho = 0.5*(drho_up+drho_dn)
        #drat = 0.5*(drat_up+drat_dn)
        dat += lineu[:5]+("%20.12f "%drho)+'\n\n'
        if (abs(drho) > 0.5):
            print "WARNING : Electron charge density is very different than expected. You should change recomputeEF and set it to 1 (switch in on)"
            
        dat += line_s1+'\n'+line_s2+'\n'+line_s3+'\n'+line_s4+'\n'+line_s5+'\n'+line_s6+'\n'+line_s7+'\n\n'
        
    
        # merging forces
        lineu, iu_ = SearchScf(0,'VALENCE-FORCE',dat_up)
        lined, id_ = SearchScf(0,'VALENCE-FORCE',dat_dn)
        if iu_<len(dat_up)-1 and id_<len(dat_dn)-1:  # Forces are present
            for KEY in [':FVA',':FSU']:
                for iat in range(nat):
                    lineu, _iu = SearchScf(iu_,KEY,dat_up)
                    lined, _id = SearchScf(id_,KEY,dat_dn)
                    for i in range(iu_,_iu): # Printing what is between two forces :FVA/:FSU
                        dat += dat_up[i]
                    iu_ = _iu+1
                    id_ = _id+1
                    dlu = lineu.split()
                    dln = lined.split()
                    # average valence force between up and dn
                    Fa = [0.5*(float(dlu[i])+float(dln[i])) for i in range(3,6)]
                    # print out average valence force
                    dat += dlu[0]+('%9s'%dlu[1])+('%17.9f'*4) % tuple([sqrt(dot(Fa,Fa))]+Fa),'\n'
    else:
        # for non-magnetic calculation it is very simple
        if os.path.exists(case+'.scf2'):
            dat += '\n------- '+case+'.scf2'+' --------\n'
            dat += open(case+'.scf2','r').readlines()
            
    if os.path.exists('Eorb_imp.dat'):
        dat += '\n---------- Eorb_imp.dat --------------\n'
        dat += open('Eorb_imp.dat','r').readlines()
        
    scfs = ['1s','2s','c']
    if wopt['updn']: scfs = ['cup','cdn']
    for i in scfs:
        if os.path.exists(case+'.scf'+i):
            dat += '\n------- '+case+'.scf'+i+' --------\n'
            dat += open(case+'.scf'+i,'r').readlines()

    fs = open(case+outfile, 'a')
    fs.writelines(dat)
    fs.close()

def scfm(case):
    dat=[]
    for i in ['m']:
        dat += '\n------- '+case+'.scf'+i+' --------\n'
        if os.path.exists(case+'.scf'+i):
            dat += open(case+'.scf'+i,'r').readlines()
            
    fs = open(case+'.scf', 'a')
    fs.writelines(dat)
    fs.close()
    
def scf1(case):
    for i in ['clmsum','vsp','vns','vrespsum','vspdn','vnsdn','vspup','vnsup','clmup','clmdn']:
        name = case+'.'+i
        if os.path.exists(name) and os.path.getsize(name)>0:
            shutil.copy2(name, name+'_old')		#save last cycle

if __name__ == '__main__':

    p={
        'cc' : 1e-3,      # the charge density precision to stop the LDA+DMFT run
        'ec' : 1e-4,      # the energy precision to stop the LDA+DMFT run
        'Nitt' : 1000     # total number of iterations
    }

    wopt={'cmplx':'', 'so':'', 'm_extn':'', 'updn':''}  
    # cmplx : no-inversion, 
    # so    : spin-orbit coupling present
    # m_extn: dmft is spin polarized
    # updn  : core-is-spin-polarized
        
    rt = os.environ['WIEN_DMFT_ROOT']
    MPI = ''
    if os.path.exists('mpi_prefix.dat'):
        MPI = open('mpi_prefix.dat', 'r').read().strip()

    fh_info = open('info.out', 'w')
    sys.stdout = fh_info
    
    case=''
    for name in glob.glob('*.struct'):
        case0 = os.path.splitext(name)[0]
        print('case0='+case0)
        if (os.path.isfile(case0+'.clmsum')):
            case = case0
            print('case='+case)
    if not case:
        print('Could not find case, hence exiting')
        sys.exit(1)
        
    cc, so = '', False
    # Checking noncentrosymmetric
    if os.path.isfile(case+".in1c") and os.path.getsize(case+".in1c")>0 :
        print('Found '+case+'.in1c file, hence assuming non-centrosymmetric structure. Switching -c switch!')
        cc = ' -c'
        wopt['cmplx']='c'
    # Check spin-orbit coupling
    if os.path.isfile(case+".inso") and os.path.getsize(case+".inso")>0 :
        print('Found '+case+'.inso file, hence assuming so-coupling exists. Switching -so switch!')
        so = True
        wopt['so']='so'
        
    for i in range(p['Nitt']):
        Step('x lapw0 -f '+case, fh_info)
        if CheckError('lapw0'): break
        Step(rt+'/x_dmft.py lapw1', fh_info)
        if CheckError('lapw1'): break
        if so:
            Step(rt+'/x_dmft.py lapwso', fh_info)
            if CheckError('lapwso'): break
        Step('x lapw2'+cc+' -f '+case, fh_info)
        if CheckError('lapw2'): break
        Step('x lcore -f '+case, fh_info)
        if CheckError('lcore'): break
        
        scf(case, wopt)
        scf1(case)
        
        Step('x mixer -f '+case, fh_info)
        if CheckError('mixer'): break
        
        scfm(case)
        drho,dene = Diff(case)
        succi = (drho<p['cc'] and dene<p['ec'])
        if succi: break
        
