#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import sys, subprocess
import re, os, glob, shutil, socket, time
from os.path import getsize
from scipy import *
from wienfile import Struct
from extractInfo import FindLatticeNd
from scf_merge import scf, scfm
import force_stop_conditon as fstc

utime = '/usr/bin/time'  # This time function is compatible with w2k timing

import utils,indmffile
import convert_processes
import oca, ctqmc

import numpy
nv = map(int,numpy.__version__.split('.'))
if (nv[0],nv[1]) < (1,6):
    loadtxt = io.read_array  # older version of numpy incompatible!
    def savetxt(filename, data):
        io.write_array(filename, data, precision=16)

                
import copy

class Params(object):
    """Class to store various control paramers which control the main flow of the program
       The class is constructed such that it can handle two types of parameters:
         - starting parameters which are read only at the instantiation
         - on fly parameters which can be refreshed (reread from the file) at any time
           by a member function refresh

       Input:
          f_sparams    -- filename of the parameters read only at the beginning
          sparams      -- a dictionary of pairs {name of the variable: default value} which are looked for in the f_sparams file
          f_params     -- filename of the parameters which can be refreshed several times by refresh member function
          params       -- a dictionary of {variable: default value} pairs which can be updated at runtime
       Output:
          None
       Variables which become members of the class are those listed in params and sparams dictionaries
    """
    def __init__(self, f_params, params, f_sparams=None):

        self.p ={}
        self.p.update(params) # default values from the __main__

        # remember params and its data file
        self.f_params = f_params

        # values from the 'params.dat' file
        # executes the parameters at startup
        if f_sparams is not None:
            execfile(f_sparams)
            sp = locals()      # stores parameters
            del sp['self']     # deletes automatic locals
            self.p.update(sp)  # new or updates variables
        
    def refresh(self, fh_info):

        # executes the parameters file
        execfile(self.f_params)
        sp = locals()   # stores locals
        del sp['self']  # deletes automatics

        # prints what was changed
        for c in sp.keys():
            if (c not in self.p.keys()) or (sp[c] != self.p[c]):
                print >>fh_info, '********------ parameter change ------***********'
                print >>fh_info, c, "=", sp[c]

        # updates the dictonary
        self.p.update(sp)
        
    def __getitem__(self, x):
        return self.p[x]

    def __setitem__(self, i, value):
        self.p[i]=value

    def keys(self):
        return self.p.keys()

def dmft1(fday, case, fh_info, extn, MPI, ROOT, m_extn=''):
    inp_name = 'indmf1'+m_extn
    name = 'dmft1'+m_extn
    if not os.path.isfile(case+'.'+inp_name):
        print '  stop error: the required input file '+case+'.'+inp_name+' for the next step could not be found!'
        print >> fday, '  stop error: the required input file '+case+'.'+inp_name+' for the next step could not be found!'
    if not os.path.isfile(name+'.def'):
        print '  stop error: the required input file '+name+'.def for the next step could not be found!'
        print >> fday, '  stop error: the required input file '+name+'.def for the next step could not be found!'
        
    print >> fday, '>%-10s' % name,  '( '+time.strftime("%H:%M:%S")+' )'
    fday.flush()

    fl = open(':log', 'a')
    print >> fl, time.strftime("%a %b %d %H:%M:%S %Z %Y")+'>     '+name
    fl.close()
    
    print >> fh_info, 'Running ---- dmft1 -----'
    cmd = utime+' '+MPI+ ' ./dmft '+name+'.def >> '+name+'_info.out '
    print >> fh_info, '#<'+name+'>: ', cmd
    fh_info.flush()
    subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)
    
    for fe in glob.glob('dmft1.error*'):
        if getsize(fe) !=0:
            print 'ERROR in dmft1 from file:', fe, open(fe,'r').read()
            #sys.exit(1)

    #if m_extn=='':
    #    shutil.copy2(case+'.cdos', case+'.cdos.'+extn)
    #    shutil.copy2(case+'.gc1', case+'.gc1.'+extn)
    #    shutil.copy2(case+'.dlt1', case+'.dlt1.'+extn)

################################################
## The charge self-consistency part
################################################

# Typical Wien2k step
def AStep(fday, case, name, inp_name, WIEN, para, fh_info, opts=''):
    if not os.path.isfile(case+'.'+inp_name) and not os.path.isfile(case+'.'+inp_name+'c'):
        print '  stop error: the required input file '+case+'.'+inp_name+' for the next step could not be found!'
        print >> fday, '  stop error: the required input file '+case+'.'+inp_name+' for the next step could not be found!'
    
    tim = time.strftime("%H:%M:%S")
    print >> fday, '>%-10s ( %s )' % (name+' '+opts, tim)
    #fday.flush()
    
    cmd = WIEN+'/x'+para+' -f '+case+' '+opts+' '+name
    print >> fh_info, ('#<'+name+'>: '), cmd
    fh_info.flush()
    
    info=subprocess.call(cmd,shell=True,stdout=fday)
    
    for fe in glob.glob(name+'.error*'):
        if getsize(fe) !=0:
            print 'ERROR in '+fe+' from file:', open(fe,'r').readlines()

def lapw0(fday, case, WIEN, para, fh_info):
    para=''
    if os.path.isfile(case+'.in0_grr'):
        AStep(fday, case, 'lapw0', 'in0_grr', WIEN, para, fh_info, '-grr')
    AStep(fday, case, 'lapw0', 'in0', WIEN, para, fh_info)

        
def lapw1(fday, case, WIEN, para, dftKS, dmfe, wopt, fh_info):
    pcmplx = '-c' if wopt['cmplx'] else ''
    if dftKS:
        name='lapw1'
        sname = 'lapw1'+wopt['cmplx']
        if not wopt['updn']:   # just one type of core file, so V_{KS} is not polarized
            print >> fday, '>%-10s' % sname, '( '+time.strftime("%H:%M:%S")+' )'
            fl = open(':log', 'a')
            print >> fl, time.strftime("%a %b %d %H:%M:%S %Z %Y")+'>     '+sname
            fl.close()
            cmd = dmfe.ROOT+'/x_dmft.py'+para+' '+pcmplx+' '+name
            print >> fh_info, ('#<'+sname+'>: '), cmd
            fh_info.flush()
            info=subprocess.call(cmd,shell=True,stdout=fh_info)
        else:                  # V_{KS} is polarized, need to lapw1 steps
            print >> fday, '>%-10s' % (sname+' --up'), '( '+time.strftime("%H:%M:%S")+' )'
            fl = open(':log', 'a')
            print >> fl, time.strftime("%a %b %d %H:%M:%S %Z %Y")+'>     '+sname+'up'
            fl.close()
            cmd = dmfe.ROOT+'/x_dmft.py'+para+' '+pcmplx+' --up '+name
            print >> fh_info, ('#<'+sname+'>: '), cmd
            fh_info.flush()
            info=subprocess.call(cmd,shell=True,stdout=fh_info)
            
            print >> fday, '>%-10s' % (sname+' --dn'), '( '+time.strftime("%H:%M:%S")+' )'
            fl = open(':log', 'a')
            print >> fl, time.strftime("%a %b %d %H:%M:%S %Z %Y")+'>     '+sname+'dn'
            fl.close()
            cmd = dmfe.ROOT+'/x_dmft.py'+para+' '+pcmplx+' --dn '+name
            print >> fh_info, ('#<'+sname+'>: '), cmd
            fh_info.flush()
            info=subprocess.call(cmd,shell=True,stdout=fh_info)
    else:
        if not wopt['updn']: # just one type of core file, so V_{KS} is not polarized
            AStep(fday, case, 'lapw1', 'in1', WIEN, para, fh_info, pcmplx)
        else:                # V_{KS} is polarized, need to lapw1 steps
            AStep(fday, case, 'lapw1', 'in1', WIEN, para, fh_info, pcmplx+' -up')
            AStep(fday, case, 'lapw1', 'in1', WIEN, para, fh_info, pcmplx+' -dn')
        
def lapwso(fday, case, WIEN, para, dftKS, dmfe, fh_info):
    if dftKS:
        name='lapwso'
        print >> fday, '>%-10s' % name, '( '+time.strftime("%H:%M:%S")+' )'
        fl = open(':log', 'a')
        print >> fl, time.strftime("%a %b %d %H:%M:%S %Z %Y")+'>     '+name
        fl.close()

        cmd = dmfe.ROOT+'/x_dmft.py'+para+' '+name
        print >> fh_info, ('#<'+name+'>: '), cmd
        fh_info.flush()
        info=subprocess.call(cmd,shell=True,stdout=fh_info)
    else:
        AStep(fday, case, 'lapwso', 'inso', WIEN, para, fh_info)

def lcore(fday, case, WIEN, wopt, fh_info):
    if not wopt['updn']:
        AStep(fday, case, 'lcore', 'inc', WIEN, '', fh_info)
    else:
        if os.path.exists(case+'.inc') : os.remove(case+'.inc')
        os.symlink(case+'.incup', case+'.inc')

        AStep(fday, case, 'lcore', 'inc', WIEN, '', fh_info, ' -up')
        os.remove(case+'.inc')
        os.symlink(case+'.incdn', case+'.inc')
        AStep(fday, case, 'lcore', 'inc', WIEN, '', fh_info, ' -dn')
        
def mixer(fday, case, WIEN, fh_info):
    AStep(fday, case, 'mixer', 'inm', WIEN, '', fh_info)

def dmft2(fday, case, fh_info, MPI, m_extn=''):
    inp_name = 'indmf2'+m_extn
    name = 'dmft2'+m_extn
    if not os.path.isfile(case+'.'+inp_name):
        print '  stop error: the required input file '+case+'.'+inp_name+' for the next step could not be found!'
        print >> fday, '  stop error: the required input file '+case+'.'+inp_name+' for the next step could not be found!'
    if not os.path.isfile(name+'.def'):
        print '  stop error: the required input file '+name+'.def for the next step could not be found!'
        print >> fday, '  stop error: the required input file '+name+'.def for the next step could not be found!'
        
    print >> fday, '>%-10s' % name, '( '+time.strftime("%H:%M:%S")+' )'
    fday.flush()
    
    fl = open(':log', 'a')
    print >> fl, time.strftime("%a %b %d %H:%M:%S %Z %Y")+'>     '+name
    fl.close()
    
    cmd = utime+' '+MPI + ' ./dmft2 '+name+'.def >> dmft2_info.out '
    
    print >> fh_info, ('#<'+name+'>: '), cmd
    fh_info.flush()
    
    info=subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)
    
    for fe in glob.glob(name+'.error*'):
        if getsize(fe) !=0:
            print 'ERROR in '+fe+' from file:', open(fe,'r').readlines()
            #sys.exit(1)

    dat = open(case+'.scf2').readlines()
    dEF=0
    for line in dat:
        if line[:5]==':DEF ':
            dEF = float(line.split()[1])
            break
    
    return dEF

    
def scf1(fday, case, WIEN):
    for i in ['clmsum','vsp','vns','vrespsum','vspdn','vnsdn','vspup','vnsup','clmup','clmdn']:
        name = case+'.'+i
        if os.path.exists(name) and os.path.getsize(name)>0:
            shutil.copy2(name, name+'_old')		#save last cycle
            
        
def Diff(fday, case, dEF):
    #
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
    if size(DChr)>1:
        drho = DChr[-1]
    else:
        drho = 1.0
    if size(DEne)>1:
        dene = abs(DEne[-1]-DEne[-2])
    else:
        dene = 1.0
    if dene==0: dene = 1.0    
    print >> fday, ':ENERGY convergence: ', dene
    print >> fday, ':CHARGE convergence: ', drho
    print >> fday, ':EF     convergence: ', dEF
    return (drho, dene)


SolveImpurity_num_calls=0
Mlist=[]
def SolveImpurity(EF, asolver, iat, extn, p, UpdateAtom, nl_imp, fh_info, fh_pinfo, fday):

    print 'Running ----- impurity solver -----'
    print >> fh_info, 'UpdateAtom=', UpdateAtom
    print >> fh_info, 'Running ----- impurity solver -----'
    fh_info.flush()

    tim = time.strftime("%H:%M:%S")
    print >> fday, '>%-10s ( %s )' % ('impurity', tim)
    fday.flush()
    fl = open(':log', 'a')
    print >> fl, time.strftime("%a %b %d %H:%M:%S %Z %Y")+'>     '+'impurity'
    fl.close()

    iprms = 'iparams'+str(iat)
    ipars = copy.deepcopy(p[iprms])
    
    global SolveImpurity_num_calls
    global Mlist
    if ipars.has_key('Mlist'):
        Mlist=ipars['Mlist'][0]
        ipars.pop('Mlist')
    if Mlist:
        #Mcurrent=ipars['M'][0]
        if len(Mlist)>SolveImpurity_num_calls:
            Mcurrent=Mlist[SolveImpurity_num_calls]
        else:
            Mcurrent=Mlist[-1]
        ipars['M'][0]=Mcurrent
    SolveImpurity_num_calls+=1
    
    
    if p['solver'] in ['OCA','NCA','CTQMC']:
        asolver._exe(ipars, p['DCs'], extn, UpdateAtom, p['wbroad'], p['kbroad'])
    else:
        print 'ERROR: Impurity solver not defined!'
    fh_info.flush()
    
    print >> fh_info, 'Taking care of high frequency'
    ntot = asolver.HighFrequency(ipars, p['DCs'], nl_imp, extn, p['wbroad'], p['kbroad'])
    
    # Reading impurity levels for printing
    fEimp = open('imp.'+str(iat)+'/Eimp.inp', 'r')
    (Eimp,Olap,Edc) = [array([float(x) for x in line.split('#')[0].split()]) for line in fEimp.readlines()]  # throw out comments
    print >> fh_info, 'Eimp=', Eimp-Edc
    print >> fh_info, 'Edc=', Edc
    print >> fh_info, 'ncorr=', ntot

    #icycle, itt = extn.split('.')
    #print >> fh_pinfo, '%3s.%3s %12.6f %12.6f %12.6f %12.6f %12.6f' % (icycle, itt, EF, Eimp[0], Eimp[-1], Edc[0], ntot)
    #fh_pinfo.flush()
    return ntot

def SimplifySiginds(siginds):
    " Takes dictionary of Sigind's and creates list of non-zero columns"

    def union(data):
        " Takes a union of array or list"
        c = []
        for d in data:
            if d not in c:
                c.append(d)
        return c
    
    cols={}
    for icix in siginds.keys():
        Sigind = siginds[icix]
        col = sorted(filter(lambda x: x>0, union(array(Sigind).flatten())))
        #col = sort(union(array(Sigind).flatten())).tolist()
        #if 0 in col: col.remove(0)
        cols[icix] = col
    return cols

def ImpurityLatticeConnection( cols, icols ):
    """Should be done better!
    We have correlated index stored in Sigind and
    impurity index stored in iSigind.
    The translation between the two is necessary when
    the impurity would like to know in which reference frame
    is impurity problem defined. Local rotation should hence
    be passed to impurity, and for that we need the translator.
    """
    #cols = SimplifySiginds(inl.siginds)
    #icols = SimplifySiginds(iSigind)

    print 'cols=', cols
    print 'icols=', icols
    
    icols_ind={}
    for i in icols.keys():
        icols_ind[i]=[]
        for j in cols.keys():
            cimp = sort(icols[i])
            catm = sort(cols[j])
            # equivalent to "catm in cimp" -> if atom j is described by impurity i
            atm_consistent_with_imp = [atm in cimp for atm in catm]
            if len(cimp)>0 and len(catm)>0 and reduce(lambda x,y: x and y, atm_consistent_with_imp) :
                icols_ind[i].append(j)
    return icols_ind

def iSigind_enumerate(iSigind,icols,fh_info):
    # iSigind starts with column number compatible with Sigind from case.indmfl.                                                                                                                                                                                                              
    # For the impurity calculation, we need to start from 1. In the past, we just shifted the index of iSigind so that the first index was 1.                                                                                                                                                 
    # But if there are negative columns (needed for static shift of some orbitals), we had a problem.                                                                                                                                                                                         
    # We recently changed sjoin.py so that columns are just enumerated in their increasing order. Here we need to be compatible with such choice.                                                                                                                                             
    icols_inv = {}
    for i,ic in enumerate(icols):
        icols_inv[ic]=i+1
    print >> fh_info, 'icols_inv=', icols_inv
    #print 'icols_inv=', icols_inv
    iSigind_new = zeros(shape(iSigind), dtype=int)
    for i in range(len(iSigind)):
        for j in range(len(iSigind[i])):
            if iSigind[i][j]>0:
                ii = iSigind[i][j]
                iSigind_new[i,j] = icols_inv[ii]
    return iSigind_new

def YLM_convert_w2k_to_standard(CF, l, qsplit):
    """WIEN2K uses a different convention for spherical harmonics than the ED code
    in the case of d-orbitals.

    Input: T(i,m) = 'Transformation matrix' printed inside case.indmfl
     - i runs over real harmonics (z^2, x^2-y^2, yz, xz, xy)
     - m runs over complex harmonics (-2, -1, 0, 1, 2)

    Output: T(i,m')

    The matrix converting between conventions is
              1  0 0 0 0
              0 -i 0 0 0
    R(m,m') = 0  0 1 0 0
              0  0 0 i 0
              0  0 0 0 1
    This function just does: T(i,m') = T(i,m) * R(m,m')
    """
    if l == 2 and (qsplit not in [1, -1, 4, 11, 13, 14]):
        # only apply conversion if this is NOT a spin-orbit run
        nspin = len(CF)/(2*l+1)  # double size of R if two spins
        R = diag([1, -1j, 1, 1j, 1]*nspin)
        CFnew = dot(CF, R)
    else:
        CFnew = CF

    return CFnew

def combineud(case, ROOT, wopt, fh_info):
    # This averaging is necessary to not double-count the exchange splitting.
    # The exchange splitting is already taken into account in DMFT, and we do
    # not want it also here on DFT level (in the valence charge)
   
    if (False):
        shutil.move(case+'.clmvalup', case+'.clmval_up')
        shutil.move(case+'.clmvaldn', case+'.clmval_dn')
        if os.path.exists(case+'.clmval'): os.remove(case+'.clmval')
        os.symlink(case+'.clmval_up', case+'.clmval')
        os.symlink(case+'.clmval_up', case+'.clmvaldn')

        cmd = ROOT+'/combineud '+case+' 0.5'
        info=subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)
        shutil.move(case+'.clmval_aver', case+'.clmvalup')

        os.remove(case+'.clmval')
        os.remove(case+'.clmvaldn')
        os.symlink(case+'.clmval_dn', case+'.clmval')
        os.symlink(case+'.clmval_dn', case+'.clmvaldn')

        cmd = ROOT+'/combineud '+case+' 0.5'
        info=subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)
        shutil.move(case+'.clmval_aver', case+'.clmvaldn')
    else:
        cmd = ROOT+'/combineud '+case
        if wopt['updn']:
            if wopt['updn']: cmd += ' 0.5'  # We need to take only 1/2 of the average charge as up and 1/2 as down charge
            if os.path.isfile(case+'.clmval'): os.remove(case+'.clmval')
            os.symlink(case+'.clmvalup', case+'.clmval')
            
        print >> fh_info, '#<combineud>: ', cmd
        fh_info.flush()
        
        info=subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)
        
        if not wopt['updn']:
            shutil.move(case+'.clmval',   case+'.clmval_up')
            shutil.move(case+'.clmvaldn', case+'.clmval_dn')
            shutil.move(case+'.clmval_aver', case+'.clmval')
        else:
            shutil.move(case+'.clmvalup', case+'.clmval_up')
            shutil.move(case+'.clmvaldn', case+'.clmval_dn')
            shutil.copy(case+'.clmval_aver', case+'.clmvaldn')
            shutil.move(case+'.clmval_aver', case+'.clmvalup')

def SSplit(ROOT, matsubara, iparams, ntail, fh_info, m_extn=''):
    cmd = ROOT+'/ssplit.py'
    if matsubara : cmd += ' -n '+str(iparams['nom'][0])+' -t '+str(ntail)
    if m_extn: cmd += ' -l '+m_extn
    print >> fh_info, '#<ssplit>: ', cmd

    info=subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)
    
def SGather(ROOT, fh_info, p, m_extn=''):
    cmd = ROOT+'/sgather.py' 
    if p['mix_sigma']<1.0 : 
        cmd += ' -m '+str(p['mix_sigma'])
    if m_extn: cmd += ' -l '+m_extn
    print >> fh_info, '#<sgather>: ', cmd
    info=subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)

def SJoin(ROOT, fh_info, p, m_extn=''):    
    cmd = ROOT+'/sjoin.py -m '+str(p['mix_delta'])
    if m_extn: cmd += ' -l '+m_extn + ' -a '+str(p['EimpAverage'])
    if p['rCF'] is not None: cmd += ' -c "'+ p['rCF']+'"'
    print >> fh_info, '#<sjoin>: ', cmd
    info=subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)

def Prepare_dmft1(dmfe, p, para, wopt, fh_info):
    __updn = '--up'  if wopt['updn'] else ''
    print >> fh_info, '--------- Preparing GF calculation ---------'
    cmd = dmfe.ROOT+'/x_dmft.py -d' + para + ' dmft1 '+__updn+' >> dmft1_info.out 2>&1'
    print >> fh_info, '#<prep-dmft1>: ', cmd
    info=subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)
    if wopt['m_extn']:
        __updn = '--dn'  if wopt['updn'] else ''
        print >> fh_info, '--------- Preparing GC calculation '+wopt['m_extn']+' ---------'
        cmd = dmfe.ROOT+'/x_dmft.py -d' + para + ' dmft1 -l dn '+__updn+' >> dmft1_info.out 2>&1'
        print >> fh_info, '#<prep-dmft1dn>: ', cmd
        info=subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)

def Prepare_dmft2(dmfe, p, para, wopt, fh_info, recomputeEF, mode='c'):
    print >> fh_info, '--------- Preparing Charge calculation ---------'
    __updn = ' --up'  if wopt['updn'] else ''
    cmd = dmfe.ROOT+'/x_dmft.py -d '+'--mode '+mode+' -m '+str(recomputeEF)+' -x '+str(p['mixEF'])+' -w '+str(p['WL'])+__updn+para+' dmft2 >> dmft2_info.out 2>&1'
    print >> fh_info, '#<prep-dmft2>: ', cmd
    info=subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)
    
    if wopt['m_extn']:
        print >> fh_info, '--------- Preparing CHR calculation '+wopt['m_extn']+'---------'
        __updn = ' --dn'  if wopt['updn'] else ''
        cmd = dmfe.ROOT+'/x_dmft.py -d '+'--mode '+mode+' -m 0'+__updn+para+' dmft2 -l dn >> dmft2_info.out 2>&1'
        print >> fh_info, '#<prep-dmft2>: ', cmd
        info=subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)

def CreateEorb(case, wopt, siginds, icols_ind):
    fsg = open('sig.inp','r')
    exec(fsg.next()[1:].lstrip())  # s_oo from sig.inp                                                                                                                   $
    exec(fsg.next()[1:].lstrip())  # Edc  from sig.inp                                                                                                                   $

    up__ = 'up' if wopt['updn'] else ''
    fsc = open(case+'.scf2'+up__,'r')
    dat = fsc.readlines()

    nd=[]
    dim=[]
    for l in range(len(dat)):
        if dat[l][:5]==':NCOR':
            dd = dat[l].split()
            nd.append( float(dd[1]) )
            dim.append( int(dd[3]) )
    if wopt['m_extn']:
        fsc = open(case+'.scf2dn','r')
        dat = fsc.readlines()
        nd_=[]
        dim_=[]
        for l in range(len(dat)):
            if dat[l][:5]==':NCOR':
                dd = dat[l].split()
                nd_.append( float(dd[1]) )
                dim_.append( int(dd[3]) )
        nd = 0.5*(array(nd)+array(nd_))
    return (nd,Edc[0])

def ConvertFromLatticeToImp(nd, icols_ind):
    nc=[]
    for i in icols_ind.keys():
        ndc=[ nd[a-1] for a in icols_ind[i] ]
        nc.append( sum(ndc)/len(ndc) )
    return nc

def FindEtot(case,wopt):
    def FindLine(dat, strng, item):
        Val = 0.0
        for line in dat:
            if re.match(strng, line) is not None:
                Val = float(line.split()[item])
                break
        return Val

    fs = open(case+'.scfm','r')
    dat = fs.readlines()[::-1]
    fs.close()
    ETOT = FindLine(dat, ':ENE', 8)

    if os.path.isfile('Eorb_imp.dat'):
        fs=open('Eorb_imp.dat','r')
        dat = fs.readlines()[::-1]
        fs.close()
        XEORB = FindLine(dat, ':XEORB',1)
        EORB  = FindLine(dat, ':EORB',1)
        ZEORB = FindLine(dat, ':ZEORB',1)
    else:
        XEORB = EORB = ZEORB = 0.
        

    up__ = 'up' if wopt['updn'] else ''
    fs = open(case+'.scf2'+up__,'r')
    dat = fs.readlines()[::-1]
    fs.close()
    SUM = FindLine(dat, ':SUM ', 6)
    XSUM = FindLine(dat, ':XSUM', 6)
    YSUM = FindLine(dat, ':YSUM', 6)
    ZSUM = FindLine(dat, ':ZSUM', 6)
    if wopt['m_extn']:
        fsn = open(case+'.scf2'+wopt['m_extn'],'r')
        datn = fsn.readlines()[::-1]
        fsn.close()
        SUM_ = FindLine(datn, ':SUM ', 6)
        XSUM_ = FindLine(datn, ':XSUM', 6)
        YSUM_ = FindLine(datn, ':YSUM', 6)
        ZSUM_ = FindLine(datn, ':ZSUM', 6)
        SUM = 0.5*(SUM+SUM_)
        XSUM = 0.5*(XSUM+XSUM_)
        YSUM = 0.5*(YSUM+YSUM_)
        ZSUM = 0.5*(ZSUM+ZSUM_)
    return (ETOT,SUM,XSUM,YSUM,ZSUM,EORB,XEORB,ZEORB)


def Connect_ImpurityProjector(icols_ind,iSigind,struct):
    "Find connection between the i-th impurity and j-th projector from projector.dat"
    #atms = sum(icols_ind.values()) # all correlated atoms
    atms = reduce(lambda x,y: x+y, icols_ind.values())
    impurity_projector=-ones(len(iSigind.keys()),dtype=int)  # the connection between the impurity and projector.dat
    latom=0
    for jatom in range(struct.nat):
        if (latom+1) in atms:           # This atom listed in indmfl file                                                                                                                                   
            impurty=0
            for imp in iSigind.keys():
                if (latom+1) in icols_ind[imp]: # This atom (jatom) is connected to this impurity (imp)
                    impurty=imp
                    break
            if impurity_projector[imp]<0: impurity_projector[imp]=jatom
        latom += struct.mult[jatom]
    return impurity_projector.tolist()

def Update_Vdc(icols,solver):
    # Finds current n_d on the atom from dmft2 calculation
    # Now transforms to format in sig.inp
    allcols = sort( reduce(lambda x,y: x+y, icols.values()) )
    noccur = zeros(max(allcols),dtype=int)
    for c in allcols: noccur[c-1]+=1
    rVdc=zeros( max(allcols) ,dtype=float)
    for imp in icols.keys():  # over all impurity problems
        Vdc = solver[imp].Lattice_Vdc(atleast_1d(nl_imp[imp]),p['iparams'+str(imp)],p['DCs'])
        mincol = min(icols[imp])
        for c in icols[imp]:
            rVdc[c-1]  += Vdc[c-mincol]*(1./noccur[c-1])
    #savetxt('Edc.dat',rVdc)
    fin = open('sig.inp','r')
    dat = fin.readlines()
    fin.close()
    dat[1]='# Edc= '+str(list(rVdc))+'\n'
    fot = open('sig.inp','w')
    fot.writelines(dat)
    fot.close()

def FermiF(dmfe, w2k, p, para, wopt, fh_info, fday, matsubara):
    # In case of ferromagnetic or ferrimagnetic metal, we need to compute common EF for up and dn spins
    case = w2k.case
    scratch = w2k.SCRATCH
    Prepare_dmft2(dmfe, p, para, wopt, fh_info, 1, 'e')
    SSplit(dmfe.ROOT, matsubara, p['iparams0'], p['ntail'], fh_info)
    dmft2(fday,case,fh_info,dmfe.MPI2)
    SSplit(dmfe.ROOT, matsubara, p['iparams0'], p['ntail'], fh_info, wopt['m_extn'])
    dmft2(fday,case,fh_info,dmfe.MPI2,wopt['m_extn'])
    cmd = dmfe.ROOT+'/fermif '+scratch+'/'+w2k.case+'.bnds '+scratch+'/'+w2k.case+'.bnds'+wopt['m_extn']
    print >> fh_info, ('#<fermif>: '), cmd
    fh_info.flush()
    tim = time.strftime("%H:%M:%S")
    print >> fday, '>%-10s ( %s )' % ('fermif', tim)
    fday.flush()
    fl = open(':log', 'a')
    print >> fl, time.strftime("%a %b %d %H:%M:%S %Z %Y")+'>     '+'fermif'
    fl.close()
    info=subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)
    if info!=0:
        print 'Problems evaluating fermif!'
        
if __name__ == '__main__':
    # -------------- Parameters which are set through input files sparams.dat and params.dat -----------
    # Default values of some parameters
    params = {'runIMP'        : True,      # Run impurity solver
              'runGF'         : True,      # Run dmft1 for G and Delta
              'UpdateAtom'    : 1,         # recompute cix-file for the impurity solver
              'DCs'           : 'nominal',    # the double counting scheme, which fixes Edc through n0
              #'DCs'          : 'FLL', # the default double counting scheme
              #'DCs'          : 'AMF',     # experimental
              #'DCs'          : 'fixn',    # old name for nominal DC
              'max_dmft_iterations': 1,    # number of iteration of the dmft-loop only
              'max_lda_iterations' : 1,    # number of iteration of the LDA-loop only
              'max_lda_iterations_optimize' : 1000,    # number of iteration of the LDA-loop when structure should be optimized
              'finish'        : 50,        # number of iterations of the charge loop
              'da'            : 0.0,       # experimental
              'wbroad'        : 0.0,       # broadening of sigma on the imaginary axis
              'kbroad'        : 0.0,       # broadening of sigma on the imaginary axis
              'solver'        : 'CTQMC',   # CTQMC solver
              #'solver'       : 'OCA',     # OCA solver
              'ntail'         : 200,       # on imaginary axis, number of points in the tail of the logarithmic mesh
              'mix_delta'     : 1.0,       # whether to mix delta, or not
              'mix_sigma'     : 1.0,       # whether to mix sigma, or not
              'riter'         : 100,        # How often to restart broyden for charge mixing
              'sleeptime'     : 2,         # If broyden file are present at the submit time, user has some time to kill the job
              'cc'            : 1e-5,      # the charge density precision to stop the LDA+DMFT run
              'ec'            : 1e-5,      # the energy precision to stop the LDA+DMFT run
              #'so'            : False,     # spin-orbit coupling
              #'c'             : False,     # non-centrosymmetric
              'rCF'           : None,      # Reduction of the crystal field splitting, if necessary
              'recomputeEF'   : 1,         # Recompute EF in dmft2 step. If recomputeEF=2, it tries to find an insulating gap.
              'mixEF'         : 1.0,       # Chemical potential can be mixed with mixEF<1.0
              'saver'         : 0.0,       # Average over the last few DMFT self-energies is used for charge. Should be between [0,1]
              'WL'            : 1.0,       # If recomputeEF=2, we look for the pole in the interval [-WL,WL]
              'EimpAverage'   : 1,         # Weather to average impurity levels over up/dn
              'dftKS'         : True,      # Use internal capabilities to solve the Kohn-Sham problem
              'remove_broyd'  : 1000,      # How many dft steps before broyden files are removed
              'Ferro'         : False,     # For ferromagnetic (or ferrimagnetic) metals it should be set to True, so that the chemical potential is properly computed
              'force_Ene_diff': 0.001,     # When checking if force is converged, we require Energy to be converged as well up to this difference
              'force_Rho_diff': 0.01,      # When checking if force is converged, we require Charge to be converged as well up to this difference
              }

    Ry2eV = 13.60569253
    p = Params('params.dat', params)

    if len(sys.argv)>1 and sys.argv[1] in ['-h', '--help']:
        help="""
        The script runs LDA+DMFT. It is a wrapper, which calls other scripts and executables.
        It executes Wien2K in parallel if .machines file exists. Otherwise only DMFT part will run in parallel.
        Usuall execution goes through the following steps in a loop:
        
          x lapw0        -- computes LDA potential on current LDA+DMFT charge
          x lapw1        -- solves LDA eigenvalue problem
          [x lapwso]     -- adds spin-orbit
          [x dmft0]      -- computes dmft eigenvalues
          [x mu]         -- computes dmft chemical potential
          x dmft1        -- computes local green's function and hybridization
          run impurity   -- runs impurity problem
          x dmft2        -- computes LDA+DMFT valence charge, and the chemical potential
          x lcore        -- computes LDA core charge
          x mixer        -- mixes total charge with the previous charge

        The most common parameters, which should be given through 'params.dat' file, include:

        name       possible values   default     help
        --------------------------------------------------------------------------------------
        solver         [OCA | CTQMC] OCA      # impurity solver
        max_iterations [int]         1        # number of iteration of the dmft-loop only
        finish         [int]         100      # number of iterations of the charge+dmft loop
        cc             [float]       1e-5,    # the charge density precision to stop the LDA+DMFT run
        ec             [float]       1e-5,    # the energy precision to stop the LDA+DMFT run
        broyd          [True| False] True     # Are we using broyden for charge mixing
        riter          [int]         99       # How often to restart broyden for charge mixing
        sleeptime      [int]         2        # If broyden file are present at the submit time, user
                                                # has some time to kill the job
        so             [True| False] False    # spin-orbit coupling
        DCs            [fixn | nominal | FLL | AMF | exact | exact_l ]
                                     fixn     # double-counting scheme
        runIMP         [True| False] True     # Run impurity solver or skip this step
        runGF          [True| False] True     # Run dmft1 for G and Delta, or skip this step
        ntail          [int]         30       # sigma.inp on imaginary axis, number of points in the
                                                # tail of the logarithmic mesh
        mix_delta      [float]       1.0      # linear mixing parameter for Delta.
        UpdateAtom     [int]         0        # How often to recompute cix-file inside the
                                              # impurity solver
        wbroad         [float]       0.0      # constant broadening of sigma on the imaginary axis
        kbroad         [float]       0.0,     # linear broadening of sigma on the imaginary axis
        recomputeEF    [0|1|2]       1        # If we want to recompute the chemical potential in dmft2 step. 2->for Mott insulator
        com            [0| 1]        0        # when computing chemical potential in [dmft0|mu] step
                                                # on imaginary axis, we can subtract Sigma_oo (com=0)
                                                # or Sigma(0) (com=1)
        mixEF          0.0-1.0                # mixing of the chemical potential
        da             [float]       0.0      # experimental - connected with DC scheme
       
        """
        print help
        sys.exit(0)
    
    # list of options for wien2k
    wopt={'cmplx':'', 'so':'', 'm_extn':'', 'updn':''}  
    # cmplx : no-inversion, 
    # so    : spin-orbit coupling present
    # m_extn: dmft is spin polarized
    # updn  : core-is-spin-polarized

    dmfe = utils.DmftEnvironment()  # DMFT paths
    w2k = utils.W2kEnvironment()    # W2k filenames and paths

    para = ''
    if p['dftKS'] and dmfe.MPI2:
        para = ' -p'  # Parallel run with internal lapw1. We could use either '-p' to have many vector files, or just '', to have single vector file. Not clear what is faster...
    if os.path.isfile('.machines') and os.path.getsize('.machines')>0 : 
        para = ' -p'              # Using w2k parallelization. For not neccessary for SO.
        p['dftKS'] = False        # Switching off internal lapw1

    # Processing 'case.indmfl' file
    inl = indmffile.Indmfl(w2k.case) # case.indmfl file
    inl.read()                       # case.indmfl read

    toclean = glob.glob('*.scf')+glob.glob('Edc.dat')+glob.glob('dmft1.error*')+glob.glob('dmft2.error*')+glob.glob('*.cdos.*')+glob.glob('*.dlt.*')+glob.glob('*.gc1.*')+glob.glob('*.outputdmf1.*')+glob.glob('*.outputdmf2.*')
    for f in toclean: os.remove(f)
    
    # Info files
    fday = open(w2k.case+'.dayfile', 'w')
    print >> fday, 'Calculating '+w2k.case+' in '+os.getcwd()+'\n'+'on '+socket.gethostname()+' with PID '+str(os.getpid())
    print >> fday, '\n\n'
    
    fh_info = open('dmft_info.out','w')
    fh_pinfo = open('info.iterate', 'a')
    
    #if len(glob.glob('dmft0_info.out'))!=0 : os.remove('dmft0_info.out')
    if len(glob.glob('dmft1_info.out'))!=0 : os.remove('dmft1_info.out')
    if len(glob.glob('dmft2_info.out'))!=0 : os.remove('dmft2_info.out')
    
    inldn = None
    if os.path.exists(w2k.case+'.indmfl'+'dn'): wopt['m_extn'] = 'dn'
    if wopt['m_extn']:
        print >> fh_info, 'INFO: case.indmfldn present => magnetic calculation with two dmft2 steps'
        inldn = indmffile.Indmfl(w2k.case, 'indmfl'+wopt['m_extn'])
        inldn.read()

    # Reading parameters from params.dat
    p.refresh(fh_info)

    # Reads structure file
    struct = Struct(w2k.case, fh_info)

    # Check spin-orbit coupling
    if os.path.isfile(w2k.case+".inso") and os.path.getsize(w2k.case+".inso")>0 :
        print >> fh_info, 'Found '+w2k.case+'.inso file, hence assuming so-coupling exists. Switching -so switch!'
        wopt['so']='so'
    if os.path.isfile(w2k.case+".in1c") and os.path.getsize(w2k.case+".in1c")>0 :
        print >> fh_info, 'Found '+w2k.case+'.in1c file, hence assuming non-centrosymmetric structure. Switching -c switch!'
        wopt['cmplx']='c'
        
    # corelated indexes
    cixs = inl.siginds.keys()        # all columns for cix

    # Nuclear charges of all correlated blocks
    Znuc={} 
    for icix in cixs:
        atm = inl.cps[icix][0][0]                # correlated atoms
        Znuc[icix] = int(struct.aZnuc[atm-1]+0.5)  # nuclear charge of each correlated block
        
    print >> fh_info, 'Znucs=', Znuc
    
    # Independent impurity problems
    iSigind = utils.ParsIndmfi(w2k.case)

    print >> fh_info, 'iSiginds=', iSigind
    
    cols = SimplifySiginds(inl.siginds)
    icols = SimplifySiginds(iSigind)
    icols_ind = ImpurityLatticeConnection( cols, icols )
    print >> fh_info, 'Impurity-lattice connection: icols_ind=', icols_ind
    
    fh_pinfo.write( '%3s %3s.%4s %12s %12s ' % ('#', '#','#','mu','Vdc') )
    fh_pinfo.write( '%15s %15s %15s ' % ('Etot', 'Ftot+T*Simp', 'Ftot+T*Simp') )
    for iat in iSigind.keys():   # Over all inequivalent impurity problems
        fh_pinfo.write( '%12s %12s %12s %12s ' % ('n_latt', 'n_imp','Eimp[0]','Eimp[-1]') )
    fh_pinfo.write('\n')
    fh_pinfo.flush()


    if p['DCs'][:5]=='exact':
        Nmax=[]
        for i in iSigind.keys():
            ic = icols_ind[i][0]
            iatom, l, qsplit = inl.cps[ic][0]  # Doesn't take into account cluster problems correctly
            Nmax.append(2*(2*l+1))

        # Finds connection between impurity and projector
        impurity_projector = Connect_ImpurityProjector(icols_ind,iSigind,struct)
        print >> fh_info, 'connection -- impurity_projector=', impurity_projector
        # Finds the screening length and corresponding J's given imput value of Coulomb U for al impurity problems
        UJ_icase={}
        for imp,icase in enumerate(impurity_projector):
            iprms = 'iparams'+str(imp)
            U=p[iprms]['U'][0]            # we always need U to determine screening length
            JH=0
            if p[iprms].has_key('J'):    # if Hunds J is awailable, we take it into account
                JH = p[iprms]['J'][0]
            lmbda=1e-6
            if p[iprms].has_key('DC_lmbda'): # if user wants to use finite lambda with exactd DC, it should set iparamsx['DC_lmbda']
                lmbda = p[iprms]['DC_lmbda'][0]
            nfraction=1.0
            if p[iprms].has_key('nf0'):
                n0 = p[iprms]['nf0'][0]
                nfraction = (Nmax[imp]-n0)/(Nmax[imp]+0.0)
            UJ_icase[icase]=[U,JH,lmbda,nfraction]       # Saves U into U_icase, which is needed below to compute lambda
        
        print 'UJ_icase=', UJ_icase
        
        import RCoulombU
        if len(p['DCs'])>5 and p['DCs'][:6]=='exactd':      # epsilon screening
            print 'Have exactd!'
            UlamJ = RCoulombU.GetDielectricFunctions(UJ_icase)  # Computes dielectric function epsilon and sets yukawa lambda=1e-6
        elif len(p['DCs'])>5 and p['DCs'][:6]=='exacty':    # yukawa screening
            for icase in UJ_icase.keys():
                UJ_icase[icase][3]=1.0 # all nfractions equal to 1.0, hence yukawa
            print 'UJ_icase=', UJ_icase
            UlamJ = RCoulombU.GetLambdas(UJ_icase) # Computes yukawa lambda and dielectric function so that they match user defined U & J.
        else :  # mixure between exactd and exacty
            print 'UJ_icase=', UJ_icase
            UlamJ = RCoulombU.GetLambdas(UJ_icase)
            
        for imp,icase in enumerate(impurity_projector): # Saves into impurity
            iprms = 'iparams'+str(imp)
            p[iprms]['icase'] = icase # Connection between impurity and projector in projector.dat'
            p[iprms]['UlamJ'] = UlamJ[icase]
        
        print >> fh_info, 'UlamJ = {iatom : (U,lambda,epsilon,[J2,J4,...]), ....}'
        print >> fh_info, 'UlamJ=', UlamJ
        
    # Solvers initialized
    if p['solver'] in ['OCA', 'NCA']:
        if inl.matsubara:
            print 'ERROR: Can not run OCA/NCA on imaginary axis! Change case.indmfl:matsubara or params.dat:solver'
        solver=[]    
        for i in iSigind.keys():
            ic = icols_ind[i][0]
            iatom, l, qsplit = inl.cps[ic][0]  # Doesn't take into account cluster problems correctly
            cftrans = inl.cftrans[ic]
            iSigind_ = iSigind_enumerate(iSigind[i],icols[i],fh_info)
            if shape(iSigind[i])[0]==2*shape(cftrans)[0] and shape(iSigind[i])[1]==2*shape(cftrans)[1]:
                zr = zeros(shape(cftrans))
                cftrans = vstack((hstack((cftrans,zr)),hstack((zr,cftrans))))
            solver.append( oca.IMP_OCA(Znuc[ic], 'imp.'+str(i), p['iparams'+str(i)], iSigind_, cftrans) )
    elif p['solver'] in ['CTQMC']:
        if not inl.matsubara:
            print 'ERROR: Can not run CTQMC on real axis! Change case.indmfl:matsubara or params.dat:solver'
        solver=[]
        for i in iSigind.keys():
            if len(icols_ind[i])==0: continue  # This atom is treated as open-core
            ic = icols_ind[i][0]
            iatom, l, qsplit = inl.cps[ic][0]  # Doesn't take into account cluster problems correctly
            cftrans = inl.cftrans[ic]
            iSigind_ = iSigind_enumerate(iSigind[i],icols[i],fh_info)
            if shape(iSigind[i])[0]==2*shape(cftrans)[0] and shape(iSigind[i])[1]==2*shape(cftrans)[1]:
                zr = zeros(shape(cftrans))
                cftrans = vstack((hstack((cftrans,zr)),hstack((zr,cftrans))))
            solver.append( ctqmc.IMP_CTQMC(Znuc[ic], 'imp.'+str(i), p['iparams'+str(i)], iSigind_, cftrans) )
    else:
        print 'ERROR: No solver specified! Bolling out.'

    shutil.copy2(dmfe.ROOT+'/dmft', '.')  # For parallel execution, executable has to be on the current directory
    shutil.copy2(dmfe.ROOT+'/dmft2', '.') # For parallel execution, executable has to be on the current directory

    if os.path.isfile(w2k.case+".incup") and os.path.isfile(w2k.case+".incdn") and os.path.getsize(w2k.case+".incup")>0 and os.path.getsize(w2k.case+".incdn")>0:
        wopt['updn']=True # Computing J_Hunds by constrained-DMFT, hence need different core density for up and dn spin.
        print >> fh_info, 'Found '+w2k.case+'.incup file, hence assuming core is spin polarized, and switching updn=up/dn on'
        if not os.path.isfile(w2k.case+'.clmup') or os.path.getsize(w2k.case+'.clmup')==0:
            print >> fh_info, 'ERROR: Missing '+w2k.case+'.clmup /clmdn files, which are needed when '+w2k.case+'.incup /incdn files are present'
            sys.exit(1)
            
    Prepare_dmft1(dmfe, p, para, wopt, fh_info)
    Prepare_dmft2(dmfe, p, para, wopt, fh_info, p['recomputeEF'])
    
    shutil.copy2('sig.inp', 'sig.inp.0')
    
    if p['finish']>1: # Charge self-consistent
        
        if os.path.isfile(w2k.case+'.clmsum_old'): shutil.copy2(w2k.case+'.clmsum_old', w2k.case+'.clmsum')
        if not os.path.isfile(w2k.case+'.clmsum'):
            print 'no '+w2k.case+'.clmsum(_old) file found, which is necessary for lapw0 !'
            print >> fday, 'no '+w2k.case+'.clmsum(_old) file found, which is necessary for lapw0 !'
            sys.exit(1)
        
        
        if os.path.isfile(w2k.case+'.broyd1'):
            print w2k.case+'.broyd* files present! I will remove previous broyd files.' 
            # print 'You have '+str(p['sleeptime'])+' seconds to kill this job ( ^C   or   kill '+str(os.getpid())+' )'
            # print 'or the script will rm *.broyd* and continue (use -NI to avoid automatic rm)'
            # time.sleep(p['sleeptime'])
            cmd = 'rm -f *.broyd*'
            subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)
            print >> fday, w2k.case+'.broyd* files removed!'
        
        print >> fday, '\n   start'+(' '*8)+ time.asctime()+' with lapw0 ('+str(1)+'/'+str(p['riter'])+' to go)'
        
        # Starting with LAPW
        print >> fday, '\n   cycle %s \t%s %s/%s to go\n' % (0,time.asctime(),p['finish'],(p['finish'])%p['riter'])
        lapw0(fday,w2k.case,w2k.WIENROOT,para, fh_info); fday.flush()
        lapw1(fday,w2k.case,w2k.WIENROOT,para,p['dftKS'],dmfe,wopt,fh_info); fday.flush()
        if wopt['so']: lapwso(fday,w2k.case,w2k.WIENROOT,para,p['dftKS'],dmfe, fh_info)
        fday.flush()
        

    if para and (not p['dftKS']): convert_processes.convert(wopt['so'], fh_info)
    
    #####################
    # Major charge loop #
    #####################
    icycle = nimp = 0
    istep0=0
    nl_imp = None
    (RelaxingStructure, mix_method) = fstc.AreWeRelaxingStructure2(w2k.case)
    print >> fh_info, 'Relaxing Structure=', RelaxingStructure
    
    while icycle < p['finish']:

        if (icycle>0): # Does not exists yet. Will find it after dmft1
            (nl_imp,nd_imp) = FindLatticeNd(w2k.case, wopt['m_extn'], inl.siginds, icols_ind, fh_info, scf_name='scf2')
        
        #####################
        # DMFT loop         #
        #####################
        sigmas=[]
        itt = 0
        n_imp = zeros(len(iSigind.keys()))
        
        if p.p.has_key('GoodGuess') and p['GoodGuess']==True:
            if icycle == 0 : 
                p['max_dmft_iterations'] = 0 # do not run impurity at the beginning, but rather converge DMFT charge on existing self-energy
                itt = 1
            if icycle == 1 : 
                p['max_dmft_iterations'] = 2 # at the next iteration, rerun impurity twice to have good histogram from the second run
                itt = 0
            if icycle > 1  : p['max_dmft_iterations'] = 1 # finally, continue as normally
            print >> fh_info, 'We have GoodGuess and max_dmft_iterations=', p['max_dmft_iterations'], 'itt=', itt
        
        while itt < p['max_dmft_iterations']:
            
            extn = str(icycle+1)+'.'+str(itt+1)  # Extension of output files
            p.refresh(fh_info)
            
            # Recomputing the chemical potential, but only if more than one DMFT step. Otherwise EF was recomputed in charge-SC step.
            if itt>1 and p['recomputeEF']:
                if p['Ferro'] and wopt['m_extn']:
                    # In case of ferromagnetic or ferrimagnetic metal, we need to compute common EF for up and dn spins
                    FermiF(dmfe, w2k, p, para, wopt, fh_info, fday, inl.matsubara)
                else:
                    Prepare_dmft2(dmfe, p, para, wopt, fh_info, 1, 'm')
                    SSplit(dmfe.ROOT, inl.matsubara, p['iparams0'], p['ntail'], fh_info)
                    dEF =  dmft2(fday, w2k.case, fh_info, dmfe.MPI2); fday.flush()

            # Computing Green's function
            if p['runGF']:
                SSplit(dmfe.ROOT, inl.matsubara, p['iparams0'], p['ntail'], fh_info)
                dmft1(fday, w2k.case, fh_info, extn, dmfe.MPI2, dmfe.ROOT)
                if wopt['m_extn']:
                    SSplit(dmfe.ROOT, inl.matsubara, p['iparams0'], p['ntail'], fh_info, wopt['m_extn'])
                    dmft1(fday, w2k.case, fh_info, extn, dmfe.MPI2, dmfe.ROOT, wopt['m_extn'])
                (nl_imp,nd_imp) = FindLatticeNd(w2k.case, wopt['m_extn'], inl.siginds, icols_ind, fh_info, scf_name='scfdmft1',cmplx=True)
            
            SJoin(dmfe.ROOT, fh_info, p, wopt['m_extn'])
            
            if os.path.isfile('EF.dat'):
                EF = float(open('EF.dat').read())
            else:
                open(w2k.case+'.indmf1')
                EF =float(open(w2k.case+'.indmf1').readlines()[1].split()[1])
            
            # ------------ Impurity solver -------------------
            if p['runIMP']:
                nimp += 1
                UpdateAtom = (p['UpdateAtom']!=0 and mod(nimp,p['UpdateAtom'])==0)
                
                for iat in iSigind.keys():   # Over all inequivalent impurity problems
                    ni=SolveImpurity(EF, solver[iat], iat, extn, p, UpdateAtom, nl_imp[iat], fh_info, fh_pinfo, fday)
                    n_imp[iat] += ni
                    
                SGather(dmfe.ROOT, fh_info, p, wopt['m_extn'])
                
                IEorb = 0.0
                XEorb = 0.0
                ZEorb = 0.0
                dat=[]
                for iat in iSigind.keys():   # Over all inequivalent impurity problems
                    if os.path.exists(solver[iat].dir+'/Eorb.dat'):
                        dat += ['--- impurity '+str(iat)+' with degeneracy '+str(len(icols_ind[iat]))+' ---\n']
                        datw = open(solver[iat].dir+'/Eorb.dat','r').readlines()
                        for line in datw:
                            if line[:6]==':IEORB':
                                IEorb += float(line[6:]) * len(icols_ind[iat]) # Takes into account the degeneracy of the impurity
                            elif line[:6]==':XEORB':
                                XEorb += float(line[6:]) * len(icols_ind[iat]) # Takes into account the degeneracy of the impurity  
                            elif line[:6]==':ZEORB':
                                ZEorb += float(line[6:]) * len(icols_ind[iat]) # Takes into account the degeneracy of the impurity  
                            else:
                                dat += [line]
                        
                dat += [':EORB   '+ str(IEorb)+'\n']
                dat += [':XEORB  '+ str(XEorb)+'\n']
                dat += [':ZEORB  '+ str(ZEorb)+'\n']
                fiorb=open('Eorb_imp.dat','w')
                fiorb.writelines(dat)  # Saves on the current directory the combination of all impurity problems
                fiorb.close()
                
                shutil.copy2('sig.inp','sig.inp.'+extn)
                sigmas.append('sig.inp.'+extn) # save the filename
                
            print >> fh_info, '------- DMFT part of the step # ', itt, 'done ---------\n'
            fh_info.flush()

            itt += 1   # END DMFT INNER LOOP
            
        n_imp *= 1./itt
        
        if int(len(sigmas)*p['saver'])>1 :
            sigmas.reverse()
            sigmas = sigmas[0:max(1,int(len(sigmas)*p['saver']))]
            print  >> fh_info, 'Averaging over the following self-energies', sigmas
            cmd = 'saverage.py -o sig.inp ' + (' '.join(sigmas))
            print >> fh_info, '#<saverage>: ', cmd
            subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)

        (RelaxingStructure, mix_method) = fstc.AreWeRelaxingStructure2(w2k.case)
        if RelaxingStructure:
            fstc.TOTtoFOR(w2k.case)
            p['max_lda_iterations']=p['max_lda_iterations_optimize']
            print >> fh_info, 'Mixing method corresponds to structural optimization.'
            print >> fh_info, ' changing TOT to FOR in case.in2 files'
            print >> fh_info, ' increasing max_lda_iterations to ', p['max_lda_iterations']
        

        for ita in range(p['max_lda_iterations']):
            if ita%2 == 1 : p.refresh(fh_info)

            recomputeEF = p['recomputeEF']

            if p['Ferro'] and wopt['m_extn'] and recomputeEF:
                # In case of ferromagnetic or ferrimagnetic metal, we need to compute common EF for up and dn spins
                FermiF(dmfe, w2k, p, para, wopt, fh_info, fday, inl.matsubara)
                recomputeEF = 0
            
            Prepare_dmft2(dmfe, p, para, wopt, fh_info, recomputeEF, 'c')
            SSplit(dmfe.ROOT, inl.matsubara, p['iparams0'], p['ntail'], fh_info)
            dEF = dmft2(fday,w2k.case,fh_info,dmfe.MPI2); fday.flush()
            if wopt['m_extn']:
                SSplit(dmfe.ROOT, inl.matsubara, p['iparams0'], p['ntail'], fh_info, wopt['m_extn'])
                dmft2(fday,w2k.case,fh_info,dmfe.MPI2,wopt['m_extn'])
                fday.flush()
                combineud(w2k.case, dmfe.ROOT, wopt, fh_info)
            
            nd,Vdc = CreateEorb(w2k.case, wopt, inl.siginds, icols_ind)

            nc = ConvertFromLatticeToImp(nd, icols_ind)

            lcore(fday,w2k.case,w2k.WIENROOT,wopt,fh_info);       fday.flush()
            scf  (w2k.case, wopt)
            scf1 (fday,w2k.case,w2k.WIENROOT);
            mixer(fday,w2k.case,w2k.WIENROOT,fh_info);       fday.flush()
            scfm (fday,w2k.case,w2k.WIENROOT);               fday.flush()
            drho, dene = Diff(fday, w2k.case, dEF);          fday.flush()
            
            # Writting to info.iterate file
            if os.path.isfile('EF.dat'): EF = float(open('EF.dat').read())
            (ETOT,SUM,XSUM,YSUM,ZSUM,EORB,XEORB,ZEORB) = FindEtot(w2k.case,wopt)
            
            fh_pinfo.write( '%3d %3s.%4s %12.6f %12.6f %15.6f %15.6f %15.6f ' % (istep0, icycle, ita, EF, Vdc, ETOT, ETOT-SUM+ZSUM-EORB+ZEORB, ETOT-SUM+YSUM-EORB+XEORB) )
            for iat in iSigind.keys():   # Over all inequivalent impurity problems
                if os.path.exists('imp.'+str(iat)+'/Eimp.inp'):
                    fEimp = open('imp.'+str(iat)+'/Eimp.inp', 'r')
                    (Eimp,Olap,Edc) = [array([float(x) for x in line.split('#')[0].split()]) for line in fEimp.readlines()]  # throw out comments
                else: # If runGF=False, this is not available
                    Eimp=array([0,0])
                fh_pinfo.write( '%12.6f %12.6f %12.6f %12.6f ' % (nc[iat], n_imp[iat], Eimp[0], Eimp[-1]) )
            fh_pinfo.write('\n')
            fh_pinfo.flush()
            istep0 += 1
            print >> fh_info, '------- LDA(charge) step', ita, icycle, 'done ---------'
            
            istep1 = icycle+1
            istep2 = ita+1

            if fstc.AreWeRelaxingStructure(w2k.case):
                succi = False # Do not stop DFT iterations when relaxing structure
                if fstc.AreForcesConverged(w2k.case, how_many_times=3, min_Ene_diff=p['force_Ene_diff'], min_Rho_diff=p['force_Rho_diff'], info=fh_info):
                    fstc.StopStructureOptimization(w2k.case,info=fh_info)
            else:
                succi = (drho<p['cc'] and dene<p['ec'])
                
            if succi or ita==(p['max_lda_iterations']-1):
                istep1 += 1
                istep2 = 1
            
            print >> fday, '\n:ITT%-2d   cycle %s \t%s %s/%s to go\n' % (istep2, istep1,time.asctime(),p['finish']-istep1,(p['finish']-istep1)%p['riter'])
            
            lapw0(fday,w2k.case,w2k.WIENROOT,para, fh_info); fday.flush()
            lapw1(fday,w2k.case,w2k.WIENROOT,para,p['dftKS'],dmfe,wopt,fh_info); fday.flush()
            if wopt['so']: lapwso(fday,w2k.case,w2k.WIENROOT,para,p['dftKS'],dmfe, fh_info)
                
            fday.flush()
            
            if succi : break
            if ita>0 and ita%p['remove_broyd']==0 : 
                for f in glob.glob('*.broyd*'): 
                    os.remove(f)

        icycle += 1   # END CHARGE LOOP
        if RelaxingStructure: # When dmft/impurity loop restarts, we should relax again
            fstc.RelaxAgain(w2k.case, mix_method)
            shutil.copy2(w2k.case+'.struct',w2k.case+'.struct.'+str(icycle-1))
        
        if (p['max_lda_iterations']>=10): 
            toclean = glob.glob('*.broyd*')
            for f in toclean: os.remove(f)
        
    toclean = glob.glob('dmft1.error*')+glob.glob('dmft2.error*')+glob.glob('*.broyd*')+glob.glob('dmft')+glob.glob('dmft2')+glob.glob('_processes*')+glob.glob('*.klist_*')+glob.glob('*.output1_*')+glob.glob('*.vector_*')+glob.glob('*.scf1_*')+glob.glob('lapw1_*.def')+glob.glob('lapw1_*.error')
    for f in toclean: os.remove(f)
