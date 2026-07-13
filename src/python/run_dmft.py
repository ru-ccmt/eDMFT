#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import sys, os, re, glob, shutil, socket, time
from os.path import getsize
from copy import deepcopy
import subprocess
from functools import reduce
from numpy import *
from itertools import chain
import fileinput
###
from wstruct import Struct
from extractInfo import FindLatticeNd
from scf_merge import scf, scfm
import force_stop_conditon as fstc
import utils
import indmffile
import convert_processes
import ctqmc, inorg
try:
    import oca
except ImportError:
    oca = None
try:
    import bhqmc
except ImportError:
    bhqmc = None
try:
    import bctqmc
except ImportError:
    bctqmc = None
from imp2lattc import ImpurityLatticeConnection, SimplifySiginds, iSigindDropNegative, ConvertFromLatticeToImp, Connect_ImpurityProjector
import builtins

NMBU_DIR = os.path.dirname(os.path.abspath(__file__))

utime = '/usr/bin/time'  # This time function is compatible with w2k timing

class Params:
    """Class to store various control paramers which control the main flow of the program
       The class is constructed such that it can handle two types of parameters:
         - starting parameters which are read only at the instantiation
         - on the fly parameters which can be refreshed (reread from the file) at any time
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
    def __init__(self, f_params, params, f_sparams=None, fh_info=sys.stdout):
        self.log = fh_info
        self.p ={}
        self.p.update(params) # default values from the __main__
        for var,val in self.p.items():
            print(var,val, file=self.log)
        # remember params and its data file
        self.f_params = f_params
        # values from the 'params.dat' file
        # executes the parameters at startup
        if f_sparams is not None:
            exec(compile(open(f_sparams, "rb").read(), f_sparams, 'exec'))
            sp = locals().copy()      # takes all variables currently defined here
            del sp['self']            # deletes some locals
            del sp['params']          #
            self.p.update(sp)         # updates variables
        for var,val in self.p.items():
            print(var,val, file=self.log)

    def refresh(self):
        # executes the parameters file
        exec(compile(open(self.f_params, "rb").read(), self.f_params, 'exec'))
        sp = locals().copy()   # stores locals
        del sp['self']  # deletes automatics
        #del sp['self.log']
        # prints what was changed
        for c in sp:
            if (c not in self.p) or (sp[c] != self.p[c]):
                print('********------ parameter change ------***********', file=self.log)
                if c in self.p:
                    if isinstance(sp[c], dict):
                        for var,val in sp[c].items():
                            print('  ',c+'['+var+']->', val, file=self.log)
                    else:
                        print(c+"->", sp[c], '(',self.p[c],')', file=self.log)
                else:
                    if isinstance(sp[c], dict):
                        for var,val in sp[c].items():
                            print('  ',c+'['+var+']->', val, file=self.log)
                    else:
                        print(c+"->", sp[c], '(None)', file=self.log)
        # updates the dictonary
        self.p.update(sp)

    def __getitem__(self, x):
        return self.p[x]

    def __setitem__(self, i, value):
        self.p[i]=value

    #def keys(self):
    #    return list(self.p.keys())

def dmft1(fday, case, fh_info, extn, dmfe, m_extn=''):
    inp_name = 'indmf1'+m_extn
    name = 'dmft1'+m_extn
    if not os.path.isfile(case+'.'+inp_name):
        print('  stop error: the required input file '+case+'.'+inp_name+' for the next step could not be found!')
        print('  stop error: the required input file '+case+'.'+inp_name+' for the next step could not be found!', file=fday)
    if not os.path.isfile(name+'.def'):
        print('  stop error: the required input file '+name+'.def for the next step could not be found!')
        print('  stop error: the required input file '+name+'.def for the next step could not be found!', file=fday)

    print('>%-10s' % name,  '( '+time.strftime("%H:%M:%S")+' )', file=fday)
    fday.flush()

    with open(':log', 'a') as fl:
        #print(time.strftime("%a %b %d %H:%M:%S %Z %Y")+'>     '+name, file=fl)
        print(time.strftime("%a %b %d %I:%M:%S %p %Z %Y")+'>     '+name, file=fl)

    print('Running ---- dmft1 -----', file=fh_info)
    cmd = utime+' '+dmfe.MPI2+' '+dmfe.ROOT+'/dmft '+name+'.def >> '+name+'_info.out '
    print('#<'+name+'>: ', cmd, file=fh_info)
    fh_info.flush()
    subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)

    for fe in glob.glob('dmft1.error*'):
        if getsize(fe) !=0:
            print('ERROR in dmft1 from file:', fe, open(fe,'r').read())
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
        print('  stop error: the required input file '+case+'.'+inp_name+' for the next step could not be found!')
        print('  stop error: the required input file '+case+'.'+inp_name+' for the next step could not be found!', file=fday)

    tim = time.strftime("%H:%M:%S")
    print('>%-10s ( %s )' % (name+' '+opts, tim), file=fday)
    #fday.flush()

    cmd = WIEN+'/x'+para+' -f '+case+' '+opts+' '+name
    print(('#<'+name+'>: '), cmd, file=fh_info)
    fh_info.flush()

    info=subprocess.call(cmd,shell=True,stdout=fday)

    for fe in glob.glob(name+'.error*'):
        if getsize(fe) !=0:
            print('ERROR in '+fe+' from file:', open(fe,'r').readlines())

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
            print('>%-10s' % sname, '( '+time.strftime("%I:%M:%S %p")+' )', file=fday)
            with open(':log', 'a') as fl:
                print(time.strftime("%a %b %d %I:%M:%S %p %Z %Y")+'>     '+sname, file=fl)
            cmd = dmfe.ROOT+'/x_dmft.py'+para+' '+pcmplx+' '+name
            print(('#<'+sname+'>: '), cmd, file=fh_info)
            fh_info.flush()
            info=subprocess.call(cmd,shell=True,stdout=fh_info)
        else:                  # V_{KS} is polarized, need to lapw1 steps
            print('>%-10s' % (sname+' --up'), '( '+time.strftime("%I:%M:%S %p")+' )', file=fday)
            with open(':log', 'a') as fl:
                print(time.strftime("%a %b %d %I:%M:%S %p %Z %Y")+'>     '+sname+'up', file=fl)
            cmd = dmfe.ROOT+'/x_dmft.py'+para+' '+pcmplx+' --up '+name
            print(('#<'+sname+'>: '), cmd, file=fh_info)
            fh_info.flush()
            info=subprocess.call(cmd,shell=True,stdout=fh_info)

            print('>%-10s' % (sname+' --dn'), '( '+time.strftime("%I:%M:%S %p")+' )', file=fday)
            with open(':log', 'a') as fl:
                print(time.strftime("%a %b %d %I:%M:%S %p %Z %Y")+'>     '+sname+'dn', file=fl)
            cmd = dmfe.ROOT+'/x_dmft.py'+para+' '+pcmplx+' --dn '+name
            print(('#<'+sname+'>: '), cmd, file=fh_info)
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
        print('>%-10s' % name, '( '+time.strftime("%I:%M:%S %p")+' )', file=fday)
        with open(':log', 'a') as fl:
            print(time.strftime("%a %b %d %I:%M:%S %p %Z %Y")+'>     '+name, file=fl)
        cmd = dmfe.ROOT+'/x_dmft.py'+para+' '+name
        print(('#<'+name+'>: '), cmd, file=fh_info)
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

def dmft2(fday, case, fh_info, dmfe, m_extn=''):
    inp_name = 'indmf2'+m_extn
    name = 'dmft2'+m_extn
    if not os.path.isfile(case+'.'+inp_name):
        print('  stop error: the required input file '+case+'.'+inp_name+' for the next step could not be found!')
        print('  stop error: the required input file '+case+'.'+inp_name+' for the next step could not be found!', file=fday)
    if not os.path.isfile(name+'.def'):
        print('  stop error: the required input file '+name+'.def for the next step could not be found!')
        print('  stop error: the required input file '+name+'.def for the next step could not be found!', file=fday)

    print('>%-10s' % name, '( '+time.strftime("%I:%M:%S %p")+' )', file=fday)
    fday.flush()

    with open(':log', 'a') as fl:
        print(time.strftime("%a %b %d %I:%M:%S %p %Z %Y")+'>     '+name, file=fl)

    cmd = utime+' '+dmfe.MPI2+' '+dmfe.ROOT+'/dmft2 '+name+'.def >> dmft2_info.out '
    print(('#<'+name+'>: '), cmd, file=fh_info)
    fh_info.flush()

    info=subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)

    for fe in glob.glob(name+'.error*'):
        if getsize(fe) !=0:
            print('ERROR in '+fe+' from file:', open(fe,'r').readlines())
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


def Diff(fday, case, dEF, diff_nimp, ddiff_nimp):
    #
    with open(case+'.scf', 'r') as fs:
        dat = fs.readlines()

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
    print(':ENERGY convergence: '+str(dene), file=fday)
    print(':CHARGE convergence: '+str(drho), file=fday)
    print(':EF     convergence: '+str(dEF),  file=fday)

    print(':ENERGY convergence: '+str(dene))
    print(':CHARGE convergence: '+str(drho))
    print(':EF     convergence: '+str(dEF))
    print(':NIMP   convergence: '+str(ddiff_nimp))
    print(':NIMP   difference : '+str(diff_nimp))
    return (drho, dene)


def FileDiffMaxAbs(f1, f2):
    if (not os.path.exists(f1)) or (not os.path.exists(f2)):
        return None
    try:
        d1 = loadtxt(f1)
        d2 = loadtxt(f2)
    except Exception:
        return None
    if shape(d1) != shape(d2):
        return None
    return float(amax(abs(d1-d2)))


def LogConvergence(iteration, sigma_diff, phi_diff, fh_info):
    print(':NAMBUDIFF iter={:d} sigma_diff={} phi_diff={}'.format(iteration, sigma_diff, phi_diff), file=fh_info)


def detect_siginp_format(sig_path):
    """Detect whether sig.inp is normal (2K+1) or Nambu (8K+1) format.

    Uses both column count AND header s_oo element count to distinguish:
      - Normal K=4: 9 cols, s_oo has 4 elements → NOT Nambu
      - Nambu K=4:  33 cols, s_oo has 16 elements → IS Nambu

    Returns:
        (is_nambu, K, ncols) where K is the number of independent normal components.
    """
    import re as _re
    with open(sig_path, 'r') as f:
        line1 = f.readline()  # s_oo header
        f.readline()  # Edc header
    # Parse s_oo element count from header
    m = _re.search(r's_oo\s*=\s*(\[.*?\])', line1)
    n_soo = len(eval(m.group(1))) if m else 0

    data = loadtxt(sig_path)
    if data.ndim == 1:
        ncols = len(data)
    else:
        ncols = data.shape[1]

    # Nambu format: s_oo has 4K elements AND ncols = 8K+1
    # Normal format: s_oo has K elements AND ncols = 2K+1
    # Disambiguation: Nambu s_oo has τ₃ symmetry (second half negated).
    # Normal K=4 and Nambu K=1 both give n_soo=4, ncols=9.
    # Use: Nambu requires ncols > 2*n_soo+1 (i.e., more cols than normal format for same K).
    # Normal: ncols = 2*n_soo + 1. Nambu: ncols = 2*n_soo + 1 only if n_soo = 4K and K=n_soo/4.
    # Simplest: Nambu needs n_soo >= 8 (at least K=2 → 4K=8 → 8K+1=17 cols).
    # For K=4 plaquette: Nambu n_soo = 16, ncols = 33.
    if n_soo >= 8 and (ncols - 1) % 8 == 0:
        K_nambu = (ncols - 1) // 8
        if K_nambu > 0 and n_soo == 4 * K_nambu:
            return (True, K_nambu, ncols)
    # Normal 2K+1 format
    if (ncols - 1) % 2 == 0:
        K = (ncols - 1) // 2
        return (False, K, ncols)
    return (False, 0, ncols)


def expand_siginp_to_nambu(sig_path, phi_seed, pairing_symmetry, fh_info):
    """Expand normal 2K+1 sig.inp to Nambu 8K+1 format in-place.

    Column layout of output (K = max normal sigind):
      col 0:       omega
      cols 1..2K:  Sigma_pp (unchanged from input)
      cols 2K+1..4K: Phi (seed with pairing symmetry)
      cols 4K+1..6K: Phi_dag (conjugate of Phi seed)
      cols 6K+1..8K: Sigma_hh (-conj(Sigma_pp) from PH symmetry)

    Header arrays extended to 4K elements:
      s_oo = [s_pp_1..K, 0..0(K), 0..0(K), -s_pp_1..-s_pp_K]
      Edc  = [E_pp_1..K, 0..0(K), 0..0(K), -E_pp_1..-E_pp_K]
    """
    import re as _re

    # Read headers
    with open(sig_path, 'r') as f:
        line1 = f.readline()
        line2 = f.readline()

    m1 = _re.search(r's_oo\s*=\s*(\[.*?\])', line1)
    m2 = _re.search(r'Edc\s*=\s*(\[.*?\])', line2)
    s_oo_pp = eval(m1.group(1)) if m1 else []
    Edc_pp = eval(m2.group(1)) if m2 else []
    K = len(s_oo_pp)

    # Read data
    data = loadtxt(sig_path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    nfreq = data.shape[0]

    # Detect and strip sentinel row (last row has omega=5000.0)
    has_sentinel = (nfreq > 1 and data[-1, 0] > 4999.0)
    if has_sentinel:
        sentinel_row = data[-1:]
        data = data[:-1]
        nfreq = data.shape[0]

    omega = data[:, 0]
    sigma_pp_cols = data[:, 1:2*K+1]  # 2K columns: Re/Im pairs

    # Generate Phi seed columns (2K columns)
    # For K=4 (d-wave): loc=1, NN_x=2, NN_y=3, NNN=4
    #   d-wave (cos kx - cos ky): NN_x and NN_y have opposite signs, loc and NNN are zero
    # For K=3 (s-wave fallback): alternate signs
    phi_cols = zeros((nfreq, 2*K))
    for i in range(K):
        if pairing_symmetry == 'd' and K == 4:
            d_wave_signs = [0.0, +1.0, -1.0, 0.0]
            sign = d_wave_signs[i]
        elif pairing_symmetry == 'd':
            # Fallback for K≠4: alternate signs (old behavior)
            sign = 1.0 if i % 2 == 0 else -1.0
        else:
            sign = 1.0
        phi_cols[:, 2*i] = sign * phi_seed      # Re(Phi_i)
        phi_cols[:, 2*i+1] = 0.0                 # Im(Phi_i)

    # Phi_dag = conj(Phi): same Re, negated Im
    phi_dag_cols = phi_cols.copy()
    phi_dag_cols[:, 1::2] = -phi_cols[:, 1::2]

    # Sigma_hh = -conj(Sigma_pp): negated Re, same Im
    sigma_hh_cols = sigma_pp_cols.copy()
    sigma_hh_cols[:, 0::2] = -sigma_pp_cols[:, 0::2]  # -Re
    # Im stays the same (since -conj(Re+iIm) = -Re+iIm)

    # Assemble 8K+1 data
    expanded = column_stack([omega.reshape(-1,1), sigma_pp_cols, phi_cols, phi_dag_cols, sigma_hh_cols])

    # Build sentinel row for 8K+1 format
    if has_sentinel:
        s_oo_4k = list(s_oo_pp) + [0.0]*K + [0.0]*K + [-x for x in s_oo_pp]
        Edc_4k = list(Edc_pp) + [0.0]*K + [0.0]*K + [-x for x in Edc_pp]
        sentinel = [5000.0]
        for i in range(4*K):
            sentinel.extend([s_oo_4k[i] - Edc_4k[i], 0.0])
        expanded = vstack([expanded, array(sentinel).reshape(1,-1)])

    # Build extended headers
    s_oo_nambu = list(s_oo_pp) + [0.0]*K + [0.0]*K + [-x for x in s_oo_pp]
    Edc_nambu = list(Edc_pp) + [0.0]*K + [0.0]*K + [-x for x in Edc_pp]

    # Write back
    with open(sig_path, 'w') as f:
        f.write('# s_oo= {}\n'.format(s_oo_nambu))
        f.write('# Edc= {}\n'.format(Edc_nambu))
    with open(sig_path, 'ab') as f:
        savetxt(f, expanded)

    print('Nambu: expanded sig.inp from {}+1 to {}+1 columns (K={})'.format(
        2*K, 8*K, K), file=fh_info)


def create_basic_arrays_nambu(case, struct_nat, fh_info):
    """Generate BasicArrays.dat for Nambu mode from normal-state data.

    Reads the normal-state indmfl backup and projector binary header
    to produce a BasicArrays.dat with normal-state orbital structure,
    which is what dmftgk_nambu.py expects (Nambu doubling is internal).

    Data sources:
      - case.indmfl.normal_bak  → Sigind, cix structure, EF
      - Udmft.0 (binary header) → norbitals, nindo
      - case.klist              → nkpt
      - case.struct (via struct_nat) → nat
    """
    import struct as pystruct

    Ry2eV = 13.60569193

    # 1. Read normal-state indmfl backup
    normal_indmfl = case + '.indmfl.normal_bak'
    if not os.path.isfile(normal_indmfl):
        raise FileNotFoundError(
            "Normal-state indmfl backup not found: {}. "
            "Run a normal DMFT iteration first to generate projectors.".format(normal_indmfl))

    inl_normal = indmffile.Indmfl(case, extn='indmfl')
    inl_normal.read(normal_indmfl)

    # Extract sigind structure
    ncix = len(inl_normal.siginds)
    sorted_icps = sorted(inl_normal.siginds.keys())
    cixdim_arr = []
    csize_arr = []
    sigind_list = []
    for icp in sorted_icps:
        sig = inl_normal.siginds[icp]
        dim = len(sig)
        cixdim_arr.append(dim)
        K = int(max(abs(sig.ravel())))
        csize_arr.append(K)
        sigind_list.append(sig)

    maxdim = max(cixdim_arr)
    maxsize = max(csize_arr)

    # 2. Read projector binary header for norbitals and nindo
    udmft_file = 'Udmft.0'
    if not os.path.isfile(udmft_file):
        raise FileNotFoundError(
            "Projector file not found: {}. "
            "Run dmftu first to generate projectors.".format(udmft_file))

    with open(udmft_file, 'rb') as fh:
        # First Fortran record: numk, nsymop, norbitals
        rec_len = pystruct.unpack('i', fh.read(4))[0]
        data = fh.read(rec_len)
        pystruct.unpack('i', fh.read(4))  # trailing length
        _numk, _nsymop, norbitals = pystruct.unpack('3i', data)

        # Second record: nindo[0..norbitals-1]
        rec_len = pystruct.unpack('i', fh.read(4))[0]
        data = fh.read(rec_len)
        pystruct.unpack('i', fh.read(4))
        nindo = list(pystruct.unpack(str(norbitals) + 'i', data))

    # 3. Count nkpt from klist
    nkpt = 0
    for fname in [case + '.klist', case + '.klist_band']:
        if os.path.isfile(fname):
            with open(fname, 'r') as f:
                for line in f:
                    if line.strip().startswith('END'):
                        break
                    nkpt += 1
            break

    # 4. Read EF
    if os.path.isfile('EF.dat'):
        EF_eV = float(open('EF.dat').readline())
    elif inl_normal.EF is not None:
        EF_eV = inl_normal.EF
    else:
        EF_eV = 0.0
    EF_Ry = EF_eV / Ry2eV

    # 5. Build per-atom orbital structure from cix
    atom_orbitals = {}  # iatom → [(L, qsplit, icp)]
    for icp, entries in inl_normal.cix.items():
        for iatom, L, qsplit in entries:
            if iatom not in atom_orbitals:
                atom_orbitals[iatom] = []
            atom_orbitals[iatom].append((L, qsplit, icp))

    natom = len(atom_orbitals)
    sorted_atoms = sorted(atom_orbitals.keys())

    nl_arr = []
    ll_list = []
    iorb_list = []
    cix_list = []
    nind_list = []

    iorb_counter = 1  # 1-based
    for icase, iatom in enumerate(sorted_atoms):
        orbs = atom_orbitals[iatom]
        nl_arr.append(len(orbs))
        ll_case, iorb_case, cix_case, nind_case = [], [], [], []
        for L, qsplit, icp in orbs:
            ll_case.append(L)
            iorb_case.append(iorb_counter)
            cix_case.append(icp)
            nind_case.append(nindo[iorb_counter - 1])
            iorb_counter += 1
        ll_list.append(ll_case)
        iorb_list.append(iorb_case)
        cix_list.append(cix_case)
        nind_list.append(nind_case)

    # 6. Build iSx: sequential positions within each cix block
    cix_orb_positions = {icp: [] for icp in sorted_icps}
    for icase in range(natom):
        for lcase in range(nl_arr[icase]):
            icp = cix_list[icase][lcase]
            iorb = iorb_list[icase][lcase] - 1  # 0-based
            cix_orb_positions[icp].append(iorb)

    iSx = [[] for _ in range(norbitals)]
    for icp in sorted_icps:
        pos = 1
        for iorb in cix_orb_positions[icp]:
            iSx[iorb] = list(range(pos, pos + nindo[iorb]))
            pos += nindo[iorb]

    # 7. Write BasicArrays.dat
    nat = struct_nat
    iso = 1          # not used by dmftgk after return
    nmat = 1000      # not used
    nume = 100       # not used
    Qcomplex = True  # not used
    lmax2 = maxdim   # not used
    maxdim2 = maxdim # read from projector header but not compared
    VOL = 1.0        # not used

    lines = []
    lines.append('# BasicArrays generated by run_dmft_nambu.py for Nambu mode')
    lines.append('{} {} {} {} {}'.format(nat, iso, norbitals, ncix, natom))
    lines.append('# nkpt, nmat, nume')
    lines.append('{} {} {}'.format(nkpt, nmat, nume))
    lines.append('# Qcomplex')
    lines.append('.TRUE.' if Qcomplex else '.FALSE.')
    lines.append('# lmax2, maxdim2, maxdim, maxsize')
    lines.append('{} {} {} {}'.format(lmax2, maxdim2, maxdim, maxsize))
    lines.append('# nindo')
    lines.append(' '.join(str(x) for x in nindo))
    lines.append('# cixdim')
    lines.append(' '.join(str(x) for x in cixdim_arr))
    lines.append('# nl')
    lines.append(' '.join(str(x) for x in nl_arr))
    lines.append('# per-atom orbital data: ll, iorbital, cix, nind')
    for icase in range(natom):
        for lcase in range(nl_arr[icase]):
            lines.append(str(ll_list[icase][lcase]))
            lines.append(str(iorb_list[icase][lcase]))
            lines.append(str(cix_list[icase][lcase]))
            lines.append(str(nind_list[icase][lcase]))
    lines.append('# csize')
    lines.append(' '.join(str(x) for x in csize_arr))
    lines.append('# iSx')
    for iorb in range(norbitals):
        for ind in range(nindo[iorb]):
            lines.append(str(iSx[iorb][ind]))
    lines.append('# Sigind')
    for icix_idx, sig in enumerate(sigind_list):
        dim = cixdim_arr[icix_idx]
        for ip in range(dim):
            for iq in range(dim):
                lines.append(str(int(sig[ip, iq])))
    lines.append('# EF VOL')
    lines.append('{} {}'.format(EF_Ry, VOL))
    # posc: one line per orbital (not used but parser reads them)
    for iorb in range(norbitals):
        lines.append('0.0 0.0 0.0')

    with open('BasicArrays.dat', 'w') as fo:
        for line in lines:
            fo.write(line + '\n')

    print('Nambu: generated BasicArrays.dat (norbitals={}, ncix={}, '
          'csize={}, nkpt={}, EF={:.6f} eV)'.format(
              norbitals, ncix, csize_arr, nkpt, EF_eV), file=fh_info)


def _extract_normal_siginp(siginp_path, K):
    """Extract normal-state portion (first 2K+1 columns) from Nambu sig.inp.

    Nambu sig.inp has 8K+1 columns (omega + 4K re/im pairs).
    This extracts the first 2K+1 columns (omega + K re/im pairs for Sigma_pp)
    and truncates headers to K values, producing a valid normal-state sig.inp
    for ssplit/dmftu.

    The file is modified in-place (caller should save a backup first).
    """
    import re as _re

    with open(siginp_path, 'r') as f:
        lines = f.readlines()

    with open(siginp_path, 'w') as f:
        for line in lines:
            if line.startswith('#'):
                # Header: truncate list to first K values
                # e.g., "# s_oo= [s1, s2, s3, 0, 0, 0, 0, 0, 0, -s1, -s2, -s3]"
                #     → "# s_oo= [s1, s2, s3]"
                try:
                    eq_pos = line.index('=')
                    key = line[:eq_pos + 1]
                    val_str = line[eq_pos + 1:].strip()
                    vals = eval(val_str)
                    f.write(key + ' ' + str(vals[:K]) + '\n')
                except (ValueError, SyntaxError):
                    f.write(line)
            else:
                parts = line.split()
                if len(parts) > 2 * K:
                    # Keep omega + first K real/imag pairs = 2K+1 columns
                    f.write('  '.join(parts[:2 * K + 1]) + '\n')
                else:
                    f.write(line)


def _join_energy_files(case, so=''):
    """Merge parallel lda.energy_N files into a single lda.energy.

    After lapw1 -p, energy data is split across per-rank files.
    dmftgk_nambu.py expects a single merged file.  Header (linearization
    energies) is taken from _0; k-point blocks are concatenated in rank order.
    """
    import re as _re
    basename = case + '.energy' + so
    split_files = sorted(
        glob.glob(basename + '_[0-9]*'),
        key=lambda f: int(f.rsplit('_', 1)[-1])
    )
    if not split_files:
        return
    if os.path.isfile(basename):
        return

    header_lines = 0
    with open(split_files[0], 'r') as f:
        for line in f:
            if _re.search(r'[0-9]E[+-]', line):
                break
            header_lines += 1

    with open(basename, 'w') as fout:
        for i, sf in enumerate(split_files):
            with open(sf, 'r') as fin:
                for j, line in enumerate(fin):
                    if i == 0:
                        fout.write(line)
                    elif j >= header_lines:
                        fout.write(line)


def generate_udmft_projectors(dmfe, w2k, para, wopt, inldn, K_normal,
                               is_already_nambu, p, fh_info):
    """Generate Udmft.N projector files by running Fortran dmftu.

    In Nambu mode, dmft1 is skipped, but dmftgk reads pre-generated binary
    projector files Udmft.N (one per MPI rank). This function calls
    x_dmft.py dmftu to produce them.

    dmftu computes projectors P(k) = <chi_m|psi_{k,n}> which depend only
    on the DFT band structure (frozen in Nambu mode), not on the self-energy.
    Therefore, the normal-state indmfl gives the correct projectors.

    For restarts where indmfl is already Nambu, we temporarily restore
    the normal-state backup before running dmftu.

    On a fresh Nambu start the DFT output files (case.vsp, case.vector*)
    may be absent because lapw0/lapw1 are skipped in Nambu mode.  If so,
    they are regenerated once from the converged charge density before
    calling dmftu.
    """
    if os.path.isfile('Udmft.0'):
        print('Nambu: Udmft.0 already exists, skipping dmftu', file=fh_info)
        return

    print('Nambu: Udmft.0 not found, generating projectors via dmftu...', file=fh_info)
    fh_info.flush()

    vsp_file = w2k.case + '.vsp'
    # dmftu needs both lda.vsp (potential) and lda.vector_* (eigenvectors).
    # lapw0 generates vsp; lapw1 generates vector files.
    # Both must exist for dmftu to run successfully.
    has_vsp = os.path.isfile(vsp_file)
    has_vectors = len(glob.glob(w2k.case + '.vector_*')) > 0
    if not has_vsp or not has_vectors:
        missing = []
        if not has_vsp:
            missing.append(vsp_file)
        if not has_vectors:
            missing.append(w2k.case + '.vector_*')
        clmsum = w2k.case + '.clmsum'
        if not os.path.isfile(clmsum):
            clmsum_old = w2k.case + '.clmsum_old'
            if os.path.isfile(clmsum_old):
                shutil.copy2(clmsum_old, clmsum)
            else:
                print('ERROR: Cannot generate DFT files - '
                      'neither {} nor {} found'.format(clmsum, clmsum_old),
                      file=fh_info)
                sys.exit(1)
        print('Nambu: missing {}, running lapw0+lapw1 to regenerate '
              'DFT files from converged charge density...'.format(
                  ', '.join(missing)),
              file=fh_info)
        fh_info.flush()
        fday_proj = open(w2k.case + '.dayfile', 'a')
        lapw0(fday_proj, w2k.case, w2k.WIENROOT, para, fh_info)
        lapw1(fday_proj, w2k.case, w2k.WIENROOT, para,
              p['dftKS'], dmfe, wopt, fh_info)
        fday_proj.close()
        _join_energy_files(w2k.case, so=('.so' if wopt['so'] else ''))
        if not os.path.isfile(vsp_file):
            print('ERROR: lapw0 failed to generate {}'.format(vsp_file),
                  file=fh_info)
            sys.exit(1)
        if len(glob.glob(w2k.case + '.vector_*')) == 0:
            print('ERROR: lapw1 failed to generate vector files', file=fh_info)
            sys.exit(1)
        print('Nambu: DFT files regenerated successfully', file=fh_info)
        fh_info.flush()

    normal_indmfl = w2k.case + '.indmfl'
    saved_files = {}

    try:
        if is_already_nambu:
            # --- Restart case: temporarily restore normal-state files ---
            bak_indmfl = normal_indmfl + '.normal_bak'
            if not os.path.isfile(bak_indmfl):
                print('ERROR: Cannot generate projectors — missing ' + bak_indmfl,
                      file=fh_info)
                sys.exit(1)

            # Save current Nambu indmfl
            saved_files['indmfl'] = normal_indmfl + '.nambu_save'
            shutil.copy2(normal_indmfl, saved_files['indmfl'])
            shutil.copy2(bak_indmfl, normal_indmfl)

            # Save current Nambu sig.inp and extract normal portion
            if os.path.isfile('sig.inp'):
                saved_files['siginp'] = 'sig.inp.nambu_save'
                shutil.copy2('sig.inp', saved_files['siginp'])
                _extract_normal_siginp('sig.inp', K_normal)

            # Rebuild indmfi for normal-state sigind structure
            inl_tmp = indmffile.Indmfl(w2k.case)
            inl_tmp.read()
            if not inl_tmp.cixgrp:
                for i, icp in enumerate(sorted(inl_tmp.siginds.keys()), 1):
                    inl_tmp.cixgrp[i] = [icp]
            indmffile.Indmfi(inl_tmp, inldn).write()

        # Ensure sig.inpN files exist — Fortran dmftu opens them with
        # status='old' even though mode='u' does not read the data.
        _inl_for_ssplit = indmffile.Indmfl(w2k.case)
        _inl_for_ssplit.read()
        if not is_already_nambu:
            # First run: rebuild indmfi from current (normal) indmfl for ssplit
            if not _inl_for_ssplit.cixgrp:
                for i, icp in enumerate(sorted(_inl_for_ssplit.siginds.keys()), 1):
                    _inl_for_ssplit.cixgrp[i] = [icp]
            indmffile.Indmfi(_inl_for_ssplit, inldn).write()
        SSplit(dmfe.ROOT, _inl_for_ssplit.matsubara,
               p['iparams0'], p['ntail'], fh_info)

        # Call x_dmft.py dmftu (prepares .def file, runs Fortran dmftu)
        __cmplx = ' -c' if wopt['cmplx'] else ''
        __updn = ' --up' if wopt['updn'] else ''
        __so = ' --so' if wopt['so'] else ''
        cmd = dmfe.ROOT + '/x_dmft.py' + para + __cmplx + __updn + __so + ' dmftu'
        print('#<dmftu>: ' + cmd, file=fh_info)
        fh_info.flush()
        ret = subprocess.call(cmd, shell=True, stdout=fh_info, stderr=fh_info)
        if ret != 0:
            print('ERROR: dmftu exited with code {}'.format(ret), file=fh_info)
            sys.exit(1)

        if not os.path.isfile('Udmft.0'):
            print('ERROR: dmftu completed but Udmft.0 was not generated', file=fh_info)
            sys.exit(1)

        n_udmft = len(glob.glob('Udmft.*'))
        print('Nambu: generated {} Udmft.N projector files'.format(n_udmft), file=fh_info)

    finally:
        # Restore saved Nambu state (even if dmftu failed)
        if 'indmfl' in saved_files:
            shutil.copy2(saved_files['indmfl'], normal_indmfl)
            os.remove(saved_files['indmfl'])
        if 'siginp' in saved_files:
            shutil.copy2(saved_files['siginp'], 'sig.inp')
            os.remove(saved_files['siginp'])

SolveImpurity_num_calls=0 # number of all impurity calls
Mlist=[]  # Mlist can appear in iparams0, and have M for different MC steps
def SolveImpurity(EF, asolver, iat, extn, p, UpdateAtom, nl_imp, cix, fh_info, fh_pinfo, fday):

    #print('Running ----- impurity solver -----')
    print('UpdateAtom=', UpdateAtom, file=fh_info)
    print('Running ----- impurity solver -----', file=fh_info)
    fh_info.flush()

    tim = time.strftime("%I:%M:%S %p")
    print('>%-10s ( %s )' % ('impurity', tim), file=fday)
    fday.flush()
    with open(':log', 'a') as fl:
        print(time.strftime("%a %b %d %I:%M:%S %p %Z %Y")+'>     '+'impurity', file=fl)

    iprms = 'iparams'+str(iat)
    ipars = deepcopy(p[iprms])

    global SolveImpurity_num_calls
    global Mlist
    if 'Mlist' in ipars:
        Mlist=ipars['Mlist'][0] # setting Mlist to its value from ipars
        ipars.pop('Mlist')      # remove from parameters, since it is set as global variable
    if Mlist: # at any other time, we check if M was a list
        # Mlist would typically contain Mlist=[1e8,1e7,5e6,5e6]
        # which means, the first time ctqmc should make 1e8 steps, next 1e7 and
        # from than on only 5e6 every other time.
        if len(Mlist)>SolveImpurity_num_calls:
            Mcurrent=Mlist[SolveImpurity_num_calls]
        else:
            Mcurrent=Mlist[-1]
        ipars['M'][0]=Mcurrent
    SolveImpurity_num_calls+=1


    if p['solver'] in ('BHQMC', 'BCTQMC'):
        # BHQMC/BCTQMC keep their native 6-arg _exe for BOTH normal and Nambu; the
        # Nambu path is selected internally by .nambu (set at construction), not by
        # the call signature.  (Without grouping BCTQMC here, a Nambu BCTQMC run
        # would fall into the 4-arg iNORG branch below and crash.)
        asolver._exe(ipars, p['DCs'], extn, UpdateAtom, p['wbroad'], p['kbroad'])
    elif p['nambu'] or p['solver'] == 'iNORG':
        asolver._exe(ipars, p['DCs'], extn, cix)
    elif p['solver'] in ['OCA','NCA','CTQMC']:
        asolver._exe(ipars, p['DCs'], extn, UpdateAtom, p['wbroad'], p['kbroad'])
    else:
        print('ERROR: Impurity solver not defined!')
    fh_info.flush()

    print('Taking care of high frequency', file=fh_info)
    ntot = asolver.HighFrequency(ipars, p['DCs'], nl_imp, extn, p['wbroad'], p['kbroad'])

    with open('imp.'+str(iat)+'/sig.out','r') as fhi:
        print(p['solver'].lower()+'.HighFrequency(s_oo) :\n', next(fhi), next(fhi), file=fh_info)
    with open('imp.'+str(iat)+'/Eorb.dat', 'r') as fhi:
        print(p['solver'].lower()+'.HighFrequency(Eorb) :\n', fhi.read(), file=fh_info)

    # Reading impurity levels for printing
    (Eimp,Olap,Edc) = loadtxt('imp.'+str(iat)+'/Eimp.inp',ndmin=2)
    #fEimp = open('imp.'+str(iat)+'/Eimp.inp', 'r')
    #(Eimp,Olap,Edc) = [array([float(x) for x in line.split('#')[0].split()]) for line in fEimp.readlines()]  # throw out comments
    print('Eimp=', Eimp-Edc, file=fh_info)
    print('Edc=', Edc, file=fh_info)
    print('ncorr=', ntot, file=fh_info)
    #icycle, itt = extn.split('.')
    #print >> fh_pinfo, '%3s.%3s %12.6f %12.6f %12.6f %12.6f %12.6f' % (icycle, itt, EF, Eimp[0], Eimp[-1], Edc[0], ntot)
    #fh_pinfo.flush()
    return ntot

def combineud(case, ROOT, wopt, fh_info):
    # This averaging is necessary to not double-count the exchange splitting.
    # The exchange splitting is already taken into account in DMFT, and we do
    # not want it also here on DFT level (in the valence charge)
    cmd = ROOT+'/combineud '+case
    if wopt['updn']:
        if wopt['updn']: cmd += ' 0.5'  # We need to take only 1/2 of the average charge as up and 1/2 as down charge
        if os.path.isfile(case+'.clmval'): os.remove(case+'.clmval')
        os.symlink(case+'.clmvalup', case+'.clmval')

    print('#<combineud>: ', cmd, file=fh_info)
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
    print('#<ssplit>: ', cmd, file=fh_info)
    info=subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)

def SGather(ROOT, fh_info, p, m_extn=''):
    cmd = ROOT+'/sgather.py'
    if p['mix_sigma']<1.0 :
        cmd += ' -m '+str(p['mix_sigma'])
    if m_extn: cmd += ' -l '+m_extn
    print('#<sgather>: ', cmd, file=fh_info)
    info=subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)
    if info != 0:
        # A failed sgather means the gathered self-energy (sig.inp) was NOT rewritten, so
        # the impurity result is silently discarded and the DMFT loop would "converge" on a
        # frozen sig.inp. Fail loudly rather than produce a meaningless self-consistent cycle.
        print('ERROR: sgather.py exited with code {} -- sig.inp was NOT updated; aborting.'.format(info), file=fh_info)
        fh_info.flush()
        sys.exit(1)

def SJoin(ROOT, fh_info, p, m_extn=''):
    cmd = ROOT+'/sjoin.py -m '+str(p['mix_delta'])
    if m_extn: cmd += ' -l '+m_extn + ' -a '+str(p['EimpAverage'])
    if p['rCF'] is not None: cmd += ' -c "'+ p['rCF']+'"'
    print('#<sjoin>: ', cmd, file=fh_info)
    info=subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)

def Prepare_dmft1(dmfe, p, para, wopt, fh_info):
    __updn = '--up'  if wopt['updn'] else ''
    print('--------- Preparing dmft1 calculation ---------', file=fh_info)
    cmd = dmfe.ROOT+'/x_dmft.py -d' + para + ' dmft1 '+__updn+' >> dmft1_info.out 2>&1'
    print('#<prep-dmft1>: ', cmd, file=fh_info)
    info=subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)
    if wopt['m_extn']:
        __updn = '--dn'  if wopt['updn'] else ''
        print('--------- Preparing dmft1 calculation '+wopt['m_extn']+' ---------', file=fh_info)
        cmd = dmfe.ROOT+'/x_dmft.py -d' + para + ' dmft1 -l dn '+__updn+' >> dmft1_info.out 2>&1'
        print('#<prep-dmft1dn>: ', cmd, file=fh_info)
        info=subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)

def Prepare_dmft2(dmfe, p, para, wopt, fh_info, recomputeEF, mode='c'):
    print('--------- Preparing dmft2 calculation ---------', file=fh_info)
    __updn = ' --up'  if wopt['updn'] else ''
    cmd = dmfe.ROOT+'/x_dmft.py -d '+'--mode '+mode+' -m '+str(recomputeEF)+' -x '+str(p['mixEF'])+' -w '+str(p['WL'])+__updn+para+' dmft2 >> dmft2_info.out 2>&1'
    print('#<prep-dmft2>: ', cmd, file=fh_info)
    info=subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)
    if wopt['m_extn']:
        print('--------- Preparing dmft2 calculation '+wopt['m_extn']+'---------', file=fh_info)
        __updn = ' --dn'  if wopt['updn'] else ''
        cmd = dmfe.ROOT+'/x_dmft.py -d '+'--mode '+mode+' -m 0'+__updn+para+' dmft2 -l dn >> dmft2_info.out 2>&1'
        print('#<prep-dmft2>: ', cmd, file=fh_info)
        info=subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)

def CreateEorb(case, wopt, siginds, imp2latt):
    """ Reads case.scf2 file (prepared by dmft2) and extracts
    density matrix, (density on each atom) which starts with :NCOR.
    Also gives the first component of double-counting
    """
    with open('sig.inp','r') as fsg:
        exec(next(fsg)[1:].lstrip(), globals())  # s_oo from sig.inp                                                                                                                   $
        exec(next(fsg)[1:].lstrip(), globals())  # Edc  from sig.inp                                                                                                                   $
    up__ = 'up' if wopt['updn'] else ''
    with open(case+'.scf2'+up__,'r') as fsc:
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

def FindEtot(case,wopt):
    def FindLine(dat, strng, item):
        Val = 0.0
        for line in dat:
            if re.match(strng, line) is not None:
                Val = float(line.split()[item])
                break
        return Val

    with open(case+'.scfm','r') as fs:
        dat = fs.readlines()[::-1]  # it sorts from the end back
    ETOT = FindLine(dat, ':ENE', 8) # last :ENE

    if os.path.isfile('Eorb_imp.dat'):
        with open('Eorb_imp.dat','r') as fs:
            dat = fs.readlines()[::-1]
        XEORB = FindLine(dat, ':XEORB',1)
        EORB  = FindLine(dat, ':EORB',1)
        ZEORB = FindLine(dat, ':ZEORB',1)
    else:
        XEORB = EORB = ZEORB = 0.

    up__ = 'up' if wopt['updn'] else ''
    with open(case+'.scf2'+up__,'r') as fs:
        dat = fs.readlines()[::-1]
    SUM  = FindLine(dat, ':SUM ', 6)
    XSUM = FindLine(dat, ':XSUM', 6)
    YSUM = FindLine(dat, ':YSUM', 6)
    ZSUM = FindLine(dat, ':ZSUM', 6)
    if wopt['m_extn']:
        with open(case+'.scf2'+wopt['m_extn'],'r') as fsn:
            datn = fsn.readlines()[::-1]
        SUM_ = FindLine(datn, ':SUM ', 6)
        XSUM_ = FindLine(datn, ':XSUM', 6)
        YSUM_ = FindLine(datn, ':YSUM', 6)
        ZSUM_ = FindLine(datn, ':ZSUM', 6)
        SUM = 0.5*(SUM+SUM_)
        XSUM = 0.5*(XSUM+XSUM_)
        YSUM = 0.5*(YSUM+YSUM_)
        ZSUM = 0.5*(ZSUM+ZSUM_)
    return (ETOT,SUM,XSUM,YSUM,ZSUM,EORB,XEORB,ZEORB)

def FermiF(dmfe, w2k, p, para, wopt, fh_info, fday, matsubara):
    # In case of ferromagnetic or ferrimagnetic metal, we need to compute common EF for up and dn spins
    case = w2k.case
    scratch = w2k.SCRATCH
    Prepare_dmft2(dmfe, p, para, wopt, fh_info, 1, 'e')
    SSplit(dmfe.ROOT, matsubara, p['iparams0'], p['ntail'], fh_info)
    dmft2(fday,case,fh_info,dmfe) #.MPI2)
    SSplit(dmfe.ROOT, matsubara, p['iparams0'], p['ntail'], fh_info, wopt['m_extn'])
    dmft2(fday,case,fh_info,dmfe,wopt['m_extn'])
    cmd = dmfe.ROOT+'/fermif '+scratch+'/'+w2k.case+'.bnds '+scratch+'/'+w2k.case+'.bnds'+wopt['m_extn']
    print(('#<fermif>: '), cmd, file=fh_info)
    fh_info.flush()
    tim = time.strftime("%I:%M:%S %p")
    print('>%-10s ( %s )' % ('fermif', tim), file=fday)
    fday.flush()
    with open(':log', 'a') as fl:
        print(time.strftime("%a %b %d %I:%M:%S %p %Z %Y")+'>     '+'fermif', file=fl)
    info=subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)
    if info!=0:
        print('Problems evaluating fermif!', file=fh_info)

if __name__ == '__main__':
    # -------------- Parameters which are set through input files sparams.dat and params.dat -----------
    # Default values of some parameters
    params = {'runIMP'        : True,      # Run impurity solver
              'runGF'         : True,      # Run dmft1 for G and Delta
              'UpdateAtom'    : 1,         # recompute cix-file for the impurity solver
              'DCs'           : 'nominal',    # the double counting scheme, which fixes Edc through n0
              'max_dmft_iterations': 1,    # number of iteration of the dmft-loop only
              'max_lda_iterations' : 1,    # number of iteration of the LDA-loop only
              'max_lda_iterations_optimize' : 1000,    # number of iteration of the LDA-loop when structure should be optimized
              'finish'        : 50,        # number of iterations of the charge loop
              'da'            : 0.0,       # experimental
              'wbroad'        : 0.0,       # broadening of sigma on the imaginary axis
              'kbroad'        : 0.0,       # broadening of sigma on the imaginary axis
              'solver'        : 'CTQMC',   # CTQMC solver
              'ntail'         : 200,       # on imaginary axis, number of points in the tail of the logarithmic mesh
              'mix_delta'     : 1.0,       # whether to mix delta, or not
              'mix_sigma'     : 1.0,       # whether to mix sigma, or not
              'riter'         : 100,        # How often to restart broyden for charge mixing
              'sleeptime'     : 2,         # If broyden file are present at the submit time, user has some time to kill the job
              'cc'            : 1e-5,      # the charge density precision to stop the LDA+DMFT run
              'ec'            : 1e-5,      # the energy precision to stop the LDA+DMFT run
              'nc'            : 1e-10,     # the impurity difference to stop LDA+DMFT run
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
              'if_details'    : False,     # If True, show detailed NORG output during particle search in iNORG solver
              'nambu'           : False,    # True: enable Nambu d-wave SC pre-flight + loop
              'phi_seed'        : 0.01,     # initial pairing seed amplitude (eV)
              'mix_phi'         : 1.0,      # mixing for anomalous self-energy
              'pairing_symmetry': 'd',      # seed symmetry: 's' or 'd'
              }

    Ry2eV = utils.Ry_in_eV

    if len(sys.argv)>1 and sys.argv[1] in ['-h', '--help']:
        help="""
        The script runs LDA+DMFT. It is a wrapper, which calls other scripts and executables.
        It executes Wien2K in parallel using dftKS program (mpi parallelization over k-points)
        or it can use Wien2k parallelization.
        The DMFT part will run in parallel even if DFT part is sequential.
        Usuall execution goes through the following steps in a loop:

          x_lapw    lapw0        -- computes LDA potential on current LDA+DMFT charge
          x_dmft.py lapw1        -- solves LDA eigenvalue problem
          [x_dmft.py lapwso]     -- adds spin-orbit
          x_dmft.py  dmft1       -- computes local green's function and hybridization
          run impurity   -- runs impurity problem
          x_dmft.py dmft2        -- computes LDA+DMFT valence charge, and the chemical potential
          x_lapw lcore        -- computes LDA core charge
          x_lape mixer        -- mixes total charge with the previous charge

        The most common parameters, which should be given through 'params.dat' file, include:

        name       possible values   default     help
        --------------------------------------------------------------------------------------
        dftKS          [True|False]  True          # uses internal lapw1 & lapwso rather than w2k equivalent
        solver         [CTQMC|OCA|INORG ] CTQMC    # impurity solver
        max_iterations [int]         1             # number of iteration of the dmft-loop only
        finish         [int]         100           # number of iterations of the charge+dmft loop
        cc             [float]       1e-5,         # the charge density precision to stop the LDA+DMFT run
        ec             [float]       1e-5,         # the energy precision to stop the LDA+DMFT run
        broyd          [True| False] True          # Are we using broyden for charge mixing
        riter          [int]         99            # How often to restart broyden for charge mixing
        sleeptime      [int]         2             # If broyden file are present at the submit time, user
                                                    # has some time to kill the job
        DCs            [exacty | exactd | exact | nominal | fixn | FLL | AMF | exact_l ]
                       exacty : dielectric+yukawa screening, such that we get predescribed U & J. Real space form is e^{-lambda*r}/(eps*r)
                       exactd : dielectric screening only, i.e., U_{1234}=<psi_1|psi_2|U|psi_3|psi_4>/eps. U is always make equal to input U, but J is not, and is typically smaller.
                       exact  : dielectric+yukawa screening, first computing lambda,epsilon from U&J, and then renormalizing
                                lambda by fractional occupancy of the shell (very little yukawa screening for late TMO and mostly yukawa for early TMOs).
                                dielectric screening is later increased to compensate for descreased lambda, to give input U. But J is probably smalller
                                than required by input.
        runIMP         [True| False] True          # Run impurity solver or skip this step
        runGF          [True| False] True          # Run dmft1 for G and Delta, or skip this step
        ntail          [int]         30            # sigma.inp on imaginary axis, number of points in the
                                                     # tail of the logarithmic mesh
        mix_delta      [float]       1.0           # linear mixing parameter for Delta.
        UpdateAtom     [int]         0             # How often to recompute cix-file inside the
                                                   # impurity solver
        wbroad         [float]       0.0           # constant broadening of sigma on the imaginary axis
        kbroad         [float]       0.0,          # linear broadening of sigma on the imaginary axis
        recomputeEF    [0|1|2]       1             # If we want to recompute the chemical potential in dmft2 step. 2->for Mott insulator
        com            [0| 1]        0             # when computing chemical potential in [dmft0|mu] step
                                                     # on imaginary axis, we can subtract Sigma_oo (com=0)
                                                     # or Sigma(0) (com=1)
        mixEF          0.0-1.0                     # mixing of the chemical potential
        da             [float]       0.0           # experimental - connected with DC scheme

        --- Nambu d-wave superconductivity (set nambu=True in params.dat) ---
        nambu            [True|False] False         If True, run d-wave Nambu CDMFT (single-shot, frozen mu)
        phi_seed         [float] 0.01               Initial pairing seed amplitude (eV)
        mix_phi          [float] 1.0                Anomalous self-energy mixing factor
        pairing_symmetry ['s'|'d'] 'd'              Pairing seed symmetry

        """
        print(help)
        sys.exit(0)

    fh_info = open('dmft_info.out','w')
    fh_pinfo = open('info.iterate', 'a')

    p = Params('params.dat', params, f_sparams=None, fh_info=fh_info)

    # list of options for wien2k
    wopt={'cmplx':'', 'so':'', 'm_extn':'', 'updn':''}
    # cmplx : no-inversion symmetry present
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

    # we want to clean certain files that might interfere with this run
    to_clean = ['*.scf','Edc.dat','dmft[1|2].error*','*.[dlt|gc1|cdos].*','*.outputdmf[1|2].*','dmft[1|2]_info.out', '*.broyd*']
    toclean = []
    for pat in to_clean:
        toclean.extend( glob.glob(pat) )
    for f in toclean: os.remove(f)

    # Info files
    fday = open(w2k.case+'.dayfile', 'w')
    print('Calculating '+w2k.case+' in '+os.getcwd()+'\n'+'on '+socket.gethostname()+' with PID '+str(os.getpid()), file=fday)
    print('\n\n', file=fday)


    # Reading 'case.indmfl' file
    inl = indmffile.Indmfl(w2k.case) # case.indmfl file
    inl.read()                       # case.indmfl read

    print('inl.emin={:s}, inl.emax={:s} inl.hybr_emin={:f} inl.hybr_emax={:f}'.format(str(inl.emin),str(inl.emax),inl.hybr_emin,inl.hybr_emax), file=fh_info)
    # If FM or AFM withouth SO, we read also case.indmfldn file
    inldn = None
    if os.path.exists(w2k.case+'.indmfl'+'dn'): wopt['m_extn'] = 'dn'
    if wopt['m_extn']:
        print('INFO: case.indmfldn present => magnetic calculation with two dmft2 steps', file=fh_info)
        inldn = indmffile.Indmfl(w2k.case, 'indmfl'+wopt['m_extn'])
        inldn.read()

    ## produce projectorw.dat if it does not yet exist.
    if inl.projector>=4 and not os.path.exists('projectorw.dat'):
        print('WARNING: projectorw.dat does not exist, hence trying to computed it. Make sure case.vsp is available', file=fh_info)
        fh_info.flush()
        if not os.path.isfile(w2k.case+'.vsp'):
            print('ERROR: projectorw.dat does not exist nor case.vsp, hence I can not compute projector. Run "init_proj.py" first.', file=fh_info)
            print('ERROR: projectorw.dat does not exist nor case.vsp, hence I can not compute projector. Run "init_proj.py" first.')
            sys.exit(1)
        inl.CmpProjector(log=fh_info)
    ## if the band range is not yet specified, do it now
    if inl.emin==inl.emax==0: # we did not yet determine band ranges
        inl.CmpBandRange(log=fh_info)
        inl.write(w2k.case+'.indmfl')
        if wopt['m_extn']:
            inldn.emin,inldn.emax = inl.emin, inl.emax
            inldn.write(w2k.case+'.indmfl'+wopt['m_extn'])

    # Reading parameters from params.dat
    p.refresh()

    # Reads structure file
    struct = Struct(inAngs=False)
    struct.ReadStruct(w2k.case+'.struct', fh_info)

    # Check spin-orbit coupling
    if os.path.isfile(w2k.case+".inso") and os.path.getsize(w2k.case+".inso")>0 :
        print('Found '+w2k.case+'.inso file, hence assuming so-coupling exists. Switching -so switch!', file=fh_info)
        wopt['so']='so'
    if os.path.isfile(w2k.case+".in1c") and os.path.getsize(w2k.case+".in1c")>0 :
        print('Found '+w2k.case+'.in1c file, hence assuming non-centrosymmetric structure. Switching -c switch!', file=fh_info)
        wopt['cmplx']='c'

    # corelated indexes
    cixs = list(inl.siginds.keys())        # all columns for cix
    ncix = len(cixs)

    nambu_cc = None
    if p['nambu']:
        p['recomputeEF'] = 0
        p['max_lda_iterations'] = 0
        p['finish'] = 1
        nambu_cc = p['cc']
        if ('m_extn' in p.p) and p['m_extn']:
            raise ValueError('Nambu mode incompatible with magnetic extension (m_extn)')
        if wopt['m_extn']:
            raise ValueError('Nambu mode incompatible with magnetic extension (case.indmfldn)')

        # --- Unified Nambu sig.inp architecture ---
        # 1. Generate Nambu indmfl (8x8 sigind) from normal-state indmfl (4x4)
        from nambu_indmfl_gen import generate_nambu_indmfl, detect_nambu_indmfl
        normal_indmfl = w2k.case + '.indmfl'
        nambu_indmfl = w2k.case + '.indmfl_nambu'

        is_already_nambu, K_normal = detect_nambu_indmfl(normal_indmfl, w2k.case)

        # 0a. Generate Udmft.N projector files if missing (BUG #5 fix)
        # MUST run with K=3 indmfl — Fortran dmftu cannot handle K=4 sigind.
        # Projectors depend only on DFT bands (frozen), not on sigind structure.
        generate_udmft_projectors(dmfe, w2k, para, wopt, inldn, K_normal,
                                   is_already_nambu, p, fh_info)

        # 0b. Auto-convert normal-state K=3 → K=4 (split NN→NN_x/NN_y for d-wave)
        # AFTER projectors generated (dmftu needs K=3), BEFORE Nambu expansion.
        # This ensures BasicArrays.dat (regenerated at line ~1524) uses K=4.
        if not is_already_nambu and K_normal == 3:
            from convert_K3_to_K4 import convert_siginp, convert_indmfl, convert_indmfi
            print('Nambu: detected K=3 normal-state, converting to K=4 (NN→NN_x/NN_y)', file=fh_info)
            if os.path.isfile('sig.inp'):
                convert_siginp('sig.inp')
            convert_indmfl(normal_indmfl)
            if os.path.isfile(w2k.case + '.indmfi'):
                convert_indmfi(w2k.case + '.indmfi')
            # Re-read converted indmfl
            inl = indmffile.Indmfl(w2k.case)
            inl.read()
            K_normal = 4
            # Force BasicArrays.dat regeneration with K=4 csize
            if os.path.isfile('BasicArrays.dat'):
                os.remove('BasicArrays.dat')
                print('Nambu: removed stale BasicArrays.dat (was csize=3, need csize=4)', file=fh_info)
            print('Nambu: K=3→K=4 conversion done', file=fh_info)

        if not is_already_nambu:
            # Save original and generate Nambu version
            shutil.copy2(normal_indmfl, normal_indmfl + '.normal_bak')
            generate_nambu_indmfl(normal_indmfl, nambu_indmfl, w2k.case)
            shutil.copy2(nambu_indmfl, normal_indmfl)
            print('Nambu: generated 8x8 indmfl from normal-state (K={})'.format(K_normal), file=fh_info)
            # Re-read the updated indmfl
            inl = indmffile.Indmfl(w2k.case)
            inl.read()
            cixs = list(inl.siginds.keys())
            ncix = len(cixs)

        # 2. Regenerate case.indmfi to match Nambu indmfl
        # Indmfl.read() does not populate cixgrp (only Indmf does),
        # so synthesize it from siginds before calling Indmfi()
        if not inl.cixgrp:
            for i, icp in enumerate(sorted(inl.siginds.keys()), 1):
                inl.cixgrp[i] = [icp]
        indmfi = indmffile.Indmfi(inl, inldn)
        indmfi.write()
        print('Nambu: regenerated case.indmfi for Nambu sigind', file=fh_info)

        # 3. Auto-convert sig.inp from 2K+1 to 8K+1 if needed
        is_nambu_sig, K_sig, ncols_sig = detect_siginp_format('sig.inp')
        if not is_nambu_sig:
            pairing_sym = p['pairing_symmetry'] if 'pairing_symmetry' in p.p else 's'
            phi_seed_val = p['phi_seed'] if 'phi_seed' in p.p else 0.01
            expand_siginp_to_nambu('sig.inp', phi_seed_val, pairing_sym, fh_info)
        else:
            print('Nambu: sig.inp already in 8K+1 format ({} cols)'.format(ncols_sig), file=fh_info)

    # Nuclear charges of all correlated blocks
    Znuc={}
    for icix in cixs:
        atm = inl.cix[icix][0][0]                # correlated atoms
        Znuc[icix] = int(round(struct.aZnuc[atm-1]))  # nuclear charge of each correlated block

    print('Znucs=', Znuc, file=fh_info)

    # Independent impurity problems
    iSigind = indmffile.ParsIndmfi(w2k.case)
    print('iSiginds=', iSigind, file=fh_info)
    cols = SimplifySiginds(inl.siginds)
    icols = SimplifySiginds(iSigind)
    imp2latt = ImpurityLatticeConnection( cols, icols, fh_info)
    print('Impurity-lattice connection: imp2latt=', imp2latt, file=fh_info)

    fh_pinfo.write( '%3s %3s.%4s %12s %12s ' % ('#', '#','#','mu','Vdc') )
    fh_pinfo.write( '%15s %15s %15s ' % ('Etot', 'Ftot+T*Simp', 'Ftot+T*Simp') )
    for iat in list(iSigind.keys()):   # Over all inequivalent impurity problems
        fh_pinfo.write( '%12s %12s %12s %12s ' % ('n_latt', 'n_imp','Eimp[0]','Eimp[-1]') )
    fh_pinfo.write('\n')
    fh_pinfo.flush()


    if p['DCs'][:5]=='exact':
        Nmax=[]
        for imp in list(iSigind.keys()):
            icx = imp2latt[imp][0]
            iatom, l, qsplit = inl.cix[icx][0]  # Doesn't take into account cluster problems correctly
            Nmax.append(2*(2*l+1))

        # Finds connection between impurity and projector
        impurity_projector = Connect_ImpurityProjector(imp2latt,inl,struct,fh_info)
        print('connection -- impurity_projector=', impurity_projector, file=fh_info)
        # Finds the screening length and corresponding J's given imput value of Coulomb U for al impurity problems
        UJ_icase={}
        for imp,icase in enumerate(impurity_projector):
            iprms = 'iparams'+str(imp)
            U=p[iprms]['U'][0]            # we always need U to determine screening length
            JH=0
            if 'J' in p[iprms]:    # if Hunds J is awailable, we take it into account
                JH = p[iprms]['J'][0]
            lmbda=1e-6
            if 'DC_lmbda' in p[iprms]: # if user wants to use finite lambda with exactd DC, it should set iparamsx['DC_lmbda']
                lmbda = p[iprms]['DC_lmbda'][0]
            nfraction=1.0
            if 'nf0' in p[iprms]:
                n0 = p[iprms]['nf0'][0]
                nfraction = (Nmax[imp]-n0)/(Nmax[imp]+0.0) # fractional occupancy
            UJ_icase[icase]=[U,JH,lmbda,nfraction]  # Saves U into U_icase, which is needed below to compute lambda

        print('UJ_icase=', UJ_icase, file=fh_info)
        import RCoulombU
        if len(p['DCs'])>5 and p['DCs'][:6]=='exactd':      # dielectric screening only
            print('Have exactd!', file=fh_info)
            UlamJ = RCoulombU.GetDielectricFunctions(UJ_icase, log=fh_info)  # Computes dielectric function epsilon and sets yukawa lambda=1e-6
        elif len(p['DCs'])>5 and p['DCs'][:6]=='exacty':    # yukawa + dielectric screening, should be best
            print('Have exacty!', file=fh_info)
            for icase in list(UJ_icase.keys()):
                UJ_icase[icase][3]=1.0 # all nfractions equal to 1.0, hence yukawa
            #print('UJ_icase=', UJ_icase)
            UlamJ = RCoulombU.GetLambdas(UJ_icase, log=fh_info) # Computes yukawa lambda and dielectric function so that they match user defined U & J.
        else :  # mixture between yukawa and dielectric. First we compute lambda&epsilon from input U&J. In the second step we reduce
                # lambda (Yukawa screening) by partial occupancy so that d^9 or d^10 system has very little or no yukawa screening, while d^0 d^1 has
                # yukawa & epsilon corresponding to input U&J. Then we increase dielectric screening to compensate for smaller Yukawa lambda.
                # This is done because in cuprates, for example, exactd seems to be the best, while in early TMO's yukawa is fine.
            print('Have exact!', file=fh_info)
            #print('UJ_icase=', UJ_icase)
            UlamJ = RCoulombU.GetLambdas(UJ_icase, log=fh_info)

        for imp,icase in enumerate(impurity_projector): # Saves into impurity
            iprms = 'iparams'+str(imp)
            p[iprms]['icase'] = [icase, '#'] # Connection between impurity and projector in projector.dat'
            p[iprms]['UlamJ'] = [UlamJ[icase], '#']

        print('UlamJ = {iatom : (U,lambda,epsilon,[J2,J4,...]), ....}', file=fh_info)
        print('UlamJ=', UlamJ, file=fh_info)

    # Solvers initialized
    if p['nambu']:
        # Unified inorg.IMP_INORG handles both normal and Nambu paths via nambu kwarg
        solver=[]
        for ic,isigind in iSigind.items():
            if len(imp2latt[ic])==0: continue
            icx = imp2latt[ic][0]
            iatom, l, qsplit = inl.cix[icx][0]
            cftrans = inl.cftrans[icx]
            if shape(isigind)[0]==2*shape(cftrans)[0] and shape(isigind)[1]==2*shape(cftrans)[1]:
                zr = zeros(shape(cftrans))
                cftrans = vstack((hstack((cftrans,zr)),hstack((zr,cftrans))))
            isigind_ = iSigindDropNegative(isigind,icols[ic],fh_info)
            Zn, imp_dir, imp_params = Znuc[icx], 'imp.'+str(ic), p['iparams'+str(ic)]
            if p['solver'] == 'BHQMC':
                # Nambu BHQMC (d-wave SC) via bold-CTQMC. isigind_ is the doubled
                # Nambu sigind (Spp|Phi|PhiD|Shh); IMP_BHQMC(nambu=True) extracts
                # the pp-block and adds d-wave anomalous baths via gen_cix(superc=1).
                if bhqmc is None:
                    raise ImportError('BHQMC (Nambu) requested but bhqmc module is not available')
                solver.append(bhqmc.IMP_BHQMC(Zn, imp_dir, imp_params, isigind_, cftrans, nambu=True))
            elif p['solver'] == 'BCTQMC':
                # K-native Nambu (d-wave SC) via Bold-CTQMC. isigind_ is the doubled
                # Nambu sigind (Spp|Phi|PhiD|Shh), already K-diagonal; IMP_BCTQMC(nambu=True)
                # reduces it inline (no site<->K Fourier) and builds the d-wave anomalous
                # baths via gen_cix(superc=1).
                if bctqmc is None:
                    raise ImportError('BCTQMC (Nambu) requested but bctqmc module is not available')
                solver.append(bctqmc.IMP_BCTQMC(Zn, imp_dir, imp_params, isigind_, cftrans, nambu=True))
            else:
                solver.append(inorg.IMP_INORG(Zn, imp_dir, imp_params, inl.matsubara, isigind_, cftrans, l, nambu=True))
    elif p['solver'] in ['OCA', 'NCA','CTQMC','iNORG','BHQMC','BCTQMC']:
        ###first information message about which frequency axis is used
        if inl.matsubara:
            if p['solver'] in ['OCA', 'NCA']:
                raise ValueError('ERROR: Can not run OCA/NCA on imaginary axis! Change case.indmfl:matsubara or params.dat:solver')
            elif p['solver'] in ['CTQMC','BHQMC','BCTQMC']:
                print('Running '+p['solver']+' solver on the imaginary axis!')
            elif p['solver'] in ['iNORG']:
                print('Running '+p['solver']+' solver on the imaginary axis!')
        if not inl.matsubara:
            if p['solver'] in ['CTQMC','BHQMC','BCTQMC']:
                raise ValueError('ERROR: Can not run CTQMC/BHQMC/BCTQMC on real axis! Change case.indmfl:matsubara or params.dat:solver')
            elif p['solver'] in ['OCA', 'NCA']:
                print('Running '+p['solver']+' solver on the imaginary axis!')
            elif p['solver'] in ['iNORG']:
                raise ValueError('Running iNORG solver on real axis is not recommended for SCF calculations, boiling out..')

        ### Then setup all unequivalents impurity problems
        solver=[]
        for ic,isigind in iSigind.items():
            if len(imp2latt[ic])==0: continue  # This atom is treated as open-core
            icx = imp2latt[ic][0]
            iatom, l, qsplit = inl.cix[icx][0]  # Doesn't take into account cluster problems correctly
            cftrans = inl.cftrans[icx]
            if shape(isigind)[0]==2*shape(cftrans)[0] and shape(isigind)[1]==2*shape(cftrans)[1]:
                zr = zeros(shape(cftrans))
                cftrans = vstack((hstack((cftrans,zr)),hstack((zr,cftrans))))
            isigind_ = iSigindDropNegative(isigind,icols[ic],fh_info)
            Zn, imp_dir, imp_params = Znuc[icx], 'imp.'+str(ic), p['iparams'+str(ic)]
            if p['solver'] in ['OCA', 'NCA']:
                if oca is None:
                    raise ImportError('OCA/NCA solver requested but oca module is not available')
                solver.append( oca.IMP_OCA(Zn, imp_dir, imp_params, isigind_, cftrans) )
            if p['solver'] in ['CTQMC']:
                solver.append( ctqmc.IMP_CTQMC(Zn, imp_dir, imp_params, isigind_, cftrans) )
            if p['solver'] in ['iNORG']:
                # add explicitely the number of max dmft iterations, useful for sector update
                imp_params['max_dmft_iterations'] = p['max_dmft_iterations']
                imp_params['if_details'] = p['if_details']
                solver.append( inorg.IMP_INORG(Zn, imp_dir, imp_params, inl.matsubara, isigind_, cftrans, l) )
            if p['solver'] in ['BHQMC']:
                if bhqmc is None:
                    raise ImportError('BHQMC solver requested but bhqmc module is not available')
                solver.append( bhqmc.IMP_BHQMC(Zn, imp_dir, imp_params, isigind_, cftrans) )
            if p['solver'] in ['BCTQMC']:
                if bctqmc is None:
                    raise ImportError('BCTQMC solver requested but bctqmc module is not available')
                solver.append( bctqmc.IMP_BCTQMC(Zn, imp_dir, imp_params, isigind_, cftrans) )

    else:
        print('ERROR: No correct solver specified! Bolling out.')
        sys.exit(1)

    # We decided not to copy executables anymore. Too much temp files
    #shutil.copy2(dmfe.ROOT+'/dmft', '.')  # For parallel execution, executable has to be on the current directory
    #shutil.copy2(dmfe.ROOT+'/dmft2', '.') # For parallel execution, executable has to be on the current directory

    if os.path.isfile(w2k.case+".incup") and os.path.isfile(w2k.case+".incdn") and os.path.getsize(w2k.case+".incup")>0 and os.path.getsize(w2k.case+".incdn")>0:
        wopt['updn']=True # Computing J_Hunds by constrained-DMFT, hence need different core density for up and dn spin.
        print('Found '+w2k.case+'.incup file, hence assuming core is spin polarized, and switching updn=up/dn on', file=fh_info)
        if not os.path.isfile(w2k.case+'.clmup') or os.path.getsize(w2k.case+'.clmup')==0:
            print('ERROR: Missing '+w2k.case+'.clmup /clmdn files, which are needed when '+w2k.case+'.incup /incdn files are present', file=fh_info)
            sys.exit(1)

    if p['recomputeEF']==0:
        # W2k now changed mixer, so that it does not normalize charge. If we have recomputeEF=0, it will result in error. Hence we need to change case.inm
        fname = w2k.case+'.inm'
        if os.path.exists(fname):
            with open(fname, "r") as fi:
                content = fi.read()
            m = re.match(r"^(\w+\s+\d+(\.\d+)?\s+)YES", content, flags=re.MULTILINE)
            if m is not None:
                content = re.sub(r"^(\w+\s+\d+(\.\d+)?\s+)YES", r"\1NO", content, flags=re.MULTILINE)
                print('since recomputeEF=0 we need to allow unrenormalized charge in mixer. Changing case.inm to', file=fh_info)
                print(content, file=fh_info)
                with open(fname, "w") as fi:
                    fi.write(content)
        if os.path.exists(w2k.case+'.scf2') and not os.path.exists('EF.dat'):
            ef = None
            with open(w2k.case+'.scf2', 'r') as fi:
                lines = fi.readlines()
                for line in lines[::-1]:
                    if line[:4]==':FER':
                        ef = float(line.split()[-1])
                        break
            if ef is not None:
                with open('EF.dat', 'w') as fo:
                    print(ef*Ry2eV, file=fo)

    if not p['nambu']:
        Prepare_dmft1(dmfe, p, para, wopt, fh_info)
        Prepare_dmft2(dmfe, p, para, wopt, fh_info, p['recomputeEF'])
    else:
        # Nambu mode: dmft1 is skipped; synthesize BasicArrays.dat for dmftgk_nambu
        regenerate = not os.path.isfile('BasicArrays.dat')
        if not regenerate:
            try:
                with open('BasicArrays.dat', 'r') as fba:
                    ba_lines = [l.strip() for l in fba]
                ba_nkpt = int(ba_lines[3].split()[0])
                klist_nkpt = 0
                for kf in [w2k.case + '.klist', w2k.case + '.klist_band']:
                    if os.path.isfile(kf):
                        with open(kf, 'r') as fk:
                            for line in fk:
                                if line.strip().startswith('END'):
                                    break
                                klist_nkpt += 1
                        break
                if klist_nkpt > 0 and ba_nkpt != klist_nkpt:
                    print('Nambu: BasicArrays.dat nkpt={} != klist nkpt={}, regenerating'.format(
                        ba_nkpt, klist_nkpt), file=fh_info)
                    regenerate = True
                else:
                    print('Nambu: BasicArrays.dat already exists (nkpt={}), reusing'.format(ba_nkpt), file=fh_info)
            except (IndexError, ValueError):
                print('Nambu: BasicArrays.dat unreadable, regenerating', file=fh_info)
                regenerate = True
        if regenerate:
            create_basic_arrays_nambu(w2k.case, struct.nat, fh_info)

    shutil.copy2('sig.inp', 'sig.inp.0')

    if (not p['nambu']) and p['finish']>1: # Charge self-consistent

        if os.path.isfile(w2k.case+'.clmsum_old'): shutil.copy2(w2k.case+'.clmsum_old', w2k.case+'.clmsum')
        if not os.path.isfile(w2k.case+'.clmsum'):
            print('no '+w2k.case+'.clmsum(_old) file found, which is necessary for lapw0 !')
            print('no '+w2k.case+'.clmsum(_old) file found, which is necessary for lapw0 !', file=fday)
            sys.exit(1)

        if os.path.isfile(w2k.case+'.broyd1'):
            print(w2k.case+'.broyd* files present! I will remove previous broyd files.')
            # print 'You have '+str(p['sleeptime'])+' seconds to kill this job ( ^C   or   kill '+str(os.getpid())+' )'
            # print 'or the script will rm *.broyd* and continue (use -NI to avoid automatic rm)'
            # time.sleep(p['sleeptime'])
            cmd = 'rm -f *.broyd*'
            subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)
            print(w2k.case+'.broyd* files removed!', file=fday)

        print('\n   start'+(' '*8)+ time.asctime()+' with lapw0 ('+str(1)+'/'+str(p['riter'])+' to go)', file=fday)

        # Starting with LAPW
        print('\n   cycle %s \t%s %s/%s to go\n' % (0,time.asctime(),p['finish'],(p['finish'])%p['riter']), file=fday)
        lapw0(fday,w2k.case,w2k.WIENROOT,para, fh_info); fday.flush()
        lapw1(fday,w2k.case,w2k.WIENROOT,para,p['dftKS'],dmfe,wopt,fh_info); fday.flush()
        if wopt['so']: lapwso(fday,w2k.case,w2k.WIENROOT,para,p['dftKS'],dmfe, fh_info)
        fday.flush()

    if (not p['nambu']) and para and (not p['dftKS']): convert_processes.convert(wopt['so'], fh_info)

    #####################
    # Major charge loop #
    #####################
    icycle = nimp = 0
    istep0=0
    EF = 0.0
    dEF = 0.0
    nl_imp = None
    diff_nimp,ddiff_nimp,pdiff_nimp=1,1,0
    (RelaxingStructure, mix_method) = fstc.AreWeRelaxingStructure2(w2k.case)
    print('Relaxing Structure=', RelaxingStructure, file=fh_info)
    sigma_prev_file = 'sig.inp.0'
    if p['nambu']:
        nl_imp = {iat: 0.0 for iat in list(iSigind.keys())}

    while icycle < p['finish']:

        if (icycle>0): # Does not exists yet. Will find it after dmft1
            (nl_imp,nd_imp) = FindLatticeNd(w2k.case, wopt['m_extn'], inl.siginds, imp2latt, fh_info, scf_name='scf2')

        #####################
        # DMFT loop         #
        #####################
        sigmas=[]
        itt = 0
        n_imp = zeros(len(list(iSigind.keys())))

        if 'GoodGuess' in p.p and p['GoodGuess']==True:
            if icycle == 0 :
                p['max_dmft_iterations'] = 0 # do not run impurity at the beginning, but rather converge DMFT charge on existing self-energy
                itt = 1
            if icycle == 1 :
                p['max_dmft_iterations'] = 2 # at the next iteration, rerun impurity twice to have good histogram from the second run
                itt = 0
            if icycle > 1  : p['max_dmft_iterations'] = 1 # finally, continue as normally
            print('We have GoodGuess and max_dmft_iterations=', p['max_dmft_iterations'], 'itt=', itt, file=fh_info)

        if p['nambu']:
            print(':NAMBU_CC convergence_threshold={:.2e}'.format(nambu_cc), file=fh_info)

        while itt < p['max_dmft_iterations']:

            extn = str(icycle+1)+'.'+str(itt+1)  # Extension for the output files
            p.refresh()                          # updating parameters in case user changes them
            if p['nambu']:
                # Re-apply Nambu overrides (refresh re-reads params.dat which has non-Nambu values)
                p['finish'] = 1
                p['recomputeEF'] = 0
                p['max_lda_iterations'] = 0

            # Recomputing the chemical potential, but only if more than one DMFT step. Otherwise EF was recomputed in charge-SC step.
            if itt>1 and p['recomputeEF']:
                if p['Ferro'] and wopt['m_extn']:
                    # In case of ferromagnetic or ferrimagnetic metal, we need to compute common EF for up and dn spins
                    FermiF(dmfe, w2k, p, para, wopt, fh_info, fday, inl.matsubara)
                else:
                    Prepare_dmft2(dmfe, p, para, wopt, fh_info, 1, 'm')
                    SSplit(dmfe.ROOT, inl.matsubara, p['iparams0'], p['ntail'], fh_info)
                    dEF =  dmft2(fday, w2k.case, fh_info, dmfe); fday.flush()

            # Computing Green's function
            if p['runGF']:
                if p['nambu']:
                    # Nambu: use dmftgk.py --nambu instead of dmft1 (unified driver)
                    SSplit(dmfe.ROOT, inl.matsubara, p['iparams0'], p['ntail'], fh_info)
                    # dmftgk is a lattice (k-point) step: launch with MPI2 (mpi_prefix.dat2) like
                    # dmft1/dmft2 and like dmftu (x_dmft.py:606) -- dmftu/dmftgk MUST share one rank
                    # count (Udmft.{rank} is written by dmftu per-rank and read per-rank by dmftgk).
                    # mpi_prefix.dat remains the impurity (boldc) launcher (official convention).
                    # Without mpi_prefix.dat2, MPI2==MPI (utils.py) -> old decks unchanged.
                    cmd_dmftgk = dmfe.MPI2 + ' ' + sys.executable + ' ' + os.path.join(NMBU_DIR, 'dmftgk.py') + ' --nambu -f'
                    print('#<dmftgk --nambu>: ', cmd_dmftgk, file=fh_info)
                    ret_dmftgk = subprocess.call(cmd_dmftgk, shell=True, stdout=fh_info, stderr=fh_info)
                    if ret_dmftgk != 0:
                        print('ERROR: dmftgk.py --nambu exited with code {}'.format(ret_dmftgk), file=fh_info)
                        sys.exit(1)
                    if nl_imp is None:
                        nl_imp = {}
                    for iat in list(iSigind.keys()):
                        nl_imp[iat] = 0.0
                else:
                    SSplit(dmfe.ROOT, inl.matsubara, p['iparams0'], p['ntail'], fh_info)
                    dmft1(fday, w2k.case, fh_info, extn, dmfe)
                    if wopt['m_extn']:
                        SSplit(dmfe.ROOT, inl.matsubara, p['iparams0'], p['ntail'], fh_info, wopt['m_extn'])
                        dmft1(fday, w2k.case, fh_info, extn, dmfe, wopt['m_extn'])
                    (nl_imp,nd_imp) = FindLatticeNd(w2k.case, wopt['m_extn'], inl.siginds, imp2latt, fh_info, scf_name='scfdmft1',cmplx=True)

            SJoin(dmfe.ROOT, fh_info, p, wopt['m_extn'])

            if os.path.isfile('EF.dat'):
                EF = float(open('EF.dat').read())
            elif not p['nambu']:
                open(w2k.case+'.indmf1')
                EF =float(open(w2k.case+'.indmf1').readlines()[1].split()[1])
            else:
                raise ValueError('Nambu mode requires EF.dat for fixed chemical potential')

            # ------------ Impurity solver -------------------
            if p['runIMP']:
                nimp += 1
                UpdateAtom = (p['UpdateAtom']!=0 and mod(nimp,p['UpdateAtom'])==0)
                if nl_imp is None:
                    nl_imp = {iat: 0.0 for iat in list(iSigind.keys())}

                for iat in list(iSigind.keys()):   # Over all inequivalent impurity problems
                    ni=SolveImpurity(EF, solver[iat], iat, extn, p, UpdateAtom, nl_imp[iat], imp2latt[iat][0], fh_info, fh_pinfo, fday)
                    n_imp[iat] += ni
                    print('STOP  IMPURITY_'+str(iat), file=sys.stderr)
                SGather(dmfe.ROOT, fh_info, p, wopt['m_extn'])
                if p['nambu']:
                    cmd_monitor = [sys.executable, os.path.join(NMBU_DIR, 'monitor_sc.py')]
                    print('#<monitor_sc>: ', ' '.join(cmd_monitor), file=fh_info)
                    subprocess.call(cmd_monitor, stdout=fh_info, stderr=fh_info)

                IEorb = 0.0
                XEorb = 0.0
                ZEorb = 0.0
                dat=[]
                for iat in list(iSigind.keys()):   # Over all inequivalent impurity problems
                    if os.path.exists(solver[iat].dir+'/Eorb.dat'):
                        dat += ['--- impurity '+str(iat)+' with degeneracy '+str(len(imp2latt[iat]))+' ---\n']
                        datw = open(solver[iat].dir+'/Eorb.dat','r').readlines()
                        for line in datw:
                            if line[:6]==':IEORB':
                                IEorb += float(line[6:]) * len(imp2latt[iat]) # Takes into account the degeneracy of the impurity
                            elif line[:6]==':XEORB':
                                XEorb += float(line[6:]) * len(imp2latt[iat]) # Takes into account the degeneracy of the impurity
                            elif line[:6]==':ZEORB':
                                ZEorb += float(line[6:]) * len(imp2latt[iat]) # Takes into account the degeneracy of the impurity
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
                if p['nambu']:
                    # Unified sig.inp contains Phi — convergence tracked through sigma_diff
                    sigma_diff = FileDiffMaxAbs(sigma_prev_file, 'sig.inp')
                    sigma_prev_file = 'sig.inp.'+extn
                    LogConvergence(itt+1, sigma_diff, sigma_diff, fh_info)
                    # itt>=1 (i.e. iteration >=2): do NOT test convergence on the FIRST DMFT iteration.
                    # sigma_prev_file is initialised to the SSplit shard 'sig.inp.0' (NOT a snapshot of
                    # the pre-loop sig.inp), so the iter-1 sigma_diff is a spurious 0.0 that would falsely
                    # trip ':CONVERGED' and exit after a single iteration -> Phi never develops for doped/SC
                    # points. From iter 2 onward, sigma_prev_file is the previous iteration's sig.inp.<extn>.
                    if sigma_diff is not None and sigma_diff < nambu_cc and itt >= 1:
                        print(':CONVERGED sigma_diff={:.2e} < {:.2e} at iteration {:d}'.format(sigma_diff, nambu_cc, itt+1), file=fh_info)
                        itt += 1
                        break

            print('------- DMFT part of the step # ', itt, 'done ---------\n', file=fh_info)
            fh_info.flush()

            itt += 1   # END DMFT INNER LOOP

        if itt > 0:
            n_imp *= 1./itt

        if int(len(sigmas)*p['saver'])>1 :
            sigmas.reverse()
            sigmas = sigmas[0:builtins.max(1,int(len(sigmas)*p['saver']))]
            print('Averaging over the following self-energies', sigmas, file=fh_info)
            cmd = 'saverage.py -o sig.inp ' + (' '.join(sigmas))
            print('#<saverage>: ', cmd, file=fh_info)
            subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)

        if p['nambu']:
            icycle += 1
            break

        (RelaxingStructure, mix_method) = fstc.AreWeRelaxingStructure2(w2k.case)
        if RelaxingStructure:
            fstc.TOTtoFOR(w2k.case)
            p['max_lda_iterations']=p['max_lda_iterations_optimize']
            print('Mixing method corresponds to structural optimization.', file=fh_info)
            print(' changing TOT to FOR in case.in2 files', file=fh_info)
            print(' increasing max_lda_iterations to ', p['max_lda_iterations'], file=fh_info)

        diff_nimp=1
        for ita in range(p['max_lda_iterations']):
            if ita%2 == 1 : p.refresh()

            recomputeEF = p['recomputeEF']

            if p['Ferro'] and wopt['m_extn'] and recomputeEF:
                # In case of ferromagnetic or ferrimagnetic metal, we need to compute common EF for up and dn spins
                FermiF(dmfe, w2k, p, para, wopt, fh_info, fday, inl.matsubara)
                recomputeEF = 0

            Prepare_dmft2(dmfe, p, para, wopt, fh_info, recomputeEF, 'c')
            SSplit(dmfe.ROOT, inl.matsubara, p['iparams0'], p['ntail'], fh_info)
            dEF = dmft2(fday,w2k.case,fh_info,dmfe); fday.flush()
            if wopt['m_extn']:
                SSplit(dmfe.ROOT, inl.matsubara, p['iparams0'], p['ntail'], fh_info, wopt['m_extn'])
                dmft2(fday,w2k.case,fh_info,dmfe,wopt['m_extn'])
                fday.flush()
                combineud(w2k.case, dmfe.ROOT, wopt, fh_info)


            if recomputeEF==0: # if we do not recomputeEF we need to make sure that density is correct, otherwise we have a doped compound!
                with open(w2k.case+'.outputdmf2') as f2:
                    lines = f2.readlines()
                drho=0
                for line in lines[::-1]:
                    if 'rho-rho_expected=' in line:
                        drho = float(line.split()[-1])
                        break
                print('rho-rho_expected=', drho, file=fh_info)
                if abs(drho)>0.1:
                    print('WARNING/ERROR: Since recomputeEF is set to 0, we expect EF in the gap.', file=fh_info)
                    print(' But density is wrong. Have rho-rho_expected=',drho, file=fh_info )
                    print(' You should adjust EF in EF.dat and see if EF falls into the gap', file=fh_info)

            nd,Vdc = CreateEorb(w2k.case, wopt, inl.siginds, imp2latt)

            nc = ConvertFromLatticeToImp(nd, imp2latt)

            lcore(fday,w2k.case,w2k.WIENROOT,wopt,fh_info);       fday.flush()
            scf  (w2k.case, wopt)
            scf1 (fday,w2k.case,w2k.WIENROOT);
            mixer(fday,w2k.case,w2k.WIENROOT,fh_info);       fday.flush()
            scfm (fday,w2k.case,w2k.WIENROOT);               fday.flush()
            drho, dene = Diff(fday, w2k.case, dEF, diff_nimp, ddiff_nimp);fday.flush()

            # Writting to info.iterate file
            if os.path.isfile('EF.dat'): EF = float(open('EF.dat').read())
            (ETOT,SUM,XSUM,YSUM,ZSUM,EORB,XEORB,ZEORB) = FindEtot(w2k.case,wopt)

            fh_pinfo.write( '%3d %3s.%4s %12.6f %12.6f %15.6f %15.6f %15.6f ' % (istep0, icycle, ita, EF, Vdc, ETOT, ETOT-SUM+ZSUM-EORB+ZEORB, ETOT-SUM+YSUM-EORB+XEORB) )
            diff_nimp=0
            for iat in list(iSigind.keys()):   # Over all inequivalent impurity problems
                if os.path.exists('imp.'+str(iat)+'/Eimp.inp'):
                    fEimp = open('imp.'+str(iat)+'/Eimp.inp', 'r')
                    (Eimp,Olap,Edc) = [array([float(x) for x in line.split('#')[0].split()]) for line in fEimp.readlines()]  # throw out comments
                else: # If runGF=False, this is not available
                    Eimp=array([0,0])
                fh_pinfo.write( '%12.6f %12.6f %12.6f %12.6f ' % (nc[iat], n_imp[iat], Eimp[0], Eimp[-1]) )
                diff_nimp += n_imp[iat]-nc[iat]
            fh_pinfo.write('\n')
            fh_pinfo.flush()
            istep0 += 1
            print('------- LDA(charge) step', ita, icycle, 'done ---------', file=fh_info)

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

            print('\n:ITT%-2d   cycle %s \t%s %s/%s to go\n' % (istep2, istep1,time.asctime(),p['finish']-istep1,(p['finish']-istep1)%p['riter']), file=fday)

            lapw0(fday,w2k.case,w2k.WIENROOT,para, fh_info); fday.flush()
            lapw1(fday,w2k.case,w2k.WIENROOT,para,p['dftKS'],dmfe,wopt,fh_info); fday.flush()
            if wopt['so']: lapwso(fday,w2k.case,w2k.WIENROOT,para,p['dftKS'],dmfe, fh_info)

            fday.flush()

            if succi : break
            if ita>0 and ita%p['remove_broyd']==0 :
                for f in glob.glob('*.broyd*'):
                    os.remove(f)

        ddiff_nimp = abs(pdiff_nimp-diff_nimp)
        pdiff_nimp = diff_nimp

        drho, dene = Diff(fday, w2k.case, dEF, diff_nimp, ddiff_nimp);fday.flush()

        icycle += 1   # END CHARGE LOOP
        if RelaxingStructure: # When dmft/impurity loop restarts, we should relax again
            fstc.RelaxAgain(w2k.case, mix_method)
            shutil.copy2(w2k.case+'.struct',w2k.case+'.struct.'+str(icycle-1))

        if (p['max_lda_iterations']>=10):
            toclean = glob.glob('*.broyd*')
            for f in toclean: os.remove(f)

        if not RelaxingStructure and ddiff_nimp<p['nc'] and drho<2*p['cc'] and dene<2*p['ec']:
            # impurity occupation is not changing much anymore, and charge and energy is converged, hence we can safely stop if structure is not relaxing.
            break

    toclean = glob.glob('dmft[1|2].error*')+glob.glob('*.broyd*')+glob.glob('dmft')+glob.glob('dmft2')+glob.glob('_processes*')+glob.glob('*.klist_*')+glob.glob('lapw1')+glob.glob('*.output1_*')+glob.glob('*.vector_*')+glob.glob('*.scf1_*')+glob.glob('lapw1_*.[def|error]')
    for f in toclean: os.remove(f)

    allFiles=glob.glob("*")
    for filename in allFiles:
        if os.stat(filename).st_size==0:
            os.remove(filename)
