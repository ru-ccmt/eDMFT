#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import sys, subprocess
import shlex, tempfile
from pathlib import Path
import re, os, glob, shutil, socket, time
from os.path import getsize
from scipy import *
import numpy
import copy
import glob
from datetime import date,datetime
from utils import shellcmd

def Step(cmd, fday):
    dt = datetime.now().strftime("%B %d  %H:%M:%S %Y")
    with open(':log', 'a') as logs:
        logs.write(dt+' '+cmd+'\n')
    #print(cmd)
    print('>%-10s' % cmd.split()[1],  '( '+dt+' )', file=fday)
    fday.flush()
    #subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)
    (out, err, retit) = shellcmd(cmd)
    print(out.decode('UTF-8'), file=sys.stdout, end='')
    print(err.decode('UTF-8'), file=sys.stderr, end='')
    if retit!=0:
        sys.exit(1)
    
def Diff(case, fday):
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
    # Equal printed SCF energies are a valid zero change at output precision.
    # Keep it zero: the caller also requires the charge residual to converge.
    # Insufficient energy history is handled by the fallback above.
    print(':ENERGY convergence: '+str(dene))
    print(':CHARGE convergence: '+str(drho))
    print(':ENERGY convergence: '+str(dene),file=fday)
    print(':CHARGE convergence: '+str(drho),file=fday)
    fday.flush()
    return (drho, dene)

def CheckError(step):
    if os.path.exists(step+'.error') and os.path.getsize(step+'.error')>0 :
        message = open(step+'.error','r').readlines()
        print(message)
        return True
    else:
        return False

def spinor_scf(up_file, down_file):
    """Trace the charge and vector observables written by nonmagnetic LAPW2."""
    # Ordinal atom labels and field widths are part of the WIEN SCF interface.
    # Band energies, entropy and occupations refer to a *single* common SOC
    # spectrum. Keep those records once. Match local quantities by atom tag.
    rows = up_file.read_text().splitlines(keepends=True)
    local_tag = re.compile(r'^:(CHA(?:\d{3}|  )|QTL\d{3}|EP[LH]\d{3}|VZZ\d{3}|F(?:VA|SU)\d{3}|STR_VAL\d{3}):')
    down = {}
    for line in down_file.read_text().splitlines(keepends=True):
        tag = local_tag.match(line)
        if tag:
            if tag[1] in down:
                raise ValueError(f'Duplicate LAPW2 record {tag[1]} in {down_file}')
            down[tag[1]] = line
    if 'CHA  ' not in down:
        raise ValueError(f'Missing total valence charge in {down_file}')
    for i, line in enumerate(rows):
        match = local_tag.match(line)
        if not match:
            continue
        key = match[1]
        other = down.pop(key)
        if line == other:
            continue
        prefix = match[0]
        if key.startswith('CHA') or (key.startswith('VZZ') and '=' in line):
            # The first '=' precedes the charge. Preserve RMT and atom labels.
            left, field = line.split('=', 1)
            a = re.match(r'(\s*[-+\d.]+)(.*)', field.rstrip('\n'))
            b = float(other.split('=', 1)[1].split()[0])
            v = (float(a[1]) + b) / 2
            prec = len(a[1].split('.')[1])
            rows[i] = f'{left}={v:{len(a[1])}.{prec}f}{a[2]}\n'
            continue
        fields = line[len(prefix):].split()
        fields_dn = other[len(prefix):].split()
        if key.startswith(('FVA', 'FSU')):
            if fields[0] != fields_dn[0]:
                raise ValueError(f'Inconsistent LAPW2 atom labels for {key}')
            force = (numpy.array(fields[2:], float) + numpy.array(fields_dn[2:], float)) / 2
            # Average the Cartesian components first: |(Fu+Fd)/2| can be zero
            # even if both one-component force magnitudes are large.
            rows[i] = prefix + f'{fields[0]:>9}' + ''.join(
                f'{v:20.12f}' for v in (numpy.linalg.norm(force), *force)) + '\n'
            continue
        a, b = numpy.array(fields, float), numpy.array(fields_dn, float)
        avg = (a + b) / 2
        if key.startswith(('EPL', 'EPH')):
            # LAPW1 can use the partial-charge energy centers. Combine their
            # first moments, then divide by the combined partial charge.
            q = a[::2] + b[::2]
            e = numpy.divide(a[::2]*a[1::2] + b[::2]*b[1::2], q,
                             out=avg[1::2].copy(), where=q != 0)
            body = '  '.join(f'{w:8.4f}{v:8.4f}' for w, v in zip(q/2, e))
        elif key.startswith('QTL'):
            body = ''.join(f'{v:7.4f}' for v in avg)
        elif key.startswith('VZZ'):
            body = ' '*8 + ''.join(f'{v:12.5f}' for v in avg[:-1]) + f'{avg[-1]:12.3f}'
        else:  # STR_VAL: Cartesian components of one tensor row
            body = '  ' + ''.join(f'{v:20.12f}' for v in avg)
        rows[i] = prefix + body + '\n'
    if down:
        raise ValueError(f'Unpaired lower-spinor SCF records: {sorted(down)}')
    return ''.join(rows)


def stage_soc_in1(in1, merged_scf, work):
    """Stage the WIEN2k global linearization-energy update for the next mixer.

    The two direct lapw2c calls bypass the postprocessing in x lapw2.
    Both spinors share one Fermi level. Use that level once, after the spin
    trace, and leave every radial, cutoff and band-selection record intact.
    """
    header, newline, basis = in1.read_bytes().partition(b'\n')
    if not newline:
        raise ValueError(f'Incomplete LAPW1 input in {in1}')
    # Older WIEN inputs may omit EF entirely. LAPW1 then uses 0.5 Ry.
    field = re.search(rb'\bEF\s*=\s*([^\s,]+)', header)
    old = float(field[1].replace(b'D', b'E').replace(b'd', b'e')) if field else 0.5
    records = [row for row in merged_scf.read_text().splitlines() if row.startswith(':FER  :')]
    if len(records) != 1 or '=' not in records[0]:
        raise ValueError(f'Expected one common :FER record in {merged_scf}')
    raw = records[0].split('=', 1)[1].split()[0]
    # Keep the x lapw2 fallback for a native output field full of asterisks.
    ef = 1.5 if raw == '**************' else float(raw.replace('D', 'E').replace('d', 'e'))
    if not numpy.isfinite([old, ef]).all():
        raise ValueError(f'Non-finite linearization reference in {in1} or {merged_scf}')

    # Normal x lapw2 step: mix two thirds of the change, bounded by 0.6 Ry.
    # MIXER reads the three numbers in order: proposed, previous, raw EF.
    mixed = old + numpy.clip((ef - old) * (2.0 / 3.0), -0.6, 0.6)
    values = f'  EF={mixed:.13f} Old {old:.13f} {ef:.13f}'.encode('ascii')
    updated = work / in1.name
    updated.write_bytes(header[:5] + values + b'  (WFFIL, WFPRI, ENFIL, SUPWF)\n' + basis)
    shutil.copymode(in1, updated)
    return updated


def complex_soc_lapw2(case, wienroot, dmftroot, fday):
    """Use both components of the same nonmagnetic, complex SOC eigensystem."""
    # In a complex cell the two spinor components at a representative k
    # need not have the same spatial density or force. LAPW2 takes one vector
    # file per invocation. Its nonmagnetic convention gives each pass the
    # full scalar degeneracy, hence the spin trace is the arithmetic mean.
    # Both passes must keep the *same* occupation, scalar potential, and
    # nonmagnetic normso convention (selecting -dn changes those inputs).
    Step(shlex.join([str(Path(wienroot)/'x'), 'lapw2', '-d', '-c', '-so', '-f', case]), fday)
    lines = Path('lapw2.def').read_text().splitlines(keepends=True)
    with tempfile.TemporaryDirectory(prefix='soc2_', dir='.') as stage:
        work = Path(stage)
        scf_files = []
        for suffix in ('', 'dn'):
            scf_file = work / (case + '.scf2' + suffix)
            scf_files.append(scf_file)
            names = {8: work / (case + '.clmval' + suffix), 21: scf_file}
            definition = work / ('lapw2' + suffix + '.def')
            updated = []
            for line in lines:
                unit, filename, rest = line.split("'", 2)
                number = int(unit.strip().rstrip(','))
                if number in (10, 30):
                    filename += suffix
                    if not Path(filename).is_file() or Path(filename).stat().st_size == 0:
                        raise FileNotFoundError(f'Missing SOC input {filename}')
                elif number in names:
                    filename = str(names[number])
                elif number == 6 and suffix:
                    filename += suffix  # keep the two detailed LAPW2 logs
                updated.append(f"{unit}'{filename}'{rest}")
            definition.write_text(''.join(updated))
            Step(shlex.join([str(Path(wienroot)/'lapw2c'), str(definition)]), fday)
            if CheckError(str(definition.with_suffix(''))):
                raise RuntimeError(f'LAPW2 failed for {definition.stem}')
        merged = work / 'spin_trace.scf2'
        merged.write_text(spinor_scf(*scf_files))
        (work / (case + '.struct')).symlink_to(Path(case + '.struct').resolve())
        subprocess.run([str(Path(dmftroot)/'combineud'), case], cwd=work, check=True)
        charge = work / (case + '.clmval_aver')
        if not charge.stat().st_size:
            raise RuntimeError('Empty SOC spin-traced valence charge')
        # Publish only after both native runs and the charge/SCF merges succeed.
        in1 = Path(case + '.in1c')
        updated_in1 = stage_soc_in1(in1, merged, work)
        charge.replace(case + '.clmval')
        merged.replace(case + '.scf2')
        updated_in1.replace(in1)


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
            print("WARNING : Electron charge density is very different than expected. You should change recomputeEF and set it to 1 (switch in on)")
            
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
            
    # Pure DFT must not inherit a stale DMFT impurity energy from a copied case.
        
    scfs = ['1s','2s','c']
    if wopt['updn']: scfs = ['cup','cdn']
    for i in scfs:
        if os.path.exists(case+'.scf'+i):
            dat += '\n------- '+case+'.scf'+i+' --------\n'
            dat += open(case+'.scf'+i,'r').readlines()

    with open(case+outfile, 'a') as fs:
        fs.writelines(dat)

def scfm(case):
    dat=[]
    for i in ['m']:
        dat += '\n------- '+case+'.scf'+i+' --------\n'
        if os.path.exists(case+'.scf'+i):
            dat += open(case+'.scf'+i,'r').readlines()
            
    with open(case+'.scf', 'a') as fs:
        fs.writelines(dat)
    
def scf1(case):
    for i in ['clmsum','vsp','vns','vrespsum','vspdn','vnsdn','vspup','vnsup','clmup','clmdn']:
        name = case+'.'+i
        if os.path.exists(name) and os.path.getsize(name)>0:
            shutil.copy2(name, name+'_old')		#save last cycle

if __name__ == '__main__':
    usage = """usage: %run_dft.py [ options ]
    Runs wien2k and dmft internal lapw1&lapwso to converge dft.
    Uses dmft mpi paralelization, which is hopefully more portable.
    """
    import optparse
    parser = optparse.OptionParser(usage)
    parser.add_option('-N', '--Nitt',   dest='Nitt', type='int', default=1000, help="maximum number of dft iterations")
    parser.add_option('-c', '--cc',     dest='cc', type='float', default=1e-3, help="charge convergence criteria")
    parser.add_option('-e', '--ec',     dest='ec', type='float', default=1e-4, help="energy convergence criteria")
    # Next, parse the arguments
    (options, args) = parser.parse_args()
    
    #p={
    #    'cc' : 1e-3,      # the charge density precision to stop the LDA+DMFT run
    #    'ec' : 1e-4,      # the energy precision to stop the LDA+DMFT run
    #    'Nitt' : 1000     # total number of iterations
    #}

    wopt={'cmplx':'', 'so':'', 'm_extn':'', 'updn':''}  
    # cmplx : no-inversion, 
    # so    : spin-orbit coupling present
    # m_extn: dmft is spin polarized
    # updn  : core-is-spin-polarized
    
    rt = os.environ['WIEN_DMFT_ROOT']
    rw = os.environ['WIENROOT']
    
    MPI = ''
    if os.path.exists('mpi_prefix.dat'):
        MPI = open('mpi_prefix.dat', 'r').readline().strip()
    
    #fh_info = open('dft_info.out', 'w')

    #sys.stdout = fh_info
    #fh_info = sys.stdout
    
    case=''
    for name in glob.glob('*.struct'):
        case0 = os.path.splitext(name)[0]
        #print(('case0='+case0))
        if (os.path.isfile(case0+'.clmsum')):
            case = case0
            #print(('case='+case))
    if not case:
        print('Could not find case, hence exiting')
        sys.exit(1)
    
    # Info file
    fday = open(case+'.dayfile', 'w')
    print('Calculating '+case+' in '+os.getcwd()+'\n'+'on '+socket.gethostname()+' with PID '+str(os.getpid()), file=fday)
    print('\n\n', file=fday)
    
    cc, so = '', False
    # Checking noncentrosymmetric
    if os.path.isfile(case+".in1c") and os.path.getsize(case+".in1c")>0 :
        print('Found '+case+'.in1c file, hence assuming non-centrosymmetric structure. Switching -c switch!', file=fday)
        cc = ' -c'
        wopt['cmplx']='c'
    # Check spin-orbit coupling
    if os.path.isfile(case+".inso") and os.path.getsize(case+".inso")>0 :
        print('Found '+case+'.inso file, hence assuming so-coupling exists. Switching -so switch!', file=fday)
        so = True
        wopt['so']='so'

    for i in range(options.Nitt):
        Step(rw+'/x lapw0 -f '+case, fday)
        if CheckError('lapw0'): break
        Step(rt+'/x_dmft.py lapw1', fday)
        if CheckError('lapw1'): break
        if so:
            Step(rt+'/x_dmft.py lapwso', fday)
            if CheckError('lapwso'): break
            if cc:
                complex_soc_lapw2(case, rw, rt, fday)
            else:
                Step(rw+'/x lapw2 -so -f '+case, fday)
        else:
            Step(rw+'/x lapw2'+cc+' -f '+case, fday)
        if not (so and cc) and CheckError('lapw2'): break
        Step(rw+'/x lcore -f '+case, fday)
        if CheckError('lcore'): break
        
        scf(case, wopt)
        scf1(case)
        
        Step(rw+'/x mixer -f '+case, fday)
        if CheckError('mixer'): break
        
        scfm(case)
        drho,dene = Diff(case, fday)
        succi = (drho<options.cc and dene<options.ec)
        if succi:
            # MIXER has just saved the density used by this last LAPW cycle as
            # case.clmsum_old.  The DFT response helper uses that density in
            # an isolated snapshot, so its stress is paired with this cycle's
            # LAPW1 vectors and LAPW2 valence stress.  Never pay for the
            # screened response during unconverged charge iterations.
            in2 = case + '.in2' + wopt['cmplx']
            if open(in2).readline().split()[0].upper() in ('STRT', 'STRFT'):
                Step(sys.executable+' '+rt+'/dft_stress_step.py --case '+case, fday)
            break
        
