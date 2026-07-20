#!/usr/bin/env python
import sys, subprocess
import re, os, glob, shutil, socket, time, json
from os.path import getsize
from numpy import *


def _json_ready(value):
    if isinstance(value, dict):
        return {str(k): _json_ready(v) for (k, v) in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(v) for v in value]
    if isinstance(value, ndarray):
        return [_json_ready(v) for v in value.tolist()]
    if isinstance(value, (bool, str)) or value is None:
        return value
    if isinstance(value, integer):
        return int(value)
    if isinstance(value, floating):
        return float(value)
    if isinstance(value, complex):
        return [float(value.real), float(value.imag)]
    return value


def _solver_name(par):
    if 'solver' in par:
        return str(par['solver']).upper()
    if 'impurity_solver' in par:
        return str(par['impurity_solver']).upper()
    print('ERROR: impurity solver is not specified in tbg_params.py.')
    print('Please set par["solver"] to "CTQMC" or "BCTQMC".')
    sys.exit(1)


def _iparam_value(iparams, key, default=None):
    if key in iparams:
        return iparams[key][0]
    return default


def _solver_output_names(solver, iparams, log=sys.stdout, warn_alias=True):
    if solver == 'CTQMC':
        return (_iparam_value(iparams, 'Sig', 'Sig.out'),
                _iparam_value(iparams, 'Gf', 'Gf.out'))
    if solver == 'BCTQMC':
        fsig = _iparam_value(iparams, 'fSig', None)
        fgf = _iparam_value(iparams, 'fGf', None)
        if fsig is None:
            fsig = _iparam_value(iparams, 'Sig', 'Sig.out')
        elif warn_alias and 'Sig' in iparams:
            print('WARNING: BCTQMC parameter "fSig" already set; ignoring alias "Sig"', file=log)
        if fgf is None:
            fgf = _iparam_value(iparams, 'Gf', 'Gf.out')
        elif warn_alias and 'Gf' in iparams:
            print('WARNING: BCTQMC parameter "fGf" already set; ignoring alias "Gf"', file=log)
        return (fsig, fgf)
    print('ERROR: Unknown impurity solver '+solver, file=log)
    sys.exit(1)


def write_ctqmc_params(dire, iparams, mu_QMC):
    with open(dire+'/PARAMS', 'w') as fp:
        for p in list(iparams.keys()):
            print(p, '  ', iparams[p][0], iparams[p][1], file=fp)
        print('mu    ', mu_QMC, file=fp)


def write_bctqmc_params(dire, iparams, mu_QMC, log=sys.stdout):
    aliases = {'M': 'Ms', 'tsample': 'tmeassure', 'Sig': 'fSig', 'Gf': 'fGf'}
    ignored = set(['mode', 'nom', 'GlobalFlip', 'svd_lmax', 'svd_x0', 'minM',
                   'minF', 'nf0', 'cix'])
    params = {}
    for key in list(iparams.keys()):
        if key in ignored:
            print('WARNING: BCTQMC ignores iparams0["'+key+'"]', file=log)
            continue
        if key in aliases:
            continue
        params[key] = iparams[key][0]

    for key in list(iparams.keys()):
        if key not in aliases:
            continue
        bkey = aliases[key]
        if bkey in params:
            print('WARNING: BCTQMC parameter "'+bkey+'" already set; ignoring alias "'+key+'"', file=log)
            continue
        params[bkey] = iparams[key][0]
        print('BCTQMC parameter alias: '+key+' -> '+bkey, file=log)

    params['mu'] = mu_QMC
    params['cix'] = 'bctqmc.cix'
    params.setdefault('Delta', 'Delta.inp')
    params['fSig'], params['fGf'] = _solver_output_names('BCTQMC', iparams, log, warn_alias=False)
    with open(dire+'/PARAMS', 'w') as fp:
        json.dump(_json_ready(params), fp, indent=2)
        fp.write('\n')


if __name__ == '__main__':
    exec(open("tbg_params.py").read())
    dirs = par['dirs']
    log = open('info.iterate','w')

    ROOT = '.'
    #ROOT = os.environ['TBG_ROOT']
    ROOT = os.environ['WIEN_DMFT_ROOT']
    
    for itt in range(par['Nitt']):
        exec(open("tbg_params.py").read())  # updating parameters if changed
        dirs = par['dirs']
        solver = _solver_name(par)
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
        
        # Preparing impurity solver input file
        for idi in range(len(dirs)):
            mu_QMC = -Eds[idi][0] + Sig_DC
            if solver == 'CTQMC':
                write_ctqmc_params(dirs[idi], iparams0, mu_QMC)
            elif solver == 'BCTQMC':
                write_bctqmc_params(dirs[idi], iparams0, mu_QMC, log)
            else:
                print('ERROR: Unknown impurity solver '+solver, file=log)
                sys.exit(1)
        
        if solver == 'CTQMC':
            sname = ROOT+'/ctqmc'
        elif solver == 'BCTQMC':
            sname = ROOT+'/bctqmc'
        fsig, fgf = _solver_output_names(solver, iparams0, log, warn_alias=False)
        print(time.strftime("%a %b %d %H:%M:%S %Z %Y")+'>     '+sname, file=log)
        for idi in range(0,1): #len(dirs)):
            cmd = 'cd '+dirs[idi]+';'+mpi_prefix+' '+sname+' >& nohup.dat '
            print(('#<'+sname+'>: '), cmd, file=log)
            log.flush()
            info=subprocess.call(cmd,shell=True,stdout=log)
        
            Sigdata = loadtxt(dirs[idi]+'/'+fsig)
            for i in range(1,Sigdata.shape[1],2):
                Sigdata[:,i] -= Sig_DC
            savetxt(dirs[idi]+'/'+par['Sigma'], Sigdata)
            
            shutil.copy2(dirs[idi]+'/'+fsig, dirs[idi]+'/Sig.out.'+str(itt))
            shutil.copy2(dirs[idi]+'/'+fgf, dirs[idi]+'/Gf.out.'+str(itt))
            shutil.copy2(dirs[idi]+'/'+par['Sigma'], dirs[idi]+'/'+par['Sigma']+'.'+str(itt))
            
