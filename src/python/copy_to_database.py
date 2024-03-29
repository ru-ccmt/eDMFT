#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
from scipy import *
import sys, re
from scipy import linalg
import optparse
import glob, os, shutil
from utils import W2kEnvironment, Ry_in_eV

def findlast(name):
    cdos = glob.glob(name+'.*')
    cdic={}
    for fil in cdos:
        m = re.match(name+'\.(\d+)\.(\d+)', fil)
        if m is not None:
            pair = m.groups()
            clast = '.'+pair[0]+'.'+pair[1]
            ckey = int(pair[0])*10000+int(pair[1])
            cdic[ckey]=clast
    sk=sorted(list(cdic.keys()))
    last = cdic[sk[-1]]
    prelast = cdic[sk[-2]]
    return (last,prelast)

if __name__ == '__main__':
    """ Will have some help
    """
    
    usage = """usage: %prog <destination-dir> -o <results-dir> -d <DOS-dir> [ options ]
    """
    parser = optparse.OptionParser(usage)
    parser.add_option("-o", "--out", dest="out", default='output', help="directory with the output of the run")
    parser.add_option("-d", "--dos", dest="DOS", default=None,    help="directory with the output DOS on the real axis")
    # Next, parse the arguments
    (options, args) = parser.parse_args()

    print('options=', options)
    print('args=', args)
    if len(args)==0:
        print('Give name of the output directory!')
        sys.exit(0)
    else:
        cdir = args[0]
        if not os.path.exists(cdir):
            print('Directory',cdir,'does not exist!')
            sys.exit(0)
    
    # Finds what is case
    w2k = W2kEnvironment()
    case=w2k.case

    wimp = 'imp.0'
    root = cdir+'/'+case
    rout = root+'/output'
    rimp = root+'/output/'+wimp
    if not os.path.exists(root):
        os.mkdir(root)
    if not os.path.exists(rout):
        os.mkdir(rout)
    if not os.path.exists(rimp):
        os.mkdir(rimp)
    
    w0in = [case+'.struct',case+'.in0',case+'.clmsum', case+'.inm', case+'.in1', case+'.klist', case+'.inso', case+'.in2', case+'.in2c', case+'.kgen', case+'.inc', case+'.scf2']
    w1in = ['params.dat', 'sig.inp', case+'.indmfl', case+'.indmfldn', case+'.indmfi', 'Sigma.000']
    
    w0out = [case+'.scf',case+'.dayfile',case+'.rotlm',':log']
    w1out = ['info.iterate','dmft0_info.out','dmft1_info.out','dmft2_info.out',case+'.outputdmf1.0',case+'.outputdmf2.0',case+'.Eimp1','Edc.dat','Eorb.dat','EF.dat']
    w2out = [case+'.gc1',case+'.dlt1',case+'.cdos']

    out = options.out.strip("/")
    last,prelast = findlast(out+'/'+case+'.cdos')
    print('last=', last)
    w3out = [case+'.gc1'+last,case+'.dlt1'+last,case+'.cdos'+last,'sig.inp'+last]
    w4out = [wimp+'/actqmc.cix', wimp+'/Gf.out', wimp+'/Delta.inp', wimp+'/PARAMS', wimp+'/Sig.out', wimp+'/nohup_imp.out'+prelast, wimp+'/Probability.dat']

    for f in w0in+w1in:
        fle = glob.glob(f)
        if fle:
            print('Input  -- Copying ... ', fle[0], 'to', root)
            shutil.copy2(fle[0], root)

    for f in w0in+w1in+w0out+w1out+w2out+w3out:
        fle = glob.glob(options.out+'/'+f)
        if fle:
            print('Output -- Copying ... ', fle[0], 'to', rout)
            shutil.copy2(fle[0], rout)

    for f in w4out:
        fle = glob.glob(options.out+'/'+f)
        if fle:
            print('Output -- Copying ... ', fle[0], 'to', rimp)
            shutil.copy2(fle[0], rimp)
            
    if options.DOS is not None:
        for f in w2out:
            fle = glob.glob(options.DOS+'/'+f)
            if fle:
                print('RealDOS-- Copying ... ', fle[0], 'to', rout)
                shutil.copy2(fle[0], rout)
        fle = glob.glob(options.DOS+'/sig.inp')
        if fle:
            print('RealDOS-- Copying ... ', fle[0], 'to', rout+'/sig.inp_real')
            shutil.copy2(fle[0], rout+'/sig.inp_real')
            
