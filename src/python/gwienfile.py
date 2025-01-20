# @Copyright 2020 Kristjan Haule
'''
Classes to handle reading of WIEN2K files.
'''
import numpy as np
import glob, os, sys
from numpy import zeros,arange,array

def Read_energy_file(filename, nat, fout, give_kname=False, give_bands=True, Extended=False):
    with open(filename, 'r') as fi:
        lines = fi.readlines()
    iln = 2*nat  # jumping for linearization energies
    if give_kname: knames=[]
    if give_bands: Ebnd=[]
    klist=[]
    wegh=[]
    hsrws=[]
    for ik in range(1000000):
        if iln >= len(lines): break
        line = lines[iln]; iln += 1
        if Extended:
            kks = [line[27*i:27*(i+1)] for i in range(3)]
            kname = line[81:91]
            hsr = [line[91+6*i:91+6*(i+1)] for i in range(2)]
            ws = line[103:108]
        else:
            kks = [line[19*i:19*(i+1)] for i in range(3)]
            kname = line[57:67]
            hsr = [line[67+6*i:67+6*(i+1)] for i in range(2)]
            ws = line[79:84]
        hsrows, Ne = int(hsr[0]), int(hsr[1])
        wgh = float(ws)
        klist.append([float(kx) for kx in kks])     # kpoint
        if give_kname: knames.append(kname)
        wegh.append(wgh)     # wk2
        hsrws.append(hsrows) # ngkir
        if give_bands:
            Ebands=[]
            for ib in range(Ne):
                line = lines[iln]; iln += 1
                dd = line.split()
                ibp1, ee = int(dd[0]), float(dd[1])
                Ebands.append(ee)
                #print ib+1, ibp1, ee
            Ebnd.append( Ebands )  # bande
        else:
            iln += Ne
    if give_kname and give_bands:
        return (klist, wegh, Ebnd, hsrws, knames)
    elif (not give_kname) and give_bands:
        return (klist, wegh, Ebnd, hsrws)
    elif give_kname and (not give_bands):
        return (klist, wegh, hsrws, knames)
    else:
        return (klist, wegh, hsrws)

def get_nat(case):
    with open(case + '.struct', 'r') as f:
        title = next(f).strip()
        second_line = next(f)
    lattice = second_line[0:4].strip()
    nat = int(second_line[27:30])
    return nat
def get_case():
    # directory name is not case (happens when submitting to cluster)
    case = None
    files = glob.glob('*.struct')
    if len(files) < 1:
        raise Exception('No struct file present.')
    elif len(files) > 1:
        # heuristic algorithm to determine case:
        # need in0 and struct present with the same case.
        candidates = [os.path.splitext(f)[0] for f in files]
        for cand in candidates:
            if os.path.isfile(cand+'.in0'):
                case = cand
                break
    else: # just one candidate exists, hence it must be it
        case, ext = os.path.splitext(os.path.basename(files[0]))
    return case
    
if __name__ == '__main__':
    case = get_case()
    print('case=', case)
    nat = get_nat(case)
    print('nat=', nat)

    so = 'so' if (os.path.isfile(case+".inso") and os.path.getsize(case+".inso")>0) else ''
    m_ext = ''
    
    energy_file = case+'.energy'+so+m_ext
    (klist, wgh, hsrws) = Read_energy_file(energy_file, nat, sys.stdout, give_kname=False, give_bands=False, Extended=False)
    print('klist=', klist)
    print('wgh=', wgh)
    print('hsrws=', hsrws)
