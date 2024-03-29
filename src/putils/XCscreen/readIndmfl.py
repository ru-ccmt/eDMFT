# @Copyright 2007 Kristjan Haule

from scipy import *
import sys

def ReadIndmfl(filename, fh_info):
    """Read the self-energy index file Sigind and the local transformation matrix CF from a file"""
    def divmodulo(x,n):
        "We want to take modulo and divide in fortran way, so that it is compatible with fortran code"
        return ( sign(x)* (abs(x)/n) , sign(x)*mod(abs(x),n))

    fh = open(filename, 'r')
    lines = [line.split('#')[0].strip() for line in fh.readlines()] # strip comments
    lines = (line for line in lines if line)  # strip blank lines & create generator expression

    hybr_emin, hybr_emax, Qrenorm, projector = [float(x) for x in lines.next().split()[:4]]
    if projector>=4:
        hybr_emin = int(hybr_emin)
        hybr_emax = int(hybr_emax)
    matsubara, broadc, broadnc, om_npts, om_emin, om_emax = [float(e) for e in lines.next().split()[:6]]
    matsubara = int(matsubara)  # recast these to integers
    om_npts   = int(om_npts) 

    atoms={}
    cps={}
    natom = int(next(lines))
    for i in range(natom):
        iatom, nL, locrot_shift = [int(x) for x in lines.next().split()]
        (shift,locrot) = divmodulo(locrot_shift,3)
        if locrot<0: locrot=3
        
        Ls, qsplits, icps = array([[int(x) for x in lines.next().split()] for i in range(nL)]).T
        new_zx = [[float(x) for x in lines.next().split()] for loro in range(abs(locrot))]
        vec_shift = [float(x) for x in lines.next().split()] if shift else None

        atoms[iatom] = (locrot, new_zx, vec_shift)
        for icp, L, qsplit in zip(icps, Ls, qsplits):
            if icp in cps:
                cps[icp] += [(iatom, L, qsplit)]
            else:
                cps[icp] = [(iatom, L, qsplit)]

    #####################################################
    # read the big block of siginds and cftrans
    ncp, maxdim, maxsize = [int(e) for e in lines.next().split()[:3]]
    legends={}
    siginds={}
    cftrans={}
    for i in range(ncp):
        icp, dim, size = [int(e) for e in lines.next().split()]
        legends[icp] = lines.next().split("'")[1::2]
        siginds[icp] = array([[int(e) for e in lines.next().split()] for row in range(dim)])
        raw_cftrans = array([[float(e) for e in lines.next().split()] for row in range(dim)])
        cftrans[icp] = raw_cftrans[:,0::2] + raw_cftrans[:,1::2]*1j

    return (siginds, cftrans, cps)

if __name__ == '__main__':

    fil='../../alpha-Ce.indmfl'
    (siginds,cftrans) = ReadIndmfl(fil,sys.stdout)
    print(siginds)
    print(cftrans)
