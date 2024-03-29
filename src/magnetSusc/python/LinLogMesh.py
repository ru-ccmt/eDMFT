#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule

from scipy import *
from scipy import interpolate
import sys
import optparse

def LinLogMeshGen(delta=0.0025,ommax=20.32,Nd=64):
    """Creates logarithmic mesh of linear meshes.
       In the first level we have Nd*2 points in linear mesh with spacing delta.
       In the second level we add Nd*2 points with spacing 2*delta.
       In each subsequent level we add Nd*2 points with doubled spacing.
    """
    #(1+2+2**2+2**3+...)=2**N-1=ommax/(delta*Nd)
    # N=log2(ommax/(delta*Nd)+1)
    mesh=[]
    Nitt = int(round(log2(1+ommax/(delta*Nd))))
    x=0.0
    for itt in range(Nitt):
        mesh += linspace(x+delta, x+delta*Nd, Nd).tolist()
        x += delta*Nd
        delta*=2
    mesh = array(mesh)
    zero_ind = len(mesh)
    mesh = hstack( ( -mesh[::-1], [0], mesh) )
    
    """Create an index array such that mesh[idx[i][j]] constitutes a linear mesh with regular spacings.
    In the first level  mesh[idx[0][:]] will give linear mesh of points in the interval [-delta*Nd,delta*Nd]
    In the second level mesh[idx[1][:]] will give linear mesh of points in the interval [-(1+2)*delta*Nd, (1+2)*delta*Nd] with spacing 2*delta
    In the k-th level   mesh[idx[k][:]] will give linear mesh of points in the interval [-(2**k-1)*delta*Nd, (2**k-1)*delta*Nd] with spacing 2**(k-1)*delta
    """
    idx=[]
    subidx=[zero_ind]
    for itt in range(Nitt):
        subidx=range(zero_ind-(itt+1)*Nd,zero_ind-itt*Nd)+subidx[::2]+range(zero_ind+itt*Nd+1,zero_ind+(itt+1)*Nd+1)
        idx.append(subidx)
    
    return (mesh, idx)

if __name__ == '__main__':
    
    usage = """usage: %prog [ options ]

    The script interpolates self-energy on log-linear mesh as appropriate for computing susceptibility.
    """

    parser = optparse.OptionParser(usage)
    parser.add_option("-i", "--sinp", dest="insig",  default='sig.inp', help="filename of the input file (sig.inp)", metavar="FILE")
    parser.add_option("-d", "--delta", dest="delta", type="float", default=0.0025,  help="The smallest frequency spacing in linear-log mesh for chi")
    parser.add_option("-m", "--ommax", dest="ommax", type="float", default=-1000,  help="The largest frequency for computing chi. Default==max(sig.inp)")
    parser.add_option("-n", "--Nd",  dest="Nd",      type="int", default=64, help="Number of points in linear mesh. Should be 2**k")
    parser.add_option("-o", "--sout", dest="outsig", default='sig.inpr', help="filename of the output file (sig.inp)")
    parser.add_option("-s", "--skipheader", dest="skipheader", type="int", default=0,  help="Should we skip reading of the header (s_oo and Edc)")
    # Next, parse the arguments
    (options, args) = parser.parse_args()

    line1=''
    line2=''
    if not options.skipheader:
        # Searching for s_oo and Edc
        fh_sig = open(options.insig, 'r')
        line1 = fh_sig.next()
        line2 = fh_sig.next()
        fh_sig.close()
    
    # Read the rest of the input file
    sigdata = loadtxt(options.insig).transpose()  # self-energy from 'sig.inp' on long mesh

    om = sigdata[0]
    ommax = 0.5*(-om[0] + om[-2])
    if options.ommax>0: ommax = options.ommax
    #print 'ommax=', ommax
    
    (oml,idxl) = LinLogMeshGen(options.delta,ommax,options.Nd)

    fmesh = open('rmesh.dat','w')
    print >> fmesh, options.delta, oml[-1], options.Nd, '# delta, ommax, Nd'
    for l in range(len(oml)):
        print >> fmesh, l, oml[l]
    fmesh.close()

    
    #omnew = hstack( (oml, [om[-1]]) )

    sigout=[oml]
    for b in range(1,len(sigdata)):
        Sg = interpolate.UnivariateSpline(om,sigdata[b],s=0)
        sigout.append( Sg(oml) )
    sigout = transpose(array(sigout))

    
    fo=open(options.outsig,'w')
    print >> fo, line1, line2,
    for i in range(len(sigout)):
        for j in range(len(sigout[i])):
            print >> fo, sigout[i][j],
        print >> fo
    fo.close()
        
    
    sys.exit(0)
    
    delta=0.0025
    ommax=20.32
    Nd=64
    
    (oml,idxl) = LinLogMeshGen(delta,ommax,Nd)
    zero_ind = oml.tolist().index(0.0)
    Oml=oml[zero_ind:]
    NOm=int(len(Oml))

    for iOm in range(1,NOm):
        level = (iOm-1)/Nd
        om_idx=idxl[level]                      # index for the linear mesh on this level
        izero= (len(om_idx)-1)/2                # zero_ind on this level
        dOm = om_idx.index(zero_ind+iOm)-izero  # om-Om in integer notation is i-dOm
        dx=float(Oml[iOm]-Oml[iOm-1])           # dx for integration
        om_idx=array(om_idx)                    #
        print 'to_integrate:', iOm, level, Oml[iOm], dOm, dOm+izero, dx, izero, oml[om_idx[izero]]
        for iom in range(izero,izero+dOm+1):
            
            i1 = om_idx[iom]
            i2 = om_idx[iom-dOm]
            print iom, i1, i2, oml[i1], oml[i2]
    


