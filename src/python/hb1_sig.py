#!/usr/bin/env python
from scipy import *
from pylab import *
import optparse
import utils

if __name__ == '__main__':
    usage = """usage: %prog [ options ]

    The script takes the self-energy file (sig.inp), which contains self-energies
    for orbitals treated dynamically, and produces a good guess for Hubbard-I type
    self-energy for orbitals treated in semicore-like approximation.
    """

    parser = optparse.OptionParser(usage)
    parser.add_option("-n", "--occ",  dest="n",    type="float", default=0.5, help="Desired occupancy per orbital")
    parser.add_option("-e", "--ef", dest="ef",  type="float", default=8,  help="Desired position of the semicore state")
    parser.add_option("-i", "--sig", dest="insig",  default='sig.inp', help="filename of the input file (sig.inp)", metavar="FILE")
    parser.add_option("-o", "--sout", dest="outsig", default='sfx.x', help="filename of the output file (sfx.x)")
    parser.add_option("-d", "--delta", dest="delta",  type="float", default=0.02,  help="small broadening of the poles")
    # Next, parse the arguments
    (options, args) = parser.parse_args()

    env = utils.W2kEnvironment()
    case = env.case
    dats = loadtxt(options.insig).transpose()

    w = dats[0]
    matsubara = False
    if w[0]>0: matsubara = True

    
    if matsubara:
        z = w*1j
    else:
        z = w + options.delta*1j
    
    n = options.n
    
    sg = z - 1./( n/(z+options.ef) + (1.-n)/(z-options.ef) )
    s_oo = options.ef*(1-2*n)
    sg -= s_oo
    
    fo = open(options.outsig, 'w')
    print('# s_oo= ['+str(s_oo)+',0]', end=' ', file=fo)
    print('     #   -n '+str(options.n)+' -e '+str(options.ef), file=fo)
    print('# Edc=  [0,0]', file=fo)

    for i in range(len(w)):
        print("%20.16f  %20.16f %20.16f   %4.2f %4.2f" % (w[i], sg[i].real, sg[i].imag, 0, 0), file=fo)
