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
    #parser.add_option("-n", "--occ",  dest="n", type="float", default=0.5, help="Desired occupancy per orbital")
    parser.add_option("-n", "--occ",  dest="n", default="[0.5]", help="Desired occupancy per orbital, default [0.5]")
    parser.add_option("-e", "--ef", dest="ef",  default="[8.0]", help="Desired position of the semicore state, default [8.0]")
    parser.add_option("-i", "--sig", dest="insig",  default='sig.inp', help="filename of the input file (sig.inp)", metavar="FILE")
    parser.add_option("-o", "--sout", dest="outsig", default='sfx.x', help="filename of the output file (sfx.x)")
    parser.add_option("-d", "--delta", dest="delta",  type="float", default=0.02,  help="small broadening of the poles")
    # Next, parse the arguments
    (options, args) = parser.parse_args()
    
    nf = eval(options.n)
    ef = eval(options.ef)
    print('n=', nf, 'ef=', ef)
    if len(nf)!=len(ef):
        print('It is required that len(n)==len(e)')
        sys.exit(0)
    
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
    
    _s_oo_=[]
    _sg_=[]
    for ia in range(len(nf)):
        n = nf[ia]
        ene = ef[ia]
        sg = z - 1./( n/(z+ene) + (1.-n)/(z-ene) )
        s_oo = ene*(1-2*n)
        sg -= s_oo
        _s_oo_.append(s_oo)
        _sg_.append(sg)
    
    _Edc_ = [0 for s in _s_oo_]
    with open(options.outsig, 'w') as fo:
        print('# s_oo= '+str(_s_oo_), end=' ', file=fo)
        print('     #   -n '+str(options.n)+' -e '+str(options.ef), file=fo)
        print('# Edc= '+str(_Edc_), file=fo)
        
        for i in range(len(w)):
            print("{:20.16f}".format(w[i]), end="  ", file=fo)
            for j in range(len(nf)):
                print("{:20.16f} {:20.16f}".format(_sg_[j][i].real, _sg_[j][i].imag), end=' ', file=fo)
            print(file=fo)
        
