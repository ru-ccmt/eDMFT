#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import re
import optparse
#from scipy import *
from scipy import interpolate
from numpy import loadtxt,savetxt,array,zeros,pi



if __name__=='__main__':
    """ Interpolates existing self-energy on matsubara mesh for different temperature
    """
    usage = """usage: %sinterp.py [ options ]

    Interpolates existing self-energy on matsubara mesh for different temperature
    """

    parser = optparse.OptionParser(usage)
    parser.add_option("-i", "--sinp", dest="insig",  default='sig.inp', help="the mesh will be used from this file.", metavar="FILE")
    parser.add_option("-o", "--sout", dest="outsig", default='sig.inp_', help="the result will be written to this file")
    parser.add_option("-b", "--beta",  dest="beta",  type="float", default=100, help="Inverse temperature")
    # Next, parse the arguments
    (options, args) = parser.parse_args()

    outsig = options.outsig+'b_'+str(options.beta)
    print('insig=%s, outsig=%s, beta=%s' %  (options.insig, outsig, options.beta))



    # Searching for s_oo and Edc
    with open(options.insig, 'r') as fh_sig:
        m = re.search(r'(s_oo\s*=\s*\[.*\])', next(fh_sig))
        if m is not None:
            exec(m.group(1), globals())
        m = re.search(r'(Edc\s*=\s*\[.*\])', next(fh_sig))
        if m is not None:
            exec(m.group(1), globals())
            
    # Read the rest of the input file
    data = loadtxt(options.insig).transpose()
    
    om0 = data[0]
    print('om(old)=', om0)
    # Creates the new Matsubara mesh
    om1 = []
    for i in range(100*len(om0)):
        omx = (2*i+1)*pi/options.beta
        if (omx<om0[-1]):
            om1.append(omx)
        else:
            break
    om1=array(om1)
    print('om(new)=', om1)

    new_data=[]
    new_data.append(om1)
    for i in range(1,len(data)):        
        tck = interpolate.splrep(om0, data[i], k=3, s=0)
        newy = interpolate.splev(om1, tck)
        new_data.append(newy)
    new_data = array(new_data)

    print('Output written to: ', outsig)
    fout = open(outsig, 'w')
    print('# s_oo=', s_oo, file=fout)
    print('# Edc=', Edc, file=fout)
    savetxt(fout, new_data.transpose())

