#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import sys,re,os
import optparse
from scipy import *
from numpy import *

def saverage_sig(fnames, osig, stdev=False, nexecute=False):
    print('files to average over:', fnames)
    print('output: ', osig)
    ws_oo=[]
    wEdc=[]
    wdata=[]
    for f in fnames:
        data = loadtxt(f)
        wdata.append(data)
        with open(f,'r') as fi:
            for i in range(2):
                line = fi.readline()
                m = re.search('#(.*)',line)
                if m is not None:
                    if not nexecute:
                        exec(m.group(1).strip(), globals())
                
            if 's_oo' in globals(): ws_oo.append(s_oo)
            if 'Edc' in globals(): wEdc.append(Edc)
    
    with open(osig, 'w') as fout:
        if len(ws_oo):
            as_oo = sum(array(ws_oo),axis=0)/len(ws_oo)
            print('s_oo=', as_oo)
            print('# s_oo=', as_oo.tolist(), file=fout)
        
        if len(wEdc):
            aEdc = sum(array(wEdc),axis=0)/len(wEdc)
            print('Edc=', aEdc)
            print('# Edc=', aEdc.tolist(), file=fout)

        wres = sum(array(wdata), axis=0)/len(wdata)
        
        if stdev:
            wstd = sqrt(abs(sum(array(wdata)**2, axis=0)/len(wdata) - wres**2))
            wres = hstack( (wres, wstd[:,1:]) )
            
        savetxt(fout,wres)

if __name__=='__main__':
    """ Takes several self-energy files and produces an average over these self-energy files
    """
    usage = """usage: %prog [ options ] argumens

    The script takes several self-energy files and produces an average self-energy

    arguments  -- all input self-energies
    option -o  -- output self-energy file
    """

    parser = optparse.OptionParser(usage)
    parser.add_option("-o", "--osig", dest="osig", default='sig.inpx', help="filename of the output self-energy file. Default: 'sig.inp'")
    parser.add_option("-l", "--lext", dest="m_extn", default='', help="For magnetic calculation, it can be 'dn'.")
    parser.add_option("-n", "--nexecute", dest="nexecute", action="store_true", default=False,  help="Do not execute the comments")
    parser.add_option("-s", dest="stdev", action='store_true', default=False, help="Computes also the standard deviation - the error of self-energy")

    # Next, parse the arguments
    (options, args) = parser.parse_args()

    saverage_sig(args, options.osig, options.stdev, options.nexecute)
