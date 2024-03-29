#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
from scipy import *
import optparse
import sys, re

def stoint(a):
    ma = re.match('.*_(\d+)$',a)
    xa=0
    if ma is not None: xa = int(ma.group(1))
    return xa

if __name__ == '__main__':
    """ Takes a list of files and concatenates them. It sortes files according to their
    suffix, which has the form <filename>_xx (xx is integer). Such files are used in Wien2k
    (for example case.energy_xxx files)
    Optionally, it cuts a few lines at the beginning of the files (in all except the first file).
    """
    usage = """usage: %prog [ options ]

       Script takes a list of files and concatenates them. It sorts files according to their
    suffix, if it has the form <filename>_xx (xx is integer). Such files are used in Wien2k.
    Optionally, it skips header consisting of a few lines (in all except the first file).
    """

    parser = optparse.OptionParser(usage)
    parser.add_option("-c", "--cut",  dest="cut",    type="int", default=0, help="Number of lines to skip")
    # Next, parse the arguments
    (options, args) = parser.parse_args()

    
    files = sorted(args, key=stoint)
    
    for i,filex in enumerate(files):
        fx = open(filex, 'r')
        for j,line in enumerate(fx):
            if i==0 or (j>=options.cut):
                print(line, end='')
        fx.close()
        

    


