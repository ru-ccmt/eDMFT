#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import re, os
from numpy import diff
from optparse import OptionParser
# Example line format:
#
# 12000000     nf=7.1504           <Sz>=0.136327         <s>=1                chiS0=4.58693          chiS=3.36033          chiD=0.0308066        Epot=172.736     Ekin=-3.05802   
#
def parse_args():
    parser = OptionParser()
    parser.add_option('-f', '--file', default='nohup_imp.out.000',
                      help='ctqmc output file to parse; default is `nohup_imp.out`')
    opts, args = parser.parse_args()
    return opts

if __name__ == '__main__':
    opts = parse_args()

    with open(opts.file, 'r') as f:
        contents = f.read()

    # Add newlines, since sometimes newline is not appended to each
    # row, probably due to simultaneous writing by parallel processes.
    lines = re.sub('(\d+\s+nf=)', r'\n\1', contents).split('\n')

    # keep only `data' rows
    data = [line for line in lines if 'Ekin' in line]

    # determine how many processes are at each stage
    cpus_at_nsteps = {}
    for line in data:
        nsteps = int(line.split()[0])
        cpus_at_nsteps.setdefault(nsteps,0)
        cpus_at_nsteps[nsteps] += 1

    # cumulative cpus: number of cpus that have reach a given stage
    stages = sorted(cpus_at_nsteps.keys())
    have_data = len(stages)
    if have_data:
        cum_cpus = [cpus_at_nsteps[k] for k in stages]
        cpus = list(-diff(cum_cpus)) + [cum_cpus[-1]]

        # ignore stages which all cpus have completed
        i = 0
        while cpus[i] == 0:
            i += 1

    # display
    print('%9s  %s' % ('qmc step', 'cpus at qmc step'))
    print('-'*40)
    if have_data:
        for s,c in zip(stages[i:],cpus[i:]):
            print('%9d %3d   %s' % (s,c,'*'*c))
        print('-'*40)
        print('%9s %3d   %s' % ('', sum(cpus), 'cpus total'))
    else:
        print(' No data yet.')
