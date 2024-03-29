#!/usr/bin/env python
from scipy import *
import utils
import sys, os
from wstruct import Struct
import re
from w2k_atpar import readpotential, readlinearizatione, atpar, rint13, rint13g
from pylab import *
import optparse

def read_core(case, struct,updn=''):
    # reading case.inc
    fc = open(case+'.inc'+updn, 'r')
    n_kappa_occup =[]
    iprint=[]
    for jatom in range(struct.nat):
        line = next(fc)
        data = line.split()
        nc = int(data[0])
        iprint.append( int(data[2]) )      # saving iprint
        #print('line=', line, 'nc=', nc)
        #nc = int(line[:3])
        nko = zeros((nc,3),dtype=int)
        for il in range(nc):
            line = next(fc)
            n, kappa, occup = (int(line[:1]), int(line[2:4]), int(line[5:6]))
            nko[il,:] = array([n,kappa,occup])
        n_kappa_occup.append(nko)
        
    print('iprint=', iprint)
    
    # reading case.corewf
    fi = open(case+'.corewf'+updn,'r')
    wavec=[]
    cType=[]
    cEne=[]
    for jatom in range(struct.nat):
        line = next(fi).split()
        nc0 = int(line[0])
        nc = len(n_kappa_occup[jatom])
        noccup0 = n_kappa_occup[jatom][0,2]
        #print noccup0, nc, nc0
        waveci=[]
        cTypei=[]
        cEnei=[]
        if (noccup0 and iprint[jatom]):
            for ic in range(nc):
                line = next(fi)
                m1 = re.search('CORE STATES = (\d\S)', line)
                m2 = re.search('CORE ENERGY=(.*) Ry', line)
                if m1 is not None and m2 is not None:
                    ctype = m1.group(1)
                    cenergy = float(m2.group(1))

                ucore = zeros((2,struct.jrj[jatom]))
                for i in range(2):
                    n=0
                    while (n<struct.jrj[jatom]):
                        line = next(fi)
                        dat = [line[3:3+19], line[3+19:3+2*19], line[3+2*19:3+3*19], line[3+3*19:3+4*19]]
                        dat = [s for s in dat if len(s)>1] # throw away some empty strings
                        dat = list(map(float,dat))
                        ucore[i,n:n+len(dat)]= dat[:]
                        n += len(dat)
                waveci.append(ucore)
                print(ctype, cenergy, shape(ucore))
                cTypei.append( ctype )
                cEnei.append(cenergy)
                
        wavec.append(waveci)
        cType.append(cTypei)
        cEne.append(cEnei)

    return (wavec, cType, cEne, iprint, n_kappa_occup)

if __name__ == '__main__':

    usage = """usage: %prog [ options ] mode
       Prints the core wave function from case.corewf[up|dn| ]
       --up  updn = up
       --dn  updn = dn
    """
    parser = optparse.OptionParser(usage)
    parser.add_option("--up", dest="up", action='store_true', default=False, help="magnetic LDA calculation with vector-up first")
    parser.add_option("--dn", dest="dn", action='store_true', default=False, help="magnetic LDA calculation with vector-dn first")
    # Next, parse the arguments
    (options, args) = parser.parse_args()
    updn=''
    if options.up:
        updn = 'up'
    if options.dn:
        updn = 'dn'

    w2k = utils.W2kEnvironment()
    #struct = struct1.Struct(w2k.case)
    struct = Struct(w2k.case)
    struct.ReadStruct(w2k.case+'.struct')
    
    (wavec, cType, cEne, iprint, n_kappa_occup) = read_core(w2k.case, struct,updn)

    for jatom in range(struct.nat):
        if iprint[jatom]:
            for ic in range(len(cType[jatom])):
                #jatom,ic = 0,4
                print(cType[jatom][ic])
                print(shape(wavec[jatom][ic]))
                
                Nr = struct.jrj[jatom]  # Number of radial points from struct file
                r0 = struct.r0[jatom]
                Rmt = struct.rmt[jatom]
                dx = log(Rmt/r0)/(Nr-1)
                Rx = r0 * exp( arange(Nr)*dx ) # radial points
                Ag = wavec[jatom][ic][0]
                Bg = wavec[jatom][ic][1]
                norm = rint13g(1.0,1.0,Ag,Bg,Ag,Bg,dx,r0)
                
                plot(Rx, Ag, label=cType[jatom][ic])
                savetxt('corewf'+updn+'_'+str(jatom)+'_'+str(ic)+'.dat', transpose(vstack((Rx,Ag,Bg))) )
                print(cType[jatom][ic], 'norm=', norm)
    title('Core wave fnunctions')
    legend(loc='best')
    show()
