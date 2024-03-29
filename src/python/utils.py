#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import os, sys, glob
from numpy import imag
import subprocess

class DmftEnvironment:
    '''This class provides the following member variables:
    ROOT    - location of EDMFTF executables
    MPI     - command to launch MPI process
    '''
    __shared_state = {}  # class has one common state, but (possibly) many instances
    def __init__(self):
        self.__dict__ = self.__shared_state
        if self.__dict__ == {}:  # only initialize once
            self.__get_root()
            self.__get_mpi()
    def __get_root(self):
        ROOT = os.environ.get('WIEN_DMFT_ROOT')
        if ROOT and os.path.isdir(ROOT):
            self.ROOT = ROOT
        else:
            errstr = "Cannot determine location of EDMFTF executables because "
            errstr += "%s does not exist." % ROOT
            errstr += "Check your WIEN_DMFT_ROOT environment variables."
            raise Exception(errstr)
        if self.ROOT not in sys.path:
            sys.path.append(self.ROOT)
    def __get_mpi(self):
        self.MPI = ''
        mpifile = 'mpi_prefix.dat'
        if os.path.isfile(mpifile):
            self.MPI = open(mpifile, 'r').readline().strip()
            print("DmftEnvironment: mpi_prefix.dat exists -- running in parallel mode.")
            print("  ", self.MPI)
        else:
            print("DmftEnvironment: mpi_prefix.dat does not exist -- running in single-processor mode.")
        mpifile2 = 'mpi_prefix.dat2'
        if os.path.isfile(mpifile2):
            self.MPI2 = open(mpifile2, 'r').readline().strip()
            print("DmftEnvironment: mpi_prefix.dat2 is different from mpi_prefix.dat -- running dmft1/dmft2 in parallel mode.")
            print("  ", self.MPI2)
        else:
            self.MPI2 = self.MPI

class W2kEnvironment:
    '''This class provides the following member variables:
    SCRATCH  - WIEN2k scratch directory
    EDITOR   - editor (vi/emacs) user chose for viewing WIEN2k files
    WIENROOT - location of WIEN2k executables
    case     - casename of current directory
    '''
    __shared_state = {}  # class has one common state, but (possibly) many instances
    def __init__(self):
        self.__dict__ = self.__shared_state
        if self.__dict__ == {}:  # only initialize once
            self.SCRATCH = '.' # os.environ.get('SCRATCH')
            self.EDITOR = os.environ.get('EDITOR')
            self.__get_wienroot()
            self.__get_case()

    def __get_wienroot(self):
        # not used at the moment
        ROOT = os.environ.get('WIENROOT')
        if ROOT:
            if os.path.isdir(ROOT):
                self.WIENROOT = ROOT
            else:
                raise Exception('WIENROOT is set, but %s is not a valid directory.' % ROOT)
        else:
            raise Exception('Cannot determine location of WIEN executables because WIENROOT is not set.')

    def __get_case(self):
        # Determine WIEN2k case name                                                                                                                                                                         
        self.case = os.path.basename(os.getcwd())
        if not (os.path.isfile(self.case+'.struct') and os.path.isfile(self.case+'.in0') ):
            # directory name is not case (happens when submitting to cluster)                                  
            files = glob.glob('*.struct')
            if len(files) < 1:
                raise Exception('No struct file present.')
            elif len(files) > 1:
                # heuristic algorithm to determine case:
                # need in0 and struct present with the same case.
                candidates = [os.path.splitext(f)[0] for f in files]
                for cand in candidates:
                    if os.path.isfile(cand+'.in0'):
                        self.case = cand
                        break
            else: # just one candidate exists, hence it must be it
                self.case, ext = os.path.splitext(os.path.basename(files[0]))
    def __repr__(self):
        string = 'W2kEnvironment: case={:s} WIENROOT={:s} SCRATCH={:s} EDITOR={:s}'.format(self.case, self.WIENROOT,self.SCRATCH, self.EDITOR)
        return string
# energy units conversion between eV and Ry

Ry_in_eV = 13.6056923
eV_in_Ry = 1./Ry_in_eV
def eV2Ry(x):
    return float(x) / Ry_in_eV
def Ry2eV(x):
    return float(x) * Ry_in_eV


def shellcmd(command):
    """ runs a shell command """
    proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = proc.communicate()
    return (out, err, proc.returncode)

# orbital angular momentum conversion between (s,p,d,f,...) and (0,1,2,3,...)

__Lnum2str = ['s', 'p', 'd', 'f', 'g', 'h']

def L2str(Lnum):
    return __Lnum2str[Lnum]

def L2num(Lstr):
    return [L for L,s in enumerate(__Lnum2str) if s==Lstr][0]

def format_matrix(mat, format='%12.8f'):
    if (imag(mat).flat == 0).all():
        # purely real matrix
        ret = '\n'.join([' '.join([format % x for x in row]) for row in mat])
    else:
        ret = '\n'.join(['   '.join([(format+' '+format) % (complex(x).real, complex(x).imag) for x in row]) for row in mat])
    return ret
