#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
#
# @file dmft_install.py
#
#  DMFT is a software package provided by Rutgers Univiversity,
#  the State University of New Jersey
#
# @version 1.0.0
# @author Kristjan Haule and Viktor Oudovenko
# @date 2016-02-15
#
###

#from .utils import writefile, shellcmd, delfiles, downloader, getURLName
from install_script.iutils import writefile, shellcmd, delfiles, downloader, getURLName
import sys
import os
#from . import framework
import install_script.iframework as framework
import shutil
import re
import subprocess

class Dmft_install(framework.Framework):
    """ This class installs DMFT packages. """
    def __init__(self, argv, config):
        framework.Framework.__init__(self, argv, config)
        self.config  = config

    def resume(self):
        # self.config  = config
        print("\n"+"="*50)
        print("  DMFT installation")
        print("="*50)

        self.write_makeinc()
        if self.src:
            self.version = 'svn'
            self.write_makeinc()
        else:
            self.down_install()

        framework.Framework.resume(self)

    def write_makeinc(self):
        """ Writing Makefile.in files for DMFT installation """
        sdir = os.getcwd()
        #print('sdir=', sdir)
        _where_ = os.path.normpath(sdir+'/../src')
        print('Writing Makefile.in to '+_where_)
        print('You can also check Makefile.in and modify it directly if you want.')
        sys.stdout.flush()

        makeinc = """
#===========================================================================
# Please set environment variable WIEN_DMFT_ROOT (in ~/.bashrc) to the same path as DESTDIR in this Makefile
DESTDIR = """+self.config.prefix+"""
#------------- SERIAL VERSION ----------------------------------------------
F90 = """+self.config.fc+"""
F77 = """+self.config.fc+"""
C++ = """+self.config.cxx+"""
CC  = """+self.config.cc+"""

OPENMP = """+self.config.ompflag+"""
FFLAGS = """+self.config.fflags+""" """+self.config.ldflags+"""
OFLAGS = """+self.config.cflags+""" 

LALIB  =  """+self.config.blaslib+""" """+self.config.lapacklib+"""
FFTLIB =  """+self.config.fftwlib+"""
FFTINC =  """+self.config.fftwinc+"""
GSLLIB =  """+self.config.gsl+"""
GSLINC =  """+self.config.gslinc+"""
PIC    =  """+self.config.PIC+"""
WFOPT  = """+self.config.fflags+""" """+self.config.ldflags+""" $(FFTINC)
WLIBS   = $(FFTLIB)  $(LALIB)
F2PL =

F2PY_LAPACK = """+self.config.f2pylapack+"""

#------------- PARALLEL VERSION --------------------------------------------
Fmpi_define =  """+self.config.mpi_define+"""

PF90 = """+self.config.pfc+""" 
PC++ = """+self.config.pcxx+"""
PCC  = """+self.config.pcc+""" 

PFLAGS   = """+self.config.mpi_define+""" -DMPICH_IGNORE_CXX_SEEK -O3 #-restrict -ipo -no-prec-div 

LLIBS = $(LALIB)
PLIBS = $(LLIBS) $(GSLLIB)

#
CMP = f2py """+self.config.f2pyflag+"""  # fortran to python converter
CMPLIBS  = """+self.config.f2pylib+""" """+self.config.f2pylapack+"""
CMPLIBS2 = --f90flags=' $(OPENMP) ' $(CMPLIBS)
#
PYBND = """+self.config.pybind11+"""
#============================================================================
"""

        writefile('../src/Makefile.in',makeinc)
        #comm = 'mv ../Makefile.in Makefile.in.`date +%s`; cp -r make.inc ../Makefile.in '
        #sys.stdout.flush()
        #(output, error, retz) = shellcmd(comm)
        #if(retz != 0):
        #    print('\n\nDMFT: Makefile.in copy failed. Aborting...')
        #    print('error is:\n','*'*50,'\n',comm,'\n',error,'\n','*'*50)
        #    sys.exit()


    def down_install(self):
        """ Download and install DMFT """
        savecwd = os.getcwd()
        #logdir = os.path.join(savecwd,'log')
        #if not os.path.isdir(logdir): os.mkdir(logdir)
        versions = self.versions
        rep_name = os.path.normpath(savecwd+'/../src')
        os.chdir(rep_name)
        
        print('You are ready to compile the DMFT code..... ')
        print('..... If you want to compile the code manually, enter ./src  and type "make". You can also type "make clean" to clean.')
        print('Please press enter to start compilation')
        if self.non_interactive:
            print('Mode is non_interactive=', self.non_interactive)
        else:
            input()
        
        comm = self.make
        #(output, error, retz) = shellcmd(comm)
        #print('current dir=', rep_name, 'comm=', comm)
        retz=subprocess.call(comm,shell=True,stdout=sys.stdout,stderr=sys.stderr)
        
        if retz:
            print('\n\nDMFT: error building DMFT')
            #print 'stderr:\n','*'*50,'\n',error,'\n','*'*50
            #writefile(os.path.join(savecwd,'log/dmft.log'), output+error)
            sys.exit()

        #liblog = os.path.join(savecwd,'log/dmft.log')
        #writefile(liblog, output+error)
        print('*'*50+'\n')
        print('Compilation of DMFT was successful.')
        print('*'*50+'\n')
        #print '(Log is in ',liblog,')'

        # Installation
        print('\nPopulatig instalation directory with executables', end=' ')
        sys.stdout.flush()

        # ADD some installation stuff here
        
        comm = 'make install'
        #(output, error, retz) = shellcmd(comm, sys.stdout)
        retz=subprocess.call(comm,shell=True,stdout=sys.stdout,stderr=sys.stderr)
        if retz:
           print('\nDMFT: make install failed... ')
           print('... Consider changing configure.py, or change some input options.')
           print('... Check "src/Makefile.in" if it is compatible with your system')
           print('... You might type make in directory src, to retry the compilation')
           #print 'stderr:\n','*'*50,'\n',error,'\n','*'*50
           #writefile(liblog, output+error)
           sys.exit()



