#! /usr/bin/env python
# -*- coding: utf-8 -*-

###
#
# @file lapack.py
#
#  DMFT is a software package provided by Rutgers Univiversity,
#  the State University of New Jersey
#
# @version 1.0.0
# @author Kristjan Haule and Viktor Oudovenko
# @date 2016-02-15
#
###


#from .utils import writefile, shellcmd, delfiles, downloader, geturl, getURLName
from install_script.iutils import writefile, shellcmd, delfiles, downloader, geturl, getURLName
import sys
import os
import urllib.request, urllib.parse, urllib.error
import shutil
#from . import framework
import install_script.iframework as framework

class Lapack(framework.Framework):
    """ This class takes care of the LAPACK installation. """
    def __init__(self, config, dmft):
        print("\n"+"="*50)
        print("  LAPACK installation/verification")
        print("="*50)

        self.config     = config
        self.downlapack = dmft.downlapack
        self.downcmd    = dmft.downcmd
#        self.prefix     = dmft.prefix
        self.dmft     = dmft
        self.lapversion = dmft.lapversion
        self.lapackurl  = "http://www.netlib.org/lapack/"+self.lapversion+".tgz"

        if self.downlapack == 2:
            self.down_install_lapack()

        ret = self.check_lapack()

        if ret != 0:
            if self.downlapack == 1:
                self.down_install_lapack()
            else:
                if not os.path.isfile(os.path.join(self.config.prefix,'lib/liblapack.a')):
                   print("""
Please provide a working LAPACK library. If a LAPACK library is not
present on the system, the netlib LAPACK library can be automatically
downloaded and installed by adding the --downlapack flag.
Most used BLAS implementations already include the LAPACK library as
MKL, ACML, Goto, Goto2 or ATLAS. If you want to use one of these
libraries, you just have to specify correctly the --blaslib option or
you can specify where is located your own LAPACK library by using the
--lapacklib option.

What do you want to do ?
    - d : download the netlib LAPACK
    - q : quit and try with another LAPACK library or define lapacklib parameter.
                   """)
                   answer = input(">[q] ")
                   if answer == "d":
                       self.down_install_lapack()
                   else:
                       sys.exit()
                else:
                  print("Netlib Lapack library is already installed at "+os.path.join(self.config.prefix,'lib/liblapack.a'))
                  print("Do you want to try it? (t)  or proceed without testing (p) or quit (q) ?")
                  answer = input(">[q] ")
                  if answer == "t":
                    self.config.lapacklib  = '-L'+os.path.join(self.config.prefix,'lib')+' -ltmg -llapack'
                    self.check_lapack()
                  elif answer == "p":
                    exit 
                  else:
                    sys.exit()


    def check_lapack(self):
        print("Checking if provided LAPACK works...", end=' ')
        # This function simply generates a C program
        # that contains few calls to LAPACK routine and then
        # checks if compilation, linking and execution are succesful

        sys.stdout.flush()
        writefile('tmpf.f',"""
      program ftest
      integer  N
      parameter (N = 1)
      double precision A(N, N), B(N)
      integer  I(N)
      integer  INFO
      B(:)   = 1
      A(:,:) = 2
      I(:)   = 0
      call cheevd( 'N', 'U', N, A, N, B, B, -1,
     $     B, -1, I, -1, INFO)
      stop
      end\n""")

        ldflg = self.config.lapacklib+' '+self.config.blaslib+'  -lm'
        ccomm = self.config.fc+' -o tmpf '+'tmpf.f '+ldflg
        (output, error, retz) = shellcmd(ccomm)

        if(retz != 0):
            if self.dmft.verbose:
                print('\n\nLAPACK: provided LAPACK cannot be used! aborting...')
                print('error is:\n','*'*50,'\n',ccomm,'\n',error.decode('UTF-8'),'\n','*'*50)
            else:
                print("no")
            return -1;

        comm = './tmpf'
        (output, error, retz) = shellcmd(comm)

        if(retz != 0):
            if self.dmft.verbose:
                print('\n\nLAPACK: provided LAPACK cannot be used! aborting...')
                print('error is:\n','*'*50,'\n',comm,'\n',error.decode('UTF-8'),'\n','*'*50)
            else:
                print("no")
            return -1;

        delfiles(['tmpf.f','tmpf'])
        print("yes")
        #print("\nChecking if provided LAPACK contains functions for test works...",)
        ## This function simply generates a C program
        ## that contains few calls to LAPACK routine and then
        ## checks if compilation, linking and execution are succesful
        #
        #sys.stdout.flush()
        #writefile('tmpf.f',"""
        #program ftest
        #double precision D(1), A(1:1), B(2)
        #integer          ISEED( 4 )
        #integer          INFO
        #B(1)   = 1
        #
        #do  I = 1, 4
        #    ISEED( I ) = 1
        #enddo
        #call dlarnv( 1, ISEED, 1, D )
        #call dlagsy( 1, 0, D, A, 1, ISEED, B, INFO )
        #stop
        #end\n""")
        #
        #ccomm = self.config.fc+' -o tmpf '+'tmpf.f '+ldflg
        #(output, error, retz) = shellcmd(ccomm)
        #
        #print('HERE output=', output, 'error=', error, 'retz=', retz)
        # 
        #if(retz != 0):
        #    print('no')
        #    self.dmft.needtmg = 1
        #else:
        #    comm = './tmpf'
        #    (output, error, retz) = shellcmd(comm)
        #    if(retz != 0):
        #        print('no')
        #        self.dmft.needtmg = 1;
        #    else:
        #        self.dmft.needtmg = 0;
        #        #   print 'yes'
        #delfiles(['tmpf.f','tmpf'])

        return 0;


    def down_install_lapack(self):

        print("""
The LAPACK library is being installed.
""")
        sys.stdout.flush()

        savecwd = os.getcwd()

        # creating the build,lib and log dirs if don't exist
        if not os.path.isdir(os.path.join(self.config.prefix,'lib')):
            os.mkdir(os.path.join(self.config.prefix,'lib'))

        if not os.path.isdir(os.path.join(os.getcwd(),'log')):
            os.mkdir(os.path.join(os.getcwd(),'log'))

        # Check if lapack.tgz is already present in the working dir
        # otherwise download it
        if not os.path.isfile(os.path.join(self.config.prefix,'lib/liblapack.a')):
#        if not os.path.isfile(os.path.join(os.getcwd(),getURLName(self.lapackurl))):


            print("Downloading LAPACK...", end=' ')
#            downloader(self.lapackurl,self.downcmd)
#            urllib.urlretrieve(self.lapackurl, "lapack.tgz")
            geturl(self.lapackurl, "lapack.tgz")
            print("Download is done")
        else:
            print("Netlib Lapack library is already installed at "+os.path.join(self.config.prefix,'lib/liblapack.a'))
            self.config.lapacklib  = '-L'+os.path.join(self.config.prefix,'lib')+' -ltmg -llapack'
            return 0;

        # unzip and untar
        os.chdir('download') 
        print('Unzip and untar Lapack...', end=' ')
        #comm = 'gunzip -f '+self.lapversion+'.tgz'
        comm = 'tar zxf  lapack.tgz'
        (output, error, retz) = shellcmd(comm)
        if retz:
            print('\n\nLAPACK: cannot unzip '+self.lapversion+'.tgz')
            print('stderr:\n','*'*50,'\n',comm,'\n',error.decode('UTF-8'),'\n','*'*50)
            sys.exit()

        #comm = 'tar xf '+self.lapversion+'.tar'
        #(output, error, retz) = shellcmd(comm)
        #if retz:
        #    print '\n\nLAPACK: cannot untar '+self.lapversion+'.tar'
        #    print 'stderr:\n','*'*50,'\n',comm,'\n',error,'\n','*'*50
        #    sys.exit()
        #os.remove(self.lapversion+'.tar')
        print('done')

        ##Apply the patch to correct [sd]lantr
        #print 'Apply patch on lapacke...',
        #comm = '(cd '+self.lapversion+' && patch -p 0 < '+(os.path.join(savecwd,'../script/patch_lantr'))+')'
        #(output, error, retz) = shellcmd(comm)
        #print 'done'

#         # Overwrite [sd]lamch.f
#         shutil.copy(os.path.join(self.dmft.installerdir,'src/dlamch.f'),
#                     os.path.join(os.getcwd(),'lapack-3.3.1/INSTALL'))
#         shutil.copy(os.path.join(self.dmft.installerdir,'src/slamch.f'),
#                     os.path.join(os.getcwd(),'lapack-3.3.1/INSTALL'))

        # change to BLAS dir
        os.chdir(os.path.join(os.getcwd(), self.lapversion))

        # Write Makefile.in
        writefile('make.inc', """
# -*- Makefile generated by DMFT installer -*-
####################################################################
#  LAPACK make include file.                                       #
#  LAPACK, Version """+self.lapversion+""""                                           #
#  April 2012                                                      #
####################################################################
#
SHELL = /bin/sh
#
#  Modify the FORTRAN and OPTS definitions to refer to the
#  compiler and desired compiler options for your machine.  NOOPT
#  refers to the compiler options desired when NO OPTIMIZATION is
#  selected.  Define LOADER and LOADOPTS to refer to the loader and
#  desired load options for your machine.
#
FORTRAN  = """+self.config.fc+"""
OPTS     = """+self.config.fflags+"""
DRVOPTS  = $(OPTS)
NOOPT    = -O0
LOADER   = """+self.config.fc+"""
LOADOPTS = 
MAKE     = make -j 8
#
# Timer for the SECOND and DSECND routines
#
# Default : SECOND and DSECND will use a call to the EXTERNAL FUNCTION ETIME
# TIMER    = EXT_ETIME
# For RS6K : SECOND and DSECND will use a call to the EXTERNAL FUNCTION ETIME_
# TIMER    = EXT_ETIME_
# For gfortran compiler: SECOND and DSECND will use a call to the INTERNAL FUNCTION ETIME
# TIMER    = INT_ETIME
# If your Fortran compiler does not provide etime (like Nag Fortran Compiler, etc...)
# SECOND and DSECND will use a call to the INTERNAL FUNCTION CPU_TIME
TIMER    = INT_CPU_TIME
# If neither of this works...you can use the NONE value... In that case, SECOND and DSECND will always return 0
# TIMER     = NONE
#
#  Configuration LAPACKE: Native C interface to LAPACK
#  To generate LAPACKE library: type 'make lapackelib'
#  Configuration file: turned off (default)
#  Complex types: C99 (default)
#  Name pattern: mixed case (default)
#  (64-bit) Data model: LP64 (default)
#
# CC is the C compiler, normally invoked with options CFLAGS.
#
CC     = """+self.config.cc+"""
CFLAGS = """+self.config.cflags+"""
#
#  The archiver and the flag(s) to use when building archive (library)
#  If you system has no ranlib, set RANLIB = echo.
#
ARCH     = ar
ARCHFLAGS= """+self.config.arflags+"""
RANLIB   = """+self.config.ranlib+"""
#
#  The location of BLAS library for linking the testing programs.
#  The target's machine-specific, optimized BLAS library should be
#  used whenever possible.
#
BLASLIB      = """+self.config.blaslib+"""
#
#  Location of the extended-precision BLAS (XBLAS) Fortran library
#  used for building and testing extended-precision routines.  The
#  relevant routines will be compiled and XBLAS will be linked only if
#  USEXBLAS is defined.
#
# USEXBLAS    = Yes
XBLASLIB     =
# XBLASLIB    = -lxblas
#
#  Names of generated libraries.
#
LAPACKLIB    = liblapack.a
TMGLIB       = libtmg.a
EIGSRCLIB    = libeigsrc.a
LINSRCLIB    = liblinsrc.a
LAPACKELIB   = liblapacke.a
""")

        # compile and generate library
        print('Compile and generate LAPACK...', end=' ')
        sys.stdout.flush()
        comm = self.make+' lapacklib tmglib'
        (output, error, retz) = shellcmd(comm)
        if retz:
            print("\n\nLAPACK: cannot compile LAPACK")
            print("stderr:\n","*"*50,"\n",comm,'\n',error.decode('UTF-8'),"\n","*"*50)
            sys.exit()


        # write the log on a file
        log = output.decode('UTF-8')+error.decode('UTF-8')
        fulllog = os.path.join(savecwd,'log/log.lapack')
        writefile(fulllog, log)
        print('Installation of liblapack.a successful.')
        print('(log is in ',fulllog,')')

        # move libcblas.a to the lib directory
        shutil.copy('liblapack.a',os.path.join(self.config.prefix,'lib/liblapack.a'))
        shutil.copy('libtmg.a',os.path.join(self.config.prefix,'lib/libtmg.a'))

        # set framework variables to point to the freshly installed BLAS library
        self.config.lapacklib  = '-L'+os.path.join(self.config.prefix,'lib')+' -ltmg -llapack'
        os.chdir(savecwd)

        self.config.lapinstalled = 1;

        # Check if the installation is successful
        self.dmft.verbose = 1
        ret = self.check_lapack()
        self.dmft.verbose = 0
        if ret != 0:
            sys.exit()
