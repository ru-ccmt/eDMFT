#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
#
# @file blas.py
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
import re

class Blas:
    blasurl      = "http://netlib.org/blas/blas.tgz"

    """ This class takes care of the BLAS installation. """
    def __init__(self, config, dmft):
        print("\n"+"="*50)
        print("  BLAS installation/verification")
        print("="*50)

        self.config  = config
#        self.downcmd = dmft.downcmd
#        self.prefix  = dmft.prefix
        downblas     = dmft.downblas;
        
#        print "BLAS LIB is ", config.blaslib
        if config.blaslib:
            self.check_blas()
        else:
            if downblas:
                self.down_install_blas()
            else:
                if not os.path.isfile(os.path.join(self.config.prefix,'lib/librefblas.a')):
                   print("""
Please provide a working BLAS library. If a BLAS library
is not present on the system, the netlib  BLAS library it can be
automatically downloaded and installed by adding the --downblas flag.
Be aware that a netlib  BLAS library will be installed with the --downblas
flag so don't expect good performance.
You can specify the optimized BLAS library **installed** on your machine using the --blaslib option
For instance
    --blaslib="-lessl" for IBM BLAS,
    --blaslib="-framework veclib" for Mac OS/X,
    --blaslib="-lgoto" for GOTO BLAS
    --blaslib="-lf77blas -lcblas -latlas" for ATLAS
    --blaslib="-lmkl_em64t -lguide" for MKL on emt64 architecture (remember to set environment variable MKL_NUM_THREADS=1)
    --blaslib="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core" for single threaded MKL on 64-bit architectures using Intel compilers
    --blaslib="-lmkl_intel -lmkl_sequential -lmkl_core" for single threaded MKL on 32-bit architectures using Intel compilers
    --blaslib="-lmkl_gf -lmkl_sequential -lmkl_core" for single threaded MKL on 32-bit architectures using GNU Fortran compilers
    etc...  .'

What do you want to do ?
    - d : download the netlib BLAS
    - q : quit and try with another BLAS library or define blaslib parameter.
                """)
                   answer = input(">[q] ")
                   if answer == "d":
                       self.down_install_blas()
                   else:
                       sys.exit()
                else:
                  print("Netlib BLAS library is already installed at "+os.path.join(self.config.prefix,'lib/liblapack.a'))
                  print("Do you want to try it? (t)  or proceed without testing (p) or quit (q) ?")
                  answer = input(">[q] ")
                  if answer == "t":
                    self.config.blaslib  = '-L'+os.path.join(self.config.prefix,'lib')+' -lrefblas'
                    self.check_blas()
                  elif answer == "p":
                    exit 
                  else:
                    sys.exit()


    def check_blas(self):

        print("Checking if provided BLAS works...", end=' ')
        # This function simply generates a FORTRAN program
        # that contains few calls to BLAS routine and then
        # checks if compilation, linking and execution are succesful

        # Try to detect which BLAS is used
        if re.search('mkl', self.config.blaslib, re.IGNORECASE):
            self.config.blasname = "mkl"
        elif re.search('acml', self.config.blaslib, re.IGNORECASE):
            self.config.blasname = "acml"
        if self.config.blasname != "Unknown":
            if self.config.compiler == "Intel":
                self.config.cflags    += '' #" -openmp"
        #        self.config.ldflags_c  += " -openmp"
        #        self.config.ldflags_fc += " -openmp"
            elif self.config.compiler == "GNU":
                self.config.cflags    += '' #" -fopenmp"
        #        self.config.ldflags_c  += " -fopenmp"
        #        self.config.ldflags_fc += " -fopenmp"
                
        sys.stdout.flush()
        writefile('tmpf.f',"""
      program ftest
      double precision da, dx(1)
      dx(1)=1
      da = 2
      call dscal(1,da,dx,1)
      stop
      end\n""")

        #fcomm = self.config.fc+' -o tmpf '+'tmpf.f '+self.config.blaslib+' '+self.config.ldflags_fc+' -lm'
        fcomm = self.config.fc+' -o tmpf '+'tmpf.f '+self.config.blaslib+' -lm'
        (output, error, retz) = shellcmd(fcomm)

        if(retz != 0):
            print('\n\nBLAS: provided BLAS cannot be used! aborting...')
            print('error is:\n','*'*50,'\n',fcomm,'\n',error.decode('UTF-8'),'\n','*'*50)
            sys.exit()

        comm = './tmpf'
        (output, error, retz) = shellcmd(comm)
        if(retz != 0):
            print('\n\nBLAS: provided BLAS cannot be used! aborting...')
            print('error is:\n','*'*50,'\n',comm,'\n',error.decode('UTF-8'),'\n','*'*50)
            sys.exit()

        delfiles(['tmpf.f','tmpf'])
        print('yes')

        return 0;


    def down_install_blas(self):
        print("""
The netlib  BLAS library is being installed.
Don't expect high performance from this netlib  library!
If you want performance, you need to use an optimized BLAS library and,
to avoid unnecessary complications, if you need to compile this optimized BLAS
library, use the same compiler you're using here.""")
        sys.stdout.flush()

        savecwd = os.getcwd()

        # creating the build,lib and log dirs if don't exist
        if not os.path.isdir(os.path.join(self.config.prefix,'lib')):
            os.mkdir(os.path.join(self.config.prefix,'lib'))

        if not os.path.isdir(os.path.join(os.getcwd(),'log')):
            os.mkdir(os.path.join(os.getcwd(),'log'))

        # Check if blas.tgz is already present in the working dir
        # otherwise download it
        if not os.path.isfile(os.path.join(self.config.prefix,'lib/librefblas.a')):
#        if not os.path.isfile(os.path.join(os.getcwd(),getURLName(self.blasurl))):
#        if not os.path.isfile(os.path.join(os.getcwd(),'BLAS')):
            print("Downloading BLAS...", end=' ')
#            downloader(self.blasurl,self.downcmd)
#            urllib.urlretrieve(self.blasurl, "blas.tgz")
            geturl(self.blasurl, "blas.tgz")
            print("Download is done")
        else:
            print("Netlib Blas library is already installed at "+os.path.join(self.config.prefix,'lib/librefblas.a'))
            self.config.blaslib  = '-L'+os.path.join(self.config.prefix,'lib')+' -lrefblas '
            return 0;

 
        # unzip and untar
        os.chdir('download')
        print('Unzip and untar netlib  BLAS...', end=' ')
#        comm = 'gunzip -f blas.tgz'
        comm = 'mkdir BLAS; tar zx --strip-components=1 -C BLAS -f blas.tgz '
        (output, error, retz) = shellcmd(comm)
        if retz:
            print('\n\nBLAS: cannot unzip blas.tgz')
            print('stderr:\n','*'*50,'\n',comm,'\n',error.decode('UTF-8'),'\n','*'*50)
            sys.exit()

#        comm = 'mkdir BLAS && tar x --strip-components=1 -C BLAS -f blas.tar'
#        (output, error, retz) = shellcmd(comm)
#        if retz:
#            print '\n\nBLAS: cannot untar blas.tgz'
#            print 'stderr:\n','*'*50,'\n',comm,'\n',error,'\n','*'*50
#            sys.exit()
#        os.remove('blas.tar')
        print('done')

        # change to BLAS dir
        os.chdir(os.path.join(os.getcwd(),'BLAS'))

        # compile and generate library
        print('Compile and generate netlib  BLAS...', end=' ')
        sys.stdout.flush()
        comm = self.config.fc+' '+self.config.fflags+" -c *.f"
        (output, error, retz) = shellcmd(comm)
        if retz:
            print("\n\nBLAS: cannot compile blas")
            print("stderr:\n","*"*50,"\n",comm,'\n',error.decode('UTF-8'),"\n","*"*50)
            sys.exit()

        log = output.decode('UTF-8')+error.decode('UTF-8')

        comm = "ar cr librefblas.a *.o"
        (output, error, retz) = shellcmd(comm)
        if retz:
            print("\n\nBLAS: cannot create blas library")
            print("stderr:\n","*"*50,"\n",comm,'\n',error.decode('UTF-8'),"\n","*"*50)
            sys.exit()
        print("done")

        log += output.decode('UTF-8')+error.decode('UTF-8')

        comm = self.config.ranlib+" librefblas.a"
        (output, error, retz) = shellcmd(comm)
        if retz:
            print("\n\nBLAS: cannot create table of contents for blas library")
            print("stderr:\n","*"*50,"\n",comm,'\n',error.decode('UTF-8'),"\n","*"*50)
            sys.exit()
        print("done")

        # write the log on a file
        log += output.decode('UTF-8')+error.decode('UTF-8')
        fulllog = os.path.join(savecwd,'log/log.blas')
        writefile(fulllog, log)
        print('Installation of netlib BLAS successful.')
        print('(log is in ',fulllog,')')

        # move librefblas.a to the lib directory
        shutil.copy('librefblas.a',os.path.join(self.config.prefix,'lib/librefblas.a'))

        # set framework variables to point to the freshly installed BLAS library
        self.config.blaslib  = '-L'+os.path.join(self.config.prefix,'lib')+' -lrefblas '
        os.chdir(savecwd)
        
if __name__ == '__main__':
    sys.path.insert(0, '../')
    import configure
    #from .dmft_install import Dmft_install
    from dmft_install import Dmft_install
    config = configure.Config((1, 0, 0))
    dmft_install = Dmft_install([], config)
    blas = Blas(config, dmft_install)
