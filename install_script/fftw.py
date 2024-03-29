#! /usr/bin/env python
# -*- coding: utf-8 -*-
" Module which checks fftw"
import re
###
#
# @file fftw.py
#
#  DMFT is a software package provided by Rutgers Univiversity,
#  the State University of New Jersey
#
# @version 1.0.0
# @author Kristjan Haule and Viktor Oudovenko
# @date 2016-02-15
#
###

from install_script.iutils import writefile, shellcmd, delfiles, downloader, geturl, getURLName, includefromlib # here removed .utils
import sys
import os
import urllib.request, urllib.parse, urllib.error
import shutil
#from . import framework
import install_script.iframework as framework

class Fftw(framework.Framework):
    """ This class takes care of the libfftw. """
    def __init__(self, config, dmft):
        print("\n"+"="*50)
        print("  FFTW installation/verification")
        print("="*50)

        self.config   = config
        self.downcmd  = dmft.downcmd
        #self.prefix   = dmft.prefix
        self.dmft     = dmft
        self.downfftw = dmft.downfftw
        self.fftwurl  = "http://fftw.org/"+self.fftwversion+".tar.gz"
        #http://fftw.org/fftw-3.3.4.tar.gz
        self.dmft.verbose = 1

        if self.downfftw == 2:
            self.down_install_fftw()

        if(self.config.fftwlib == ""):
          if (os.path.isfile(os.path.join(self.config.prefix,'fftw/lib/libfftw3.a')) or\
             os.path.isfile(os.path.join(self.config.prefix,'fftw/lib64/libfftw3.a'))):
                self.set_fftlib()

        ret = self.check_fftw()
        if ret != 0:
            if self.downfftw == 1:
                self.down_install_fftw()
            else:
                if not (os.path.isfile(os.path.join(self.config.prefix,'fftw/lib/libfftw3.a')) or os.path.isfile(os.path.join(self.config.prefix,'fftw/lib64/libfftw3.a'))):
                   print(""" Please provide a working FFTW library using --fftwlib. 
 If the FFTW library is not awailable in the system, the FFTW library can be 
 automatically downloaded and installed by adding the --downfftw flag.
                    
 What do you want to do ?
   - s : specify the path and library if you have it
   - d : download and install the FFTW library.
   - q : quit to download and install manually the FFTW.
   - i : ignore and proceed
""")
                   answer = input(">[q] ")
                   if answer == "d":
                       self.down_install_fftw()
                   elif answer =="s":
                       self.config.fftwlib = input("> ")
                       ret = self.check_fftw()
                       if ret!=0:
                           sys.exit()
                   elif answer =="i":
                       pass
                   else:
                       sys.exit()
                else:
                   print("FFTW library is available at "+os.path.join(self.config.prefix,'fftw/lib/libfftw3.a'))
                   print("Do you want to try it? (t)  or proceed without testing (p) or quit (q) ?")
                   answer = input(">[q] ")
                   if answer == "t":
                       self.set_fftlib()
                       self.check_fftw()
                   elif answer == "p":
                       exit
                   else:
                       sys.exit()

    def check_fftw(self):
        """ This function simply generates a C program
            that contains few FFTW calls routine and then
            checks if compilation, linking and execution are succesful"""
        
        sys.stdout.flush()
        code="""
        #include<fftw3.h>
        #include<complex.h>
        int main(void)
        {
            int N;
            int i;
            N=1;
            fftw_complex *in, *out;
            fftw_plan my_plan;
            in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
            out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
            my_plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(my_plan);
            fftw_destroy_plan(my_plan);
            fftw_free(in);
            fftw_free(out);
            return 0;
         }
        """
        writefile('tmpc.c',code)

        if self.config.fftwlib == "" or self.config.fftwlib is None:
            self.config.fftwlib='-lfftw3'  # just trying default"

        #print('config.fftwinc=', self.config.fftwinc, 'config.fftwlib=', self.config.fftwlib)
        if self.config.fftwinc == "" or self.config.fftwinc is None:
            self.config.fftwinc = includefromlib(self.config.fftwlib)
        ccomm = self.config.cc+' '+self.config.fftwinc+' -o tmpc  tmpc.c '+self.config.fftwlib
        print('checking with:', ccomm)
        
        (output, error, retz) = shellcmd(ccomm)
        
        print("Checking if provided FFTW works...", end=' ')

        if(retz != 0):
            if self.dmft.verbose:
                print('\n\nlibfftw: provided FFTW cannot be used! aborting...')
                print('error is:\n','*'*50,'\n',ccomm,'\n',error.decode('UTF-8'),'\n','*'*50)
            else:
                print("no")
            return 1

        comm = './tmpc'
        (output, error, retz) = shellcmd(comm)
        if(retz != 0):
            if self.dmft.verbose:
                print('\n\nlibfftw: provided FFTW cannot be used! aborting...')
                print('error is:\n','*'*50,'\n',comm,'\n',error.decode('UTF-8'),'\n','*'*50)
                print(retz)
            else:
                print("no")
            return 1

        if self.config.fftwinc=='':
            # No include file given, hence checking which includes were used when compilation succeded
            ccomm = self.config.cc+' -E '+self.config.fftwinc+' tmpc.c |grep fftw | grep include'
            (output, error, retz) = shellcmd(ccomm)
            # compiler output in lines
            lines = output.decode('UTF-8').split('\n')
            incl={}
            for i,line in enumerate(lines):
                dat=line.split()
                for d in dat:
                    m = re.search('fftw',d) # take out the directory
                    if m is not None:
                        incl[os.path.dirname(d[1:-1])]=True # path has extra "xxx" so take them out
        
            for inc in list(incl.keys()):
                self.config.fftwinc += ' -I'+inc 
        delfiles(['tmpc.c','tmpc'])
        print('yes')
        return 0;


    def down_install_fftw(self):
        print("""
        The FFTW library is being installed.
        """)
        sys.stdout.flush()

        savecwd = os.getcwd()

        # creating the build,lib and log dirs if don't exist
        if not os.path.isdir(os.path.join(self.config.prefix,'fftw')):
            os.mkdir(os.path.join(self.config.prefix,'fftw'))

        if not os.path.isdir(os.path.join(os.getcwd(),'log')):
            os.mkdir(os.path.join(os.getcwd(),'log'))

        # Check if fftw.tgz is already present in the working dir
        # otherwise download it
        if not os.path.isfile(os.path.join(os.getcwd(),"fftw.tgz")):
            print("Downloading FFTW ...", end=' ')
            #downloader(self.lapackurl,self.downcmd)
            #urllib.urlretrieve(self.fftwurl, "fftw.tgz")
            geturl(self.fftwurl, "fftw.tgz")
            print("done")

        # unzip and untar
        os.chdir('download')
        print('Unzip and untar FFTW...', end=' ')
        comm = 'tar zxf fftw.tgz '
        (output, error, retz) = shellcmd(comm)
        if retz:
            print('\n\nlibfftw: cannot unzip '+self.fftwversion+'.tgz')
            print('stderr:\n','*'*50,'\n',comm,'\n',error.decode('UTF-8'),'\n','*'*50)
            sys.exit()

        print('done')

        # change to FFTW dir
        os.chdir(os.path.join(os.getcwd(), self.fftwversion))

        # compile and generate library
        print('Configure  FFTW...', end=' ')
        sys.stdout.flush()
        if(self.config.ompflag == "" ):
           comm = './configure MPICC='+self.config.pcc+' CC='+self.config.cc+' F77='+self.config.fc+\
           ' --enable-mpi  --enable-threads  --enable-shared --prefix='+os.path.join(self.config.prefix,'fftw')
        else:
           comm = './configure MPICC='+self.config.pcc+' CC='+self.config.cc+' F77='+self.config.fc+\
           ' --enable-mpi  --enable-threads  --enable-shared  -enable-openmp --prefix='+os.path.join(self.config.prefix,'fftw')

        (output, error, retz) = shellcmd(comm)
        if retz:
            print("\n\nlinfftw: cannot configure FFTW")
            print("stderr:\n","*"*50,"\n",comm,'\n',error.decode('UTF-8'),"\n","*"*50)
            sys.exit()


        # write log on a file
        log = output.decode('UTF-8')+error.decode('UTF-8')
        fulllog = os.path.join(savecwd,'log/log.fftw')
        writefile(fulllog, log)
        print('Configuration of FFTW  successful.')
        print('(log is in ',fulllog,')')

        # compile and generate library
        print('Compile and generate FFTW...', end=' ')
        sys.stdout.flush()
        comm = self.make+' -j4; '+self.make+' install'
        (output, error, retz) = shellcmd(comm)
        if retz:
            print("\n\nlinfftw: cannot compile FFTW")
            print("stderr:\n","*"*50,"\n",comm,'\n',error.decode('UTF-8'),"\n","*"*50)
            sys.exit()


        # write the log on a file
        log = output.decode('UTF-8')+error.decode('UTF-8')
        fulllog = os.path.join(savecwd,'log/log.fftw')
        writefile(fulllog, log)
        print('Installation of FFTW successful.')
        print('(log is in ',fulllog,')')

        # move libcblas.a to the lib directory
        #shutil.copy('libtmg.a',os.path.join(self.config.prefix,'fftw/libtmg.a'))

        self.set_fftlib()

        os.chdir(savecwd)

        # Check if the installation is successful
        self.dmft.verbose = 1
        self.check_fftw()
        self.dmft.verbose = 0


    def set_fftlib(self):

        # set framework variables to point to the freshly installed FFTW library
        
        fftwlibpath = os.path.join(self.config.prefix,'fftw')+"/lib"
        if(os.path.isdir(fftwlibpath)):
            if(self.config.ompflag == ''):
              self.config.fftwlib = '-L'+fftwlibpath+'             -lfftw3 '
            else:
              self.config.fftwlib = '-L'+fftwlibpath+' -lfftw3_omp -lfftw3 '
        else:
            fftwlibpath = os.path.join(self.config.prefix,'fftw')+"/lib64"
            if(self.config.ompflag == ''):
              self.config.fftwlib = '-L'+fftwlibpath+'             -lfftw3 '
            else:
              self.config.fftwlib = '-L'+fftwlibpath+' -lfftw3_omp -lfftw3 '

if __name__ == '__main__':
  sys.path.insert(0, '../')
  import configure
  from dmft_install import Dmft_install
  config = configure.Config((1, 0, 0))
  dmft_install = Dmft_install([], config)
  fft = Fftw(config, dmft_install)
  #fft.check_fftw()
