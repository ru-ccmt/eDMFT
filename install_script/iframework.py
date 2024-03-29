#! /usr/bin/env python
# -*- coding: utf-8 -*-

###
#
# @file framework.py
#
#  DMFT is a software package provided by Rutgers Univiversity,
#  the State University of New Jersey
#
# @version 1.0.0
# @author Kristjan Haule and Viktor Oudovenko
# @date 2016-02-15
#
###

#from .utils import shellcmd, writefile, delfiles, fixpaths, includefromlib
from install_script.iutils import shellcmd, writefile, delfiles, fixpaths, includefromlib
import sys
import os
import getopt
import string
import subprocess

class Framework:
    """ This class is framework for DMFT package installation. """
    #set default values
    prefix        = None                      # The install directory
    build         = "./build"                 # The build directory
    make          = "make"                    # the "make" command
    downcmd       = ""                        # the command used to download stuff
    ranlib        = ""                        # the ranlib command
    downblas      = 0                         # whether or not to download reference Blas
    downfftw      = 0                         # whether or not to download reference FFTW
    downgsl       = 0                         # whether or not to download reference GSL
    downlapack    = 0                         # whether or not to download reference Lapack

    lapversion   = "lapack-3.6.0"
    fftwversion   = "fftw-3.3.4"
    gslversion   = "gsl-1.16"

    versions     = ("1.0.0")
    clean        = 0
    src          = 0
    installerdir = ""
    verbose      = 0
    #nbcores      = 2

    def __init__(self, argv, config):
        #print "*-+-*"*16
        print("Setting up the framework")

        self.config = config

        self.config.PIC = '-fPIC'

        if config.prefix==None:
            if 'WIEN_DMFT_ROOT' in os.environ and os.environ['WIEN_DMFT_ROOT']!='':
                config.prefix=os.environ['WIEN_DMFT_ROOT']
            else:
                config.prefix=os.getcwd()+'/install/'
                
        self.config.prefix = os.path.abspath(config.prefix)
        
        # parse input arguments
        self.parse_args(argv)

        if self.clean:
            self.cleanup()
            sys.exit()

        if ((str.upper(self.config.compiler) == 'INTEL')):
            if self.config.cc==None or self.config.cc == '':
              self.config.cc = 'icc'
              print("""
You specified """+self.config.compiler+""" as your main compiler 
and left  --cc  parameter empty, assuming it equal to """+self.config.cc+".\n")
            
            if self.config.fc==None or self.config.fc == '':
              self.config.fc = 'ifort'
              print("""
You specified """+self.config.compiler+""" as your main compiler 
and left  --fc  parameter empty, assuming it equal to """+self.config.fc+".\n")
            
            if self.config.cxx==None or self.config.cxx == '':
              self.config.cxx = 'icpc'
              print("""
You specified """+self.config.compiler+""" as your main compiler 
and left  --cxx  parameter empty, assuming it equal to """+self.config.cxx+".\n")
        
            if self.config.preproc==None or self.config.preproc == '':
              self.config.preproc = 'fpp'
              print("""
You specified """+self.config.compiler+""" as your main compiler 
and left  --preproc  parameter empty, assuming it equal to """+self.config.preproc+".\n")


        if ((str.upper(self.config.compiler) == 'GNU')):
            if self.config.cc==None or self.config.cc == '':
                self.config.cc = 'gcc'
                print("""
You specified """+self.config.compiler+""" as your main compiler 
and left  --cc  parameter empty, assuming it equal to """+self.config.cc+".\n")

            if self.config.fc==None or self.config.fc == '':
              self.config.fc = 'gfortran'
              print("""
You specified """+self.config.compiler+""" as your main compiler 
and left  --fc  parameter empty, assuming it equal to """+self.config.fc+".\n")
            
            if self.config.cxx==None or self.config.cxx == '':
              self.config.cxx = 'g++'
              print("""
You specified """+self.config.compiler+""" as your main compiler 
and left  --cxx  parameter empty, assuming it equal to """+self.config.cxx+".\n")
            
            if self.config.preproc==None or self.config.preproc == '':
              self.config.preproc = 'cpp -traditional-cpp'
              print("""
You specified """+self.config.compiler+""" as your main compiler 
and left  --preproc  parameter empty, assuming it equal to """+self.config.preproc+".\n")

        if self.config.cc == '':
            #  check if no C compiler is provided
            print("""
C compiler is required to compile DMFT software.  
Please use  --cc flag or edit configure.py file.""")
            sys.exit()

        if self.config.fc == "":
            #  chekc if Fortran compiler is provided
            print("""
Fortran compiler is required to compile DMFT software.
Please use  --fc flag or edit configure.py file.""")
            sys.exit()

        if self.config.cxx == "":
            #  check if  C++ compiler is provided 
            print("""
C++ compiler is required to compile DMFT software.
Please use  --cxx flag or edit configure.py file.""")
            sys.exit()

        if((self.config.mpi_define != "") and ( (self.config.pcc =="")  or (self.config.pfc == "") or (self.config.pcxx == "") )):
            #  check if  MPI compilers are provided 
            print("""
MPI Fortran, C and C++  compilers are required to compile DMFT software 
with non empty mpi_define="""+self.config.mpi_define+""" flag.
Please use  --pfc, --pcc and  --pcxx flags or edit configure.py file 
and specify missing parameters.\n""")
            sys.exit()

        if(self.config.mpi_define == ""):
          self.config.pcc = self.config.cc 
          self.config.pfc = self.config.fc 
          self.config.pcxx = self.config.cxx 

        if self.config.prefix == "" :
            self.config.prefix = "./install"
        self.config.prefix = fixpaths(self.config.prefix)
        if(not os.path.isdir(self.config.prefix)):
            print("Creating directory", self.config.prefix)
            try:
                os.mkdir(self.config.prefix)                
                print("Directory created successfully")
            except OSError as e:
                print("Failed to create directory:", e)
        print('Install directory is...', self.config.prefix)

        if self.build == "":
            self.build = "./"
        self.build = fixpaths(self.build)

        #print('HERE: self.build=', self.build)
        if(not os.path.isdir(self.build)):
            print("Creating directory",self.build)
            try:
                os.mkdir(self.build)
                print("Directory created successfully")
            except OSError as e:
                print("Failed to create directory:", e)
        print('Build directory is...', self.build)

        self.installerdir = os.getcwd()
        os.chdir(self.build)

        # CUSTOM CHECKS
        self.check_cc()
        self.check_fc()
        self.set_ranlib()
        #self.set_download()
        self.detect_compilers()
        self.detect_blaslibs()
        #self.check_linking()
        #if self.testing:
        #self.set_nbcores()
        
        print('C compiler is...       ', self.config.cc)
        print('C flags are...         ', self.config.cflags)
        print('Fortran compiler is... ', self.config.fc)
        print('Fortran flags are...   ', self.config.fflags)
        #print 'Ar flags are...        ', self.config.arflags
        
        if (self.downblas > 1) :
            print('BLAS library will to be downloaded and installed )')
        else:
            print('BLAS library is... ', self.config.blaslib)
        
        if (self.downfftw > 1) :
            print('FFTW library will be downloaded and installed ')
        else :
            print('FFTW library is... ', self.config.fftwlib)
        
        if (self.downgsl > 1) :
            print('GSL library will be downloaded and installed ')
        else :
            print('GSL library is...  ', self.config.gsl)

        if (self.downlapack == 2) :
            print('Lapack library will to be downloaded and installed )')
        elif (self.config.lapacklib == ""):
            if (self.downlapack == 1):
                print('LAPACK library is... will check if it is part of Blas Library call and download if it is not')
            else:
                print('LAPACK library is... will check if it is part of Blas library call')
        else:
            self.downlapack = -1
            print('LAPACK library is...', self.config.lapacklib)

        return

    def usage(self):
          print("*-+-*"*16)
          print("""
   DMFT configuration script version %d.%d.%d.
   The script will help you to create makefile & compile DMFT code.
   Please provide as much information as you can using flags listed below.
   Alternatively, you can edit "configure.py" file in the current directory.

   -h or --help        : Display this help and exit

   --prefix=[DIR]      : Install files in DIR [%s]
   --build=[DIR]       : Software building DIR [%s]. Contains logs, downloads and builds.

   --fc=[CMD]          : Fortran compiler. [%s]
   --cc=[CMD]          : C compiler. [%s]
   --cxx=[CMD]         : C++ compiler. [%s]

   --cflags=[FLAGS]    : Flags for the C compiler [%s]
   --fflags=[FLAGS]    : Flags for the Fortran compiler [%s]
   --ompflag=[FLAGS]   : Flags for openmp compiler [%s]

   --pfc=[CMD]         : MPI Fortran compiler. [%s]
   --pcc=[CMD]         : MPI C compiler. [%s]
   --pcxx=[CMD]        : MPI C++ compiler. [%s]

   --blaslib=[LIB]     : BLAS library [%s]
   --lapacklib=[LIB]   : Lapack library (if it's not included in BLAS) [%s]
   --fftwlib=[LIB]     : FFTW library [%s]
   --gsl=[LIB]         : GSL library [%s]

   --downblas          : Download and install BLAS.
   --downlapack        : Download and install LAPACK.
   --downfftw          : Download and install FFTW.
   --downgsl           : Download and install GSL.
   --downall           : Download and install all missing external libraries.
                         If you don't have access to network please provide 
                         the following packages in directory %s/download:
        		
			 http://netlib.org/blas/blas.tgz		->	[blas.tgz]
        		 http://www.netlib.org/lapack/lapack-%s.tgz	->	[lapack.tgz]
        		 http://fftw.org/fftw-%s.tar.gz		->	[fftw.tgz]
         		 ftp://ftp.gnu.org/gnu/gsl/gsl-%s.tar.gz	->	[gsl.tgz]

   --clean             : cleans up the installation directory.
   """ % (self.config.version[0], self.config.version[1], self.config.version[2], 
            self.config.prefix,os.path.abspath(self.build),
            self.config.fc, self.config.cc,self.config.cxx,
            self.config.cflags,self.config.fflags,self.config.ompflag,
            self.config.pfc, self.config.pcc,self.config.pcxx,
            self.config.blaslib, self.config.lapacklib,self.config.fftwlib,self.config.gsl, 
            self.build, "3.6.0","3.3.4","1.16"))
          print("*-+-*"*16)

    def parse_args(self, argv):
        """ Parse input argument to get compilers etc. from command line. """

        #if len(argv) == 1: print("Let's go")
        try:
          opts, args = getopt.getopt(argv[1:], "?hvp:b:n:",
                                         ["help", "compiler=", "prefix=", "build=",
                                          "cc=", "fc=", "cxx=", 
                                          "cflags=", "fflags=", "oflags=", "gflags=",
                                          "pcc=", "pfc=", "pcxx=", 
                                          "blaslib=", "lapacklib=",  "fttwlib=", "gsl=",
                                          "mpi_define=", "make=",
                                          "downblas", "downlapack","downgsl","downfftw",
                                          "downall", "verbose", "clean", "src"])

        except getopt.error as msg:
          print(msg)
          print("for help use --help")
          sys.exit(2)

        if len(args) > 0 :
            print('Unknown arguments : ', args)
            print("for help use --help")
            sys.exit(2);

        # process options
        for o, a in opts:
            if o in ("-h", "--help"):
                self.usage()
                sys.exit(0)
            else:
                if o == '--clean':
                    self.clean = 1
                    return
                elif o in ('-p', '--prefix'):
                    self.prefix = a
                elif o in ('-b', '--build'):
                    self.build = a
                elif o == '--cflags':
                    self.config.cflags = a
                elif o=='--fflags':
                    self.config.fflags = a
                elif o=='--oflags':
                    self.config.oflags = a
                elif o=='--gflags':
                    self.config.gflags = a
                elif o=='--make':
                    self.make = a
                elif o=='--compiler':
                    self.config.ompiler = a
                elif o=='--cc':
                    self.config.cc = a
                elif o=='--fc':
                    self.config.fc = a
                elif o=='--cxx':
                    self.config.cxx = a
                elif o=='--pcc':
                    self.config.cc = a
                elif o=='--pfc':
                    self.config.fc = a
                elif o=='--pcxx':
                    self.config.cxx = a
                elif o == '--downblas':
                    self.downblas = 2
                elif o == '--downlapack':
                    self.downlapack = 2
                elif o == '--downfftw':
                    self.downfftw = 2
                elif o == '--downgsl':
                    self.downgsl = 2
                elif o == '--downall':
                    self.downblas   = max(1, self.downblas  )
                    self.downlapack = max(1, self.downlapack)
                    self.downfftw   = max(1, self.downfftw  )
                elif o == '--blaslib':
                    self.config.blaslib = fixpaths(a)
                elif o == '--gsl':
                    self.config.gsl = fixpaths(a)
                elif o == '--lapacklib':
                    self.config.lapacklib = fixpaths(a)
                elif o == '--fttwlibs':
                    self.config.fttwlibs = fixpaths(a)
                elif o == '--mpi_define':
                    self.config.mpi_define = fixpaths(a)
                elif (o in ('-v', '--verbose')):
                    self.verbose = 1
                else :
                    print("Unknown option : ", o)
                    sys.exit()

        # Set correctly downloads
        if (((self.config.blaslib == "") and (self.downblas > 0))
            or (self.config.blaslib == "download") ):
            self.config.blaslib = ""
            self.downblas = max(1, self.downblas)
        else :
            self.downblas = 0

        if (((self.config.fftwlib == "") and (self.downfftw > 0))
            or (self.config.fftwlib == "download" )):
            self.config.fftwlib = ""
            self.downfftw = max(1, self.downfftw)
        else :
            self.downfftw = 0

        if (((self.config.lapacklib == "") and (self.downlapack > 0))
            or (self.config.lapacklib == "download" )):
            self.config.lapacklib = ""
            self.downlapack = max(1, self.downlapack)
        else :
            self.downlapack = 0

        if (((self.config.gsl == "") and (self.downgsl > 0))
            or (self.config.gsl == "download" )):
            self.config.gsl = ""
            self.downgsl = max(1, self.downgsl)
        else :
            self.downgsl = 0


    #if (self.config.ldflags_fc == "") and (self.config.ldflags_c):
    #   self.config.ldflags_fc = self.config.ldflags_c

    def check_cc(self):
        """ checking if cc works """
        # simply generates a C program containing a couple of calls
        # to MPI routines and checks if the compilation and execution
        # are succesful
        print('Checking if cc works...', end=' ')
        sys.stdout.flush()
        # generate
        writefile('tmpc.c',"""
            #include <stdio.h>
            int main(int argc, char **argv){
            int iam;
            fprintf(stdout, \"success\" );fflush(stdout);
            return 0;
            }\n""")

        # compile
        #ccomm = self.config.cc+" "+self.config.cflags+" "+self.config.ldflags_c+" -o tmpc "+os.path.join(os.getcwd(),"tmpc.c")
        import re
        cflags_ = re.sub('-O3', '-O0', self.config.cflags) # -O3 would take too much time. Do not need optimization here....
        ccomm = self.config.cc+" "+cflags_+"  -o tmpc "+os.path.join(os.getcwd(),"tmpc.c")
        (output, error, retz) = shellcmd(ccomm)

        if retz:
            print('\n\nCOMMON: C compiler not working! aborting...')
            print('stderr:\n','*'*50,'\n',error.decode('UTF-8'),'\n','*'*50)
            sys.exit()

        # run
        comm = './tmpc'
        (output, error, retz) = shellcmd(comm)
        if retz:
            print('\n\nCOMMON: cc not working! aborting...')
            print('error is:\n','*'*50,'\n',error.decode('UTF-8'),'\n','*'*50)
            sys.exit()

        # cleanup
        delfiles(['tmpc.c','tmpc'])
        print('yes')
        return 0;

    def check_fc(self):
        """ check if the Fortran compiler works """
        # simply generates a F77 program and checks if the compilation and execution
        # are succesful
        print("Checking if the Fortran compiler works...", end=' ')
        sys.stdout.flush()
        # generate
        writefile("tmpf.f","""
      program ftest
      integer i
      print*,'success'
      stop
      end\n""")

        # compile
        #fcomm = self.config.fc+' '+self.config.fcflags+" "+self.config.ldflags_fc+' -o tmpf '+'tmpf.f'
        fcomm = self.config.fc+' '+self.config.fflags+'  -o tmpf '+'tmpf.f'
        (output, error, retz) = shellcmd(fcomm)

        if retz:
            print('\n\nCOMMON: the Fortran compiler is not working! aborting...')
            print('error is:\n','*'*50,'\n',error.decode('UTF-8'),'\n','*'*50)
            sys.exit()

        # run
        comm = './tmpf'
        (output, error, retz) = shellcmd(comm)
        if retz:
            print('\n\nCOMMON: the Fortran compiler is not working! aborting...')
            print('error is:\n','*'*50,'\n',error.decode('UTF-8'),'\n','*'*50)
            sys.exit()

        # cleanup
        delfiles(['tmpf.f','tmpf','tmpf.o'])
        print('yes')

        return 0;

    def set_mangling(self):
        """ Sets the INTFACE variable in Bmake.inc """
        # This one generates a program equivalent to that in BLACS/INSTALL
        # that checks the mangling in FORTRAN function symbols
        print('Setting Fortran mangling...', end=' ')
        sys.stdout.flush()
        writefile('tmpf.f',"""
      program intface
      external c_intface
      integer i
      call c_intface(i)
      stop
      end\n""")
        writefile('tmpc.c',"""
      #include <stdio.h>
      void c_intface_(int *i){fprintf(stdout, \"-DADD_\");fflush(stdout);}
      void c_intface(int *i){fprintf(stdout, \"-DNOCHANGE\");fflush(stdout);}
      void c_intface__(int *i){fprintf(stdout, \"-DfcIsF2C\");fflush(stdout);}
      void C_INTFACE(int *i){fprintf(stdout, \"-DUPCASE\");fflush(stdout);}\n""")

        ccomm = self.config.cc+' '+self.config.cflags+' -c tmpc.c -o tmpc.o'
        fcomm = self.config.fc+' '+self.config.fflags+'  tmpf.f tmpc.o -o xintface'

        (output, error, retz) = shellcmd(ccomm)
        if retz:
            print('\n\nCOMMON: in set_mangling: cannot compile')
            print('error is:\n','*'*50,'\n',error.decode('UTF-8'),'\n','*'*50)
            #sys.exit()

        (output, error, retz) = shellcmd(fcomm)
        if retz:
            print('\n\nCOMMON: in set_mangling: cannot compile')
            print('error is:\n','*'*50,'\n',error.decode('UTF-8'),'\n','*'*50)
            #sys.exit()

        comm = os.path.join(os.getcwd(),'xintface')
        (output, error, retz) = shellcmd(comm)
        if retz:
            print('\n\nCOMMON: in set_mangling: cannot run xintface')
            print('error is:\n','*'*50,'\n',error.decode('UTF-8'),'\n','*'*50)
            #sys.exit()

        self.mangling = output.decode('UTF-8')
        delfiles(['xintface', 'tmpf.f', 'tmpf.o', 'tmpc.c', 'tmpc.o'])

        print(self.mangling)
        return 1;


    def check_linking(self):
        """ Check if C main can be linked to Fortran subroutine """

        # This one checks if the linking command works out of the box or
        # if any specific flag is required. For example if the linker if the
        # Intel FORTRAN compiler, then the "-nofor_main" is usually required.
        # This function only checks if linker works but does not automatically
        # detect the required flags
        print('Checking loader...', end=' ')
        sys.stdout.flush()
        writefile('tmpf.f',"""
      subroutine fsub()
      write(*,*)'success'
      stop
      end\n""")
        writefile('tmpc.c',"""
      #if defined ADD_
      #define fsub fsub_
      #elif defined NOCHANGE
      #define fsub fsub
      #elif defined fcIsF2C
      #define fsub fsub_
      #elif defined UPCASE
      #define fsub FSUB
      #endif
      void main(){
      fsub();}\n""")

        #ccomm = self.config.cc+' '+self.config.ccflags+' '+self.mangling+' -c -o tmpc.o tmpc.c'
        #fcomm = self.config.fc+' '+self.config.fcflags+' -c -o tmpf.o tmpf.f'
        #lcomm = self.config.fc+' '+self.config.ldflags_fc+' '+self.config.ld_fcmain+' -o lnk tmpf.o tmpc.o'
        ccomm = self.config.cc+' '+self.config.cflags+' -c -o tmpc.o tmpc.c'
        fcomm = self.config.fc+' '+self.config.fflags+' -c -o tmpf.o tmpf.f'
        lcomm = self.config.fc+'   -o lnk tmpf.o tmpc.o'

        (output, error, retz) = shellcmd(ccomm)
        if retz:
            print('\n\nCOMMON: in check_linking: cannot compile')
            print('command is: ',ccomm)
            print('error is:\n','*'*50,'\n',error.decode('UTF-8'),'\n','*'*50)
            #sys.exit()

        (output, error, retz) = shellcmd(fcomm)
        if retz:
            print('\n\nCOMMON: in check_linking: cannot compile')
            print('command is: ',fcomm)
            print('error is:\n','*'*50,'\n',error.decode('UTF-8'),'\n','*'*50)
            #sys.exit()

        (output, error, retz) = shellcmd(lcomm)
        if retz:
            print("""\n\nCOMMON: in check_linking: cannot link
            Cannot link a C main program to a Fortran77 subroutine
            Make sure that the appropriate flags are passed to the linker.""")
            print('command is: ',lcomm)
            print('error is:\n','*'*50,'\n',error.decode('UTF-8'),'\n','*'*50)
            #sys.exit()


        delfiles(['lnk', 'tmpf.f', 'tmpf.o', 'tmpc.c', 'tmpc.o'])

        print('works')
        return 1;



    def set_download(self):
        """ Figures out how to download files """
        print('Setting download command...')
        wget = 0
        urllib = 0
        if urllib == 0:
            # if urllib2 is not present checks if wget is present
            # in the PATH and if yes it sets the download command
            # to be wget
            print("Checking availablility of wget...", end=' ')
            path=str(os.getenv('PATH')).split(os.pathsep)
            for i in path:
                if (os.path.isfile(os.path.join(i,'wget'))):
                    print("available")
                    wget = 1
                    break
            if wget:
                # test wget
                print("Testing wget...", end=' ')
                comm = 'wget --tries=2 --timeout=5 http://www.netlib.org/lapack/index'
                (output, error, retz) = shellcmd(comm)
                if(retz != 0):
                    print('not working.')
                    wget = -1
                else:
                    print("working")
                    self.downcmd="wget"
                    os.remove("index")
                    return
            else:
                # wget not available
                print("not available")
                wget=0


    def set_ranlib(self):
        """ Sets the ranlib command """
        # Some systems don't have the ranlib command (e.g. SGIs).
        # In the case where ranlib is not present in the PATH,
        # echo is used instead of ranlib
        print("Setting ranlib command...", end=' ')

        path=str(os.getenv('PATH')).split(os.pathsep)
        for i in path:
            if os.path.isfile(os.path.join(i,'ranlib')):
                self.config.ranlib=os.path.join(i,'ranlib')
                print(self.config.ranlib)
                return

        for i in path:
            if os.path.isfile(os.path.join(i,'echo')):
                self.config.ranlib=os.path.join(i,'echo')
                print(self.config.ranlib)
                return

    def detect_compilers(self):
        """ Tries to detect the compilers type """
        # By users experience it is known which compiler flags are required
        # in some cases. This function tries to detect which compilers are used
        # and sets the flags accordingly

        print('Detecting Fortran compiler...', end=' ')
        if self.fc_is_intel():
            # The Intel FORTRAN compiler requires -nofor_main flag
            # for the linking and the -mp flag to maintain the
            # floating-point precision
            self.config.ldflags    += ' -free -no-prec-div -pc80 $(OPENMP) '
            # self.config.fflags    += ' -diag-disable vec -fltconsistency -fp_port'
            # self.config.ldflags_c  += ' '   # used to link
            # self.config.ldflags_fc += ' '
            self.config.ld_fcmain   = ' -nofor_main'
            #self.config.noopt      += ' -mp'
            self.config.f2pyflag   += ' --fcompiler=intelem '
            # self.testing = 0; # Cannot compile lintest with fc_main option
            print('Intel')
        elif self.fc_is_gnu():
            print('GNU')
            self.config.ldflags    += ' -ffree-form -ffree-line-length-none $(OPENMP) '
            self.config.f2pyflag   += ' --fcompiler=gnu95 '
            # self.config.ld_fcmain   = ''
        elif self.fc_is_xlf():
            self.config.fflags    += ' -qstrict -qthreaded'
            # self.config.ld_fcmain   = ''
            print('IBM')
        elif self.fc_is_pgi():
            #  self.config.ldflags_c  += ''
            self.config.ldflags_fc += ''
            #  self.config.ld_fcmain   = ' -Mnomain'
            # # self.testing = 0; # Cannot compile lintest with fc_main option
        else:
            self.config.compiler = "Unknown"
            print('unknown')

        print('Detecting C compiler...', end=' ')
        if self.cc_is_intel():
            self.config.compiler = "Intel"
            self.config.cflags    += ''#' -diag-disable vec'
            print('Intel')
        elif self.cc_is_gnu():
            self.config.compiler = "GNU"
            print('GNU')
        elif self.cc_is_xlc():
            self.config.compiler = "XLC"
            self.config.cflags    += ''#' -qstrict -qthreaded'
            print('IBM')
        elif self.cc_is_pgi():
            self.config.compiler = "PGI"
            print('PGI')
        else:
            print('unknown')

        print('Selected C compiler flags: '+self.config.cflags)
        print('Selected Fortran compiler flags: '+self.config.fflags)
        #print 'Selected loader flags (C main): '+self.config.ldflags_c
        #print 'Selected loader flags (Fortran main): '+self.config.ldflags_fc
        return


    def fc_is_intel(self):
        comm = self.config.fc+' -V'
        (output, error, retz) = shellcmd(comm)
        #isifort = string.find(error,'Intel(R) Fortran')
        isifort = error.find(b'Intel(R) Fortran')
        if isifort != -1:
            return 1

        return 0

    def fc_is_intel_new(self):
        if subprocess.call([self.config.fc,"-V"]) == 0:
             print("Found C")
             return 1

        return 0


    def fc_is_gnu(self):
        comm = self.config.fc+' -v'
        (output, error, retz) = shellcmd(comm)
        #isifort = string.find(error,'gcc')
        isifort = error.find(b'gcc')
        if isifort != -1:
            return 1

        return 0

    def fc_is_gnu_new(self):
        if subprocess.call([self.config.fc,"--version"]) == 0:
            print("\n Found GNU C")
            return 1

        return 0


    def fc_is_xlf(self):
        comm = self.config.fc
        (output, error, retz) = shellcmd(comm)
        #isifort = string.find(output,'xlf')
        isifort = output.find(b'xlf')
        if isifort != -1:
            return 1

        return 0
    def fc_is_pgi(self):
        comm = self.config.fc+' -V'
        (output, error, retz) = shellcmd(comm)
        #isifort = string.find(error,'pgifc')
        isifort = error.find(b'pgifc')
        if isifort != -1:
            return 1
        #isifort = string.find(error,'pgif95')
        isifort = error.find(b'pgif95')
        if isifort != -1:
            return 1
        #isifort = string.find(error,'Portland')
        isifort = error.find(b'Portland')
        if isifort != -1:
            return 1
        #isifort = string.find(output,'pgifc')
        isifort = output.find(b'pgifc')
        if isifort != -1:
            return 1
        #isifort = string.find(output,'pgif95')
        isifort = output.find(b'pgif95')
        if isifort != -1:
            return 1
        #isifort = string.find(output,'Portland')
        isifort = output.find(b'Portland')
        if isifort != -1:
            return 1

        return 0

    def cc_is_intel(self):
        comm = self.config.cc+' -V'
        (output, error, retz) = shellcmd(comm)
        #isifort = string.find(error,'Intel(R) C')
        isifort = error.find(b'Intel(R) C')
        if isifort != -1:
            return 1

        return 0

    def cc_is_gnu(self):
        comm = self.config.cc+' -v'
        (output, error, retz) = shellcmd(comm)
        #isifort = string.find(error,'gcc')
        isifort = error.find(b'gcc')
        if isifort != -1:
            return 1

        return 0

    def cc_is_xlc(self):
        comm = self.config.cc+' -qversion'
        (output, error, retz) = shellcmd(comm)
        #isifort = string.find(output,'XL')
        isifort = output.find(b'XL')
        if isifort != -1:
            return 1

        return 0

    def cc_is_pgi(self):
        comm = self.config.cc+' -V'
        (output, error, retz) = shellcmd(comm)
        #isifort = string.find(error,'pgicc')
        isifort = error.find(b'pgicc')
        if isifort != -1:
            return 1
        #isifort = string.find(error,'Portland')
        isifort = error.find(b'Portland')
        if isifort != -1:
            return 1
        #isifort = string.find(output,'pgicc')
        isifort = output.find(b'pgicc')
        if isifort != -1:
            return 1
        #isifort = string.find(output,'Portland')
        isifort = output.find(b'Portland')
        if isifort != -1:
            return 1

        return 0

    def set_nbcores(self):
        print('Setting environment variable...', end=' ')
        os.environ['OMP_NUM_THREADS']='1'
        os.environ['GOTO_NUM_THREADS']='1'
        os.environ['MKL_NUM_THREADS']='1'
        print('done.')
        if "SC_NPROCESSORS_ONLN" in os.sysconf_names:
           ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
           print('Checking number of cores available on the machine: ', ncpus)
           if int(ncpus)<int(self.nbcores):
              print('---> INFO: you do not have' , self.nbcores , 'cores on that machine! (will use the default 2)')
              self.nbcores=2
           print('Number of cores to be used for DMFT testing: ',  self.nbcores)
           if int(self.nbcores)==1:
              print('WARNING: DMFT is meant to run on more than one core')

        return 0



    def detect_blaslibs(self):
        """ Tries to detect type of Blas libraris installed in the system"""

        print('\n Detecting MKL library...', end=' ')
        if self.fc_is_intel():
            self.config.fflags    += '' 
            #self.config.fflags    += ' -diag-disable vec -fltconsistency -fp_port'
            self.config.ld_fcmain   = ''# '-nofor_main'
            #self.config.noopt      += '' #' -mp'
            self.config.f2pylib    += " --opt='-fast' " 
           # self.testing = 0; # Cannot compile lintest with fc_main option
            print('Intel')
        elif self.fc_is_gnu():
            print('GNU')
            self.config.ld_fcmain   = ''
        elif self.fc_is_xlf():
            self.config.fcflags    += ' ' #'-qstrict -qthreaded'
            self.config.ld_fcmain   = ''
            print('IBM')
        elif self.fc_is_pgi():
            self.config.ldflags_c  += ''
            self.config.ldflags_fc += ''
            self.config.ld_fcmain   = '' #' -Mnomain'
           # self.testing = 0; # Cannot compile lintest with fc_main option
        else:
            self.config.compiler = "Unknown"
            print('unknown')
        if self.config.f2pylapack==None or self.config.f2pylapack=='':
            self.config.f2pylapack='--link-lapack_opt'
            
    def resume(self):
        print("""

CONGRATULATIONS!
Installation of DMFT software has been completed.
Results are stored in : """+self.config.prefix+"""
******************************************************
WARNING: Copy the bin directory ("""+self.config.prefix+""")
         to desired location, and set "WIEN_DMFT_ROOT" 
         environment variable to point to this location.

         Please make sure that WIEN2K program is installed
         and corresponding eviroment variables are set.
         Also make sure that numpy and scipy Python 
         packages are available in the system.
******************************************************
""")
        return 0


    def cleanup(self):
        " Cleans up the installer directory "

        print("Cleaning up...", end=' ')
        sys.stdout.flush()
        
        DIR = os.path.abspath(os.path.dirname(__file__))
        maindir = os.path.normpath(os.path.join(DIR, "../src"))
        builddir= os.path.normpath(os.path.join(DIR, "../build"))
        
        print('builddir=', builddir)

        os.chdir(maindir) 
        comm = 'make clean'
        retz=subprocess.call(comm,shell=True,stdout=sys.stdout,stderr=sys.stderr)
        #(output, error, retz) = shellcmd(comm)
        
        comm = 'rm -rf '+builddir
        retz=subprocess.call(comm,shell=True,stdout=sys.stdout,stderr=sys.stderr)
        
        comm = 'rm -rf '+self.config.prefix
        retz=subprocess.call(comm,shell=True,stdout=sys.stdout,stderr=sys.stderr)
        
        print("done.")
