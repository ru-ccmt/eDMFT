#! /usr/bin/env python
# -*- coding: utf-8 -*-

###
#
# @file gsl.py
#
#  DMFT is a software package provided by Rutgers Univiversity,
#  the State University of New Jersey
#
# @version 1.0.0
# @author Kristjan Haule and Viktor Oudovenko
# @date 2016-02-15
#
###
from install_script.iutils import writefile, shellcmd, delfiles, downloader, geturl, getURLName,includefromlib
import sys
import os
import urllib.request, urllib.parse, urllib.error
import shutil
#from . import framework
import install_script.iframework as framework
import re

class Gsl(framework.Framework):

    """ This class takes care of the libgsl. """
    def __init__(self, config, dmft):
        print("\n"+"="*50)
        print("  GSL installation/verification")
        print("="*50)

        self.config   = config
        self.downcmd  = dmft.downcmd
        #self.prefix   = dmft.prefix
        self.dmft     = dmft
        self.downgsl = dmft.downgsl
        self.gslurl  = "ftp://ftp.gnu.org/gnu/gsl/"+self.gslversion+".tar.gz"
                        #ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz

        self.dmft.verbose = 1

        if self.downgsl == 2:
            self.down_install_gsl()

        if(self.config.gsl == ""):
          if (os.path.isfile(os.path.join(self.config.prefix,'gsl/lib/libgsl.a')) or\
              os.path.isfile(os.path.join(self.config.prefix,'gsl/lib64/libgsl.a'))):
                self.set_gsl()

        ret = self.check_gsl()

        if ret != 0:
            if self.downgsl == 1:
                self.down_install_gsl()
            else:
                if not os.path.isfile(os.path.join(self.config.prefix,'gsl/lib/libgsl.a')):
                   print("""
 Please provide a working GSL library using --gsl. 
 If the GSL library is not awailable in the system, the GSL library can be 
 automatically downloaded and installed by adding the --downgsl flag.

 What do you want to do ?
    - s : specify the path and library if you have it
    - d : download and install the GSL library.
    - q : quit to download and install manually the GSL.
    - i : ignore are proceed
                   """)
                   answer = input(">[q] ")
                   if answer == "d":
                       self.down_install_gsl()
                   elif answer == "s":
                       self.config.gsl = input("> ")
                       ret = self.check_gsl()
                       if ret!=0:
                           sys.exit()
                   elif answer == "i":
                       pass
                   else:
                       sys.exit()
                else:
                   print("Netlib Lapack library is already installed at "+os.path.join(self.config.prefix,'lib/liblapack.a'))
                   print("Do you want to try it? (t)  or proceed without testing (p) or quit (q) ?")
                   answer = input(">[q] ")
                   if answer == "t":
                      self.config.gsl = '-L'+os.path.join(self.config.prefix,'gsl/lib')+' -lgsl '
                      self.check_gsl()
                   elif answer == "p":
                      exit
                   else:
                      sys.exit()

    def check_gsl(self):
        """ This function simply generates a C program
            that contains few GSL calls routine and then
            checks if compilation, linking and execution are succesful"""
        

        sys.stdout.flush()
        code="""#include <gsl/gsl_rng.h>
         const gsl_rng_type * T = gsl_rng_default; 
         gsl_rng * r = gsl_rng_alloc (T); 
         int main(void)
         {
           gsl_rng_env_setup();
           gsl_rng_set (r,1); 
           return 0;
         }
        """
        writefile('tmpc.cc',code)
        
        if self.config.gsl == "" or self.config.gsl is None:  # just trying default == -lgsl
            self.config.gsl='-lgsl'

        if self.config.gslinc == "" or self.config.gslinc is None:
            self.config.gslinc = includefromlib(self.config.gsl)
        ccomm = self.config.cxx+' '+self.config.gslinc+ ' -o tmpc '+'tmpc.cc '+self.config.gsl
        print('checking with:', ccomm)
        (output, error, retz) = shellcmd(ccomm)
        
        print("Checking if provided GSL works...", end=' ')
        if(retz != 0):
            if self.dmft.verbose:
                print('\n\nlibgsl: provided GSL cannot be used! aborting...')
                print('error is:\n','*'*50,'\n',ccomm,'\n',error.decode('UTF-8'),'\n','*'*50)
                return -1   
            else:
                print("no")
                return -1   
                #sys.exit()

        comm = './tmpc'
        (output, error, retz) = shellcmd(comm)
        if(retz != 0):
            if self.dmft.verbose:
                print('\n\nlibgsl: provided GSL cannot be used! aborting...')
                print('error is:\n','*'*50,'\n',comm,'\n',error.decode('UTF-8'),'\n','*'*50)
                print(retz)
                return -1   
            else:
                print("no")
                return -1   
            # sys.exit()

        if self.config.gslinc=='':
            # It worked, but we do not know the include files, hence figuring it out
            ccomm = self.config.cc+' -E  tmpc.cc |grep gsl | grep include'
            (output, error, retz) = shellcmd(ccomm)
            #print 'output=', output
            # compiler output in lines
            lines = output.decode('UTF-8').split('\n')
            incl={}
            for i,line in enumerate(lines):
                dat=line.split()
                for d in dat:
                    m = re.search('gsl',d) # take out the directory
                    if m is not None:
                        incl[os.path.dirname(d[1:-1])]=True
                    
            for inc in list(incl.keys()):
                self.config.gslinc += ' -I'+inc[:-4] # path has extra gsl. Take it out
        
        delfiles(['tmpc.cc','tmpc'])
        print('yes')

        return 0;


    def down_install_gsl(self):
        print("The GSL library is being installed.")
        sys.stdout.flush()

        savecwd = os.getcwd()

        # creating the build,lib and log dirs if don't exist
        if not os.path.isdir(os.path.join(self.config.prefix,'gsl')):
            os.mkdir(os.path.join(self.config.prefix,'gsl'))

        if not os.path.isdir(os.path.join(os.getcwd(),'log')):
            os.mkdir(os.path.join(os.getcwd(),'log'))

        # Check if gsl.tgz is already present in the working dir
        # otherwise download it
        if not os.path.isfile(os.path.join(os.getcwd(),getURLName(self.gslurl))):
            print("Downloading GSL ...", end=' ')
            #downloader(self.lapackurl,self.downcmd)
            #urllib.urlretrieve(self.gslurl, "gsl.tgz")
            geturl(self.gslurl, "gsl.tgz")

            print("done")

        # unzip and untar
        os.chdir('download')
        print('Unzip and untar GSL...', end=' ')
        comm = 'tar zxf gsl.tgz '
        (output, error, retz) = shellcmd(comm)
        if retz:
            print('\n\nlibgsl: cannot unzip '+self.gslversion+'.tgz')
            print('stderr:\n','*'*50,'\n',comm,'\n',error.decode('UTF-8'),'\n','*'*50)
            sys.exit()

        print('done')

        # change to GSL dir
        os.chdir(os.path.join(os.getcwd(), self.gslversion))

        # compile and generate library
        print('Configure  GSL...', end=' ')
        sys.stdout.flush()
        comm = './configure CC='+self.config.cc+'  --prefix='+os.path.join(self.config.prefix,'gsl')
        (output, error, retz) = shellcmd(comm)
        if retz:
            print("\n\nlingsl: cannot configure GSL")
            print("stderr:\n","*"*50,"\n",comm,'\n',error.decode('UTF-8'),"\n","*"*50)
            sys.exit()

        log = output.decode('UTF-8')+error.decode('UTF-8')

        # write log on a file
        #log = log+output+error
        fulllog = os.path.join(savecwd,'log/log.gsl')
        writefile(fulllog, log)
        print('Configuration of GSL  successful.')
        print('(log is in ',fulllog,')')

        # compile and generate library
        print('Compile and generate GSL...', end=' ')
        sys.stdout.flush()
        comm = self.make+' -j4; '+self.make+' install'
        (output, error, retz) = shellcmd(comm)
        if retz:
            print("\n\nlingsl: cannot compile GSL")
            print("stderr:\n","*"*50,"\n",comm,'\n',error.decode('UTF-8'),"\n","*"*50)
            sys.exit()

        log = output.decode('UTF-8')+error.decode('UTF-8')

        # write the log on a file
        #log = log+output+error
        fulllog = os.path.join(savecwd,'log/log.gsl')
        writefile(fulllog, log)
        print('Installation of GSL successful.')
        print('(log is in ',fulllog,')')

        # move libcblas.a to the lib directory
        #shutil.copy('libtmg.a',os.path.join(self.config.prefix,'gsl/libtmg.a'))

        # set framework variables to point to the freshly installed GSL library
        self.config.gsl = '-L'+os.path.join(self.config.prefix,'gsl')+' -lgsl '

        os.chdir(savecwd)

        # Check if the installation is successful
        self.dmft.verbose = 1
       # self.check_gsl()
        self.dmft.verbose = 0


    def set_gsl(self):
        # set framework variables to point to  installed GSL library
        
        gsllibpath = os.path.join(self.config.prefix,'gsl')+"/lib"
        if(os.path.isdir(gsllibpath)):
              self.config.gsl = '-L'+gsllibpath+' -lgsl -lgslcblas '
        else:
            gsllibpath = os.path.join(self.config.prefix,'gsl')+"/lib64"
            self.config.gsl = '-L'+gsllibpath+' -lgsl -lgslcblas '


            
if __name__ == '__main__':
    sys.path.insert(0, '../')
    import configure
    from install_script.dmft_install import Dmft_install
    config = configure.Config((1, 0, 0))
    dmft_install = Dmft_install([], config)
    gsl = Gsl(config, dmft_install)
