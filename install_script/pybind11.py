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
import subprocess

class Pybind11(framework.Framework):

    """ This class takes care of the libgsl. """
    def __init__(self, config, dmft):
        print("\n"+"="*50)
        print("  Pybind11 installation/verification")
        print("="*50)

        self.config   = config
        self.downcmd  = dmft.downcmd
        self.dmft     = dmft
        self.dmft.verbose = 1
        
        #if self.config.pybind11 == "" or self.config.pybind11 is None:  # just trying default:
            #self.config.pybind11='`python -m pybind11 --includes` -shared -std=c++11 -fPIC'
        from sysconfig import get_paths
        DIR = os.path.abspath(os.path.dirname(__file__))
        include_dir = os.path.normpath(os.path.join(DIR, "../src/includes"))
        extra_config = '-I'+include_dir             # pybind11 includes
        extra_config+= ' -I'+get_paths()['include'] # python includes
        extra_config+= ' -shared -std=c++11 -fPIC '
        if self.config.pybind11 is not None:
            self.config.pybind11 = extra_config+self.config.pybind11
        else:
            self.config.pybind11 = extra_config
        
        ret = self.check_pybind()
        
        if ret != 0:
            print("""
 Please install pybind11 module for your python installation.
 It could be done through 
            `conda install -c conda-forge pybind11`
 or
            `pip install pybind11`
        """)
            sys.exit(1)
            
    def check_pybind(self):
        """ This function simply generates a C program
            that contains few GSL calls routine and then
            checks if compilation, linking and execution are succesful"""
        

        sys.stdout.flush()
        code="""#include "pybind11/pybind11.h"
        #include "pybind11/numpy.h"
        #include "pybind11/stl.h"
        #include <cstdint>
        #include <vector>
        
        namespace py = pybind11;
        void st(py::array_t<double>& data, const std::vector<int>& ext)
        {
          auto dat = data.mutable_unchecked<2>();
          int i0 = ext[0], i1=ext[1];
        }
        
        PYBIND11_MODULE(simple,m){
          m.doc() = "pybind11 wrap checking";
          m.def("simple", &st);
        }
        """
        writefile('simple.cc',code)



        ccomm = self.config.cxx+' '+self.config.pybind11+ ' simple.cc -o simple.so'
        #print('checking with:', ccomm)
        print("Checking if pybind11 works...", end=' ')
        #(output, error, retz) = shellcmd(ccomm)
        retz=subprocess.call(ccomm,shell=True,stdout=sys.stdout,stderr=sys.stderr)         
        if(retz != 0):
            if self.dmft.verbose:
                print('\n\npybind11: provided pybind11 cannot be used! aborting...')
                print('configuration for pybind11 is:',self.config.pybind11)
                print('The following command is used to compile')
                print(ccomm)
                #print('error is:\n'+'*'*50,'\n',ccomm,'\n',error,'\n'+'*'*50)
                return -1   
            else:
                print("no")
                return -1   
                #sys.exit()
        print(" cmpile yes; ", end='')
        comm = 'python -c "import simple;"'
        #(output, error, retz) = shellcmd(comm)
        retz=subprocess.call(comm,shell=True,stdout=sys.stdout,stderr=sys.stderr)
        if(retz != 0):
            if self.dmft.verbose:
                print('\n\npybind11: provided pybind11 cannot be used! aborting...')
                print('execution command:', comm)
                print('Did not succeed with return status=', retz)
                print('compilation for pybind11 is:',self.config.pybind11)
                print('The following command is used to compile')
                print(ccomm)
                print()
                print('You might try to compile end test the file "build/simple.cc" modifying compilation.')
                return -1
            else:
                print("no")
                return -1   
            # sys.exit()
        print(" exec yes; ", end='')
        
        delfiles(['simple.cc', 'simple.so'])
        print('yes')

        return 0

            
if __name__ == '__main__':
    sys.path.insert(0, '../')
    import configure
    from install_script.dmft_install import Dmft_install
    config = configure.Config((1, 0, 0))
    dmft_install = Dmft_install([], config)
    py11 = Pybind11(config, dmft_install)
