#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
#
# @file utils.py
#
#  DMFT is a software package provided by Rutgers Univiversity,
#  the State University of New Jersey
#
# @version 1.0.0
# @author Kristjan Haule and Viktor Oudovenko
# @date 2016-02-15
#
###

import os
import re
import shutil
import subprocess
import select


def writefile(fname, fill):
    """ writes the file fname with content fill """
    fp = open(fname,'w')
    fp.write(fill)
    fp.close()

def delfiles(lst):
    """ deletes a list of files """
    for i in lst:
        if(os.path.isfile(i)):
            os.remove(i)

def openPipe(command):

    pipe = None
    if hasattr(popen2, 'Popen3'):
        pipe   = popen2.Popen3(command, 1)
        input  = pipe.tochild
        output = pipe.fromchild
        err    = pipe.childerr
    else:
        (input, output, err) = os.popen3(command)

    return (input, output, err, pipe)

class NullFile:
    def write(self, s): pass
    def flush(self): pass

def shellcmd(command, outfile=None):
    """ runs a shell command """
    if not outfile:
      outfile = NullFile()

    proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = proc.communicate()
    return (out, err, proc.returncode)

def getURLName(url):
    directory=os.curdir

    name="%s%s%s" % (
        directory,
        os.sep,
        url.split("/")[-1]
    )

    return name



def downloader(uri,cmd):
    """ downloads the content of an URL """

    savedir = os.getcwd()
    downdir = savedir+"/download"
    if(not os.path.isdir(downdir)):
        print("Creating directory", downdir)
        os.mkdir(downdir)
    os.chdir(downdir)

    name = getURLName(uri)
    try:
        if(os.path.isfile(downdir+"/"+name)):
            print("Package "+name+" already downloaded")
        elif(cmd == 'urllib2'):
            import urllib.request, urllib.error, urllib.parse
            url = urllib.request.urlopen(uri)
            f = open(name,'w')
            for line in url.readlines():
                f.write(line)
            url.close()
            f.close()
        elif(cmd == 'wget'):
            comm = 'wget '+uri
            (output, error, retz) = shellcmd(comm)
        else:
            raise
    except:
        print(" ")
        print("=================================================================================")
        print("ERROR: Cannot download "+name)
        print("Make sure the network is reachable.")
        print("If you still have troubles, you can manually download "+name+" from this URL:")
        print(uri)
        print("into download directory:")
        print(os.getcwd())
        print("")

    os.chdir(savedir)
    shutil.copy('download/'+name, './')

def geturl(uri,packagename):
    """ downloads the content of an URL """

    savedir = os.getcwd()
    downdir = savedir+"/download"
    if(not os.path.isdir(downdir)):
        print("Creating directory", downdir)
        os.mkdir(downdir)

    os.chdir(downdir)
    #name = getURLName(uri)
    name = packagename
    try:
        if(os.path.isfile(downdir+"/"+name)):
            print("Package "+name+" already downloaded")
        else:
            import urllib.request, urllib.parse, urllib.error
            urllib.request.urlretrieve(uri,name)
            #url = urllib2.urlopen(uri)
    except:
        print(" ")
        print("=================================================================================")
        print("ERROR: Cannot download "+name)
        print("Make sure the network is reachable.")
        print("If you still have troubles, you can manually download "+name+" from this URL:")
        print(uri)
        print("into download directory:")
        print(os.getcwd())
        print("")

    os.chdir(savedir)
    # shutil.copy('download/'+name, './')



def fixpaths(inpath):
    lst = inpath.split(" ")

    outpath = ""

    if ((len(lst) == 1) and (inpath[0] != "download") and (inpath[0] != '-') and (inpath[0] != '')):
        outpath = os.path.abspath(inpath)
        return outpath
    else:
        for i in lst:
            if re.search("^-L",i):
                p = "-L"+os.path.abspath(i[2:])
            else:
                p = i

            outpath = outpath+p+" "	

    return outpath

def includefromlib(inpath):
    lst = inpath.split(" ")
    outpath = ""
    #print('includefromlib=', lst)
    for i in lst:
        if re.search("^-L",i):
            #print('YES: ', i)
            # gets absolute path
            abspath = os.path.abspath(i[2:])
            # trying to remove trailing '/lib' if possible
            m = re.match('(.*)\/lib(64)?$',abspath)
            if m is not None:
                p = m.group(1)+'/include'
            else:
                p = abspath+"/../include"
            outpath += ' -I'+p
    return outpath

#def includefromlib(inpath):
#
#    lst = inpath.split(" ")
#    outpath = ""
#
#    for i in lst:
#        if re.search("^-L",i):
#           p = os.path.abspath(i[2:])+"/../include"
#           outpath = p	
#
#    return outpath
