#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
"""Script cleans eDMFT runs.
 The temporary and unimportant files for analisys are removed (check fremove and cremove lists below)
 Most of the files that might be relevant but are not very important in analisys are kept in a tarball, 
 for example, in imp? directories we produce impfiles.tgz and in other directories we produce wfiles.tgz.

"""
import re,glob,os,sys,copy,shutil
import subprocess
from utils import shellcmd

# Careful: these names to be matched are regular expressions, not wildcarts used in terminal
# files to be removed 
fremove= ['^dmft$', '^lapw1$', '^lapw2$', '^lapwso$', 'broyd.*', 'dmft\d?.error.*', 'lapw\d?.error','lcore.error','mixer.error']
cremove=['.energy_\d*', '.vectorso[nd]?[up]?[.gz]?', '.vector[nd]?[up]?[.gz]?', '.vns_old']
# files to be tarballed
f_tar=['dmft\d+_info.out', 'sig.inp\.\d+.\d+']
c_tar=['.outputdmf\d+']

# for impurity directory, which files should be tarballed.
impfiles = ['ctqmc', 'Probability\.dat\..*', 'status\..*', 'check.dat.*', 'g_hb\d\..*', 'new.cix', 'nohup_imp.out.*', 'nohup.out', 'Delta.tau\..*', 'Sigma\..*', 'oca_log\..*', 'Ac\.inp\..*', 'Aw.out.*', 'Gcoeff.dat', 'Sw.dat', 'gs_qmc.dat', 'Gt.dat', 'Gw.dat']


if len(sys.argv)<2:
    print('Give directory name')
    sys.exit(1)
else:
    drn = sys.argv[1]
    print('Directory is', drn)


case = os.path.splitext(os.path.split(glob.glob(drn+"/*.struct")[0])[1])[0]
print('case=', case)

# Files to remove 
fremovea = fremove + [case+fname for fname in cremove]
# Files to tar
to_tar = f_tar + [case+fname for fname in c_tar]

# First, find all files on the directory tree
impzip={}
remove = {}
tarz={}
for root, dirs, files in os.walk(drn):
    #print('root=', root, 'dirs=', dirs)
    #print('files=', files)
    
    if re.match(drn+'/?imp.', root) is not None:
        for name in files:
            #print( root, name )
            for imp in impfiles:
                if re.match(imp, name):  # should be zipped
                    if root in impzip:
                        impzip[root].append(name)
                    else:
                        impzip[root]=[name]
                    break    
    else: # not impurity directory, could be current or subdirectory
        for name in files: # all existing files
            for fname in fremovea: # files to remove
                #print('checking for match1', fname, name)
                if re.match(fname, name):  # should be removed
                    if root in remove:
                        remove[root].append(name)
                    else:
                        remove[root]=[name]
                    break
            for fname in to_tar:
                #print('checking for match2', fname, name)
                if re.match(fname, name):  # should be tagz
                    if root in tarz:
                        tarz[root].append(name)
                    else:
                        tarz[root]=[name]
                    break
            longfname = os.path.join(root, name)
            size = os.stat(longfname).st_size
            if size==0:
                if root in remove:
                    remove[root].append(name)
                else:
                    remove[root]=[name]

## this just prints what will happen
#for root,files in impzip.items():
#    for i,name in enumerate(files):
#        print(i, 'itar', root, name)
#for root,files in tarz.items():
#    for i,name in enumerate(files):
#        print('tar',i, root, name)
#for root,files in remove.items():
#    for i,name in enumerate(files):
#        print('rm', i, root, name)
#sys.exit(0)

N=20  # note that system command has problems with too many arguments, hence we rm 20 files at the time.
curd = os.getcwd()
for root,files in impzip.items():
    os.chdir(root)
    (out, err, retz) = shellcmd('tar -czvf impfiles.tgz '+' '.join(files))
    print(err.decode('UTF-8'),end='')
    print(out.decode('UTF-8'),end='')
    #print(proc.returncode)
    for i in range(int(len(files)/N)+1):
        pack = files[i*N:(i+1)*N]
        (out, err, retz) = shellcmd('rm '+' '.join(pack))
        print(err.decode('UTF-8'),end='')
        print(out.decode('UTF-8'),end='')
    os.chdir(curd)
    
for root,files in tarz.items():
    os.chdir(root)
    (out, err, retz) = shellcmd('tar -czvf wfiles.tgz '+' '.join(files))
    print(err.decode('UTF-8'),end='')
    print(out.decode('UTF-8'),end='')
    #print(proc.returncode)
    for i in range(int(len(files)/N)+1):
        pack = files[i*N:(i+1)*N]
        (out, err, retz) = shellcmd('rm '+' '.join(pack))
        print(err.decode('UTF-8'),end='')
        print(out.decode('UTF-8'),end='')
    os.chdir(curd)

for root,files in remove.items():
    os.chdir(root)
    for i in range(int(len(files)/N)+1):
        pack = files[i*N:(i+1)*N]
        (out, err, retz) = shellcmd('rm '+' '.join(pack))
        print(err.decode('UTF-8'),end='')
        print(out.decode('UTF-8'),end='')
    os.chdir(curd)
