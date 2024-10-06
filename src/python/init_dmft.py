#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import optparse
import os, shutil
import subprocess
from numpy import zeros, log
import indmffile
from utils import DmftEnvironment, W2kEnvironment
import sys
from sigen import CreateIndmflFromIndmf

def print_prompt(text):
    print("----->", text)

if __name__ == '__main__':
    usage = """Program : inid_dmft.py
    
    Initializes the EDMFTF run by producing case.indmfl and case indmfi files.
    The first connects DMFT problem with the Kohn-Sham orbitals, 
          and the second connects DMFT with the impurity solver.
    
    Usage init_dmft.py  [-option val]
         
     where option's are:
       -h, --help     -- displays this messgae
       -ca  1,2,3     -- list of correlated atoms (which atoms from 
                                     case.struct are treated by DMFT)
       -ot  d,d,d     -- orbital type being correlated 
                         (which orbital is treated by DMFT, i.e, d or f)
       -qs  7,7,7     -- q-split described the symmetry of the orbital. 
                          (7 cubic, 2 real harmonics,..). 
                          For help execute -h qs
       -p  5          -- projector type (For help execute -h p )
       -hr -10,10     -- hybridization window in eV
       -a  i          -- real (r) or imaginary (i) axis calculation
       -us "1,2 3,4"  -- unique set of orbitals, i.e., atom 1 and 2 map to
                         an equivalent imputy problem, and 3,4 map to the
                         second equivalent impurity problem.
       -cl "1,2 3,4"  -- grouping of atoms into cluster for cluster-DMFT 
                         by default atoms are treated by single-site DMFT
    """
        
    inpt={}  # input from command line
    if len(sys.argv)>1:
        if sys.argv[1] in ['-h', '--help']:
            print(usage)
            if len(sys.argv)>2 and sys.argv[2] == 'qs':
                print(indmffile.qsplit_doc)
            if len(sys.argv)>2 and sys.argv[2] == 'p':
                print(indmffile.projector_doc)
            sys.exit(0)
        else:        
            for i in range(int((len(sys.argv)-1)/2)):
                inpt[sys.argv[2*i+1][1:]] = sys.argv[2*i+2]
            print('given options are', inpt)
    
    dmfenv = DmftEnvironment()
    w2kenv = W2kEnvironment()
    indmf = indmffile.Indmf(w2kenv.case)
    
    # First, check if case.indmf exists
    create_indmf = True
    if indmf.file_exists():
        prompt = "File `%s` already exists.  Do you want to create new inputs? (y/n): " % indmf.filename()
        userin = input(prompt).lower().strip()
        if userin.lower() == 'n':
            create_indmf = False
    
    # create new case.indmf (if desired)
    if create_indmf:
        indmf.user_input(inpt)
        indmf.write()
        fname = indmf.filename()
        fnametemp = fname + '_st'
    
        if inpt:
            finished_editing = True  # no editing when using command line
        else:
            finished_editing = False
        while not finished_editing:
            shutil.copy(fname, fnametemp)  # copy to temp file for user to edit
            print_prompt("Edit %s file." % (fnametemp,))
            subprocess.call(w2kenv.EDITOR.split() + [fnametemp])
    
            # read user-edited file in just to check for any syntax errors
            indmf_st = indmffile.Indmf(w2kenv.case)
            indmf_st.extn = 'indmf_st'
            try:
                indmf_st.read()
                finished_editing = True
                shutil.copy(fnametemp, fname)  # move it back            
            except:
                print_prompt("Your edits have introduced syntax errors in %s." % (fnametemp,))
                prompt = "Do you want to edit the file again? (y/n): "
                userin = input(prompt).lower().strip()
                if userin.lower() == 'n':
                    finished_editing = True
    
    # ask user if this is a SO run
    so_run = False
    if os.path.isfile(w2kenv.case+".inso") and os.path.getsize(w2kenv.case+".inso")>0 :
        print('Found '+w2kenv.case+'.inso file, hence assuming so-coupling exists. Switching -so switch!')
        so_run = True

    # parse main input file case.indmf
    try:
        indmf = indmffile.Indmf(w2kenv.case)
        indmf.read()
    except Exception as e:
        print("ERROR: reading file `%s.indmf` failed." % w2kenv.case) 
        print(e)
        exit(1)
        
    cmpProjector = os.path.exists(w2kenv.case+'.vsp') # w2k was already done
    inl = CreateIndmflFromIndmf(indmf, so_run, cmpProjector=cmpProjector)
    inl.write()
    
    ## run sigen to calculate siginds, cftrans and write output to case.indmfl file
    #sigen = os.path.join( dmfenv.ROOT, 'sigen.py' )
    ##sigen = './sigen.py'
    #args = ["--so"] if so_run else []
    #ireturn = subprocess.call([sigen] + args, shell=True,stdout=sys.stdout,stderr=sys.stdout)
    #if ireturn:
    #    print("Error executing sigen.py.")
    #    exit(1)
    
    # copy case.indmfl --> case.indmfl_st to allow user to make changes
    findmfl = w2kenv.case + '.indmfl'
    findmfltemp = findmfl + '_st'
    
    if inpt:
        finished_editing = True  # no editing when using command line
    else:
        finished_editing = False
    while not finished_editing:
        shutil.copy(findmfl, findmfltemp)
        print_prompt("Edit %s file." % (findmfltemp,))
        subprocess.call(w2kenv.EDITOR.split() + [findmfltemp])
    
        # read user-edited file in just to check for any syntax errors
        inl = indmffile.Indmfl(w2kenv.case)
        inl.extn = 'indmfl_st'
        try:
            inl.read()
            finished_editing = True
            shutil.copy(findmfltemp, findmfl)  # move it back
        except:
            print_prompt("Your edits have introduced syntax errors in %s." % (findmfltemp,))
            prompt = "Do you want to edit the file again? (y/n): "
            userin = input(prompt).lower().strip()
            if userin.lower() == 'n':
                finished_editing = True
    
    # Generate case.indmfi file -----------------------------------
    indmf.read()
    inl = indmffile.Indmfl(w2kenv.case)
    inl.read()
    
    f = open(w2kenv.case + '.indmfi', 'w')
    print(len(indmf.cixgrp), "  # number of sigind blocks", file=f)
    for iucp,icps in indmf.cixgrp.items():
        icp = icps[0]  # just take sigind of first correlated problem  (TODO: modify appropriately for broken symmetry runs)
        sigind = inl.siginds[icp]
        print(sigind.shape[0], "  # dimension of this sigind block", file=f)
        max_sigfig = 1 + int(log(max(sigind.flat))/log(10))
        format = '%' + str(max_sigfig) + 'd'
        for row in sigind:
            print(' '.join([format % elem for elem in row]), file=f)
    f.close()
    

