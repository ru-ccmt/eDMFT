#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule  
import glob
import os, sys

def RemoveSomeFiles(fls):
    if fls:
        prompt = "Remove files (y/n): ", fls
        userin = input(prompt).lower().strip()
        if userin.lower() == 'y':
            for file in fls: os.remove(file)


# Removes all empty files
for root,dirs,files in os.walk('.'):
    for f in files:
        file_path = os.path.join(root, f) #Full path to the file
        size = os.path.getsize(file_path) #pass the full path to getsize()
        if size == 0:
            print('removing: ', file_path)
            os.remove(file_path)

fls = glob.glob('dmft')+glob.glob('dmft2')+glob.glob('imp.0/ctqmc')
RemoveSomeFiles(fls)

fls = glob.glob('*.vector_*')
RemoveSomeFiles(fls)

fls = glob.glob('imp.0/status.*')
RemoveSomeFiles(fls)

fls = glob.glob('_processes_*')
RemoveSomeFiles(fls)

fls = glob.glob('lapw1_*')
RemoveSomeFiles(fls)

fls = glob.glob('*.broyd*')
RemoveSomeFiles(fls)

fls = glob.glob('*.output1_*')
RemoveSomeFiles(fls)
