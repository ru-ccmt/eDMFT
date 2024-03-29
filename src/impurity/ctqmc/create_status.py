#!/usr/bin/env python
import re, glob, sys
import subprocess

def str3(i):
    if i<10:
        return '00'+str(i)
    elif i<100:
        return '0'+str(i)
    else:
        return str(i)
    
if len(sys.argv)<2:
    print 'Give number of needed status files!';
    sys.exit(1)
    
N = int(sys.argv[1])
olds = sorted(glob.glob('status.*'))

if len(olds)==0:
    print 'No status file found. Can not multiply if no status file is present!'
    sys.exit(1)


for i in range(N):
    old_file = olds[i % len(olds)]
    new_file = 'status.'+str3(i)
    if old_file!=new_file:
        cmd = 'cp '+old_file+ ' '+new_file
        print cmd,
        out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        print out, err
print
    
    
    


