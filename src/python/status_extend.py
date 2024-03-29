#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import glob, os, sys

def giveString(i, clen=3):
    s = str(i)
    s = '0'*(3-len(s)) + s
    return s

fname='status.'

if len(sys.argv)<=1:
    print('Give desired number of status files!')
    sys.exit(1)

desired=int(sys.argv[1])
sln = len(fname)
status = glob.glob(fname+'*')

print('current status files: ', status)

ilast= int(status[-1][sln:])
for i in range(ilast+1,desired+1):
    ii = (ilast-(i-ilast-1))%(ilast+1)
    cmd = 'cp '+fname+giveString(ii)+' '+fname+giveString(i)
    print(cmd)
    stdin, stdout, stderr = os.popen3( cmd )
    print(stdout.read(), stderr.read(), end=' ')





