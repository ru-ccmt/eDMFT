#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import re, sys
import os, glob

def convert(so, fh_info):
    if not os.path.exists('.processes'):
        print('Can not find .processes. This must be serial run. I will not create _processes_x', file=fh_info)
        return 
    
    f = open('.processes')
    lines = f.readlines()
    f.close()
    
    machines=[]
    work={}
    nkstart=0
    for line in lines:
        dat = line.split()
        if dat[0]=='init:':
            machines.append(dat[1])
        if re.match('\s*\d+',dat[0]):
            VECN=['','','','']
    
            (i,nkp,nprc) = (int(dat[0]),int(dat[4]),int(dat[8]))
    
            if not so:
                fdef = open('lapw1_'+str(i)+'.def','r')
                for lin in fdef:
                    dat = lin.split(',')
                    if int(dat[0])==10:
                        m=re.search('.*[\'|\"](.*)_(\d+)',dat[1])
                        if m is None:
                            print('ERROR: Could not find vector file to macth in lapw1.def file')
                        VECN[0] = m.group(1)+'_'+m.group(2)
                        VECN[1] = m.group(1)+'dn_'+m.group(2)
                    if int(dat[0])==11:
                        m=re.search('.*[\'|\"](.*)_(\d+)',dat[1])
                        if m is None:
                            print('ERROR: Could not find vector file to macth in lapw1.def file')
                        VECN[2] = m.group(1)+'_'+m.group(2)
                        VECN[3] = m.group(1)+'dn_'+m.group(2)
                fdef.close()
            else:
                print('lapwso_'+str(i)+'.def', file=fh_info)
                fdef = open('lapwso_'+str(i)+'.def','r')
                for lin in fdef:
                    dat = lin.split(',')
                    if int(dat[0])==42: VECN[0]=dat[1].split("'")[1]
                    if int(dat[0])==41: VECN[1]=dat[1].split("'")[1]
                    if int(dat[0])==52: VECN[2]=dat[1].split("'")[1]
                    if int(dat[0])==51: VECN[3]=dat[1].split("'")[1]
                fdef.close()
            
            if nprc in work:
                work[nprc].append( (i,nkp, nkstart, VECN) )
            else:
                work[nprc]=[ (i,nkp,nkstart, VECN) ]
            nkstart = nkstart + nkp


    #for prc in sorted(work.keys()):
    #    for (i,nkp,nkstart,VECN) in work[prc]:
    #        print i, VECN[0]
    #        print i, VECN[1]
    #        print i, VECN[2]
    #        print i, VECN[3]
    #print '------------'

                   
    for prc in sorted(work.keys()):
        fo = open('_processes_'+str(prc-1),'w')
        for (i,nkp,nkstart,VECN) in work[prc]:
            print(i, nkp, nkstart, '"'+VECN[0]+'"', '"'+VECN[1]+'"', '"'+VECN[2]+'"', '"'+VECN[3]+'"', file=fo)
            #print  i, nkp, nkstart, '"'+VECN[0]+'"', '"'+VECN[1]+'"', '"'+VECN[2]+'"', '"'+VECN[3]+'"'
        
        

if __name__ == '__main__':
    
    so = False
    
    inso=glob.glob('*.inso')
    if len(inso) and os.path.getsize(inso[0]):
        so = True
    
    convert(so, sys.stdout)
    
