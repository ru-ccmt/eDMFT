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
        if not dat:
            continue
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

                   
    if not work:
        return

    # Install complete rank files before retiring any previous layout.
    # dftKS uses the common _processes_ file, which otherwise takes
    # precedence over these newly converted WIEN2k rank files in dmft2.
    from pathlib import Path
    from tempfile import TemporaryDirectory
    written = set()
    with TemporaryDirectory(prefix='.dmft-vector-maps-', dir='.') as staging:
        for prc, rows in sorted(work.items()):
            if prc<1:
                raise ValueError('Invalid processor index in WIEN2k .processes')
            dest = Path('_processes_'+str(prc-1))
            with (Path(staging)/dest).open('w') as output:
                for i, nkp, start, names in rows:
                    if nkp<0 or not all(names):
                        raise ValueError('Incomplete vector assignment in WIEN2k job '+str(i))
                    # Older readers consume only the first seven fields.
                    print(i, nkp, start, *('"'+name+'"' for name in names),
                          'nk_total', nkstart, file=output)
            written.add(dest)
        for dest in sorted(written):
            os.replace(Path(staging)/dest, dest)
    for old in Path('.').glob('_processes_*'):
        if re.fullmatch(r'_processes_[0-9]+',old.name) and old not in written:
            old.unlink()
    Path('_processes_').unlink(missing_ok=True)

if __name__ == '__main__':
    
    so = False
    
    inso=glob.glob('*.inso')
    if len(inso) and os.path.getsize(inso[0]):
        so = True
    
    convert(so, sys.stdout)
    
