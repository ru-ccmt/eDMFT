#!/usr/bin/env python
import sys

def OneLine(f1, f2):
    linea = f1.readline().split()
    lineb = f2.readline().split()

    if (linea != lineb):
        print('>', linea)
        print('<', lineb)
    return linea

def NumLine(linea):
    sdat=[]
    for c in range(4):
        w = linea[3+19*c:3+19*(c+1)].strip()
        if (len(w)>0): sdat.append(float(w))
    return sdat

def NumLine2(linea):
    fi = [0,0,0]
    for c in range(3):
        fi[c] = int(linea[3+5*c:3+5*(c+1)])
    fn = [float(linea[18:18+19]),float(linea[18+19:18+2*19])]
    return (fi, fn)


files = []
for arg in sys.argv:
    files.append(arg)

f1 = open(files[1], 'r')
f2 = open(files[2], 'r')


f1.readline()
f2.readline()
f1.readline()
f2.readline()
f1.readline()
f2.readline()

line = OneLine(f1, f2) # Atomnumber
#print('line=', line)
jatom = int(line[2])


for iatom in range(100):
    
    line = OneLine(f1, f2) # Number of LM
    nlm = int(line[3])
    print('jatom=', jatom, 'nlm=', nlm)
    
    
    for ilm in range(nlm):
    
        for i in range(2):
            f1.readline()
            f2.readline()
            
        line = OneLine(f1, f2)  # CLM(R) for L, M
        
        f1.readline()
        f2.readline()
        
        dsum=0.
        for i in range(1000):
            
            linea = f1.readline()
            lineb = f2.readline()
    
            numa = NumLine(linea)
            numb = NumLine(lineb)
            #print 'a=', numa
            #print 'b=', numb
            
            for i in range(len(numa)):
                dsum += abs(numa[i]-numb[i])
            #print dsum, '   ', linea
                
            if (len(numa)<4 and len(numb)<4): break
            
        
        #print "at=%d ilm=%d diff=%g" % (jatom, ilm, dsum)
        if (dsum>1e-16):
            print("at=%d ilm=%d diff=%g" % (jatom, ilm, dsum))
    
    
    for i in range(6):
        f1.readline()
        f2.readline()
        
    line = OneLine(f1, f2)
    if (line[0]!='ATOMNUMBER'): break
    jatom = int(line[2])

#print line

f1.readline()
f2.readline()

line = OneLine(f1,f2)
NPW = int(line[0])

wisum=0
wdsum=0
for i in range(NPW):
    linea = f1.readline()
    lineb = f2.readline()
    #2071 FORMAT(3X,3I5,2E19.12)

    (fia, fna) = NumLine2(linea)
    (fib, fnb) = NumLine2(lineb)

    isum=0
    for l in range(3): isum += abs(fia[l]-fib[l])
    dsum=0
    for l in range(2): dsum += abs(fna[l]-fnb[l])

    wisum += isum
    wdsum += dsum
    #print isum, dsum

print('Interstitsials', wisum, wdsum)


    


    

