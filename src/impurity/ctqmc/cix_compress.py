#!/usr/bin/env python

from scipy import *
import sys, glob
import operator


def cmp_small_large(wlarge_small):
    nsize = len(wlarge_small)
    large_small = []
    for i in range(nsize): large_small.append(-1)
    ii=0
    for i in range(nsize):
        if len(wlarge_small[i])>0:
            large_small[i] = ii
            ii+=1
    
    
    small_large=[]
    for i in range(nsize):
        if len(wlarge_small[i])>0:
            small_large.append((i,wlarge_small[i]))
            
    return (small_large, large_small)

def union(data):
    " Takes a union of array or list"
    c = []
    for d in data:
        if d not in c:
            c.append(d)
    return c


def PrintHeader(small_large, large_small, start, sindex, Ene, Spin):
    for i in range(len(small_large)):
        ii = small_large[i][0]

        print "%3s %3s %3s %4s" % (i+1, start[ii][1], start[ii][2], start[ii][3]), 
        print "%3d   " % len(small_large[i][1]), #sdim[ii],
        for j in sindex[ii]:
            if j>=0: print "%3d " % (large_small[j]+1),
            else: print "%3d " % (j+1),
            
        for j in small_large[i][1]: print Ene[ii][j],
        for j in small_large[i][1]: print Spin[ii][j],
        print

def PrintHeader1(small_large, large_small, start, sindex, findex, Ene, Spin):
    for i in range(len(small_large)):
        ii = small_large[i][0]

        print "%3s %3s %3s %3s %4s" % (i+1, findex[i]+1, start[ii][2], start[ii][3], start[ii][4]), 
        print "%3d   " % len(small_large[i][1]), #sdim[ii],
        for j in sindex[ii]:
            if j>=0: print "%3d " % (large_small[j]+1),
            else: print "%3d " % (j+1),
            
        for j in small_large[i][1]: print Ene[ii][j],
        for j in small_large[i][1]: print Spin[ii][j],
        print

if __name__ == '__main__':

    cutof = 1e-6
    maxDim = 1000
    
    files=[]
    par=[]
    for a in sys.argv[1:]:
        if len(glob.glob(a))>0:
            files.append(a)
        else:
            par.append(a)

    if len(files)<2:
        print 'Need input parameters: <cix file> <probability>'
        sys.exit(0)
    if len(par)>0:
        cutof = float(par[0])

    
        
    fcix = files[0]
    fProb = files[1]

    
    fc = open(fcix, 'r')

    lines = fc.readlines()
    count = 2

    (Nm, nsize, N_ifl, max_size) = map(int, lines[count].split())
    
    count += 2
    Nbaths=0
    for b in range(N_ifl):
        bath = map(int, lines[count].split())
        if bath[0] != b : print 'something wrong in reading cix file ib!=ibath'
        dimbath = bath[1]
        Nbaths += dimbath
        count +=1

    #print Nbaths
    
    count += 3

    bath_header = count

    
    state=[]
    start=[]
    Ene=[]
    Spin=[]
    sdim=[]
    sindex=[]
    
    for i in range(nsize):
        state.append(lines[count])
        slines = lines[count].split()    

        start.append(slines[0:5])
        wsdim = int(slines[4])
        wsindex = map(int, slines[5:5+Nbaths])
        
        sdim.append(wsdim)

        #print 'b=', wsindex
        for iw in range(len(wsindex)): wsindex[iw]-=1
        #print 'a=', wsindex
        
        sindex.append(wsindex)
        #print wsdim, sindex

        Ene.append(slines[5+Nbaths:5+Nbaths+wsdim])
        Spin.append(slines[5+Nbaths+wsdim:5+Nbaths+2*wsdim])

        #print Ene
        #print Spin
        count += 1

    count += 1

    FM=[]
    for i in range(nsize):
        bFM=[]
        for b in range(Nbaths):
            (ii, il, s1, s2)  = map(int, lines[count].split()[0:4])
            #print lines[count], len(sdim)
            if (ii!=i+1): print 'something wrong reading cix two'
            if (il!=sindex[i][b]+1): print 'something wrong reading cix three', il, sindex[i][b]
            if il>0:
                if (s1!=sdim[i]): print 'something wrong reading cix four', s1, sdim[i]
                if (s2!=sdim[il-1]): print 'something wrong reading cix five', s2, sdim[il-1]

            clines = lines[count].split()[4:]
            qFM = zeros((s1,s2), dtype=float)
            for i0 in range(s1):
                for i1 in range(s2):
                    qFM[i0,i1] = float(clines[i0*s2 + i1])
            bFM.append(qFM)
                        
            count += 1

        FM.append(bFM)

    # HERE WE NEED TO CHANGE IF OPERATORS ARE PRESENT
    count += 4

    # Reading the HB1 part from the cix file
    (Nm, wnsize, wN_ifl, wmax_size) = map(int, lines[count].split())
    
    count += 2

    wstate=[]
    wstart=[]
    wEne=[]
    wSpin=[]
    wsdim=[]
    wsindex=[]
    wrindex=[]
    for i in range(wnsize):
        
        wstate.append(lines[count])
        slines = lines[count].split()    

        wstart.append(slines[0:6])
        wrindex.append(int(slines[1])-1)
        twsdim = int(slines[5])
        twsindex = map(int, slines[6:6+Nbaths])
        
        wsdim.append(twsdim)

        for iw in range(len(twsindex)): twsindex[iw]-=1
        
        wsindex.append(twsindex)
        
        wEne.append(slines[6+Nbaths:6+Nbaths+twsdim])
        wSpin.append(slines[6+Nbaths+twsdim:6+Nbaths+2*twsdim])

        count += 1

    count += 1

    wFM=[]
    for i in range(wnsize):
        bFM=[]
        for b in range(Nbaths):
            (ii, il, s1, s2)  = map(int, lines[count].split()[0:4])
                        
            if (ii!=i+1): print 'something wrong reading cix two 2'
            if (il!=wsindex[i][b]+1): print 'something wrong reading cix three 2', il, wsindex[i][b]
            if il>0:
                if (s1!=wsdim[i]): print 'something wrong reading cix four 2', s1, wsdim[i]
                if (s2!=wsdim[il-1]): print 'something wrong reading cix five 2', s2, wsdim[il-1]

            clines = lines[count].split()[4:]
            qFM = zeros((s1,s2), dtype=float)
            for i0 in range(s1):
                for i1 in range(s2):
                    qFM[i0,i1] = float(clines[i0*s2 + i1])
            bFM.append(qFM)
            
            count += 1
            
        wFM.append(bFM)
        
    # Creating inverse index between the short and long list of states in cix file
    inv_wrindex=zeros(nsize,dtype=int); inv_wrindex -= 1
    for i in range(wnsize):
        if wrindex[i]>=0:
            inv_wrindex[wrindex[i]]=i

    #print inv_wrindex







    # DEBUGGING HERE
    # Reading probability
    fp = open(fProb, 'r')
    plines = fp.readlines()
    pcount = 1
    Prob = []
    for i in range(nsize):
        inw = inv_wrindex[i]
        #print 'i=', i, sdim[i]
        for j in range(sdim[i]):
            splines = plines[pcount].split()
            (ii, ij) = map(int, splines[:2])
            P = float(splines[2])
            if i+1 != ii : print 'something wrong reading probability 1'
            if j != ij : print 'something wrong reading probability 2'

            Prob.append([P,i,j])
            pcount += 1

    Prob.sort(lambda x, y: cmp(abs(y[0]),abs(x[0])))


    # Creates wlarge for dynamic treatment
    wlarge_small=[]
    for i in range(nsize): wlarge_small.append([])
    for ip,pr in enumerate(Prob):
        if (pr[0]>cutof):
            if len(wlarge_small[pr[1]])<maxDim:
                wlarge_small[pr[1]].append(pr[2])
    
    #print 'wlarge_small', wlarge_small


    # creates small_large and large_small index
    (small_large, large_small) = cmp_small_large(wlarge_small)
    

    #print 'small_large=', small_large
    #for i in range(len(small_large)):
    #    ii = small_large[i][0]
    #    for j in small_large[i][1]:
    #        print 'ii+1=', ii+1, 'j=', j, 'len(Ene[ii])=', len(Ene[ii]), 'True=', j<len(Ene[ii]), '  ', Ene[ii]



    max_size_new = max(map(lambda x: len(x[1]), small_large))
    
    # Start printing !!!
    for l in lines[0:2]: print l,
    print Nm, len(small_large), N_ifl, max_size_new
    for l in lines[3:bath_header]: print l,

    PrintHeader(small_large, large_small, start, sindex, Ene, Spin)
    print "# matrix elements"
    
    for i in range(len(small_large)):
        
        ii = small_large[i][0]
        for b in range(Nbaths):
            ifi = sindex[ii][b]
            
            if ifi>=0 and large_small[ifi]>=0:
                sifi = large_small[ifi]
                ws0 = small_large[i][1]
                ws1 = small_large[sifi][1]
            
                print "%3s %3s %3s %3s  " % (i+1, sifi+1, len(ws0), len(ws1)),
                for i0 in ws0:
                    for i1 in ws1:
                        print FM[ii][b][i0,i1],
                print
            else:
                print "%3s %3s %3s %3s  " % (i+1, 0, 0, 0)
            





    
    
    # Creates wlarge for static treatment
    # First creates static list of states
    dynamic=[]
    for sm in small_large: dynamic.append(inv_wrindex[sm[0]])

    static=dynamic[:]

    for i in dynamic:
        for j in wsindex[i]:
            if (j>=0): static.append(j)

    for i in range(wnsize):
        for j in wsindex[i]:
            if j in dynamic:
                static.append(i)
    
    static = union(static)
    static.sort()


    # Finally creates wlarge index
    wlarge_small=[]
    for i in range(wnsize): wlarge_small.append([])
    for i in range(len(static)):
        ii = static[i]
        wlarge_small[ii] = range(wsdim[ii])

    # creates small_large and large_small index
    (qsmall_large, qlarge_small) = cmp_small_large(wlarge_small)

    
    findex = zeros(len(qsmall_large),dtype=int)
    findex -= 1
    for i in range(len(qsmall_large)):
        #print qsmall_large[i][0] , wrindex[qsmall_large[i][0]]
        findex[i] = large_small[wrindex[qsmall_large[i][0]]]
        
    #print findex


    print "HB1"
    print "# number of operators needed"
    print "0"
    print "# Data for HB1"
    
    max_size_new1 = max(map(lambda x: len(x[1]), qsmall_large))
    
    print Nm, len(qsmall_large), N_ifl, max_size_new1
    print "# ind   N   K   Jz size"
    
    PrintHeader1(qsmall_large, qlarge_small, wstart, wsindex, findex, wEne, wSpin)

    print "# matrix elements"
    
    for i in range(len(qsmall_large)):
        
        ii = qsmall_large[i][0]
        for b in range(Nbaths):
            ifi = wsindex[ii][b]
            
            if ifi>=0 and qlarge_small[ifi]>=0:
                sifi = qlarge_small[ifi]
                ws0 = qsmall_large[i][1]
                ws1 = qsmall_large[sifi][1]
            
                print "%3s %3s %3s %3s  " % (i+1, sifi+1, len(ws0), len(ws1)),
                for i0 in ws0:
                    for i1 in ws1:
                        print wFM[ii][b][i0,i1],
                print
            else:
                print "%3s %3s %3s %3s  " % (i+1, 0, 0, 0)
