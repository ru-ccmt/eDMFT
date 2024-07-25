#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import sys,re,os
import optparse
import utils
from indmffile import Indmfl, ParsIndmfi
from numpy import *
#import numpy as np
from functools import reduce
from numpy import zeros,array,loadtxt,savetxt,atleast_1d

def SimplifySiginds(siginds):
    """Takes dictionary of Sigind's and creates list of non-zero values
    that appear in Siginds, cols[icix]=[1,2,3,...]
    """
    cols={}
    for icix in siginds:
        cols[icix] = sorted([x for x in set(siginds[icix].ravel()) if x>0])
        # should remove all negative indices!!!
    return cols

#def SimplifySiginds(siginds):
#    " Takes dictionary of Sigind's and creates list or non-zero columns"
#    def union(data):
#        " Takes a union of array or list"
#        c = []
#        for d in data:
#            if d not in c:
#                c.append(d)
#        return c
#    cols={}
#    for icix in siginds:
#        Sigind = siginds[icix]
#        col = sorted([x for x in union(array(Sigind).flatten()) if x>0])
#        cols[icix] = col
#    return cols


if __name__=='__main__':
    """ Takes the delta files, created by dmft1, and combines them into the input file for the
    impurity solver. There can be more output delta's from dmft1, then there are impurity problem
    to be solved.
    """
    usage = """usage: %prog [ options ]

    The script takes the delta files, created by dmft1, and combines them into
    the input file for the impurity solver. There can be more output delta's
    from dmft1, then there are impurity problem to be solved.

    To give filename, you can use the following expressions:
      - word 'case', which will be replaced by current case, determined by
                    the presence of struct file.
      - '?', which will be replaced by icix.
    """

    parser = optparse.OptionParser(usage)
    parser.add_option("-i", "--delta", dest="indlt",  default='case.dlt?', help="filename of the input delta file. Default: 'case.dlt?'")
    parser.add_option("-o", "--Delta", dest="outdlt", default='imp.?/Delta.imp', help="filename of the output Delta to be used by the impurity solver. Default: 'imp.?/Delta'")
    parser.add_option("-p", "--pDelta", dest="pdelta", default=None, help="If need to print the intermediate file which contains all deltas in one file")
    
    parser.add_option("-e", "--Eimp", dest="inE",  default='case.Eimp?', help="filename of the input Eimp file. Default: 'case.Eimp?'")
    parser.add_option("-u", "--Eout", dest="outE", default='imp.?/Eimp.inp', help="filename of the output Eimp to be used by the impurity solver. Default: 'imp.?/Eimp.inp'")
    parser.add_option("-a", "--EimpAverage", dest="EimpAverage", default=1, help="By default we average over up/dn impurity levels. With EimpAverage=0, we skip that.")
    parser.add_option("-d", "--inDC", dest="inDC", default='Edc.dat', help="filename of the input DC-file (Edc.dat)")
    parser.add_option("-m", "--mix", dest="mix", type=float, default=1.0, help="Mixing parameter for hybridization. Default=1 -- no mixing")
    parser.add_option("-c", "--cfield", dest="cf", default=None, help="If CF needs to be reduced. Default=1 -- no reduction")
    parser.add_option("-l", "--lext", dest="m_extn", default='', help="For magnetic calculation, it can be 'dn'.")

    parser.add_option("-g", "--ginp", dest="ginp",  default='case.gc?', help="filename of the input lattice green's function file. Default: 'case.gc?'")
    parser.add_option("-j", "--gout", dest="gout", default='imp.?/Glatt.imp', help="filename of the output lattice Green's function to be used by the impurity solver. Default: 'imp.?/Glatt.imp'")
    
    # Next, parse the arguments
    (options, args) = parser.parse_args()

    env = utils.W2kEnvironment()
    case = env.case

    options.indlt = re.sub(r'case', case, options.indlt)
    options.outdlt = re.sub(r'case', case, options.outdlt)
    options.inE = re.sub(r'case', case, options.inE)
    options.outE = re.sub(r'case', case, options.outE)
    options.ginp = re.sub(r'case', case, options.ginp)
    options.gout = re.sub(r'case', case, options.gout)
    if options.pdelta is not None:
        options.pdelta = re.sub(r'case', case, options.pdelta)
    
    print('case=%s, indlt=%s, outdlt=%s, inE=%s, outE=%s, m_extn=%s' %  (case, options.indlt, options.outdlt, options.inE, options.outE, options.m_extn))
    
    cf_sets=[]
    if options.cf is not None:
        cf_sets0 = options.cf.split(';')
        for c in cf_sets0:
            wcols = eval(c.split()[0])
            wfact = float(c.split()[1])
            cf_sets.append([wcols, wfact])
    
    inl = Indmfl(case)
    inl.read()
    if options.m_extn:
        inldn = Indmfl(case, 'indmfl'+options.m_extn)
        inldn.read()
    
    iSiginds = ParsIndmfi(case)
    icols = SimplifySiginds(iSiginds)
    print('icols=', icols)
    
    cols = SimplifySiginds(inl.siginds)
    print('cols=', cols)
    if options.m_extn:
        colsdn = SimplifySiginds(inldn.siginds)
        print('colsdn=', colsdn)
    
    allcols_ = list(cols.values())
    if options.m_extn: allcols_ += list(colsdn.values())
    allcols = sorted( reduce(lambda x,y: x+y, allcols_) )
    #allcols = sort(array(allcols).flatten())
    
    print('allcols=', allcols)

    if len(allcols)==0: sys.exit(0)

    print('max(allcols)=', max(allcols))
    noccur = zeros(max(allcols),dtype=int)
    #print 'len(noccur)=', len(noccur)
    for c in allcols:
        noccur[c-1]+=1

    print('noccur=', noccur)
    
    # to check number of frequency points
    filename = re.sub(r'\?', str(1), options.indlt)
    print('filename=', filename)
    data = loadtxt( filename, ndmin=2 ).transpose()
    om = data[0]
    
    # Array of all Delta's, after averaging over all columns and over all impurity problems.
    rDeltas=zeros(( max(allcols)*2+1, len(om) ),dtype=float)
    rDeltas[0] = om
    rGlatt=zeros(( max(allcols)*2+1, len(om) ),dtype=float)
    rGlatt[0] = om
    rEimps = zeros( max(allcols), dtype=float )
    rOlaps = zeros( max(allcols), dtype=float )
    
    # Reading output of dmft1
    for icix in list(cols.keys()):
        if len(cols[icix])==0: continue
        filename = re.sub(r'\?', str(icix), options.indlt)
        data = loadtxt( filename ).transpose()
        mincol = min(cols[icix])
        print(('icix=%d reading from: %s' % (icix, filename)), 'cols=', cols[icix])
        for c in cols[icix]:
            #print '%d Taking from col #(%d,%d) and putting into #(%d,%d)' % (c, 2*(c-mincol+1)-1, 2*(c-mincol+1), 2*c-1, 2*c )
            rDeltas[2*c-1] += data[2*(c-mincol+1)-1]*(1./noccur[c-1])
            rDeltas[2*c]   += data[2*(c-mincol+1)  ]*(1./noccur[c-1])

        filename = re.sub(r'\?', str(icix), options.ginp)
        data = loadtxt( filename ).transpose()
        mincol = min(cols[icix])
        print(('icix=%d reading from: %s' % (icix, filename)), 'cols=', cols[icix])
        for c in cols[icix]:
            #print '%d Taking from col #(%d,%d) and putting into #(%d,%d)' % (c, 2*(c-mincol+1)-1, 2*(c-mincol+1), 2*c-1, 2*c )
            rGlatt[2*c-1] += data[2*(c-mincol+1)-1]*(1./noccur[c-1])
            rGlatt[2*c]   += data[2*(c-mincol+1)  ]*(1./noccur[c-1])
            
        # Processing impurity levels
        Eimpfile = re.sub(r'\?', str(icix), options.inE)
        exec(compile(open(Eimpfile, "rb").read(), Eimpfile, 'exec'))
        for c in cols[icix]:
            #print 'Eimpfile=%s  %d Taking from col #(%d) and putting into #(%d)' % (Eimpfile, c, c-mincol, c-1 )
            rEimps[c-1] += Ed[c-mincol]*(1./noccur[c-1])
            rOlaps[c-1] += Olap[c-mincol]*(1./noccur[c-1])

    if options.m_extn:
        # Reading output of dmft1
        for icix in list(colsdn.keys()):
            if len(colsdn[icix])==0: continue
            filename = re.sub(r'\?', str(icix), options.indlt)
            filename += options.m_extn
            data = loadtxt( filename ).transpose()
            
            mincol = min(colsdn[icix])
            print(('icix=%d reading from: %s' % (icix, filename)), 'cols=', colsdn[icix])
            for c in colsdn[icix]:
                #print '%d Taking from col #(%d,%d) and putting into #(%d,%d)' % (c, 2*(c-mincol+1)-1, 2*(c-mincol+1), 2*c-1, 2*c )
                #print 'c=', c, 'mincol=', mincol, 'noccur=', noccur[c-1], 'ind=', 2*(c-mincol+1)-1
                #print len(data[2*(c-mincol+1)-1]*(1./noccur[c-1]))
                #print len(rDeltas[2*c-1])
                rDeltas[2*c-1] += data[2*(c-mincol+1)-1]*(1./noccur[c-1])
                rDeltas[2*c]   += data[2*(c-mincol+1)  ]*(1./noccur[c-1])

            filename = re.sub(r'\?', str(icix), options.ginp)
            filename += options.m_extn
            data = loadtxt( filename ).transpose()
            
            print(('icix=%d reading from: %s' % (icix, filename)), 'cols=', colsdn[icix])
            for c in colsdn[icix]:
                #print '%d Taking from col #(%d,%d) and putting into #(%d,%d)' % (c, 2*(c-mincol+1)-1, 2*(c-mincol+1), 2*c-1, 2*c )
                rGlatt[2*c-1] += data[2*(c-mincol+1)-1]*(1./noccur[c-1])
                rGlatt[2*c]   += data[2*(c-mincol+1)  ]*(1./noccur[c-1])
            
                
            # Processing impurity levels
            Eimpfile = re.sub(r'\?', str(icix), options.inE)
            Eimpfile += options.m_extn
            exec(compile(open(Eimpfile, "rb").read(), Eimpfile, 'exec'))
            for c in colsdn[icix]:
                #print 'Eimpfile=%s  %d Taking from col #(%d) and putting into #(%d)' % (Eimpfile, c, c-mincol, c-1 )
                rEimps[c-1] += Ed[c-mincol]*(1./noccur[c-1])
                rOlaps[c-1] += Olap[c-mincol]*(1./noccur[c-1])

        print('options.EimpAverage=', options.EimpAverage)
        if int(options.EimpAverage):
            print('Averaging!')
            # Averaging impurity leveles between up and down. Should be the same, except for the numerical error
            for icix in list(colsdn.keys()):
                for ic in range(len(cols[icix])):
                    cu = cols[icix][ic]
                    cd = colsdn[icix][ic]
                    Eaver = 0.5*(rEimps[cu-1]+rEimps[cd-1])
                    Ediff = 0.5*(rEimps[cu-1]-rEimps[cd-1])
                    #print 'averaging over up and down', ic, cu, cd, rEimps[cu-1], rEimps[cd-1], Eaver, Ediff
                    rEimps[cu-1] = Eaver
                    rEimps[cd-1] = Eaver


    
    # Storing intermediate data, if necessary
    if options.pdelta is not None:
        savetxt(options.pdelta, rDeltas.transpose())
    #if options.hlpE is not None:
    #    fE = open(options.hlpE, 'w')
    #    print >> fE, 'Ed=',rEimps.tolist(), '\n', 'Olap=',rOlaps.tolist()


    Edc = atleast_1d(loadtxt(options.inDC))

    # Reducing CF if necessary
    for cf in cf_sets:
        factor = cf[1]
        mean=zeros((2,len(om)), dtype=float)
        Emean=0
        for c in cf[0]:
            mean[0] += rDeltas[2*c-1]
            mean[1] += rDeltas[2*c]
            Emean += rEimps[c-1]
        mean *= 1./len(cf[0])
        Emean *= 1./len(cf[0])
        for c in cf[0]:
            rDeltas[2*c-1] = mean[0] + (rDeltas[2*c-1]-mean[0])*factor
            rDeltas[2*c]   = mean[1] + (rDeltas[2*c]-mean[1])*factor
            rEimps[c-1] = Emean + (rEimps[c-1]-Emean)*factor
            
    
    
    # Arranging the Delta's for the impurity problems.
    for icix in list(icols.keys()):
        # Processing Delta
        mincol = min(icols[icix])
        maxcol = max(icols[icix])
        filename = re.sub(r'\?', str(icix), options.outdlt)
        fileGout = re.sub(r'\?', str(icix), options.gout)
        print(('icix=%d saving into: %s %s' % (icix, filename,fileGout)), 'cols=', icols[icix])

        #data=zeros((2*(maxcol-mincol+1)+1,len(om)), dtype=float)
        data=zeros((2*len(icols[icix])+1, len(om)), dtype=float)
        data[0] = om
        #datG=zeros((2*(maxcol-mincol+1)+1,len(om)), dtype=float)
        datG=zeros((2*len(icols[icix])+1, len(om)), dtype=float)
        datG[0] = om
        for ic,c in enumerate(icols[icix]):
            #data[2*(c-mincol+1)-1] = rDeltas[2*c-1]
            #data[2*(c-mincol+1)]   = rDeltas[2*c]
            #datG[2*(c-mincol+1)-1] = rGlatt[2*c-1]
            #datG[2*(c-mincol+1)]   = rGlatt[2*c]
            #print '%d Taking from col #(%d,%d) and putting into #(%d,%d)' % (c, 2*c-1, 2*c, 2*(c-mincol+1)-1, 2*(c-mincol+1) )
            data[2*ic+1] = rDeltas[2*c-1]
            data[2*ic+2] = rDeltas[2*c]
            datG[2*ic+1] = rGlatt[2*c-1]
            datG[2*ic+2] = rGlatt[2*c]
            print('%d Taking from col #(%d,%d) and putting into #(%d,%d)' % (c, 2*c-1, 2*c, 2*ic+1, 2*ic+2 ))
            
        cdir = os.path.split(filename)[0]
        if cdir and not os.path.exists(cdir):
            print('Output directory '+cdir+' does not exist! Creating it....!')
            os.makedirs(cdir)

        if options.mix!=1.0 and os.path.isfile(filename) and os.path.getsize(filename)>0:
            print('Mixing delta with mix=', options.mix)
            data0 = loadtxt(filename).transpose()
            data = options.mix*data + (1-options.mix)*data0
            
        # Writting Delta
        savetxt(filename, data.transpose())
        # Writting Glatt
        savetxt(fileGout, datG.transpose())
        
        # Processing Impurity levels
        Edat=zeros(len(icols[icix]), dtype=float)
        Odat=zeros(len(icols[icix]), dtype=float)
        Edca=zeros(len(icols[icix]), dtype=float)
        for ic,c in enumerate(icols[icix]):
            Edat[ic] += rEimps[c-1]
            Odat[ic] += rOlaps[c-1]
            Edca[ic] += Edc[c-1]

        # Writting impurity levels to the file
        filename = re.sub(r'\?', str(icix), options.outE)
        
        if options.mix!=1.0 and os.path.isfile(filename) and os.path.getsize(filename)>0:
            print('Mixing impurity levels with mix=', options.mix)
            data0 = loadtxt(filename)
            
            print('shape1', shape(data0[0]))
            print('shape2', shape(Edat))
            
            Edat = options.mix*Edat + (1-options.mix)*data0[0]
            Odat = options.mix*Odat + (1-options.mix)*data0[1]
            Edca = options.mix*Edca + (1-options.mix)*data0[2]
        
        fE = open(filename, 'w')
        print("%15.10f "*len(Edat) % tuple(Edat), ' # Eimp', file=fE)
        print("%15.10f "*len(Odat) % tuple(Odat), ' # Overlap', file=fE)
        print("%15.10f "*len(Edca) % tuple(Edca), ' # Edc', file=fE)
