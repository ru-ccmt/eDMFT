#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import sys,re,os
import optparse
import utils
from indmffile import Indmfl, ParsIndmfi
#from scipy import *
from numpy import array,zeros,loadtxt,savetxt,shape
from functools import reduce

#def union(data):
#    " Takes a union of array or list"
#    c = []
#    for d in data:
#        if d not in c:
#            c.append(d)
#    return c
#
#def SimplifySiginds(siginds):
#    " Takes dictionary of Sigind's and creates list or non-zero columns"
#    colsp={}
#    colsm={}
#    for icix in list(siginds.keys()):
#        Sigind = siginds[icix]
#        cols_all = union(array(Sigind).flatten())
#        cols_all = sorted(cols_all,key=lambda x: abs(x))
#        colp=[x for x in cols_all if x>0]
#        colm=[x for x in cols_all if x<0]
#        #col = sort(union(array(Sigind).flatten())).tolist()
#        #if 0 in col: col.remove(0)
#        colsp[icix] = colp
#        colsm[icix] = colm
#    return (colsp,colsm)

def SimplifySiginds(siginds):
    """Takes dictionary of Sigind's and creates list of non-zero values
    that appear in Siginds, cols[icix]=[1,2,3,...]
    """
    colsp={}
    colsm={}
    for icix in siginds:
        unique_values = set(siginds[icix].ravel())
        colsp[icix] = sorted([x for x in unique_values if x>0])
        colsm[icix] = sorted([x for x in unique_values if x<0], reverse=True)
        # should remove all negative indices!!!
    return (colsp,colsm)

if __name__=='__main__':
    """ Takes the self-energy files from all impurity problems, and combines
    them into a single self-energy file.
    """
    usage = """usage: %prog [ options ]

    The script takes the self-energy files from all impurity problems,
    and combines them into a single self-energy file.

    To give filename, you can use the following expressions:
      - word 'case', which will be replaced by current case, determined by
                    the presence of struct file.
      - '?', which will be replaced by icix.
    """

    parser = optparse.OptionParser(usage)
    parser.add_option("-o", "--osig", dest="osig", default='sig.inp', help="filename of the output self-energy file. Default: 'sig.inp'")
    parser.add_option("-i", "--isig", dest="isig", default='imp.?/sig.out', help="filename of the input self-energy from all impurity problems. Default: 'imp.?/sig.out'")
    parser.add_option("-l", "--lext", dest="m_extn", default='', help="For magnetic calculation, it can be 'dn'.")
    parser.add_option("-m", "--mix", dest="mix", type=float, default=1.0, help="Mixing parameter for self-energy. Default=1 -- no mixing")

    # Next, parse the arguments
    (options, args) = parser.parse_args()
    
    env = utils.W2kEnvironment()
    case = env.case
    
    options.osig = re.sub(r'case', case, options.osig)
    options.isig = re.sub(r'case', case, options.isig)
    
    print('sgather.py: case=%s, isig=%s, osig=%s' %  (case, options.isig, options.osig))
    

    inl = Indmfl(case)
    inl.read()
    if options.m_extn:
        inldn = Indmfl(case, 'indmfl'+options.m_extn)
        inldn.read()

    # Impurity quantities
    iSiginds = ParsIndmfi(case)
    icols,icolsm = SimplifySiginds(iSiginds)
    print('sgather: icols=', icols)

    # lattice quantities
    colsp,colsm = SimplifySiginds(inl.siginds)
    print('sgather: colsp=', colsp, 'colsm=', colsm)

    if options.m_extn:
        colspdn, colsmdn = SimplifySiginds(inldn.siginds)
        print('sgather: colspdn=', colspdn, 'colsmdn=', colsmdn)


    # these are columns we can fill in whit our impurity problems
    allcols = sorted( reduce(lambda x,y: x+y, list(icols.values())) )
    print('sgather: allcols=', allcols)
    #print('len(allcols)=', len(allcols))
    noccur = zeros(max(allcols),dtype=int)
    #print('len(noccur)=', len(noccur))
    # these columns can not be obtained from impurity problems. Need another way.
    missing_columns={}
    for icix in list(colsm.keys()):
        # This needs particular attention. Note that there are two types of negative cix indices :
        # a) If both positive and negative indices appear for the same orbital (same icix) then we expect
        #    to just merely shift the spectra that corresponds to the negative columns, and in this case
        #    we add such terms in s_oo and Edc.
        # b) If all indices of certain orbitals are negative, such orbital is treated as "open core" orbital,
        #    and is not computed by solving an impurity problem. Rather it is split by a fixed self-energy,
        #    stored in the file sfx.[icix]. This case needs particular attention in ssplit.py part, but
        #    not in sgather.py step, as such "open core" orbitals do not appear in sig.inp, and hence can
        #    be safely ignored in this step.
        if len(colsp[icix]) > 0:  # only if some orbitals are positive in this icix, such case is not "open core" case, and we need to treat it.
            for i in colsm[icix]: # missing columns do not appear in self-energy
                missing_columns[abs(i)]=1
    print('sgather: missing_columns=', missing_columns)

    for c in allcols:
        noccur[c-1]+=1
    print('sgather: noccur=', noccur)
    
    filename = re.sub(r'\?', str(0), options.isig)
    om = loadtxt(filename).transpose()[0]
    
    # Array of all sigmas
    rSigma=zeros(( len(allcols)*2+1, len(om) ),dtype=float)
    rs_oo=zeros( len(allcols)+len(missing_columns) ,dtype=float)
    rEdc=zeros( len(allcols)+len(missing_columns) ,dtype=float)
    rSigma[0] = om
    print('sgather: shape(rSigma)', shape(rSigma))
    print('sgather: shape(rs_oo)=', shape(rs_oo))

    # Reading self-energies from all impurity problems
    for icix in list(icols.keys()):
        # Processing Delta
        #mincol = min(icols[icix])
        #maxcol = max(icols[icix])
        
        filename = re.sub(r'\?', str(icix), options.isig)
        print(('sgather: icix=%d reading from: %s' % (icix, filename)), 'cols=', icols[icix])
        
        # Searching for s_oo and Edc
        with open(filename, 'r') as fh_sig:
            line1 = next(fh_sig)[1:]
            line2 = next(fh_sig)[1:]
            print('sgather: line1=', line1,end='')
            print('sgather: line2=', line2,end='')
            exec(line1.strip(), globals())
            exec(line2.strip(), globals())
        # reading sigma
        data = loadtxt(filename).transpose()
        
        for iic,c in enumerate(icols[icix]):
            imiss = len([x for x in missing_columns if x<=c])
            ic = c-imiss
            #print 'writting into ', 2*ic-1, 'and', 2*ic, 'from data at ', 2*(c-mincol+1)-1, 'and', 2*(c-mincol+1)
            print('sgather: writting into ', 2*ic-1, 'and', 2*ic, 'from data at ', 2*iic+1, 'and', 2*iic+2)
            #rSigma[2*ic-1] += data[2*(c-mincol+1)-1]*(1./noccur[c-1])
            #rSigma[2*ic]   += data[2*(c-mincol+1)  ]*(1./noccur[c-1])
            rSigma[2*ic-1] += data[2*iic+1]*(1./noccur[c-1])
            rSigma[2*ic]   += data[2*iic+2]*(1./noccur[c-1])
        
        for iic,c in enumerate(icols[icix]):
            #print 'writting s_oo into ', c-1, 'from', c-mincol
            print('writting s_oo into ', c-1, 'from', iic)
            rs_oo[c-1] += s_oo[iic]*(1./noccur[c-1])
            rEdc[c-1]  += Edc[iic]*(1./noccur[c-1])
        
    print('sgather: rs_oo=', rs_oo, 'rEdc=', rEdc)
    # The impurity self-energy has a static contribution
    # equal to s_oo-Edc
    # This self-energy will however go to zero in infinity. No static contribution
    for c in allcols: 
        imiss = len([x for x in missing_columns if x<=c])
        ic = c-imiss
        rSigma[2*ic-1] -= (rs_oo[c-1]-rEdc[c-1])


    #cm = colsm.values()
    #if options.m_extn:
    #    cm += colsmdn.values()
    #cm=array(union(cm)).flatten()
    #print 'cm=', cm

    if len(missing_columns)>0: # Some columns contain only s_oo and Edc (no dynamic component)
        # old-s_oo and old-Edc
        with open(options.osig,'r') as fi: # checking old s_oo and Edc from sig.inp
            exec(next(fi)[1:].strip(), globals())
            exec(next(fi)[1:].strip(), globals())
            
        for i in list(missing_columns.keys()):
            rs_oo[i-1] = s_oo[i-1]
            rEdc[i-1] = Edc[i-1]

    if options.mix!=1.0 and os.path.isfile(options.osig) and os.path.getsize(options.osig)>0:
        print('Mixing self-energy with mix=', options.mix)
        with open(options.osig,'r') as fi: # checking old s_oo and Edc
            exec(next(fi)[1:].strip(), globals())
            exec(next(fi)[1:].strip(), globals())
            
        rSigma_old = loadtxt(options.osig).transpose()
        
        # The number of frequencies should be the same, but sometimes some frequencies are mising in the new iteration. Than mix just existing frequencies
        nom_new = shape(rSigma)[1]
        print('shape sigma old', shape(rSigma_old))
        print('shape sigma new', shape(rSigma))
        
        rs_oo[:] = options.mix*rs_oo[:]  + (1-options.mix)*array(s_oo)
        rEdc[:]  = options.mix*rEdc[:]   + (1-options.mix)*array(Edc)
        rSigma = options.mix*rSigma[:,:] + (1-options.mix)*rSigma_old[:,:nom_new]

    
    print('osig=', options.osig)
    with open(options.osig, 'w') as fout:
        print('# s_oo=', rs_oo.tolist(), file=fout)
        print('# Edc=', rEdc.tolist(), file=fout)
        savetxt(fout,rSigma.transpose())
