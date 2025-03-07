#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule

from scipy import *
import sys
import scipy

from distutils.version import StrictVersion
if StrictVersion(scipy.__version__) > StrictVersion('0.19.0'):
    import weave
else:
    import scipy.weave as weave


def ReadKlist(fklist,nkp):
    #nk = [line[:3]=='END' for line in data].index(True)
    fk = open(fklist,'r')
    data = fk.readlines()
    if data[nkp][:3]!='END': 
        print 'wrong klist ', fklist
    kp=[]
    for i in range(nkp):
        kp.append( map(int, [data[i][10:15], data[i][15:20], data[i][20:25], data[i][25:30]]) )

    #BS = [map(float,line.split()) for line in data[nkp+1:nkp+4]]
    return array(kp) #, array(BS))

def ReadN(fi,size):
    a = next(fi).split()
    while len(a)<size:
        a += next(fi).split()
    return a
    
def ReadBasicArrays():
    fi = open('BasicArrays.dat','r')
    next(fi)
    nat,iso,norbitals,ncix,natom = map(int,next(fi).split())
    next(fi)
    nkpt,nmat,nume = map(int,next(fi).split())
    next(fi)
    next(fi)
    next(fi)
    lmax2,maxdim2,maxdim,maxsize = map(int,next(fi).split())
    next(fi)

    nindo = map(int,ReadN(fi,norbitals))
    #print 'nindo=', nindo
    #nindo = map(int,fi.next().split())
    next(fi)
    cixdim = map(int,ReadN(fi,ncix))
    #print 'cixdim=', cixdim
    #cixdim = map(int,fi.next().split())
    next(fi)
    nl = map(int,ReadN(fi,natom))
    #print 'nl=', nl
    #nl = map(int,fi.next().split())
    next(fi)
    cix_orb=zeros(norbitals,dtype=int)
    ll=[]
    iorbital=[]
    cix=[]
    nind=[]
    for icase in range(natom):
        ll_=[]
        iorbital_=[]
        cix_=[]
        nind_=[]
        for lcase in range(nl[icase]):
            llt = int(next(fi))
            iorb = int(next(fi))
            icix = int(next(fi))
            nindt = int(next(fi))
            ll_.append( llt )
            iorbital_.append( iorb )
            cix_.append( icix )
            nind_.append( nindt )
            cix_orb[iorb-1] = icix-1
            #print llt, iorb, icix, nindt

        ll.append( ll_ )
        iorbital.append( iorbital_ )
        cix.append( cix_ )
        nind.append( nind_ )

    print 'cix_orb=', cix_orb, 'ncix=', ncix
    return (cix_orb, ncix)


if __name__ == '__main__' :
    
    if len(sys.argv)<2:
        print 'ERROR : need input filename. Use the same input as for dmftgk.'
        print '[0] mode           # mode=e/g'
        print '....'
        print '[4]  case.klist    #  klist in line 4'
        print '....'
        print '[10] G_k1          # name of the output gk'
        print '[11] G_local1           # name of the output glocal'
        sys.exit(1)
    
    fin = open(sys.argv[1], 'r')

    mode=next(fin).split()[0]
    for i in range(3): next(fin)
    fklist = next(fin).split()[0]
    for i in range(5): next(fin)
    
    print 'mode=', mode
    print 'fklist=', fklist

    if mode=='g':
        filegk = next(fin).split()[0]
        fileglc = next(fin).split()[0]
        print 'filegk=', filegk
        print 'fileglc=', fileglc
    elif mode=='e':
        next(fin)
        fileUR = next(fin).split()[0]
        fileUL = next(fin).split()[0]
        print 'fileUL=', fileUR
        print 'fileUR=', fileUL
        
    if not (mode=='e' or mode=='g'):
        print 'mode should be one of e or g. Currently it is ', mode
    

    if mode=='g':
        
        # Reading some basic information from G_k
        fg = open(filegk,'r')
        first_line = next(fg)
        nkp,nsymop,nom,cixdm,norbitals = map(int,first_line.split()[1:6])  # dimensions
        second_line = next(fg)
        R_a = array(map(float,second_line.split()[1:1+3*norbitals]))  # actual positions of atoms
        R_a = R_a.reshape(norbitals,3)
        
        cdm1 = cixdm/len(R_a) # dimension of the unfolded block
        
        # correct the first line for upfolded case
        dfirst_line = first_line.split()
        dfirst_line[4] = str(cdm1)
        dfirst_line[5] = str(1)
        first_line = reduce(lambda x,y: x+' '+y, dfirst_line)
        
        # get k-points from klist
        kps = ReadKlist(fklist,nkp)
        
        
        fn = open(filegk+'_', 'w')
        print >> fn, first_line
        print >> fn, second_line.strip()
        
        
        for ik in range(nkp):
            k = kps[ik][:3]/float(kps[ik][3])
            for isym in range(nsymop):
                for iom in range(nom):
                    gg = array(map(float,next(fg).split()))
                    om = gg[0]
                    gkc = (gg[1::2]+gg[2::2]*1j).reshape(cixdm,cixdm)
                    gs = zeros((cdm1,cdm1),dtype=complex)
                    code="""
                    #line 60 "unfold2.py"
                    using namespace std;
                    using namespace blitz;
                    int dim = R_a.extent(0);
                    for (int i=0; i<dim; i++){
                       for (int j=0; j<dim; j++){
                          double phase  = (R_a(j,0)-R_a(i,0))*k(0)+(R_a(j,1)-R_a(i,1))*k(1)+(R_a(j,2)-R_a(i,2))*k(2);
                          complex<double> expi(cos(2*M_PI*phase),sin(2*M_PI*phase));
                          gs += gkc( Range(i*cdm1,(i+1)*cdm1), Range(j*cdm1,(j+1)*cdm1) ) * (expi/static_cast<double>(dim));
                       }
                    }
                    """
                    weave.inline(code, ['R_a', 'k', 'gs', 'gkc', 'cdm1'],type_converters=weave.converters.blitz, compiler='gcc')
                    line = ["%14.8f "%om] + ["%14.8f %14.8f "%(g.real, g.imag) for g in gs.ravel()]
                    print >> fn, reduce(lambda x,y: x+' '+y, line)
            print 'Finished k-point', ik
        fg.close()
        fn.close()
        
        fi = open(fileglc, 'r')
        fn = open(fileglc+'_', 'w')
        
        for iom in range(nom):
            gg = array(map(float,next(fi).split()))
            om = gg[0]
            glc = (gg[1::2]+gg[2::2]*1j).reshape(cixdm,cixdm)
        
            gs = zeros((cdm1,cdm1),dtype=complex)
            for i in range(len(R_a)):
                gs += glc[i*cdm1:(i+1)*cdm1, i*cdm1:(i+1)*cdm1]
            gs *= 1./len(R_a)
            
            line = ["%14.8f "%om] + ["%14.8f %14.8f "%(g.real, g.imag) for g in gs.ravel()]
            print >> fn, reduce(lambda x,y: x+' '+y, line)
        fi.close()
        fn.close()
    elif mode=='e':
        #import utils, indmffile
        ## UL and UR contain all orbitals and atoms, not just for single cix. We need to find which atoms should be unfolded.
        ## First we need to find what contains single icix
        #w2k = utils.W2kEnvironment()    # W2k filenames and paths
        ## Processing 'case.indmfl' file
        #inl = indmffile.Indmfl(w2k.case) # case.indmfl file
        #inl.read()                       # case.indmfl read
        #cixs = inl.siginds.keys()        # all columns for cix
        #atom={}   # atom[icix] contains atoms which need to be considered for unfolding.
        #for icix in cixs:
        #    atoms = [alq[0] for alq in inl.cps[icix]]
        #    atom[icix]=atoms
        #ax=ravel(atom)
        #???sorted(range(len(ax)),key=lambda i:ax[i])
        
        (cix_orb, ncix)=ReadBasicArrays()
        
        # Reading some basic information from U*AR and AL*U
        fr = open(fileUR,'r')
        fl = open(fileUL,'r')

        first_line = next(fr)
        dat = first_line.split()
        nkp,nsymop,nom,norbitals = map(int,dat[1:5])  # dimensions
        dims = map(int,dat[5:5+norbitals])
        ind_orb=[0]
        for idm in dims:
            ind_orb.append(ind_orb[-1]+idm)
        ind_orb=array(ind_orb)
        print 'cix_orb=', cix_orb, 'dims=', dims
        print 'ind_orb=', ind_orb
        
        second_line = next(fr)
        R_a = array(map(float,second_line.split()[1:1+3*norbitals]))  # actual positions of atoms
        R_a = R_a.reshape(norbitals,3)

        cdims=zeros(ncix,dtype=int)
        for icix in range(ncix):
            for iorb in range(norbitals):
                if cix_orb[iorb]==icix: break
            cdims[icix]=dims[iorb]

        cdm1 = dims[0]
        print 'cdm1=', cdm1
        
        next(fl)
        next(fl)
        
        
        #cdm1 = cixdm/len(R_a) # dimension of the unfolded block
        
        # correct the first line for upfolded case
        dfirst_line = first_line.split()
        dfirst_line[4] = str( ncix )
        for i in range(ncix):
            dfirst_line[5+i] = str(cdims[i])
        for i in range(5+ncix,5+norbitals):
            dfirst_line[i]=''
        first_line = reduce(lambda x,y: x+' '+y, dfirst_line)
        
        # get k-points from klist
        kps = ReadKlist(fklist,nkp)
        
        
        frn = open(fileUR+'_', 'w')
        fln = open(fileUL+'_', 'w')
        print >> frn, first_line
        print >> frn, second_line.strip()
        print >> fln, first_line
        print >> fln, second_line.strip()
        
        for ik in range(nkp):
            k = kps[ik][:3]/float(kps[ik][3])
            for isym in range(nsymop):
                for iom in range(nom):
                    line1 = next(fr)
                    dat = line1.split()
                    omega = float(dat[0])
                    nbands = int(dat[1])
                    print >> frn, line1,
                    print >> fln, line1,
                    next(fl)
                    for ibnd in range(nbands):
                        dat = next(fr).split()
                        ii = int(dat[0])
                        gr = array(map(float,dat[1:]))
                        URX = gr[::2]+gr[1::2]*1j

                        dat = next(fl).split()
                        ii = int(dat[0])
                        gl = array(map(float,dat[1:]))
                        ULX = gl[::2]+gl[1::2]*1j
                        
                        #if len(URX)!=len(R_a)*cdm1 or len(ULX)!=len(R_a)*cdm1:
                        if len(URX)!=sum(dims) or len(ULX)!=sum(dims):
                            print 'ERROR: dimension of URX or ULX are not correct', len(URX), len(ULX), sum(dims)
                        
                        liner = ["%5d "%ii]
                        linel = ["%5d "%ii] 
                        for icix in range(ncix):
                            cdim1 = cdims[icix]  # all correlated atoms have to have the same dimension for this to work!!!!
                            URY = zeros(cdim1,dtype=complex)
                            ULY = zeros(cdim1,dtype=complex)
                            dim=norbitals/ncix
                            
                            code="""
                            #line 182 "unfold2.py"
                            using namespace std;
                            using namespace blitz;
                            double nn = 1/sqrt(dim);
                            for (int iorb=0; iorb<norbitals; iorb++){
                               if (cix_orb(iorb)==icix){
                                  double phase  = R_a(iorb,0)*k(0)+R_a(iorb,1)*k(1)+R_a(iorb,2)*k(2);
                                  complex<double> expi(cos(2*M_PI*phase),sin(2*M_PI*phase));
                                  //URY += URX( Range(iorb*cdm1,(iorb+1)*cdm1) ) * (conj(expi) * nn);
                                  //ULY += ULX( Range(iorb*cdm1,(iorb+1)*cdm1) ) * (expi * nn);
                                  URY += URX( Range(ind_orb(iorb),ind_orb(iorb+1)) ) * (conj(expi) * nn);
                                  ULY += ULX( Range(ind_orb(iorb),ind_orb(iorb+1)) ) * (expi * nn);
                               }
                            }
                            """
                            weave.inline(code, ['R_a', 'k', 'URX', 'ULX', 'URY', 'ULY', 'cdm1', 'dim', 'norbitals','cix_orb','icix','ind_orb'],type_converters=weave.converters.blitz, compiler='gcc')
                            
                            liner += ["%14.8f %14.8f "%(g.real, g.imag) for g in URY]
                            
                            linel += ["%14.8f %14.8f "%(g.real, g.imag) for g in ULY]

                        print >> frn, reduce(lambda x,y: x+' '+y, liner)
                        print >> fln, reduce(lambda x,y: x+' '+y, linel)
                    
                        
            print 'Finished k-point', ik
        fr.close()
        fl.close()
        frn.close()
        fln.close()
