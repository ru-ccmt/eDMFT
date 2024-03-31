from numpy import *
from itertools import chain
import sys

def SimplifySiginds(siginds):
    """Takes dictionary of Sigind's and creates list of non-zero values
    that appear in Siginds, cols[icix]=[1,2,3,...]
    """
    cols={}
    for icix in siginds:
        # should remove all negative indices, because they are treated as static self-energy!!!
        cols[icix] = [x for x in set(siginds[icix].ravel()) if x>0] 
    return cols

def ImpurityLatticeConnection( cols, icols, log ):
    """We have correlated index stored in Sigind[icix] and
    impurity index stored in iSigind[iimp].
    The translation between the two is necessary when
    the impurity would like to know in which reference frame
    is impurity problem defined. Local rotation should hence
    be passed to impurity, and for that we need the translator.
    """
    imp2latt={}
    for i in icols: # all impurity problems
        imp2latt[i]=[]
        for j in cols: # all correlated atoms
            dif = set(icols[i])-set(cols[j]) # checks if the same indices
            atm_consistent_with_imp = len(dif)==0 # true when all indices are equal and set diff is empty
            if atm_consistent_with_imp and len(icols[i])==len(cols[j]) and len(icols[i])>0:
                imp2latt[i].append(j)
        if len(imp2latt)==0:
            print('WARNING: in ImpurityLatticeConnection impurity['+str(i)+'] has no connection with correlated atom',file=log)

    all_atoms_in_impurities = sorted(chain(*imp2latt.values()))
    all_atoms = sorted([j for j in cols])
    if all_atoms_in_impurities != all_atoms:
        print('ERROR: in ImpurityLatticeConnection all_atoms=',str(all_atoms), file=log)
        print(' while atoms found in impurities', all_atoms_in_impurities, 'are different', file=log)
        sys.exit(1)
    return imp2latt

def iSigindDropNegative(iSigind,icols,fh_info):
    """ iSigind starts with column number compatible with Sigind from case.indmfl.
    For the impurity calculation, we need to start from 1. In the past,
    we just shifted the index of iSigind so that the first index was 1.
    But if there are negative columns (needed for static shift of some orbitals), this does not work.
    We recently changed sjoin.py so that columns are just enumerated in their increasing order.
    Here we need to be compatible with such choice.
    """
    icols_inv = {}
    for i,ic in enumerate(icols):
        icols_inv[ic]=i+1
    print('icols_inv=', icols_inv, 'icols=', icols, file=fh_info)
    iSigind_new = zeros(shape(iSigind), dtype=int)
    for i in range(shape(iSigind)[0]):
        for j in range(shape(iSigind)[1]):
            if iSigind[i,j]>0:
                iSigind_new[i,j] = icols_inv[ iSigind[i,j] ]
    return iSigind_new

def ConvertFromLatticeToImp(nd, imp2latt):
    nc=[]
    for ic in imp2latt:
        ndc=[ nd[a-1] for a in imp2latt[ic] ]
        nc.append( sum(ndc)/len(ndc) )
    return nc

def Connect_ImpurityProjector(imp2latt,indmfl,struct,log):
    """Find connection between the i-th impurity and j-th projector from projector.dat
    returns list impurity_projector[imp]=icase
    Here imp is unique impurity problem (not cix, which enumates all impurity problems)
    and icase is successive index of projector in projectw.dat
    """
    latt2imp={} # this is inverse index of imp2latt. 
    for imp,cixs in imp2latt.items(): # imp2latt[imp]=[cix_1,cix_2,...]
        for icix in cixs:
            latt2imp[icix]=imp       # latt2imp[cix_i]=imp
    #print('latt2imp=', latt2imp, file=log)
    atms = list(indmfl.atoms.keys())
    impurity_projector=-ones(len(imp2latt),dtype=int)  # the connection between the impurity and projector.dat
    first_atom, icase = 0, 0
    for jatom in range(struct.nat):          # all atom types in structure, just like when building projector
        if (first_atom+1) in atms:           # This atom listed in indmfl file and projector exists
            ia = atms.index(first_atom+1)    # It appears as the ia consequitive atom in indmfl file
            for icix in indmfl.icpsa[ia]:    # all correlated indices for this atom
                if icix>0:                   # icix>0, hence it is correlated
                    imp = latt2imp[icix]   # and icase is the index of projector
                    if impurity_projector[imp]<0:
                        # impurity_projector[imp] should always be negative here
                        # because icix should normally be connected with single
                        # correlated atom with projector.
                        # However, there are exceptions, which are not yet handled
                        # correctly. For example, in Ta2NiSe5 we have a cluster between
                        # two different atoms, Ta-Ni, and in this case one cix is
                        # connected to both Ta and Ni, i.e., two projectors.
                        # We can not assign a single icase to such impurity. We would
                        # need to specify a list of icase in general. But for now we
                        # just take the first one, with understanding that exact-DC
                        # does not work when cluster is between two atoms of different type
                        impurity_projector[imp]=icase # and imp is index of impurity
                    print('ImpurityProjector: Found icix=', icix, 'jatom=', jatom, 'first_atom=', first_atom,
                              'imp=', imp, 'icase=', icase, file=log)
                    icase += 1
        first_atom += struct.mult[jatom]
    return impurity_projector.tolist()

if __name__ == '__main__':
    import utils
    import sys
    import indmffile
    from wstruct import Struct
    
    w2k = utils.W2kEnvironment()    # W2k filenames and paths
    fh_info = sys.stdout
    # Reading 'case.indmfl' file
    inl = indmffile.Indmfl(w2k.case) # case.indmfl file
    inl.read()                       # case.indmfl read

    iSigind = indmffile.ParsIndmfi(w2k.case)
    cols = SimplifySiginds(inl.siginds)
    icols = SimplifySiginds(iSigind)
    imp2latt = ImpurityLatticeConnection( cols, icols, fh_info)
    print('Impurity-lattice connection: imp2latt=', imp2latt, file=fh_info)

    struct = Struct(inAngs=False)
    struct.ReadStruct(w2k.case+'.struct', fh_info)

    Connect_ImpurityProjector(imp2latt,inl,struct,fh_info)
