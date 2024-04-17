#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import utils,indmffile,sys,re,os,shutil
import optparse #, subprocess
#from scipy import *
from scipy import optimize
import numpy
from numpy import *

nv = list(map(int,numpy.__version__.split('.')))
if (nv[0],nv[1]) < (1,6):
    loadtxt = io.read_array
    #
    def savetxt(filename, data):
        io.write_array(filename, data, precision=16)

def union(data):
    " Takes a union of array or list"
    c = []
    for d in data:
        if d not in c:
            c.append(d)
    return c

def GiveTanMesh(x0,L,Nw):
    def fun(x,x0,L,Nw):
        "x[0]=d, x[1]=w"
        d=x[0]
        w=x[1]
        #print 'd=', d, 'w=', w
        return array([L-w/tan(d), x0-w*tan(pi/(2*Nw)-d/Nw) ])
    
    xi=x0/L
    d0 = Nw/2.*(tan(pi/(2*Nw))-sqrt(tan(pi/(2*Nw))**2 - 4*xi/Nw))
    w0 = L*d0

    sol=optimize.root(fun, [d0,w0], args=(x0,L,Nw) )
    (d,w) = sol.x
    om = w*tan(linspace(0,1,2*Nw+1)*(pi-2*d) -pi/2+d)
    return om
    
#def default_realaxis_mesh():
#    # Probably should allow user to control the mesh generation
#    # to some extent view program options...
#    
#    # use gaumesh.py to generate standard mesh
#    gaumesh = os.path.join( utils.DmftEnvironment().ROOT, 'gaumesh.py' )
#    args = [
#        'x0=[0,0]',
#        'dx0=[0.25,0.01]',
#        'fwhm=[30, 3]',
#        'xmin=-20',
#        'xmax=20'
#        ]
#    stdoutdata, stderrdata = subprocess.Popen([gaumesh] + args, stdout=subprocess.PIPE).communicate()
#    return [float(x) for x in stdoutdata.split('\n') if x.strip() and not x.strip().startswith('#')]


if __name__=='__main__':
    """ Create a zero input self-energy for dmft0, dmft1 and dmft0
    """
    usage = """usage: %prog [ options ]

    Create a zero input self-energy for dmft0, dmft1 and dmft0
    """
    # n==(200/T-1)/2.
    parser = optparse.OptionParser(usage)
    parser.add_option("-e", "--Edc",    dest="Edc",      type="float", default=0.0, help="Starting double counting and Hartree value of the self-energy. Should be close to U(n-1/2)")
    parser.add_option("-i", "--sinp",   dest="insig",   default='sig.inp', help="the mesh will be used from this file.", metavar="FILE")
    parser.add_option("-o", "--sout",   dest="outsig",  default='sig.inp', help="the result will be written to this file")
    parser.add_option("-T", "--Temp",   dest="T",        type="float", default=0.0, help="Temperature")
    parser.add_option("-n", "--nom",    dest="nom",      type="int", default=None, help="Number of frequency points")
    parser.add_option("-L", "--range",  dest="L",        type="float", default=20., help="energy range on real axis")
    parser.add_option("-x", "--x0",     dest="x",        type="float", default=0.05, help="energy range on real axis")
    parser.add_option("-N", "--Nom",    dest="Nom",      type="int", default=400, help="Number of frequency points on real axis")
    parser.add_option("-m", "--magnet", dest="Qmagnetic", action='store_true', default=False, help="should split up/dn self-energy to allows broken symmetry")
    parser.add_option("-B", "--B",      dest="B",        type="float", default=1.0, help="energy shift for magnetic splitting in s_oo")
    
    # Next, parse the arguments
    (options, args) = parser.parse_args()
    
    if options.nom==None and options.T>0:
        options.nom = (300./options.T-1)/(2.*pi)
        
    env = utils.W2kEnvironment()
    case = env.case

    print('Edc=%s case=%s, insig=%s, outsig=%s, nom=%s T=%s' %  (options.Edc, case, options.insig, options.outsig, options.nom, options.T))

    inl = indmffile.Indmfl(case)
    inl.read()
    m_extn = 'dn' if os.path.exists(case+'.indmfl'+'dn') else ''
    if m_extn:
        print('INFO: case.indmfldn present => magnetic calculation with two dmft2 steps')
        inldn = indmffile.Indmfl(case, 'indmfl'+m_extn)
        inldn.read()

    if options.T>0 and inl.matsubara:
        print('..... creating new matsubara mesh of size '+str(options.nom)+' of T='+str(options.T))
        omega = (arange(1,options.nom,1)*2-1)*pi*options.T
    elif os.path.isfile(options.insig) and os.path.getsize(options.insig)>0:
        # Read the input file
        sigdata = loadtxt(options.insig)  # self-energy from 'sig.inp' on long mesh
        if len(shape(sigdata))==1:
            omega = sigdata
        else:
            omega = (sigdata.transpose())[0]
    else:
        if not inl.matsubara:
            print('..... Could not find '+options.insig+'. Generating default real axis mesh.')
            omega = GiveTanMesh(options.x, options.L,int(options.Nom/2))
        else:
            Found=False
            if os.path.isfile('params.dat'):
                exec(compile(open('params.dat', "rb").read(), 'params.dat', 'exec'))
                if 'beta' in iparams0: 
                    Found=True
                    beta = iparams0['beta'][0]
                    if options.nom==None: options.nom = (300*beta-1)/(2.*pi)
                    print('..... creating new matsubara mesh of size '+str(options.nom)+' of T='+str(1./beta))
                    omega = (arange(1,options.nom,1)*2-1)*pi/beta
            if not Found:
                print('..... Could not find '+options.insig+'. Do not know the temperature. Can not create self-energy!')
                print('..... Boiling out.....')
                sys.exit(1)
        
    print('..... Going over all correlated blocks')
    cols=[]
    for icix in inl.siginds:    # over all imp.problems, even those that are equivalent
        Sigind = inl.siginds[icix]
        cols = sort(union(array(Sigind).flatten().tolist()+cols)).tolist()
    if m_extn:
        for icix in inldn.siginds:    # over all imp.problems, even those that are equivalent
            Sigind = inldn.siginds[icix]
            cols = sort(union(array(Sigind).flatten().tolist()+cols)).tolist()
    if 0 in cols: cols.remove(0)

    print('cols=', cols)
    Nc = cols[-1]
    
    if options.Qmagnetic and m_extn:
        cols_up, cols_dn = set(), set()
        for icix in inl.siginds:
            dcols_up = set(  inl.siginds[icix].ravel())-{0}
            dcols_dn = set(inldn.siginds[icix].ravel())-{0}
            if len( dcols_up & cols_dn )==0:
                cols_up.update(dcols_up)
            if len( dcols_dn & cols_up )==0:
                cols_dn.update(dcols_dn)
        cols_up = list(cols_up)
        cols_dn = list(cols_dn)
        print('cols_up=', cols_up, 'cols_dn=', cols_dn)
    # saving the original self-energy if necessary
    if options.insig==options.outsig and os.path.isfile(options.insig) and os.path.getsize(options.insig)>0:
        shutil.copy2(options.insig, options.insig+'.bak')


    #print 'om=', omega

    imp_cols=[]
    findmfi = case+'.indmfi'
    if os.path.exists(findmfi):
        fi = open(findmfi)
        N_imp = int(next(fi).split()[0])
        for i in range(N_imp):
            dim = int(next(fi).split()[0])
            Sigi=[]
            for j in range(dim):
                Sigi.append( list(map(int,next(fi).split()[:dim])) )
            #print 'Sigi=', Sigi
            cols = sort(union(array(Sigi).flatten().tolist())).tolist()
            if 0 in cols: cols.remove(0)
            imp_cols.append(cols)
    else:
        imp_cols=[list(range(1,Nc+1))]
    #print 'imp_cols=', imp_cols
    
    Edc = zeros(Nc).tolist()
    if options.Edc==0 and os.path.isfile('params.dat'):
        exec(compile(open('params.dat', "rb").read(), 'params.dat', 'exec'))
        for imp in range(len(imp_cols)):
            ii = str(imp)
            cols = imp_cols[imp]
            if eval('"U" in iparams'+ii+' and "J" in iparams'+ii+' and "nf0" in iparams'+ii):
                U =  eval('iparams'+ii+'["U"][0]')
                J =  eval('iparams'+ii+'["J"][0]')
                nf = eval('iparams'+ii+'["nf0"][0]')
                edc = U*(nf-0.5)-0.5*J*(nf-1.)
                for ic in cols:
                    Edc[ic-1] = edc
    else:
        Edc = (ones(Nc)*options.Edc).tolist()
    print('Edc=', Edc)
        
    # writing to the file
    fo = open(options.outsig, 'w')
    soo = Edc.copy()

    if options.Qmagnetic:
        for i in cols_up: soo[i-1] -= options.B
        for i in cols_dn: soo[i-1] += options.B
        
    
    print('# s_oo=', soo, file=fo)
    print('#  Edc=', Edc, file=fo)
    for iom,om in enumerate(omega):
        print(("%20.15f "%om), ("0.0 "*(2*Nc)), file=fo)
        
    print(options.outsig, 'written to the disc.')
