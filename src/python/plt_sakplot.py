#!/usr/bin/env python
# @Copyright 2024 Kristjan Haule
import numpy as np
from pylab import *
import glob, os, sys
import scipy
import argparse
import operator
from functools import reduce

def expand_intlist(input):
    '''Expand out any ranges in user input of integer lists.
    Example: input  = "1,2,4-6"
             output = [1, 2, 4, 5, 6]'''
    def parse1(x):
        y = x.split('-')
        return [int(x)] if len(y) == 1 else list(range(int(y[0]), int(y[1])+1))
    return reduce(operator.add, [parse1(x) for x in input.split(',')])

mingle_names={'W':'$W$','L':'$L$','LAMBDA':'$\Lambda$','GAMMA':'$\Gamma$','DELTA':'$\Delta$','X':'$X$','Z':'$Z$','W':'$W$','K':'$K$',
              'R': '$R$', 'S':'$S$', 'T':'$T$', 'U':'$U$', 'Y':'$Y$'}

def ReadDataFile(fname, nkp):
    large_number = 1e6
    try:
        ekw = np.loadtxt(fname)
    except ValueError:
        # Fallback: parse line by line
        data_list = []
        with open(fname, "r") as f:
            for line in f:
                # Strip whitespace and skip empty or comment lines if needed
                line = line.strip()
                if not line or line.startswith(('#','!')):
                    continue
                # Split and convert each field to float
                data_list.append(np.fromstring(line, sep=' '))
        max_columns = np.max([len(line) for line in data_list])
        #print('Fallback max_columns=', max_columns, [len(line) for line in data_list])
        # some data has less bands, but we are not allowed to cut them. We rather increase the bands where they are missing
        data_list = [hstack( (line, ones(max_columns-len(line))*large_number) ) for line in data_list]
        ekw = np.array(data_list)
    return ekw
    
if __name__ == '__main__':
    usage = 'Plots the spectral furnction after dmftk.py or dmftk step has been executed'
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('fname', nargs='?', default='eigenvalues.dat', type=str, help='filename eigvals.dat or eigenvalues.dat. Default=eigvals.dat')
    parser.add_argument('-i', type=float, default=0.97, help='color intensity, a number slightly smaller than 1.(default 0.97). Percentage of points being used to find maximum value.')
    parser.add_argument('-b', type=float, default=1e-5, help='small broadening in calculating A(k,w), default 1e-5')
    parser.add_argument('-d', type=float, default=0, help='shift of zero in y axis when there is a gap and fermi level can be moved away from zero')
    parser.add_argument('-c', type=str, default='cm.coolwarm', help='color map, default is cm.coolwarm but could be changed to cm.hot, cm.Purples or any other matplotlib color map')
    parser.add_argument('-l', type=str, default='k', help='color of the lines.default w for white')
    parser.add_argument('-g', default=False, action='store_true', help='add color bar')
    parser.add_argument('-o', type=str, default=None, help='which orbitals to use for color, example "[\'41,43\',\'42,44\',\'1-40\']" for [R,G,B]')
    parser.add_argument('-a', type=float, default=0.8, help='aspect ratio of the plot, default 0.8')
    parser.add_argument('-f', default=True, action='store_false', help='ignore coherence factors')
    parser.add_argument('-r', type=str, default=None, help='color_intensity default=[1,1,1,1] for [R,G,B,A]. Can be increased/decreased.')
    args = parser.parse_args()
    
    fname = args.fname
    intensity = args.i
    small = args.b
    DY = args.d
    _cmap_ = eval(args.c)
    _col_ = args.l
    color_intensity = [1,1,1,1]
    if args.r is not None:
        color_intensity = eval(args.r)
    
    print('-r color_intensity=', color_intensity)
        
    with open('EF.dat', 'r') as fEF:
        mu = float(fEF.read())
    
    print('using fname=', fname, 'and chemical potential', mu)

    wg = glob.glob('*.klist_band')
    if len(wg)>0:
        fg = open(wg[0], 'r')
        wkpointi=[]
        wkpoints=[]
        for il,line in enumerate(fg):
            if line[:3]=='END': break
            com = line[:10].split()
            if com:
                legnd=line.split()[0]
                if legnd in mingle_names:
                    legnd = mingle_names[legnd]
                else:
                    legnd = '$'+legnd+'$'
                wkpoints.append(legnd)
                wkpointi.append(il)
        print('high-symmetry points found', [(wkpointi[i],wkpoints[i]) for i in range(len(wkpoints))] )
    
    nkp = wkpointi[-1]+1
    print('nkp=', nkp, 'args.f=', args.f)

    if os.path.isfile('A_kom.npy') and os.path.getsize('A_kom.npy')>0:
        COHERENCE=args.f
        import plt_auxiliary as au
        case = au.get_case()
        inl = au.Indmfl(case)
        inl.read()
        iorbs=[]
        ii=0
        for icix in inl.siginds:
            sigind = inl.siginds[icix]
            a = sigind.ravel()
            positives = a[a>0]
            imin = positives.min()
            unique_positive = np.unique(positives)
            inl.legends[icix]
            iatm = inl.cix[icix][0][0]
            #atm = atom_names[iatm-1]
            for it in unique_positive:
                if it!=0:
                    print(f'ii={ii+1:-3d} iatm={iatm:-3d} icix={icix:-3d} it={it:-3d} {inl.legends[icix][it-imin]}')
                    iorbs.append([icix,it,iatm,inl.legends[icix][it-imin]])
                    ii += 1
        
        A_kom = load('A_kom.npy')   # [numk,nom,nii+1]
        om = load('omega.npy')
        (nkp,nom,nii) = shape(A_kom)
        nii-=1
        print('nkp=', nkp, 'nom=', nom, 'nii=', nii, 'shape(A_kom)=', shape(A_kom))
        #print('(icix,it,iatm,leg)=', iorbs)
        if COHERENCE:
            if args.o==None:
                indices=[]
                print("Coherence factors to color. Given a list of "+str(nii)+" entries from Sigind give RGB for all orbitals:")
                for ic,col in enumerate(["red  ","green","blue "]):
                    prompt = col+" (ex: 1,2 3,4 5-8): "
                    userin = input(prompt).lower().strip()
                    expanded = [expand_intlist(group) for group in userin.split()]
                    expandedflat = reduce(operator.add, expanded)
                    indices.append(expandedflat)
            else:
                tcs = eval(args.o)
                if len(tcs)!=3:
                    print('len(args.o) should be exactly 3. It is', len(tcs))
                    sys.exit(0)
                cols=["red  ", "green", "blue "]
                indices=[]
                for ic,tc in enumerate(tcs):
                    expanded = [expand_intlist(group) for group in tc.split()]
                    expandedflat = reduce(operator.add, expanded)
                    indices.append(expandedflat)
                    print(cols[ic], ':', expandedflat)
                #orb_plot=eval(args.o)
                #orb_plot = array(orb_plot)
                #print('orb_plot=', orb_plot, 'and color_intensity=', color_intensity)
                #indices = [np.where(orb_plot == 0)[0], np.where(orb_plot == 1)[0], np.where(orb_plot == 2)[0]]
                #for ind in indices:
                #    if len(ind): ind+=1
                
            print('indices=', indices)
            Akom = zeros((nkp,nom,4))
            A_for_h = A_kom[:,:,0]
            for j in range(3):
                Akom[:,:,j] = sum(A_kom[:,:,indices[j]],axis=2)
            Akom[:,:,3] = A_kom[:,:,0]
            print('max=', [np.max(Akom[:,:,i]) for i in range(4)])
        else:
            Akom = A_kom[:,:,0]
            A_for_h = Akom
    else:
        #print('cmap=', args.c, 'color=', args.l, 'intensity=', args.i, 'small=', args.b, 'DY=', args.d)
        COHERENCE=False
        if args.f and os.path.isfile('UL.dat') and os.path.getsize('UL.dat')>0 and os.path.isfile('UR.dat') and os.path.getsize('UR.dat')>0:
            COHERENCE=True
            print('Using coherence factors')
            fL = open('UL.dat', 'r')
            fR = open('UR.dat', 'r')
            dat = fL.readline().split()
            fR.readline()
            nkp2,nsymop,nom,norbitals = list(map(int,dat[1:5]))  # dimensions
            dims = [int(x) for x in dat[5:5+norbitals]]
            n_all_orbitals = sum(dims)
            fL.readline()
            fR.readline()
            if nkp!=nkp2:
                print('ERROR: UL.dat does not have the same number of k-points as '+fname)
                sys.exit(1)
            if nsymop>1:
                print('ERROR: More than one group operation. We do not want to symmetrize over group operations when plotting spectra!')
                sys.exit(1)
            if args.o==None:
                prompt = "Coherence factors from UL.dat used. Give a list with "+str(n_all_orbitals)+" entries with 0=R,1=G,2=B,3=None > "
                userin = input(prompt).lower().strip()
                orb_plot=eval(userin)
            else:
                orb_plot=eval(args.o)
            if len(orb_plot)<n_all_orbitals:
                norb_plot = ones(n_all_orbitals,dtype=int)*3
                norb_plot[:len(orb_plot)] = orb_plot[:]
                orb_plot = norb_plot
            if len(orb_plot)>n_all_orbitals:
                orb_plot = orb_plot[:n_all_orbitals]
            for i,x in enumerate(orb_plot):
                if x not in [0,1,2]:
                    orb_plot[i]=4
            print('orb_plot=', orb_plot, 'and color_intensity=', color_intensity)
            
        with open(fname, 'r') as fi:
            first_line = fi.readline()
        if first_line.startswith('#'):
            dmft1 = True
        elif first_line.startswith('!'):
            dmft1 = False
        else:
            print('WARNING: Dont recognize file', fname,'. It should start with # or !')
            sys.exit(0)
            
            
        ekw = ReadDataFile(fname,nkp)
        nom = int(len(ekw)/nkp)
        if nkp*nom != len(ekw):
            print('ERROR data length in', fname, 'seems to be incompatible with case.klist_band. We have nkp=', nkp, 'and num lines=', len(ekw), 'which is incomensurate with number of k-points')
            sys.exit(1)
        om = ekw[:nom,0]
        if dmft1:
            zekw = ekw[:,1::2]+ekw[:,2::2]*1j
        else:
            zekw = ekw[:,3::2]+ekw[:,4::2]*1j
        
        nbnd = shape(zekw)[1]
        zekw = reshape(zekw, (nkp,nom,nbnd))
        
        print('nom=', nom, 'shape(zekw)=', shape(zekw),'=(nkp,nom,nbnd)')
        # ensure causality
        zekw = np.where(zekw.imag < -small, zekw, zekw.real - small * 1j)
        
        if COHERENCE:
            #alpha_set=False
            #if 3 in orb_plot: alpha_set=True
            orb_plot = array(orb_plot)
            Akom = zeros((nkp,nom,4))
            for ik in range(nkp):
                cohi = zeros((nom,nbnd,4),dtype=float)
                cohi[:,:,3] = 1.0  # transparency is always set to 1/(om-ek)
                for iw in range(nom):
                    dat = fL.readline().split()
                    omega2 = float(dat[0])
                    if abs(omega2-om[iw])>1e-5:
                        print('WARNING: Frequency in UL.dat '+str(omega2)+' and in '+fname+' '+str(om[iw])+' are not equal!')
                    nbands2 = int(dat[1])
                    #if nbands!=nbands2:
                    #    print('ERROR: Number of bands in UL.dat and '+fname+' is not consistent')
                    fR.readline()
                    
                    for ibnd in range(nbands2):
                        datL = np.fromstring(fL.readline(), sep=' ')
                        datR = np.fromstring(fR.readline(), sep=' ')
                        UL = datL[1::2]+datL[2::2]*1j
                        UR = datR[1::2]+datR[2::2]*1j
                        coh = abs(UL*UR)
                        if ibnd<nbnd:
                            #for iorb in range(len(coh)):
                            #    cohi[orb_plot[iorb],ibnd] += coh[iorb]
                            np.add.at(cohi, (iw, ibnd, orb_plot), coh)
                # Here we normalize weight for each color, because some colors are used for multiple bands.
                #for icolor in range(3):
                #    how_many = np.count_nonzero(orb_plot == icolor)
                #    if how_many>0:
                #        cohi[:,:,icolor] *= color_intensity[icolor]/how_many
                for icolor in range(4):
                    cohi[:,:,icolor] *= color_intensity[icolor]
                #print('cohi=', cohi)
                Akom[ik,:,:] = np.sum( (cohi[:,:,:]/(om[:,None,None] + mu - zekw[ik,:,:,None])).imag, axis=1)*(-1/pi)
        
                #print(Akom[0,:,:])
                #sys.exit(0)
                
                if ik%10==9:
                    print('reading k-pnt', ik+1)
            A_for_h = Akom[:,:,3]
        else:
            # Akom[nkp,nom]
            Akom = np.sum( (1.0/(om[np.newaxis,:,np.newaxis] + mu - zekw)).imag, axis=2)*(-1/pi)
            A_for_h = Akom

    nc = 4 if COHERENCE else 1
    for ic in range(4):
        # finds how to scale the image, and what is the maximum value for density plot
        A_for_h = Akom[:,:,ic].ravel() if COHERENCE else Akom
        ht, bin_edges = histogram(A_for_h,bins=5000)
        xh = 0.5*(bin_edges[1:]+bin_edges[:-1])
        cums = cumsum(ht)/sum(ht)
        i = searchsorted(cums, intensity)
        print('min,max values['+str(ic)+']=', xh[0],xh[-1])
        print('with intensity=', intensity, 'we determine cutoff at max(Ak[:,:,'+str(ic)+'])=', xh[i])
        if ic<3:
            cc = color_intensity[ic]
        else:
            if len(color_intensity)>3:
                cc = color_intensity[ic]*1.8
            else:
                cc = 1.8
        if COHERENCE:
            Akom[:,:,ic] *= cc/xh[i]
        else:
            Akom *= cc
            break
        
    vmm = [0,xh[i]]

    (ymin,ymax) = (om[0]+DY,om[-1]+DY)
    (xmin,xmax) = (0, nkp-1)

    #plot(xh, ht)
    #show()
    
    #print('xmin,xmax,ymin,ymax=', xmin, xmax, ymin, ymax)
    if COHERENCE:
        Akom = np.transpose(Akom, axes=(1, 0, 2))
        imshow(Akom, interpolation='bilinear', origin='lower', vmin=vmm[0], vmax=vmm[1], extent=[xmin,xmax,ymin,ymax], aspect=(xmax-xmin)*args.a/(ymax-ymin) )
        save('Akom', Akom)
        save('rest', hstack( (vmm,[xmin,xmax,ymin,ymax]) ) )
    else:
        imshow(Akom.T, interpolation='bilinear', cmap=_cmap_, origin='lower', vmin=vmm[0], vmax=vmm[1], extent=[xmin,xmax,ymin,ymax], aspect=(xmax-xmin)*args.a/(ymax-ymin) )
    
    for i in range(len(wkpointi)):
        plot([wkpointi[i],wkpointi[i]], [ymin,ymax], _col_+'-')
        
    plot([xmin,xmax],[0,0], _col_+':')

    dytck=0.005
    Ntck=5
    for j in range(len(wkpointi)-1):
        for ix in range(1,Ntck):
            x = wkpointi[j]+(wkpointi[j+1]-wkpointi[j])*ix/float(Ntck)
            plot([x,x],[-dytck,dytck], _col_+'-')
        
    axis([xmin,xmax,ymin,ymax])
    xticks( wkpointi, wkpoints, fontsize='x-large' )
    if args.g:
        colorbar()
    show()
    
