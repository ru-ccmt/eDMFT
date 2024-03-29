#!/usr/bin/env python
import yaml
import numpy as np
import pylab as plt
import os, sys

def plot_phonon_bands(params):
    #def plot_phonon_band(idx,path,conv,labl):
    path, plot_total, Emin, Emax, ylab, yconv, savename, figsize, atom_type, atom_name, atom_type_u, natom_u, atom_counts, atom_indx_u, atom_name_u = read_params(params)
    data = yaml.safe_load(open(os.path.join(path, 'band.yaml')))

    # read band energy for each q point and each band
    nqpoint = data['nqpoint']; 
    #print (nqpoint)
    nband = len(data['phonon'][0]['band'] )

    band_ene = np.zeros( (nband, nqpoint), dtype=float)
    qq = np.zeros(nqpoint, dtype=float)

    for iq in range(nqpoint):
        q = data['phonon'][iq]['q-position']
        q_dist = data['phonon'][iq]['distance']
        qq[iq] = q_dist
        for ib in range(nband):
            bands = data['phonon'][iq]['band'][ib]
            freq = bands['frequency']
            band_ene[ib, iq] = freq

    clist=['blue','red','green','orange','cyan']
    lws = [2, 2, 2, 2, 2, 2 ]

    # sort the energy values
    for iq in range(nqpoint):
        ene_q = band_ene[:,iq]
        ene_q_sorted= np.sort(ene_q)
        band_ene[:,iq] = ene_q_sorted

    #plt.figure(1, figsize=(8,4))
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    for ib in range(nband):
        ax.plot( qq, yconv*band_ene[ib], color=clist[0], lw=lws[0])
        if ib==nband-1:
            ax.plot( qq, yconv*band_ene[ib], color=clist[0], lw=lws[0], label=ylab)
    ax.set_xlim(min(qq), max(qq))
    ax.set_axhline(y=0, ls='--', color='k', lw=0.5)

    segs = np.cumsum( data['segment_nqpoint'] )

    #print (segs)

    qticks=[]
    for iseg in segs:
        plt.axvline(x=qq[iseg-1], lw=0.5, color='k', ls='--')
        qticks.append(qq[iseg-1])
    qticks.insert(0, 0)

    # read Q point labels and Segmentations
    if 'labels' in list(data.keys()):
        qlabels = data['labels']
        qlabs=[qlabels[0][0]]
        for qlab in qlabels:
            qlabs.append( qlab[1] )
    else:
        print ('Reading the band path labels from band.conf file.' )
        band_conf = open(os.path.join(path, 'band.conf'), 'r')
        for line in band_conf.readlines():
            if 'BAND_LABELS' in line:
                qlabs = line.split('=')[1].split('#')[0].split()

    for i in range(len(qlabs)):
        if qlabs[i] in ['G', 'Gamma', 'GAMMA', 'Gam', 'gamma']:
            qlabs[i] = '$\Gamma$'

    ax.set_xticks(qticks)
    ax.set_xticklabels(qlabs, fontsize='large', weight = 'bold')
    ax.set_ylabel(ylab, fontsize='large')#, fontweight='bold')
    plt.subplots_adjust(hspace=2, bottom=0.2)
    plt.savefig(savename, dpi=300, bbox_inches='tight')
    plt.show()

def plot_phonon_dos(params):
    path, plot_total_dos, dos_only, band_dos, Emin, Emax, ylab, yconv, savename, figsize, atom_type, atom_name, atom_type_u, natom_u, atom_counts, atom_indx_u, atom_name_u = read_params(params)

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    dos=np.loadtxt(path+'/projected_dos.dat').transpose()

    if plot_total_dos:
        dos_tot=np.loadtxt(path+'/total_dos.dat').transpose()
        ax.plot(         dos_tot[0], 0.5*dos_tot[1], color='lightgray')
        ax.fill_between( dos_tot[0], 0.5*dos_tot[1], color='lightgray')

    for i, iat in enumerate(atom_indx_u):
        clr = [0., 0., 0.]
        clr[i] = 1.
        ax.plot( yconv*dos[0], dos[iat+1],  color=clr, label=atom_name[iat], lw=2.)

    ax.legend(frameon=False, fontsize='large', prop={'weight':'bold'})
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_xlabel(ylab, fontsize='large')#, fontweight='bold')
    ax.set_xlim([Emin, Emax])
    ax.set_ylim(bottom=0)

    plt.subplots_adjust(hspace=2, bottom=0.2)
    plt.savefig(savename, dpi=300, bbox_inches='tight')
    plt.show()

def plot_phonon_bands_character(params):
    path, plot_total_dos, dos_only, band_dos, Emin, Emax, ylab, yconv, savename, figsize, atom_type, atom_name, atom_type_u, natom_u, atom_counts, atom_indx_u, atom_name_u = read_params(params)

    data = yaml.safe_load(open(os.path.join(path, 'band.yaml')))
    # read band energy for each q point and each band
    nqpoint = data['nqpoint']; 
    natom = data['natom']
    nband = len(data['phonon'][0]['band'] )

    print(("nqpoint=%1d, nband=%1d, natom=%1d, natom_unique=%1d" %(nqpoint, nband, natom, natom_u)))

    band_ene = np.zeros( (nband, nqpoint), dtype=float)
    band_char = np.zeros( (nband, nqpoint, natom_u), dtype=float) 


    qq = np.zeros(nqpoint, dtype=float)
    for iq in range(nqpoint):
        q = data['phonon'][iq]['q-position']
        q_dist = data['phonon'][iq]['distance']
        qq[iq] = q_dist
        for ib in range(nband):
            freq  = data['phonon'][iq]['band'][ib]['frequency']
            evecs = data['phonon'][iq]['band'][ib]['eigenvector']
            for iat in range(natom):
                evec_iat = evecs[iat]
                iat_u = atom_type[iat]
                # Sum the like atoms contribution, weighted by the each atoms counts
                band_char[ib, iq, iat_u] += 1./(atom_counts[iat_u] ) * np.sum(np.array(evec_iat)**2)
            band_ene[ib, iq] = freq


    print(('band_ene.shape:', band_ene.shape))
    print(('band_char.shape:', band_char.shape))

    #band_ene.shape  = (nband, nqpoint)
    #band_char.shape = (nband, nqpoint, natom_u)
    # Sort the bands 
    for iq in range(nqpoint):
        indx = np.argsort(band_ene[:, iq])
        band_ene[:, iq] = band_ene[indx, iq]
        band_char[:, iq, :] = band_char[indx, iq, :]

    if band_dos:
        fig, ax = plt.subplots(1, 2, figsize=figsize, gridspec_kw={'width_ratios': [4, 1]})
        ax0 = ax[0]
        ax1 = ax[1]
    else:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax0 = ax
        ax1 = None

    # plot the bands with character
    for ib in range(nband):
        xs = list(zip( qq[:-1], qq[1:] ) )
        ys = list(zip( band_ene[ib, :-1], band_ene[ib, 1:] ) )

        # sum over the band character at a q point and the next q point for the color for a line segment
        chars = [ (band_char[ib, :-1, i] + band_char[ib, 1:, i] )  for i in range(natom_u)]

        rr = chars[0]
        bb = chars[1]
        gg = np.zeros(nqpoint, dtype=float)
        if natom_u>=3:
            gg = chars[2]

        for iq in range(nqpoint-1):
            r = rr[iq]
            b = bb[iq]
            g = gg[iq]

            if r>1.: r = 1.
            if b>1.: b = 1.
            if g>1.: g = 1.

            ax0.plot(xs[iq], yconv*np.array(ys[iq]), ls='-', color=(r, b, g), lw=2.)

    ax0.set_xlim(min(qq), max(qq))
    ax0.axhline(y=0, ls='--', color='k', lw=0.5)

    segs = np.cumsum( data['segment_nqpoint'] )
    qticks=[]
    for iseg in segs:
        ax0.axvline(x=qq[iseg-1], lw=0.5, color='k', ls='--')
        qticks.append(qq[iseg-1])
    qticks.insert(0, 0)
    
    # read Q point labels and Segmentations
    if 'labels' in list(data.keys()):
        labels = data['labels']
        qlabs=[labels[0][0]]
        for lab in labels:
            qlabs.append( lab[1] )
    else:
        print ('Reading the band path labels from band.conf file.')
        band_conf = open(os.path.join(path, 'band.conf'), 'r')
        for line in band_conf.readlines():
            if 'BAND_LABELS' in line:
                qlabs = line.split('=')[1].split('#')[0].split()

    for i in range(len(qlabs)):
        if qlabs[i] in ['G', 'Gamma', 'GAMMA', 'Gam', 'gamma']:
            #qlabs[i] = '$\mathrm{\Gamma}$'
            qlabs[i] = '$\mathbf{\Gamma}$'

    ax0.set_xticks(qticks)
    ax0.set_xticklabels(qlabs, fontsize='large', weight = 'bold')
    ax0.set_ylabel(ylab, fontsize='large')#, fontweight='bold')
    ax0.set_ylim([Emin, Emax])

    # DOS plot
    if ax1:
        dos=np.loadtxt(path+'/projected_dos.dat').transpose()

        if plot_total_dos:
            dos_tot=np.loadtxt(path+'/total_dos.dat').transpose()
            ax1.plot(         0.5*dos_tot[1], dos_tot[0], color='lightgray')
            ax1.fill_between( 0.5*dos_tot[1], dos_tot[0], color='lightgray')
    
        for i, iat in enumerate(atom_indx_u):
            clr = [0., 0., 0.]
            clr[i] = 1.
            ax1.plot( dos[iat+1] , yconv*dos[0], color=clr, label=atom_name[iat], lw=2.)
    
        ax1.legend(frameon=False, fontsize='large', prop={'weight':'bold'})
    
        ax1.set_ylim([Emin, Emax])
    
        ax1.set_xlim(left=Emin)
        ax1.set_xticks([])
        ax1.set_yticklabels([])
        ax1.set_xlabel('DOS', fontsize='large', fontweight='bold')
    
    plt.subplots_adjust(wspace=0.01, bottom=0.15)
    plt.savefig(savename, dpi=300, bbox_inches='tight')
    plt.show()

def check_files(params):
    path = read_params(params)[0]
    files = ['band.yaml', 'projected_dos.dat']
    for f in files:
        if not os.path.exists(path+'/'+f):
            print(path+'/'+f, 'is not present in the given path. Please generate this file.\n\nExitting...')
            sys.exit()

def read_params(params):
    path, plot_total_dos    = params['path'], params['total_dos']
    dos_only, band_dos      = params['dos_only'], params['band_dos']
    Emin, Emax, ylab, yconv = params['Emin'], params['Emax'], params['ylabel'], params['yconv']
    savename, figsize       = params['figname'], params['figsize']
    atom_type, atom_name    = params['atom_type'], params['atom_name']

    # Find the unique atoms and their indices
    atom_type_u = set(atom_type)
    natom_u     = len(atom_type_u)
    atom_counts = [atom_type.count(i) for i in atom_type_u]
    atom_indx_u = [atom_type.index(i) for i in atom_type_u]
    atom_name_u = [atom_name[i]       for i in atom_indx_u]

    return path, plot_total_dos, dos_only, band_dos, Emin, Emax, ylab, yconv, savename, figsize, atom_type, atom_name, atom_type_u, natom_u, atom_counts, atom_indx_u, atom_name_u

if __name__=="__main__":
    #----------Parameters----------#
    p = {'path'    : './FeSeee/',    # directory where the phonopy output files are present
         'total_dos' : False,        # If you want to plot total dos in the DOS plot
         'dos_only'  : False,        # Plot DOS only
         'band_dos'  : False,        # plot bands along with dos in the same figure
         'Emin'      : 0.,           # minimum value of the energy axis
         'Emax'      : 40,           # maximum value of the energy axis
         'ylabel'    :'Energy (meV)',# label for y axis 
         'yconv'     : 1,            # unit conversion for energy axis for meV 33.35643 for cm^-1
         'figname'   : 'fig.pdf',    # figure name to save the figure
         'figsize'   : (8,4),        # figure size to save the plots (8,4) for bands and (8,3) for dos look nice
         'atom_type' : [0, 0, 1, 1], # indices for each atom in the system 
         'atom_name' : ['Fe', 'Fe', 'Se', 'Se'] # name of the atoms in the same order as atom_type
        }
    #-------------------------------#
    
    if  os.path.exists("params.dat"):
        with open("params.dat","r") as fo:
            exec(fo.read(), globals())
        sp = globals()['params']       # stores parameters
        p.update(sp)
    else:
        print ("params.dat was not found. Boiling out...")
        sys.exit(1)
        
    check_files(p)

    if p['dos_only']:
        plot_phonon_dos(p)
    else:
        plot_phonon_bands_character(p)
