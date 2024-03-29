#!/usr/bin/env python
import yaml
import numpy as np
import pylab as plt
import optparse
import os

def plot_phonon_band(idx,path,conv,labl):
    data = yaml.load(open(os.path.join(path, 'band.yaml')))

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

    plt.figure(1, figsize=(8,4))
    for ib in range(nband):
        plt.plot( qq, conv*band_ene[ib], color=clist[idx], lw=lws[idx])
        if ib==nband-1:
            plt.plot( qq, conv*band_ene[ib], color=clist[idx], lw=lws[idx], label=labl)
    plt.xlim(min(qq), max(qq))
    plt.axhline(y=0, ls='--', color='k', lw=0.5)

    segs = np.cumsum( data['segment_nqpoint'] )

    #print (segs)

    qticks=[]
    for iseg in segs:
        plt.axvline(x=qq[iseg-1], lw=0.5, color='k', ls='--')
        qticks.append(qq[iseg-1])
    qticks.insert(0, 0)

    # read Q point labels and Segmentations
    if 'labels' in list(data.keys()):
        labels = data['labels']
        labs=[labels[0][0]]
        for lab in labels:
            labs.append( lab[1] )
    else:
        print('Reading the band path labels from band.conf file.')
        band_conf = open(os.path.join(path, 'band.conf'), 'r')
        for line in band_conf.readlines():
            if 'BAND_LABELS' in line:
                labs = line.split('=')[1].split('#')[0].split()

    for i in range(len(labs)):
        if labs[i] in ['G', 'Gamma', 'GAMMA', 'Gam', 'gamma']:
            labs[i] = '$\Gamma$'

    plt.xticks(qticks, labs, fontsize='large')

def plot_phonon_dos(ip,path,conv,labl):
    dos=np.loadtxt(os.path.join(path,'total_dos.dat')).transpose()
    #pdos=np.loadtxt(os.path.join(path,'projected_dos.dat')).transpose()
    #pdos_fil = 'projected_dos.dat'
    #dos_fil  = 'total_dos.dat'
    #pdos = np.loadtxt(pdos_fil).transpose()
    #dos  = np.loadtxt( dos_fil).transpose()
    
    plt.plot(dos[0], dos[1], label=labl, lw=2.)
    plt.legend(loc='upper left')
    plt.ylim(bottom=0.)
    plt.xlim(left=0.)
    plt.grid(True, axis='x')

if __name__=="__main__":
    usage="""
    Plots phonon calculations.
    """
    parser = optparse.OptionParser(usage)
    parser.add_option("-x", dest="xlim",  type="str",  help="x range")
    parser.add_option("-y", dest="ylim",  type="str",  help="y range")
    parser.add_option("-p", dest="plotting",  type="str", default='band', help="what to plot? options: dos, band, pdos")
    parser.add_option("-u", dest="unit",  type="str", default='mev',help="unit for plotting. cm, mev, thz")
    parser.add_option("-l", dest="lab",  type="str", help="labels for paths.")
    parser.add_option("-d", dest="datdir",  type="str", help="data file")
    parser.add_option("-s", dest="savename",  type="str", default='fig.pdf', help="name to save the file.")
    parser.add_option("-t", dest="title",  type="str", default='title', help="title of the plot")

    (options, args) = parser.parse_args()
    
    xl=None;yl=None;labels=None;conv=1.0; ylab='THz'
    if options.xlim:
        xl=options.xlim.split(':')
        if xl[0]=='':    xl=[0,float(xl[1])]
        else:  xl=list(map(float, xl))
    if options.ylim:
        yl=options.ylim.split(':')
        if yl[0]=='':    yl=[0,float(yl[1])]
        else:  yl=list(map(float, yl))

    if options.lab: labels=options.lab.split(',')
    if options.unit=='thz':     (conv, ylab)=(1., 'ThZ')
    elif  options.unit=='mev':     (conv, ylab)=(1, 'meV')#(4.13567, 'meV')
    elif options.unit=='cm':    (conv, ylab)=(33.35641, 'cm^-1')
    savename=options.savename

    if len(args)==0: args=['./']
    toplot=options.plotting

    for ip,path in enumerate(args):
        if labels is None: labl='phonon_dir'+str(ip)
        else:              labl=labels[ip]
        if toplot=='band':
            plot_phonon_band(ip,path,conv,labl)
        elif toplot=='dos':
            plot_phonon_dos(ip,path,conv,labl)
        else:
            print("Not implemented!")

    if xl: plt.xlim([xl[0],xl[1]])
    if yl: plt.ylim([yl[0],yl[1]])
    plt.legend(loc='best',frameon=False)
    plt.ylabel('Energy ['+ylab+']')

    if options.title !='title':
        plt.title(options.title, fontsize=16)

    plt.savefig(savename)
    plt.show()
