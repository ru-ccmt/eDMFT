#!/usr/bin/env python
from pylab import *
import numpy as np
import re, os
import glob
import argparse
import utils
from indmffile import Indmfl, ParsIndmfi
from imp2lattc import ImpurityLatticeConnection, SimplifySiginds
from wstruct import Struct

class InteractiveLegend:
    def __init__(self,lines,leg,fig):
        self.map_legend_to_ax = {}  # Will map legend lines to original lines.
        pickradius = 5  # Points (Pt). How close the click needs to be to trigger an event.
        for legend_line, ax_line in zip(leg.get_lines(), lines):
            legend_line.set_picker(pickradius)  # Enable picking on the legend line.
            self.map_legend_to_ax[legend_line] = ax_line
        # Works even if the legend is draggable. This is independent from picking legend lines.
        self.leg = leg
        self.leg.set_loc('best')
        self.leg.set_draggable(True)
        self.fig = fig
        
    def __call__(self, event):
        # On the pick event, find the original line corresponding to the legend
        # proxy line, and toggle its visibility.
        legend_line = event.artist
    
        # Do nothing if the source of the event is not a legend line.
        if legend_line not in self.map_legend_to_ax:
            return
        
        ax_line = self.map_legend_to_ax[legend_line]
        visible = not ax_line.get_visible()
        ax_line.set_visible(visible)
        # Change the alpha on the line in the legend, so we can see what lines
        # have been toggled.
        legend_line.set_alpha(1.0 if visible else 0.2)
        self.fig.canvas.draw()

if __name__=='__main__':
    usage = 'Plots Green\'s function or self-energy by reading imp.?/Gf.out.* and case.gc? files'
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('fname', nargs='?', default='gc', type=str, help='filename either gc or sig.inp. Default=gc')
    parser.add_argument('-x', type=str, help='x range of the plot in the form -x0:10')
    parser.add_argument('-y', type=str, help='y range of the plot in the form -y0:10')
    parser.add_argument('-l', type=str, default=None, help='lattice lcols to show in the form of dic like {0:[1,2,3], 1:[1,2,3]}')
    parser.add_argument('-n', type=str, default=None, help='lattice lcols_dn to show in the form of dic like {0:[1,2,3], 1:[1,2,3]}')
    parser.add_argument('-g', default=False, action='store_true', help='grid')
    parser.add_argument('-a', default='1', type=str, help='how many impurity steps displayed. Default one.')
    args = parser.parse_args()

    xrng = [None,None]
    yrng = [None,None]
    if args.x!=None:
        w = args.x.split(':')
        xrng = [float(w[i]) if w[i]!='' else None for i in range(2)]
    if args.y!=None:
        w = args.y.split(':')
        yrng = [float(w[i]) if w[i]!='' else None for i in range(2)]
    lastn = int(args.a)

    env = utils.W2kEnvironment()

    dat = loadtxt(env.case+'.cdos').T
    with open(env.case+'.cdos') as fi:
        line = fi.readline()
        leg = re.findall(r"total|a=\s*\d+\s+L=\s*\d+", line) # legends parsing
        #print(leg)


    lines=[]
    # printing case.dos
    fig, axs = plt.subplots(1,1)
    om = dat[0]
    DOS = dat[1:,:]
    il=0
    for l in range(len(DOS)):
        linei,= axs.plot(om,DOS[l,:], 'C'+str(il%10), label=leg[l])
        lines.append(linei)
        il+=1
    
    inl = Indmfl(env.case)
    inl.read()

    # below is Green's function plotting, to extract partial dos in MT-spheres
    struct = Struct()
    struct.ReadStruct(env.case+'.struct')
    atom_names0 = [[struct.aname[iat]]*struct.mult[iat] for iat in range(struct.nat)]
    atom_names = [item for sublist in atom_names0 for item in sublist]
    print('atom_names=', atom_names)
    
    iSiginds = ParsIndmfi(env.case)
    _icols = SimplifySiginds(iSiginds)
    _lcols = SimplifySiginds(inl.siginds)
    _lcols_dn=None
    if os.path.isfile(env.case+'.indmfldn'):
        inldn = Indmfl(env.case, 'indmfldn')
        inldn.read()
        _lcols_dn = SimplifySiginds(inl.siginds)

    so_present = (os.path.isfile(env.case+'.inso') and os.path.getsize(env.case+'.inso')>0)
    #nspin=2
    nspin=1
    if _lcols_dn is not None or so_present:
        nspin=1
        
    imp2latt = ImpurityLatticeConnection(_lcols, _icols, sys.stdout)
    # these are columns we can fill in with our impurity problems
    print('_icols=', _icols)
    print('_lcols=',_lcols)
    # First create simple columns for communicating with the user.
    icols={}
    for ii in _icols:
        icols[ii+1] = list(range(len(_icols[ii])))
    lcols={}
    degs={}
    for icix in _lcols:
        lcols[icix] = (array(_lcols[icix])-min(_lcols[icix])).tolist()
        # degeneracies of each column
        degs[icix]  = [np.count_nonzero(inl.siginds[icix] == ic) for ic in _lcols[icix]]
    lcols_dn=None
    degsdn={}
    if _lcols_dn is not None:
        lcols_dn={}
        for icix in _lcols_dn:
            lcols_dn[icix] = (array(_lcols_dn[icix])-min(_lcols_dn[icix])).tolist()
            degsdn[icix]  = [np.count_nonzero(inldn.siginds[icix] == ic) for ic in _lcols_dn[icix]]
    
    print('imp2lattice connection=', imp2latt)
    for ii in imp2latt:
        print(' imp.'+str(ii)+'/Gf.out icols['+str(ii+1)+']=', icols[ii+1])
        print(' connected to ')
        for icix in imp2latt[ii]:
            print('  gc'+str(icix)+'   lcols['+str(icix)+']=',array(lcols[icix]))
            if lcols_dn is not None:
                print('  gc'+str(icix)+'dn lcols['+str(icix)+']=',array(lcols_dn[icix]))
                
        
    print('All available columns:')
    print(' icols=', icols)
    print(' lcols=', lcols)
    if lcols_dn is not None:
        print(' lcols_dn=', lcols_dn)
    changed=False
    if args.l is not None:
        lcols = eval(args.l)
        changed = True
    if lcols_dn is not None and args.n is not None:
        lcols_dn = eval(args.n)
        changed = True
    if changed:
        print('user changed to')
        print(' icols=', icols)
        print(' lcols=', lcols)
        if lcols_dn is not None:
            print(' lcols_dn=', lcols_dn)

    for ii in imp2latt:
        il=0
        foundOneYet,foundOneYetdn=False,False
        for icix in imp2latt[ii]:
            if icix in lcols and not foundOneYet:
                fname = env.case+'.gc'+str(icix)
                if not os.path.exists(fname):
                    continue
                else:
                    foundOneYet=True
                iatm = inl.cix[icix][0][0]
                atm = atom_names[iatm-1]
                print('atm=', atm)
                cls = array(lcols[icix])
                print('reading', fname, 'cols=', cls)
                dat = loadtxt(fname).T
                w = dat[0]
                dos = array([-dat[2+2*j]*degs[icix][j]*nspin*len(imp2latt[ii])/pi for j in cls])
                for i,j in enumerate(cls):
                    linei,= axs.plot(w,dos[i], 'C'+str(il%10)+':', label=atm+'['+str(icix)+',$'+inl.legends[icix][j]+'$]')
                    lines.append(linei)
                    il+=1
            if lcols_dn is not None and icix in lcols_dn and not foundOneYetdn:
                fname = env.case+'.gc'+str(icix)+'dn'
                if not os.path.exists(fname):
                    continue
                else:
                   foundOneYetdn=True
                cls = array(lcols_dn[icix])
                print('reading', fname, 'cols=', cls)
                dat = loadtxt(fname).T
                w = dat[0]
                Glat = array([dat[1+2*j]+dat[2+2*j]*1j for j in cls])
                for i,j in enumerate(cls):
                    linei,=axs.plot(w,-Glat[i].imag/pi, 'C'+str(il%10)+':', label='glat['+str(icix)+'dn,$'+inl.legends[icix][j]+'$]')
                    lines.append(linei)
                    il+=1
            
    if xrng[0]!=None and xrng[1]!=None:
        axs.set_xlim(xrng)
    elif xrng[0]!=None and xrng[1]==None:
        axs.set_xlim([xrng[0],om[-1]])
    elif xrng[1]!=None and xrng[0]==None:
        axs.set_xlim([om[0],xrng[1]])

    if yrng[0]!=None and yrng[1]!=None:
        axs.set_ylim(yrng)
    elif yrng[0]!=None and yrng[1]==None:
        axs.set_ylim(bottom=yrng[0])
    elif yrng[0]==None and yrng[1]!=None:
        axs.set_ylim(top=yrng[1])
    
    int_leg = InteractiveLegend(lines,axs.legend(),fig)
    fig.canvas.mpl_connect('pick_event', int_leg)
    
    #axs.legend(loc='best')
    #if args.g: axs.grid()
    
    show()
