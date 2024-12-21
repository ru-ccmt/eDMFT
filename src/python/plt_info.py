#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import glob
import utils
Ry2eV = 13.60569193

if __name__ == '__main__':
    
    dat = np.loadtxt('info.iterate')
    ind=[]
    istr=0
    for j in range(len(dat)-1):
        if dat[j+1][2]<=dat[j][2]:
            ind.append(j)
        if dat[j+1][0]<dat[j][0]:
            ind.clear()
            istr=j+1
    ind.append(len(dat)-1)
    print('# ind=', ind)
    
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, sharex=False)
    
    ax1.plot(dat[istr:,0],dat[istr:,8], 'o-', label='nlat')
    ax1.plot(dat[istr:,0],dat[istr:,9], 'o-', label='nimp')
    ax1.legend(loc='best')
    ax1.set_xlabel('dft+dmft iterations')
    ax1.xaxis.set_label_position('top')
    ax1.xaxis.tick_top()
    
    cdat = np.array([dat[j] for j in ind])
    itt = np.array(cdat[:,1],dtype=int)
    nlat = cdat[:,8]
    nimp = cdat[:,9]
    dEtot = (cdat[:,5]-cdat[-1,5])*Ry2eV
    dFtot = (cdat[:,7]-cdat[-1,7])*Ry2eV

    #ax2.sharex(ax3)
    ax2.set_xticks(itt,[])
    ax3.set_xticks(itt,itt)
    #print('itt=', itt)
    
    ax2.plot(itt,nlat, 'o-', label='nlat')
    ax2.plot(itt,nimp, 'o-', label='nimp')
    ax2.legend(loc='best')

    ax3.plot(itt,dEtot, 'o-', label='Etot[eV]')
    ax3.plot(itt,dFtot, 'o-', label='Ftot+T*Simp[eV]')

    print('Etot=', (cdat[:,5]*Ry2eV).tolist())
    print('Ftot+T*Simp=', (cdat[:,7]*Ry2eV).tolist())
    
    ax3.legend(loc='best')
    ax3.set_ylabel('E[eV]')
    ax3.set_xlabel('dmft_iterations')
    plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0.02)
    plt.show()
    
