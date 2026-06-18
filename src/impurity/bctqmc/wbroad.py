from wbroad1 import FindMaxDeriv, FindRangeToRemove, BroadenFunctions, GiveDh
from pylab import *

def fermi(x):
    return 1.0/(exp(x)+1.0)

p={"brd_exclude": 1.0, "brd_max_der": 2.0, "brd_beta": 2.0, "broad" : 2.0, "Wtypical": 1.0};

#dat = loadtxt('Stval_.4').T
dat = loadtxt('g00_.5').T

om = array(dat[0])
dh = GiveDh(om)
fl = np.ascontiguousarray(array(dat[1:,:]), dtype=np.float64)


md = FindMaxDeriv(om, fl, p['Wtypical'])
if md<0.2: md=0.2
ky = 0.015*md*p['broad']
kwidth = ky/(1 + ky)

iistart=0
x_remove = FindRangeToRemove(om, fl, p['brd_exclude'], p['brd_max_der'])
artificial_beta = p['brd_beta']
    
nplt=2
for l in range(1):
    istart,iend = iistart+l*6, iistart+(l+1)*nplt
    
    fig, axs = subplots(nplt//2,2, figsize=(12,8), squeeze=False)
    for i,ii in enumerate(range(istart,iend)):
        fi = fl[ii:ii+1,:]
        
        fo = BroadenFunctions(om,fi, kwidth, True, 200)
        
        if (x_remove[ii,1] < om[-1]):
            x_start = x_remove[ii,1]-(1.0+0.15*x_remove[ii,1])/artificial_beta
            print('x_remove=', x_remove[ii,1], 'x_start=', x_start)
            fo *= fermi((om[:]-x_start)*artificial_beta)
        if (x_remove[ii,0] > om[0]):
            x_start = x_remove[ii,0]+(1.0+0.15*abs(x_remove[ii,0]))/artificial_beta
            print('x_remove=', x_remove[ii,0], 'x_start=', x_start)
            fo *= 1-fermi((om[:]-x_start)*artificial_beta)
        
        
        print(f'Int0[{ii}]={sum(fi[0,:]*dh):10.7f} Int1[{ii}]={sum(fo[0,:]*dh):10.7f} x_rm=[{x_remove[ii,0]},{x_remove[ii,1]}]')
        ax = axs[i//2,i%2]
        
        line2, = ax.plot(om, fo[0], 'C'+str(ii%10), label='J['+str(ii)+']')
        line1, = ax.plot(om, fi[0], 'C'+str(ii%10), lw=0.5, label='I['+str(ii)+']')
        #ax.set_xlim(np.min(om), np.max(om))
        #ax.set_ylim(np.min(fo[0])*2, np.max(fo[0])*2)

        ax.legend(loc='best')
show()
