#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
# 
"""Module analytically continues self-energy from imaginary axis to real axis
   Input:
       Mandatory:
          -sig   filename     -- input self-energy on imaginary axis
          -nom   number       -- number of matsubara points used 
          -beta  float       -- inverse temperature

       Optional (poles of self-energy determined from spectral function)
          -poles  [[y0,A0,B0],[y1,A1,B1],...]  --  poles will be forced at y0,y1,...
                       This is achieved by contribution to chi2+= Al/((x_i-yl)**2+w_i*Bl)
                       where x_i are positions and w_i weights of lorentzians to be determined
                       
       Optional (basis functions):   
          -b    float:[0-1]  -- b parameter to determin family of basis functions (0.8)
          -Ng   number       -- number of basis functions (12)
          -wexp float:[1-2]  -- wexp^i is position where basis functions are peaked (1.6)
          
       Optional (low energy fit):
          -ifit number[3-..] -- number of Matsubara points used to fit the lowe energy expansion
          
       Optional weights for chi2 computations:
          -alpha3 float[0-1] -- weight for the normalization in functional to be minimized (0.01)
          -alpha4 float[0-1] -- weight for the low energy expansion in the functional to be minimized (0.001)

   Output:
           Siom.nn  -- current function on imaginary axis (nn counts every 1000 steps)
           Sres.nn  -- current analytic continuation to real axis (nn counts every 1000 steps)
           
"""

import sys
from pylab import *
from scipy import *
import pmesh

import chi2f

def fparab(par, x, data):
    chi2=0;
    for i in range(len(x)):
        chi2 += ((par[0] + x[i]*par[1] + x[i]**2*par[2])-data[i])**2
    return chi2

class lowEnergyFermiLiquid(object):
    
    def __init__(self, ifit, iom0, isig, p0, valueAtUnity=0.05):
        xfit = iom0[0:ifit]
        yfit_i = array(isig)[0:ifit,1]
        yfit_r = array(isig)[0:ifit,0]
        self.expan_i = optimize.fmin_powell(fparab, [0,0,0], args=(xfit,yfit_i))
        self.expan_r = optimize.fmin_powell(fparab, [0,0,0], args=(xfit,yfit_r))
        print('estimated derivatives (real): ', self.expan_r)
        print('estimated derivatives (imag): ', self.expan_i)
        self.expan = self.expan_r.tolist() + self.expan_i.tolist()
        
        self.a0 = abs(self.expan_i[0])
        self.a1 = self.expan_r[1]
        self.a2 = abs(self.expan_i[2])
        if abs(self.a1) > sqrt(2*self.a0*self.a2): # parabola would otherwise becomes negative at some point.
            self.a1 = sqrt(2*self.a0*self.a2)*sign(self.a1)
        #valueAtUnity = 0.05   # function should becomes smaller than valueAtUnity at x=1.0
        p0n = ((self.a0 + self.a1 + 0.5*self.a2)/valueAtUnity - 1)*(2/self.a2)**4
        if p0>p0n: self.p0 = p0
        else: self.p0 = p0n # function should be smaller than valueAtUnity at x=1.0    
        #self.p0 = p0
        self.pf = sqrt(2.)*self.a2**2*self.p0**(3/4.)/(2*pi*(2+self.a0*self.a2*sqrt(self.p0))) # Normalization for the function
        
        #a1max = sqrt(2*self.a0*self.a2)
        #if (abs(self.a1)>0.5*a1max) : self.a1 = 0.5*a1max*(self.a1/abs(self.a1))
        
        print(self.expan)

    def Fun(self, x):
        return self.pf*(self.a0 + self.a1*x + 0.5*self.a2*x**2)/(1 + self.p0*(self.a2*x/2)**4)

    def Func(self, om):
        return array([self.Fun(x) for x in om])
        
    def Matsubara(self, iom):
        #gm = 1/self.expan_i[2]
        gm = 1/self.a2
        (x0, dh0) = pmesh.MakeLogTanMesh(500, 1e-5*gm, 300*gm, gm)
        F0 = array([self.Fun(x) for x in x0])
        F00 = self.Fun(0.0)
        weigh0 = abs(self.expan_i[0])/(pi*F00)

        datai=zeros(len(x0),dtype=float)
        datar=zeros(len(x0),dtype=float)
        F0i=[]
        F0r=[]
        for n in range(len(iom)):
            omn = iom[n]
            if (omn<0.3): subtract=1
            else: subtract=0
            for i in range(len(F0)):
                datai[i] = (F0[i]-F00*subtract)/(omn**2+x0[i]**2)
                datar[i] = F0[i]*x0[i]/(omn**2+x0[i]**2)
            wi = -(omn*integrate.simpson(datai, x0) + F00*pi*subtract)
            wr = -integrate.simpson(datar, x0)
            F0i.append(wi)
            F0r.append(wr)
            
        #F0r = zeros(len(iom), dtype=float)
        return (array(F0r), array(F0i), weigh0)
    
    def RealParts(self, om):
        gm = 1/self.expan_i[2]
        (x0, dh0) = pmesh.MakeLogTanMesh(500, 1e-5*gm, 300*gm, gm)
        F0 = array([self.Fun(x) for x in x0])
        F0j = array([self.Fun(x) for x in om])

        Frc = chi2f.kramskron(om, F0j, F0, x0, dh0, 0.0)

        
        spl = interpolate.splrep(om, real(Frc), k=3, s=0.0)
        dersr = [interpolate.splev(0.0, spl, der=0),
                interpolate.splev(0.0, spl, der=1),
                interpolate.splev(0.0, spl, der=2)]
        spl = interpolate.splrep(om, imag(Frc), k=3, s=0.0)
        dersi = [interpolate.splev(0.0, spl, der=0),
                interpolate.splev(0.0, spl, der=1),
                interpolate.splev(0.0, spl, der=2)]
        
        ders = array(dersr + dersi)
        print('der0=', ders.tolist())
        return (Frc, ders)


class lowEnergyInsulator(object):
    
    def __init__(self, ifit, iom0, isig, width):
        xfit = iom0[0:ifit]
        yfit_i = array(isig)[0:ifit,1]
        yfit_r = array(isig)[0:ifit,0]
        self.expan_i = optimize.fmin_powell(fparab, [0,0,0], args=(xfit,yfit_i))
        self.expan_r = optimize.fmin_powell(fparab, [0,0,0], args=(xfit,yfit_r))
        print('estimated derivatives (real): ', self.expan_r)
        print('estimated derivatives (imag): ', self.expan_i)
        self.expan = self.expan_r.tolist() + self.expan_i.tolist()
        
        self.a0 = abs(self.expan_i[0])
        self.width = width

    def Fun(self, x):
        return (self.width/pi)/(x**2 + self.width**2)
    
    def Funr(self, x):
        return 1/(x + self.width*1j)
    
    def Func(self, om):
        return array([self.Fun(x) for x in om])
        
    def Matsubara(self, iom):
        F00 = self.Fun(0.0)
        weigh0 = abs(self.expan_i[0])/(pi*F00)

        F0i=[]
        for n in range(len(iom)): F0i.append(-1/(iom[n]+self.width))

        F0r = zeros(len(iom), dtype=float)
        return (F0r, array(F0i), weigh0)
    
    def RealParts(self, om):
        Frc = array([self.Funr(x) for x in om])
        dersr = [0.0, 1/self.width**2, 0.0]
        dersi = [-1/self.width, 0.0, 2/self.width**3]
        ders = array(dersr + dersi)
        print('der0=', ders.tolist())
        return (Frc, ders)
        

class highEnergy(object):
    def __init__(self, b_, pos_):
        self.b = b_
        self.pos = pos_

    def F0(self, om, En):
        if (En*om > 0):
            return 1./(self.b*abs(En)*sqrt(pi))*exp(-self.b**2/4.-(log(om/En)/self.b)**2)    
        else:
            return 0.0
    def Func(self, om):
        F=[]
        for i in range(len(self.pos)):
            En = self.pos[i]
            F.append(array([self.F0(x, En) for x in om]))
        return F
        
    def Matsubara(self, iom):
        # Functions on imaginary axis
        Funcr=[]
        Funci=[]
        for i in range(len(self.pos)):
            En = self.pos[i]

            (x0, dh0) = pmesh.MakeLogTanMesh(300, 1e-2*abs(En), 2*abs(En), 7*abs(En))

            wb = []
            for x in x0:
                wb.append(self.F0(x+En, En))
            wb = array(wb)
            
            (gr, gi) = chi2f.matsum(En, iom, x0, dh0, wb)
            
            Funcr.append(array(gr))
            Funci.append(array(gi))
            
        return (Funcr, Funci)

    def RealParts(self, om):
        Func=[]
        derivs=[]
        for i in range(len(self.pos)):
            En = self.pos[i]

            (x0, dh0) = pmesh.MakeLogTanMesh(300, 1e-2*abs(En), 2*abs(En), 7*abs(En))

            wb = array([self.F0(x+En, En) for x in x0])
            
            F0j = array([self.F0(x, En) for x in om])

            Frc = chi2f.kramskron(om, F0j, wb, x0, dh0, En)


            
#            datai=zeros(len(x0), dtype=float)
#            Fre = zeros(len(om), dtype=float)
#            Fri = zeros(len(om), dtype=float)
#            nonzero=[]
#            for j in range(len(om)):
#                if (abs(om[j])<cutoff):
#                    nonzero.append(j)
#                    omj = om[j]
#                    wj = F0j[j] #wj = self.F0(omj, En)
#                    for k in range(len(x0)):
#                        datai[k] = (wb[k]-wj)/(omj-(x0[k]+En))
#  
#                    Fre[j] = integrate.simps(datai, x0) - wj*log(abs((x0[-1]+En-omj)/(omj-x0[0]-En)))
#                    Fri[j] = -pi*F0j[j] #self.F0(om[j], En)
  
  
#            if (i==10):
#                plot(om, real(Frc), 'ro-')
#                plot(om, imag(Frc), 'gs-')
#                show()
#                sys.exit(0)

#            spl = interpolate.splrep(om[nonzero[0]:nonzero[-1]+1], Fre[nonzero[0]:nonzero[-1]+1], k=3, s=0.0)
#            dersr = [interpolate.splev(0.0, spl, der=0),
#                     interpolate.splev(0.0, spl, der=1),
#                     interpolate.splev(0.0, spl, der=2)]
#
#            spl = interpolate.splrep(om[nonzero[0]:nonzero[-1]+1], Fri[nonzero[0]:nonzero[-1]+1], k=3, s=0.0)
#            dersi = [interpolate.splev(0.0, spl, der=0),
#                     interpolate.splev(0.0, spl, der=1),
#                     interpolate.splev(0.0, spl, der=2)]
            


            spl = interpolate.splrep(om, real(Frc), k=3, s=0.0)            
            dersr = [interpolate.splev(0.0, spl, der=0),
                     interpolate.splev(0.0, spl, der=1),
                     interpolate.splev(0.0, spl, der=2)]

            spl = interpolate.splrep(om, imag(Frc), k=3, s=0.0)
            dersi = [interpolate.splev(0.0, spl, der=0),
                     interpolate.splev(0.0, spl, der=1),
                     interpolate.splev(0.0, spl, der=2)]


            Func.append(Frc)
            ders = array(dersr + dersi)
            print('derx['+str(En)+']=', ders.tolist())
            derivs.append(array(ders))
        return (Func, derivs)

        

globalc = 0
global_print = 5000
def chi2(gweigh, vary, gwfix, fixed, sqmc, ifunr, ifuni, iom, intg, om, rfun, expand, ders, alphas):
    
    expand_sig = ones(len(expand), dtype=float)
    expand_sig[2] = 10
    
    sig0 = sqrt(0.5*(sqmc[0,2]**2+sqmc[0,3]**2))
    chi = 0
    for im in range(len(iom)):
        gi = 0.0
        gr = 0.0
        for i in range(len(gweigh)):
            gi += ifuni[vary[i],im]*gweigh[i]
            gr += ifunr[vary[i],im]*gweigh[i]
        for i in range(len(fixed)):
            gi += ifuni[fixed[i],im]*gwfix[i]
            gr += ifunr[fixed[i],im]*gwfix[i]
        
        chi += ((sqmc[im,0]-gr)*sig0/sqmc[im,2])**2+((sqmc[im,1]-gi)*sig0/sqmc[im,3])**2

    nrm = (sum(gweigh)+sum(gwfix))/intg
    chi3 = (nrm-1.)**2

    tders = zeros(len(ders[0]), dtype=float)
    for i in range(len(gweigh)):
        tders += ders[vary[i]]*gweigh[i]
    for i in range(len(fixed)):
        tders += ders[fixed[i]]*gwfix[i]

    chi4 = 0
    for i in range(len(tders)):
        chi4 += ((tders[i]-expand[i])/expand_sig[i])**2


    print('chi2=', chi, 'chi3=', chi3, 'chi4=', chi4)
    
    global globalc
    global global_print
    globalc = globalc + 1
    if (globalc%global_print == 0):
        ii = globalc/global_print
        print(ii, 'chi2=', chi, 'chi3=', alphas[0]*chi3, 'chi4=', alphas[1]*chi4, 'nrm=', nrm)
        print('der=', tders.tolist())
        print('exp=', expand)
        
        fh = open('Sres.'+str(ii), 'w')
        for im in range(len(om)):
            csum=0
            for i in range(len(gweigh)):
                csum += rfun[vary[i],im]*gweigh[i]
            for i in range(len(fixed)):
                csum += rfun[fixed[i],im]*gwfix[i]
            print(om[im], csum, file=fh)

        fn = open('Siom.'+str(ii), 'w')
        for im in range(len(iom)):
            gc = 0j
            for i in range(len(gweigh)):
                gc += (ifunr[vary[i],im]+1j*ifuni[vary[i],im])*gweigh[i]
            for i in range(len(fixed)):
                gc += (ifunr[fixed[i],im]+1j*ifuni[fixed[i],im])*gwfix[i]
            print(iom[im], gc.real, gc.imag, sqmc[im,0], sqmc[im,1], file=fn)

    
    return chi + alphas[0]*chi3 + alphas[1]*chi4


def chi2n(gweigh, vary, gwfix, fixed, sqmc, ifunr, ifuni, iom, intg, om, rfun, rfunc, expand, ders, alphas, gpos, poles):

    expand_sig = ones(len(expand), dtype=float)
    expand_sig[2] = 100
    expand_sig[4] = 100
    expand_sig[5] = 100
    
    (chi, nrm, chi4) = chi2f.chi2(gweigh, vary, gwfix, fixed, sqmc, ifunr, ifuni, ders, expand, expand_sig)

    nrm /= intg
    chi3 = (nrm-1.)**2

    # chi5 += Al/((x_i-yl)**2+w_i*Bl)
    chi5 = 0
    for i in range(len(gpos)):
        for l in range(len(poles)):
            chi5 += poles[l][1]/(((gpos[i]-poles[l][0])/poles[l][2])**4+gweigh[i])

        #print gpos[i], gweigh[i], poles[l][1]/(((gpos[i]-poles[l][0])/poles[l][2])**4+gweigh[i])
        #print gpos[i],  poles[0][1]/((gpos[i]-poles[0][0])**2+gweigh[i]*poles[0][2])

    tot_chi2 = chi + alphas[0]*chi3 + alphas[1]*chi4 + chi5
#    print 'chi2=', chi, 'chi3=', chi3, 'chi4=', chi4
    
    global globalc
    global global_print
    globalc = globalc + 1
    if ((globalc-1)%global_print == 0):
        ii = globalc/global_print
        print(ii, 'chi2=', chi, 'chi3=', alphas[0]*chi3, 'chi4=', alphas[1]*chi4, 'chi5=', chi5, 'nrm=', nrm)


        tders = zeros(len(ders[0]), dtype=float)
        for i in range(len(gweigh)):
            tders += ders[vary[i]-1]*gweigh[i]
        for i in range(len(fixed)):
            tders += ders[fixed[i]-1]*gwfix[i]
        

        
        print('der=', tders.tolist())
        print('exp=', expand.tolist())
        
        fh = open('Sres.'+str(ii), 'w')
        for im in range(len(om)):
            csum=0
            zsum=0
            for i in range(len(gweigh)):
                csum += rfun[vary[i]-1,im]*gweigh[i]
                zsum += rfunc[vary[i]-1,im]*gweigh[i]
            for i in range(len(fixed)):
                csum += rfun[fixed[i]-1,im]*gwfix[i]
                zsum += rfunc[fixed[i]-1,im]*gwfix[i]
            print(om[im], csum, zsum.real+sinfty, zsum.imag, file=fh)

        fn = open('Siom.'+str(ii), 'w')
        for im in range(len(iom)):
            gc = 0j
            for i in range(len(gweigh)):
                gc += (ifunr[vary[i]-1,im]+1j*ifuni[vary[i]-1,im])*gweigh[i]
            for i in range(len(fixed)):
                gc += (ifunr[fixed[i]-1,im]+1j*ifuni[fixed[i]-1,im])*gwfix[i]
            print(iom[im], gc.real, gc.imag, sqmc[im,0], sqmc[im,1], file=fn)

        fm = open('coeff.'+str(ii), 'w')
        print('nom=', params['nom'][0], end=' ', file=fm)
        print('beta=', params['beta'][0], end=' ', file=fm)
        print('b=', params['b'][0], end=' ', file=fm)
        print('Ng=', params['Ng'][0], end=' ', file=fm)
        print('wexp=', params['wexp'][0], end=' ', file=fm)
        print('ifit=', params['ifit'][0], end=' ', file=fm)
        print('alpha3=', params['alpha3'][0], end=' ', file=fm)
        print('alpha4=', params['alpha4'][0], file=fm)
        for i in range(len(gweigh)):
            print(i, gpos[i], gweigh[i], file=fm)
            
 

    return tot_chi2
#    return chi + alphas[0]*chi3 + alphas[1]*chi4

def compare(gweigh, sqmc, ifunr, ifuni, iom, intg):
    sig0 = sqrt(0.5*(sqmc[0,2]**2+sqmc[0,3]))
    chi = 0
    fgr=[]
    fgi=[]
    for im in range(len(iom)):
        gi = 0.0
        gr = 0.0
        for i in range(len(gweigh)):
            gi += ifuni[i,im]*gweigh[i]
            gr += ifunr[i,im]*gweigh[i]
        fgr.append(gr)
        fgi.append(gi)
#        print iom[im], gr, gi, sqmc[im,0], sqmc[im,1], sqmc[im,2], sqmc[im,3]

    plot(iom, fgr, iom, fgi)
    plot(iom, sqmc[:len(iom),0], iom, sqmc[:len(iom),1])
    nrm = sum(gweigh)/intg
    print('nrm = ', nrm)
    show()


def log_functions(om, rfun):
    f = open('bfunc.dat', 'w')
    for i in range(len(om)):
        print(om[i], end=' ', file=f)
        for j in range(len(rfun)):
            print(rfun[j,i], end=' ', file=f)
        print(file=f)
    f.close()

def log_functions_im(iom, ifunr, ifuni):
    f = open('ifunc.dat', 'w')
    for i in range(len(iom)):
        print(iom[i], end=' ', file=f)
        for j in range(len(ifunr)):
            print(ifunr[j,i], ifuni[j,i], end=' ', file=f)
        print(file=f)
    f.close()

    
def log_history(filename, argv):
    f = open(filename, 'a')
    for a in argv:
        print(a, end=' ', file=f)
    print(file=f)

        

if __name__ == '__main__':

    # from command line
    arguments = sys.argv[1:]
    # parameters in this module
    params = {'sig'   : ['Sig.out', '# input self-energy on imaginary axis'],
              'FermiLiquid' : [True, '# Low energy expansion of a Fermi liquid of Mott insulator'],
              'nom'   : [100,       '# number of matsubara points used '],
              'lcut'  : [0.0,       '# the lowest frequency lorentzian position'],
              'beta'  : [100.,      '# inverse temperature'],
              'b'     : [0.8,       '# b parameter to determin family of basis functions'],
              'Ng'    : [15,        '# number of basis functions'],
              'wexp'  : [1.5,       '# wexp^i is position where basis functions are peaked'],
              'ifit'  : [4,         '# number of Matsubara points used to fit the lowe energy expansion'],
              'alpha3': [0.01,      '# weight for the normalization in functional to be minimized'],
              'alpha4': [0.001,     '# weight for the low energy expansion in the functional to be minimized'],
              'p0'    : [1,         '# how fast should the low energy function fall off'],
              'poles' : [[],        '# additional poles'],
              'vunity': [0.05,      '# value of low-energy function at unity'],
              'L0'    : [15,        '# High energy cut-off for plotting'],
              'N0'    : [400,       '# Number of frequency points for plotting'],
              'a0'    : [5e-3,      '# Frequency mesh start'],
              'b0'    : [0.5,       '# Frequency mesh parameter'],
              'rps'   : [1,         '# The lowest lorentzian is usualy at pi*T. This factor can bring it closer to zero.'],
              'maxsteps':[7500,     '# Maximum number of function evaluations in minimization routine'],
              }
    
    if (len(arguments)==0 or arguments[0]=="-h" or arguments[0]=="--help"): # if no arguments is given
        print(__doc__)
        sys.exit(0)
    log_history('history.dat', sys.argv)
    
    for i in range(len(arguments)/2):
        if (arguments[2*i][1:] in list(params.keys())):
            var = arguments[2*i][1:]             # variable name
            if (type(params[var][0]) is str):
                params[var][0] = arguments[2*i+1]    # new value is assigned to string parameters
            else:
                params[var][0] = eval(arguments[2*i+1]) # new value is assigned to non-string parameters
            
    for var in list(params.keys()):                    # We go over all parameters and evaluate them
        print('%10s  %-10s  %s' % ('-'+var, params[var][0], params[var][1]))
        val = params[var][0]
        if (type(val) is str): val = '"'+str(val)+'"'
        else: val = str(val)
        exec(var + '=' + val)     # variable is set to its current value

    # to compute derivatives (not important)
#    cutoff = 0.1

    # reading qmc data
    f = open(sig)
    iom0 = []
    sqmc0 = []
    for line in f:
        spl = list(map(float,line.split()))
        iom0.append(spl[0])
        sqmc0.append([spl[1],spl[2],spl[3],spl[4]])
    f.close()

    # subtracting Sigma_{infty}
    sinfty = sqmc0[-1][0]
    for i in range(len(sqmc0)):
        sqmc0[i][0] -= sinfty
    sqmc0 = array(sqmc0)
    # integral of self-energy from sig-infty

    Nstart = int(size(iom0)*0.1)
    if Nstart<nom : Nstart = nom+1
    Nend = int(size(iom0)*0.7)
    htail = [-sqmc0[i][1]*iom0[i] for i in range(Nstart,Nend)]
    intg = sum(htail)/len(htail)    
    #intg = -sqmc0[-1][1]*iom0[-1]
    print('integral of sigma is ', intg, 'omstart,omend=', iom0[Nstart], iom0[Nend])

    T = 1./beta

    # Position of gaussians
    gpos=[]
    for i in range(Ng):
        x0 = -rps*pi*T*wexp**(Ng-1-i)
        if abs(x0)>params['lcut'][0] : 
            gpos.append(-rps*pi*T*wexp**(Ng-1-i))
        #gpos.append(-pi*T*wexp**(Ng-1-i))
    for i in range(Ng):
        x0 = rps*pi*T*wexp**i
        if abs(x0)>params['lcut'][0] :
            gpos.append(rps*pi*T*wexp**i)
#        gpos.append(pi*T*wexp**i)

    print('positions=', gpos)
    
    # Creat real-frequency mesh
    (om, dh) = pmesh.MakeLogTanMesh(N0, 5e-3, 0.5, L0)
    (om, dh) = pmesh.MakeLogTanMesh(N0, a0, b0, L0)

    # Create Matsubara mesh
    iom = [(2*n+1)*pi/beta for n in range(nom)]
    if (abs(iom[0]-iom0[0])<T):  ### HERE WAS A BUG
        sqmc = array(sqmc0[0:len(iom),:], order='F')
    else:
        sqmc = array(sqmc0[1:len(iom)+1,:], order='F')
        
    # Fitting low energy on imaginary axis
    # The broad function at low frequency
    # is prepared on real axis
    if (FermiLiquid):
        lwe = lowEnergyFermiLiquid(ifit, iom0, sqmc0, p0, vunity)
    else:
        lwe = lowEnergyInsulator(ifit, iom0, sqmc0, 2*T)

    expandr = [lwe.expan_r[0], lwe.expan_i[1], -lwe.expan_r[2]]
    expandi = [lwe.expan_i[0], -lwe.expan_r[1], -lwe.expan_i[2]]
    print('Expecting real-axis derivatives (real):', expandr)
    print('Expecting real-axis derivatives (imag):', expandi)
    expand = array(expandr + expandi, order='F')

    # The low energy basis function
    A0 = lwe.Func(om)

    print('INTEGRAL of A0 IS: ', integrate.simpson(A0,om))
    
    # The low energy function on Matsubara axis
    (F0r, F0i, weight0) = lwe.Matsubara(iom)
    weight0 = weight0
    
    # Real parts
#    (F0_real, ders0) = lwe.RealParts(om, cutoff)
    (Gw0, ders0) = lwe.RealParts(om) #, cutoff)

    # all gaussians 
    hge = highEnergy(b, gpos)
    # basis functions
    Arest = hge.Func(om)
    # rest of the functions on matsubara axis
    (Frest_r, Frest_i) = hge.Matsubara(iom)
    # Real parts
#    (Frest_real, ders_rest) = hge.RealParts(om, cutoff)
    (Gw_rest, ders_rest) = hge.RealParts(om) #, cutoff)


    # Assigning weights to each basis functions
    gweigh = []
    for i in range(len(gpos)):
        gweigh.append( abs(gpos[i])**1 )
    gweigh = array(gweigh)
    sm = sum(gweigh)
    gweigh *= (intg-weight0)/sm  # normalizing to start with the correct weight

    # Putting them rogether
    ifunr = array(Frest_r   + [F0r], order='F')
    ifuni = array(Frest_i   + [F0i], order='F')
    rfun  = array(Arest     + [A0])
    rfunc = array(Gw_rest   + [Gw0])
    ders  = array(ders_rest + [ders0], order='F')

    # Preparing index arrays for minimization
    gwfix= array([weight0], order='F')
    abounds = [(1e-12,None)] * len(gpos)
    alphas = [alpha3,alpha4]

    vary=array(list(range(1,len(gpos)+1)), order='F', dtype=int)
    fixed=array([len(gpos)+1], order='F', dtype=int)

    log_functions(om, rfun)
    log_functions_im(iom, ifunr, ifuni)

    # Minimization routine is called
    (gweigh, fmin, dinf) = optimize.fmin_l_bfgs_b(chi2n, gweigh, approx_grad=1,
                                args=(vary, gwfix, fixed, sqmc, ifunr, ifuni, iom, intg, om, rfun, rfunc, expand, ders, alphas, gpos, poles),
                                                  bounds=abounds, maxfun=maxsteps)

    for i in range(len(gpos)):
        print(i, gpos[i], gweigh[i])

    f1 = open("gim.dat", 'w')
    for im in range(len(iom)):
        csumr=0
        csumi=0
        for i in range(len(vary)):
            csumr += ifunr[vary[i]-1,im]*gweigh[i]
            csumi += ifuni[vary[i]-1,im]*gweigh[i]
        for i in range(len(fixed)):
            csumr += ifunr[fixed[i]-1,im]*gwfix[i]
            csumi += ifuni[fixed[i]-1,im]*gwfix[i]
        print(iom[im], csumr, csumi, file=f1)

    f2 = open("gre.dat", 'w')
    for im in range(len(om)):
        csum=0
        for i in range(len(vary)):
            csum += rfun[vary[i]-1,im]*gweigh[i]
        for i in range(len(fixed)):
            csum += rfun[fixed[i]-1,im]*gwfix[i]
        print(om[im], csum, file=f2)


    f3 = open('Sig.out', 'w')
    for im in range(len(om)):
        zsum=0
        for i in range(len(gweigh)):
            zsum += rfunc[vary[i]-1,im]*gweigh[i]
        for i in range(len(fixed)):
            zsum += rfunc[fixed[i]-1,im]*gwfix[i]
        print(om[im], zsum.real+sinfty, zsum.imag, file=f3)
