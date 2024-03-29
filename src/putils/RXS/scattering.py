from scipy import *
from scipy import optimize

class Scattering:
    def __init__(self,hkl,abc,lmbda,Sprint=False):
        self.h, self.k, self.l = hkl
        self.a, self.b, self.c = abc
        self.lmbda = lmbda
        self.Sprint = Sprint
        
        self.q = 2*pi*array([self.h/self.a, self.k/self.b, self.l/self.c])
        self.xi = 0.5*self.lmbda*((self.h/self.a)**2+(self.k/self.b)**2+(self.l/self.c)**2)
        
        self.theta0 = self.Get_Azimuthal_theta0()
    
    def ScatteringVector(self,theta):
        # solving the equation k*q = q**2/2
        # with k= 2*pi/lmbda * [sin(t)*cos(fi),sin(t)*sin(fi),cos(t)]
        if self.h==0 and self.k==0:
            if (cos(theta) - self.lmbda/(2*self.c))>1e-4:
                if self.Sprint: print 'If you scatter in z direction, there is only one solution, namely cos(theta)=lmbda/(2*c), and you did not satisfy it'
                return [[],[],[],[]]
            else:
                if self.Sprint: print 'phi is still undertermined'
                return [[],[],[],[]]
            
        if (theta==0):
            if ( abs(self.xi*self.c/self.l - 1)>1e-6):
                if self.Sprint: print 'No solution 1'
                return [[],[],[],[]]
            else:
                ek  = array([0,0,1.])
                ekp = array([-self.h/self.a,-self.k/self.b,1/self.lmbda-self.l/self.c])
                ekp *= 1/sqrt(sum(ekp**2))
                sig = cross(ek,ekp)
                sig *= 1.0/sqrt(sum(sig**2))
                p_in = cross(sig,ek)
                p_in *= 1.0/sqrt(sum(p_in**2))
                p_out = cross(sig,ekp)
                p_out *= 1.0/sqrt(sum(p_out**2))
                return [sig,p_in,p_out]
        
        sol = self.Find_fi(theta)
        if len(sol)==0:
            return [[],[],[],[]]
            
        pi_in=[]
        pi_out=[]
        sigm=[]
        eks=[]
        for i,(cos_fi,sin_fi) in enumerate(sol):
            ek=array([sin(theta)*cos_fi,sin(theta)*sin_fi,cos(theta)])
            ekp = (2*pi/self.lmbda*ek - self.q)/(2*pi/self.lmbda)
    
            sig = cross(ek,ekp)
            sig *= 1.0/sqrt(sum(sig**2))
    
            p_in = cross(sig,ek)
            p_in *= 1.0/sqrt(sum(p_in**2))
            
            p_out = cross(sig,ekp)
            p_out *= 1.0/sqrt(sum(p_out**2))
            
            pi_in.append(p_in)
            pi_out.append(p_out)
            sigm.append(sig)
            eks.append(ek)
            
            if self.Sprint: 
                print 'k_in = ', ek*2*pi/self.lmbda
                print 'k_out= ', ekp*2*pi/self.lmbda
                print 'k_in-k_out=', (ek-ekp)*2*pi/self.lmbda
                print 'q=', self.q
                print 'norm(ek)=', sum(ek**2)
                print 'norm(ekp)=', sum(ekp**2)
                print 'sig=', sig, sum(sig**2)
                print 'pi_in=', p_in, sum(p_in**2)
                print 'pi_out=', p_out, sum(p_out**2)
        return (sigm,pi_in,pi_out,eks)
        
    def Get_Azimuthal_theta0(self):
        which_axis=1
        self.psi_is_zero = zeros(3)
        if self.h==0:
            self.psi_is_zero[0]=1  # azimuthal angle is measured from the x-axis
            return arccos(self.lmbda*self.l/(2*self.c))
        elif self.k==0:
            self.psi_is_zero[1]=1  # azimuthal angle is measured from the y-axis
            return arccos(self.lmbda*self.l/(2*self.c))
        else:
            # should_vanish when  dot(ek,cross(eq,psi_is_zero))=0
            # will find solution numerically
            psi = [1,0,0]
            qr = [self.h/self.a,self.k/self.b,self.l/self.c]
            self.cr = cross(qr,psi)              #  q x psi
            self.psi_is_zero = cross(self.cr,qr) #  psi_is_zero = (q x psi) x q
            self.psi_is_zero *= 1./sqrt(sum(self.psi_is_zero**2))
            xs=[]
            ys=[]
            for it in linspace(0.00001,pi,300):
                y = self.what_theta(it)
                #print it/pi, y
                if y is not None:
                    ys.append(y)
                    xs.append(it)
            
            
            for i in range(len(xs)-1):
                if ys[i]*ys[i+1]<0:
                    thta = optimize.brentq(self.what_theta,a=xs[i],b=xs[i+1])
                    return thta

    def what_theta(self, theta):
        sol = self.Find_fi(theta)
        if len(sol)>0:
            _cos_fi = sol[0][0]
            _sin_fi = sol[0][1]
            _sin_th = sin(theta)
            return dot(self.cr, [_sin_th*_cos_fi,_sin_th*_sin_fi,cos(theta)])
        else:
            return None
            
    def Find_fi(self, theta):
        #xi=0.5*self.lmbda*((self.h/self.a)**2+(self.k/self.b)**2+(self.l/self.c)**2)
        eta=(self.xi-self.l/self.c*cos(theta))/sin(theta)
        ax = (self.h/self.a)**2 + (self.k/self.b)**2

        if self.h!=0:
            bx = eta*self.k/self.b
            cx = eta**2-(self.h/self.a)**2
            if bx**2 < ax*cx : return []
            sinf1 = (bx + sqrt(bx**2-ax*cx))/ax
            sinf2 = (bx - sqrt(bx**2-ax*cx))/ax
            sinfi = [ sinf1, sinf1, sinf2, sinf2]
            cosfi = [ sqrt(1-sinf1**2),-sqrt(1-sinf1**2), sqrt(1-sinf2**2),-sqrt(1-sinf2**2)]
        else:
            bx = eta*self.h/self.a
            cx = eta**2-(self.k/self.b)**2
            if bx**2 < ax*cx : return []
            cosf1 = (bx + sqrt(bx**2-ax*cx))/ax
            cosf2 = (bx - sqrt(bx**2-ax*cx))/ax
            cosfi = [ cosf1, cosf1, cosf2, cosf2]
            sinfi = [ sqrt(1-cosf1**2),-sqrt(1-cosf1**2), sqrt(1-cosf2**2),-sqrt(1-cosf2**2)]

        if eta > max(abs(self.h/self.a),abs(self.k/self.b)) or eta < -max(abs(self.h/self.a),abs(self.k/self.b)):
            if self.Sprint : print 'No solution 2', eta, self.h/self.a, self.k/self.b
            return []

        # is this solution good?
        good = []
        for i in range(4):
            if ( abs(self.h/self.a*cosfi[i] + self.k/self.b*sinfi[i]-eta)<1e-6): good.append(i)
        return [(cosfi[i],sinfi[i]) for i in good]
    
        
    def Azimuthal_single_valued(self, theta):
        if theta<self.theta0:
            return -1
        else:
            return 1
