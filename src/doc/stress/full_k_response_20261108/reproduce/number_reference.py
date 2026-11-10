"""Independent nonlinear trace for fixed-mu self-energy endpoint root finding.

At fixed Sigma(iw), all mu dependence is a scalar shift. Cache eigenVALUES of
H_W + U_W Sigma(iw) U_W^dagger once, then trace each resolvent as a scalar sum.
This is an endpoint reference, not the analytic density-response engine. No
linear-response dmu, endpoint eigenvector matching or screened information
enters this calculation. Validate against dense inverse traces before use.
"""
import numpy as np
from scipy.special import expit

class NumberReference:
    def __init__(self,model):
        self.meta=model.meta;self.kw=model.kw.copy()
        W=model.W;Q=np.setdiff1d(np.arange(model.eps.shape[1]),W)
        self.outside=model.eps[:,Q].copy();self.static=[];self.dynamic=[]
        self.omega=model.meta.omega.copy();self.weights=model.meta.weights.copy()
        for e,u in zip(model.eps,model.uc):
            u=u[W];h=np.diag(e[W]);ud=u.conj().T
            self.static.append(np.linalg.eigvalsh(h+u@model.sinf@ud))
            parts=[]
            for j in range(0,len(self.omega),32):
                parts.append(np.linalg.eigvals(h[None]+(u@model.sigma[j:j+32])@ud))
            self.dynamic.append(np.concatenate(parts))
        self.static=np.array(self.static);self.dynamic=np.array(self.dynamic)

    def __call__(self,mu):
        T=self.meta.temperature;z=mu+1j*self.omega
        number=expit((mu-self.static)/T).sum(axis=1)+expit((mu-self.outside)/T).sum(axis=1)
        dg=(1/(z[None,:,None]-self.dynamic)).sum(axis=-1)-(1/(z[None,:,None]-self.static[:,None])).sum(axis=-1)
        number+=2*T*(dg.real@self.weights)
        return float(2*self.kw@number)
