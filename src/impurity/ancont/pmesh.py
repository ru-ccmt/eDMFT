# @Copyright 2007 Kristjan Haule
# 

from scipy import *

def MeshJoin(u, dwt, xt):
    tg = tan(u)
    return u-arctan(tg/xt)-dwt*xt*tg/(xt*xt+tg*tg)

def MakePositiveLogTanMesh(N, x0, x1, x2, alpha=0):
    eta = log(x1/x0)/(x2/x1-1)
    N1_min = int((1+eta*N)/(1+eta)+0.5)
    N1 = int((1+alpha)*N1_min)
    if (N1>N-2): N1=N-2
    N2 = N-N1
    xt = x2/x1
    dwt = N2*(log(x1)-log(x0))/(N1-1)

    if (alpha==0):
        ut=1e-5
    else:
        ut = optimize.brentq(MeshJoin, 1e-5, pi/2.-1e-5, args=(dwt,xt))
  
    a = arctan(tan(ut)/xt)
    b = dwt*sin(a)*cos(a)
    w = x1/tan(a)

    om=[]
    for i in range(N1):
        om.append(exp(log(x0)+i*(log(x1)-log(x0))/(N1-1)))
    for i in range(N2):
        om.append(w*tan(a+(i+1)*b/N2))

    return om

def MakeLogTanMesh(N, x0, x1, x2, alpha=0):
    # For building mesh which is logarithmic at small frequency and tan at large frequency
    # Mesh is symmetric around zero frequency
    N2 = N/2
    om1 = MakePositiveLogTanMesh(N2,x0,x1,x2,alpha)

    om=[]
    for i in range(N2):
        om.append(-om1[N2-i-1])
    
    for i in range(N2):
        om.append(om1[i])

    dh=[]
    dh.append(0.5*(om[1]-om[0]))
    for i in range(1,len(om)-1):
        dh.append(0.5*(om[i+1]-om[i-1]))
    dh.append(0.5*(om[-1]-om[-2]))
        
    return (array(om), array(dh))


if __name__ == '__main__':
    
    om, dh = MakeLogTanMesh(200, 1e-4, 1., 20.)

    for i in range(len(om)):
        print(i, om[i], " ", dh[i])

