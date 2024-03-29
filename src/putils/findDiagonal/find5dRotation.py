#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
from numpy import *
from scipy import linalg
import copy
import sys

def mprint(Us):
    for i in range(shape(Us)[0]):
        for j in range(shape(Us)[1]):
            print("%11.8f %11.8f " % (real(Us[i,j]), imag(Us[i,j])), end=' ')
        print()

def MakeOrthogonal(a, b, ii):
    a -= (a[ii]/b[ii])*b
    a *= 1/sqrt(dot(a,a.conj()))
    b -= dot(b,a.conj())*a
    b *= 1/sqrt(dot(b,b.conj()))
    return (a,b)

def StringToMatrix(cfstr):
    mm=[]
    for line in cfstr.split('\n'):
        line = line.strip()
        if line:
            data = array(list(map(float,line.split())))
            mm.append( data[0::2]+data[1::2]*1j )
    mm=matrix(mm)
    return mm

def RealPhase(vec):
    for j in range(len(vec)):
        v = vec[j]
        imax = 0
        vmax = abs(v[imax])
        for i in range(len(v)):
            if abs(v[i])>vmax:
                vmax=abs(v[i])
                imax = i
        vec[j,:] = array(v)*abs(v[imax])/v[imax]
    return vec

def to_normalize(a):
    return 1./sqrt(abs(dot(conj(a), a)))

def swap(a,b):
    an = copy.deepcopy(a)
    bn = copy.deepcopy(b)
    return (bn,an)

def findMax(v):
    av=abs(v)
    return sorted(range(len(v)), key=lambda a: -av[a])

def GiveNewT2C(Hc, T2C, KeepOrder=False):
    
    ee = linalg.eigh(Hc)
    Es = ee[0]
    Us = matrix(ee[1])
    Es = Es[::-1]
    Us = Us[:,::-1]
    #print 'In Eigensystem:'
    #mprint(Us.H * Hc * Us)
    # Us.H * Hc * Us  === diagonal
    dim = len(Hc)
    
    print('Eigenvalues=', Es.tolist())

    if KeepOrder: # If we want the output to be very similar to input
        # Next we will sort eigenvectors so that the matrix of eigenvectors is as close as possible to unity
        N = len(Us)
        R = abs(Us)
        mprint(R)
        ind = zeros(N, dtype=int)
        for l in range(N):
            ii = argmax(R)
            i,j = int(ii/N), ii%N  # maximum index in the matrix
            #print(ii, int(ii/N), ii%N, R[i, j])
            ind[j] = i
            #print('setting ind['+str(j)+']=',i)
            R[:,j]=0
            #print('ind=', ind)
        print('ind to rearange eigenvectors=', ind)
        if sum(arange(N)-array(sorted(ind)))==0:
            # means that ind is a permutation, hence its sorted version is just arange
            Usn = Us.copy()
            for j in range(N):
                Usn[:,ind[j]] = Us[:,j]
            Esn = Es.copy()
            for j in range(N):
                Esn[ind[j]] = Es[j]
            #print('Us=')
            #mprint(Us)
            Us = Usn
            Es = Esn
            #print('Us=')
            #mprint(Us)
        else:
            print('WARNING: Rearangement did not succeed ind=',ind)
    mprint(Us)
    print('Eigenvalues=', Es.tolist())
    print()

    #checkDiag = Us.H * Hc * Us
    #mprint(checkDiag)

    
    for i0 in range(0,dim,2):
        vects = Us[:,i0:(i0+2)]
        print('vects=')
        mprint(vects)
        vtc = vects.T.conj()
        #print 'vects^H='
        #mprint(vtc)
        ind0=findMax(vects[:,0])  # Finds which components of the eigenvector are large for the first atomic state
        ind1=findMax(vects[:,1])  # Finds which components of the eigenvector are large for the second atomic state
        print('ind0=', ind0, 'ind1=', ind1)
        
        # The two largest components will be analized
        if ind0[0]!=ind1[0]: 
            j0,j1 = min(ind0[0],ind1[0]), max(ind0[0],ind1[0])
        else:
            j0,j1 = min(ind0[1],ind1[0]), max(ind0[1],ind1[0])
            
        # We will make sure that the largest components of the two eigenvectors are maximally orthogonal
        O = hstack((vtc[:,j0],vtc[:,j1]))
        print('O=')
        mprint(O)
        
        (u_,s_,v_) = linalg.svd(O)
        print('S=', s_.tolist())
        
        m = min(shape(u_)[1],shape(v_)[0])
        R = dot(u_[:,:m],v_[:m,:])
        #print 'R=', R
        
        vectn = dot(vects,R)
        #print 'vectn^H='
        #mprint(transpose(conj(vectn)))
        
        Us[:,i0:(i0+2)] = vectn[:,:]
        
    #Us = u_ * s_ * v_
    #print 'Eigenvalues'
    #print "%10.5f "*len(Es) % tuple(Es)
    
    print('Transformation in crystal harmonics=')
    mprint(Us)
    print()
    #print 'shape(T2C)=', shape(T2C)
    #mprint(T2C)
    #print
    final = Us.T*T2C
    
    final = array(final)
    final2 = RealPhase(final)
    final=copy.deepcopy(final2)
    return final

def Check(final, T2C, Hc):
    # the modified final transofrmation is rotated back to t2g-eg base to see how weell diagonal remains
    Us_new = transpose(matrix(final)*T2C.H)
    print('Check-diagonal:')
    mprint(Us_new.H * Hc * Us_new)
    
    print('Check unitary:')
    mprint( matrix(final) * matrix(final).H )
    print()

def CheckDet(final, T2Crest):
    totalfinal = vstack((final,T2Crest))
    Det = linalg.det(totalfinal)
    print('Determinant=', Det)
    if abs(Det+1)<1e-3:
        print('Determinant is -1, you need to change an eigenvector, to make the rotation proper!')
    return Det



if __name__ == '__main__':
    
    #strHc1="""
   #-1.52945958     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000       0.06742707     0.01545415       0.00000000     0.00000000       0.00000000     0.00000000   
    #0.00000000     0.00000000      -1.52945958     0.00000000       0.06742705    -0.01545415       0.00000000     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000   
    #0.00000000     0.00000000       0.06742705     0.01545415      -2.30650667     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000   
    #0.06742707    -0.01545415       0.00000000     0.00000000       0.00000000     0.00000000      -2.30650668     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000   
    #0.00000000     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000      -2.28241642     0.00000000       0.00000000     0.00000000   
    #0.00000000     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000      -2.28241642     0.00000000   
    #"""
    #
    #strHc2="""
    #0.13290679     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000       0.05659197    -0.02929068       0.00000000     0.00000000       0.00000000     0.00000000   
    #0.00000000     0.00000000       0.13290691     0.00000000       0.05659196     0.02929068       0.00000000     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000   
    #0.00000000     0.00000000       0.05659196    -0.02929068      -1.08239101     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000   
    #0.05659197     0.02929068       0.00000000     0.00000000       0.00000000     0.00000000      -1.08239098     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000   
    #0.00000000     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000      -0.70798000     0.00000000       0.00000000     0.00000000   
    #0.00000000     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000      -0.70797997     0.00000000   
    #"""
    #
    #
    #strT2C="""
    #-0.00000000 -0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.87675550  0.00000000    0.00000000  0.00000000   -0.33968059 -0.01634050   -0.00000000 -0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.33968059  0.01634050  
    # 0.33968059 -0.01634050   -0.00000000 -0.00000000   -0.00000000 -0.00000000    0.00000000 -0.00000000   -0.33968059  0.01634050   -0.00000000  0.00000000    0.87675550 -0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000    0.00000000 -0.00000000  
    # 0.61995975  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000 -0.00000000   -0.61995975  0.00000000    0.00000000 -0.00000000   -0.48038090 -0.02310896    0.00000000  0.00000000    0.00000000  0.00000000   -0.00000000  0.00000000  
    # 0.00000000  0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000   -0.48038090  0.02310896   -0.00000000 -0.00000000   -0.61995976 -0.00000000    0.00000000  0.00000000   -0.00000000 -0.00000000   -0.00000000 -0.00000000    0.61995976 -0.00000000  
    # 0.00000000  0.00000000    1.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000  
    # 0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    1.00000000  0.00000000    0.00000000  0.00000000  
    #"""
    #strT2Crest="""
    # 0.00000000  0.00000000    0.00000000  0.00000000    1.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000  
    # 0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    1.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000  
    # 0.70710679  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.70710679  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000  
    # 0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.70710679  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000    0.70710679  0.00000000  
    #"""
   
    if len(sys.argv)<2:
        print()
        print("Give input file which conatins impurity levels (strHc) and transformation (strT2C)")
        print()
        sys.exit(0)
    
    fpar=sys.argv[1]
    exec(compile(open(fpar, "rb").read(), fpar, 'exec'))
    
    Hc = StringToMatrix(strHc)
    print('shape(Hc)=', shape(Hc))
    T2C0=StringToMatrix(strT2C)
    print('shape(T2C0)=', shape(T2C0))
    T2C = T2C0[:len(Hc),:]
    print('shape(T2C)=', shape(T2C))
    T2Crest = T2C0[len(Hc):,:]
    print('shape(T2Crest)=', shape(T2Crest))
    
    final = GiveNewT2C(Hc, T2C)
    
    print('Rotation to input : ')
    mprint( final )
    mprint( T2Crest )
