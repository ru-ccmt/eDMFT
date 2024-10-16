#!/usr/bin/env python
from numpy import *
from numpy import linalg
#from scipy.linalg import *
import utils
import os, re, sys

def Print(Us):
    for i in range(shape(Us)[0]):
        print('  ', end=' ')
        for j in range(shape(Us)[1]):
            print("%11.8f " % Us[i,j], end=' ')
        print()

def Get_BR1_DIR(case):
    file = case+'.outputd'
    found = False
    if os.path.exists(file) and os.path.getsize(file)>0:
        fi = open(file,'r')
        lines = fi.readlines()
        for il,line in enumerate(lines):
            if re.search('BR1_DIR',line):
                found = True
                break
        if found:
            a1 = list(map(float,lines[il+1].split()))
            a2 = list(map(float,lines[il+2].split()))
            a3 = list(map(float,lines[il+3].split()))
            return vstack( (a1,a2,a3) ).transpose()
            
    file = case+'.rotlm'
    if os.path.exists(file) and os.path.getsize(file)>0:
        fi = open(file,'r')
        next(fi)
        b1 = list(map(float, next(fi).split()))
        b2 = list(map(float, next(fi).split()))
        b3 = list(map(float, next(fi).split()))
        BR1 = vstack( (b1,b2,b3) )
        S2C = linalg.inv(BR1)*2*pi
        return S2C

def read_POSCAR(file):
    fi = open(file, 'r')
    next(fi)
    a = float(next(fi).split()[0])
    a1 = list(map(float, next(fi).split()))
    a2 = list(map(float, next(fi).split()))
    a3 = list(map(float, next(fi).split()))
    vaspbasis = vstack( (a1,a2,a3) ).transpose()
    return vaspbasis

def read_w2k_disp(strfile):
    fi = open(strfile, 'r')
    for i in range(4):
        next(fi)
    line = next(fi)  # should be first atoms
    # Here we assume that the first atom is displaced away from (0,0,0).
    # This should be generalized!
    ii = int(line[4:8])
    disp = zeros(3)
    for i,col in enumerate(12+13*arange(3)):
        disp[i] = float( line[col:col+10])
    return disp

def read_disp_yaml(filename):
    fi = open(filename, 'r')
    lines = fi.readlines()
    nat = int(lines[0].split()[1])
    
    disp=[]
    iat_disp=[]
    c=2
    for i in range(nat):
        iat = lines[c].split()[2]
        vdir = eval(lines[c+2])
        vdis = eval(lines[c+4])
        disp.append( vdis )
        iat_disp.append(iat)
        c += 5
        if lines[c][:6] != '- atom': break
    if lines[c][:7] != 'lattice':
        print('Something wrong reading', filename)
    c+=1
    vaspbasis=zeros((3,3))
    for i in range(3):
        ai = eval(lines[c+i][2:])
        vaspbasis[i,:] = array(ai)
    vaspbasis = vaspbasis.transpose()
    return (vaspbasis, disp)
    
def grepForce(filename):
    fi = open(filename, 'r')
    lines = fi.readlines()

    if lines[-1][:7]!=':FCHECK':
        print('ERROR: ', filename, 'does not contain forces at the end')
    frl = []  # forces in lattice basis
    c=2
    for i in range(2,1000):
        c+=1
        if lines[-i][:4]!=':FGL':
            break
        frl.append( list(map(float,lines[-i].split()[2:5])) )
    frl = frl[::-1]
    frc=[]  # forces in cartesian basis
    for i in range(c,1000):
        if lines[-i][:4]!=':FCA':
            break 
        frc.append( list(map(float,lines[-i].split()[2:5])) )
    frc = frc[::-1]
    return [array(frc), array(frl)]
    
au_2_angs = 0.52917721067
Ry2eV = 13.605693009

if __name__ == '__main__':
    POSCAR_file = 'POSCAR-001'
    
    (vaspbasis, dispVasp) = read_disp_yaml('disp.yaml')

    w2k = utils.W2kEnvironment()
    S2C = Get_BR1_DIR(w2k.case)
    print('S2C=')
    Print(S2C)
    S2C *= au_2_angs
    
    #vaspbasis = read_POSCAR(POSCAR_file)

    Vw2k = linalg.det(S2C)
    Vvasp = linalg.det(vaspbasis)
    print('Volume of wien2k primitive cell is %10.4f' % Vw2k)
    print('Volume of vasp   primitive cell is %10.4f' % Vvasp)
    if abs(Vw2k/Vvasp-1.)/Vw2k > 1e-4:
        print('ERROR : Volumes of the two unit cells are very different')
        sys.exit(1)

    # This transformation takes a vector expressed in cartesian coordinates in w2k and transforms it to
    # cartesian coordinates in Vasp.
    # Note that if vector is expressed in lattice vectors (call it r), than we just apply vaspbasis*r
    M = dot(vaspbasis, linalg.inv(S2C) )
    print("rotation matrix from w2kbasis_BR1 to vaspbasis is")
    Print(M)
    detM = linalg.det(M)
    print('and its determinant is', detM)
    print('check of unitarity : ')
    Print( dot(M.T,M) )

    force,forcel = grepForce(w2k.case+'.scf')
    print('force in cartesian coordinates extracted from case.scf file :')
    Print(force)
    
    disp = read_w2k_disp(w2k.case+'.struct')  # This must be generalized
    print('Displacement extracted from case.struct file', disp)
    disp_w2k = dot(S2C, disp)   # to cartesian coordinates
    dispVasp0 = dot(M,disp_w2k) # to Vasp lattice coordinates
    #print "Unrotated displacement vector in w2k in Ang"
    #print disp_w2k, 'norm is', norm(disp_w2k)
    print('Converted w2k displacement to Vasp basis is [Ang]:', "%12.8f "*3 % tuple(dispVasp0))
    print('Correct Vasp displacement from disp.yaml         :', "%12.8f "*3 % tuple(dispVasp[0]))

    #w2k_cartesian_force  = zeros( shape(forcel) )
    #a=S2C[2,2]
    #S2C_normalized = S2C/a
    #for i in range(len(forcel)):
    #    w2k_cartesian_force[i] = dot(S2C_normalized, forcel[i])
    #print 'w2k fore converted to cartesian coordinates:'
    #Print(w2k_cartesian_force)  # should be equal to force
    
    # change to Vasp units
    mRy_au_2_ev_ang = Ry2eV/1000./(au_2_angs)
    force *= mRy_au_2_ev_ang
    Vasp_force = zeros( shape(force) )
    for i in range(len(force)):
        Vasp_force[i] = dot(M, force[i])
    print('w2k forces in vaspbasis are in ev/Ang')
    Print(Vasp_force)
    
    fo = open('FORCE_SETS', 'w')
    print(len(force), file=fo)
    print('1', file=fo)
    print(file=fo)
    print(len(dispVasp), file=fo)
    for i in range(len(dispVasp)):
        print('%15.8f '*3 % tuple(dispVasp[i]), file=fo)
    for i in range(len(force)):
        print('%15.8f '*3 % tuple(Vasp_force[i]), file=fo)
    fo.close()
