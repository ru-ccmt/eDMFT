#!/usr/bin/env python
""" Given a set of neighbors around the central atom it finds the best 
   local coordinate axis for hybridization, i.e., the axis in which 
   hybridization is maximally diagonal.
   
  Most algorithm is inside function
     FindCageBasis(neighs, matrix, log, max_bond_variance=0.7, max_angle_variance=1.0)
  which needs list of neighbors and the lattice matrix, which  transform 
  fractional coordinates to cartesian basis.
  
  If run from command line, it reads w2k outputnn file. If run from cif2indmfl, it can be
  use pymatgen to determine neighbors and determine polyhedra.
"""
from numpy import * 
import sys
from numpy import linalg
import numpy as np
import heapq
    
def FindCageBasis(neighs, matrix, log, max_bond_variance=0.85, max_angle_variance=1.0):
    """Find the best directions for the local coordinate system, which is non-orthogonal.
    """
    def Find_two_closest_smallest_values(ph0, where):
        """This is used for trigonal prism. 
        Given all angles, we want to find the vectors to the triangle of the trigonal prism. 
        We assume first neigbor (closes neighbor) is in one triangle, and we are searching 
        for the remaining two neighbors (call them j1,j2) in the triangle. 
        We search for j1 and j2 such that phi(0,j1) ~ phi(0,j2) and also want 
        phi(0,j1) and phi(0,j2) to be small, because the opposite triangle will have 
        180-phi(0,j1) and 180-phi(0,j2), which we want to avoid.
        """
        cm=2*pi
        j1=-1
        #print('starting with ph=', ph0*180/pi)
        wh = np.where(ph0<pi)[0]
        for j in wh:
            pj = ph0[j]
            ph0[j]=4*pi
            #print('ph0=', ph0, 'pj=', pj)
            mm = min(abs(ph0-pj)) + abs(pj)/10.
            ph0[j]=pj
            #print(j, 'pj=',pj*180/pi, 'mm=', mm*180/pi)
            if mm<cm:
                cm=mm
                j1=j
        #print('cm=', cm, 'j1=', j1)
        pj=ph0[j1]
        ph0[j1]=3*pi
        j2 = argmin(abs(ph0-pj))
        ph0[j1]=pj
        return (j1,j2)
    # supported cages which are currently detected
    # (type, number of vertices, angles)
    cases = [('cube',8, [arccos(1/3.)]*3+[arccos(-1/3.)]*3+[pi]),
            ('octahedra',6, [pi/2]*4+[pi]),
            ('tetrahedra',4, [arccos(-1/3.)]*3),
            ('square-piramid',5, [pi/2]*4),
            ('cuboctahedron',12,[pi/3]*4+[pi/2]*2+[2*pi/3]*4+[pi]),
            ('triangular-bipyramid', 5, [pi/2]*2+[2/3*pi]*2),
            ('truncated tetrahedron',12,[arccos(7/11)]*3+[arccos(-1/11)]*4+[arccos(-5/11)]*2+[arccos(-9/11)]*2)]
            # 
        # NEED:   square antiprism
    cases_outside=[('peak-of-tetrahedron',3, [pi/2]*2),
                   ('peak-of-square-piramid',4,[pi/3]*3),
                   ('peak-of-hexagonal-piramid',6,[pi/6]*5)    ]

    #for case in cases:
    #    print(case[0], case[1], array(case[2])*180/pi)
    
    # how many neighbors with the same distance we have
    grps = array( [i+1 for i in range(len(neighs)-1) if (neighs[i+1][0]-neighs[i][0])>1e-5] + [len(neighs)], dtype=int)
    polyhedra_sizes = set([poly[1] for poly in cases] + [poly[1] for poly in cases_outside])
    print('poly_sizes=', polyhedra_sizes, 'group size=', grps, file=log)
    while len(grps) and not (set(grps) & polyhedra_sizes):
        # all neighbor groups are of such sizes that we can not try any of our polyhedra. So the algorithm will fail.
        # In this case we just remove anough groups at the beginning, and try to get something with the second or third set of neighbors.
        print('WARN: The group of neighbors with groups=',grps,' are not compatible with any polyhedra. We are eliminating first group and are left with ', (grps[1:]-grps[0]), file=log)
        neighs = neighs[grps[0]:] # We take out the first group of neighbors
        grps = grps[1:]-grps[0] # We take out the first group
    if len(grps) == 0: # Could not find any shape
        print('Could not detect the type of environment, Boiling out.', file=log)
        return None,None
    
    tobohr = 1/0.5291772083
    
    dist0 = [n[0]*tobohr for n in neighs]
    rij_fract = array([n[2] for n in neighs])
    names = [n[1] for n in neighs]
    # cartesian=fractional*M
    # [fa,fb,fc]*[[ax,ay,az],
    #             [bx,by,bz],
    #             [cx,cy,cz]]
    rij = rij_fract @ matrix  
    dist = linalg.norm(rij,axis=1)    # precompute all distances, which should match dist0
    rij_norm = (rij.T * 1/dist).T     # precompute all unit vectors

    # Compute all angles between all unit vectors vectors
    phi = arccos( clip( rij_norm @ rij_norm.T, -1, 1) )  # we remove possible values beyond 1, so that arccos works
    w = exp(1-(dist/dist[0])**6)
    _lav_ = dist * w
    print('phi=', file=log)
    N = len(phi)
    for i in range(N):
        print(('{:6.1f} '*N).format(*(phi[i,:]*180/pi)), file=log)

    
    criteria=[]
    # Now go over all atoms in the cage, and guess what type of cage it is
    for (ctype,N,angles) in cases: # over possible environments: cube, octahedra, tetrahedra
        if N>len(rij): continue
        if len(dist)>N and abs(dist[N-1]-dist[N])<1e-4: continue # we should not cut neighbors of the next one is equally far away.
        
        # First, compute effective coordination number (see VESTA manual)
        lav = sum(_lav_[:N])/sum(w[:N])
        ws = sum(exp(1-(dist[:N]/lav)**6))
        # Second, compare all angles to those expected
        angs = array([sorted(phi[i,:N])[1:] for i in range(N)])    # remove the smallest angle, because it is zero on diagonal
        chi1 = abs(ws-N)/2.
        
        if ctype=='square-piramid':
            var0 = np.var(angs, axis=1)
            itop = argmin(var0) # index of the top atom
            base = list(range(N))
            base.remove(itop)   # index to the points in the base of piramid
            r_base = rij[base]  # these vectors should be in the same plane. Let's check
            r_inplane = r_base[1:]-r_base[0]
            vol = abs(cross(r_inplane[0],r_inplane[1]) @ r_inplane[2])/6. # vol should be zero
            top_to_base = mean(angs[itop]) 
            cs = cos(top_to_base)    # 
            two_base = arccos(cs**2) # 
            one_base = arccos(2*cs**2-1) # one_base = 2*top_to_base
            if one_base>top_to_base:
                angles = array([two_base,two_base,top_to_base,one_base])
            else:
                angles = array([two_base,two_base,one_base,top_to_base])
            sm=sum((angs[base]-angles)**2)
            #print('Vol in square piramid=', vol, file=log)
            chi2 = (sm/(N-1)**2 + var0[itop])*2 + vol/3
        elif ctype=='triangular-bipyramid':
            max_ang = [max(angs[i]) for i in range(N)]
            bipyramid_ind = sorted(range(N), key=lambda i: -max_ang[i])
            top_bot = bipyramid_ind[0:2]
            base = bipyramid_ind[2:N]
            angles = [pi/2]*3 + [pi]
            diff = sum((angs[top_bot]-array(angles))**2)
            #print('angs[top_bot]=', angs[top_bot]*180/pi, 'diff=', diff, file=log)
            chi2 = sqrt(sum((angs[top_bot]-array(angles))**2)/(N*(N-1.)))*2
            angles = [pi/2]*2 + [2*pi/3]*2            
            diff = sum((angs[base]-array(angles))**2)
            #print('angs[base]=', angs[base]*180/pi, 'diff=', diff, file=log)
            chi2 += sqrt(sum((angs[base]-array(angles))**2)/(N*(N-1.)))*2
            #print('chi2=', chi2)
        else:
            chi2 = sqrt( sum((angs-angles)**2)/angs.size )*2
        
        #print('chi2(angle)=', chi2, file=log)
        criteria.append( (ctype, chi1, chi2 ) ) # now remember how good is this guess in terms of coordination number and angle variance
        print('trying {:15s}: accuracy {:6.3f} <w>-N={:6.3f} <phi>-<phi0>={:6.3f}'.format(ctype,chi1**2+chi2**2,chi1,chi2), file=log)

    # It has n neighbors if we find n such atoms, or, first n atoms have certain distance, and n+1 atom has larger distance
    have_n_neigbors = [len(neighs)==n or (len(neighs)>n and neighs[n][0]-neighs[n-1][0]>1e-4) for n in range(16)]
    _cases_ = [case for case in cases_outside if have_n_neigbors[case[1]]]
    for (ctype,N,angles) in _cases_:
        rjk0 = rij[1:]-rij[0]
        normal_plane = cross(rjk0[0],rjk0[1]) # normal to the plane formed by rij[1]-rij[0] and rij[2]-rij[0]
        normal_plane *= 1/linalg.norm(normal_plane)
        #print('normal_plane=', normal_plane, file=log)
        voln=0
        for l in range(2,N-1):
            voln += abs(rjk0[l] @ normal_plane)/6.
        #print('voln[N='+str(N)+']=', voln, 'normal_plane=', normal_plane, file=log)
        if voln>1e-4: break  # these points are not on the same plane
        # check if all equal distances.
        chi1=0
        for j in range(N):
            for k in range(j+1,N):
                chi1 += (dist[j]-dist[k])**2
        chi1 *= 1/(N*(N-1)/2)
        chi2=0
        if N==3:
            #Vol = abs(dot(cross(rij[0],rij[1]),rij[2]))/6.
            #chi1 += 2*exp(-Vol*20.)  # we don't want volume to be small, because than it is different shape
            chi2 = (phi[0,1]-phi[0,2])**2
            n_down = rij[0]+rij[1]+rij[2]
            n_down /= linalg.norm(n_down)
            n_plane = cross(rij[1]-rij[0],rij[2]-rij[0])
            n_plane /= linalg.norm(n_plane)
            is_in_center = (1-abs(dot(n_down,n_plane)))**2 * 3.
            #print('rij=', rij[:3,:], file=log)
            #print('    is_in_center of peak-of-tetrahedron=', is_in_center, 'n_down=', n_down, 'n_plane=', n_plane, file=log)
            chi2 += 0.26+is_in_center # so that we take tetrahedra when chi2 of tetrahedra is 0.08 or smaller
        if N==4:
            n_down = rij[0]+rij[1]+rij[2]+rij[3]
            #print('n_down=', n_down, file=log)
            if linalg.norm(n_down)>1e-6:
                n_down /= linalg.norm(n_down)
                in_plane_=[]
                in_plane_.append( cross(rij[1]-rij[0],rij[2]-rij[0]) )
                in_plane_.append( cross(rij[1]-rij[0],rij[3]-rij[0]) )
                in_plane_.append( cross(rij[2]-rij[0],rij[3]-rij[0]) )
                i_which = argmax([linalg.norm(in_plane_[i]) for i in range(3)])
                in_plane = in_plane_[i_which]
                in_plane /= linalg.norm(in_plane)
                is_in_center = (1-abs(dot(n_down,in_plane)))**2 * 3.
                chi2 += is_in_center
                #print('is_in_center=', is_in_center, file=log)
            chi2 += 0.2 # so that we take octahedra is chi2 of octahedar is smaller than 0.05
        criteria.append( (ctype, chi1, chi2 ) ) # now remember how good is this guess in terms of coordination number and angle variance
        print('trying {:15s}: accuracy {:6.3f} <w>-N={:6.3f} <phi>-<phi0>={:6.3f}'.format(ctype,chi1**2+chi2**2,chi1,chi2), file=log)
    
    # On the basis of these criteria, we should be able to tell which environment we have
    criteria = sorted(criteria, key=lambda ct: ct[1]**2 + ct[2]**2 ) # sort according to the  (angle-variance*2*pi)^2+(bond-variance)^2
    #best_so_far = criteria[0][1]**2+criteria[0][2]**2
    best_so_far = criteria[0][2]
    cutoff_for_best = 0.2
    #print('best_so_far=', criteria[0][2],file=log)
    have_four_neighbors = have_n_neigbors[4] # len(neighs) == 4 or (len(neighs)>4 and neighs[4][0]-neighs[3][0]>1e-4)
    if best_so_far>cutoff_for_best and have_four_neighbors:
        # It is likely not a regular tetrahedron, hence checking for planar quadrilateral
        N=4
        dr = rij[1:N]-rij[0]
        #print('dr=', dr)
        # volome^2 of the tetrahedron composed of points [r0,r1,r2,r3]
        V2 = abs(cross(dr[0], dr[1]) @ dr[2])/6.
        # volome^2 of the tetrahedron composed of points [center,r0,r1,r2]
        V3 = abs(cross(rij[0],rij[1]) @ rij[2])/6.
        if (V2+V3<best_so_far):
            criteria.insert(0, ('planar quadrilateral', V2, V3) )
    have_five_neigbors = have_n_neigbors[5] # len(neighs) == 5 or (len(neighs)>5 and neighs[5][0]-neighs[4][0]>1e-4)
    if best_so_far>cutoff_for_best and have_five_neigbors:
        # trying square piramid in which origin is in the basal plane of the piramid
        N=5
        min_angle_diff = 0.5
        j1=-1
        for i in range(N):
            phi[i,i]=pi/2
            all_ninety = sum(abs(phi[i,:]-pi/2))
            phi[i,i]=0
            if all_ninety < min_angle_diff:
                j1 = i
                break
        if j1>=0:
            chi2=all_ninety
            #print('chi2=', chi2, 'j1=', j1)
            for i in range(N):
                if i==j1: continue
                phi[i,i]=pi/2
                i_max = argmax(phi[i,:])
                irest = list(range(i))+list(range(i+1,N))
                irest.remove(i_max)
                dchi = abs(phi[i,i_max]-pi) + sum([abs(phi[i,j]-pi/2) for j in irest])
                chi2 += abs(phi[i,i_max]-pi) + sum([abs(phi[i,j]-pi/2) for j in irest])
                #print('dchi=', dchi, 'chi2=', chi2, 'im=', i_max, 'irest=', irest)
                phi[i,i]=0
            chi2 *= 1./(N*(N-1))
            if (chi2<max_angle_variance and chi2 < best_so_far):
                criteria.insert(0, ('square prism with base', chi2, 0) )
                nz = rij_norm[j1]
                iwhich = list(range(j1))+list(range(j1+1,N))
                #print('nz=', nz, 'iwhich=', iwhich)
                nis=[]
                for ii,i in enumerate(iwhich):
                    angs=[]
                    _nis_=[]
                    for j in iwhich:  # find one with 90 degree, i.e., max scalar product
                        if i==j: continue
                        angs.append( rij_norm[i] @ rij_norm[j] )
                        ni = rij[i]+rij[j]           # new unit vectors are between each pair of vertices
                        ni *=  1./linalg.norm(ni)    # and normalize
                        _nis_.append(ni)
                        #print('i,j=', i,j, 'angs=', angs[-1], 'ni=', ni)
                    icolin = argmin(angs)
                    del _nis_[icolin]
                    nis.extend(_nis_)
                    #print('nis=', _nis_)
                #print('nis=', nis)
                R0=[]
                for i in range(2):
                    nu = nis.pop(0) # first unit vector
                    which = argmin([dot(nu,ni) for ni in nis] ) # which one has most negative dot product is the opposite
                    nv1 = nis.pop(which)  # nv is opposite to nu
                    which = argmin([dot(nu,ni) for ni in nis] ) # which one has most negative dot product is the opposite
                    nv2 = nis.pop(which)
                    which = argmax([dot(nu,ni) for ni in nis] ) # which one has most negative dot product is the opposite
                    nv3 = nis.pop(which)
                    n1 = nu+nv3-(nv1+nv2)           # now take the average of the two, which are opposite
                    n1 *= 1./linalg.norm(n1) # and normalize
                    R0.append(n1)       # we have a good unit vector
                    #print('nu=', nu, 'nv1=', nv1, 'nv2=', nv2, 'nv3=', nv3, 'n1=', n1)
                R0.append(nz)
                R0_ready = R0

    #best_so_far = criteria[0][2]
    have_six_neighbors = have_n_neigbors[6] # len(neighs)==6 or (len(neighs)>6 and (neighs[6][0]-neighs[5][0]>1e-3))
    if best_so_far>cutoff_for_best and have_six_neighbors:
        # likely not octahedra, hence checking for trigonal prism
        # Could not find something in literature, hence long elaborate algorithm.
        N=6
        rjk = zeros((N,N,3)) # vectors between all neighbors
        djk = zeros((N,N))   # and their distances
        for j in range(N):
            for k in range(j+1,N):
                rjk[j,k] = rij[j]-rij[k]
                rjk[k,j] = rij[k]-rij[j]
                djk[j,k] = djk[k,j] = linalg.norm(rjk[j,k])
        #phi = arccos( clip( rij_norm @ rij_norm.T, -1, 1) ) # all angles again
        ph0 = copy(phi[0,:N]) # we assume the first atom (index 0) is in the triangle of trigonal prism.
        ph0[0]=3*pi          # for searching for the minimum we eliminate first atom.
        (j1,j2) = Find_two_closest_smallest_values(ph0, range(1,N)) # Now finding two other neighbors in triangle with 0.
        #print('ph['+str(j1)+']=',ph0[j1]*180/pi,'ph['+str(j2)+']=',ph0[j2]*180/pi)
        n1 = (rij[0]+rij[j1]+rij[j2]) # unit vector from the origin to the center of the triangle. It is along the trigonal prism.
        n1 *= 1./linalg.norm(n1)
        #print('n1=', n1)
        triangle_1 = rij[j1]-rij[0]
        triangle_2 = rij[j2]-rij[0]
        # Now we search for the neighbors on the opposite triangle.
        # The first triangle has neigbors [0,j1,j2]
        # On the opposite end we will have [j3,j4,j5]
        # hence vectors from j3 to 0 and j4 to j1 and j5 to j2 are along n1.
        in_direction_n1 = zeros(N)
        for j in range(1,N):
            if (j==j1 or j==j2): continue # searching for index not in first triangle
            n03 = rjk[0,j]/djk[0,j]  #== (rij[0]-rij[j])/linalg.norm(rij[0]-rij[j])
            in_direction_n1[j] = dot(n03,n1)  # is vector from j to 0 in direction of n1?
            #print('j=', j, 'n03=', n03, 'in_dir=', in_direction_n1[j])
        #print('in_direction_n1=', in_direction_n1, file=log)
        j3 = argmax(in_direction_n1) # find the best one, for which vector 0->j3 is along n1.
        #print('j3=', j3)
        in_direction_n1 = zeros(N)
        for j in range(1,N): # searching for neighbor j4, which is on opposite side of j1.
            if (j==j1 or j==j2 or j==j3): continue
            n04 = rjk[j1,j]/djk[j1,j] #(rij[j1]-rij[j])/linalg.norm(rij[j1]-rij[j])
            in_direction_n1[j] = dot(n04,n1) # vector from j->j1
        j4 = argmax(in_direction_n1) # and best vector parallel no n1 gives j4
        #print('j4=', j4, 'in_dir=', in_direction_n1)
        # The last remaining index is j5, and corresponds to the second triangle.
        j5 = [x for x in range(1,6) if x not in [j1,j2,j3,j4]][0]

        # The x will be choosen along the trigonal prism. It is best approximated
        # by the direction to both triangles formed by [0,j1,j2] and [j3,j4,j5]
        nx = (rij[0]+rij[j1]+rij[j2])-(rij[j3]+rij[j4]+rij[j5])
        nx *= 1/linalg.norm(nx)
        #print('nx=', nx, 'js=', [0,j1,j2,j3,j4,j5], rij[[0,j1,j2,j3,j4,j5],:], file=log)

        prism_vol = abs(cross(triangle_1,triangle_2) @ nx)
        #print('prism_vol=', prism_vol, file=log)
        # We compute the chi^2 for trigonal prism by checking how parallel are
        # the intersections of three planes, i.e., edges of trigonal prism.
        # These are the vectors 0->j3, j1->j4, and j2->j5. They should all be along nx,
        # as computed above. 
        paral_0_j3 = dot(nx,rjk[0,j3]/djk[0,j3])
        paral_j1_j4 = dot(nx,rjk[j1,j4]/djk[j1,j4])
        paral_j2_j5 = dot(nx,rjk[j2,j5]/djk[j2,j5])
        #print('paral=', paral_0_j3, paral_j1_j4, paral_j2_j5, 'vol=', prism_vol, file=log)
        
        chi2 = (abs(paral_0_j3-1)+abs(paral_j1_j4-1)+abs(paral_j2_j5-1))*1.5
        chi2r = 2*exp(-prism_vol*20) # volume should not be zero
        # We accept trigonal prism with triangle which is not necessary equilateral, but can have
        # two sides with equal length, hence only two angles are equal. We find the minimal
        # difference in angles, and add this is penalty.
        chi3 = min(abs(phi[0,j1]-phi[0,j2]),abs(phi[0,j1]-phi[j1,j2]),abs(phi[0,j2]-phi[j1,j2]))
        chi3 +=min(abs(phi[j3,j4]-phi[j3,j5]),abs(phi[j3,j4]-phi[j4,j5]),abs(phi[j3,j5]-phi[j4,j5]))
        chi3 *= 2*pi
        #print('chi2=', chi2, 'chi3=', chi3, 'chi2r=', chi2r, file=log)
        chi2 += chi2r
        print('trying {:15s}: accuracy {:6.3f} <w>-N={:6.3f} <phi>-<phi0>={:6.3f}'.format('trigonal prism',chi2**2+chi3**2,chi2,chi3), file=log)
        if (chi2<max_bond_variance and chi3<max_angle_variance and chi2<best_so_far):
            # Since we already have all ingredients here, we just compute where unit vectors should point.
            # nx is along the trigonal prism. The xy plane is perpendicular to the triangles and
            # we take the largest side of the triangle to lay in the xy plane.
            cases = (djk[0,j1],djk[j1,j2],djk[0,j2])
            if abs(cases[0]-cases[1])<1e-3 and abs(cases[0]-cases[2])<1e-3: # degenerate case
                ny1 = rjk[j1,0]+rjk[j4,j3]
                ny1 *= 1/linalg.norm(ny1)
                nz1 = -(rij[0]+rij[j1]+rij[j3]+rij[j4])
                #print('ny1=', ny1, 'nz1=', nz1, file=log)
                lnz1 = linalg.norm(nz1)
                if lnz1<1e-4:
                    nz1 = cross(nx, ny1)
                    print('   It seems tensegrity-prism, choosing nz1=', nz1, file=log)
                else:
                    nz1 /= lnz1
                    if linalg.norm(cross(nz1,ny1))<1e-4:
                        nz1 = cross(nx, ny1)
                        print('   It seems tensegrity-prism, choosing nz1=', nz1, file=log)
                ny2 = rjk[j1,j2]+rjk[j4,j5]
                ny2 *= 1/linalg.norm(ny2)
                nz2 = -(rij[j1]+rij[j2]+rij[j4]+rij[j5])
                #print('ny2=', ny2, 'nz2=', nz2, file=log)
                lnz2 = linalg.norm(nz2)
                if lnz2<1e-4:
                    nz2 = cross(nx, ny2)
                    print('   It seems tensegrity-prism, choosing nz2=', nz2, file=log)
                else:
                    nz2 /= lnz2
                    if linalg.norm(cross(nz2,ny2))<1e-4:
                        nz2 = cross(nx, ny2)
                        print('   It seems tensegrity-prism, choosing nz2=', nz2, file=log)
                ny3 = rjk[j2,0]+rjk[j5,j3]
                ny3 *= 1/linalg.norm(ny3)
                nz3 = -(rij[0]+rij[j2]+rij[j3]+rij[j5])
                #print('ny3=', ny3, 'nz3=', nz3, file=log)
                lnz3 = linalg.norm(nz3)
                if lnz3<1e-4:
                    nz3 = cross(nx,ny3)
                    print('   It seems tensegrity-prism, choosing nz3=', nz3, file=log)
                else:
                    nz3 /= lnz3
                    if linalg.norm(cross(nz3,ny3))<1e-4:
                        nz3 = cross(nx, ny3)
                        print('   It seems tensegrity-prism, choosing nz3=', nz3, file=log)
                #print('nx=', nx, 'ny1=', ny1, 'nz1=', nz1, file=log)
                #print('nx=', nx, 'ny1=', ny2, 'nz1=', nz2, file=log)
                #print('nx=', nx, 'ny1=', ny3, 'nz1=', nz3, file=log)
                
                Rs = [ResortToDiagonal(array([nx,ny1,nz1])),
                      ResortToDiagonal(array([nx,ny2,nz2])),
                      ResortToDiagonal(array([nx,ny3,nz3]))]
                I = identity(3)
                iwhich = argmin( [sum(abs(Rs[0]-I)), sum(abs(Rs[1]-I)), sum(abs(Rs[2]-I))] )
                R0_ready = Rs[iwhich]
                #print('Rs[0]=', Rs[0], 'Rs[1]=', Rs[1], 'Rs[2]=', Rs[2], 'choosing', iwhich)
            else:
                iwhich = argmin(cases)
                if iwhich==0:
                    nz1 = -(rij[0]+rij[j1]+rij[j3]+rij[j4])
                    nz1 *= 1/linalg.norm(nz1)
                    ny1 = rjk[j1,0]+rjk[j4,j3]
                    ny1 *= 1/linalg.norm(ny1)
                    R0_ready = [nx,ny1,nz1]
                elif iwhich==1:
                    nz2 = -(rij[j1]+rij[j2]+rij[j4]+rij[j5])
                    nz2 *= 1/linalg.norm(nz2)
                    ny2 = rjk[j1,j2]+rjk[j4,j5]
                    ny2 *= 1/linalg.norm(ny2)
                    R0_ready = [nx,ny2,nz2]
                else:
                    nz3 = -(rij[0]+rij[j2]+rij[j3]+rij[j5])
                    nz3 *= 1/linalg.norm(nz3)
                    ny3 = rjk[j2,0]+rjk[j5,j3]
                    ny3 *= 1/linalg.norm(ny3)
                    R0_ready = [nx,ny3,nz3]
            #print('R0_ready=', R0_ready, file=log)
            criteria.insert(0, ('trigonal prism', chi2, chi3) )
        
    # First take the cage, which has most similar angles.
    #print('We are trying', criteria[0][0], 'with c1=', criteria[0][1],'<', max_bond_variance,
    #        'and c2=', criteria[0][2],'<',max_angle_variance, file=log)
    if criteria[0][1]<max_bond_variance and criteria[0][2]<max_angle_variance: # If the bond variance is not too bad, i.e., |<w>-N|< 0.1, we found the type of cage.
        ctype = criteria[0]
    elif len(criteria)>1 and criteria[1][1]<max_bond_variance and criteria[1][2]<max_angle_variance: # If the bond variance was very bad for the first case, we go down the list
        ctype = criteria[1]
    else: # If the second is not OK, we boil out at the moment
        print('Could not detect the type of environment, Boiling out.', file=log)
        return None,grps[0]
    print('Found the environment is {:s}  (<w-n>= {:8.3f},<phi>= {:8.3f})'.format(*ctype), file=log)
    
    if ctype[0]=='octahedra':
        N=6
        # Now that we know it is octahedra, we take all vectors to atoms, and construct the best coordinate system that
        # goes through these atoms
        nis = [rij_norm[i] for i in range(N)] # all unit vectors, but list
        R0=[]
        for i in range(3):
            nu = nis.pop(0)      # first unit vector
            which = argmin([dot(nu,ni) for ni in nis]) # which one has most negative dot product is the opposite
            nv = nis.pop(which)  # nv is opposite to nu
            n1 = (nu-nv)/linalg.norm(nu-nv) # now take the average of the two, which are opposite
            R0.append(n1)        # we have a good unit vector
            #print( n1, nu, nv, arccos(cosp[which])*180/pi )
    elif ctype[0]=='tetrahedra':
        # Now that we know it is tetrahedra, we take vectors which go through the middle of the two atoms, and construct
        # the best coordinate system with such vectors.
        N=4
        n0 = rij[:N] # all unit vectors
        nis=[]
        for i in range(N):
            for j in range(i+1,N):
                ni = n0[i]+n0[j]           # new unit vectors are between each pair of vertices
                ni *=  1./linalg.norm(ni)  # and normalize
                nis.append(ni)
        R0=[]
        for i in range(3):
            nu = nis.pop(0) # first unit vector
            which = argmin([dot(nu,ni) for ni in nis] ) # which one has most negative dot product is the opposite
            nv = nis.pop(which)  # nv is opposite to nu
            n1 = nu-nv           # now take the average of the two, which are opposite
            n1 *= 1./linalg.norm(n1) # and normalize
            R0.append(n1)       # we have a good unit vector
            #print n1, nu, nv, arccos(cosp[which])*180/pi
    elif ctype[0]=='cube':
        N=8
        # For each vector to atom, we should find 3 other atoms which have smallest angle between them
        # Now we have four vectors with small angle, which define the phase of a cube.
        # We then sum this four vectors, and get unit vector along x,y, or z
        def get_sorted_index(phi0):
            return sorted(range(len(phi)),key=lambda i:phi0[i])
        indx = get_sorted_index(phi[0])
        j1,j2,j3 = indx[1:4]
        j4,j5,j6 = indx[4:7]
        j7 = indx[7]
        n1 = rij[j1]-rij[0]
        n2 = rij[j2]-rij[0]
        n3 = rij[j3]-rij[0]
        nr = [rij[j4]-rij[j7], rij[j5]-rij[j7], rij[j6]-rij[j7]]
        k1 = argmin([nr[i]@n1 for i in range(3)])
        k2 = argmin([nr[i]@n2 for i in range(3)])
        k3 = argmin([nr[i]@n3 for i in range(3)])
        nx = n1-nr[k1]
        nx *= 1/linalg.norm(nx)
        ny = n2-nr[k2]
        ny *= 1/linalg.norm(ny)
        nz = n3-nr[k3]
        nz *= 1/linalg.norm(nz)
        R0 = [nx,ny,nz]
        R0 = [nx,ny,nz]
    elif ctype[0]=='planar quadrilateral':
        N=4
        n0 = rij[:N]
        olap=zeros((N,N))
        for i in range(N):
            for j in range(i,N):
                olap[i,j] = n0[i] @ n0[j]
        i_opposite = argmin(olap[0,:])
        nx = n0[0]-n0[i_opposite]
        nx *= 1/linalg.norm(nx)
        remain = [i for i in range(1,4) if i!=i_opposite]
        j_opposite = argmin(olap[remain[0],:])
        ny = n0[remain[0]]-n0[j_opposite]
        ny *= 1/linalg.norm(ny)
        R0 = [nx,ny,cross(nx,ny)]
        #
        ##n0 = rij_norm[:N] # all unit vectors
        ## We again take vectors which go through the middle of the two atoms, and construct
        ## the best coordinate system with such vectors.
        #nis=[]
        #for i in range(N):
        #    for j in range(i+1,N):
        #        ni = n0[i]+n0[j]           # new unit vectors are between each pair of vertices
        #        d = linalg.norm(ni)
        #        if d>1e-10:
        #            ni *=  1./linalg.norm(ni)  # and normalize
        #        nis.append(ni)
        #nrm = linalg.norm(nis,axis=1)
        #indx = sorted(range(len(nis)), key=lambda i:-nrm[i])
        ## the smallest two nis are essentially zero, because the system is complanar. We throw away the two zero vectors.
        #nis_ = [nis[indx[i]] for i in range(0,4)]
        #
        #print('rij=', n0, file=log)
        #print('nis_=', nis_, file=log)
        #
        #R0=[]
        #for i in range(2):
        #    nu = nis_.pop(0) # first unit vector
        #    olap = [dot(nu,ni) for ni in nis_]
        #    which = argmin(olap) # which one has most negative dot product is the opposite
        #    nv = nis_.pop(which)  # nv is opposite to nu
        #    n1 = nu-nv           # now take the average of the two, which are opposite
        #    n1 *= 1./linalg.norm(n1) # and normalize
        #    print('')
        #    print('nu=', nu, 'olap=', olap, 'which=', which, 'nv=', nv, 'n1=', n1, file=log)
        #    print('nis_=', nis_, file=log)
        #    R0.append(n1)       # we have a good unit vector
        #R0.append(cross(R0[0],R0[1])) # because is planar, the z should be out of the plane
    elif ctype[0]=='peak-of-tetrahedron':
        N=3
        vol = cross(rij[0],rij[1]) @ rij[2]
        n_base = cross(rij[1]-rij[0], rij[2]-rij[0])
        n_base *= 1./linalg.norm(n_base)
        hv = (rij[0]+rij[1]+rij[2])/3.
        h = abs(hv @ n_base)              # height
        ra = (dist[0]+dist[1]+dist[2])/3. # average distance
        #print('vol=', vol, 'h/(ra/3)=', h/(ra/3), file=log)
        if vol>10 and (h > 0.5*ra/3 and h < 2*ra/3):
            # if tetrahedron has sufficient volume so that it is not just in-plane
            R0= [rij_norm[0], rij_norm[1], rij_norm[2]]
        else:
            #nz = -hv
            nz_length = linalg.norm(hv)
            if (nz_length > 1e-5):
                nz = -hv/nz_length
                n1 = rij[0]-(rij[0]@nz)*nz
                n1 *= 1./linalg.norm(n1)
                n2 = rij[1]+rij[2]-((rij[1]+rij[2])@nz)*nz
                n2 *= 1./linalg.norm(n2)
                nx = (n1-n2)/linalg.norm(n1-n2)
                ny = cross(nz,nx)
            else:
                n1 = rij_norm[0]
                n2 = rij[1]+rij[2]
                n2 *= 1/linalg.norm(n2)
                nx = (n1-n2)/linalg.norm(n1-n2)
                ny = rij[1]-rij[2]
                ny = ny - (ny@nx)*nx
                ny *= 1/linalg.norm(ny)
                nz = cross(nx,ny)
            R0=[nx,ny,nz]
    elif ctype[0]=='peak-of-square-piramid':
        N=4
        djk0 = linalg.norm(rij[1:N]-rij[0],axis=1)
        i_diagonal = argmax(djk0)+1               # this one is on diagonal of the square
        which = [i for i in range(1,4) if i!=i_diagonal] # these two are nearest neighbors in a square
        n1 = rij[which[0]]-rij[0] # 
        n2 = rij[which[1]]-rij[0] #
        n3 = rij[which[0]]-rij[i_diagonal]
        n4 = rij[which[1]]-rij[i_diagonal]
        nx = n1-n4
        nx *= 1/linalg.norm(nx)
        ny = n2-n3
        ny *= 1/linalg.norm(ny)

        vol = abs(cross(rij[0],rij[which[0]]) @ rij[i_diagonal])+(cross(rij[0],rij[which[1]]) @ rij[i_diagonal])
        #print('vol2=', vol, file=log)
        #print('i_diagonal=', i_diagonal, 'which=', which, file=log)
        #print('n1=', n1, 'n2=', n2, 'n3=', n3, 'n4=', n4, file=log)
        #print('nx=', nx, 'ny=', ny, file=log)
        if (vol>0.5):
            nz = -(rij[0]+rij[1]+rij[2]+rij[3])/4.
            nz_length = linalg.norm(nz)
            if (nz_length > 1e-10):
                nz *= 1./nz_length
            #print('nz=', nz, file=log)
        else:
            # if the piramid is completely squashed, we can not
            nz = cross(nx,ny)
        R0=[(nx+ny)/sqrt(2),(nx-ny)/sqrt(2),nz]
    elif ctype[0]=='peak-of-hexagonal-piramid':
        N=6
        ix = argmax(rij_norm[:6,0])
        nx = rij_norm[ix]
        pr = rij_norm[:6,:] @ nx
        indx = sorted(range(6), key=lambda i: pr[i])
        nx_better = rij[ix] - rij[indx[0]]
        nx_better *= 1/linalg.norm(nx_better)
        if linalg.norm(rij[1]-rij[3])<linalg.norm(rij[1]-rij[4]):
            ny = rij[indx[1]]+rij[indx[3]] - (rij[indx[2]]+rij[indx[4]])
        else:
            ny = rij[indx[1]]+rij[indx[4]] - (rij[indx[2]]+rij[indx[3]])
        ny *= 1/linalg.norm(ny)
        nz = cross(nx,ny)
        R0 = [nx,ny,nz]
    elif ctype[0]=='cuboctahedron':
        N=12
        djk = zeros((N,N))   # and their distances
        for j in range(N):
            for k in range(j+1,N):
                djk[j,k] = djk[k,j] = linalg.norm(rij[j]-rij[k])
        indx = sorted(range(N),key=lambda i:djk[0,i])
        # nearest neighbors are indx[1:5], and next-nn are indx[5],indx[6]
        n1 = rij[indx[5]]-rij[0] # next-nearest neighbor is in one of orthogonal directions
        n1 *= 1/linalg.norm(n1)
        n2 = rij[indx[6]]-rij[0] # next-nearest neighbor is in one of orthogonal directions
        n2 *= 1/linalg.norm(n2)
        # Now we want to get four vertices of the square to improve the precision
        # of the unit vector in case of distortions.
        in_plane1 = abs((rij[indx[1:5]]-rij[0]) @ n1)
        two_largest = heapq.nlargest(2, range(4), key=lambda i:in_plane1[i])
        indx1 = [indx[i+1] for i in two_largest]
        # four points of the square give better approximation
        n1_better = rij[0]+rij[indx[5]]+rij[indx1[0]]+rij[indx1[1]]
        n1_better *= 1/linalg.norm(n1_better)
        # Now the same for the other direction
        in_plane2 = abs((rij[indx[1:5]]-rij[0]) @ n2)
        two_largest = heapq.nlargest(2, range(4), key=lambda i:in_plane2[i])
        indx2 = [indx[i+1] for i in two_largest]
        # and four points of the square for better approximation
        n2_better = rij[0]+rij[indx[6]]+rij[indx2[0]]+rij[indx2[1]]
        n2_better *= 1/linalg.norm(n2_better)
        # the third direction is of course cross produce
        n3 = cross(n1_better,n2_better)
        # but we want to have a better approximation for distorted object
        # nearest neighbors of rij[0] and vectors from them to rij[0]
        vc3 = [rij[indx[i]]-rij[0] for i in range(1,5)]
        # two of these nearest neigbors should be on the square in positive n3 direction
        pos_square = [indx[i+1] for i in range(4) if ( vc3[i] @ n3 > 0) ]
        # and the other two on the square in the opposite n3 direction
        #neg_square = [indx[i+1] for i in range(4) if ( vc3[i] @ n3 < 0) ]
        # now we have two out of four points on each of the two squares
        pos_square2=[] # finding the rest two points for positive square
        for i in range(2):
            djkp0 = linalg.norm(rij - rij[pos_square[i]], axis=1)
            indxz = sorted(range(N), key=lambda i:djkp0[i])
            poss = set(indxz[1:5])-set([0,*pos_square,indx[5],indx[6]])
            pos_square2.append(poss.pop())
        pos_square.extend(pos_square2)
        # Finally pos_square should contain all four points of a square
        n3_better = rij[pos_square[0]]+rij[pos_square[1]]+rij[pos_square[2]]+rij[pos_square[3]]
        n3_better *= 1/linalg.norm(n3_better)
        R0 = [n1_better, n2_better, n3_better]
    elif ctype[0]=='truncated tetrahedron':
        N=12
        #rt = array([[3, 1, 1],[1, 3, 1],[1, 1, 3],[-3, -1, 1],[-1, -3, 1],[-1, -1, 3],
        #        [-3, 1, -1],[-1, 3, -1],[-1, 1, -3],[3, -1, -1],[1, -3, -1],[1, -1, -3]])
        print('WARNING:ERROR: truncated tetrahedron not yet implemented! You should work on this....', file=log)
        R0=identity(3)
    elif ctype[0]=='trigonal prism' or ctype[0]=='square prism with base':
        R0 = R0_ready
        N=6 if ctype[0]=='trigonal prism' else 5
    elif ctype[0]=='square-piramid':
        N=5
        var0 = np.var(angs, axis=1)
        itop = argmin(var0) # index of the top atom
        base = list(range(N))
        base.remove(itop)   # index to the points in the base of piramid
        r_base = rij[base]  # these vectors should be in the same plane.
        n2 = sum(r_base,axis=0)
        n2 *= 1/linalg.norm(n2)
        n1 = rij_norm[itop]
        nz = n1-n2 if dot(n1,n2)<0 else n1+n2
        nz *= 1/linalg.norm(nz)
        j1 = base[0]            # first atom in the base
        phi[j1,j1]=2*pi
        phi[j1,itop]=2*pi
        j2 = argmin(phi[j1,:5]) # j2 is closest neigbor of j1 in the base
        j3,j4 = list( set(range(5))-set([itop,j1,j2]) )
        nx = rij[j1]+rij[j2]-(rij[j3]+rij[j4])
        nx -= (nx @ nz)* nz
        nx *= 1/linalg.norm(nx)
        ny = rij[j1]-rij[j2]
        if (ny @ (rij[j3]-rij[j4]))>0:
            ny += rij[j3]-rij[j4]
        else:
            ny -= rij[j3]-rij[j4]
        ny -= (ny @ nz)* nz
        ny *= 1/linalg.norm(ny)
        #print('n2=', n2, 'n1=', n1, 'nz=', nz, file=log)
        #print('base=', base, 'itop=', itop, 'j1=', j1, 'j2=', j2, file=log)
        #print('nx=', nx, 'ny=', ny, file=log)
        #R0 = [nx,ny,nz]
        # We decided to rotate for pi/2 around z axis
        R0 = [(nx+ny)/sqrt(2),(nx-ny)/sqrt(2),nz]
    elif ctype[0]=='triangular-bipyramid':
        N=5
        top_bot = bipyramid_ind[0:2]
        base = bipyramid_ind[2:N]
        #print('base=', base, 'top_bot=', top_bot, file=log)
        #print('r_base=', rij[base], file=log)
        #print('r_tb = ', rij[top_bot], file=log)
        nz = rij[top_bot[0]] - rij[top_bot[1]]
        nz *= 1/linalg.norm(nz)
        #print('nz=', nz, file=log)
        #imin = argmin(dist[base])
        imax = argmax(dist[base])
        imax = base[imax]
        #print('imax=', imax, file=log)
        ny = rij[imax]
        ny *= 1/linalg.norm(ny)
        #print('ny=', ny, file=log)
        base.remove(imax)
        nx = (rij[base[0]]-rij[base[1]])
        nx *= 1/linalg.norm(nx)
        #print('nx=', nx, file=log)
        R0 = [nx,ny,nz]
    else:
        print('Noy yet implemented', file=log)
        return None,grps[0]
    
    #for i in range(3): print( ('{:12.8f} '*3).format(*R0[i]))
    # Now orthogonalizing the set of vectors
    U,S,V = linalg.svd(R0)

    if (min(S)<0.3 or max(S)>2.):
        print('WARN: Since singular values are far from unity (S=',S, ') we decided that this polyhedron is a bad choice. Resorting to global coordinate axis.', file=log)
        return None,N
    
    R = dot(U,V)
    R = ResortToDiagonal(R)
    print('singular values=', S, file=log)
    return R,N

def ResortToDiagonal(R):
    # We now resort, so that the rotation is close to identity
    # This is not necessary, but is convenient
    permutations=[(0,1,2),(0,2,1),(1,0,2),(1,2,0),(2,1,0),(2,0,1)]
    ii, wi = 0, 1000
    for ip,p in enumerate(permutations):
        Rt = array( [R[p[0]], R[p[1]], R[p[2]]] )
        wj = sum(abs(abs(diag(Rt))-1.0))
        if wj<wi:
            ii=ip
            wi=wj
    p=permutations[ii]
    Rt = array( [R[p[0],:], R[p[1],:], R[p[2],:]] )
    for i in range(3):
        if Rt[i,i]<0 : Rt[i,:] *= -1
    Rt *= linalg.det(Rt)
    return Rt


if __name__ == '__main__':
    import glob, os
    import w2k_nn
    tobohr = 1/0.5291772083
    
    files = glob.glob('*.struct')
    if len(files) < 1:
        print('ERROR No struct file present.')
        sys.exit(1)
    # get case
    case, ext = os.path.splitext(os.path.basename(files[0]))
    # reading w2k output about the real space vectors a,b,c

    Using_w2k_nn = False
    if (Using_w2k_nn):
        with open(case+'.outputd', 'r') as g:
            for n,line in enumerate(g):
                if line[13:20]=='BR1_DIR':
                    break
            S2C=[]
            for i,line in enumerate(g):
                S2C.append( [float(x) for x in line.split()] )
                if i>=2: break
        S2C = array(S2C).T
        # reading w2k output about the neighbors
        with open(case+'.outputnn', 'r') as f:
            lines = f.readlines()
    else:
        # executes nn in python module to obtain case.outputnn_ and case.rotlm_
        w2k_nn.w2knn(case)
        g = open(case+'.outputnn_', 'r')
        with open(case+'.rotlm_', 'r') as fi:
            next(fi)
            BR1=[]
            for i in range(3):
                BR1.append( [float(x) for x in next(fi).split()] )
            BR1=array(BR1)
        S2C = linalg.inv(BR1).T*2*pi
        # reading w2k output about the neighbors
        with open(case+'.outputnn_', 'r') as f:
            lines = f.readlines()

    print(('w2k_conventional=\n'+('\n'.join(['{:9.5f} '*3+' ==a'+str(i+1) for i in range(3)]))).format(*ravel(S2C)))
    to_frac = linalg.inv(S2C)
    
    lines = [line.strip() for line in lines if (not line.startswith(' RMT(')) and (not line.startswith(' SUMS TO'))]
    headers = [(n,line) for n,line in enumerate(lines) if ('ATOM:' in line) and ('EQUIV.' in line)]
    
    print('Atoms found in ',case+'.outputnn')
    for i,(n,line) in enumerate(headers):
        dat = line.split()
        print('{:3d} {:3s}{:4s} {:8s} {:8s} {:8s}'.format(i+1,dat[1],dat[4],dat[6],dat[7],dat[8]))

    _iatom_ = input('Enter the atom for which you want to determine local axes: ')
    iatom = int(_iatom_)
    
    startline, text = headers[iatom-1]
    x,y,z = [float(coord) for coord in text.split()[-3:]]  # coordinates of central atom
    name = text[20:30].strip()
    print(name, x,y,z)

    d = 8 if Using_w2k_nn else 14
    p = [23+i*d for i in range(4)]
    p1 = [49,59] if Using_w2k_nn else [81,91]
    # Reads all atoms around the central atom
    neighbrs = [] # these are in lattice units
    for line in lines[startline+1:]:
        if not line: break # empty line stops the environment of this atom
        neigh = [float(line[p[0]:p[1]])-x, float(line[p[1]:p[2]])-y, float(line[p[2]:p[3]])-z]
        name = line[10:20].strip()
        dst=float(line[p1[0]:p1[1]])
        cneigh =  neigh @ S2C,neigh # convert coordinates into orthonormal basis
        neighbrs.append([dst,name,neigh,0])
    # Main part of the algorithm
    log = sys.stdout
    R,N = FindCageBasis(neighbrs, S2C, log)
    if R is not None:
        Rf = R @ to_frac
        print('Rotation to input into case.indmfl by locrot=-1 : ', file=log)
        print(file=log)
        for i in range(3): print( ('{:12.8f} '*3).format(*R[i,:]), file=log)
        print(file=log)
        print('Rotation in fractional coords : ', file=log)
        print(file=log)
        for i in range(3): print( ('{:12.8f} '*3).format(*Rf[i,:]*tobohr), file=log)
        print(file=log)
