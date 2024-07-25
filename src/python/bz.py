#!/usr/bin/env python
import numpy as np
from matplotlib.patches import Circle, PathPatch, Rectangle, FancyArrowPatch
import mpl_toolkits.mplot3d.art3d as art3d
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.proj3d import proj_transform
import sys
from numpy import *

from utils import W2kEnvironment
from wlatgen import Latgen
from wstruct import Struct

def rotation_matrix(d):
    """
    Calculates a rotation matrix given a vector d. The direction of d
    corresponds to the rotation axis. The length of d corresponds to 
    the sin of the angle of rotation.

    Variant of: http://mail.scipy.org/pipermail/numpy-discussion/2009-March/040806.html
    """
    sin_angle = np.linalg.norm(d)

    if sin_angle == 0:
        return np.identity(3)

    d /= sin_angle

    eye = np.eye(3)
    ddt = np.outer(d, d)
    skew = np.array([[    0,  d[2],  -d[1]],
                  [-d[2],     0,  d[0]],
                  [d[1], -d[0],    0]], dtype=np.float64)

    M = ddt + np.sqrt(1 - sin_angle**2) * (eye - ddt) + sin_angle * skew
    return M

def pathpatch_2d_to_3d(pathpatch, z = 0, normal = 'z'):
    """
    Transforms a 2D Patch to a 3D patch using the given normal vector.

    The patch is projected into they XY plane, rotated about the origin
    and finally translated by z.
    """
    if type(normal) is str: #Translate strings to normal vectors
        index = "xyz".index(normal)
        normal = np.roll((1.0,0,0), index)

    normal /= np.linalg.norm(normal) #Make sure the vector is normalised

    path = pathpatch.get_path() #Get the path and the associated transform
    trans = pathpatch.get_patch_transform()

    path = trans.transform_path(path) #Apply the transform

    pathpatch.__class__ = art3d.PathPatch3D #Change the class
    pathpatch._code3d = path.codes #Copy the codes
    pathpatch._facecolor3d = pathpatch.get_facecolor #Get the face color    

    verts = path.vertices #Get the vertices in 2D

    if normal[0]!=0:
        cn = normal @ array([0,sqrt(normal[0]**2+normal[1]**2),normal[2]])
        #cn = cos(60/180*pi)
        sn = sqrt(1-cn**2)
        M = array([[cn,-sn],[sn,cn]])
        verts = [ M @ array([x, y]) for x, y in verts]
        verts = array(verts)
        
    #print('verts=', verts)
    d = np.cross(normal, (0, 0, 1)) #Obtain the rotation vector    
    M = rotation_matrix(d) #Get the rotation matrix

    rs = array([[x, y, 0] for x, y in verts])
    rs = array([M @ rs[i] + array([0,0,z]) for i in range(len(rs))])
    
    #pathpatch._segment3d = array([ M @ array([x, y, 0]) + array([0, 0, z]) for x, y in verts])
    pathpatch._segment3d = rs

def pathpatch_translate(pathpatch, delta):
    """
    Translates the 3D pathpatch by the amount delta.
    """
    pathpatch._segment3d += delta

class Arrow3D(FancyArrowPatch):
    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)
    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)
    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return np.min(zs)
def _arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    '''Add an 3d arrow to an `Axes3D` instance.'''
    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)
setattr(Axes3D, 'arrow3D', _arrow3D)

def get_brillouin_zone_3d(cell):
    """
    Generate the Brillouin Zone of a given cell. The BZ is the Wigner-Seitz cell
    of the reciprocal lattice, which can be constructed by Voronoi decomposition
    to the reciprocal lattice.  A Voronoi diagram is a subdivision of the space
    into the nearest neighborhoods of a given set of points. 

    https://en.wikipedia.org/wiki/Wigner%E2%80%93Seitz_cell
    https://docs.scipy.org/doc/scipy/reference/tutorial/spatial.html#voronoi-diagrams
    """
    cell = np.asarray(cell, dtype=float)
    assert cell.shape == (3, 3)

    px, py, pz = np.tensordot(cell, np.mgrid[-1:2, -1:2, -1:2], axes=[0, 0])
    points = np.c_[px.ravel(), py.ravel(), pz.ravel()]

    from scipy.spatial import Voronoi
    vor = Voronoi(points)

    bz_facets = []
    bz_ridges = []
    bz_vertices = []

    # for rid in vor.ridge_vertices:
    #     if( np.all(np.array(rid) >= 0) ):
    #         bz_ridges.append(vor.vertices[np.r_[rid, [rid[0]]]])
    #         bz_facets.append(vor.vertices[rid])

    for pid, rid in zip(vor.ridge_points, vor.ridge_vertices):
        # WHY 13 ????
        # The Voronoi ridges/facets are perpendicular to the lines drawn between the
        # input points. The 14th input point is [0, 0, 0].
        if(pid[0] == 13 or pid[1] == 13):
            bz_ridges.append(vor.vertices[np.r_[rid, [rid[0]]]])
            bz_facets.append(vor.vertices[rid])
            bz_vertices += rid

    bz_vertices = list(set(bz_vertices))

    return vor.vertices[bz_vertices], bz_ridges, bz_facets



if __name__ == '__main__':
    log = sys.stdout
    w2k = W2kEnvironment()
    strc = Struct()
    strc.ReadStruct(w2k.case+'.struct', log)
    latgen = Latgen(strc, w2k.case, log)

    #k2cartes = linalg.inv(latgen.br2)
    #c2f = latgen.k2icartes @ linag.inv(latgen.br2) @ diag(latgen.pia)
    #c2f_inf = diag(1/latgen.pia) @ latgen.br2 @ linalg.inv(latgen.k2icartes)
    #print('c2f_inf=', c2f_inf, '')
    
    cell = latgen.br2.T
    # length of the basis vectors
    b1, b2, b3 = np.linalg.norm(cell, axis=1)

    v, e, f = get_brillouin_zone_3d(cell)
    
    fig = plt.figure(figsize=(3, 3), dpi=300)
    ax = plt.subplot(111, projection='3d')
    ax.view_init(elev=35, azim=45, roll=0)
    
    for ii,xx in enumerate(e):
        ax.plot(xx[:, 0], xx[:, 1], xx[:, 2], color='k', lw=1.0)
        print('bz_line'+str(ii))
        for i in range(len(xx)):
            #print('({:9.6f},{:9.6f},{:9.6f})'.format(xx[i,0]/latgen.pia[0],xx[i,1]/latgen.pia[1],xx[i,2]/latgen.pia[2]), end=',')
            print('  ({:13.10f},{:13.10f},{:13.10f})'.format(*(xx[i]/latgen.pia)), end='\n')

    # Draw a rectangle
    #p = Rectangle([-2/3.*b1,-2/3.*b2], 4/3.*b1, 4/3*b2, alpha=0.5)
    #p = Rectangle([-2/3.*b1,-0.5*b3], 4/3.*b1, b3, alpha=0.5)

    #a=0.532728
    #x0=latgen.pia[0]*sqrt(3)/2.*a
    #y0=latgen.pia[1]*a
    #z0=latgen.pia[2]*1.5
    #p = Rectangle([-x0,-y0], 2*x0, 2*y0, alpha=0.5)
    #print([-x0,-y0], 2*x0, 2*y0)
    kpth2=None
    kpth3=None
    exec(compile(open('2D_params.py', 'rb').read(), '2D_params.py', 'exec'))
    print('kpth=', kpth)
    corners = [[-1,-1],[1,-1],[1,1],[-1,1],[-1,-1]]
    kc = [array(eval(kpth))*latgen.pia for x,y in corners]
    kc = array(kc)
    ax.plot(kc[:,0],kc[:,1],kc[:,2], color='b', lw=1)
    XYZ = [array([ [kc[0,d],kc[1,d]],[kc[3,d],kc[2,d]] ]) for d in range(3)]
    ax.plot_surface(XYZ[0],XYZ[1],XYZ[2],alpha=0.5)

    if kpth2 is not None:
        print('kpth2=', kpth2)
        corners = [[-1,-1],[1,-1],[1,1],[-1,1],[-1,-1]]
        kc = [array(eval(kpth2))*latgen.pia for x,y in corners]
        kc = array(kc)
        ax.plot(kc[:,0],kc[:,1],kc[:,2], color='r', lw=1, alpha=0.3)
        XYZ = [array([ [kc[0,d],kc[1,d]],[kc[3,d],kc[2,d]] ]) for d in range(3)]
        ax.plot_surface(XYZ[0],XYZ[1],XYZ[2],alpha=0.3)
    if kpth3 is not None:
        print('kpth3=', kpth2)
        corners = [[-1,-1],[1,-1],[1,1],[-1,1],[-1,-1]]
        kc = [array(eval(kpth3))*latgen.pia for x,y in corners]
        kc = array(kc)
        ax.plot(kc[:,0],kc[:,1],kc[:,2], color='g', lw=1, alpha=0.3)
        XYZ = [array([ [kc[0,d],kc[1,d]],[kc[3,d],kc[2,d]] ]) for d in range(3)]
        ax.plot_surface(XYZ[0],XYZ[1],XYZ[2],alpha=0.3)
        

    #p = Rectangle([-0.5*b1,-0.5*b3], b1, b3, alpha=0.5)
    #ax.add_patch(p)
    #pathpatch_2d_to_3d(p, z = 0.0, normal=(-sqrt(3)/2.,0.5,0))

            
    ax.arrow3D(0,0,0,
               0.5,0,0,
               mutation_scale=20,
               arrowstyle="->",
               fc ='red',
               lw=0.2,
               linestyle='dashed')
    ax.arrow3D(0,0,0,
               0,0.5,0,
               mutation_scale=20,
               arrowstyle="->",
               fc ='green',
               lw=0.2,
               linestyle='dashed')
    ax.arrow3D(0,0,0,
               0,0,0.5,
               mutation_scale=20,
               arrowstyle="->",
               fc ='blue',
               lw=0.2,
               linestyle='dashed')
    ax.arrow3D(0,0,0,
               *cell[0],
               fc ='red',
               mutation_scale=20,
               arrowstyle="-|>")
    ax.arrow3D(0,0,0,
               *cell[1],
               fc ='green',
               mutation_scale=20,
                arrowstyle="-|>")
    ax.arrow3D(0,0,0,
               *cell[2],
               fc ='blue',
               mutation_scale=20,
               arrowstyle="-|>")

    # make the panes transparent
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # make the grid lines transparent
    ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    
    plt.show()
    
