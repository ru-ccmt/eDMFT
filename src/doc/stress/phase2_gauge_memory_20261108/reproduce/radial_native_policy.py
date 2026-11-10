from pathlib import Path
import sys, ctypes as ct, json
import numpy as np
root=Path(__file__).resolve().parent
sys.path[:0]=[str(root/'stage'),'/Users/haule/projects/eDMFT_new/src/python']
from lapw_valence_state import _snapshot_inputs
from lapw_radial_native import NativeRadialContext
case=sys.argv[1];source=Path('/private/tmp/edmft_full_k_response_20261108')/case/'charge/0'
mat=_snapshot_inputs(case,source,case=='SrVO3');s=mat.s
lib=ct.CDLL(str(root/'native/libnative_quadrature.dylib'));fn=lib.native_outwin
arr=np.ctypeslib.ndpointer(dtype=np.float64,flags='F_CONTIGUOUS');fn.argtypes=[ct.c_int,ct.c_int,arr,arr,ct.c_double,ct.c_double,ct.c_int,ct.c_double,arr];fn.restype=None
charge=lib.native_charge
charge.argtypes=[ct.c_int,arr,ct.c_double,arr,ct.POINTER(ct.c_double)];charge.restype=None
rows=[]
for a in range(s.nat):
 for l in range(11):
    e=mat.header[a][l];lapw=e<150;e=e if lapw else e-200
    n=s.jrj[a];dx=np.log(s.rmt[a]/s.r0[a])/(n-1)
    elo=mat.lo[a][l,:mat.ilo[a][l]] if l<4 else np.array([])
    with NativeRadialContext(mat.v[a][0,0],r0=s.r0[a],dx=dx,rmt=s.rmt[a],energy=e,lo_energies=elo,ell=l,charge=s.Znuc[a],relativistic='RELA' in s.mode,lapw=lapw,quadrature='wien2k24',library=str(root/'build/liblapw_response.so')) as c:
        def integrate(a,b):
            factor=1/137.0359895**2 if 'RELA' in s.mode else 1e-22
            fields=np.asfortranarray(a[:,0]*b[:,0]+factor*a[:,1]*b[:,1])
            result=ct.c_double();charge(n,np.asfortranarray(c.radius),dx,fields,ct.byref(result))
            return result.value
        def native(en):
            out=np.empty((n,4),order='F');fn(n,int('RELA' in s.mode),np.asfortranarray(mat.v[a][0,0]/2),np.asfortranarray(c.radius),dx,en/2,l,s.Znuc[a],out)
            return out/np.sqrt(integrate(out,out))
        u=native(e);ud=(native(e+.004)-native(e-.004))/.008;ud-=integrate(u,ud)*u
        row=dict(a=a,l=l,energy=e,central=float(np.linalg.norm(u-c.primitive[:,:,0])/np.linalg.norm(u)),udot=float(np.linalg.norm(ud-c.primitive[:,:,1])/np.linalg.norm(ud)),native_norm=float(integrate(ud,ud)),context_norm=float(np.sum(c.overlap_metric*c.primitive[:,:,1]**2)))
        rows.append(row);print(row,flush=True)
(root/(case+'_radial_native_policy.json')).write_text(json.dumps(rows,indent=2))
