from pathlib import Path
import sys, json, resource, collections, types, time, os
import numpy as np
root=Path(__file__).resolve().parent
sys.path[:0]=[str(root/'stage'),'/Users/haule/projects/eDMFT_new/src/python']
from lapw_valence_state import load_native_valence_state
case=sys.argv[1]
source=Path('/private/tmp/edmft_full_k_response_20261108')/case/'charge/0'
if os.environ.get('EDMFT_FULL_MESH')=='1':source=source.parents[2]/'full_bz'/case/'charge/0'
with load_native_valence_state(source,case,complex_vectors=case=='SrVO3',energy_policy='fixed',matrix_quadrature='wien2k24',library='/private/tmp/edmft_phase2_gauge_memory_20261108/build/liblapw_response.so') as f:
    seen=set();roots={};groups=collections.Counter();counts=collections.Counter()
    def walk(x,path):
        if id(x) in seen:return
        seen.add(id(x))
        if isinstance(x,np.ndarray):
            base=x
            while isinstance(base.base,np.ndarray):base=base.base
            if id(base) not in roots:
                roots[id(base)]=dict(path=path,shape=base.shape,bytes=base.nbytes,dtype=str(base.dtype))
            return
        if isinstance(x,dict):
            for k,v in x.items():walk(v,path+'.'+str(k))
        elif isinstance(x,(list,tuple)):
            for i,v in enumerate(x):walk(v,path+f'[{i}]')
        elif not isinstance(x,(type,types.ModuleType,types.FunctionType,types.MethodType)) and hasattr(x,'__dict__'):
            for k,v in vars(x).items():walk(v,path+'.'+k)
    walk(f,'state')
    import re
    for a in roots.values():
        key=re.sub(r'\[\d+\]','[*]',a['path'])
        groups[key]+=a['bytes'];counts[key]+=1
    out=dict(case=case,preparation_seconds=f.preparation_seconds,rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,owned_bytes=sum(a['bytes'] for a in roots.values()),
             groups=[dict(path=k,bytes=v,count=counts[k]) for k,v in groups.most_common()],
             largest=sorted(roots.values(),key=lambda x:x['bytes'],reverse=True)[:30])
    (root/(case+'_memory'+('_new' if os.environ.get('EDMFT_PROFILE_NEW')=='1' else '')+'.json')).write_text(json.dumps(out,indent=2))
    print(case,'seconds',f.preparation_seconds,'peak GB',out['rss']/1e9,'numpy GB',out['owned_bytes']/1e9,flush=True)
    for k,v in groups.most_common(25):print(f'{v/1e6:.3f} MB x{counts[k]} {k}',flush=True)
