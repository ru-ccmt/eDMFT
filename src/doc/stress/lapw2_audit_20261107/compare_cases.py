from pathlib import Path
import re,struct,hashlib,json
import numpy as np

root=Path(__file__).resolve().parent
FAC=14710.5164057
def tag(path,key):
    return re.findall(r'^'+re.escape(key)+r'\s*(.*)$',path.read_text(),re.M)
def number(path,key):return float(tag(path,key)[-1].split('=')[-1].split()[0])
def tensor(path,key):
    return np.array([[float(x) for x in tag(path,f':{key}{i:03d}:')[-1].split()[:3]] for i in range(1,4)])
def energy(path):return number(path,':ENE  :')
def vectors(folder,case):
    nat=int((folder/(case+'.struct')).read_text().splitlines()[1][27:30])
    nlos=set(map(int,re.findall(r'number of local orbitals, nlo \(hamilt\)\s+(\d+)',(folder/(case+'.output1')).read_text())))
    assert len(nlos)==1,nlos
    nlo=nlos.pop()
    result=[]
    with (folder/(case+'.vector')).open('rb') as f:
        def rec():
            raw=f.read(4)
            if not raw:return None
            n=struct.unpack('<i',raw)[0];data=f.read(n)
            assert len(data)==n and f.read(4)==raw
            return data
        for _ in range(2*nat):rec()
        while (h:=rec()) is not None:
            x,y,z,name,n,nb,w=struct.unpack('<3d10s2id',h)
            g=rec();assert len(g)==12*n
            pw=np.frombuffer(g,dtype='<i4').reshape(n,3)[:n-nlo]
            rows=sorted(map(tuple,pw.tolist()))
            digest=hashlib.sha256(np.array(rows,dtype='<i4').tobytes()).hexdigest()
            for _ in range(nb):rec();rec()
            result.append(dict(k=[x,y,z],weight=w,basis=n,bands=nb,pw_hash=digest))
    return result

report={}
for case,dirs in [('FeSe',['DFT_RELA_m001','DFT_RELA','DFT_RELA_p001']),
                  ('SrVO3',['DFT_m001','DFT_center','DFT_p001'])]:
    fs=[root/'runs'/(case+'_'+d) for d in dirs]
    if not all((f/'audit/corrected_total'/(case+'.scfm')).exists() for f in fs):continue
    vec=[vectors(f,case) for f in fs]
    keys=['k','weight','basis','bands','pw_hash']
    consistency={key:all([v[key] for v in vec[j]]==[v[key] for v in vec[1]] for j in [0,2]) for key in keys}
    bands=sorted(set(x['bands'] for x in vec[1]));basis=[x['basis'] for x in vec[1]]
    vols=[number(f/(case+'.scf0'),':VOL  :') for f in fs]
    es=[energy(f/'audit/corrected_total'/(case+'.scfm')) for f in fs]
    es_old=[energy(f/'audit/original'/(case+'.scfm')) for f in fs]
    # Isotropic a->a(1+t): dE/dt=Omega*Tr(sigma), hence divide by 3
    # for the mean normal stress. Endpoint lattice constants are rounded;
    # also report the slope wrt volume as a closely related diagnostic.
    num=(es[2]-es[0])/.002/(3*vols[1])*FAC
    num_vol=(es[2]-es[0])/(vols[2]-vols[0])*FAC
    center=fs[1]
    den=tensor(center/(case+'.scf0'),'STR_DENS')
    core=tensor(center/(case+'.scfc'),'STR_CORE')
    modes={}
    for mode in ['before','corrected_trace','corrected_full']:
        p=center/'audit'/mode/(case+'.scf2')
        if p.exists():
            # Match mixer's physical convention: sigma = DENS-CORE-VAL.
            stress=(den-core-tensor(p,'STR_VAL'))*FAC
            modes[mode]=dict(diagonal=np.diag(stress).tolist(),mean=float(np.trace(stress)/3),
                             difference_from_fd=float(np.trace(stress)/3-num))
    mixer=tensor(center/'audit/corrected_total'/(case+'.scfm'),'STRESS_GPa')
    assert np.max(abs(np.diag(mixer)-modes['corrected_trace']['diagonal']))<1e-5
    clm_before=center/'audit/before'/(case+'.clmval')
    clm_after=center/'audit/corrected_trace'/(case+'.clmval')
    full=center/'audit/corrected_full'/(case+'.clmval')
    def ordinary_records(path):
        return [line for line in path.read_text().splitlines()
                if line.startswith(':') and not line.startswith((':STR',':DBG_STR'))]
    records_equal=ordinary_records(center/'audit/before'/(case+'.scf2'))==ordinary_records(center/'audit/corrected_trace'/(case+'.scf2'))
    report[case]=dict(directories=[str(f) for f in fs],basis_consistency=consistency,
        nk=len(vec[1]),bands=bands,basis_min=min(basis),basis_max=max(basis),
        volume=vols,energy_before=es_old,energy_after=es,
        energy_change_ry=[a-b for a,b in zip(es,es_old)],
        DIS=[tag(f/'audit/corrected_total'/(case+'.scfm'),':DIS  :')[-1] for f in fs],
        numerical_mean_gpa=num,numerical_volume_slope_gpa=num_vol,
        analytic=modes,clmval_identical=clm_before.read_bytes()==clm_after.read_bytes(),
        ordinary_scf2_records_identical=records_equal,
        full_clmval_identical=full.exists() and full.read_bytes()==clm_after.read_bytes())
    (root/(case+'_vector_manifest.json')).write_text(json.dumps(vec,indent=2)+'\n')
(root/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
