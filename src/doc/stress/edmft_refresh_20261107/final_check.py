from pathlib import Path
import sys,importlib.util,json,re,hashlib,datetime,shutil
sys.dont_write_bytecode=True
spec=importlib.util.spec_from_file_location('fsk','/Users/haule/dbin/frozen_sigma_klein.py');fsk=importlib.util.module_from_spec(spec);sys.modules['fsk']=fsk;spec.loader.exec_module(fsk)
ROOT=Path('/private/tmp/edmft_refresh_check_20261107'); old=json.loads((ROOT/'results.json').read_text());out={};fac=14710.5164057

def text(p):return p.read_text(errors='replace') if p.exists() else ''
def matches(s,key):return re.findall(r'^'+re.escape(key)+r'\s*(.*)$',s,re.M)
def tensor(s,key):return [[float(z) for z in matches(s,f':{key}{i:03d}:')[-1].split()[:3]] for i in [1,2,3]]
def number(s,key):return float(matches(s,key)[-1].split('=')[-1].split()[0])
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None

def energy(d,case):
 nat=int((d/(case+'.struct')).read_text().splitlines()[1][27:30])
 cand=[p for p in d.glob(case+'.energy*') if re.fullmatch(re.escape(case)+r'\.energy(?:_\d+)?',p.name) and p.stat().st_size]
 newest=max(p.stat().st_mtime for p in cand)
 files=sorted([p for p in cand if newest-p.stat().st_mtime<30],key=lambda p:p.name)
 if any(re.search(r'_\d+$',p.name) for p in files):files=[p for p in files if re.search(r'_\d+$',p.name)]
 rows=[]
 for p in files:
  ls=p.read_text().splitlines();i=2*nat
  while i<len(ls):
   h=ls[i];xyz=[float(h[j:j+19]) for j in [0,19,38]];a=h[57:].split();n=int(a[1]);nb=int(a[2]);w=float(a[3]);ev=[float(ls[i+1+j].split()[1]) for j in range(nb)]
   rows.append({'k':xyz,'n':n,'nb':nb,'w':w,'ev':ev});i+=nb+1
 rows.sort(key=lambda r:r['k']);return files,rows
for mat,a in old.items():
 dirs=[Path(r['directory']) for r in a['runs']];d=dirs[1];case=next(k.split('.scfm')[0] for k in a['runs'][1]['files'] if k.endswith('.scfm'))
 snapshot=ROOT/'snapshots'/mat;snapshot.mkdir(parents=True,exist_ok=True)
 for i,dd in enumerate(dirs):
  dest=snapshot/['minus','center','plus'][i];dest.mkdir(exist_ok=True)
  for suffix in ['scf','scfm','scf0','scfc','scf2','scf2dn','scf1','in0','in1','in1c','in2','in2c','indmfl','indmfldn','indmf2','indmf2dn','struct']:
   p=dd/(case+'.'+suffix)
   if p.exists():shutil.copy2(p,dest/p.name)
  for name in ['info.iterate','EF.dat','params.dat','Eorb_imp.dat','dmft2.error.0','dmft2dn.error.0','dmft2_info.out']:
   p=dd/name
   if p.exists():shutil.copy2(p,dest/name)
 s=text(snapshot/'center'/(case+'.scf'));sm=text(snapshot/'center'/(case+'.scfm'))
 ten=tensor(s,'STRESS_GPa');den=tensor(s,'STR_DENS');core=tensor(s,'STR_CORE');val=tensor(s,'STR_VAL')
 ass=[[fac*(den[i][j]-core[i][j]-val[i][j]) for j in range(3)] for i in range(3)]
 mean=sum(ten[i][i] for i in range(3))/3;vol=number(s,':VOL  :')
 terms=[fsk.load_terms(dd,case) for dd in [dirs[0],dirs[2]]]
 fd=(terms[1].frozen_reduced-terms[0].frozen_reduced)/(.006*vol)*fsk.RY_PER_BOHR3_TO_GPA
 res={'dirs':list(map(str,dirs)),'case':case,'analytic_mean_GPa':mean,'analytic_diagonal_GPa':[ten[i][i] for i in range(3)],'numerical_mean_GPa':fd,'difference_GPa':mean-fd,'assembly_max_deviation_GPa':max(abs(ten[i][j]-ass[i][j]) for i in range(3) for j in range(3)),'volume_bohr3':vol,'Kminus_Ry':terms[0].frozen_reduced,'Kplus_Ry':terms[1].frozen_reduced,'frozen_checks':fsk.frozen_file_checks(dirs),'central_iteration':matches(s,':ITE')[-1:] if ':ITE' in s else [],'energies':[]}
 allrows=[]
 for i,dd in enumerate(dirs):
  files,rows=energy(dd,case);allrows.append(rows)
  ss=text(snapshot/['minus','center','plus'][i]/(case+'.scf'));p1=dd/(case+'.in1');p1=p1 if len(text(p1).splitlines())>1 else dd/(case+'.in1c')
  r={'files':[p.name for p in files],'nk':len(rows),'bands':sorted(set(r['nb'] for r in rows)),'dimension_range':[min(r['n'] for r in rows),max(r['n'] for r in rows)],'rkmax':float(text(p1).splitlines()[1].split()[0]),'mode':text(dd/(case+'.struct')).splitlines()[2],'distance':matches(ss,':DIS  :')[-1:],'charge':matches(ss,':CHA  :')[-1:],'pressure_history':matches(ss,':PRESSURE:')[-6:],'mu_input':text(dd/(case+'.indmf2')).splitlines()[-1], 'input_hashes':{name:sha(dd/name) for name in ['sig.inp1','sig.inp2','sig.inp1dn','sig.inp2dn','projectorw.dat',case+'.indmfl',case+'.indmfldn',case+'.klist'] if (dd/name).exists()}}
  ind=text(dd/(case+'.indmfl')).splitlines()[0].split();lo,hi=int(ind[0]),int(ind[1]);r['projector_window']=[lo,hi]
  if hi>lo>0:
   r['window_min_gaps_Ry']=[min(rr['ev'][lo-1]-rr['ev'][lo-2] for rr in rows),min(rr['ev'][hi]-rr['ev'][hi-1] for rr in rows if rr['nb']>hi)]
  res['energies'].append(r)
 res['energy_metadata_equal']={key:all([r[key] for r in allrows[i]]==[r[key] for r in allrows[1]] for i in [0,2]) for key in ['k','w','n','nb']}
 res['input_hashes_equal']={key:len({r['input_hashes'].get(key) for r in res['energies']})==1 for key in res['energies'][1]['input_hashes']}
 res['iteration_tag']=re.findall(r'^:ITE.*$',s,re.M)[-1]
 res['central_scf_time']=datetime.datetime.fromtimestamp((d/(case+'.scf')).stat().st_mtime).isoformat(timespec='seconds')
 res['standalone_scf2_complete']=all(':VSUM' in text(d/name) for name in fsk.dmft2_files(d,case)) if False else {name:':VSUM' in text(d/name) for name in [case+'.scf2']+([case+'.scf2dn'] if mat=='MnTe' else [])}
 if mat=='MnTe' and all(res['standalone_scf2_complete'].values()):
  vals=[tensor(text(d/(case+'.'+suffix)),'STR_VAL') for suffix in ['scf2','scf2dn']]
  res['raw_vs_merged_valence_max_difference_Ry_per_bohr3']=max(abs((vals[0][i][j]+vals[1][i][j])/2-val[i][j]) for i in range(3) for j in range(3))
  res['individual_mapping_means_GPa']=[sum(fac*(den[i][i]-core[i][i]-v[i][i]) for i in range(3))/3 for v in vals]
 out[mat]=res
 print(mat,mean,fd,mean-fd,'diag',res['analytic_diagonal_GPa'],'assembly error',res['assembly_max_deviation_GPa'])
 print(' basis',res['energy_metadata_equal'],'inputs',res['input_hashes_equal'])
 print(' setup',[(r['rkmax'],r['nk'],r['bands'],r['projector_window'],r.get('window_min_gaps_Ry')) for r in res['energies']])
 print('DIS',[r['distance'] for r in res['energies']],'last cycle',res['iteration_tag'])
 print('raw complete',res['standalone_scf2_complete'],'magnetic means',res.get('individual_mapping_means_GPa'))
(ROOT/'comparison.json').write_text(json.dumps(out,indent=2)+'\n')
