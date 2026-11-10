from pathlib import Path
import json,re,decimal,sys,importlib.util
sys.dont_write_bytecode=True
spec=importlib.util.spec_from_file_location('fk','/Users/haule/dbin/frozen_sigma_klein.py');fk=importlib.util.module_from_spec(spec);sys.modules['fk']=fk;spec.loader.exec_module(fk)
base=Path('/Users/haule/projects/eDMFT_new/src/doc/stress/edmft_refresh_20261107');data=json.loads((base/'comparison.json').read_text());root=Path('/private/tmp/numeric_stress_audit_20261107');out={};D=decimal.Decimal

def tag(s,key):
 a=re.findall(r'^'+re.escape(key)+r'\s*(.*)$',s,re.M)
 if not a:raise ValueError(key)
 return D(a[-1].split('=')[-1].split()[0].replace('D','E'))
def blocks(s):return re.split(r'(?=^:ITE\d+:)',s,flags=re.M)[1:]
def vals(s):
 d={k:tag(s,k) for k in [':ENE  :',':SUM  :',':VSUM :',':YSUM :',':EORB',':VOL  :',':DEN  :']};d['K']=d[':ENE  :']-d[':SUM  :']+d[':VSUM :']-d[':EORB'];return d
for mat,a in data.items():
 case=a['case'];tr=[]
 for name,p in zip(['minus','center','plus'],a['dirs']):
  path=Path(p);s=(path/(case+'.scf')).read_text();seq=[]
  for b in blocks(s):
   try:v=vals(b)
   except ValueError:continue
   seq.append(v)
  v=seq[-1];external=None
  try:external=fk.load_terms(path,case)
  except Exception:pass
  tr.append(v)
  item={'last_complete':{k:float(z) for k,z in v.items()},'recent_K_Ry':[float(z['K']) for z in seq[-6:]],'K_range_recent_Ry':float(max(z['K'] for z in seq[-3:])-min(z['K'] for z in seq[-3:])), 'external_K_Ry':external.frozen_reduced if external else None,'external_minus_samecycle_K_Ry':float(D(str(external.frozen_reduced))-v['K']) if external else None}
  out.setdefault(mat,{})[name]=item
 vol=tr[1][':VOL  :'];fac=D(str(fk.RY_PER_BOHR3_TO_GPA))/(6*D('.001')*vol)
 contributions={k:float((tr[2][k]-tr[0][k])*fac*sgn) for k,sgn in [(':ENE  :',1),(':SUM  :',-1),(':VSUM :',1),(':EORB',-1)]}
 contributions['density_DFT_part_ENE_minus_SUM_minus_EORB']=float(((tr[2][':ENE  :']-tr[2][':SUM  :']-tr[2][':EORB'])-(tr[0][':ENE  :']-tr[0][':SUM  :']-tr[0][':EORB']))*fac)
 slopes=[float((tr[1]['K']-tr[0]['K'])*fac*2),float((tr[2]['K']-tr[1]['K'])*fac*2)]
 shapes=[]
 for p in a['dirs']:
  ls=(Path(p)/(case+'.struct')).read_text().splitlines();shapes.append([float(ls[3][10*i:10*(i+1)]) for i in range(6)])
 strain=[[shapes[j][i]/shapes[1][i]-1 for i in range(3)] for j in [0,2]]
 hm,hp=-sum(strain[0])/3,sum(strain[1])/3
 # Three-point derivative for rounded, slightly unequal strain steps.
 fd_nonuniform=(-hp/(hm*(hm+hp))*float(tr[0]['K'])+(hp-hm)/(hm*hp)*float(tr[1]['K'])+hm/(hp*(hm+hp))*float(tr[2]['K']))/float(3*vol)*fk.RY_PER_BOHR3_TO_GPA
 out[mat].update(mean_stress_decimal_GPa=float((tr[2]['K']-tr[0]['K'])*fac),contributions_GPa=contributions,one_sided_GPa=slopes,strain_actual=strain,actual_step_mean_GPa=fd_nonuniform,printed_energy_rounding_bound_GPa=float(fac*D('1.2e-8')))
 print(mat,'FD',out[mat]['mean_stress_decimal_GPa'],'actualstrain',strain,'correctedFD',fd_nonuniform)
 print('contributions',contributions,'one-sided',slopes)
 print('crossfile', [out[mat][n]['external_minus_samecycle_K_Ry'] for n in ['minus','center','plus']],'recentKrange', [out[mat][n]['K_range_recent_Ry'] for n in ['minus','center','plus']])
old=json.loads(Path('/Users/haule/projects/eDMFT_new/src/doc/stress/g1_g2_second_pass_20261107/endpoint_references/SrVO3_existing_fd.json').read_text());a=out['SrVO3'];fac=D(str(fk.RY_PER_BOHR3_TO_GPA))/(6*D('.001')*D(str(old['reference_volume_bohr3'])))
a['RK7_contributions_GPa']={k:float((D(str(old['plus'][k]))-D(str(old['minus'][k])))*fac*sgn) for k,sgn in [('ene',1),('band_sum',-1),('vsum',1),('eorb',-1)]}
a['RK7_sigma_hashes_match_current']=old['frozen_input_checks'][0]['records'][0]['sha256']==data['SrVO3']['frozen_checks'][0]['records'][0]['sha256']
print('Sr old',a['RK7_contributions_GPa'],'sameSigma',a['RK7_sigma_hashes_match_current'])
(root/'bookkeeping.json').write_text(json.dumps(out,indent=2)+'\n')
