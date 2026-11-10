from pathlib import Path
import json,re,hashlib
import numpy as np
root=Path(__file__).resolve().parent;fac=14710.5164057
# Reuse the read-only vector record parser from the prior DFT audit.
space={'__file__':'/private/tmp/lapw2_stress_audit_20261107/compare_cases.py'}
exec(Path('/private/tmp/lapw2_stress_audit_20261107/compare_cases.py').read_text().split('report={}')[0],space)
vectors=space['vectors']
def tag(path,key):return re.findall(r'^'+re.escape(key)+r'\s*(.*)$',path.read_text(),re.M)[-1]
def number(path,key):return float(tag(path,key).split('=')[-1].split()[0])
def ten(path,key):return np.array([[float(z) for z in tag(path,f':{key}{i:03d}:').split()[:3]] for i in [1,2,3]])
old=json.loads(Path('/Users/haule/projects/eDMFT_new/src/doc/stress/g1_g2_second_pass_20261107/material_comparison.json').read_text())
out={};details={}
for label,d in old.items():out[label]={'analytic_mean_gpa':d['mean_gpa'],'numerical_mean_gpa':d['existing_fd_mean_gpa'],'difference_gpa':d['analytic_minus_fd_gpa'],'status':'previous second-pass row; unchanged'}
for mat,case,names in [('MnO','MnO',['MnO_DFT_m_0.001','MnO_DFT_center','MnO_DFT_p_0.001']),('MnTe','0_800_MnTe',['0_800_MnTe_MnTe_m_0.001','0_800_MnTe_0_800_MnTe','0_800_MnTe_MnTe_p_0.001'])]:
 dirs=[root/'DFT_runs'/n for n in names]
 if not all((f/'mixer_current.log').exists() for f in dirs):continue
 en=[number(f/(case+'.scfm'),':ENE  :') for f in dirs];en_temp=en.copy()
 if mat=='MnTe' and (root/'temps_check.json').exists():en=[x['energy_TEMPS_Ry'] for x in json.loads((root/'temps_check.json').read_text())]
 vol=number(dirs[1]/(case+'.scf0'),':VOL  :');fd=(en[2]-en[0])/.002/(3*vol)*fac
 center=dirs[1];stress=(ten(center/(case+'.scf0'),'STR_DENS')-ten(center/(case+'.scfc'),'STR_CORE')-ten(center/(case+'.scf2'),'STR_VAL'))*fac
 mix=ten(center/(case+'.scfm'),'STRESS_GPa');assert abs(mix-stress).max()<2e-5
 ve=[vectors(f,case) for f in dirs]
 consistency={key:all([x[key] for x in v]==[x[key] for x in ve[1]] for v in [ve[0],ve[2]]) for key in ['k','weight','basis','bands','pw_hash']}
 distances=[float(tag(f/(case+'.scfm'),':DIS  :').split()[-1]) for f in dirs]
 details[mat]={'directories':[str(p) for p in dirs],'energies_Ry':en,'TEMP_energies_Ry':en_temp,'energy_convention':'full Mermin TEMPS' if mat=='MnTe' and (root/'temps_check.json').exists() else 'ENE','volume_bohr3':vol,'basis_consistency':consistency,'nk':len(ve[1]),'bands':sorted({x['bands'] for x in ve[1]}),'charge_distances':distances,'diagonal_gpa':stress.diagonal().tolist(),'in2_integration':(center/(case+'.in2')).read_text().splitlines()[2]}
 out['DFT_'+mat]={'analytic_mean_gpa':float(np.trace(stress)/3),'numerical_mean_gpa':fd,'difference_gpa':float(np.trace(stress)/3-fd),'status':'one current full DFT evaluation of all three saved densities; see convergence/basis checks'}
for mat,case,folders,suffixes in [('MnO','MnO',['DMFT_MnO'],['scf2']),('MnTe','0_800_MnTe',['DMFT_MnTe_up','DMFT_MnTe_dn'],['scf2','scf2dn'])]:
 dirs=[root/p for p in folders]
 if not all('DMFT2 END' in (f/'current_run.log').read_text() for f in dirs):continue
 val=np.mean([ten(f/(case+'.'+s),'STR_VAL') for f,s in zip(dirs,suffixes)],axis=0)
 den=ten(dirs[0]/(case+'.scf0'),'STR_DENS');core=ten(dirs[0]/(case+'.scfc'),'STR_CORE');stress=(den-core-val)*fac
 numerical=json.loads((root/(mat+'_fd.json')).read_text());fd=numerical['stress_derivative_gpa']/3
 details['DMFT_'+mat]={'directories':[str(f) for f in dirs],'spin_mapping_count':len(dirs),'diagonal_gpa':stress.diagonal().tolist(),'individual_mapping_means_gpa':[float(np.trace((den-core-ten(f/(case+'.'+s),'STR_VAL'))*fac)/3) for f,s in zip(dirs,suffixes)],'frozen_input_checks':numerical['frozen_input_checks'],'endpoint_charges':[[tag(Path(numerical[key]['directory'])/(case+'.'+s),':CHA  :') for s in suffixes] for key in ['minus','plus']]}
 out['DMFT_'+mat]={'analytic_mean_gpa':float(np.trace(stress)/3),'numerical_mean_gpa':fd,'difference_gpa':float(np.trace(stress)/3-fd),'status':'current binary; matched prior controlled inputs; fixed-mu FD caveat' if mat=='MnO' else 'current binary; both magnetic mappings averaged once; endpoint convergence/window caveat'}
(root/'expanded_comparison.json').write_text(json.dumps({'rows':out,'details':details},indent=2)+'\n')
for k,x in out.items():print(k,f"{x['analytic_mean_gpa']:.9f} {x['numerical_mean_gpa']:.9f} {x['difference_gpa']:+.9f}")
for k,x in details.items():
 if 'basis_consistency' in x:print(k,'basis',x['basis_consistency'],'DIS',x['charge_distances'],'nk',x['nk'],'bands',x['bands'])
