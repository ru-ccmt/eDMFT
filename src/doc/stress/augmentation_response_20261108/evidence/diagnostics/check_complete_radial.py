from material_augmentation import *
from scipy.special import eval_legendre
from lapw_response_assembly import NativeBandProjection
case=sys.argv[1];mat=Material(case);s=mat.s;ev,c,labels=mat.cv[0];npw=len(c)-mat.nlo_total
with (mat.source/'native.k00000001.bin').open('rb') as f:f.seek(44);kf=np.fromfile(f,'f8',3)
q=2*np.pi*(labels.T+kf)/np.array([s.a,s.b,s.c]);qn=np.linalg.norm(q,axis=1);qu=np.divide(q,qn[:,None],out=np.zeros_like(q),where=qn[:,None]>0)
ltop=int((mat.source/(case+('.in1c' if mat.complex else '.in1'))).read_text().splitlines()[1].split()[1])
# Use the same actual LTOP as native HAMILTON, not the larger radial-header size.
h=np.zeros((len(ev),len(ev)),complex);o=h.copy();h_reverse=h.copy();latom=0
for ia in range(s.nat):
 n=s.jrj[ia];r=s.r0[ia]*np.exp(np.log(s.rmt[ia]/s.r0[ia])/(n-1)*np.arange(n));dv=.07*r*np.sin(1.3*r)
 for ell in range(ltop+1):
  ctx=mat.radial(ia,ell,dv)
  inp=mat.work/f'radial_a{ia}_l{ell}.input';lines=inp.read_text().splitlines();lines[0]+=' '+str(npw)
  b=np.array([spherical_jn(ell,qn[:npw]*s.rmt[ia]),qn[:npw]*spherical_jn(ell,qn[:npw]*s.rmt[ia],derivative=True)])
  source=inp.with_suffix('.spherical_input');target=inp.with_suffix('.spherical_bin')
  source.write_text('\n'.join(lines)+'\n'+' '.join(f'{x:.17e}' for x in b.ravel(order='F'))+'\n')
  subprocess.run([str(R/'probe_spherical'),str(source),str(target)],check=True)
  nr=ctx['lo'].shape[1];m=npw+nr;rad=np.fromfile(target,'f8').reshape(2,m,m).transpose(0,2,1)
  active=list(range(npw));proxy=list(range(npw));factor=list(np.ones(npw))
  for j in range(nr):
   group=npw+mat.lo_indices[ia,ell,j];active.extend(group);proxy.extend([npw+j]*len(group));factor.extend([s.rmt[ia]**2]*len(group))
  active=np.array(active);proxy=np.array(proxy);factor=np.array(factor)
  qs=qu[active];g=labels[:,active].T
  phase=sum(np.exp(2j*np.pi*((g@pos)[None,:]-(g@pos)[:,None])) for pos in s.pos[latom:latom+s.mult[ia]])
  angular=4*np.pi/(s.a*s.b*s.c)*(2*ell+1)*eval_legendre(ell,np.clip(qs@qs.T,-1,1))*phase*factor[:,None]*factor[None,:]
  projection=NativeBandProjection(angular,c[active]);hr=rad[0][np.ix_(proxy,proxy)];orr=rad[1][np.ix_(proxy,proxy)]
  h+=projection.project(hr);o+=projection.project(orr);h_reverse+=projection.project(hr.T)
 latom+=s.mult[ia]
 print(case,'atom',ia,'assembled',flush=True)
root=R/case/'moving_probe';ns=np.load(root/'central_response.npz')['h_nonspherical'];h+=ns;h_reverse+=ns
rows=[]
for step in [2e-3,1e-3,5e-4]:
 hs=[];os=[]
 for sign in ['m','p']:
  d=root/(f'full_{sign}_{step:.0e}');hs.append(read_basis_matrix(d/'native.ham.k00000001.bin',b'EDMFTHV1').matrix);os.append(read_basis_matrix(d/'native.k00000001.bin',b'EDMFTOV1').matrix)
 dh=c.conj().T@((hs[1]-hs[0])/(2*step))@c;do=c.conj().T@((os[1]-os[0])/(2*step))@c
 row=dict(step=step,H_scale=float(np.linalg.norm(dh)),O_scale=float(np.linalg.norm(do)),H_relative=float(np.linalg.norm(dh-h)/np.linalg.norm(dh)),H_max=float(np.max(abs(dh-h))),O_relative=float(np.linalg.norm(do-o)/np.linalg.norm(do)),O_max=float(np.max(abs(do-o))),transposed_radial_H_relative=float(np.linalg.norm(dh-h_reverse)/np.linalg.norm(dh)))
 rows.append(row);print(case,row,flush=True)
(root/'complete_matrices.json').write_text(json.dumps(rows,indent=2)+'\n')
np.savez(root/'complete_matrices.npz',H=h,O=o,central_energies=ev)
