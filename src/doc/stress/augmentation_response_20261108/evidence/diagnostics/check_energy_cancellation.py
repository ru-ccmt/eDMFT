from material_augmentation import *
reports={}
for case in ['FeSe','SrVO3']:
 mat=Material(case);s=mat.s;root=R/case/'moving_probe';data=np.load(root/'complete_matrices.npz');e=data['central_energies'];h,o=data['H'],data['O']
 predicted=(np.diag(h)-e*np.diag(o)).real;rows=[]
 for step in [2e-3,1e-3,5e-4]:
  energies=[next(vectors(root/(f'full_{sign}_{step:.0e}')/(case+'.vector'),s.nat,mat.complex))[0] for sign in ['m','p']]
  fd=(energies[1]-energies[0])/(2*step);rows.append(dict(step=step,max_absolute=float(np.max(abs(fd-predicted))),relative=float(np.linalg.norm(fd-predicted)/np.linalg.norm(fd)),max_derivative=float(np.max(abs(fd)))))
 reports[case]=rows;print(case,rows,flush=True)
(R/'energy_cancellation.json').write_text(json.dumps(reports,indent=2)+'\n')
