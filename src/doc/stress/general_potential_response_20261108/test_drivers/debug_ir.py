from pathlib import Path
import sys,numpy as np
ROOT=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT/'stage'),'/Users/haule/projects/eDMFT_new/src/python']
from lapw_potential_response import NativeInterstitialPotential
from frozen_sigma_tangent import CachedFrozenSigmaResponse
from check_oneparticle_stress_fd import read_basis_matrix
from wstruct import Struct
case='SrVO3';source=Path('/private/tmp/edmft_coupled_valence_response_20261108/SrVO3_b703_one/central')
m,records=CachedFrozenSigmaResponse.from_export(source/'response');rec=records[0];npw=rec.nbasis-rec.nlo
s=Struct();s.ReadStruct(str(source/(case+'.struct')))
lines=(source/(case+'.vns')).read_text().splitlines();j=next(i for i,l in enumerate(lines) if 'NUMBER OF PW' in l);n=int(lines[j].split()[0]);reps=[];pw=[]
for line in lines[j+1:j+1+n]:
 reps.append([int(line[i:i+5]) for i in [3,8,13]]);pw.append(complex(float(line[18:37]),float(line[37:56])))
reps=np.array(reps);pw=np.array(pw)*.1
ends=[read_basis_matrix(ROOT/case/'interstitial'/(sign+'_0.001')/'native.ham.k00000001.bin',b'EDMFTHV1').matrix for sign in ['minus','plus']]
fd=(ends[1]-ends[0])/.002
print('npw',npw,'nbasis',rec.nbasis,'nonPW_fd',np.linalg.norm(fd[npw:]),'fdNorm',np.linalg.norm(fd))
for transpose in [False,True]:
 for signtau in [-1,1]:
  for signg in [-1,1]:
   op=NativeInterstitialPotential(np.eye(npw),signg*rec.labels[:,:npw].T,reps,s.timat.transpose(0,2,1) if transpose else s.timat,signtau*s.tau,complex_values=True)
   expect=op.coefficient_matrix(pw);d=fd[:npw,:npw]-expect
   print(transpose,signtau,signg,np.linalg.norm(d)/np.linalg.norm(fd))
   if transpose and signtau==signg==1:
    idx=np.unravel_index(np.argsort(abs(d).ravel())[-12:],d.shape)
    for i,j in zip(*idx):print(i,j,rec.labels[:,i]-rec.labels[:,j],expect[i,j],fd[i,j])
