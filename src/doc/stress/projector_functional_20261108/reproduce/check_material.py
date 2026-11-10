"""Native endpoint projector scalar and physical-source adjoint checks.

Uses archived, independently run native endpoints. No material run is changed.
The scalar oracle uses full native eigenvector/metric transport and a low-rank
log-determinant identity, independently of the new analytic gradient chain.
"""
from pathlib import Path
import sys,time,json,resource
import numpy as np
from scipy.linalg import eigh,eigvalsh,lu_factor
root=Path(__file__).resolve().parent
if root.name=='reproduce':root=root.parent
sys.path[:0]=[str(root/'stage'),'/Users/haule/projects/eDMFT_new/src/python']
from frozen_sigma_response import FrozenSigmaResponse,herm
from lapw_valence_state import load_native_valence_state
from lapw_native_snapshot import read_basis_matrix
from lapw_potential_core_state import load_potential_core_state
from lapw_screened_state import NativeScreenedState
from lapw_projector_functional import PreparedProjectorFunctional
from wien_response_density import Density
from stress_response_adjoint import PairedOperator,ScreenedResponse


def spectral_reference(valence,source):
    records=[]
    for k,b in enumerate(valence.bases):
        eq,uq=eigh(b.solver.h)
        c=np.column_stack((b.retained,b.complement@uq))
        o=read_basis_matrix(source/f'native.k{k+1:08d}.bin',b'EDMFTOV1').matrix
        left=c.conj().T@o
        assert np.linalg.norm(left@c-np.eye(len(c)))<2e-9
        records.append(dict(energies=np.r_[valence.model.eps[k],eq],basis=c,left=left))
    return records


def native_projectors(folder,spectral):
    m,records=FrozenSigmaResponse.from_export(folder/'response')
    out=[]
    for k,(rec,ref) in enumerate(zip(records,spectral)):
        oo=read_basis_matrix(folder/f'native.k{k+1:08d}.bin',b'EDMFTOV1')
        assert np.array_equal(oo.labels,rec.labels)
        metric=herm(ref['basis'].conj().T@oo.matrix@ref['basis'])
        e,v=eigh(metric);assert min(e)>0
        root_metric=(v*np.sqrt(e))@v.conj().T
        # C_native U_native is invariant under arbitrary endpoint band phases.
        # The symmetric metric root places it in the common orthonormal gauge
        # used by the analytic generalized derivative. No band matching occurs.
        out.append(root_metric@ref['left']@rec.coefficients[:,:,0]@m.uc[k])
    return out


def logabsdet(matrices):
    # Independent scalar LU reference. On this host, NumPy's batched slogdet
    # warned even for a finite 5x5 batch with condition <=2.76 and singular
    # values >=0.425. Its numbers agreed with this route to 4.45e-16. Check
    # the actual pivots rather than masking floating-point warnings.
    values=[]
    for matrix in matrices:
        lu,_=lu_factor(matrix,check_finite=True)
        pivots=abs(np.diag(lu))
        if not np.all(np.isfinite(pivots)) or np.any(pivots==0):
            raise ValueError('singular or nonfinite scalar-oracle determinant')
        values.append(np.log(pivots).sum())
    return np.array(values)


def lattice_scalar(model,spectral,projectors):
    """Independent full-space scalar at FIXED central physical H and mu.

    Only U varies. Evaluate dynamic/static determinant ratios in orbital
    space by det(I-Sigma U^dagger(z-H0)^-1 U), but the static free energy uses
    the full finite native basis. This is the same exported quadrature, not
    a new tail approximation or the full Klein functional.
    """
    result=0.;T=model.meta.temperature;mu=model.meta.mu
    for k,(ref,u) in enumerate(zip(spectral,projectors)):
        eps=ref['energies'];norb=u.shape[1]
        static=np.diag(eps-mu).astype(complex)+(u@model.sinf)@u.conj().T
        eigenvalues=eigvalsh(herm(static),check_finite=False)
        value=-T*np.logaddexp(0,-eigenvalues/T).sum()
        for j in range(0,len(model.meta.omega),64):
            sl=slice(j,j+64);z=mu+1j*model.meta.omega[sl]
            m=np.einsum('ni,wn,nj->wij',u.conj(),1/(z[:,None]-eps),u,optimize=True)
            ls=logabsdet(np.eye(norb)-model.sinf@m)
            ld=logabsdet(np.eye(norb)-model.sigma[sl]@m)
            value+=2*T*np.dot(model.meta.weights[sl],ls-ld)
        result+=2*model.kw[k]*value
    return float(result)


def run(case):
    source=Path('/private/tmp/edmft_full_k_response_20261108')/case/'charge/0'
    output=root/case;output.mkdir(exist_ok=True)
    previous=Path('/private/tmp/edmft_screened_native_20261108/three_k')/case
    with load_native_valence_state(source,case,complex_vectors=(case=='SrVO3'),energy_policy='fixed',
                                  matrix_quadrature='wien2k24',library='/Users/haule/dbin/liblapw_response.so') as v:
        pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
        core=Density(Path('/private/tmp/edmft_xc_core_response_20261108')/case/'hxc_core/0'/(case+'.clmcor'),source/(case+'.struct'))
        state=NativeScreenedState(v,pc,core,target_charge={'SrVO3':85,'FeSe':120}[case],fixed_number=True)
        start=time.monotonic();obs=PreparedProjectorFunctional(state)
        report=dict(case=case,kpoints=len(v.records),retained=v.model.eps.shape[1],
                    scope='projector-only regularized lattice scalar; not total stress',
                    gradient_preparation_seconds=time.monotonic()-start,central=state.metadata)
        direction=np.load('/private/tmp/edmft_hartree_periodic_20261108/'+case+'_direction.npz')['density']
        potential=pc.potential_operator.matvec(direction)
        joint=state.chain.basis.matvec(potential)
        analytic=obs.potential_contraction(potential)
        report['analytic_density_direction_derivative_Ry']=analytic
        report['joint_vs_native_pullback_error']=float(abs(analytic-obs.joint_gradient@joint))
        h,o,p,hq,oq=obs.coordinates.unpack(joint)
        gh,go,gp,ghq,goq=obs.coordinates.unpack(obs.joint_gradient)
        inner=lambda a,b:float(np.vdot(a,b).real)
        report['contributions_Ry']=dict(H_RR=inner(gh,h),O_RR=inner(go,o),raw=inner(gp,p),
                                      H_QR=2*sum(inner(a,b) for a,b in zip(ghq,hq)),
                                      O_QR=2*sum(inner(a,b) for a,b in zip(goq,oq)))
        assert abs(sum(report['contributions_Ry'].values())-analytic)<1e-10
        spectral=spectral_reference(v,source)
        center=native_projectors(source,spectral)
        s0=lattice_scalar(v.model,spectral,center)
        report['central_scalar_full_minus_retained_Ry']=s0-v.model.lattice_scalar()
        report['native_steps']={}
        for step in [.001,.0005]:
            values=[]
            for sign in [1,-1]:
                folder=source.parent/str(sign*step)
                values.append(lattice_scalar(v.model,spectral,native_projectors(folder,spectral)))
            fd=(values[0]-values[1])/(2*step)
            row=dict(fd_Ry=fd,absolute_error_Ry=abs(fd-analytic),relative_error=abs(fd-analytic)/max(abs(analytic),1e-30))
            report['native_steps'][step]=row
            print(case,'native scalar',step,row,flush=True)
        fd1=report['native_steps'][.001]['fd_Ry'];fd2=report['native_steps'][.0005]['fd_Ry']
        report['richardson']=dict(fd_Ry=(4*fd2-fd1)/3,relative_error=abs((4*fd2-fd1)/3-analytic)/max(abs(analytic),1e-30))
        # Existing hydrostatic q uses the identical prepared charge Jacobian.
        # This first readout is the INDUCED contribution only; no explicit
        # strain potential, geometric band derivative or baseline is added.
        b=np.load(previous/'source_0.000125.npy');q=np.load(previous/'analytic_density.npy')
        direct=float(obs.charge_gradient@q)
        start=time.monotonic();lam,adj=state.screened.solve(obs.charge_gradient,adjoint=True,rtol=2e-9,atol=1e-12)
        reverse=float(lam@b)
        old=json.loads((previous/'screening_report.json').read_text())
        bound=float(np.linalg.norm(lam)*old['forward']['residual_norm']+np.linalg.norm(q)*adj['residual_norm'])
        report['screened_induced']=dict(adjoint=adj,forward_Ry=direct,adjoint_Ry=reverse,
                                       difference_Ry=abs(direct-reverse),residual_bound_Ry=bound)
        assert abs(direct-reverse)<=bound*1.001+1e-12
        np.savez(output/'projector_covectors.npz',joint=obs.joint_gradient,potential=obs.potential_gradient,
                 charge=obs.charge_gradient,adjoint=lam,direction=direction)
        report['peak_process_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        report['source_fingerprints']=v.fingerprints
        (output/'REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
        print(case,'done',json.dumps({k:x for k,x in report.items() if k not in ['source_fingerprints','central']}),flush=True)

if __name__=='__main__':run(sys.argv[1])
