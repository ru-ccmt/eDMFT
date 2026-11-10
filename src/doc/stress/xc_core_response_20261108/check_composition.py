"""First-k unmixed charge derivative; no normalization or screening claim."""

from pathlib import Path
import sys, time, json
import numpy as np

root = Path(__file__).resolve().parent
sys.path[:0] = [
    str(root / "stage"),
    str(root),
    "/Users/haule/projects/eDMFT_new/src/python",
]
from lapw_potential_response import NativeInterstitialPotential, NativePotentialBands
from lapw_complement_response import complement_electronic_operator
from stress_response_adjoint import PairedOperator
from check_xc_periodic import build as build_xc
from check_hxc_core import build_core

sys.path.insert(0, "/private/tmp/edmft_complement_response_20261108")
from build_reduced_material import build

case = sys.argv[1]
nr = {"FeSe": 66, "SrVO3": 43}[case]
source = root / case / "charge" / "0"
t = time.perf_counter()
f = build(case, source, nr)
s = f["mat"].s
b = f["bands"]
lines = (source / (case + ".vns")).read_text().splitlines()
j = next(i for i, l in enumerate(lines) if "NUMBER OF PW" in l)
n = int(lines[j].split()[0])
reps = np.array(
    [[int(l[a : a + 5]) for a in [3, 8, 13]] for l in lines[j + 1 : j + 1 + n]]
)
irs = []
for record in f["records"]:
    npw = record.nbasis - record.nlo
    irs.append(
        NativeInterstitialPotential(
            record.coefficients[:npw, :, 0],
            record.labels[:, :npw].T,
            reps,
            s.timat.transpose(0, 2, 1),
            s.tau,
            complex_values=f["mat"].complex,
        )
    )
bands = NativePotentialBands(b, irs, retained=nr)
e, inputs, outputs = complement_electronic_operator(f["complement"], fixed_number=False)
shape = b.coordinates.band.shape


def forward(x):
    rr, qr = outputs.unpack(x)
    rho = np.zeros(shape, complex)
    rho[:, :nr, :nr] = rr
    rho[:, nr:, :nr] = qr
    rho[:, :nr, nr:] = np.array(qr).swapaxes(1, 2).conj()
    return f["spatial"].action(rho)


def reverse(x):
    rho, _, _ = f["spatial"].pullback(x)
    return outputs.pack(rho[:, :nr, :nr], rho[:, nr:, :nr])


r = PairedOperator(
    outputs.size, len(f["reference"].vector()), forward, reverse, "spatial valence"
)
valence = r @ e @ bands.operator + f["direct"].operator @ bands.spherical_selector
xc, h, density = build_xc(case)
core, _ = build_core(case, h)
select, coords = h.native_potential_output(
    tuple(b.slices),
    {a: part.layout.harmonics for a, (part, ks) in bands.ns.items()},
    reps,
    complex_interstitial=f["mat"].complex,
)
assert (
    coords.nonspherical_shapes == bands.input.nonspherical_shapes
    and coords.size == bands.input.size
)
assert tuple(f["reference"].rad) == h.keys
# C acts in the full potential coordinates; valence uses the explicit native
# LAPW1 selection. This avoids silently dropping core on atoms with no
# correlated orbital or assigning a second spin/k weight to either map.
charge_from_potential = valence @ select + core.operator
combined = charge_from_potential @ xc.hxc_operator
d = np.load("/private/tmp/edmft_hartree_periodic_20261108/" + case + "_direction.npz")[
    "density"
]
result = dict(
    scope="consistent native central state; fixed-mu first-k complete-complement valence + core + Hartree/XC; unnormalized; no screened stress",
    preparation_seconds=time.perf_counter() - t,
)
t = time.perf_counter()
y = combined.matvec(d)
result["forward_seconds"] = time.perf_counter() - t
t = time.perf_counter()
result["transpose"] = combined.check_transpose(samples=1, tolerance=1e-10)
result["transpose_seconds"] = time.perf_counter() - t
np.savez(
    root / (case + "_hxc_valence_core.npz"),
    density=y,
    potential=xc.hxc_operator.matvec(d),
)

en,_,_=complement_electronic_operator(f['complement'],fixed_number=True)
valence_n=r@en@bands.operator+f['direct'].operator@bands.spherical_selector
combined_n=(valence_n@select+core.operator)@xc.hxc_operator
yn=combined_n.matvec(d)
result['fixed_n_transpose']=combined_n.check_transpose(samples=1,tolerance=1e-10)
np.savez(root/(case+'_hxc_valence_core_fixed_n.npz'),density=yn)

from wien_response_density import Density
from check_core import read_states

result["steps"] = {}
result["fixed_n_steps"] = {}
for step in [0.01, 0.005, 0.001, 0.0005]:
    end = []
    for q in [step, -step]:
        folder = root / case / "charge" / str(q)
        val = Density(folder / (case + ".clmval"), folder / (case + ".struct")).vector()
        native = read_states(root / case / "hxc_core" / str(q) / "core_native.bin")
        i = 0
        for a, atom in core.atoms.items():
            for state in atom.states:
                nn = native[i]
                i += 1
                val[core.slices[a]] += (
                    state.occ * (nn["p"] ** 2 + nn["q"] ** 2)[: atom.n]
                )
        end.append(val)
    numeric = (end[0] - end[1]) / (2 * step)
    result["steps"][step] = dict(
        relative=float(np.linalg.norm(y - numeric) / np.linalg.norm(numeric)),
        maximum=float(np.max(abs(y - numeric))),
    )
for step in [.001,.0005]:
    end=[]
    for q in [step,-step]:
        folder=root/case/'fixed_n'/str(q)
        val=Density(folder/(case+'.clmval'),folder/(case+'.struct')).vector()
        native=read_states(root/case/'hxc_core'/str(q)/'core_native.bin');i=0
        for a,atom in core.atoms.items():
            for state in atom.states:
                nn=native[i];i+=1
                val[core.slices[a]]+=state.occ*(nn['p']**2+nn['q']**2)[:atom.n]
        end.append(val)
    numeric=(end[0]-end[1])/(2*step)
    result['fixed_n_steps'][step]=dict(relative=float(np.linalg.norm(yn-numeric)/np.linalg.norm(numeric)),maximum=float(np.max(abs(yn-numeric))))
from check_oneparticle_stress_fd import read_basis_matrix
v=select.matvec(xc.hxc_operator.matvec(d))
ah,ao,ap,ahq,aoq=bands.output.unpack(bands.action(v))
result['matrix_steps']={}
for step in [.001,.0005]:
    hh=[];oo=[];hq=[];oq=[]
    for k,record in enumerate(f['records']):
        matrices=[]
        for q in [step,-step]:
            folder=root/case/'charge'/str(q)
            hm=read_basis_matrix(folder/f'native.ham.k{k+1:08d}.bin',b'EDMFTHV1')
            om=read_basis_matrix(folder/f'native.k{k+1:08d}.bin',b'EDMFTOV1')
            np.testing.assert_array_equal(hm.labels,record.labels)
            matrices.append((hm.matrix,om.matrix))
        c=record.coefficients[:,:,0]
        dh=c.conj().T@((matrices[0][0]-matrices[1][0])/(2*step)@c[:,:nr])
        do=c.conj().T@((matrices[0][1]-matrices[1][1])/(2*step)@c[:,:nr])
        hh.append(dh[:nr]);hq.append(dh[nr:]);oo.append(do[:nr]);oq.append(do[nr:])
    dh=bands.output.matrix.pack(np.array(hh),hq);do=bands.output.matrix.pack(np.array(oo),oq)
    result['matrix_steps'][step]=dict(H_relative=float(np.linalg.norm(dh-bands.output.matrix.pack(ah,ahq))/np.linalg.norm(dh)),O_relative=float(np.linalg.norm(do-bands.output.matrix.pack(ao,aoq))/np.linalg.norm(do)))
    np.savez(root/(case+'_hxc_matrix_fd_'+str(step)+'.npz'),analytic=bands.action(v),numeric=bands.output.pack(np.array(hh),np.array(oo),ap,hq,oq))
result['central_density_relative']=float(np.linalg.norm(f['spatial'].value()-f['reference'].vector())/np.linalg.norm(f['reference'].vector()))

from material_augmentation import read_vns
native_inputs=[]
for q in [.0005,-.0005]:
    folder=root/case/'charge'/str(q)
    sphere=read_vns(folder/(case+'.vsp'));nonsphere=read_vns(folder/(case+'.vns'))
    sph=b.pack({a:fields[0,0][:s.jrj[a]] for a,fields in enumerate(sphere)})
    ns={a:np.array([nonsphere[a][lm][:s.jrj[a]] for lm in part.layout.harmonics]).T for a,(part,ks) in bands.ns.items()}
    lines=(folder/(case+'.vns')).read_text().splitlines();i=next(i for i,line in enumerate(lines) if 'NUMBER OF PW' in line)
    pw=np.array([complex(float(line[18:37]),float(line[37:56])) for line in lines[i+1:i+1+len(reps)]])
    if not f['mat'].complex:pw=pw.real
    native_inputs.append(bands.input.pack(sph,ns,pw))
fdv=(native_inputs[0]-native_inputs[1])/.001
result['native_potential_steps']={}
for name,sl in [('spherical',slice(0,bands.input.spherical_size)),('nonspherical',slice(bands.input.spherical_size,bands.input.pw_offset)),('interstitial',slice(bands.input.pw_offset,None))]:
    result['native_potential_steps'][name]=dict(relative=float(np.linalg.norm(v[sl]-fdv[sl])/np.linalg.norm(fdv[sl])),maximum=float(np.max(abs(v[sl]-fdv[sl]))))
np.savez(root/(case+'_native_potential_comparison.npz'),analytic=v,numeric=fdv)
fh,fo,fp,fhq,foq=bands.output.unpack(bands.action(fdv))
result['H_with_native_potential_relative']=float(np.linalg.norm(dh-bands.output.matrix.pack(fh,fhq))/np.linalg.norm(dh))
# Substitute the independently measured H/O derivative, leaving the analytic
# projector and direct basis response intact. This is diagnostic only.
fd_matrices=bands.output.pack(np.array(hh),np.array(oo),ap,hq,oq)
response_fdho=r.matvec(en.matvec(fd_matrices))+f['direct'].operator.matvec(bands.spherical_selector.matvec(v))
response_fdho+=core.action(xc.hxc_operator.matvec(d))
result['fixed_n_with_native_HO_relative']=float(np.linalg.norm(response_fdho-numeric)/np.linalg.norm(numeric))

(root / (case + "_composition.json")).write_text(json.dumps(result, indent=2) + "\n")
print(case, json.dumps(result, indent=2), flush=True)
