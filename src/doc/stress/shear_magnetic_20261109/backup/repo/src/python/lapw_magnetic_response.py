"""Collinear DMFT spins coupled to ONE unpolarized KS potential.

Each existing electronic branch includes the native nonmagnetic factor 2.
The physical density/scalar are therefore (branch_up + branch_down)/2.
Projector normalization belongs to each branch's actual orbital mapping;
charge normalization, HXC, core and screening belong to their sum only.

No spin-dependent DFT potential, spin-orbit coupling or impurity vertex is
implied. Common mu is checked at the snapshot boundary, never fitted per spin.
"""
from types import SimpleNamespace
import numpy as np
from lapw_density_response import checked
from response_parallel import parallel_of
from stress_response_adjoint import PairedOperator


class BlockCoordinates:
    def __init__(self, blocks):
        self.blocks=tuple(blocks)
        edges=np.r_[0,np.cumsum([b.size for b in blocks])]
        self.slices=tuple(slice(int(a),int(b)) for a,b in zip(edges[:-1],edges[1:]))
        self.size=int(edges[-1])

    def pack(self, *vectors):
        if len(vectors)!=len(self.blocks): raise ValueError("missing spin block")
        return np.concatenate([checked(v,(b.size,),"spin block",real=True) for v,b in zip(vectors,self.blocks)])

    def unpack(self, vector):
        vector=checked(vector,(self.size,),"spin coordinates",real=True)
        return tuple(vector[s] for s in self.slices)


def block_operator(operators):
    inputs=BlockCoordinates([SimpleNamespace(size=o.input_size) for o in operators])
    outputs=BlockCoordinates([SimpleNamespace(size=o.output_size) for o in operators])
    return PairedOperator(inputs.size,outputs.size,
        lambda x:outputs.pack(*(o.matvec(v) for o,v in zip(operators,inputs.unpack(x)))),
        lambda y:inputs.pack(*(o.rmatvec(v) for o,v in zip(operators,outputs.unpack(y)))),
        "independent spin electronic blocks"),inputs,outputs


def common_mu_operator(fixed_mu, mu_density, number_gradient, compressibility, *, parallel):
    """One rank-one fixed-N correction, including the moving overlap metric.

    E_N x = E_mu x - z * (b.x)/chi, where b.x is dN at fixed mu.
    z is d(rho_coeff)/dmu, while b includes Tr(rho O') as well as the
    coefficient-density trace. A local trace alone is not the total number.
    MPI reductions belong to b.x and z.y; chi is already globally reduced.
    """
    z=checked(mu_density,(fixed_mu.output_size,),"mu density",real=True)
    b=checked(number_gradient,(fixed_mu.input_size,),"number gradient",real=True)
    chi=float(compressibility)
    if not np.isfinite(chi) or abs(chi)<1e-12:
        raise ValueError("common fixed-N constraint is ill conditioned; use explicit fixed-mu policy")
    return PairedOperator(fixed_mu.input_size,fixed_mu.output_size,
        lambda x:fixed_mu.matvec(x)-z*(parallel.sum(b@x)/chi),
        lambda y:fixed_mu.rmatvec(y)-b*(parallel.sum(z@y)/chi),
        "two spins with one chemical-potential constraint")


class CollinearComplement:
    def __init__(self, up, down):
        self.branches=(up,down)
        self.parallel=parallel_of(up.model)
        for name in ('mu','temperature'):
            if not np.isclose(getattr(up.model.meta,name),getattr(down.model.meta,name),rtol=0,atol=1e-13):
                raise ValueError("spin branches do not share "+name)
        if parallel_of(down.model).size!=self.parallel.size or parallel_of(down.model).rank!=self.parallel.rank:
            raise ValueError("spin MPI partitions differ")

    def electronic_operator(self, *, fixed_number):
        from lapw_complement_response import complement_electronic_operator
        pieces=[complement_electronic_operator(b,fixed_number=False) for b in self.branches]
        independent,_,_=block_operator([x[0] for x in pieces])
        inputs=BlockCoordinates([x[1] for x in pieces]);outputs=BlockCoordinates([x[2] for x in pieces])
        if not fixed_number: return independent,inputs,outputs
        z=[]; gradients=[]; chi=0.
        for branch,(op,ins,outs) in zip(self.branches,pieces):
            m=branch.model; m.prepare()
            zero_q=[np.zeros((s.nq,s.nr),complex) for s in branch.solvers]
            z.append(outs.pack(m.mu_response,zero_q))
            identity=np.broadcast_to(np.eye(m.eps.shape[1]),m.eps.shape+(m.eps.shape[1],))
            number=op.rmatvec(outs.pack(identity,zero_q))
            # rho_coeff is contravariant. dN includes the metric derivative.
            number+=ins.pack(np.zeros_like(branch.rho),branch.rho,np.zeros_like(m.raw),zero_q,zero_q)
            gradients.append(number);chi+=m.compressibility
        # The factors 1/2 in physical N and chi cancel in delta mu. Keep the
        # densities in native branch units until the spatial reconstruction.
        constrained=common_mu_operator(independent,outputs.pack(*z),inputs.pack(*gradients),chi,parallel=self.parallel)
        return constrained,inputs,outputs

    def functional_gradients(self):
        from lapw_strain_functional import lattice_functional_gradients
        pieces=[lattice_functional_gradients(b) for b in self.branches]
        coordinates=BlockCoordinates([x[1] for x in pieces])
        gradients={key:.5*coordinates.pack(*(x[0][key] for x in pieces)) for key in ('projector','direct','lattice')}
        return gradients,coordinates


def combine_valence(up,down):
    """Combine fully validated native spin branches; no density solve occurs.

    The shared KS geometry, vector/potential fingerprints and density layouts
    must agree. Orbital mappings, Sigma and normalized projectors may differ.
    Branch contexts remain owned by the caller (typically an ExitStack).
    """
    from lapw_screened_state import _same_structure,_same_density_layout
    _same_structure(up.structure,down.structure)
    _same_density_layout(up.reference,down.reference)
    for key,value in up.fingerprints.items():
        if key.endswith(('.struct','.vector','.vsp','.vns','.in1','.in1c','.inm')) or key.startswith('native.'):
            if down.fingerprints.get(key)!=value: raise ValueError("spin KS snapshots differ: "+key)
    for name in ('size','spherical_size','nonspherical_shapes','nstar','complex_interstitial'):
        if getattr(up.bands.input,name)!=getattr(down.bands.input,name):
            raise ValueError("spin potential layout differs")
    # Native DMFT uses this precise Ry/eV conversion. Both exports must refer
    # to the common EF.dat, rather than silently retaining independent mus.
    ef_up = getattr(up, "fingerprint_paths", {}).get("EF.dat")
    ef_down = getattr(down, "fingerprint_paths", {}).get("EF.dat")
    if ef_up is None or ef_down is None or up.fingerprints["EF.dat"] != down.fingerprints["EF.dat"]:
        raise ValueError("magnetic response requires the same recorded EF.dat")
    mu = float(ef_up.read_text().split()[0]) / 13.60569253
    if not np.isfinite(mu) or not np.isclose(mu, up.model.meta.mu, rtol=0, atol=1e-10):
        raise ValueError("exported common chemical potential differs from EF.dat")
    complement=CollinearComplement(up.complement,down.complement)
    _,coordinates,density_coordinates=complement.electronic_operator(fixed_number=False)
    basis=PairedOperator(up.bands.operator.input_size,coordinates.size,
        lambda x:coordinates.pack(up.bands.operator.matvec(x),down.bands.operator.matvec(x)),
        lambda y:sum(v.bands.operator.rmatvec(z) for v,z in zip((up,down),coordinates.unpack(y))),
        "shared potential into two spin band maps")
    spatial=PairedOperator(density_coordinates.size,len(up.reference.vector()),
        lambda x:.5*sum(v.spatial.operator.matvec(z) for v,z in zip((up,down),density_coordinates.unpack(x))),
        lambda y:density_coordinates.pack(*(.5*v.spatial.operator.rmatvec(y) for v in (up,down))),
        "spin-averaged native spatial density")
    direct=PairedOperator(up.direct.operator.input_size,up.direct.operator.output_size,
        lambda x:.5*(up.direct.operator.matvec(x)+down.direct.operator.matvec(x)),
        lambda y:.5*(up.direct.operator.rmatvec(y)+down.direct.operator.rmatvec(y)),
        "spin-averaged moving-basis density")
    reference=up.reference.copy().add(down.reference).scale(.5)
    result=SimpleNamespace(branches=(up,down),complement=complement,model=up.model,
        structure=up.structure,reference=reference,records=up.records,
        bands=SimpleNamespace(operator=basis,output=coordinates,input=up.bands.input,
          spherical=up.bands.spherical,spherical_selector=up.bands.spherical_selector,ns=up.bands.ns,ir=up.bands.ir),
        spatial=SimpleNamespace(operator=spatial,coordinates=density_coordinates,
          value=lambda:.5*(up.spatial.value()+down.spatial.value())),
        direct=SimpleNamespace(operator=direct),fingerprints=up.fingerprints,
        fingerprint_paths=up.fingerprint_paths,
        mixer_normalization_enabled=getattr(up, "mixer_normalization_enabled", None),
        preparation_seconds=up.preparation_seconds+down.preparation_seconds,
        matrix_quadrature=up.matrix_quadrature)
    return result
