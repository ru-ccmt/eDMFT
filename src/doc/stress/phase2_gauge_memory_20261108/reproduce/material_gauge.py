"""Native material gauge control with all reference energies shifted together.

This is NOT a claim that a fixed-linearization-energy policy is gauge invariant.
For a true gauge shift, V, E_l, E_LO, mu and the core reference shift together.
The core's shifted potential is then unchanged. Frozen Sigma(iw) is unchanged.
"""
import numpy as np
from lapw_complement_response import complement_electronic_operator


def gauge_check(f,state,select):
    nk,nr=f.model.eps.shape
    h=np.zeros((nk,nr,nr),complex);o=np.zeros_like(h)
    hq=[np.zeros((b.solver.nq,nr),complex) for b in f.bases]
    oq=[np.zeros_like(q) for q in hq]
    raw=np.zeros_like(f.model.raw)
    radial={};max_radial=0.
    for a,channels in f.contexts.items():
        radial[a]={}
        for l,c in channels.items():
            dp,dm=c.action(c.radius,1.,np.ones(c.nlo))
            radial[a][l]=(dp,dm)
            max_radial=max(max_radial,float(np.max(abs(dp))),float(np.max(abs(dm),initial=0)))
    for a,l,k,op in f.bands.spherical.spherical:
        c=f.contexts[a][l]
        dh,do=op.action(*radial[a][l],energy_dot=1.,lo_energies_dot=np.ones(c.nlo))
        h[k]+=dh[:nr];hq[k]+=dh[nr:];o[k]+=do[:nr];oq[k]+=do[nr:]
    for a,op,ks in f.bands.spherical.nonspherical:
        for k,dh in zip(ks,op.action({l:radial[a][l][0] for l in op.phi},{l:radial[a][l][1] for l in op.phi},retained=nr)):
            h[k]+=dh[:nr];hq[k]+=dh[nr:]
    for a,l,k,ic,cols,op in f.bands.spherical.projectors:
        raw[k,:,cols,ic]+=op.action(*radial[a][l])[:nr]
    field=np.zeros(state.hartree.size)
    pw=np.zeros(state.hartree.stars.nstar,complex)
    reps=state.hartree.stars.representatives
    pw[np.flatnonzero(np.all(reps==0,axis=1))[0]]=1.
    masked=state.hartree.reanalysis.action(pw)
    field[state.hartree.pw_offset:]=masked.view(float)
    native=select.matvec(field);_,_,p=f.bands.input.unpack(native)
    for k,op in enumerate(f.bands.ir):
        rr,qr=op.action(p,retained=nr);h[k]+=rr;hq[k]+=qr
    y,cross,info=f.complement.response(h,o,raw,hq,oq,fixed_number=True)
    spatial=f.spatial.action(f.spatial.coordinates.pack(y,cross))
    # Native spherical harmonics are truncated at LMAX. Thus H'=I may leave
    # small native basis-completeness errors; record them rather than fitting.
    return dict(radial_maximum=max_radial,
        H_identity_RR_norm=float(np.linalg.norm(h-np.eye(nr))),
        H_QR_norm=float(np.sqrt(sum(np.linalg.norm(q)**2 for q in hq))),
        O_derivative_norm=float(np.sqrt(np.linalg.norm(o)**2+sum(np.linalg.norm(q)**2 for q in oq))),
        raw_projector_derivative_norm=float(np.linalg.norm(raw)),
        dmu=float(info['dmu']),band_number_derivative=float(info['dN']),
        spatial_density_norm=float(np.linalg.norm(spatial)),
        density_coefficient_norm=float(np.sqrt(np.linalg.norm(y)**2+2*sum(np.linalg.norm(q)**2 for q in cross))),
        policy='constant +1 Ry potential; all valence linearization/LO energies and core shift +1 Ry; fixed N; not fixed-energy-policy probe')


def gauge_cutoff_check(f,state,select,source,l_extra=24):
    """Independently resolve the constant-potential residual into cutoffs.

    Test whether an omitted high-l plane-wave continuation or IR mask
    truncation explains the observed H'-O. Do not assume either is missing
    from the native finite basis. This hypothesis is REJECTED in the recorded
    FeSe/SrVO3 probes: adding the candidate increases the residual. These
    diagnostic integrals are not added to the source or stress calculation.
    """
    from scipy.special import spherical_jn,eval_legendre
    from numpy.polynomial.legendre import leggauss
    from pathlib import Path
    s=f.structure;volume=s.a*s.b*s.c;nr=f.model.eps.shape[1]
    pw=np.zeros(state.hartree.stars.nstar,complex)
    pw[np.flatnonzero(np.all(state.hartree.stars.representatives==0,axis=1))[0]]=1.
    field=np.zeros(state.hartree.size)
    field[state.hartree.pw_offset:]=state.hartree.reanalysis.action(pw).view(float)
    _,_,p=f.bands.input.unpack(select.matvec(field))
    gx,gw=leggauss(48);actual=[];corrected=[];tails=[];irs=[]
    for k,(record,base,ir) in enumerate(zip(f.records,f.bases,f.bands.ir)):
        c=np.column_stack((base.retained,base.complement));npw=record.nbasis-record.nlo
        g=record.labels[:,:npw].T
        with (Path(source)/f'native.k{k+1:08d}.bin').open('rb') as file:
            file.seek(44);kf=np.fromfile(file,'f8',3)
        q=2*np.pi*(g+kf)/np.array([s.a,s.b,s.c]);qn=np.linalg.norm(q,axis=1)
        unit=np.divide(q,qn[:,None],out=np.zeros_like(q),where=qn[:,None]!=0)
        angle=np.clip(unit@unit.T,-1,1);difference=np.linalg.norm(q[:,None]-q[None,:],axis=-1)
        theta=np.eye(npw,dtype=complex);tail=np.zeros_like(theta)
        site=0
        for a,mult in enumerate(s.mult):
            radius=s.rmt[a];x=difference*radius
            sphere=4*np.pi*radius**3/volume*np.divide(spherical_jn(1,x),x,out=np.full_like(x,1/3),where=x!=0)
            rr=(gx+1)*radius/2;ww=gw*radius/2*rr**2
            localtail=np.zeros((npw,npw))
            for l in range(max(f.contexts[a])+1,l_extra+1):
                jl=spherical_jn(l,qn[:,None]*rr)
                localtail+=(4*np.pi/volume)*(2*l+1)*eval_legendre(l,angle)*((jl*ww)@jl.T)
            for pos in s.pos[site:site+mult]:
                phase=np.exp(2j*np.pi*(g@pos));phase=phase.conj()[:,None]*phase[None,:]
                theta-=sphere*phase;tail+=localtail*phase
            site+=mult
        exact_ir=c[:npw].conj().T@(theta@c[:npw,:nr])
        high=c[:npw].conj().T@(tail@c[:npw,:nr])
        a,b=ir.action(p,retained=nr);native_ir=np.concatenate((a,b))
        mt=np.zeros_like(exact_ir)
        for a,l,kk,op in f.bands.spherical.spherical:
            if kk!=k:continue
            ctx=f.contexts[a][l]
            dh,_=op.action(np.zeros_like(ctx.primitive),np.zeros_like(ctx.lo_coefficients),energy_dot=1.,lo_energies_dot=np.ones(ctx.nlo))
            mt+=dh
        identity=np.zeros_like(mt);identity[:nr]=np.eye(nr)
        actual.append(mt+native_ir-identity)
        corrected.append(mt+exact_ir+high-identity)
        tails.append(high);irs.append(native_ir-exact_ir)
    norm=lambda xs:float(np.sqrt(sum(np.linalg.norm(x)**2 for x in xs)))
    return dict(high_l_max=l_extra,original_H_minus_O_columns_norm=norm(actual),
                missing_high_l_columns_norm=norm(tails),native_IR_minus_exact_mask_norm=norm(irs),
                residual_after_independent_cutoff_completion=norm(corrected))
