"""Test stationarity with an eigenvector-dependent finite-band projector.

This is a finite Hilbert-space oracle, not a material stress correction.
The basis, orbital shapes, overlap and self-energy are fixed; only the
one-body potential changes. The density is set exactly equal to the
lattice density at t=0. Thus any residual comes from the window projector,
not from incomplete SCF convergence or omitted radial derivatives.
"""
import contextlib,io,json
with contextlib.redirect_stdout(io.StringIO()):
    import joint_projector_test as m
import numpy as np
m.dynamic=True
m.do=[np.zeros((m.n,m.n),complex) for _ in range(m.nk)]
m.dq=[np.zeros((m.n,m.nc),complex) for _ in range(m.nk)]
out={}
for full in [False,True]:
    if full:
        m.lo,m.hi=0,m.n;m.nw=m.n;m.W=slice(0,m.n)
        m.maps=[m.I([v]) for v in [1,m.n,1,m.n,m.n,m.n,m.nc,1,m.nc,1]]
    a=m.analytic(True,True)
    bare=a['pieces']['oneparticle']
    fd=(m.endpoint(1e-5,True)-m.endpoint(-1e-5,True))/(2e-5)
    # Gamma(V)=F_sp(V)-Tr[n_ref(V-V_ref)]+constant.
    out['complete_space' if full else 'finite_window']={
        'dFsp':fd,'Tr_rho_dV':bare,'dGamma_dV_at_n_equal_rho':fd-bare,
        'predicted_window_chain':a['pieces']['window']+a['pieces']['norm_window']}
    assert abs(fd-a['total'])<2e-9
    if full:assert abs(fd-bare)<2e-9
print(json.dumps(out,indent=2))
(m.root/'potential_stationarity_results.json').write_text(json.dumps(out,indent=2)+'\n')
