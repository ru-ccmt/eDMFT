"""WIEN cubic harmonic coordinate maps (SUMFAC, SROLYL, CUB_XC_BACK).

Stored density channels remain separate real harmonics. Native XC first
contracts each cubic multiplet with c, evaluates the combined harmonic,
and stores the result in its FIRST channel divided by c_first. ENERGY
pairs both fields after contraction: the metric on each group is c c^T.
This metric is not the identity and the output adapter is not C^T.
"""
import numpy as np

# Normalized cubic-harmonic coefficients, using WIEN signed-l labels.
_SQUARED = {
 (0,0):(1,1.), (3,2):(1,1.),
 (4,0):(1,7/12), (4,4):(1,5/12),
 (6,0):(1,1/8), (6,2):(1,11/16), (6,4):(-1,7/8), (6,6):(-1,5/16),
 (7,2):(1,13/24), (7,6):(1,11/24),
 (8,0):(1,33/64), (8,4):(1,7/48), (8,8):(1,65/192),
 (9,2):(1,3/16), (9,4):(1,17/24), (9,6):(-1,13/16), (9,8):(-1,7/24),
 (10,0):(1,65/384), (10,2):(1,247/384), (10,4):(-1,11/32),
 (10,6):(1,19/768), (10,8):(-1,187/384), (10,10):(-1,85/256),
}

def coefficient(l,m):
    try:
        sign, square = _SQUARED[abs(l),m]
    except KeyError as exc:
        raise ValueError("unsupported cubic harmonic") from exc
    return sign*np.sqrt(square)


def cubic_maps(harmonics):
    """Return C (stored -> compressed) and D (compressed -> stored XC)."""
    lm = list(map(tuple, harmonics))
    groups=[]; i=0
    while i<len(lm):
        l,m=lm[i]
        count=1 if (l,m) in [(0,0),(-3,2)] else 2 if l in (4,6,-7,-9) else 3 if l in (8,10) else 0
        if not count or lm[i:i+count] != [(l,m+4*j) for j in range(count)]:
            raise ValueError("invalid or incomplete native cubic multiplet")
        groups.append((i,count)); i+=count
    c=np.zeros((len(groups),len(lm))); d=np.zeros(c.T.shape)
    for row,(i,count) in enumerate(groups):
        for j in range(count): c[row,i+j]=coefficient(*lm[i+j])
        d[i,row]=1/c[row,i]
    return c,d


def boundary_factor(l,m):
    """POTFAC's cubic boundary conversion; radial volume term is unchanged."""
    return 1. if (l,m) in [(0,0),(-3,2)] else .5/coefficient(l,m)**2


def cubic_potential_gaunt(radial_counts, harmonics):
    """CBCOMB followed by cubic_atom_gaunt in ordinary real-harmonic units.

    CBCOMB's 1/sqrt(2) for m!=0 is already in real_harmonic_gaunt. The
    remaining operation is C^T C on stored potential channels. Replacing it
    by identity gives a different H_NS and breaks the potential response.
    """
    from lapw_augmentation_response import real_harmonic_gaunt
    c, _ = cubic_maps(harmonics)
    angular = np.array([real_harmonic_gaunt(max(radial_counts), l, m)
                        for l, m in harmonics])
    return np.einsum("ab,aij->bij", c.T @ c, angular)
