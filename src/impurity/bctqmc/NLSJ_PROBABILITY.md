# Atomic NLSJ labels in Probability2.dat

`Probability.dat` remains unchanged. At each probability update, the master
rank also writes `Probability2.dat` with one row per atomic-Hamiltonian
eigenstate:

```text
state block eigen P E_cix E_F0 E_gc N L S J label_weight status
```

`state` is a one-based running number, while `block` is one-based and
`eigen` is zero-based, matching `Probability.dat`. `P` is the same
normalized `Prob_eig` entry as in that file. The `N` column is the
expectation value `⟨N⟩` (an integer in fixed-`N` atomic blocks).
`E_gc` is the eigenvalue of the solver Hamiltonian after the `F0` and
chemical-potential shifts in `boldc.cc`, but before its pseudo-particle
`lambda0` shift. `E_F0=E_gc+mu*⟨N⟩` and, for a single-site atom,
`E_cix=E_F0-F0*⟨N(N-1)/2⟩`; for a cluster `E_cix=E_F0`.
In a fixed-`N` atomic block, `E_cix` is precisely the eigenvalue of the
Hamiltonian stored in the CIX. In a mixed-`N` block, it is that Hamiltonian's
expectation value in the solver eigenstate.

For an f shell, `l=3` describes each electron, **not** the many-electron
`L` printed here. `L,S,J` are the dominant many-electron subspace of
the CIX Hamiltonian eigenstate. `label_weight` is the sum of squared
overlaps with all basis vectors carrying that same `(N,L,S,J)` label.
`status=mixed` when the dominant weight is below 0.7. A genuinely
degenerate Hamiltonian subspace may have basis-dependent per-eigenstate
labels; consider its members together.

The labels are optional. The generator `atomh.py` writes a `<cix>.nlsj`
sidecar whenever it generates a new CIX. Its CIX-style V3 layout starts with
`nstates nblocks norbitals cix_fingerprint` and the dimensions of every block.
The first section then has exactly one row per **Hamiltonian energy
eigenstate**:

```text
#  idx blk  N   L   S   J   Jz              E_CIX label_weight
```

`L,S,J` are the dominant label of that eigenstate, not the label of an
intermediate NLSJ basis vector. `label_weight` exposes mixing; exactly
degenerate eigenvectors can rotate, so an individual state's label can be
ambiguous. `-1 -1 -1` and weight zero mean that a block cannot be labeled.
The `state` column is a successive eigenstate number; it is not an occupation
bitstring or a claim that this vector equals the same-numbered direct state.
`E_CIX` comes from the Hamiltonian actually stored in the CIX, before the
solver's separate `F0` and chemical-potential corrections.

The second section has exactly one full eigenvector row per state, in the
same order, with components in CIX block order. Numeric coefficients have six
significant digits. There are no occupation-bitstring annotations; the full
vectors remain available for solver matching. The solver also accepts the
older V1 and V2 sidecars. For an existing
UTe2 CIX, without changing any CIX bytes or atomic energies, use:

```sh
python /path/to/eDMFT_new/src/impurity/atom_ctqmc/atomh.py -l 3 --label-existing /path/to/bctqmc.cix
```

Copy `bctqmc.cix` **and** `bctqmc.cix.nlsj` to the impurity working
directory on the cluster. The solver verifies a whole-file CIX fingerprint,
block dimensions, energies, and eigenvector orthonormality before attaching
labels. A missing
sidecar is harmless: `Probability2.dat` still contains `P`, the three
energy columns, and `N`, with `L,S,J` marked unavailable. An existing
but mismatched sidecar is an error rather than silently giving wrong labels.
Artificially split approximate-valence blocks are likewise left unlabeled.
