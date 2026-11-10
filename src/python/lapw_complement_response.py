"""Complete an empty LAPW response complement without high-state eigenvectors.

The retained space R must contain the DMFT window and every appreciably occupied
KS state. Q is the entire remaining finite coefficient space, not a second band
window. All central representations are O-orthonormal. Energies/potentials are Ry.

For moving LAPW coefficients, the missing mixing T solves
    H_Q T - T eps_R = -(delta H_QR - delta O_QR eps_R).
The retained electronic problem sees delta phi_R + T^dagger phi_Q. Its density
derivative is supplemented by delta rho_QR = T rho_R. This includes the moving
window/normalization feedback through the retained electronic response exactly
once. It is not an additional independent stress term.

Thermal Q occupations are neglected only after explicit occupation AND slope
bounds pass. This module supplies a reusable analytic response component; it
does not change production dmft2 or assemble a screened stress tensor.
"""

from dataclasses import dataclass
import time
import numpy as np
from scipy.linalg import cholesky, eigh_tridiagonal, hessenberg, qr, solve_triangular
from scipy.special import expit

from lapw_density_response import checked, herm
from stress_response_adjoint import HermitianCoordinates, PairedOperator


class ShiftedComplement:
    """Prepare all shifted Hermitian solves using one tridiagonal reduction.

    H_Q = U tridiag(b*,d,b) U^dagger. Each retained epsilon gives a positive
    definite tridiagonal system, so its LDL^dagger factors cost O(Nq), not a
    separate O(Nq^3) factorization. Forward and transpose reuse the SAME factors.
    The only complement eigenvalue computed is its lower spectral endpoint,
    used for validation. No complement eigenvectors or Green functions are used.

    Preparation remains dense and cubic in Nq; actions are O(Nq^2 Nr), with
    O(Nq^2 + Nq Nr) stored arrays. This is a direct baseline, not a claim of
    linear scaling or an optimal large-basis iterative implementation.
    """

    def __init__(
        self,
        hamiltonian,
        energies,
        *,
        mu,
        temperature,
        gap_tolerance=1e-8,
        charge_tolerance=1e-12,
        slope_tolerance=1e-12,
        active_columns=None,
    ):
        started = time.perf_counter()
        h = np.asarray(hamiltonian)
        e = np.asarray(energies)
        if h.ndim != 2 or h.shape[0] != h.shape[1] or not len(h):
            raise ValueError("empty or nonsquare complement Hamiltonian")
        if e.ndim != 1 or not len(e):
            raise ValueError("empty retained spectrum")
        self.h = np.array(
            checked(h, h.shape, "complement Hamiltonian", hermitian=True), complex
        )
        self.energies = np.array(
            checked(e, e.shape, "retained energies", real=True), float
        )
        values = np.asarray(
            [mu, temperature, gap_tolerance, charge_tolerance, slope_tolerance]
        )
        if (
            np.iscomplexobj(values)
            or not np.all(np.isfinite(values))
            or np.any(values[1:] <= 0)
        ):
            raise ValueError("invalid complement tolerances or temperature")
        self.mu, self.temperature = float(mu), float(temperature)
        self.nq, self.nr = len(h), len(e)
        reduced, self.rotation = hessenberg(self.h, calc_q=True, check_finite=False)
        d = np.diag(reduced).real.copy()
        b = np.diag(reduced, -1).copy()
        tri = np.diag(d).astype(complex) + np.diag(b, -1) + np.diag(b.conj(), 1)
        self.reduction_error = float(
            np.linalg.norm(reduced - tri) / max(1.0, np.linalg.norm(h))
        )
        if self.reduction_error > 1e-11:
            raise ValueError("Hermitian complement did not reduce to tridiagonal form")
        low = float(
            d[0]
            if self.nq == 1
            else eigh_tridiagonal(
                d,
                abs(b),
                eigvals_only=True,
                select="i",
                select_range=(0, 0),
                check_finite=False,
            )[0]
        )
        self.minimum_energy = low
        # Mixing two thermally empty states outside the DMFT window has no
        # density/projector readout. Its eigenvector derivative is undefined
        # at a degeneracy, but the physical derivative is regular. The caller
        # identifies columns irrelevant to the projector; independently check
        # their thermal occupation and slope before excluding any solve.
        self.active_columns = np.ones(len(e), bool) if active_columns is None else np.array(active_columns, copy=True)
        if self.active_columns.shape != e.shape or self.active_columns.dtype != bool or not np.any(self.active_columns):
            raise ValueError("invalid active complement columns")
        inactive = ~self.active_columns
        f_inactive = expit((self.mu-e[inactive])/self.temperature)
        self.inactive_charge_bound = float(2*np.sum(f_inactive))
        if self.inactive_charge_bound > charge_tolerance or self.inactive_charge_bound/self.temperature > slope_tolerance:
            raise ValueError("inactive retained columns are not thermally empty")
        self.minimum_gap_all = low - float(np.max(e))
        self.minimum_gap = low - float(np.max(e[self.active_columns]))
        if self.minimum_gap <= gap_tolerance:
            raise ValueError(
                "retained/complement boundary cuts a degeneracy or misses lower states; retain more bands"
            )
        # The factor two is the nonmagnetic spin convention of the electronic
        # model. Weighted sums over normalized k weights obey the same bound.
        f = float(expit((self.mu - low) / self.temperature))
        self.charge_bound = 2 * self.nq * f
        self.slope_bound = self.charge_bound / self.temperature
        if self.charge_bound > charge_tolerance or self.slope_bound > slope_tolerance:
            raise ValueError("complement is not thermally empty; retain more bands")
        safe_energies = np.where(self.active_columns, e, np.min(e[self.active_columns]))
        self.diagonal = d[:, None] - safe_energies[None, :]
        self.lower = np.empty((self.nq - 1, self.nr), complex)
        for i in range(self.nq - 1):
            if np.any(self.diagonal[i] <= 0):
                raise ValueError("nonpositive shifted complement pivot")
            self.lower[i] = b[i] / self.diagonal[i]
            self.diagonal[i + 1] -= abs(b[i]) ** 2 / self.diagonal[i]
        if np.any(self.diagonal[-1] <= 0):
            raise ValueError("nonpositive shifted complement pivot")
        self.preparation_seconds = time.perf_counter() - started
        for array in (self.h, self.energies, self.rotation, self.diagonal, self.lower, self.active_columns):
            array.setflags(write=False)

    def action(self, rhs):
        """Return T with H_Q T-T eps_R=-rhs; also the real-Frobenius transpose."""
        rhs = checked(rhs, (self.nq, self.nr), "complement right-hand side")
        x = -(self.rotation.conj().T @ (rhs * self.active_columns))
        for i in range(1, self.nq):
            x[i] -= self.lower[i - 1] * x[i - 1]
        x /= self.diagonal
        for i in range(self.nq - 2, -1, -1):
            x[i] -= self.lower[i].conj() * x[i + 1]
        return self.rotation @ x

    def residual(self, tangent, rhs):
        checked(tangent, (self.nq, self.nr), "complement tangent")
        checked(rhs, (self.nq, self.nr), "complement right-hand side")
        return float(
            np.linalg.norm((self.h @ tangent - tangent * self.energies + rhs) * self.active_columns)
            / max(np.linalg.norm(rhs * self.active_columns), np.finfo(float).tiny)
        )

    def info(self):
        return dict(
            complement_size=self.nq,
            retained_size=self.nr,
            minimum_gap_Ry=self.minimum_gap,
            retained_endpoint_gap_Ry=self.minimum_gap_all,
            inactive_retained_columns=int(np.sum(~self.active_columns)),
            inactive_retained_charge_bound=self.inactive_charge_bound,
            omitted_charge_bound=self.charge_bound,
            omitted_fermi_slope_bound_per_Ry=self.slope_bound,
            reduction_relative_error=self.reduction_error,
            preparation_seconds=self.preparation_seconds,
            stored_bytes=sum(
                a.nbytes for a in (self.h, self.rotation, self.diagonal, self.lower)
            ),
        )


@dataclass
class LAPWComplementBasis:
    """O-orthogonal complement from native H/O and ONLY retained eigenvectors.

    The native coefficient order must agree in H, O and C. Cholesky whitening
    and a complete QR supply Q; they do not solve for unoccupied KS vectors.
    Project raw overlaps by C^dagger phi_basis in this same gauge. The basis
    object's retained/complement arrays are fixed central data.
    """

    retained: np.ndarray
    complement: np.ndarray
    solver: ShiftedComplement
    checks: dict

    @classmethod
    def prepare(
        cls,
        hamiltonian,
        overlap,
        retained,
        energies,
        *,
        mu,
        temperature,
        eigen_tolerance=1e-8,
        **controls
    ):
        h, o, c, e = map(np.asarray, (hamiltonian, overlap, retained, energies))
        if (
            h.ndim != 2
            or h.shape[0] != h.shape[1]
            or c.ndim != 2
            or c.shape[0] != len(h)
        ):
            raise ValueError("invalid native coefficient dimensions")
        if not 0 < c.shape[1] < len(h):
            raise ValueError("need a nonempty retained space and complement")
        checked(h, h.shape, "native Hamiltonian", hermitian=True)
        checked(o, h.shape, "native overlap", hermitian=True)
        checked(c, c.shape, "retained coefficient matrix")
        checked(e, (c.shape[1],), "retained energies", real=True)
        if not np.isfinite(eigen_tolerance) or eigen_tolerance <= 0:
            raise ValueError("invalid eigenpair tolerance")
        lower = cholesky(o, lower=True)
        y = lower.conj().T @ c
        gram_error = float(np.linalg.norm(y.conj().T @ y - np.eye(len(e))))
        if gram_error > eigen_tolerance:
            raise ValueError("retained coefficients are not O-orthonormal")
        unitary, _ = qr(y, mode="full", check_finite=False)
        q = solve_triangular(
            lower.conj().T, unitary[:, len(e) :], lower=False, check_finite=False
        )
        # Check the full eigenpair residual in whitened coordinates, rather
        # than relative to a huge nuclear/high-energy Hamiltonian norm.
        residual = solve_triangular(
            lower, h @ c - (o @ c) * e, lower=True, check_finite=False
        )
        eigen_error = float(np.linalg.norm(residual) / max(1.0, np.linalg.norm(e)))
        if eigen_error > eigen_tolerance:
            raise ValueError("native H/O and retained eigenpairs are inconsistent")
        solver = ShiftedComplement(
            herm(q.conj().T @ h @ q), e, mu=mu, temperature=temperature, **controls
        )
        c = np.array(c, complex, copy=True)
        c.setflags(write=False)
        q.setflags(write=False)
        return cls(
            c,
            q,
            solver,
            dict(orthonormality_error=gram_error, eigenpair_relative_error=eigen_error),
        )

    def project_matrix(self, matrix):
        """Return retained RR and complement-retained QR blocks."""
        m = checked(
            matrix, (len(self.retained),) * 2, "coefficient matrix", hermitian=True
        )
        mc = m @ self.retained
        return herm(self.retained.conj().T @ mc), self.complement.conj().T @ mc

    def project_raw(self, raw):
        raw = np.asarray(raw)
        if (
            raw.ndim < 2
            or raw.shape[0] != len(self.retained)
            or not np.all(np.isfinite(raw))
        ):
            raise ValueError("invalid coefficient-space raw overlaps")
        shape = raw.shape[1:]
        flat = raw.reshape(len(raw), -1)
        return (
            (self.retained.conj().T @ flat).reshape(self.retained.shape[1], *shape),
            (self.complement.conj().T @ flat).reshape(self.complement.shape[1], *shape),
        )

    def density(self, rr, qr_block):
        """Contravariant density derivative in native LAPW coefficients."""
        checked(rr, (self.solver.nr,) * 2, "retained density", hermitian=True)
        checked(qr_block, (self.solver.nq, self.solver.nr), "cross density")
        cross = (self.complement @ qr_block) @ self.retained.conj().T
        return herm(
            self.retained @ rr @ self.retained.conj().T + cross + cross.conj().T
        )


class EmptyComplementResponse:
    """Add the static complement to a retained CachedFrozenSigmaResponse.

    raw_complement[k] is LOCAL impurity-block storage (Q,orbital,block), matching
    model.raw; it is not normalized. All supplied central data must share a
    gauge. The model includes spin and k weights; none are reapplied here.

    Output/input matrix block pairing is
       <Y,dR> = sum_k Re Tr(Y_RR^dagger dR_RR)
                         + 2 Re Tr(Y_QR^dagger dR_QR).
    The factor two matters in the transpose; raw overlaps have no such factor.
    QQ density and its derivative vanish under the checked empty-Q restriction.
    """

    def __init__(self, model, complements, raw_complement):
        self.model = model
        self.solvers = tuple(complements)
        nk, nr = model.eps.shape
        if len(self.solvers) != nk or len(raw_complement) != nk:
            raise ValueError("missing k-point complement")
        self.raw_q = []
        for k, solver in enumerate(self.solvers):
            if (
                solver.nr != nr
                or not np.array_equal(solver.energies, model.eps[k])
                or solver.mu != model.meta.mu
                or solver.temperature != model.meta.temperature
            ):
                raise ValueError("complement and retained central state differ")
            raw = np.array(
                checked(
                    raw_complement[k],
                    (solver.nq, *model.raw.shape[2:]),
                    "raw complement projector",
                ),
                complex,
                copy=True,
            )
            raw.setflags(write=False)
            self.raw_q.append(raw)
        for solver in self.solvers:
            if not np.all(solver.active_columns[model.W]):
                raise ValueError("DMFT window columns cannot be omitted from complement response")
        self.rho = model.central_density()
        self.rho.setflags(write=False)

    def response(self, hrr, orr, raw_r, hqr, oqr, *, fixed_number=True):
        if len(hqr) != len(self.solvers) or len(oqr) != len(self.solvers):
            raise ValueError("missing complement perturbation blocks")
        raw = np.array(
            checked(raw_r, self.model.raw.shape, "retained raw derivative"),
            complex,
            copy=True,
        )
        mixed = []
        for k, solver in enumerate(self.solvers):
            h = checked(hqr[k], (solver.nq, solver.nr), "cross Hamiltonian derivative")
            o = checked(oqr[k], h.shape, "cross overlap derivative")
            tangent = solver.action(h - o * solver.energies)
            mixed.append(tangent @ self.rho[k])
            raw[k] += (tangent.conj().T @ self.raw_q[k].reshape(solver.nq, -1)).reshape(
                raw[k].shape
            )
        rr, info = self.model.generalized_response(
            hrr, orr, raw, fixed_number=fixed_number
        )
        info = dict(info, complement=[solver.info() for solver in self.solvers])
        return rr, mixed, info

    def adjoint(self, yrr, yqr, *, fixed_number=True):
        if len(yqr) != len(self.solvers):
            raise ValueError("missing complement density covectors")
        hb, ob, rb = self.model.generalized_adjoint(yrr, fixed_number=fixed_number)
        hq, oq = [], []
        for k, solver in enumerate(self.solvers):
            y = checked(yqr[k], (solver.nq, solver.nr), "cross density covector")
            # Reverse delta phi_R += T^dagger phi_Q and delta rho_QR=T rho_R.
            tb = (
                self.raw_q[k].reshape(solver.nq, -1)
                @ rb[k].reshape(solver.nr, -1).conj().T
            )
            tb += 2 * y @ self.rho[k]
            b = solver.action(tb)
            # Recombine QR and RQ in a Hermitian full-matrix pairing.
            hq.append(b / 2)
            oq.append(-b * solver.energies / 2)
        return hb, ob, rb, hq, oq


class ComplementMatrixCoordinates:
    """Isometric real coordinates for Hermitian RR and rectangular QR blocks.

    QQ is identically zero in the density output and irrelevant in the H/O
    input for an empty complement. sqrt(2) on BOTH real/imaginary QR parts
    accounts for the conjugate RQ block. Complement sizes may vary with k.
    """

    def __init__(self, retained_shape, complement_sizes):
        self.retained = HermitianCoordinates(retained_shape)
        nk, nr, _ = self.retained.shape
        if len(complement_sizes) != nk or any(
            int(n) != n or n < 1 for n in complement_sizes
        ):
            raise ValueError("invalid k-resolved complement sizes")
        self.shapes = [(int(n), nr) for n in complement_sizes]
        self.size = self.retained.size + 2 * sum(n * nr for n in complement_sizes)

    def pack(self, rr, qr_blocks):
        if len(qr_blocks) != len(self.shapes):
            raise ValueError("missing cross blocks")
        result = [self.retained.pack(rr)]
        for value, shape in zip(qr_blocks, self.shapes):
            value = checked(value, shape, "cross matrix coordinates")
            result.extend(
                [np.sqrt(2) * value.real.ravel(), np.sqrt(2) * value.imag.ravel()]
            )
        return np.concatenate(result)

    def unpack(self, vector):
        v = checked(vector, (self.size,), "complement matrix coordinates", real=True)
        pos = self.retained.size
        rr = self.retained.unpack(v[:pos])
        blocks = []
        for shape in self.shapes:
            n = int(np.prod(shape))
            blocks.append(
                (
                    (v[pos : pos + n] + 1j * v[pos + n : pos + 2 * n]) / np.sqrt(2)
                ).reshape(shape)
            )
            pos += 2 * n
        return rr, blocks


class ComplementPerturbationCoordinates:
    """Joint H/O/raw-R coordinates, with complete rectangular complement legs."""

    def __init__(self, matrix_coordinates, raw_shape):
        self.matrix = matrix_coordinates
        self.raw_shape = tuple(raw_shape)
        if (
            len(self.raw_shape) != 4
            or min(self.raw_shape) < 1
            or self.raw_shape[:2] != self.matrix.retained.shape[:2]
        ):
            raise ValueError("invalid retained raw-overlap shape")
        self.raw_size = int(np.prod(raw_shape))
        self.size = 2 * self.matrix.size + 2 * self.raw_size

    def pack(self, hrr, orr, raw_r, hqr, oqr):
        raw = checked(raw_r, self.raw_shape, "raw retained coordinates")
        return np.concatenate(
            [
                self.matrix.pack(hrr, hqr),
                self.matrix.pack(orr, oqr),
                raw.real.ravel(),
                raw.imag.ravel(),
            ]
        )

    def unpack(self, vector):
        v = checked(vector, (self.size,), "joint complement coordinates", real=True)
        n, r = self.matrix.size, self.raw_size
        hrr, hqr = self.matrix.unpack(v[:n])
        orr, oqr = self.matrix.unpack(v[n : 2 * n])
        raw = (v[2 * n : 2 * n + r] + 1j * v[2 * n + r :]).reshape(self.raw_shape)
        return hrr, orr, raw, hqr, oqr


def complement_electronic_operator(model, *, fixed_number=True):
    """Return (paired E, input coordinates, density coordinates) for composition.

    Unlike a full-square representation, these interfaces carry no QQ block.
    A spatial-density adapter must reconstruct BOTH RR and QR/RQ contributions,
    and add the direct moving-basis map separately as in R E B + D.
    """
    if hasattr(model, "electronic_operator"):
        return model.electronic_operator(fixed_number=fixed_number)
    nk, nr = model.model.eps.shape
    output = ComplementMatrixCoordinates((nk, nr, nr), [s.nq for s in model.solvers])
    inputs = ComplementPerturbationCoordinates(output, model.model.raw.shape)

    def forward(vector):
        rr, qr_blocks, _ = model.response(
            *inputs.unpack(vector), fixed_number=fixed_number
        )
        return output.pack(rr, qr_blocks)

    def reverse(vector):
        return inputs.pack(
            *model.adjoint(*output.unpack(vector), fixed_number=fixed_number)
        )

    return (
        PairedOperator(
            inputs.size,
            output.size,
            forward,
            reverse,
            "retained DMFT response with empty LAPW complement",
        ),
        inputs,
        output,
    )
