#!/usr/bin/env python
"""Shared scalar cluster-Hubbard interaction helpers.

The eDMFT impurity parameter dictionaries store values as ``[value, comment]``.
This module resolves the solver-specific legacy spellings without letting those
spellings change the physical model.  In particular, scalar ``Uc`` means the
literal active-model interaction ``Uc * sum_i n_i_up n_i_down``; it is not a
Slater ``F0`` and is never combined with ``U`` or ``J``.
"""

from dataclasses import dataclass
import hashlib
import json
import math
import numbers
import warnings

import numpy


_MISSING = object()
_DENSITY_DENSITY_FORMS = frozenset(("Ising", "IsingS", "IsingK"))


def parse_parameter_value(params, key, default=_MISSING):
    """Return the value component of an eDMFT impurity parameter.

    Parameters are normally stored as ``[value, comment]``.  A bare scalar is
    accepted as well so the helpers remain convenient in focused unit tests.
    """
    if key not in params:
        if default is _MISSING:
            raise KeyError(key)
        return default

    entry = params[key]
    if isinstance(entry, (list, tuple)):
        if not entry:
            raise ValueError("empty eDMFT parameter entry for %s" % key)
        return entry[0]
    return entry


def _finite_real(value, key):
    if isinstance(value, (bool, numpy.bool_)) or not isinstance(value, numbers.Real):
        raise ValueError("%s must be a finite real scalar, got %r" % (key, value))
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("%s must be a finite real scalar, got %r" % (key, value))
    # Canonicalize signed zero for stable logging and hashing.
    return 0.0 if result == 0.0 else result


def _coulomb_form(params):
    value = parse_parameter_value(params, "CoulombF", "Ising")
    return str(value).strip().strip("'").strip('"')


@dataclass(frozen=True)
class InteractionSpec:
    """Resolved impurity interaction and its provenance."""

    mode: str
    value: object
    source: object
    coulomb_form: str
    solver: str
    cluster: bool
    legacy: bool = False

    @property
    def Uc(self):
        """Compatibility spelling used by the impurity wrappers."""
        return self.value

    @property
    def model(self):
        return self.mode

    def physical_payload(self):
        """Return only fields that determine the scalar interaction itself."""
        payload = {"schema": 1, "mode": self.mode}
        if self.mode == "single_band_hubbard":
            payload["Uc"] = self.value
            # All accepted Ising-family labels serialize to this same model.
            payload["interaction_form"] = "density_density"
        elif self.coulomb_form:
            payload["interaction_form"] = self.coulomb_form
        return payload

    def to_payload(self):
        """Return the physical fields together with diagnostic provenance."""
        payload = self.physical_payload()
        payload.update({
            "source": self.source,
            "solver": self.solver,
            "cluster": self.cluster,
            "legacy": self.legacy,
            "requested_coulomb_form": self.coulomb_form,
        })
        return payload


def resolve_interaction_spec(params, solver, cluster):
    """Resolve ``Uc`` and solver-specific legacy behavior.

    ``Uc`` always selects the literal single-band cluster-Hubbard model.
    BCTQMC additionally accepts legacy ``Ud`` when ``Uc`` is absent.  iNORG
    without ``Uc`` keeps its existing Slater construction unchanged.
    """
    solver_name = str(solver)
    solver_key = solver_name.lower()
    is_cluster = bool(cluster)
    form = _coulomb_form(params)

    if "Uc" in params:
        if not is_cluster:
            raise ValueError(
                "%s: scalar Uc is supported only for a cluster impurity"
                % solver_name)
        value = _finite_real(parse_parameter_value(params, "Uc"), "Uc")
        if form not in _DENSITY_DENSITY_FORMS:
            raise ValueError(
                "%s: scalar Uc supports only CoulombF=%s; got %r"
                % (solver_name, sorted(_DENSITY_DENSITY_FORMS), form))
        return InteractionSpec(
            mode="single_band_hubbard",
            value=value,
            source="Uc",
            coulomb_form=form,
            solver=solver_name,
            cluster=True,
            legacy=False,
        )

    if solver_key == "bctqmc" and is_cluster:
        if "Ud" not in params:
            raise ValueError(
                "BCTQMC cluster impurity requires Uc "
                "(or legacy Ud for backward compatibility)")
        value = _finite_real(parse_parameter_value(params, "Ud"), "Ud")
        warnings.warn(
            "BCTQMC cluster interaction from Ud is deprecated; use Uc. "
            "Ud retains its existing double-counting override meaning.",
            FutureWarning,
            stacklevel=2,
        )
        return InteractionSpec(
            mode="single_band_hubbard",
            value=value,
            source="Ud",
            coulomb_form=form,
            solver=solver_name,
            cluster=True,
            legacy=True,
        )

    if solver_key == "inorg":
        return InteractionSpec(
            mode="legacy_slater",
            value=None,
            source="U/J",
            coulomb_form=form,
            solver=solver_name,
            cluster=is_cluster,
            legacy=True,
        )

    # The BCTQMC single-site path continues to use its existing cix/runtime-F0
    # division.  It has no scalar cluster interaction to resolve.
    if solver_key == "bctqmc":
        return InteractionSpec(
            mode="legacy_solver",
            value=None,
            source="U/J",
            coulomb_form=form,
            solver=solver_name,
            cluster=False,
            legacy=True,
        )

    raise ValueError("unsupported impurity solver %r" % solver_name)


def build_single_band_cluster_tensor(spec, Sigind, mem_indx, nat, spin_deg):
    """Build a literal onsite Hubbard tensor in the physical site basis.

    Exactly one correlated physical orbital must be present in every atom's
    global-index mapping.  For doubled spin/Nambu layouts (``spin_deg == 2``),
    only the first, physical block is returned; the serializer handles the
    particle-hole representation downstream.
    """
    if spec.mode != "single_band_hubbard" or spec.value is None:
        raise ValueError("a resolved single_band_hubbard InteractionSpec is required")

    sigind = numpy.asarray(Sigind)
    if sigind.ndim != 2 or sigind.shape[0] != sigind.shape[1]:
        raise ValueError("Sigind must be a square matrix")
    if spin_deg not in (1, 2):
        raise ValueError("spin_deg must be 1 or 2, got %r" % (spin_deg,))
    if spin_deg == 2 and sigind.shape[0] % 2:
        raise ValueError("a doubled Sigind must have even dimension")

    nat = int(nat)
    if nat <= 0:
        raise ValueError("cluster must contain at least one atom")
    nphys = sigind.shape[0] // 2 if spin_deg == 2 else sigind.shape[0]
    tensor = numpy.zeros((nphys, nphys, nphys, nphys), dtype=complex)
    active_global = []

    for iat in range(nat):
        if iat not in mem_indx:
            raise ValueError("missing orbital mapping for cluster atom %d" % iat)
        candidates = []
        for raw_index in mem_indx[iat]:
            index = int(raw_index)
            if 0 <= index < nphys and sigind[index, index] != 0:
                if index not in candidates:
                    candidates.append(index)
        if len(candidates) != 1:
            raise ValueError(
                "scalar Uc requires exactly one active physical orbital per "
                "cluster atom; atom %d has %d (%r)"
                % (iat, len(candidates), candidates))
        active_global.append(candidates[0])

    if len(set(active_global)) != nat:
        raise ValueError(
            "cluster atoms must map to distinct active physical orbitals: %r"
            % active_global)

    for index in active_global:
        tensor[index, index, index, index] = spec.value
    return tensor


def _canonicalize(value):
    if isinstance(value, InteractionSpec):
        return _canonicalize(value.physical_payload())
    if isinstance(value, numpy.ndarray):
        return _canonicalize(value.tolist())
    if isinstance(value, numpy.generic):
        return _canonicalize(value.item())
    if isinstance(value, dict):
        return {
            str(key): _canonicalize(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (list, tuple)):
        return [_canonicalize(item) for item in value]
    if isinstance(value, (set, frozenset)):
        items = [_canonicalize(item) for item in value]
        return sorted(items, key=lambda item: json.dumps(item, sort_keys=True))
    if isinstance(value, complex):
        real = _finite_real(value.real, "fingerprint complex real part")
        imag = _finite_real(value.imag, "fingerprint complex imaginary part")
        return {"real": real, "imag": imag}
    if isinstance(value, numbers.Real) and not isinstance(value, bool):
        if isinstance(value, numbers.Integral):
            return int(value)
        return _finite_real(value, "fingerprint float")
    if value is None or isinstance(value, (str, bool)):
        return value
    raise TypeError("unsupported fingerprint value %r" % (value,))


def interaction_fingerprint(spec, extra=None):
    """Return a stable SHA-256 digest of physical semantics and structure.

    Provenance such as whether the value was spelled ``Uc`` or legacy ``Ud``
    is deliberately excluded.  Wrappers add cluster shape, basis, and mode via
    ``extra``.
    """
    payload = {"interaction": spec.physical_payload()}
    if extra is not None:
        payload["extra"] = extra
    encoded = json.dumps(
        _canonicalize(payload),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


__all__ = [
    "InteractionSpec",
    "build_single_band_cluster_tensor",
    "interaction_fingerprint",
    "parse_parameter_value",
    "resolve_interaction_spec",
]
