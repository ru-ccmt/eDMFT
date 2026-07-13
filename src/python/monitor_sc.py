#!/usr/bin/env python
"""
monitor_sc.py - Read-only monitoring tool for superconducting convergence diagnostics.

Reads a unified Nambu sig.inp file (8K+1 columns) and computes SC convergence quantities:
- Renormalization factor Z
- Anomalous self-energy at lowest Matsubara frequency
- Estimated SC gap
- Maximum anomalous self-energy magnitude
- SC convergence status
- Particle-hole symmetry violation

Column layout of unified sig.inp (K = number of independent normal components):
  col 0:       omega (Matsubara frequency)
  cols 1..2K:  Sigma_pp (normal particle self-energy)
  cols 2K+1..4K: Phi (anomalous self-energy, upper off-diagonal)
  cols 4K+1..6K: Phi_dag (anomalous self-energy, lower off-diagonal)
  cols 6K+1..8K: Sigma_hh (hole self-energy, = -conj(Sigma_pp) by symmetry)
"""

import argparse
import numpy as np
import sys


def read_unified_siginp(filename):
    """
    Read unified Nambu sig.inp file (8K+1 columns).

    Returns:
        omega: 1D array of Matsubara frequencies
        K: number of independent normal components
        sigma_pp: complex array (n_freq x K) — normal self-energy
        phi: complex array (n_freq x K) — anomalous self-energy
        sigma_hh: complex array (n_freq x K) — hole self-energy
    """
    try:
        raw_data = np.loadtxt(filename, skiprows=2)
    except IOError as e:
        print(f"Error: Cannot read file '{filename}': {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: Invalid data format in '{filename}': {e}", file=sys.stderr)
        sys.exit(1)

    if raw_data.ndim == 1:
        raw_data = raw_data.reshape(1, -1)

    ncols = raw_data.shape[1]
    K = (ncols - 1) // 8

    if 8 * K + 1 != ncols:
        print(f"Error: sig.inp has {ncols} columns, not a valid 8K+1 Nambu format.", file=sys.stderr)
        print(f"  Expected 8K+1 columns (K must be integer). Got remainder={(ncols-1)%8}", file=sys.stderr)
        sys.exit(1)

    omega = raw_data[:, 0]
    data = raw_data[:, 1:]  # 8K data columns

    def extract_complex(block_data):
        """Extract complex array from Re/Im column pairs."""
        re_part = block_data[:, 0::2]
        im_part = block_data[:, 1::2]
        return re_part + 1j * im_part

    sigma_pp = extract_complex(data[:, 0:2*K])
    phi = extract_complex(data[:, 2*K:4*K])
    phi_dag = extract_complex(data[:, 4*K:6*K])
    sigma_hh = extract_complex(data[:, 6*K:8*K])

    return omega, K, sigma_pp, phi, phi_dag, sigma_hh


def compute_sc_diagnostics(sig_file, threshold):
    """
    Compute superconducting convergence diagnostics from unified sig.inp.

    Args:
        sig_file: Path to unified Nambu sig.inp file (8K+1 columns)
        threshold: SC detection threshold (eV)

    Returns:
        Dictionary with diagnostic quantities
    """
    omega, K, sigma_pp, phi, phi_dag, sigma_hh = read_unified_siginp(sig_file)

    iw0 = 0

    z_factor = 1.0 - sigma_pp[iw0, 0].imag / omega[iw0]

    phi_0_re = phi[iw0, 0].real
    phi_0_im = phi[iw0, 0].imag

    if abs(z_factor) > 1e-10:
        gap_estimate = phi_0_re / z_factor
    else:
        gap_estimate = 0.0

    phi_magnitude = np.abs(phi)
    max_phi_magnitude = np.max(phi_magnitude)
    max_phi_re = np.max(np.abs(phi.real))

    # Particle-hole symmetry check: Sigma_hh should be -conj(Sigma_pp)
    ph_violation = np.max(np.abs(sigma_hh + np.conj(sigma_pp)))

    # Hermitian conjugate check: Phi_dag should be conj(Phi)
    hc_violation = np.max(np.abs(phi_dag - np.conj(phi)))

    if max_phi_re < 1e-5:
        status = "No superconductivity (Phi -> 0)"
    elif max_phi_re > threshold:
        status = "Superconducting state found!"
    else:
        status = "Still converging..."

    return {
        'K': K,
        'z_factor': z_factor,
        'phi_0_re': phi_0_re,
        'phi_0_im': phi_0_im,
        'gap_estimate': gap_estimate,
        'max_phi_magnitude': max_phi_magnitude,
        'max_phi_re': max_phi_re,
        'ph_violation': ph_violation,
        'hc_violation': hc_violation,
        'status': status,
        'omega_0': omega[iw0]
    }


def print_diagnostics(diagnostics):
    """Print formatted diagnostic output."""
    print("\n" + "="*60)
    print("SUPERCONDUCTING CONVERGENCE DIAGNOSTICS")
    print("="*60)
    print(f"Unified sig.inp: K={diagnostics['K']} ({8*diagnostics['K']+1} columns)")
    print(f"Lowest Matsubara frequency: {diagnostics['omega_0']:.6f} eV")
    print(f"Renormalization factor Z: {diagnostics['z_factor']:.6f}")
    print(f"Anomalous self-energy at omega_0:")
    print(f"  Phi(omega_0) = {diagnostics['phi_0_re']:+.6f} + {diagnostics['phi_0_im']:+.6f}i eV")
    print(f"Estimated SC gap Delta: {diagnostics['gap_estimate']:.6f} eV")
    print(f"Maximum |Phi(omega)| across all frequencies: {diagnostics['max_phi_magnitude']:.6f} eV")
    print(f"Maximum |Phi_re(omega)|: {diagnostics['max_phi_re']:.6f} eV")
    print(f"Symmetry checks:")
    print(f"  PH violation (Sigma_hh + conj(Sigma_pp)): {diagnostics['ph_violation']:.2e}")
    print(f"  HC violation (Phi_dag - conj(Phi)):       {diagnostics['hc_violation']:.2e}")
    print(f"\nSTATUS: {diagnostics['status']}")
    print("="*60 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description='Monitor superconducting convergence from unified Nambu sig.inp.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python monitor_sc.py
  python monitor_sc.py --sig sig.inp --threshold 0.005
        """
    )

    parser.add_argument(
        '-s', '--sig',
        default='sig.inp',
        help='Path to unified Nambu sig.inp file (default: sig.inp)'
    )
    parser.add_argument(
        '--threshold',
        type=float,
        default=0.001,
        help='SC detection threshold in eV (default: 0.001)'
    )

    args = parser.parse_args()

    diagnostics = compute_sc_diagnostics(args.sig, args.threshold)
    print_diagnostics(diagnostics)


if __name__ == '__main__':
    main()
