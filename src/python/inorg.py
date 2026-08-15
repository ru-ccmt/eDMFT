#!/usr/bin/env python
#
import sys, re, os, glob, shutil, types
from scipy import interpolate
from utils import DmftEnvironment
from cmpEimp3 import GiveImpFunctional
import brod
import subprocess
import numpy
import dmft_DC
import gaunt
import scipy
import dpybind
from numpy import zeros
from inorg_edmft_io import plaquette_unitary
from cluster_interaction import (
    build_single_band_cluster_tensor,
    resolve_interaction_spec,
)


def _parse_nmat_txt(nmat_path, nambu=False, fh_info=None):
    """Parse nmat.txt from iNORG solver and return occupations.

    Normal-state keys: up1, up2, ..., dw1, dw2, ..., sum
    Nambu keys: K0p, K1p, ..., K0h, K1h, ..., s_p, s_h, sum

    Args:
        nmat_path: Path to nmat.txt file.
        nambu: If True, parse Nambu-format keys and convert hole occupancies.
        fh_info: File handle for logging (optional).

    Returns:
        (nf, mom, s_p, s_h) where:
            nf   — total physical electron count (float)
            mom  — list of per-K (or per-orbital) total occupancies
            s_p  — sum of particle occupancies (None if normal-state)
            s_h  — sum of hole occupancies (None if normal-state)
    """
    mom_dict = {}
    nf = 0.0
    s_p = None
    s_h = None

    with open(nmat_path, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) < 2 or line.strip().startswith('#'):
                continue

            key, val_str = parts[0], parts[1]

            if key == 'sum':
                nf = float(val_str)
            # Normal-state keys: up1, dw1, up2, dw2, ...
            elif (key.startswith('up') or key.startswith('dw')) and key[2:].isdigit():
                spin = 0 if key.startswith('up') else 1
                orb_idx = int(key[2:]) - 1
                mom_dict[(spin, orb_idx)] = float(val_str)
            # Nambu keys: K0p, K1p, ..., K0h, K1h, ...
            elif len(key) >= 3 and key[0] == 'K' and key[-1] in ('p', 'h') and key[1:-1].isdigit():
                orb_idx = int(key[1:-1])
                spin = 0 if key[-1] == 'p' else 1  # p=particle(up), h=hole(down)
                mom_dict[(spin, orb_idx)] = float(val_str)
            elif key == 's_p':
                s_p = float(val_str)
            elif key == 's_h':
                s_h = float(val_str)

    mom = []
    if mom_dict:
        orb_indices = set(orb_idx for (spin, orb_idx) in mom_dict.keys())
        max_orb = max(orb_indices) if orb_indices else -1
        norb_K = max_orb + 1

        if nambu and s_p is not None and s_h is not None:
            # Nambu: convert hole occupancies to physical electron occupancies
            # n_{K,down} = 1 - ñ_{K,down}  where ñ is the hole occupation (Kh value)
            for orb_idx in range(norb_K):
                if (1, orb_idx) in mom_dict:
                    mom_dict[(1, orb_idx)] = 1.0 - mom_dict[(1, orb_idx)]
            # Physical nf = s_p + (norb_K - s_h)
            nf = s_p + (norb_K - s_h)

        # Build mom list: total occupation per orbital = up + down
        for orb_idx in range(norb_K):
            up_occ = mom_dict.get((0, orb_idx), 0.0)
            down_occ = mom_dict.get((1, orb_idx), 0.0)
            mom.append(up_occ + down_occ)

    return (nf, mom, s_p, s_h)



def merge_nambu_sigout(sig_path, phi_path, fh_info=None):
    """Merge sig.out + phi.out into unified Nambu sig.out in sigind order.

    Reads:
      sig.out: 2K+1 cols — ω + K×(Re,Im) for Σ_pp (Edc-subtracted)
      phi.out: 2K+1 cols — ω + K×(Re,Im) for Φ (no Edc)

    Writes unified sig_path with 8K+1 cols in sigind order:
      ω | Σ_pp_1 Σ_pp_2 .. Σ_pp_K | Φ_1 .. Φ_K | Φ†_1 .. Φ†_K | Σ_hh_1 .. Σ_hh_K

    The sigind numbering matches the Nambu indmfl:
      sigind 1..K     = Σ_pp components
      sigind K+1..2K  = Φ components
      sigind 2K+1..3K = Φ† components
      sigind 3K+1..4K = Σ_hh components

    Headers s_oo/Edc have 4K elements in the same sigind order.
    sgather.py reads columns sequentially by sigind, so this layout
    ensures correct mapping without modifying sgather.

    Sentinel row at ω=5000.0 with high-frequency limits.
    If phi_path does not exist, returns without action (non-Nambu).
    """
    if not os.path.exists(phi_path):
        return

    # --- Read sig.out header and data ---
    with open(sig_path, 'r') as f:
        hdr1 = f.readline()   # # s_oo= [...]
        hdr2 = f.readline()   # # Edc= [...]
    s_oo_pp = eval(hdr1.split('=', 1)[1].strip())
    Edc_pp  = eval(hdr2.split('=', 1)[1].strip())
    K = len(s_oo_pp)
    # Safety: Edc may have more elements than s_oo (e.g., 12 Nambu Edc from Eimp.inp
    # vs 3 pp-block data components). Truncate to match K = len(s_oo_pp).
    Edc_pp = list(Edc_pp)[:K]

    sig_data = numpy.loadtxt(sig_path, comments='#')   # (nfreq, 2K+1)
    phi_data = numpy.loadtxt(phi_path, comments='#')   # (nfreq, 2K+1)

    omega    = sig_data[:, 0]
    sigma_pp = sig_data[:, 1:]     # (nfreq, 2K) — Edc-subtracted Re
    phi_cols = phi_data[:, 1:]     # (nfreq, 2K) — no Edc

    # --- Construct Φ† and Σ_hh blocks ---
    # Φ†(iω) = Φ*(iω): Re unchanged, Im negated
    phi_dag = phi_cols.copy()
    phi_dag[:, 1::2] *= -1

    # Σ_hh(iω) = -conj(Σ_pp(iω))
    # Stored format: Re col = Re(Σ)-Edc, Im col = Im(Σ)
    # With Edc_hh = -Edc_pp:
    #   Re(Σ_hh - Edc_hh) = -(Re(Σ_pp) - Edc_pp)  → negate Re cols
    #   Im(Σ_hh) = Im(Σ_pp)                         → keep Im cols
    sigma_hh = sigma_pp.copy()
    sigma_hh[:, 0::2] *= -1

    # --- Physics verification: Σ_hh should be -conj(Σ_pp) by construction ---
    pp_cplx = sigma_pp[:, 0::2] + 1j * sigma_pp[:, 1::2]
    hh_cplx = sigma_hh[:, 0::2] + 1j * sigma_hh[:, 1::2]
    ph_viol = numpy.max(numpy.abs(hh_cplx + numpy.conj(pp_cplx)))
    if ph_viol > 1e-10 and fh_info is not None:
        print(f"WARNING: merge_nambu_sigout PH violation = {ph_viol:.2e}", file=fh_info)

    # --- Unified headers in sigind order (4K elements) ---
    # sigind 1..K: Σ_pp,  K+1..2K: Φ,  2K+1..3K: Φ†,  3K+1..4K: Σ_hh
    s_oo = list(s_oo_pp) + [0.0]*K + [0.0]*K + [-x for x in s_oo_pp]
    Edc  = list(Edc_pp)  + [0.0]*K + [0.0]*K + [-x for x in Edc_pp]

    # --- Interleave blocks in sigind order ---
    # Each block has K components × 2 columns (Re, Im).
    # Output order: pp_1 pp_2 .. pp_K  phi_1 .. phi_K  phiD_1 .. phiD_K  hh_1 .. hh_K
    # This matches sigind numbering so sgather reads columns sequentially.
    merged_data = numpy.column_stack([sigma_pp, phi_cols, phi_dag, sigma_hh])
    # merged_data shape: (nfreq, 8K)

    # --- Write unified sig.out (8K+1 cols, 4K header elements) ---
    with open(sig_path, 'w') as f:
        print('# s_oo=', s_oo, file=f)
        print('# Edc=', Edc, file=f)
        for iom in range(len(omega)):
            print(("%20.15f " % omega[iom]), end=' ', file=f)
            for i in range(4 * K):
                print(("%20.15f %20.15f  " % (merged_data[iom, 2*i],
                                               merged_data[iom, 2*i + 1])),
                      end=' ', file=f)
            print(file=f)
        # Sentinel row: high-frequency limit (s_oo[i]-Edc[i] for Re, 0 for Im)
        print(("%20.15f " % 5000.0), end=' ', file=f)
        for i in range(4 * K):
            print(("%20.15f %20.15f  " % (s_oo[i] - Edc[i], 0.0)),
                  end=' ', file=f)
        print(file=f)

    # Remove phi.out — now embedded in unified sig.out
    os.remove(phi_path)

    if fh_info is not None:
        print(f'Nambu merge: sig.out + phi.out -> unified sig.out '
              f'(K={K}, {8*K+1} cols, 4K={4*K} sigind components, '
              f'{len(omega)} freqs + sentinel)',
              file=fh_info)



def extract_impurity_matrix_from_outputdmf1(outputdmf1_file):
    """Extract impurity levels matrix from case.outputdmf1 file with adaptive dimensions"""
    impurity_matrices = {}
    
    try:
        with open(outputdmf1_file, 'r') as f:
            # Search for "Full matrix of impurity levels follows" marker
            found = False
            for line in f:
                if line.strip() == 'Full matrix of impurity levels follows':
                    found = True
                    break
            
            if not found:
                return None
            
            # Read matrices for each icix
            while True:
                line = f.readline()
                if not line:
                    break
                
                # Look for icix= line
                if 'icix=' in line:
                    parts = line.split()
                    icix = int(parts[1])
                    
                    # Read all matrix rows for this icix
                    matrix_rows = []
                    while True:
                        line = f.readline()
                        if not line:
                            break
                        
                        # Skip empty lines
                        line = line.strip()
                        if not line:
                            continue
                            
                        # Check if we reached the next section
                        if 'icix=' in line or line.startswith('#') or 'omega' in line:
                            break
                        
                        try:
                            row_data = [float(x) for x in line.split()]
                            if len(row_data) > 0:  # Valid data row
                                matrix_rows.append(row_data)
                        except (ValueError, IndexError):
                            continue  # Skip unparseable lines
                    
                    if matrix_rows:
                        # Determine matrix dimensions from the data
                        first_row = matrix_rows[0]
                        dim = len(first_row) // 2  # Each complex number takes 2 floats
                        expected_rows = dim
                        
                        if len(matrix_rows) >= expected_rows:
                            # Initialize matrix
                            matrix = zeros((dim, dim), dtype=complex)
                            
                            # Fill matrix from all rows
                            for i in range(min(expected_rows, len(matrix_rows))):
                                row_data = matrix_rows[i]
                                if len(row_data) >= 2 * dim:
                                    for j in range(dim):
                                        matrix[i, j] = row_data[2*j] + 1j * row_data[2*j + 1]
                            
                            impurity_matrices[icix] = matrix
                            print(f"Successfully extracted {dim}x{dim} matrix for icix={icix}")
    
    except (IOError, OSError) as e:
        print(f"Warning: Could not read {outputdmf1_file}: {e}")
        return None
    
    return impurity_matrices

Default_CoulombF = 'Ising'


def _resolve_inorg_sym(nambu, params):
    """Return a solver symmetry token consistent with the execution mode."""
    requested = None
    if 'sym' in params:
        requested = str(params['sym'][0]).strip("'\"")

    if nambu:
        allowed = {'psc', 'dsc', 'K-dsc'}
        if requested is None:
            return 'K-dsc'
        if requested not in allowed:
            raise ValueError(
                "iNORG: Nambu mode requires sym in "
                f"{sorted(allowed)}, got {requested!r}")
        return requested

    allowed = {'normal-pm', 'normal-af'}
    if requested is None:
        return 'normal-pm'
    if requested not in allowed:
        raise ValueError(
            "iNORG: normal mode requires sym in "
            f"{sorted(allowed)}, got {requested!r}")
    return requested


def _validate_inorg_interaction_spec(spec, coulomb_form):
    """Reject interaction models that the iNORG serializer cannot represent."""
    if (spec.mode == 'single_band_hubbard'
            and coulomb_form not in ['Ising', 'IsingS', 'IsingK']):
        raise ValueError(
            "iNORG: Uc selects the single-band cluster Hubbard model, "
            "which supports only the Ising CoulombF family; got "
            + repr(coulomb_form))

class IMP_INORG(object):
    """
    Wrapper for inorg
    """
    def __init__(self, Znuc, dire, params, matsubara, Sigind, CF, l, nambu=False):
        # Nambu flag — stored before anything else so it is always available
        self.nambu = nambu

        # Environment
        self.env = DmftEnvironment()
        self.mpi_prefix = self.env.MPI

        # impurity directory
        self.dir = dire+'/'

        # create impurity directory if it does not yet exist
        if not os.path.isdir(dire):
            os.mkdir( dire )

        # log-file is opened
        self.fh_info = open(self.dir + 'info.out', 'w')

        # Check type of interaction
        if 'CoulombF' in params:
            self.CoulombF = params['CoulombF'][0].strip("'").strip('"')
        else:
            params['CoulombF']=["'"+Default_CoulombF+"'", "# Full Coulomb repulsion"]
            self.CoulombF = Default_CoulombF

        # Check nbr of impurity problem
        if 'icase' in params:
            self.icase = params['icase'][0]
        # Load UlamJ values if there are some
        if 'UlamJ' in params:
            self.UlamJ = params['UlamJ'][0]

        # Copy the correlated indices
        self.Sigind = numpy.array(Sigind)
        # Crytal field splitting
        self.CF = CF
        # Real or imaginary axis
        self.matsubara = matsubara
        # Spin polarized or not
        if 'spin_dependent' in params:
            self.spin_dep = params['spin_dependent'][0]
        else:
            self.spin_dep = False

        # find smallest sigind to shift siginds so they start from 1
        minsigind = min([s for s in self.Sigind.flat if s != 0])
        self.Sigind_shifted = numpy.where(self.Sigind != 0, self.Sigind-minsigind+1, 0)

        # Transformation between spheric harmonics and local coordinates needs to be
        # used in exact diagonalization of the atom.
        with open(self.dir + 'Trans.dat', 'w') as fh_T:
            print(len(Sigind), len(Sigind[0]), '#  size of Sigind and CF', file=fh_T)
            print('#---- Sigind follows', file=fh_T)
            max_sigfig = 1 + int(numpy.log(max(self.Sigind_shifted.flat))/numpy.log(10))
            fmt = '{:'+str(max_sigfig)+'d} '
            for row in self.Sigind_shifted:
                print((fmt*len(row)).format(*row), file=fh_T)
            print('#---- CF follows', file=fh_T)
            for row in CF:
                for elem in row:
                    print('{:12.8f} {:12.8f} '.format(elem.real,elem.imag), end=' ', file=fh_T)
                print(file=fh_T)

        # parameters, which will be given to INORG and will be written to PARAMS
        # Need to make difference between real or imaginary frequency
        if self.matsubara:
            _nom_ = 100
            if 'beta' in params: _nom_ = int(3*params['beta'][0])
            self.sparams={'exe':'inorg_run.py', 'Sig':'Sig.out', 'Delta':'Delta.inp', 'Gf':'Gf.out', 'EDinput':'input.tsv',
                          'old_bathsec':'old_bathsec.out', 'nom':_nom_, 'cluster':False, 'force_szero':None, 'hyb_optmet':'PRAXIS'}
        else:
            raise ValueError("INORG: Do not use INORG solver on the real frequency axis in the self-consistent"
                             +" calculation, please switch to the Matsubara axis.")
        self.PARAMS = 'PARAMS.norg'


        ####### defining l and J
        self.l = l                    # -> for constructing Coulomb matrix
        # Quit if f ORBITALS -> not doable with ED SOLVER
        if (self.l >= 3):
            raise ValueError("INORG: INORG solver is not suited for f shell and beyond. Consider using CTQMC..")
            sys.exit(1)

        # UlamJ= {0: (U, lmbda, epsilon, [J2, J4])}
        self.J = params['J'][0]
        if 'UlamJ' in params and self.UlamJ[1]>1e-5:  # DCs='exact'
            # If we use Yukawa screening, we need to change Jhunds to be compatible with Yukawa form of screening
            # For dielectric screening we do not change Jhunds from user specified values
            Jh = self.UlamJ[3]
            self.J = Jh

        if self.nambu:
            print('IMP_INORG: Nambu superconducting mode enabled (sym=psc)', file=self.fh_info)

        #########################



    def _exe(self, params, DCs, extn, cix, maxAc=200.):
        """
        Executes the INORG impurity solver.

        Nambu additions:
          - Forces sym='psc' in PARAMS.norg when self.nambu is True.
          - Injects bin-nambu/ into PYTHONPATH so inorg_run.py can import
            the Nambu-aware inorg_edmft_io.
          - Verifies that Phi.out appears after the solver returns.
        """

        ##### Reading input, preparing impurity levels  ##########################################
        REMOVE_STATUS=False
        Qaverage = False
        # Reading impurity levels
        (Eimp,Olap,Edc) = numpy.loadtxt(self.dir+'/Eimp.inp',ndmin=2)
        # Reading Delta
        Delta_dat = numpy.loadtxt(self.dir+'/Delta.imp').T
        om = Delta_dat[0]
        Delta = Delta_dat[1::2] + Delta_dat[2::2]*1j

        # Strip duplicate/non-increasing frequency entries (e.g. sentinel 5000.0 repeated at end)
        _mask = numpy.concatenate(([True], numpy.diff(om) > 0))
        om = om[_mask]
        Delta = Delta[:, _mask]

        # Strip sentinel row (high-frequency tail marker at ω >> physical mesh).
        # merge_nambu_sigout appends a sentinel at ω=5000 for ssplit/sgather.
        # If kept, the interpolation at line 578 expands to ~40000 Matsubara points
        # instead of the physical ~500, causing bath fitting instability.
        if len(om) > 1:
            dw = om[1] - om[0]  # Matsubara spacing = 2π/β
            gap = om[-1] - om[-2]
            if gap > 10 * dw:  # sentinel has a huge gap from previous point
                print(f'Stripping sentinel frequency: om[-1]={om[-1]:.1f} (gap={gap:.1f} >> dw={dw:.4f})')
                om = om[:-1]
                Delta = Delta[:, :-1]

        ### If option to force zero self-energy for some orbitals
        # we have to be careful with Delta/Eimp and eventually replace
        # by the ones defined in some other impurity problem
        # as for instance in case we mix cluster and single-site impurities
        if 'force_szero' in params and params['force_szero']:
            # loop over indices registered in params['force_szero']
            # and replace Delta and Eimp by the one from the given other
            # impurity problem
            # Edc is not replaced however!
            for b in params['force_szero'][0]:
                # check provided target, if nothing then do nothing
                if isinstance(b,list) and len(b)>1:
                    imp, Dind2 = b[1].strip().split(':')
                    # read Delta from the given impurity problem
                    Ddat_tmp = numpy.loadtxt(self.dir+'../imp.'+imp
                                        +'/Delta.imp').transpose()
                    Dtmp = Ddat_tmp[1::2] + Ddat_tmp[2::2]*1j
                    # and now replace
                    Delta[b[0]] = Dtmp[int(Dind2)]
                    # same for Eimp now
                    Edat_tmp = numpy.loadtxt(self.dir+'../imp.'+imp
                                        +'/Eimp.inp',ndmin=2)[0]
                    Eimp[b[0]] = Edat_tmp[int(Dind2)]

        # subtracting Edc, because the ksum used (Hk+s_oo-Edc)
        # In Nambu mode, Edc in Eimp.inp ALREADY has τ₃ structure inherited from
        # sig.inp (Edc_hh = -Edc_pp). No additional sign flip needed here.
        # Previously line 452 applied an extra τ₃ flip, causing a double-flip
        # that made T_hh ≠ -T_pp (T_pp + T_hh = -2*Edc instead of 0).
        self.Eimps = Eimp - Edc

        # In case of fixn0 double counting, impurity levels shoud be the same
        if DCs=='fixn0':
            self.Eimps = numpy.ones(len(self.Eimps))*numpy.sum(self.Eimps)/len(self.Eimps)
            print('Eimps(fixn0)=', self.Eimps, file=self.fh_info)

        # iNORG Ed parameter: Edc from Eimp.inp already has τ₃ structure
        norg_eimps = Eimp - (Edc * Olap)
        formatted_eimps = "[{}]".format(" ".join(["{:.8e}".format(v) for v in norg_eimps.flat]))
        params['Ed'] = [formatted_eimps, "# Impurity levels"]


        ##### Prepare Uijkl interaction matrix  ###################################################
        is_cluster = bool('cluster' in params and params['cluster'][0])
        interaction_spec = resolve_interaction_spec(
            params, solver='iNORG', cluster=is_cluster)
        if is_cluster:
            print("Impurity "+self.dir+" is a cluster impurity, we reconstruct"
                  +" the interaction matrix accordingly", file=self.fh_info)
            _validate_inorg_interaction_spec(
                interaction_spec, self.CoulombF)
            U_tensor = ClusterCoulombMatrix(
                params, self.dir, self.fh_info, cix,
                interaction_spec=interaction_spec)
        else:
            print("Impurity "+self.dir+" is a single-site impurity, setting up"
                  +" interaction matrix normally", file=self.fh_info)
            U_tensor = CoulombMatrix(
                params, self.dir+'/Trans.dat', self.l, self.fh_info)
        self.fh_info.flush()

        #### Prepare executable  ##################################################################
        self.q_exe = self.env.ROOT+'/'+params['exe'][0]

        #### Preparing parameters file ############################################################
        # Parameters to skip (handled by solver or inorg_run.py)
        skip_params = {'Ed', 'U', 'J', 'fit_nbaths', 'norg_sets'}
        if self.nambu:
            skip_params.add('ful_pcl_sch')  # forced to 0 explicitly below
        # Parameters whose list values should be written without [] and commas
        list_params = {'fit_points', 'restrain', 'distribute', 'nb'}
        with open(self.dir+self.PARAMS, 'w') as fp:
            # Write if_details as 0/1 integer for C++ solver
            if 'if_details' in params:
                val = params['if_details']
                if isinstance(val, list):
                    val = val[0]
                int_val = 1 if val else 0
                print('{:10s}  {:10s}'.format('if_details', str(int_val)), file=fp)
            for p in ['if_spin_incl', 'orbital', 'if_mat_type', 'if_norg_degenerate', 'CoulombF', 'beta', 'basis', 'fit_range', 'fit_points', 'NOOC', 'restrain', 'distribute', 'nb', 'Minsulator', 'ful_pcl_sch', 'weight_nooc', 'weight_freze']:
                if p in params and p not in skip_params:
                    val = params[p][0]
                    if p in list_params:
                        val_str = str(val).replace('[', '').replace(']', '').replace(',', '')
                    else:
                        val_str = str(val)
                    print('{:10s}  {:10s}  {:s}'.format(p, val_str, params[p][1]), file=fp)

            # ---------------------------------------------------------------- #
            # Write sym parameter.  The bridge identifies normal mode from a   #
            # ``normal-*`` token, so accepting a superconducting token in a    #
            # normal calculation (or vice versa) changes interaction signs.    #
            # Fail instead of silently changing the requested model.           #
            # ---------------------------------------------------------------- #
            sym_value = _resolve_inorg_sym(self.nambu, params)

            print('{:10s}  {:10s}'.format('sym', sym_value), file=fp)

            # ---------------------------------------------------------------- #
            # Nambu mode: force ful_pcl_sch=0                                   #
            # ---------------------------------------------------------------- #
            if self.nambu:
                print('{:10s}  {:10s}'.format('ful_pcl_sch', '0'), file=fp)
                print('Nambu mode: forcing ful_pcl_sch=0 in PARAMS.norg',
                      file=self.fh_info)

            # ---- DISABLED (no longer needed): extract impurity-levels matrix from lda.outputdmf1 ----
            # outputdmf1_file = self.dir+'../lda.outputdmf1'  # Relative path from imp.x/ directory
            # if os.path.exists(outputdmf1_file):
            #     matrices = extract_impurity_matrix_from_outputdmf1(outputdmf1_file)
            #     if matrices:
            #         for icix, matrix in matrices.items():
            #             dim = matrix.shape[0]
            #             print("impurity levels =", file=fp)
            #             for i in range(dim):
            #                 print("  ", end="", file=fp)
            #                 for j in range(dim):
            #                     real_part = matrix[i, j].real
            #                     imag_part = matrix[i, j].imag
            #                     print(f" {real_part:14.8f}     {imag_part:14.8f}", end="", file=fp)
            #                 print("", file=fp)
            #             print(f"# {dim}x{dim} complex impurity levels matrix for icix={icix}", file=fp)


        # ---------------------------------------------------------------- #
        # nb provenance, and a tripwire against a silent drop.               #
        # 'nb' is the per-K bath-size vector.  A copy of this file predating  #
        # the key drops it without a word and the solver falls back to        #
        # AUTOMATIC bath splitting, which converges cleanly to a DIFFERENT    #
        # physical setup -- an unequal X/Y split feeds spurious X/Y asymmetry #
        # straight into Phi_d = (Phi_Y - Phi_X)/2.  Read PARAMS.norg back and #
        # record which source the solver will actually use, so the run log    #
        # settles the question with no guesswork.                             #
        # ---------------------------------------------------------------- #
        nb_written = None
        with open(self.dir+self.PARAMS, 'r') as fp:
            for line in fp:
                tok = line.split('#')[0].split()
                if tok and tok[0] == 'nb':
                    nb_written = ' '.join(tok[1:])
        if 'nb' in params:
            if nb_written is None:
                raise ValueError(
                    "iNORG: requested nb=" + str(params['nb'][0]) + " was DROPPED -- no nb line in "
                    + self.dir+self.PARAMS + ". This copy of inorg.py cannot pass nb through, so the "
                    "solver would silently fall back to automatic bath splitting (a different "
                    "physical setup, and an unequal X/Y split contaminates Phi_d). Refusing to run.")
            print('iNORG: nb=[' + nb_written + '] source=' + self.PARAMS, file=self.fh_info)
        else:
            print('iNORG: nb source=auto (not requested in impurity params; solver splits baths automatically)',
                  file=self.fh_info)

        Sigind = numpy.array(self.Sigind_shifted)
        Diagonal={}
        for i in range(len(Sigind)):
            for j in range(len(Sigind)):
                if Sigind[i,j]>0: Diagonal[Sigind[i,j]-1]= (i==j)

        print('Diagonal=', Diagonal, file=self.fh_info)


        #### Preparing Delta.inp file  ############################################################
        fp = open(self.dir+self.sparams['Delta'], 'w')             #Open Delta.inp
        # Checking that Delta is causal and not regular
        for ib in range(numpy.shape(Delta)[0]):
            if Diagonal[ib]:
                for im in range(numpy.shape(Delta)[1]):
                    if abs(Delta[ib,im])>maxAc*numpy.pi : Delta[ib,im] = -1e-10j
                    if Delta[ib,im].imag>0 : Delta[ib,im] = Delta[ib,im].real-1e-10j


        # creating big mesh to interpolate Delta
        T=1./params['beta'][0]
        m_max = int(om[-1]/(numpy.pi*T)+1e-6)         # n of the last mastubara point
        lom = numpy.arange(1,m_max+1,2)*numpy.pi*T          # all Matsubata points included

        # give some information on new mesh
        print('shape(Delta)=', numpy.shape(Delta), file=self.fh_info)
        print('shape(om)=', numpy.shape(om), file=self.fh_info)
        print('shape(lom)=', numpy.shape(lom), file=self.fh_info)
        print('ss=', numpy.shape(Delta[0].real), file=self.fh_info)
        print('om[0]=', om[0], 'lom[0]=', lom[0], 'om[-1]=', om[-1], 'lom[-1]=', lom[-1], file=self.fh_info)

        # interpolate on big mesh
        lDelta=[]
        for ib in range(len(Delta)):
            fDr = interpolate.UnivariateSpline(om, Delta[ib].real, s=0)
            fDi = interpolate.UnivariateSpline(om, Delta[ib].imag, s=0)
            lDelta.append( fDr(lom) + fDi(lom)*1j )

        for im,ome in enumerate(lom):
            print(ome, end=' ', file=fp)
            for b in range(len(lDelta)):
                print("{:19.12f}  {:19.12f} ".format(lDelta[b][im].real, lDelta[b][im].imag), end=' ', file=fp)
            print(file=fp)
        fp.close()
        # Now copy the Delta.inp to Delta.inp.globaliter.dmftiter
        shutil.copy2(self.dir+self.sparams['Delta'], self.dir+self.sparams['Delta']+'.'+extn)

        #### Prepare the input file with Eimp, Delta and Uijkl formatted for the ED solver
        prepare_EDinput(self.dir+'input.tsv', self.Sigind, self.Eimps,
                        lDelta, lom, U_tensor, self.CoulombF,
                        self.spin_dep, self.fh_info, nambu=self.nambu)


        # Write mpi_prefix.dat so inorg_run.py can use it for the C++ binary
        with open(self.dir + 'mpi_prefix.dat', 'w') as mpi_f:
            mpi_f.write(self.mpi_prefix)

        # Solver command — inorg_run.py imports the (now unified) inorg_edmft_io
        # which auto-detects Nambu vs normal via PARAMS.norg `sym` field.
        # The bootstrap-script injection used previously (when inorg_edmft_io_nambu.py
        # was separate) is no longer needed: a single inorg_edmft_io.py handles both.
        cmd = 'cd '+self.dir+'; python3 '+self.q_exe+' > nohup_imp.out 2>&1 '

        print(cmd, file=self.fh_info)
        print('Running ---- inorg -----', file=self.fh_info)
        if subprocess.call(cmd,shell=True,stdout=self.fh_info,stderr=self.fh_info) !=0:
            print("ERROR IN inorg SOLVER: check nohup_imp file. Might be an error due to lanczos not converging, modyfing these global parameters might help: \n"
            + "max_iter_lanczos, Ground_state_method\n"
            + "The bath model can also be at fault...")
            sys.exit(0)

        # ------------------------------------------------------------------ #
        # Nambu: verify that Phi.out was produced by the solver.              #
        # iNORG writes Phi.out only when sym='psc' is active.                 #
        # ------------------------------------------------------------------ #
        if self.nambu:
            phi_out_path = self.dir + 'Phi.out'
            if not os.path.exists(phi_out_path):
                raise FileNotFoundError(
                    f"Nambu mode: Phi.out not found at {phi_out_path!r} after "
                    "iNORG solver run.  Check that the C++ iNORG binary was "
                    "compiled with Nambu support and that PARAMS.norg has a "
                    "Nambu-compatible sym (psc/dsc/K-dsc).  The anomalous self-energy output file (Phi.out) "
                    "is required for phi.out to be written in HighFrequency()."
                )
            print(f'Nambu mode: verified Phi.out exists at {phi_out_path}', file=self.fh_info)


    def HighFrequency(self, params, DCs, nl_imp, extn, wbroad=0.0, kbroad=0.0):
        """ Reads impurity output and prepares for further execution
          Output:
               Sig[b][iom]  -- DMFT dynamic self-energy which vanishes at infinity
               sinf[b]      -- DMFT self-energy at infinity
               Edc[b]       -- double counting using DMFT occupancies

          Nambu extension:
               When self.nambu is True, additionally reads Phi.out (anomalous
               self-energy) and writes phi.out with headers
                   # s_oo= [0.0, ...]
                   # Edc=  [0.0, ...]
               The anomalous self-energy Φ decays to zero at high frequency
               (s_oo=0) and has no double-counting correction (Edc=0).
               The phi.out values are written *without* any DC subtraction.
        """
        # Following edqcm.py, read occupations from impurity solver output
        mom = []
        nf = 0.0
        nmat_path = self.dir + 'nmat.txt'
        if os.path.exists(nmat_path):
            print('Reading impurity occupation from', nmat_path, file=self.fh_info)
            (nf, mom, s_p, s_h) = _parse_nmat_txt(nmat_path, nambu=self.nambu,
                                                    fh_info=self.fh_info)

        if self.nambu:
            print('From nmat.txt (Nambu K-basis): nf=', nf, 'mom=', mom,
                  's_p=', s_p, 's_h=', s_h, file=self.fh_info)
        else:
            print('From nmat.txt (K-basis): nf=', nf, 'mom=', mom, file=self.fh_info)

        # K-native discrimination via PP-BLOCK Sigind diagonality (matches
        # inorg_run.py; works for normal AND Nambu): a DIAGONAL pp-block means the
        # solver worked in the K basis, so nmat.txt occupations are per-K and must
        # NOT be rotated back to site.  Nambu sigind is 2N x 2N with the anomalous
        # Phi off the main diagonal; test the pp block sig[:N,:N].
        _sig_full = numpy.array(self.Sigind)
        _half = _sig_full.shape[0] // 2
        _sig_pp = _sig_full[:_half, :_half] if _half >= 1 else _sig_full
        is_k_native = (_sig_pp.size > 0 and numpy.count_nonzero(
            _sig_pp - numpy.diag(numpy.diag(_sig_pp))) == 0)
        kspace_numsites = 0

        # Transform K-basis occupations to site-basis for eDMFT double-counting
        # C++ solver works in K-basis, but ComputeExactDC expects site-basis occupations
        # n_site[r] = sum_K |U[K,r]|^2 * n_K[K]
        # For 2x2 plaquette: U is real, |U[K,r]|^2 = 1/4 for all K,r
        if mom:
            nband = len(mom)
            Sigind_block = numpy.array(self.Sigind)[:nband, :nband]
            max_sig = int(numpy.max(Sigind_block)) if nband > 0 else 0
            if is_k_native and max_sig > 0:
                # K-native: mom from nmat.txt is already K-basis per-orbital
                # [n_Γ, n_X, n_Y, n_M].  Compress by SUMMING per diagonal-Sigind
                # component (NO K->site rotation):  -> [n_Γ, n_X+n_Y, n_M].
                # This applies to Nambu as well (the discriminator above already
                # tests the pp block): Sigind_block is then sig[:nband,:nband], so
                # mom keeps the same K-resolved meaning as in the normal run
                # instead of being rotated into (by symmetry degenerate) sites.
                kspace_numsites = nband
                mom_Kbasis = list(mom)
                mom = [0.0] * max_sig
                for k in range(nband):
                    s = int(Sigind_block[k, k])
                    if s > 0:
                        mom[s-1] += mom_Kbasis[k]
                print('K-native Sigind-compressed mom (sum per K-component):',
                      mom_Kbasis, '->', mom, file=self.fh_info)
            else:
                U = plaquette_unitary()  # 4x4 real unitary
                if nband == U.shape[0]:
                    nK = numpy.diag(mom)  # K-basis occupation matrix (diagonal)
                    n_site_matrix = U.conj().T @ nK @ U  # K->site transform
                    print('K->site n_site_matrix diagonal=',
                          [n_site_matrix[i, i].real for i in range(nband)], file=self.fh_info)

                    # Compress site-basis density matrix to sigind-compressed format
                    # Sum (not average) n_site_matrix elements grouped by sigind component
                    # Convention: mom[s-1] = total occupancy of sigind component s
                    # DC code divides by n_per_elem_ (or deg) to get per-orbital occupancy
                    if max_sig > 0:
                        mom_Kbasis = list(mom)
                        mom = [0.0] * max_sig
                        for s in range(1, max_sig + 1):
                            positions = [(i, j) for i in range(nband) for j in range(nband)
                                         if int(Sigind_block[i, j]) == s]
                            if positions:
                                mom[s-1] = sum(n_site_matrix[i, j].real for i, j in positions)
                        print('Sigind-compressed mom:', mom_Kbasis, '->', mom, file=self.fh_info)
                else:
                    print('WARNING: K->site transform skipped (nband=%d != U.shape=%d)' %
                          (nband, U.shape[0]), file=self.fh_info)

        # Prepending header to Gf.out and Sig.out
        if nf is not None and mom is not None:
            mom_str = f"[{','.join(map(str, mom))}]"
            header = f"# nf={nf} mom={mom_str}\n"
            files_to_header = [
                os.path.join(self.dir, self.sparams['Gf']),
                os.path.join(self.dir, self.sparams['Sig'])
            ]
            for file_path in files_to_header:
                if os.path.exists(file_path):
                    try:
                        with open(file_path, 'r+') as f:
                            content = f.read()
                            f.seek(0, 0)
                            f.write(header + content)
                    except IOError as e:
                        print(f"Could not write header to {file_path}: {e}", file=self.fh_info)


        if kbroad>0.0 or wbroad>0.0:
            bexe = self.env.ROOT + '/broad'
            shutil.move(self.dir+self.sparams['Sig'], self.dir+self.sparams['Sig']+'b')
            broad_cmd = 'cd '+self.dir+'; '+bexe+' -w '+str(wbroad)+' -k '+str(kbroad)+' '+self.sparams['Sig']+'b  > '+self.sparams['Sig']
            print(broad_cmd, file=self.fh_info)
            ret = subprocess.call(broad_cmd,shell=True,stdout=self.fh_info,stderr=self.fh_info)

        OffDiagonalExist = False
        Sigind = numpy.array(self.Sigind_shifted)
        self.Diagonal={}
        for i in range(len(Sigind)):
            for j in range(len(Sigind)):
                if Sigind[i,j]>0:
                    self.Diagonal[Sigind[i,j]-1]= (i==j)
                    if i!=j: OffDiagonalExist = True


        # Reading Self-energy of the impurity
        Sg_data = numpy.loadtxt(self.dir + self.sparams['Sig'], comments='#').T
        om = Sg_data[0]
        Sigo = Sg_data[1::2,:]+Sg_data[2::2,:]*1j


        # Reading Gf.out which was produced by inorg
        try:
            with open(self.dir + self.sparams['Gf'],'r') as fS:
                first_line = next(fS).strip()
            adat = first_line.split()
            header_data = {}
            for par in adat:
                if '=' in par:
                    key, val_str = par.split('=', 1)
                    key = key.strip()
                    try:
                        header_data[key] = eval(val_str)
                    except (NameError, SyntaxError):
                        header_data[key] = val_str
            
            if 'nf' in header_data:
                nf = header_data['nf']
            if 'mom' in header_data:
                mom = header_data['mom']
            TrSigmaG = header_data.get('TrSigmaG', 0.0)
            ntot = nf
            print('From impurity Gf.out: nf=', nf, 'TrSigmaG=', TrSigmaG, 'mom=', mom, file=self.fh_info)
        except (IOError, StopIteration):
            print("Could not read header from Gf.out, using values from nmat.txt", file=self.fh_info)
            TrSigmaG = 0.0
            ntot = nf

        for b in range(numpy.shape(Sigo)[0]):   #print b, 'This is diagonal bath, correcting it'
            for iom in range(numpy.shape(Sigo)[1]):
                if self.Diagonal[b] and Sigo[b,iom].imag>0:
                    Sigo[b,iom] = Sigo[b,iom].real
                    print('Disregarding non-causal self-energy at Sig['+str(b)+','+str(iom)+']=',Sigo[b,iom], file=self.fh_info)

        sinf = numpy.copy(Sigo[:,-1]).real

        print('sinf=', sinf, file=self.fh_info)


        # Reading old Edc and impurity levels
        (Eimp_old,Olap_old,Edc_old) = numpy.loadtxt(self.dir+'/Eimp.inp',ndmin=2)

        print('Diagonal=', self.Diagonal, file=self.fh_info)
        print('len(Edc_old)=', len(Edc_old), 'Eimp_old=', Eimp_old, file=self.fh_info)


        # Here we compute self-energy(infinity) and
        # double-counting using impurity occupancy
        print('DCs=', DCs, file=self.fh_info)
        Ud = params['U'][0]
        if type(self.J) is list:
            Jd = sum(self.J)/len(self.J)
        else:
            Jd = self.J
        # If we want to have U different than the U used in double counting
        # (useful for cluster-DMFT when U can not be simply added in ctqmc)
        if 'Ud' in params: Ud = params['Ud'][0]
        if 'Jd' in params: Jd  = params['Jd'][0]
        DC_ones = numpy.array([int(self.Diagonal[i]) for i in range(len(Edc_old))])

        # Handle missing UlamJ gracefully — fall back to FLL if exact DC not available
        _DCs_effective = DCs
        if DCs[:5]=='exact' and getattr(self, 'UlamJ', None) is None:
            print(f'WARNING: DCs={DCs!r} requires UlamJ in params but it was not set. '
                  f'Falling back to FLL double-counting.',
                  file=self.fh_info)
            _DCs_effective = 'FLL'

        _nambu_freeze_dc = bool(params['nambu_freeze_dc'][0]) if 'nambu_freeze_dc' in params else True

        if getattr(self, 'nambu', False) and _nambu_freeze_dc:
            nf_ = mom
            Edc0 = numpy.array(Edc_old, dtype=float)
            Phidc0 = 0.0
            print('Nambu: DC frozen at the Eimp.inp value (normal-state reference); Vdc=',
                  Edc0, file=self.fh_info)
        elif _DCs_effective[:5]=='exact':
            # We allow DCs = ['exact', 'exact_l', 'exactd', 'exactdl']
            # exact and exactd use impurity occupation. exact uses yukawa&dielectric, while exactd uses dielectric screening only
            # exact_l and exactdl use lattice occupations. exact_l uses yukawa&dielectric, while exactdl uses dielectric screening only
            lmbda = self.UlamJ[1]
            epsilon = self.UlamJ[2]
            nf_ = mom                         # 'exact': Density determined from impurity occupation
            print('Here we start n_imp=', nf_, 'isscalar(nlattice)=', numpy.isscalar(nl_imp), file=self.fh_info)
            if len(DCs)>6 and DCs[6]=='l' and not numpy.isscalar(nl_imp):
                nf_ = nl_imp
            print('Here we have nl_imp=', nf_, file=self.fh_info)
            print('lmbda=', lmbda, 'epsilon=', epsilon, file=self.fh_info)

            print('Here we have nl_imp=', mom)

            if is_k_native:
                # K-native exact DC (mirror bctqmc.py): exact DC is a LOCAL atomic
                # quantity.  Reduce the K-decomposed occupations to the per-site
                # occupation n_local = Tr(n)/numsites and call the SINGLE-SITE
                # ComputeExactDC with the invariant local x2y2 d-orbital T2C.  Vdc
                # is K-independent (proven K-invariant) -> broadcast to all comps.
                numsites = kspace_numsites if kspace_numsites > 0 else len(nf_)
                n_local = sum(nf_) / float(numsites)   # uniform plaquette |U[K,r]|^2=1/numsites
                _s = numpy.sqrt(0.5)
                _T2C = numpy.array([[   _s,     0, 0,     0,    _s],
                                    [    0,     0, 1,     0,     0],
                                    [    0,    _s, 0,   -_s,     0],
                                    [    0, 1j*_s, 0, 1j*_s,     0],
                                    [-1j*_s,     0, 0,     0, 1j*_s]], dtype=complex)
                _pp_trans = self.dir + 'Trans_loc_pp.dat'
                with open(_pp_trans, 'w') as _ft:
                    print('5 5', file=_ft)
                    print('#---- Sigind (single-site pp-block for exact DC; x2y2 correlated)', file=_ft)
                    for _i in range(5):
                        print(' '.join(('1' if (_i == 0 and _j == 0) else '0') for _j in range(5)), file=_ft)
                    print('#---- CF', file=_ft)
                    for _i in range(5):
                        for _j in range(5):
                            print('%12.8f %12.8f' % (_T2C[_i, _j].real, _T2C[_i, _j].imag), end='  ', file=_ft)
                        print(file=_ft)
                (Edc_loc, Phidc0_site) = dmft_DC.ComputeExactDC(
                    [n_local], lmbda, epsilon, self.icase, self.fh_info,
                    projector='projectorw.dat', trans=_pp_trans)
                Vdc = Edc_loc[0]
                if not (1.0 < Vdc < 10.0):
                    print('WARNING exact-DC: Vdc=%g outside sane [1,10] eV (n_local=%g)'
                          % (Vdc, n_local), file=self.fh_info)
                if self.nambu:
                    # Broadcast the scalar normal-Hartree Vdc into the Nambu tau3
                    # structure: pp gets +Vdc, hh gets -Vdc, the anomalous Phi/PhiD
                    # blocks get 0 (no DC for the anomalous channel).  This mirrors
                    # run_dmft.py:455 Edc_nambu = Edc_pp + [0]*K + [0]*K + [-Edc_pp]
                    # and the sig.out tau3 header.  (DC_ones*Vdc would put +Vdc on the
                    # hh block too -- wrong sign.)
                    # Length of ONE tau3 block = number of independent pp-block
                    # components = len(nf_).  That is NOT the site count whenever
                    # cluster momenta are degenerate; kspace_numsites is the site
                    # count and belongs to numsites above, not here.
                    K_pp = len(nf_)
                    Edc0 = numpy.array([Vdc]*K_pp + [0.0]*K_pp + [0.0]*K_pp + [-Vdc]*K_pp)
                    assert len(Edc0) == len(Edc_old), \
                        'Nambu tau3 Edc0 len %d != Edc_old len %d' % (len(Edc0), len(Edc_old))
                else:
                    Edc0 = DC_ones * Vdc
                Phidc0 = numsites * Phidc0_site        # cluster total (equivalent sites)
                print('K-native Exact DC: Vdc=', Edc0, 'PhiDC=', Phidc0,
                      '(local n=%g, numsites=%d)' % (n_local, numsites), file=self.fh_info)

            # Have to check here if it is cluster impurity or not
            # If not, we copmute DC in the standard way
            # If cluster, we treat each correlated atom separately
            # and we ignore the inter-atom off-diagonal terms
            elif 'cluster' in params and params['cluster'][0]:
                print('Dealing with cluster impurity: exact DC is calculated'
                      +' for each atom individually, and information'
                      + ' is recombined in Vdc and Phi_DC. WARNING: inter-atomic elements are ignored!',
                      file=self.fh_info)

                #### initialize
                Phidc0 = 0
                # first find all independent elements in total Sigind: this will be the shape of final output
                elem_ = numpy.unique(self.Sigind)[1:]     # remove non-correlated elements
                n_per_elem_ = []
                for ind in range(len(elem_)):
                    n_per_elem_.append((self.Sigind).size - numpy.count_nonzero(Sigind-ind-1))
                Edc0 = numpy.zeros(len(elem_))
                divide_ = numpy.zeros(len(elem_),dtype=int)

                #### loop over all atoms, copmute exact DC for each, add everything to PhiDC and Edc0
                # In Nambu mode, nf_ only has pp-block sigind components (1..K_pp).
                # Nambu-extended sigind values (K_pp+1..4K_pp) must be skipped.
                # Also, ComputeExactDC reads Trans_loc internally and crashes on Nambu
                # sigind (10×10 with indices up to 12, but nocc has only K_pp elements).
                # Fix: write a temporary pp-block-only Trans file for ComputeExactDC.
                max_nf = len(nf_)  # number of pp-block sigind components
                iat = 0
                while (os.path.isfile(self.dir+'Trans_loc'+str(iat)+'.dat')):   # we search all the local Trans files
                    # read the local sigind
                    Sigind_loc, CF_loc = ReadTrans(self.dir+'Trans_loc'+str(iat)+'.dat', self.fh_info)

                    # adapt the nf_ table to the local atom
                    # Skip Nambu-extended sigind values (icor > max_nf)
                    nf_loc_ = []
                    for icor in numpy.unique(Sigind_loc)[1:]:
                        if icor <= max_nf:
                            nf_loc_.append(nf_[icor-1]/n_per_elem_[icor-1])              # should divide by number of atoms having this icor

                    # In Nambu mode, write a pp-block-only Trans file for ComputeExactDC.
                    # ComputeExactDC reads the Trans file internally and crashes on the
                    # Nambu-doubled Sigind (10×10 with sigind up to 12, but nocc has K_pp elements).
                    trans_file = self.dir+'Trans_loc'+str(iat)+'.dat'
                    if self.nambu:
                        ndim = Sigind_loc.shape[0]
                        ndim_pp = ndim // 2
                        Sigind_pp = Sigind_loc[:ndim_pp, :ndim_pp]
                        CF_pp = CF_loc[:ndim_pp, :ndim_pp] if CF_loc is not None else numpy.eye(ndim_pp)
                        trans_file = self.dir+'Trans_loc'+str(iat)+'_pp.dat'
                        with open(trans_file, 'w') as ft:
                            # ReadTrans format: line0=dims, line1=comment, lines2..=Sigind, comment, CF
                            print(f'{ndim_pp} {ndim_pp}', file=ft)
                            print('#---- Sigind (pp-block only for Nambu DC)', file=ft)
                            for i in range(ndim_pp):
                                print(' '.join(str(int(Sigind_pp[i,j])) for j in range(ndim_pp)), file=ft)
                            print('#---- CF', file=ft)
                            for i in range(ndim_pp):
                                for j in range(ndim_pp):
                                    print(f'{CF_pp[i,j].real:12.8f} {CF_pp[i,j].imag:12.8f}', end='  ', file=ft)
                                print(file=ft)
                        print(f'Nambu: wrote pp-block Trans_loc{iat} ({ndim_pp}x{ndim_pp}) for ComputeExactDC', file=self.fh_info)

                    # get the local DC
                    (Edc0_loc,Phidc0_loc)=dmft_DC.ComputeExactDC(nf_loc_, lmbda, epsilon, self.icase, self.fh_info, projector='projectorw.dat',trans=trans_file)

                    # scatter the data to the right elements
                    Phidc0 += Phidc0_loc
                    ii = 0
                    for icor in numpy.unique(Sigind_loc)[1:]:
                        if icor > max_nf:
                            continue  # skip Nambu-extended sigind (Phi/PhiD/hh blocks)
                        # only diagonal elements of Sigind_loc are taken into account in ComputeExactDc
                        if icor in numpy.diagonal(Sigind_loc):
                            Edc0[icor-1]+=Edc0_loc[ii]
                            divide_[icor-1]+= 1
                            ii+=1
                        # intra-atomic off-diagonal terms also have zero double-counting
                        else:
                            Edc0[icor-1]+= 0.0
                            divide_[icor-1]+= 1
                    iat += 1

                # Vdc should be quantity per orbital -> need to divide
                for iel in range(len(elem_)):
                    # if divide is zero it means it is off-diagonal -> we let it to zero
                    if divide_[iel] !=0:
                        Edc0[iel] = Edc0[iel]/divide_[iel]

            else:
                # Even if not cluster, stil have to check if non-zero off diagonal Sigind
                elem_ = numpy.unique(self.Sigind)[1:]     # remove non-correlated elements
                n_per_elem_ = []
                for ind in range(len(elem_)):
                    n_per_elem_.append((self.Sigind).size - numpy.count_nonzero(Sigind-ind-1))
                Edc0 = numpy.zeros(len(elem_))

                # get double-counting
                (Edc0_loc,Phidc0)=dmft_DC.ComputeExactDC(nf_, lmbda, epsilon, self.icase, self.fh_info, projector='projectorw.dat',trans=self.dir+'Trans.dat')

                # scatter the Edc0 data tot ake care of possible intra-atomic off-diagonals
                ii = 0
                for icor in numpy.unique(self.Sigind)[1:]:
                    # only diagonal elements of Sigind_loc are taken into account in ComputeExactDc
                    if icor in numpy.diagonal(self.Sigind):
                        Edc0[icor-1]=Edc0_loc[ii]
                        ii+=1
                    # intra-atomic off-diagonal terms also have zero double-counting
                    else:
                        Edc0[icor-1]= 0.0

            print('The Exact Double-counting is: Vdc=', Edc0, 'PhiDC=', Phidc0, file=self.fh_info)

        elif _DCs_effective=='FLL' or _DCs_effective=='default' or _DCs_effective=='FLL_l':
            nf_ = nf
            if DCs=='FLL_l' and nf_:
                nf_ = sum(nl_imp)
            # This is the fully localized limit double-counting by Anisimov
            Edc0 = DC_ones*( Ud*(nf_-0.5) - Jd*(nf_/2.-0.5) )
            Phidc0 = Ud*0.5*nf_*(nf_-1.) - Jd*0.25*nf_*(nf_-2.)
        elif _DCs_effective=='AMF':
            nf_ = nf
            if DCs=='AMF_l' and nf_:
                nf_ = sum(nl_imp)
            Ueff = Ud*(1-1./(2.*(2.*self.l+1.))) - Jd*self.l/(2*self.l+1.)
            Edc0 = DC_ones*Ueff*nf_
            Phidc0 = Ueff*nf_**2/2.
        elif _DCs_effective=='fixn' or _DCs_effective=='nominal':         # If the scheme if fixn, we keep double-counting potential fixed
            nf0 = params['nf0'][0]
            Vdc = Ud*(nf0-0.5) - Jd*(nf0/2.-0.5)
            Edc0 = DC_ones * Vdc
            Phidc0 = ntot*Vdc
        else:
            print('ERROR: Unkown double-counting', DCs)
            sys.exit(1)

        if 'dDC' in params:
            Edc0 += numpy.array(params['dDC'][0])
            Phidc0 += numpy.sum([params['dDC'][0][i]*nf_[i] for i in range(len(nf_))])
            print('The Double-counting corrected to: Vdc=', Edc0, 'PhiDC=', Phidc0, file=self.fh_info)


        Edc = Edc_old # we want to keep Edc<-1000 for orbitals which are projected out!
        for i in range(len(Edc)):
            if (Edc_old[i]>-1000.): # this is true for all orbitals which are kept!
                Edc[i] = Edc0[i]
            else:                   # these orbitals are projected out -> we want to keep sinf=Edc for orbitals projected out!
                sinf[i] = Edc[i]

        # Copying data at each iteration to follow the evolution
        shutil.copy2(self.dir+self.sparams['Sig'], self.dir+self.sparams['Sig']+'.'+extn)
        shutil.copy2(self.dir+self.sparams['Gf'], self.dir+self.sparams['Gf']+'.'+extn)
        shutil.copy2(self.dir+'PARAMS.norg', self.dir+'PARAMS.norg'+'.'+extn)
        shutil.copy2(self.dir+'nohup_imp.out', self.dir+'norg.log'+'.'+extn)
        if os.path.exists(self.dir+'hb_fit.txt'): shutil.copy2(self.dir+'hb_fit.txt', self.dir+'hb_fit'+'.'+extn)
        if os.path.exists(self.dir+'ose_hop'): shutil.copy2(self.dir+'ose_hop', self.dir+'ose_hop'+'.'+extn)
        if os.path.exists(self.dir+'nmat.txt'): shutil.copy2(self.dir+'nmat.txt', self.dir+'nmat'+'.'+extn)
        if os.path.exists(self.dir+'h0.txt'): shutil.copy2(self.dir+'h0.txt', self.dir+'h0'+'.'+extn)
        if os.path.exists(self.dir+'solver_input.tsv'): shutil.copy2(self.dir+'solver_input.tsv', self.dir+'solver_input.tsv'+'.'+extn)
        if os.path.exists(self.dir+'F.out'): shutil.copy2(self.dir+'F.out', self.dir+'F.out'+'.'+extn)

        self.fh_info.flush()
        
        Ry2eV = 13.60569193
        Epotential = TrSigmaG
        
        with open(self.dir+'/Eorb.dat', 'w') as fE:
            print('#  Phidc=', 0.0, file=fE)
            print('#  Tr(Sigma*G)/2=', 0.0, file=fE)
            print('#  Tr(Sigma*G)/2-Phidc=', 0.0, file=fE)
            print('#  Phi_DMFT=', 0.0, file=fE)
            print('#  Phi_DMFT-Phi_DC=', 0.0, file=fE)
            print('#  Tr(log(-Gimp))=', 0.0, file=fE)
            print('#  Fimpurity+TS=', 0.0, file=fE)
            print('#  E_pot+E_kin-Tr(w*dD/dw)=', 0.0, file=fE)
            print(':IEORB ', 0.0, file=fE)
            print(':XEORB ', 0.0, file=fE)
            print(':ZEORB ', 0.0, file=fE)

        ### Option here to force certain orbitals
        # to have zero self-energy
        if 'force_szero' in params.keys() and params['force_szero']:
            # loop over indices registered in params['force_szero']
            # put the corresponding Sigo and sinf to EDc
            for b in params['force_szero'][0]:
                # quick check that b is list or not
                # whther there is replacement of Delta or not
                if isinstance(b,list):
                    ind = b[0]
                else:
                    ind = b
                sinf[ind] = Edc[ind]
                Sigo[ind,:] = sinf[ind]

        # ------------------------------------------------------------------ #
        # Write sig.out — normal self-energy (identical to upstream)          #
        # ------------------------------------------------------------------ #
        fE = open(self.dir+'/sig.out', 'w')
        # Truncate Edc to match data dimensions: Sigo has len(sinf) components
        # but Edc may have more (e.g., 12 Nambu elements from Eimp.inp vs 3 pp-block data cols).
        # merge_nambu_sigout reconstructs the full 4K Edc with τ₃ signs from the pp-block Edc.
        n_data_comp = len(sinf)
        print('# s_oo=', sinf.tolist(), file=fE)
        print('# Edc=', Edc[:n_data_comp].tolist(), file=fE)

        for iom in range(len(om)):
            print(("%20.15f " % om[iom]), end=' ', file=fE)
            for b in range(len(Sigo)):
                print(("%20.15f %20.15f  " % (Sigo[b,iom].real-Edc[b], Sigo[b,iom].imag)), end=' ', file=fE)
            print(file=fE)
        fE.close()

        # ------------------------------------------------------------------ #
        # Nambu: read Phi.out → write phi.out                                  #
        #                                                                       #
        # Phi.out is the anomalous self-energy produced by iNORG's sym='psc'   #
        # mode.  It has the same column layout as Sig.out:                     #
        #   col 0: omega_n                                                      #
        #   col 2b+1: Re(Phi_b)                                                 #
        #   col 2b+2: Im(Phi_b)                                                 #
        #                                                                       #
        # Key differences from sig.out:                                         #
        #   - s_oo = 0  (Φ decays to zero at high frequency)                   #
        #   - Edc  = 0  (no double-counting for anomalous component)            #
        #   - NO subtraction of Edc[b] when writing data                        #
        # ------------------------------------------------------------------ #
        if self.nambu:
            phi_out_solver = self.dir + 'Phi.out'  # output from iNORG C++ solver

            if not os.path.exists(phi_out_solver):
                # _exe() should have caught this already, but guard here too.
                print('WARNING: Nambu mode active but Phi.out not found; '
                      'phi.out will NOT be written.', file=self.fh_info)
            else:
                # Read anomalous self-energy — same format as Sig.out
                Phi_data = numpy.loadtxt(phi_out_solver, comments='#').T
                phi_om = Phi_data[0]
                Phio = Phi_data[1::2, :] + Phi_data[2::2, :] * 1j

                num_components = len(Phio)
                phi_s_oo = [0.0] * num_components
                phi_Edc  = [0.0] * num_components

                print(f'Nambu: writing phi.out ({num_components} components, '
                      f'{len(phi_om)} frequencies)', file=self.fh_info)

                fPhi = open(self.dir + '/phi.out', 'w')
                print('# s_oo=', phi_s_oo, file=fPhi)
                print('# Edc=', phi_Edc, file=fPhi)

                for iom in range(len(phi_om)):
                    print(("%20.15f " % phi_om[iom]), end=' ', file=fPhi)
                    for b in range(num_components):
                        # NO Edc subtraction — anomalous self-energy has Edc=0
                        print(("%20.15f %20.15f  " % (Phio[b, iom].real,
                                                        Phio[b, iom].imag)),
                              end=' ', file=fPhi)
                    print(file=fPhi)
                fPhi.close()

                # Archive Phi.out and F.out for iteration tracking
                shutil.copy2(phi_out_solver, self.dir + 'Phi.out.' + extn)
                f_out_solver = self.dir + 'F.out'
                if os.path.exists(f_out_solver):
                    shutil.copy2(f_out_solver, self.dir + 'F.out.' + extn)
                print(f'Nambu: phi.out written, Phi.out/F.out archived as .{extn}',
                      file=self.fh_info)

                # Merge sig.out + phi.out → unified sig.out (8K+1 cols)
                merge_nambu_sigout(
                    self.dir + '/sig.out',
                    self.dir + '/phi.out',
                    fh_info=self.fh_info
                )

        return ntot


def SlaterF(U, J, l):
    Fk = numpy.zeros((4,4), dtype=float)
    if l==0:
        # F0 for s-electrons
        Fk[0,0] = U
    elif l==1:
        # F2 for p-electrons
        Fk[0,1] = U
        if type(J) is list:
            Fk[1,1] = 5*J[0]
        else:
            Fk[1,1] = 5*J
    elif l==2:
        # F2 and F4 for d-electrons
        Fk[0,2] = U
        if type(J) is list:
            Fk[1,2] = 14./1.625 * J[0]
            Fk[2,2] = 14.*0.625/1.625 * J[1]
        else:
            Fk[1,2] = 14./1.625 * J
            Fk[2,2] = 14.*0.625/1.625 * J
    elif l==3:
        # F2, F4 and F6 for f-electrons
        Fk[0,3] = U
        if type(J) is list:
            Fk[1,3] = 6435./(286+195*0.668+250*0.494) * J[0]
            Fk[2,3] = 0.668*6435./539.76 * J[1]
            Fk[3,3] = 0.494*6435./539.76 * J[2]
        else:
            Fk[1,3] = 6435./(286+195*0.668+250*0.494) * J
            Fk[2,3] = 0.668*6435./539.76 * J
            Fk[3,3] = 0.494*6435./539.76 * J
    return Fk

def Spheric2Cubic(L):
    """Transformation matrix from spherical to cubic harmonics:
    From |l,m> (l=0,1,2,3) basis to
    basis of triplets, doublets and singlets real functions
    (see appendices 2.1 and 4.2 of H. Watanabe 'Operator Methods in Ligand Field Theory'
    Prentice Hall, 1966.)
    """
    s2 = numpy.sqrt(2)
    f=1/s2
    g=1j/s2
    if L == 0:
        return numpy.array([[1]])
    elif L == 1:
        T = [[ f, 0,-f],  # x
             [ g, 0, g],  # y
             [ 0, 1, 0]]  # z
        return numpy.array(T)
    elif L == 2:
        T = [[ 0, 0, 1, 0, 0],  # z^2
             [ f, 0, 0, 0, f],  # x^2-y^2
             [ 0, f, 0,-f, 0],  # xz
             [ 0, g, 0, g, 0],  # yz
             [ g, 0, 0, 0,-g]]  # xy
        return numpy.array(T)
    elif L == 3:
        s3 = numpy.sqrt(3)
        s5 = numpy.sqrt(5)
        s8 = numpy.sqrt(8)
        a=s5/4.
        b=s3/4.
        c=1j*a
        d=1j*b
        T=[[ 0, g, 0, 0, 0,-g, 0], # singlet
           [ a, 0,-b, 0, b, 0,-a], # x
           [-c, 0,-d, 0,-d, 0,-c], # y
           [ 0, 0, 0, 1, 0, 0, 0], # z
           [-b, 0,-a, 0, a, 0, b], # ksi
           [-d, 0, c, 0, c, 0,-d], # eta
           [ 0, f, 0, 0, 0, f, 0]] # zeta
        return numpy.array(T)

def CoulUsC2(l, T2C):
    # Precomputes gaunt coefficients for speed
    # shape(gck)=(4,7,7,4). It contains gck(l, m4, m1, k)
    gck = gaunt.cmp_all_gaunt()

    gck = numpy.array(gck, order='C')   # This is essential for pybind11 code to work with fortran gck
    T2C = numpy.array(T2C, order='C')
    nw = len(T2C)
    UC = numpy.zeros((l+1, nw, nw, nw, nw), dtype=complex)
    dpybind.FromSlaterToMatrixU(UC, gck, l, T2C)     # uses the dpybind routines for atom_d
    return UC

def ReadTrans(filename, fh_info):
    """Read the self-energy index file Sigind and the local transformation matrix CF from a file"""
    with open(filename, 'r') as fh:
        data = fh.readlines()

    (n1,n2) = list(map(int, data[0].split()[:2]))

    Sigind=[]
    for i in range(n1):
        Sigind.append( list(map(int, data[i+2].split()[:n2])) )
    Sigind = numpy.array(Sigind)

    print('len(data)', len(data), file=fh_info)
    print('n1=', n1, file=fh_info)
    if len(data) >= n1+n1+3:
        n2 = n1
        CF=[]
        for i in range(n2):
            cl = numpy.array(list(map(float, data[n1+3+i].split())))
            CF.append( cl[0::2]+cl[1::2]*1j )
        CF = numpy.array(CF)
    elif len(data)>=n1+n1/2+3:
        n2 = int(n1/2)
        CF=[]
        for i in range(n2):
            cl = numpy.array(list(map(float, data[n1+3+i].split())))
            CF.append( cl[0::2]+cl[1::2]*1j )
        CF = numpy.array(CF)
        CFN = numpy.zeros((2*n2,2*n2), dtype=complex)
        CFN[:n2,:n2] = CF
        CFN[n2:,n2:] = CF
        CF = CFN
    else:
        CF = numpy.identity(n1)
    print('shape(CF)=', numpy.shape(CF), 'shape(Sigind)=', numpy.shape(Sigind), file=fh_info)

    return (Sigind, CF)

def Check_T2C_Real(T2C, l, fh_info, small):
    """
     Here we added a routine which checks that cubic harmonics are real.
     Only in this case operators F^+ and F used in ctqmc will be real.
     Otherwise these matrix elements might be complex.
     The condition for cubic harmonics to be real is:

      Imag( T2C[m,i] + (-1)**m * T2C[-m,i] ) == 0
        and
      Real( T2C[m,i] - (-1)**m * T2C[-m,i] ) ==0

     which follows from the requirement: \sum_m T2C[m,i]*exp(i*m*phi) is real for any phi

     We can also derive this from the definition
        T2C[m,i]==<Y_{lm}|Phi_i>
     and the fact that localized orbital |Phi> is real, hence
       T2C[m,i]^* = <Y_{lm}|Phi_i>^* = <Y_{lm}^*|Phi_i> = (-1)^m <Y_{l,-m}|Phi_i> = (-1)^m T2C[l,-m]

     We are free to add any phase to cubic harmonics, hence T2C[m,i] -> T2C[m,i]*exp(i*phi_i)
       with phi_i being arbitrary

     This leads to the following 2x2 system of equations:
      ( Rp[m,i], Qp[m,i] ) ( sin(phi_i) )
      ( Qm[m,i], Rm[m,i] ) ( cos(phi_i) ) = 0

     where
          Qp[m,i] = Imag( T2C[m,i] + (-1)**m * T2C[-m,i] )
          Rm[m,i] = Real( T2C[m,i] - (-1)**m * T2C[-m,i] )
          Rp[m,i] = Real( T2C[m,i] + (-1)**m * T2C[-m,i] )
          Qm[m,i] = Imag(-T2C[m,i] + (-1)**m * T2C[-m,i] )
    """

    for i in range(2*l+1):
        ctg=None
        indexm = list(range(l+1))
        indexm = sorted(indexm, key=lambda m: abs(T2C[m+l,i])+abs(T2C[-m+l,i]), reverse=True)

        print('sorted by magnitude', file=fh_info)
        for m in indexm:
            print(m, T2C[m+l,i],T2C[-m+l,i], file=fh_info)

        for m in indexm[:1]:  # We will do it such that only the largest component is OK
            Qp = T2C[m+l,i].imag + (-1)**m * T2C[-m+l,i].imag
            Rm = T2C[m+l,i].real - (-1)**m * T2C[-m+l,i].real
            Rp = T2C[m+l,i].real + (-1)**m * T2C[-m+l,i].real
            Qm =-T2C[m+l,i].imag + (-1)**m * T2C[-m+l,i].imag
            #print 'Qp=', Qp, 'Rm=', Rm
            if abs(Qp) > small or abs(Rm) > small:
                if abs(Qp) > small :
                    ctg = -Rp/Qp
                    xb = -Rp
                    xa = Qp
                if abs(Rm) > small :
                    ctg = -Qm/Rm
                    xb = -Qm
                    xa = Rm
                #break

        if ctg is not None:
            phi = numpy.arctan2(xa, xb)
            print('Correcting T2C because original cubic harmonics were not real. Vector=', i, 'phase=', phi, file=fh_info)
            print('T2C^T before correction:', file=fh_info)
            cprint2(fh_info, T2C.transpose())
            T2C[:,i] = T2C[:,i] * numpy.exp(phi*1j)
            print('T2C^T after correction:', file=fh_info)
            cprint2(fh_info, T2C.transpose())

            for m in range(0,l+1):
                Rm = T2C[m+l,i].real - (-1)**m * T2C[-m+l,i].real
                Qp = T2C[m+l,i].imag + (-1)**m * T2C[-m+l,i].imag
                #Rp = T2C[m+l,i].real + (-1)**m * T2C[-m+l,i].real
                #Qm =-T2C[m+l,i].imag + (-1)**m * T2C[-m+l,i].imag
                #print 'm=', m, 'Rm=', Rm, 'Qp=', Qp

                if abs(Rm)>10*small or abs(Qp)>10*small:
                    print('ERROR: Could not find an angle to make all cubic harmonics real for vector=', i, 'and m=', m, 'D[Real]=', Rm, 'D[Imag]=', Qp)
                    print('ERROR: Could not find an angle to make all cubic harmonics real for vector=', i, 'and m=', m, 'D[Real]=', Rm, 'D[Imag]=', Qp, file=fh_info)

def cprint2(fh, U):
    for i in range(numpy.shape(U)[0]):
        for j in range(numpy.shape(U)[1]):
            f = U[i,j]
            if abs(f)<1e-10: f = 0j
            print("%12.8f %12.8f*i " % (f.real, f.imag), end=' ', file=fh)
        print(file=fh)
    print(file=fh)

def CoulombMatrix(params,ftrans,l,fh_info):
    """
    Generates the general Coulomb interaction matrix-elements U_{m1,m2,m3,m4}
    in the spin independent basis. Contains all the terms, the difference
    between full U or Ising type is taken care of in the EDQCM solver.

    We use here the implementation of atom_d.py.
    """

    Q3d=None
    small_t2c=1e-5

    ### IF Trans.dat exists read it an obtain transformation
    ### to the local basis, otherwise construct the transformation
    ### here.
    if (len(glob.glob(ftrans))!=0):
        print('Reading file', ftrans, file=fh_info)
        (Sigind, CF) = ReadTrans(ftrans, fh_info)

        if (len(Sigind)!= 2*(2*l+1) and len(Sigind)!= (2*l+1)):
            print('INORG ERROR: Sigind from file '+ftrans+' does not have correct dimension!')
            sys.exit(1)

        if len(CF)==2*l+1:
            if Q3d==None: Q3d=True # this must be 3d orbital
        elif len(CF)==2*(2*l+1):
            dn=2*l+1
            off_diag = CF[dn:,:dn]
            off = sum(sum(abs(off_diag)))
            if Q3d==None:
                if off>1e-5:
                    Q3d=False
                else:
                    Q3d=True
        else:
            print('INORG ERROR: Transformation CF=T2C does not have correct dimension')
            sys.exit(1)

        # Matrix contained in Trans.dat should be T(i,m):
        #   - i runs over real harmonics (z^2, x^2-y^2, yz, xz, xy)
        #   - m runs over complex harmonics (-2, -1, 0, 1, 2)
        # The matrix in Trans.dat, read into CF, is the transpose of
        # what we want in T2C, hence the T2C = transpose(CF) statement
        if Q3d:
            ##### change 2013
            T2C = numpy.transpose(CF[:(2*l+1),:(2*l+1)])
            Check_T2C_Real(T2C, l, fh_info, small_t2c)
        else:
            #CFN = zeros(shape(CF),dtype=complex)
            #CFN[:,:(2*l+1)] = CF[:,(2*l+1):]
            #CFN[:,(2*l+1):] = CF[:,:(2*l+1)]
            CFN = CF
            T2C = numpy.transpose(CFN)
        ##### change 2013, you will need to generalize this

    else:
        """
        Cubic harmonics have the order compatible with the Wien2K package:
        z^2, x^2-y^2, yz, zx, xy
        """
        print(ftrans, 'file not found; generating Sigind & complex-to-real spherical harmonics transformation inside atom_d.')
        T2C = numpy.transpose(Spheric2Cubic(l))
        if Q3d==None: Q3d=True

        Sigind = numpy.zeros((2*(2*l+1),2*(2*l+1)), dtype=int)
        if Q3d:
            for i in range(2*l+1):
                Sigind[i,i] = i+1
                Sigind[i+(2*l+1),i+(2*l+1)] = i+1
        else:
            zr = numpy.zeros(2*l+1)
            CF = numpy.transpose(T2C)
            CFn=[]
            for i in range(2*l+1):
                CFn.append( CF[i].tolist()+zr.tolist()  )
                CFn.append( zr.tolist()+CF[i].tolist()  )
                Sigind[2*i,2*i] = i+1
                Sigind[2*i+1,2*i+1] = i+1
            CFn = numpy.array(CFn)
            T2C = numpy.transpose(CFn)


    ### Print some information to fh_info file
    print('T2C=', file=fh_info)
    for i in range(len(T2C)):
        for j in range(len(T2C)):
            print("%6.3f %6.3f   " % (T2C[i,j].real, T2C[i,j].imag), end=' ', file=fh_info)
        print(file=fh_info)


    print('Impurity level structure Sigind is:', file=fh_info)
    print(Sigind, file=fh_info)

    print('T2C^T follows:', file=fh_info)
    print('\n'.join('   '.join('%10f %10f' % (x.real, x.imag) for x in row) for row in T2C.transpose()), file=fh_info)
    print('shape(T2C)=', numpy.shape(T2C), file=fh_info)
    print('T2C is Unitary=', sum(abs(numpy.matrix(T2C) * numpy.matrix(T2C).H - numpy.identity(len(T2C)))), file=fh_info)

    ### get slater integrals
    Fk = SlaterF(params['U'][0], params['J'][0], l)
    print('Slater integrals F^k are ', Fk[:,l], file=fh_info)

    ### Get first Coulomb matrix elements
    # These are possible values for CoulombF:
    #   'Ising'  : density-density interaction of Slater type
    #   'Full'   : Fully rotationally invariant interaction of Slater type, but with SU(N) symmetry [like in Kanamori]
    #   'FullS'  : Fully rotationally invariant interaction of Slater type, but with SU(3) symmetry [usual Slater interaction]
    #   'FullK'  : Kanamori type interaction, but U and J are computed so that within t2g block the interaction has the same form as Slater interaction

    if params['CoulombF'][0].replace("'","") in ['Full','Ising','FullS', 'IsingS']:      # Slater parametrization
        UC = CoulUsC2(l,T2C)              # Coulomb repulsion matrix -> in this case we keep all terms, no need to care about sign problem
    elif params['CoulombF'][0].replace("'","") in ['FullK','IsingK']:  # Kanamori parametrization
        UC = numpy.zeros((l+1, 2*l+1, 2*l+1, 2*l+1, 2*l+1), dtype=complex)
        for i in range(5):
            for j in range(5):
                if i==j:
                    UC[0,i,j,j,i]=1.
                    UC[1,i,j,j,i]=4./49.
                    UC[2,i,j,j,i]=4./49.
                else:
                    UC[0,i,j,j,i]=1.      # F0
                    UC[1,i,j,j,i]=-2./49. # F2 exact for t2g's
                    UC[2,i,j,j,i]=-4./441.# F4 exact for t2g's
                    UC[1,i,j,i,j]= 3./49. # F2 exact for t2g's
                    UC[2,i,j,i,j]=20./441.# F4 exact for t2g's
                    UC[1,i,i,j,j]= 3./49. # F2 exact for t2g's
                    UC[2,i,i,j,j]=20./441.# F4 exact for t2g's
    else:
        print('ERROR: CoulombF form ', params['CoulombF'][0], ' not yet implemented!')
        sys.exit(0)

    ###Get interaction Matrix and write it to file
    Uc = numpy.zeros(numpy.shape(UC)[1:],dtype=complex)
    for m1 in range(numpy.shape(Uc)[0]):
        for m2 in range(numpy.shape(Uc)[1]):
            for m3 in range(numpy.shape(Uc)[2]):
                for m4 in range(numpy.shape(Uc)[3]):
                    # note that we don't need the 1/2 factor in front due to the
                    # convention used in pyqcm
                    Uc[m1,m2,m3,m4] = numpy.sum(numpy.array(UC[:,m1,m2,m3,m4])*numpy.array(Fk[:l+1,l]))

    return Uc

def ClusterCoulombMatrix(params,dir,fh_info,cix,interaction_spec=None):
    """
    Reconstructs the Coulomb matrix for a cluster calculation
    """

    # First read the case.indmfl file
    case = get_case()
    try:
        f = open(case+'.indmfl')
        indmfl = f.readlines()
        f.close()
    except ValueError:
        print("File "+case+'.indmfl'+" not found!")

    # Extract info from atoms of the cluster cix
    nat, at_l = GetClusterAtoml(indmfl,cix)
    print("We extracted the following l for the cluster atoms: ", at_l,file=fh_info)

    # Check that Sigind's size corresponds
    (Sigind, CF) = ReadTrans(dir+'Trans.dat', fh_info)
    if len(Sigind)!=numpy.sum([2*l+1 for l in at_l]) and len(Sigind)!=numpy.sum([2*(2*l+1) for l in at_l]):
        raise ValueError("Expected size for Sigind is different from the one in Trans.dat file!")
    spin_deg=1
    if len(Sigind)==numpy.sum([2*(2*l+1) for l in at_l]):   # if spin_dep
        spin_deg=2

    # Construct each Trans_local_i.dat and remember its global-orbital mapping.
    # The legacy path also builds a Slater tensor for each atom.  The explicit
    # Uc path must bypass SlaterF: SlaterF(Uc,J) gives
    # U_x2-y2=Uc+8J/7 rather than the requested literal onsite Uc.
    if interaction_spec is None:
        interaction_spec = resolve_interaction_spec(
            params, solver='iNORG', cluster=True)
    Uc_cluster = {}
    mem_indx = {}
    for iat in range(nat):
        print("Constructing Coulomb matrix for atom ", iat, file=fh_info)
        ftrans, mem_indx[iat], spin_deg = ConstructLocalTrans(iat,at_l,dir,Sigind,CF,spin_deg)
        if interaction_spec.mode == 'legacy_slater':
            Uc_cluster[iat] = CoulombMatrix(
                params, ftrans, at_l[iat], fh_info)

    # Reconstruct the full interaction matrix.
    print("Reconstructing Coulomb matrix for cluster ", iat, file=fh_info)
    if interaction_spec.mode == 'legacy_slater':
        if spin_deg==2:
            U_tensor = ReconstructCoulombMatrix(
                Uc_cluster, mem_indx, int(len(Sigind)/2))
        else:
            U_tensor = ReconstructCoulombMatrix(
                Uc_cluster, mem_indx, len(Sigind))
    elif interaction_spec.mode == 'single_band_hubbard':
        U_tensor = build_single_band_cluster_tensor(
            interaction_spec, Sigind, mem_indx, nat, spin_deg)
        print('Using literal single-band cluster interaction Uc=',
              interaction_spec.value, 'source=', interaction_spec.source,
              file=fh_info)
    else:
        raise ValueError(
            'iNORG: unsupported cluster interaction mode '
            + repr(interaction_spec.mode))

    return U_tensor

def ReconstructCoulombMatrix(Uc_cluster,mem_indx,size):
    """
    From each local Uc matrix, reconstruct the full matrix for the cluster.
    We take care of the reordering of each orbital in Siging/CF with the mem_indx dict.
    """
    Uc = numpy.zeros([size,size,size,size],dtype=complex)

    def find_at(iorb):
        indx = numpy.unravel_index(numpy.argmin(numpy.abs(numpy.array([*mem_indx.values()])-iorb)),numpy.shape(numpy.array([*mem_indx.values()])))
        return indx[0],indx[1]

    for iorb1 in range(size):
        at1,ic1 = find_at(iorb1)
        for iorb2 in range(size):
            at2,ic2 = find_at(iorb2)
            if at1==at2:                 #Only local Coulomb interaction
                for iorb3 in range(size):
                    at3,ic3 = find_at(iorb3)
                    if at2==at3:                 #Only local Coulomb interaction
                        for iorb4 in range(size):
                            at4,ic4 = find_at(iorb4)
                            if at3==at4:                 #Only local Coulomb interaction
                                Uc[iorb1,iorb2,iorb3,iorb4] = Uc_cluster[at1][ic1,ic2,ic3,ic4]
    return Uc

def ConstructLocalTrans(iat,at_l,dir,Sigind,CF,spin_deg):
    """
    Constructs the local Trans file for atom iat in the cluster,
    which will enable to use the CoulombMatrix routine
    """

    n_elem = spin_deg*(2*at_l[iat]+1)                   # target number of lines for local trans
    n_tot = spin_deg*int(numpy.sum([(2*l+1) for l in at_l]))     # total number of lines
    CF_loc = numpy.zeros([n_elem,n_elem],dtype=complex)
    Sigind_loc = numpy.zeros([n_elem,n_elem],dtype=int)

    #### Find corresponding CF elements
    indx_min = int(numpy.sum([(2*l+1) for l in at_l[:iat]]))    #need to force int here
    indx_max = indx_min+int(n_elem/spin_deg)

    # we check for each CF raw if it corresponds to atom iat
    cnt=0
    mem_indx = []

    # K-native cluster fix: the cftrans is the DELOCALIZED plaquette-Fourier(x)T2C
    # (each K-row carries weight on ALL atoms), so the per-atom row extraction below
    # -- which assumes a site-localized, block-diagonal cftrans -- overflows CF_loc
    # (cnt >> n_elem).  The on-site Coulomb is LOCAL: build it from the standard
    # d-shell T2C (the Kronecker factor of build_plaquette_cftrans, dx2y2 first),
    # with x2-y2 correlated and a site-ordered mem_indx.  The resulting site-basis
    # Vijkl is mapped to K downstream by transform_interaction_to_K_vijkl (kept ON
    # for K-native), giving the SAME solver interaction as the legacy site path.
    _sigA = numpy.array(Sigind)
    _hh = _sigA.shape[0] // 2
    _ppA = _sigA[:_hh, :_hh] if _hh >= 1 else _sigA
    _is_k_native = (numpy.count_nonzero(_ppA - numpy.diag(numpy.diag(_ppA))) == 0)
    if _is_k_native and spin_deg == 1 and at_l[iat] == 2:
        _s = numpy.sqrt(0.5)
        CF_loc = numpy.array([[   _s,     0, 0,     0,    _s],
                              [    0,     0, 1,     0,     0],
                              [    0,    _s, 0,   -_s,     0],
                              [    0, 1j*_s, 0, 1j*_s,     0],
                              [-1j*_s,     0, 0,     0, 1j*_s]], dtype=complex)
        Sigind_loc = numpy.zeros([n_elem, n_elem], dtype=int)
        Sigind_loc[0, 0] = 1                          # x2-y2 correlated (matches d08 Trans_loc)
        nloc = int(n_elem / spin_deg)
        nat = len(at_l)                               # number of cluster atoms (= N_K, 4 for the plaquette)
        # Interleave so the 4 K-correlated rows (diagonal Sigind 0..3, the 4 K-points of
        # x2-y2) map to 4 DISTINCT atoms, each at local ic=0 (x2-y2) -- reproducing d08's
        # distinct-atom routing.  The buggy [iat*nloc+k] routed all 4 correlated rows into
        # atom-0's d-shell, so ReconstructCoulombMatrix.find_at pulled atom-0's INTER-ORBITAL
        # Coulomb (U-2J) instead of the 4 sites' on-site U -> full vij_abba instead of diag(U)
        # -> Vijkl ~U (256 entries) instead of the correct non-diagonal K interaction ~U/N_K (64).
        mem_indx = [iat + nat * k for k in range(nloc)]   # atom iat -> [iat, iat+nat, ...]
    elif _is_k_native and spin_deg == 2 and at_l[iat] == 2:
        # K-native NAMBU (spin_deg==2): the 40x40 doubled Sigind is pp(rows 0-19) +
        # hh(rows 20-39); the delocalized cftrans is block_diag(plaq(x)T2C, plaq(x)T2C),
        # so the else-branch row-scan over-extracts (same bug as normal).  Build the
        # 10x10 local CF = block_diag(T2C, T2C) (verified: d08 Trans_loc hh block == T2C,
        # NOT conj), and interleave each atom's d-shell within BOTH the pp and hh blocks
        # so the 4 atoms' x2-y2 land at the pp correlated rows 0-3 and hh rows N_pp+0..3.
        # Sigind_loc (derived from mem_indx) then carries pp/hh/Phi -> CoulombMatrix +
        # prepare_EDinput yield the 8x8 site vij_abba=diag([U,U,U,U,0,0,0,0]); the io-layer
        # Nambu transform (block_diag(U,U), sign=-1) maps it to the 64-entry -U/N_K K interaction.
        _s = numpy.sqrt(0.5)
        _T2C = numpy.array([[   _s,     0, 0,     0,    _s],
                            [    0,     0, 1,     0,     0],
                            [    0,    _s, 0,   -_s,     0],
                            [    0, 1j*_s, 0, 1j*_s,     0],
                            [-1j*_s,     0, 0,     0, 1j*_s]], dtype=complex)
        CF_loc[:5, :5] = _T2C
        CF_loc[5:, 5:] = _T2C                         # block_diag(T2C, T2C): pp + hh d-shells
        nloc = int(n_elem / spin_deg)                 # 5 (per pp/hh block)
        nat = len(at_l)                               # 4 cluster atoms
        N_pp = nat * nloc                             # 20 (pp block size in the 40x40 Sigind)
        mem_indx = ([iat + nat * k for k in range(nloc)]              # pp d-shell, interleaved
                    + [N_pp + iat + nat * k for k in range(nloc)])     # hh d-shell, interleaved
        for i1 in range(n_elem):
            for i2 in range(n_elem):
                Sigind_loc[i1, i2] = Sigind[mem_indx[i1], mem_indx[i2]]
    else:
        for i in range(len(CF)):
            # if spin polarized, need to shift the indx_min/max when
            # changing from spin up lines to spin down
            if (spin_deg==2 and i>=int(n_tot/2)):
                shift = int(n_tot/2)
            else:
                shift=0

            if numpy.sum([abs(x) for x in CF[i,indx_min+shift:indx_max+shift]])>1e-8:   # we know that iat is the one with nonzero elements in these columns
                if spin_deg==2:
                    CF_loc[cnt,:] = numpy.concatenate((CF[i,indx_min:indx_max],CF[i,indx_min+int(n_tot/2):indx_max+int(n_tot/2)]))
                else:
                    CF_loc[cnt,:] = CF[i,indx_min:indx_max]
                mem_indx.append(i)
                cnt+=1

        #### Now construct local Sigind
        for i1 in range(n_elem):
            for i2 in range(n_elem):
                Sigind_loc[i1,i2] = Sigind[mem_indx[i1],mem_indx[i2]]

    #### We left with creating file
    ftrans = dir+'Trans_loc'+str(iat)+'.dat'
    f = open(ftrans,'w')
    # first size of Sigind and CF
    print(len(Sigind_loc), len(CF_loc), " #  size of Sigind and CF", file=f)
    # Then Sigind
    print("#---- Sigind follows", file=f)
    for i1 in range(n_elem):
        for i2 in range(n_elem):
            f.write(' '+str(Sigind_loc[i1,i2]))
        f.write('\n')
    # Then CF
    print("#---- CF follows",file=f)
    for i1 in range(n_elem):
        for i2 in range(n_elem):
            f.write(' '+str(CF_loc[i1,i2].real))
            f.write(' '+str(CF_loc[i1,i2].imag))
        f.write('\n')
    f.close()

    return ftrans, mem_indx, spin_deg

def GetClusterAtoml(indmfl,cix):
    """
    Extract the l of each atom in cluster
    """
    def NextIs(line):
        locrot = int(line.strip().split()[2])
        nl = int(line.strip().split()[1])-1
        if locrot==0:
            return 2+nl
        elif locrot==-1:
            return 5+nl
        elif locrot==-4:
            return 6+nl
        else:
            raise ValueError("locrot "+str(locrot)+" not recognized when building cluster interaction matrix!")
    # total number of atoms
    natom = int(indmfl[2].strip().split()[0])

    at_l = []
    lindx=3
    # Loop over atoms and extract those with right cix
    for i in range(natom):
        # check each l of this atom
        for il in range(int(indmfl[lindx].strip().split()[1])):
            # check cix
            if int(indmfl[lindx+1+il].strip().split()[2])==cix:
                at_l.append(int(indmfl[lindx+1+il].strip().split()[0]))
        # calculate next line index
        lindx+=NextIs(indmfl[lindx])

    return len(at_l), at_l

def prepare_EDinput(f, Sigind, Eimps, Delta, lom, U_tensor, CoulombF, spin_dep, fh_info, nambu=False):
    """
    Function that writes the impurity energies, hybridization function,
    and interaction matrix in the ED solver format.

    In Nambu mode, Eimp and Delta use the full Nambu Sigind, but U_tensor remains
    a physical (normal-state) interaction and is indexed only with the pp block.
    """
    norb = numpy.count_nonzero(numpy.diag(Sigind))  #nbr of correlated orbitals (spin included if spin dependent)
    ntot = len(Sigind)

    # Nambu: U_tensor is a normal-state (pp-block only) interaction even though
    # Sigind/Eimp/Delta are 2x larger.
    Sigind_uc = Sigind
    ntot_uc = ntot
    if nambu:
        ntot_uc = ntot // 2
        Sigind_uc = Sigind[:ntot_uc, :ntot_uc]

    ##### Info
    print("Writing the input for the ED solver in "+f, file=fh_info)

    ##### Open file
    f = open(f,'w')

    ##### Start with the impurity on-site energies
    f.write("Tij\n")
    # construct numpy array
    to_save = numpy.zeros([norb,norb],dtype=float)     # Force float type, though it should already be the case
    i1=0
    for iorb1 in range(ntot):
        i2=0
        if Sigind[iorb1,iorb1]!=0:     # keep only the non-zero terms
            for iorb2 in range(ntot):
                if Sigind[iorb1,iorb2]!=0:     # keep only the non-zero terms
                    to_save[i1,i2] = Eimps[Sigind[iorb1,iorb2]-1]
                    i2+=1
                elif Sigind[iorb2,iorb2]!=0:  #if off diagonal zero but iorb2 is correlated, still have to count the index i2
                    i2+=1
            i1+=1

    # save it
    numpy.savetxt(f,to_save)
    f.write("\n")


    ##### Then hybridization function
    f.write("Delta\n")
    # construct the right shape
    to_save = numpy.zeros((len(lom), (norb*norb + 1)), dtype=numpy.complex128)
    to_save[:,0] = lom
    for i in range(len(lom)):
        l = 1
        for iorb1 in range(ntot):
            for iorb2 in range(ntot):
                # keep only the terms allowed by Sigind
                if Sigind[iorb1,iorb2]!=0:
                    to_save[i,l] = Delta[Sigind[iorb1,iorb2]-1][i]
                    l += 1
                elif (Sigind[iorb2,iorb2]!=0 and Sigind[iorb1,iorb1]!=0):    #if off diagonal zero but iorb2 and iorb1 are correlated I have to count
                    l += 1
    # save it
    numpy.savetxt(f,to_save)
    f.write("\n")


    ##### Then the interaction tensor
    # We have to care here in case of spin-polarized calculations
    divider=1
    if spin_dep:
        divider=2
    print(CoulombF, file=fh_info)
    # if Ising, prepare the Vij format
    Vij_abba = numpy.zeros([int(norb/divider),int(norb/divider)],dtype=float)
    Vij_abab = numpy.zeros([int(norb/divider),int(norb/divider)],dtype=float)
    if CoulombF in ['Ising','IsingS', 'IsingK']:
        i1=0
        for m1 in range(int(ntot_uc/divider)):
            i2=0
            if Sigind_uc[m1,m1]!=0:
                for m2 in range(int(ntot_uc/divider)):
                    i3=0
                    if Sigind_uc[m2,m2]!=0:
                        for m3 in range(int(ntot_uc/divider)):
                            i4=0
                            if Sigind_uc[m3,m3]!=0:
                                for m4 in range(int(ntot_uc/divider)):
                                    if Sigind_uc[m4,m4]!=0:
                                        if (m1==m4 and m2==m3):
                                            if numpy.abs(U_tensor[m1,m2,m3,m4])>1e-8:
                                                Vij_abba[i1,i2] = U_tensor[m1,m2,m3,m4].real
                                        if (m1==m3 and m2==m4 and m1!=m2):
                                            if numpy.abs(U_tensor[m1,m2,m3,m4])>1e-8:
                                                Vij_abab[i1,i2] = U_tensor[m1,m2,m3,m4].real
                                        i4+=1
                                i3+=1
                        i2+=1
                i1+=1
        f.write("Vij_abba\n")
        numpy.savetxt(f, Vij_abba)
        f.write("\n")
        f.write("Vij_abab\n")
        numpy.savetxt(f, Vij_abab)
        f.write("\n")
    # if full, prepare Vijkl
    else:
        f.write("Vijkl\n")
        i1=0
        for m1 in range(int(ntot_uc/divider)):
            i2=0
            if Sigind_uc[m1,m1]!=0:
                for m2 in range(int(ntot_uc/divider)):
                    i3=0
                    if Sigind_uc[m2,m2]!=0:
                        for m3 in range(int(ntot_uc/divider)):
                            i4=0
                            if Sigind_uc[m3,m3]!=0:
                                for m4 in range(int(ntot_uc/divider)):
                                    if Sigind_uc[m4,m4]!=0:
                                        if numpy.abs(U_tensor[m1,m2,m3,m4])>1e-8:
                                            print(i1+1, i2+1, i3+1, i4+1, U_tensor[m1,m2,m3,m4].real, file=f)    #only real interaction terms
                                        i4+=1
                                i3+=1
                        i2+=1
                i1+=1

    #close file
    f.close()

def get_case():
    """
    Finds the casename of th calculation in a general way
    """
    case= os.path.basename(os.getcwd())
    if not (os.path.isfile(case+'.struct') and os.path.isfile(case+'.in0') ):  # directory name is not case (happens when submitting to cluster)
        files = glob.glob('*.struct')
        if len(files) < 1:
            raise Exception('No struct file present.')
        elif len(files) > 1:
            # heuristic algorithm to determine case:
            # need in0 and struct present with the same case.
            candidates = [os.path.splitext(f)[0] for f in files]
            for cand in candidates:
                if os.path.isfile(cand+'.in0'):
                    case = cand
                    break
        else: # just one candidate exists, hence it must be it
            case, ext = os.path.splitext(os.path.basename(files[0]))
    return case
