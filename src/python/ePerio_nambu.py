#!/usr/bin/env python

"""
Script ePerio_nambu.py

Extended version of ePerio.py supporting Nambu Green's functions for
superconducting cluster DMFT periodization.

Without --nambu: behaves identically to ePerio.py.
With --nambu:    reads Nambu G_k files (cixdim = 2*norb in header),
                 extracts normal G and anomalous F blocks, applies
                 periodize() independently to each, and writes both
                 to case.gunfold.

Wrapper handling the periodization of the 
correlated Green's function:
    1- call x_dmft.py routine for producing the momentum-resolved
    Green's function for the correlated subspace. It is an external
    call, the parallelization settings should be provided in 
    mpi_prefix.dat* files in the jobfile. 
    2-  reading the Green's function from the previous step, and
     periodizing it following Eq. (7) of [D. Sénéchal, arXiv:0806.2690].

MANDATORY ARGUMENTS:
--------------------

    CASE : str
            Name of the case, to find case.energy* files, etc.

OPTIONAL ARGUMENTS:
-------------------

    --gamma : float (default: 0.0025)
                Broadening for non-correlated bands.
    --gammac : float (default: 0.0025)
                Broadening for correlated bands.        
    --up : bool (default: False)
                If option is used, calculations carried for up spins
                only.   
    --down : bool (default: False)
                If option is used, calculations carried for down spins
                only.
    --so : bool (default: False)
                If option is used, spin-orbit is activated.
    -c : bool (default: False)
                If option is used, complex version is used for calls to 
                x_dmft.py.
    --no-dmft : bool (default: True)
                If option is used, the x_dmft.py routines will not be 
                called, only the calculation of hybridization function
                is carried. 
    -C : int (default: 0)
                Size of k-point chunks on which dmftgk is run, if 0 no 
                chunking.
    --nambu : bool (default: False)
                Enable Nambu mode for superconducting periodization.
                G_k files must have cixdim = 2*norb (Nambu doubled matrix).
                Periodizes normal (G) and anomalous (F) blocks independently.
                Only compatible with --no-dmft; run dmftgk_nambu.py first.

"""

import os
import sys
import argparse
import numpy as np
import subprocess
import json
import shutil
import glob


def RunXdmft(cmd: str, up: str, dn: str, so: str, c:str
             ) -> None:
    """
    General function to run x_dmft.py routine for the specific
    cmd task.

    Parameters:
    -----------
        cmd : str
            Command to run, can be: 'lapw1', 'dmft1', 'dmftu', 'dmftgk'.
        up : str
            If not empty, calculations carried for up spins only.
        dn : str
            If not empty, calculations carried for down spins only.
        so : str
            If not empty, spin-orbit is activated.  
        c : str
            If not empty, complex version is used.

    Raises:
    -------
        ValueError: if the cmd argument is not valid.
        Exit error if the subprocess call to x_dmft.py fails. 
    """

    # check that the command is valid
    if cmd not in ["lapw1", "dmft1", "dmftu -g", "dmftgk"]:
        raise ValueError(
            "ERROR: command "
            + cmd
            + " not valid."
            + " Should be either 'lapw1' or 'dmft1',"
            + " 'dmftu -g', 'dmftgk'."
        )

    # call x_dmft.py
    exe = ""
    exe += "x_dmft.py " + cmd + " -p"
    if so:
        exe += " --so"
    if up:
        exe += " --up"
    if dn:
        exe += " --dn"
    if c:
        exe += " -c"
    if subprocess.call(exe, shell=True) != 0:
        print("ERROR when calling 'x_dmft.py " + cmd + " -p'! Aborting..")
        sys.exit(0)


def JoinEnergyFiles(case, so="", up="", dn=""):
    """Join split energy files from parallel lapw1 into a single file.

    After lapw1 -p, energy files are split as case.energy_0, _1, etc.
    dmftgk expects a single case.energy file. This function joins them
    by keeping the header from _0 and appending k-point data from all
    ranks in order.
    """
    import re
    basename = case + ".energy" + so + up + dn
    split_files = sorted(
        glob.glob(basename + "_[0-9]*"),
        key=lambda f: int(f.rsplit("_", 1)[-1])
    )
    if not split_files:
        return  # no split files, single energy file should already exist

    # detect header length from the first file:
    # header = linearization energy lines (before first k-point header)
    # k-point header lines contain scientific notation (e.g. '0.000E+00')
    header_lines = 0
    with open(split_files[0], "r") as f:
        for line in f:
            if re.search(r'[0-9]E[+-]', line):
                break
            header_lines += 1

    # join: header from file_0 + k-point data from all files in rank order
    with open(basename, "w") as fout:
        for i, sf in enumerate(split_files):
            with open(sf, "r") as fin:
                for j, line in enumerate(fin):
                    if i == 0:
                        fout.write(line)
                    else:
                        if j >= header_lines:
                            fout.write(line)

    print("Joined", len(split_files), "split energy files into", basename)


def PrepareDmftgkInput(
    case: str, ncix: int, gamma: float, 
    gammac: float, up: str, dn: str, so: str
) -> None:
    """
    Prepares the input file dmftgk.in for the x_dmft.py dmftgk call.
    The number of sig.inp* files is given by ncix. The number of 
    output Gf files is set accordingly.

    Parameters:
    -----------
        case : str
            Name of the case.
        ncix : int
            Number of independent kcix blocks.
        up : str
            If not empty, calculations carried for up spins only.
        dn : str
            If not empty, calculations carried for down spins only.
        so : str
            If not empty, spin-orbit is activated.
    
    """

    # write dmftgk.in
    # we have to take care of the ordering of parameters
    # as it is note flxible within the fortran routine
    space = "     "
    with open("dmftgk.in", "w") as f:
        # mode
        f.write("g" + space + "# mode [g/e]: greens function or eigen values\n")
        # matsubara mode, always one in our case
        f.write("0" + space + "# matsubara mode\n")
        # LDA-energy-file
        f.write(
            case
            + ".energy"
            + so
            + up
            + dn
            + space
            + "# LDA energy filename, case.energy(so)(updn)\n"
        )
        # klist file
        f.write(case + ".klist" + space + "# klist filename\n")
        # reciprocal vectors
        f.write(case + ".rotlm" + space + "# reciprocal vectors filename\n")
        # projector filename
        f.write("Udmft." + space + "# projector filename\n")
        # gamma for non-correlated
        f.write(str(gamma) + space + "# gamma for non-correlated\n")
        # gamma for correlated
        f.write(str(gammac) + space + "# gamma for correlated\n")
        # self-energy filenames
        for ix in range(ncix):
            f.write("sig.inp" + str(ix + 1) + " ")
        f.write(space + "# self-energy filenames\n")
        # output gk
        for ix in range(ncix):
            f.write("G_k" + str(ix + 1) + " ")
        f.write(space + "# output gk files\n")
        # output glocal
        for ix in range(ncix):
            f.write("G_local" + str(ix + 1) + " ")
        f.write(space + "# output G_local files\n")
        # output glocal-diagonal
        for ix in range(ncix):
            f.write("g_local" + str(ix + 1) + " ")
        f.write(space + "# output diagonal g_local files\n")
    
    # here is a quick fix to avoid the problem of sig.inp19 which somehow
    # disapears sometimes. We copy it from sig.inp20 if it exists
    if os.path.isfile("sig.inp20") and not os.path.isfile("sig.inp19"):
        shutil.copy("sig.inp20","sig.inp19")


def Getncix(indmfl: list) -> int:
    """
    Extract the number of independent kcix blocks from the
    indmfl file.

    Parameters:
    -----------
        indmfl : list[str]
                 List of lines of the indmfl file.

    Returns:
    --------
        ncix : int
                Number of independent kcix blocks.
                 
    Raises:
    -------
        ValueError: if the a value of locrot is not recognized.
    """

    def NextIs(line: str) -> int:
        """
        Quick routine that provides the number of lines to skip
        given the locrot read in line.

        Parameters:
        -----------
            line : str
                    Line of the indmfl file.
        Raises:
        -------
            ValueError: if the a value of locrot is not recognized.
        """
        locrot = int(line.strip().split()[2])
        nl = int(line.strip().split()[1]) - 1
        if locrot == 0:
            return 2 + nl
        elif locrot == -1:
            return 5 + nl
        elif locrot == -4:
            return 6 + nl
        else:
            raise ValueError(
                "locrot "
                + str(locrot)
                + " not recognized when building cluster"
                + " interaction matrix!"
            )

    # total number of atoms
    natom = int(indmfl[2].strip().split()[0])

    lindx = 3
    # Loop over atoms to get final index
    for i in range(natom):
        lindx += NextIs(indmfl[lindx])

    # the ncix is located on line below
    ncix = int(indmfl[lindx + 1].strip().split()[0])

    return ncix


def SearchSigind(case: str, ncix: int) -> list:
    """
    Search all sigind tables in case.indmfl file. Number of 
    sigind tables should match with ncix. 

    Parameters:
    -----------
        case : str
            Name of the case.
        ncix : int
            Number of independent kcix blocks.
    
    Returns:
    --------
        sigind : list
            List of sigind tables. The tables are truncated to the
            non-zero (correlated) elements only.
        sigind_noshift : list
            List of sigind tables. The tables are truncated to the
            non-zero (correlated) elements only. The indices are
            not shifted to start from zero. 
        ucix : list
            List of indices of unique sigind tables.
        nucix : int
            Number of unique sigind tables.
    
    Raises:
    -------
        ValueError: if the number of sigind tables does not match
                    with ncix.
    """

    sigind = []
    sigind_noshift = []
    with open(case + ".indmfl", "r") as f:
        lines = f.readlines()
        # search all sigind tables
        # they are signaled by ... --- # Sigind follows -- statement
        for iline in range(len(lines)):
            if "Sigind follows" in lines[iline]:
                # read the sigind table
                sigind.append(ReadSigind(lines, iline)[0])
                sigind_noshift.append(ReadSigind(lines, iline)[1])
    
    # check nbr of sigind tables
    if len(sigind) != ncix:
        raise ValueError("number of sigind tables does not match with ncix!")
    
    # some sigind are repeated, find index of unique ones
    ucix = [0]
    for i in range(1,len(sigind_noshift)):
        equal = False
        for u in ucix:
            if np.array_equal(sigind_noshift[i], sigind_noshift[u]):
                equal = True
                break
        if not equal:
            ucix.append(i)    
    
    return sigind, sigind_noshift, ucix, len(ucix)


def ReadSigind(lines: list, iline: int) -> np.ndarray:
    """
    Reads sigind table located at iline in lines.

    Parameters:
    -----------
        lines : list[str]
            List of lines of the indmfl file.
        iline : int
            Line index where the sigind table starts.
    
    Returns:
    --------
        sigind : np.ndarray
            Sigind table truncated to the non-zero (correlated) elements
            only.
        sigind_noshift : np.ndarray
            Sigind table truncated to the non-zero (correlated) elements
            only. The indices are not shifted to start from zero.
    """

    # size of sigind table
    size = int(lines[iline-3].strip().split()[1])
    
    # read sigind table
    tmp_sigind = np.zeros((size, size), dtype=int)
    for i in range(size):
        tmp_sigind[i,:] = [eval(j) for j in lines[iline+i+1].strip().split()]
    
    #### restrict to non-zero elements
    #### this takes a few steps: look at diagonal elements
    #### and find non-zero elements
    #### and then assign values to new sigind table
    diag = np.diag(tmp_sigind)

    # find indices which are non-zero
    nonzer_indx = []
    for iorb in range(len(diag)):
        if diag[iorb]!=0:
            nonzer_indx.append(iorb)
    
    sigind = np.zeros([len(nonzer_indx),len(nonzer_indx)],dtype=int)
    # now assign values
    i1=0
    for iorb1 in nonzer_indx:
        i2=0
        for iorb2 in nonzer_indx:
            sigind[i1,i2]=tmp_sigind[iorb1,iorb2]
            i2+=1
        i1+=1

    # non-zero indices in these sigind tables should start
    # by zero since we will read data from sig.inp* files
    # generated from ssplit.py
    sigind_noshift = sigind.copy()
    sigind -= np.min(sigind[np.nonzero(sigind)])-1
    # replace all negative values by 0
    sigind[sigind<0] = 0

    return sigind, sigind_noshift


def ReadGf(filename: str, sigind: np.ndarray, nambu: bool = False) -> np.ndarray:
    """
    Read green's function from filename. We follow the sigind
    table and neglect potential non-zero (but small) off-diag 
    elements corresponding to zero sigind indices. 

    Parameters:
    -----------
        filename : str
            Name of the file containing the green's function.
        sigind : np.ndarray
            Sigind table (shape: norb x norb, where norb is the
            cluster size without Nambu doubling).
        nambu : bool (default: False)
            If True, the G_k file contains a Nambu-doubled matrix
            with cixdim = 2*norb in the header. Returns the full
            (2*norb, 2*norb) matrix; the caller extracts normal and
            anomalous blocks.
    
    Returns:
    --------
        gf : np.ndarray
            Green's function array.
            - Normal mode:  shape (nkpt, nfreq, norb, norb)
            - Nambu mode:   shape (nkpt, nfreq, 2*norb, 2*norb)
              where the upper-left block is G (normal) and
              upper-right block is F (anomalous).
    
    Raises:
    -------
        ValueError: if the number of columns do not correspond
                    with norb (normal) or 2*norb (Nambu). 

    Notes:
    ------

        We enforce the symmetry of the matrix to compensate small
        errors in the calculation of G(k,w). 

        In Nambu mode the file header contains cixdim = 2*norb.
        The data lines have 1 + 2*(2*norb)^2 entries per (ik,isymop,iw)
        triplet, encoding the full Nambu matrix in row-major order.

    """

    with open(filename, "r") as f:
        lines = f.readlines()
        # get nbr of k-points, frequencies, symop and cixdims
        nk = int(lines[0].strip().split()[1])
        nsymop = int(lines[0].strip().split()[2])
        nw = int(lines[0].strip().split()[3])
        cixdim = int(lines[0].strip().split()[4])
        norb = len(sigind)

        # In Nambu mode the file was written with cixdim = 2*norb.
        # In normal mode cixdim must equal norb.
        expected_cixdim = 2 * norb if nambu else norb
        if expected_cixdim != cixdim:
            if nambu:
                raise ValueError(
                    "Nambu mode: expected cixdim = 2*norb = "
                    + str(2 * norb)
                    + " in green's function file "
                    + filename
                    + ", but got cixdim = "
                    + str(cixdim)
                    + "."
                )
            else:
                raise ValueError(
                    "number of columns in green's function file"
                    + " does not match with nbr of columns in siding!"
                )

        # dim is the actual matrix dimension stored in the file
        # (= norb in normal mode, = 2*norb in Nambu mode)
        dim = cixdim

        # initialize gf with file matrix dimensions
        gf = np.zeros((nk, nsymop, nw, dim, dim), dtype=complex)
        # read gf
        line_cnt = 2
        for ik in range(nk):
            for isymop in range(nsymop):
                for iw in range(nw):
                    for iorb1 in range(dim):
                        for iorb2 in range(dim):
                            # we keep all terms, even those non-zero
                            # although corresponding sigind is non-zero
                            # for the renormalization factor calculated
                            # in CmpHybkw
                            ind = 2*(iorb1*dim+iorb2)+1
                            gf[ik,isymop,iw,iorb1,iorb2] = complex(
                                    eval(lines[line_cnt].strip().split()[ind]),
                                    eval(lines[line_cnt].strip().split()[ind+1])
                                    )
                    line_cnt += 1
    # sum over symetry operations
    gf = np.sum(gf,axis=1)/nsymop

    return gf


def ReadKlist(filename: str) -> np.ndarray:
    """
    Read the case.klist file to obtain the list of k-points.
    Returns the k-points in units of 1/a. 

    PARAMETERS:
    -----------
        case : str
            case name of the calculation. 
    
    RETURNS:
    --------
        klist : np.ndarray (nkpt, 4)
            List of k-points in units of 1/a, with relative weights. 

    RAISES:
    -------
        FileNotFoundError: if the case.klist file does not exist.
        
    """
    # check that the file exists
    if not os.path.isfile(filename):
        raise FileNotFoundError("ERROR: file " + filename
                                + " does not exist.")

    with open(filename, "r") as f:
        lines = f.readlines()
        # get nbr of kpoints
        for iline in range(len(lines)):
            if "END" in lines[len(lines)-1-iline]:
                # then las k-point is just above
                nkpt = len(lines)-1-iline
                break

        # initialize klist and read
        klist = np.zeros((nkpt, 4), dtype=float)
        klabel = []
        for ik in range(nkpt):
            # have to check if name if defined or not
            if len(lines[ik].strip().split())==6:
                shift=1
                klabel.append(lines[ik].strip().split()[0])
            else:
                shift=0
                klabel.append("")
            # norm & weight
            norm, klist[ik,0] = [eval(i) for i in 
                                 lines[ik].strip().split()[3+shift:5+shift]]
            # kpts
            klist[ik, 1:] = [2*np.pi*eval(i)/norm for i in 
                             lines[ik].strip().split()[0+shift:3+shift]]
            
    return klist, klabel


def periodize(gf: np.ndarray, k_unfold: list, cluster: list) -> np.ndarray:
    """
    Periodize the provided Green's function following Eq. (7) of
    [D. Sénéchal, arXiv:0806.2690].
    """
    # initialize Gper
    Gper = np.zeros(np.shape(gf)[:2],dtype=complex)
    # loop over kpts
    for ik in range(len(k_unfold)):
        for iR1 in range(len(cluster)):
            for iR2 in range(len(cluster)):
                # compute the components of k*(R1-R2)
                krx = k_unfold[ik,1]*(cluster[iR1][0]-cluster[iR2][0])
                kry = k_unfold[ik,2]*(cluster[iR1][1]-cluster[iR2][1])
                krz = k_unfold[ik,3]*(cluster[iR1][2]-cluster[iR2][2])
                # compute the phase factor
                phase = np.exp(-1j*(krx+kry+krz))
                # add to Gper
                Gper[ik,:] += phase*gf[ik,:,iR1,iR2]

    return Gper/len(cluster)


def save_Gper(ucix, nucix, Gper, wmesh, case, nkpt_tab, Fper=None):
    """
    Saves the unfolded Green's functions into case.gunfold file.

    Parameters:
    -----------
        ucix : list
            List of unique cix indices.
        nucix : int
            Number of unique cix blocks.
        Gper : dict
            Dictionary of periodized normal Green's functions,
            keyed by cix index.
        wmesh : np.ndarray
            Frequency mesh array.
        case : str
            Name of the case.
        nkpt_tab : int
            Number of k-points in the tabulation (from case.inperio).
        Fper : dict or None (default: None)
            If not None, dictionary of periodized anomalous (Nambu F)
            Green's functions, keyed by cix index. Written after Gper
            sections.
    """
    with open(case+".gunfold","w") as f:

        # first write 2 lines of information
        f.write(str(nucix)
                +' '+str(np.shape(Gper[ucix[0]])[0])
                +' '+str(np.shape(Gper[ucix[0]])[1])
                + "\t# nucix, nkpt, nfreq\n")
        f.write(str(nkpt_tab)+"\t# nkpt_tab\n")
        f.write("\n")

        # Then write down the frequency mesh
        f.write("# wmesh\n")
        np.savetxt(f,wmesh)
        f.write("\n")

        # then write down the normal Green's function
        for icix in range(nucix):
            f.write("# Gper, icix="+str(ucix[icix])+"\n")
            np.savetxt(f,Gper[ucix[icix]])
            f.write("\n")

        # if Nambu mode, also write anomalous (F) block
        if Fper is not None:
            for icix in range(nucix):
                f.write("# Fper (anomalous), icix="+str(ucix[icix])+"\n")
                np.savetxt(f,Fper[ucix[icix]])
                f.write("\n")


def Prepare_klist_chunk(case: str, chunks: list, ichunk: int
                        ) -> None:
    """
    Prepare the klist file for the specific chunk.

    Parameters:
    -----------
        case : str
            Name of the case.
        chunks : list
            List of chunks sizes.
        ichunk : int
            Index of the chunk to prepare.
    
    Raises:
    -------
        FileNotFoundError: if the case.klist_save file does not exist.
    """

    # check that the file exists
    if not os.path.isfile(case+".klist_save"):
        raise FileNotFoundError("ERROR: file " + case+".klist_save"
                                + " does not exist.")

    # we simply copy the lines corresponding to the chunk from
    # file_klist to case.klist
    with open(case+".klist_save", "r") as f:
        lines = f.readlines()

    # start index
    if ichunk == 0:
        start = 0
    else:
        start = np.sum(chunks[:ichunk])
    
    # open case.klist
    with open(case+".klist","w") as f:
        # write the chunk
        for i in range(chunks[ichunk]):
            f.write(lines[start+i])
        # add "END"
        f.write("END\n")


if __name__ == "__main__":
    """
    Main part of the script.

    Handles the reading of arguments, calls to the necessary
    routines to periodize the cluster Green's function.

    In Nambu mode (--nambu):
        - Only --no-dmft is supported (raise error if --dodmft used).
        - G_k files must be pre-computed by dmftgk_nambu.py with
          cixdim = 2*norb in their headers.
        - ReadGf returns (nkpt, nfreq, 2*norb, 2*norb).
        - Upper-left block  [:norb, :norb]  -> normal G.
        - Upper-right block [:norb, norb:]  -> anomalous F.
        - periodize() applied independently to each block.
        - Both Gper and Fper written to case.gunfold.

    Raises:
    -------
        FileNotFoundError: if the CASE.indmfl file does not exist.
        FileNotFoundError: in --no-dmft mode, if the G_k files do 
                            not exist.
        ValueError: if --nambu is combined with --dodmft (not supported).
    """

    #### Check arguments provided  ###################################
    parser = argparse.ArgumentParser(
        prog="ePerio_nambu",
        usage="%(prog)s CASE [options]",
        description=(
            "Periodization of cluster DMFT Green's functions. "
            "Extends ePerio.py with Nambu (superconducting) support."
        ),
    )
    # mandatory arguments
    parser.add_argument("CASE", type=str, help="case name of the calculation")
    # optional arguments
    parser.add_argument(
        "--gamma",
        type=float,
        default=0.0025,
        help="broadening for non-correlated bands",
    )
    parser.add_argument(
        "--gammac",
        type=float,
        default=0.0025,
        help="brodening for correlated bands",
    )
    parser.add_argument(
        "--so",
        dest="so",
        action="store_true",
        default=False,
        help="switches on spin-orbit coupling",
    )
    parser.add_argument(
        "--up",
        dest="up",
        action="store_true",
        default=False,
        help="magnetic LDA calculation with vector-up first",
    )
    parser.add_argument(
        "--dn",
        dest="dn",
        action="store_true",
        default=False,
        help="magnetic LDA calculation with vector-dn first",
    )
    parser.add_argument(
        "-c",
        dest="c",
        action="store_true",
        default=False,
        help="Complex version",
    )
    parser.add_argument(
        "--no-dmft",
        dest="dodmft",
        action="store_false",
        default=True,
        help="If used, x_dmft routines are not called",
    )
    parser.add_argument(
        "-C",
        dest="chunks",
        type=int,
        default=0,
        help="Size of k-point chunks on which dmftgk is run, if 0 no chunking",
    )
    parser.add_argument(
        "--nambu",
        dest="nambu",
        action="store_true",
        default=False,
        help=(
            "Enable Nambu mode for superconducting periodization. "
            "G_k files must have cixdim = 2*norb (pre-computed by "
            "dmftgk_nambu.py). Periodizes normal G and anomalous F "
            "blocks independently and writes both to case.gunfold. "
            "Only compatible with --no-dmft."
        ),
    )
    args = parser.parse_args()

    # Nambu mode is incompatible with running dmftgk internally.
    # The user must pre-compute G_k files with dmftgk_nambu.py.
    if args.nambu and args.dodmft:
        raise ValueError(
            "ERROR: --nambu mode only supports --no-dmft. "
            "Run dmftgk_nambu.py separately first to produce "
            "the Nambu G_k files, then re-run with --no-dmft --nambu."
        )

    # check that if chunks is not zero, the do-dmft option is not used
    if args.chunks != 0 and not args.dodmft:
        raise ValueError("ERROR: --no-dmft option cannot be used"
                         + " with --chunks option.")
    

    #### Check that the case exists ##################################
    if not os.path.isfile(args.CASE + ".indmfl"):
        raise FileNotFoundError("ERROR: case " + args.CASE 
                                + ".indmfl does not exist.")
    # check in CASE.indmfl that matsubara mode is set to 1
    with open(args.CASE + ".indmfl", "r") as f:
        if int(f.readlines()[1].strip().split()[0]) != 0:
            raise ValueError(
                "ERROR: matsubara mode is not set to 0 in "
                + args.CASE
                + ".indmfl file."
                + " Periodization only meaningful on real frequency"
                + " axis, i.e. in matsubara mode 0."
            )
    
    # search cix in case.indmfl
    with open(args.CASE + ".indmfl", "r") as f:
        lines = f.readlines()
        # extract ncix
        ncix = Getncix(lines)

    # search sigind tables and unique cix
    sigind, sigind_noshift, ucix, nucix = SearchSigind(args.CASE, ncix)

    # load json case.inperio file for cluster atomic positions
    # and klist file
    if not os.path.isfile(args.CASE + ".inperio"):
        raise FileNotFoundError("ERROR: case " + args.CASE 
                                + ".inperio does not exist.")
    with open(args.CASE + ".inperio", "r") as f:
        inperio = json.load(f)
    # load klist
    k_unfold, klabel = ReadKlist(inperio['k_unfold_file'])
    # copy the klist file to klist_save if the latter does not exist
    if not os.path.isfile(args.CASE+".klist_save"):
        shutil.copy(args.CASE+".klist",args.CASE+".klist_save")

    #### Set up, dn and so strings ###################################
    up = ""
    if args.up:
        up = "up"
    dn = ""
    if args.dn:
        dn = "dn"
    so = ""
    if args.so:
        so = "so"
    c = ""
    if args.c:
        c = "c"

    if args.dodmft:
        # Normal (non-Nambu) dodmft path — identical to ePerio.py.
        # (args.nambu and args.dodmft together would have raised above.)

        # find the number of chunks and size of last one if args.chunks
        # is not a multple of number of kpts
        if args.chunks != 0:
            # number of chunks
            nchunks = int(len(k_unfold)/args.chunks)
            chunks = [args.chunks for i in range(nchunks)]
            # size of last chunk
            last_chunk = len(k_unfold)%args.chunks
            # if last chunk is non-zero, we add it to the list
            if last_chunk != 0:
                chunks.append(last_chunk)
                nchunks += 1
        else:
            nchunks = 1
            chunks = [len(k_unfold)]

        # we also have to initialize the Gf arrays because we will
        # read it chunk by chunk
        # We need the number of frequency points, so we read sig.inp
        wmesh = np.loadtxt("sig.inp", usecols=(0,), dtype=float)[:]
        gf = {}
        for icix in ucix:
            # have to add +1 to len(wmesh) because ssplit
            # will add one additional freq far away
            gf[icix] = np.zeros((len(k_unfold),len(wmesh)+1,len(sigind[icix]),
                                len(sigind[icix])),dtype=complex)

        # loop over chunks and run x_dmft.py routines for each 
        for ichunk in range(nchunks):

            #### Prepare klist file ###########################################
            Prepare_klist_chunk(args.CASE, chunks, ichunk)

            #### Call x_dmft.py lapw1 #########################################
            RunXdmft("lapw1", up, dn, so, c)

            #### Join split energy files for dmftgk ############################
            JoinEnergyFiles(args.CASE, so, up, dn)

            #### Call x_dmft.py dmft1 #########################################
            RunXdmft("dmft1", up, dn, so, c)

            #### Call x_dmft.py dmftu #########################################
            RunXdmft("dmftu -g", up, dn, so, c)

            #### Call x_dmft.py dmftgk ########################################
            # we first need preparation of input file
            PrepareDmftgkInput(args.CASE, ncix, args.gamma, args.gammac, up,
                                dn, so)

            # then we can call dmftgk
            RunXdmft("dmftgk", up, dn, so, c)

            #### Read Green's function ########################
            for icix in ucix:
                gf[icix][ichunk*args.chunks:ichunk*args.chunks+chunks[ichunk],:,:,:] = \
                    ReadGf("G_k"+str(icix+1), sigind[icix])
    
    # if no dmft run, we have to check that G_k files exists
    # and then we read them directly
    else:
        for icix in ucix:
            if not os.path.isfile("G_k"+str(icix+1)):
                raise FileNotFoundError("ERROR: G_k"+str(icix+1)
                                        + " file does not exist."
                                        + " You should not use --no-dmft"
                                        + " option without running the"
                                        + " necessary x_dmft.py steps first.")

        if args.nambu:
            # --nambu --no-dmft path:
            # G_k files were produced by dmftgk_nambu.py with
            # cixdim = 2*norb.  ReadGf returns (nkpt, nfreq, 2*norb, 2*norb).
            # We extract normal [:norb, :norb] and anomalous [:norb, norb:].
            gf = {}
            gf_anom = {}
            for icix in ucix:
                norb = len(sigind[icix])
                gf_full = ReadGf("G_k"+str(icix+1), sigind[icix], nambu=True)
                gf[icix]      = gf_full[:, :, :norb, :norb]
                gf_anom[icix] = gf_full[:, :, :norb, norb:]

        else:
            # Normal --no-dmft path — identical to ePerio.py.
            gf = {}
            for icix in ucix:
                gf[icix] = ReadGf("G_k"+str(icix+1), sigind[icix])
        
        # also read frequency mesh from sig.inp1
        wmesh = np.loadtxt("sig.inp", usecols=(0,), dtype=float)[:]
        # check that frequency mesh is consistent with other self-energy files
        for icix in ucix:
            if np.shape(wmesh)[0] != np.shape(gf[icix])[1]:
                raise ValueError("frequency mesh in self-energy"
                                    + " files do not match!")
        
    #### Periodization ###################################################
    Gper = {}
    Fper = {} if args.nambu else None

    for icix in ucix:
        # check if this correlated block is a cluster by checking
        # dimensions of Gper compared to cluster size
        if np.shape(gf[icix])[2] == len(inperio['cluster']):
            Gper[icix] = periodize(gf[icix], k_unfold, inperio['cluster'])
            if args.nambu:
                # Anomalous block has same cluster geometry — periodize
                # with the same cluster positions (Eq. 7 applies to F too)
                Fper[icix] = periodize(gf_anom[icix], k_unfold, inperio['cluster'])
        else:
            # if not just create zero array
            Gper[icix] = np.zeros(np.shape(gf[icix])[:2],dtype=complex)
            if args.nambu:
                Fper[icix] = np.zeros(np.shape(gf_anom[icix])[:2],dtype=complex)

    # and finally save
    save_Gper(ucix, nucix, Gper, wmesh, args.CASE, inperio['nkpt'], Fper=Fper)

    # clean up
    # vector files
    for f in glob.glob("./"+args.CASE+".vector_*"):
        os.remove(f)
    # Udmft* files
    for f in glob.glob("./Udmft.*"):
        os.remove(f)
    # G_* files
    for f in glob.glob("./G_*"):
        os.remove(f)
    # g_* files
    for f in glob.glob("./g_*"):
        os.remove(f)
