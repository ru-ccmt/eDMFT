"""Shared chemical-potential sequencing for DMFT-only collinear magnetism.

Native calls keep one unpolarized KS potential. The two spin densities are
independent; only their total electron number is constrained. Temporary DEF
and control files preserve the user's prepared inputs, including live edits.
"""
from contextlib import contextmanager
from pathlib import Path
import math
import subprocess
import tempfile
from dmft2_step import entries, replace_filenames, execution_definition
from lapw_response_policy import control_fields


def chemical_potential(directory='.'):
    mu = float((Path(directory)/'EF.dat').read_text().split()[0])
    if not math.isfinite(mu):
        raise ValueError('finite common chemical potential required in EF.dat')
    return mu


@contextmanager
def measurement_definition(definition, *, mode, recompute, mu):
    """Change only the three measurement fields in a disposable control copy."""
    with execution_definition(definition) as native:
        text = Path(native).read_text()
        source = Path(entries(text)[2])
        lines, i, words = control_fields(source)
        words[:3] = [mode, format(mu, '.17g'), str(recompute)]
        lines[i] = ' '.join(words)
        with tempfile.TemporaryDirectory(prefix='.dmft2-spin-input-', dir='.') as tmp:
            folder = Path(tmp)
            control = folder/source.name
            control.write_text('\n'.join(lines)+'\n')
            df = folder/'dmft2.def'
            df.write_text(replace_filenames(text, {2: control}))
            yield str(df)


def run_channels(definitions, policy, run, fermi, *, export=False, search_only=False):
    """Invoke run(label, DEF) in the proper common-mu order.

    Ordinary compensated magnets adjust mu in the first charge call. Ferro
    states require BOTH eigenvalue files before fermif and both charge calls.
    export=True reuses the established mu in every channel and skips searching.
    That is a measurement choice, not a change of the response ensemble.
    """
    if len(definitions) != len(policy.channels):
        raise ValueError('missing DMFT spin definition')
    if not export and policy.magnetic and policy.ferro and policy.fixed_number:
        for channel, definition in zip(policy.channels, definitions):
            with measurement_definition(definition, mode='e', recompute=1,
                                        mu=chemical_potential()) as df:
                run('dmft2_eigen'+channel, df)
        fermi([entries(Path(d).read_text())[15] for d in definitions])
    if search_only and (not policy.fixed_number or policy.ferro):
        return
    for index, (channel, definition) in enumerate(zip(policy.channels, definitions)):
        search = int(not export and index == 0 and policy.fixed_number and not policy.ferro)
        with measurement_definition(definition, mode='c', recompute=search,
                                    mu=chemical_potential()) as df:
            run('dmft2'+channel, df)
        if search_only:
            return


def run_fermif(executable, bands, output):
    """Find one chemical potential from both spin spectra; fail before charge."""
    if len(bands) != 2 or any(not Path(p).is_file() for p in bands):
        raise ValueError('fermif requires both freshly computed band files')
    command = [str(executable), *map(str, bands)]
    print('FERMIF ..... running:', command, file=output, flush=True)
    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True) as child:
        for line in child.stdout:
            print(line, end='', file=output, flush=True)
        code = child.wait()
    if code:
        raise RuntimeError('fermif failed with status '+str(code))
    chemical_potential()


def average_charge(case, root, output):
    """Native branch outputs contain a factor two; combineud takes their mean.

    Keep the raw branches under the existing run_dmft backup names. This is
    called only after BOTH charge calls, never once per spin or before fermif.
    Core density and stress assembly remain shared downstream operations.
    """
    with subprocess.Popen([str(Path(root)/'combineud'), case], stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True) as child:
        for line in child.stdout:
            print(line, end='', file=output, flush=True)
        code = child.wait()
    average = Path(case+'.clmval_aver')
    if code or not average.is_file() or average.stat().st_size == 0:
        raise RuntimeError('combineud did not produce the averaged valence charge')
    Path(case+'.clmval').replace(case+'.clmval_up')
    Path(case+'.clmvaldn').replace(case+'.clmval_dn')
    average.replace(case+'.clmval')
