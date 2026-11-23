"""Physical ensemble, magnetic channels and MIXER policy for native response.

The flags describe different choices. INM controls charge normalization;
recomputeEF selects the reference Fermi search. Ferro controls how the two
spin channels obtain their common chemical potential. After fermif each
Ferro charge call uses recomputeEF=0 even when the common mu was searched.
Screened stress records the derivative constraint separately in response_policy;
its default conserves the central particle number even for a supplied EF.dat.
"""
from dataclasses import asdict, dataclass
from pathlib import Path
import ast


@dataclass(frozen=True)
class ResponsePolicy:
    magnetic: bool
    ferro: bool
    fixed_number: bool
    normalization_enabled: bool

    @property
    def channels(self):
        return ("", "dn") if self.magnetic else ("",)

    def record(self):
        return asdict(self)


def control_fields(path):
    lines = Path(path).read_text().splitlines()
    for i, line in enumerate(lines):
        words = line.split('#', 1)[0].split()
        if words:
            if len(words) < 3 or words[0] not in ('c', 'e', 'm'):
                raise ValueError('invalid dmft2 charge/eigenvalue control: ' + str(path))
            return lines, i, words
    raise ValueError('empty dmft2 control: ' + str(path))


def literal_parameters(path):
    """Read the two simple policy assignments without executing params.dat.

    run_dmft may supply already evaluated values for computed parameters.
    Standalone callers must use literal Ferro/recomputeEF assignments, rather
    than accidentally executing imports or impurity setup a second time.
    """
    path = Path(path)
    if not path.is_file():
        return {}
    tree = ast.parse(path.read_text(), filename=str(path))
    out = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            names = [t.id for t in node.targets if isinstance(t, ast.Name)]
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names = [node.target.id]
        else:
            continue
        for name in set(names) & {'Ferro', 'recomputeEF'}:
            if node not in tree.body:
                raise ValueError('computed policy requires an explicit override: ' + name)
            try:
                out[name] = ast.literal_eval(node.value)
            except (ValueError, TypeError):
                raise ValueError('computed policy requires an explicit override: ' + name) from None
    return out


def detect_policy(directory, case, *, control=None, override=None, electronic_model="dmft", native_soc=False):
    """Detect the second DMFT branch; it does not imply a spin-polarized KS map."""
    directory = Path(directory)
    nonempty = lambda suffix: (directory/(case+suffix)).is_file() and (directory/(case+suffix)).stat().st_size > 0
    if (nonempty('.inso') and not native_soc) or any(nonempty(s) for s in ('.clmup', '.clmdn', '.incup', '.incdn')):
        raise NotImplementedError('response requires no SOC and one unpolarized KS/core potential')
    if electronic_model not in ('dmft', 'dft'):
        raise ValueError('unknown electronic response model')
    magnetic = bool(nonempty('.indmfldn')) if electronic_model == 'dmft' else False
    if electronic_model == 'dmft':
        _, _, words = control_fields(control or directory/(case+'.indmf2'))
        mode = int(words[2])
        if mode not in (0, 1):
            raise NotImplementedError('response ensemble for gap-search recomputeEF=2 is not defined')
    else:
        if nonempty('.indmfldn'):
            raise ValueError('DFT response cannot use DMFT magnetic channels')
        mode = 1
    override = {} if override is None else dict(override)
    if set(override) - {'magnetic', 'ferro', 'fixed_number', 'normalization_enabled'}:
        raise ValueError('unknown physical response policy')
    if any(not isinstance(value, bool) for value in override.values()):
        raise ValueError('physical policy flags must be booleans')
    params = {} if electronic_model == 'dft' or {'ferro', 'fixed_number'} <= set(override) else literal_parameters(directory/'params.dat')
    ferro = override.get('ferro', params.get('Ferro', False))
    if not isinstance(ferro, bool):
        raise ValueError('Ferro must be True or False')
    # Without Ferro, the live prepared primary input is authoritative. With
    # Ferro, params carries the physical policy hidden by temporary m=0 calls.
    physical_mode = params.get('recomputeEF', mode) if ferro else mode
    if 'fixed_number' not in override and (type(physical_mode) is not int or physical_mode not in (0, 1)):
        raise ValueError('physical recomputeEF must be 0 or 1 for response')
    fixed = override.get('fixed_number', bool(physical_mode))
    if not isinstance(fixed, bool):
        raise ValueError('fixed_number must be a boolean')
    if ferro and not magnetic:
        raise ValueError('Ferro=True requires the second DMFT channel')
    if electronic_model == 'dft' and (ferro or not fixed):
        raise ValueError('Mermin DFT response requires a nonmagnetic fixed-N ensemble')
    fields = (directory/(case+'.inm')).read_text().splitlines()[0].split()
    if len(fields) < 3 or fields[2].upper() not in ('YES', 'NO'):
        raise ValueError('explicit INM YES/NO normalization required')
    if float(fields[1]) != 0:
        raise NotImplementedError('charged-background INM policy requires a separate response')
    normalize = fields[2].upper() == 'YES'
    if 'magnetic' in override and override['magnetic'] != magnetic:
        raise ValueError('recorded spin policy disagrees with case.indmfldn')
    if 'normalization_enabled' in override and override['normalization_enabled'] != normalize:
        raise ValueError('recorded charge normalization disagrees with case.inm')
    return ResponsePolicy(magnetic, ferro, fixed, normalize)
