"""
cli.py

Command-line tools for afrc.

``afrc-ensemble`` generates a 3D conformational ensemble (one bead per residue)
for a sequence from any of the package's polymer models, writes it as a PDB/XTC
pair (or a single multi-model PDB), and reports how well the ensemble reproduces
the model's statistics.

Copyright Alex Holehouse 2018-2026 (holehouselab.com).

"""

from __future__ import annotations

import argparse
import os
import sys
import textwrap
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import numpy as np

from .afrc import AnalyticalFRC
from .ensemble import save_conformations, validate_sequence
from .exceptions import AFRCException
from .polymer_models.fjc import FJCException, FreelyJointedChain
from .polymer_models.frc import FRCException, FreelyRotatingChain
from .polymer_models.nudep_saw import NuDepSAW, NuDepSAWException
from .polymer_models.saw import SAW, SAWException
from .polymer_models.wlc import WLCException, WormLikeChain
from .polymer_models.wlc2 import WLC2Exception, WormLikeChain2

# every exception a model can raise for bad input; these become a clean error
_MODEL_ERRORS: tuple[type[Exception], ...] = (AFRCException, FJCException, FRCException, WLCException,
                                               WLC2Exception, SAWException, NuDepSAWException)

# the PDB format numbers models with four digits
_MAX_PDB_MODELS: int = 9999


# .....................................................................................
#
@dataclass(frozen=True)
class EnsembleModel:
    """
    How ``afrc-ensemble`` builds and describes one polymer model.

    Attributes
    ----------
    name : str
        Display name.

    parameters : tuple of str
        The command-line parameter options (as argparse destinations, e.g.
        ``'lp'``) this model accepts.

    build : callable
        Function ``(sequence, args) -> (model, per_call_kwargs, description)``.
        ``per_call_kwargs`` are passed to the model's ``sample_conformations()``
        and ``check_ensemble()``; ``description`` summarizes the parameters.

    """

    name: str
    parameters: tuple[str, ...]
    build: Callable[[str, argparse.Namespace], tuple[Any, dict[str, float], str]]


def _value(args: argparse.Namespace, name: str, default: float) -> float:
    """
    Return a parameter from the command line, or its default if it was not given.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments.

    name : str
        The parameter's argparse destination.

    default : float
        The model's default value.

    Returns
    -------
    float
        The value to use.

    """
    value = getattr(args, name)
    return float(default if value is None else value)


def _build_afrc(seq: str, args: argparse.Namespace) -> tuple[Any, dict[str, float], str]:
    """Build the AFRC (no parameters)."""
    return AnalyticalFRC(seq), {}, 'AFRC'


def _build_fjc(seq: str, args: argparse.Namespace) -> tuple[Any, dict[str, float], str]:
    """Build a freely jointed chain."""
    b = _value(args, 'segment_length', 3.8)
    return FreelyJointedChain(seq, b=b), {}, f'freely jointed chain (b = {b:g} A)'


def _build_frc(seq: str, args: argparse.Namespace) -> tuple[Any, dict[str, float], str]:
    """Build a freely rotating chain."""
    b, c_inf = _value(args, 'segment_length', 3.8), _value(args, 'c_inf', 2.0)
    return FreelyRotatingChain(seq, b=b, c_inf=c_inf), {}, f'freely rotating chain (b = {b:g} A, c_inf = {c_inf:g})'


def _build_wlc(seq: str, args: argparse.Namespace) -> tuple[Any, dict[str, float], str]:
    """Build a worm-like chain, reported against the Zhou (2004) form."""
    b, lp = _value(args, 'segment_length', 3.8), _value(args, 'lp', 3.0)
    return WormLikeChain(seq, lp=lp, aa_size=b), {}, f'worm-like chain (lp = {lp:g} A, aa_size = {b:g} A; Zhou)'


def _build_wlc2(seq: str, args: argparse.Namespace) -> tuple[Any, dict[str, float], str]:
    """Build a worm-like chain, reported against the O'Brien (2009) form."""
    b, lp = _value(args, 'segment_length', 3.8), _value(args, 'lp', 3.0)
    return WormLikeChain2(seq, lp=lp, aa_size=b), {}, f"worm-like chain (lp = {lp:g} A, aa_size = {b:g} A; O'Brien)"


def _build_saw(seq: str, args: argparse.Namespace) -> tuple[Any, dict[str, float], str]:
    """Build a self-avoiding walk (Gaussian approximation)."""
    prefactor = _value(args, 'prefactor', 5.5)
    return SAW(seq), {'prefactor': prefactor}, f'self-avoiding walk (prefactor = {prefactor:g} A, nu = 0.598; Gaussian approximation)'


def _build_saw_nu(seq: str, args: argparse.Namespace) -> tuple[Any, dict[str, float], str]:
    """Build a nu-dependent self-avoiding walk (Gaussian approximation)."""
    nu, prefactor = _value(args, 'nu', 0.5), _value(args, 'prefactor', 5.5)
    return NuDepSAW(seq), {'nu': nu, 'prefactor': prefactor}, f'nu-dependent SAW (nu = {nu:g}, prefactor = {prefactor:g} A; Gaussian approximation)'


# the polymer models afrc-ensemble can draw from
ENSEMBLE_MODELS: dict[str, EnsembleModel] = {
    'afrc': EnsembleModel('AFRC', (), _build_afrc),
    'fjc': EnsembleModel('freely jointed chain', ('segment_length',), _build_fjc),
    'frc': EnsembleModel('freely rotating chain', ('segment_length', 'c_inf'), _build_frc),
    'wlc': EnsembleModel('worm-like chain (Zhou)', ('segment_length', 'lp'), _build_wlc),
    'wlc2': EnsembleModel("worm-like chain (O'Brien)", ('segment_length', 'lp'), _build_wlc2),
    'saw': EnsembleModel('self-avoiding walk', ('prefactor',), _build_saw),
    'saw-nu': EnsembleModel('nu-dependent self-avoiding walk', ('nu', 'prefactor'), _build_saw_nu),
}

# every model parameter option, mapped to its flag for error messages
_PARAMETER_FLAGS: dict[str, str] = {'segment_length': '--segment-length', 'c_inf': '--c-inf', 'lp': '--lp',
                                    'prefactor': '--prefactor', 'nu': '--nu'}

_EXAMPLES = """\
models (-m):
  afrc     Analytical Flory Random Coil (exact)
  fjc      freely jointed chain (exact)                   --segment-length
  frc      freely rotating chain (exact)                  --segment-length, --c-inf
  wlc      worm-like chain, compared with Zhou's P(r)      --segment-length, --lp
  wlc2     worm-like chain, compared with O'Brien's P(r)   --segment-length, --lp
  saw      self-avoiding walk (Gaussian approximation)    --prefactor
  saw-nu   nu-dependent SAW (Gaussian approximation)      --nu, --prefactor

examples:
  afrc-ensemble -s MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDI -n 5000 -o p53
      writes p53.pdb, p53.xtc and p53_report.txt from the AFRC

  afrc-ensemble -s MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDI -m wlc --lp 4 -o p53_wlc
      a worm-like chain with a 4 A persistence length

  afrc-ensemble -s MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDI -m saw-nu --nu 0.55 -n 500 --pdb-only --seed 1
      a single multi-model PDB (no mdtraj needed) from the nu-dependent SAW
"""


# .....................................................................................
#
def _positive_int(value: str) -> int:
    """
    argparse type for a strictly positive integer.

    Parameters
    ----------
    value : str
        The command-line value.

    Returns
    -------
    int
        The value as an int.

    Raises
    ------
    argparse.ArgumentTypeError
        If the value is not an integer greater than zero.

    """
    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f'must be a positive integer (got {value!r})')
    if number < 1:
        raise argparse.ArgumentTypeError(f'must be a positive integer (got {value!r})')
    return number


# .....................................................................................
#
def _positive_float(value: str) -> float:
    """
    argparse type for a strictly positive number.

    Parameters
    ----------
    value : str
        The command-line value.

    Returns
    -------
    float
        The value as a float.

    Raises
    ------
    argparse.ArgumentTypeError
        If the value is not a number greater than zero.

    """
    try:
        number = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f'must be a positive number (got {value!r})')
    if not number > 0:
        raise argparse.ArgumentTypeError(f'must be a positive number (got {value!r})')
    return number


# .....................................................................................
#
def _exponent(value: str) -> float:
    """
    argparse type for a scaling exponent strictly between 0 and 1.

    Parameters
    ----------
    value : str
        The command-line value.

    Returns
    -------
    float
        The value as a float.

    Raises
    ------
    argparse.ArgumentTypeError
        If the value is not a number strictly between 0 and 1.

    """
    try:
        number = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f'must be a number between 0 and 1 (got {value!r})')
    if not 0 < number < 1:
        raise argparse.ArgumentTypeError(f'must be a number between 0 and 1 (got {value!r})')
    return number


# .....................................................................................
#
def build_ensemble_parser() -> argparse.ArgumentParser:
    """
    Build the argument parser for ``afrc-ensemble``.

    Returns
    -------
    argparse.ArgumentParser
        The parser.

    """

    parser = argparse.ArgumentParser(
        prog='afrc-ensemble',
        description=('Generate a 3D conformational ensemble (one bead per residue) for a sequence from one of the '
                     'afrc polymer models, write it as a PDB/XTC pair, and report how well it reproduces the '
                     'model\'s statistics.'),
        epilog=_EXAMPLES,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-s', '--sequence', required=True,
                        help='amino acid sequence (one-letter codes, case insensitive)')
    parser.add_argument('-n', '--number-of-conformers', type=_positive_int, default=1000, metavar='N',
                        help='number of conformations to generate (default: 1000)')
    parser.add_argument('-m', '--model', type=str.lower, choices=list(ENSEMBLE_MODELS), default='afrc',
                        help='polymer model to draw the ensemble from (default: afrc; see the list below)')
    parser.add_argument('-o', '--out', default='out', metavar='PREFIX',
                        help='output prefix: writes PREFIX.pdb, PREFIX.xtc and PREFIX_report.txt (default: out)')
    parser.add_argument('--pdb-only', action='store_true',
                        help='write every conformation to a single multi-model PDB file instead of a PDB/XTC pair '
                             '(no mdtraj needed; at most 9999 conformations)')
    parser.add_argument('--seed', type=int, default=None,
                        help='random seed, for a reproducible ensemble (default: a fresh seed, recorded in the report)')

    parameters = parser.add_argument_group('model parameters (each defaults to the model\'s own default)')
    parameters.add_argument('--segment-length', type=_positive_float, default=None, metavar='A',
                            help='length per residue in Angstroms: the bond length b (fjc, frc) or aa_size (wlc, wlc2); default 3.8')
    parameters.add_argument('--c-inf', type=_positive_float, default=None, metavar='C',
                            help='characteristic ratio (frc); default 2.0')
    parameters.add_argument('--lp', type=_positive_float, default=None, metavar='A',
                            help='persistence length in Angstroms (wlc, wlc2); default 3.0')
    parameters.add_argument('--prefactor', type=_positive_float, default=None, metavar='A',
                            help='size prefactor in Angstroms (saw, saw-nu); default 5.5')
    parameters.add_argument('--nu', type=_exponent, default=None, metavar='NU',
                            help='Flory scaling exponent, between 0 and 1 (saw-nu); default 0.5')
    return parser


# .....................................................................................
#
def ensemble_main(argv: list[str] | None = None) -> int:
    """
    Entry point for the ``afrc-ensemble`` command.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (without the program name). Defaults to
        ``sys.argv[1:]``.

    Returns
    -------
    int
        Exit status: 0 on success and 1 if the ensemble could not be generated
        or written. (Invalid arguments exit with status 2, via argparse.)

    """

    parser = build_ensemble_parser()
    args = parser.parse_args(argv)

    # a parameter option that the chosen model does not use is a usage error
    model = ENSEMBLE_MODELS[args.model]
    for parameter, flag in _PARAMETER_FLAGS.items():
        if getattr(args, parameter) is not None and parameter not in model.parameters:
            users = ', '.join(key for key, entry in ENSEMBLE_MODELS.items() if parameter in entry.parameters)
            parser.error(f'{flag} does not apply to the {args.model} model (it is used by: {users})')

    try:
        return _run_ensemble(args, model)
    except _MODEL_ERRORS as error:
        print(f'afrc-ensemble: error: {error}', file=sys.stderr)
        return 1


# .....................................................................................
#
def _run_ensemble(args: argparse.Namespace, model_entry: EnsembleModel) -> int:
    """
    Generate, write and check an ensemble from parsed ``afrc-ensemble`` arguments.

    Everything that can be checked up front (the sequence, the model
    parameters, the output directory, mdtraj for XTC output, the PDB model
    limit) is checked before the ensemble is generated, so a bad run fails fast
    and writes nothing.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    model_entry : EnsembleModel
        The chosen model.

    Returns
    -------
    int
        Exit status (0).

    Raises
    ------
    AFRCException
        If the inputs are invalid or the ensemble cannot be written (or a
        model-specific exception for invalid model parameters).

    """

    from . import __version__

    # output files: tolerate a prefix given with a .pdb or .xtc extension
    prefix = str(args.out)
    root, extension = os.path.splitext(prefix)
    if extension.lower() in ('.pdb', '.xtc'):
        prefix = root
    directory = os.path.dirname(os.path.abspath(prefix))
    if not os.path.isdir(directory):
        raise AFRCException(f'Output directory {directory} does not exist')

    pdb_file, xtc_file, report_file = f'{prefix}.pdb', f'{prefix}.xtc', f'{prefix}_report.txt'
    n_conformations = int(args.number_of_conformers)

    if args.pdb_only:
        if n_conformations > _MAX_PDB_MODELS:
            raise AFRCException(f'A multi-model PDB holds at most {_MAX_PDB_MODELS} conformations (asked for {n_conformations}); drop --pdb-only to write an XTC file')
        written = [pdb_file]
    else:
        try:
            import mdtraj  # type: ignore[import-untyped]  # noqa: F401
        except ImportError:
            raise AFRCException('Writing an XTC file needs mdtraj; install it ("pip install mdtraj") or use --pdb-only')
        written = [pdb_file, xtc_file]

    # the sequence and the model (which validates its own parameters)
    sequence = validate_sequence(args.sequence)
    model, model_kwargs, description = model_entry.build(sequence, args)

    # always use an explicit seed, so any run can be reproduced from its report
    seed = int(args.seed) if args.seed is not None else int(np.random.default_rng().integers(2**31))

    conformations = model.sample_conformations(n=n_conformations, seed=seed, **model_kwargs)
    save_conformations(conformations, sequence, prefix, pdb_only=args.pdb_only,
                       remark=f'{description} ensemble from afrc {__version__} (seed {seed}): {n_conformations} conformations')

    report = model.check_ensemble(conformations, **model_kwargs)

    sequence_lines = textwrap.wrap(sequence, width=60)
    header = [f'afrc version  : {__version__}',
              f'Model         : {description}',
              f'Sequence      : {sequence_lines[0]}']
    header += [f'                {line}' for line in sequence_lines[1:]]
    header += [f'Seed          : {seed} (rerun with --seed {seed} to reproduce this ensemble)',
               f'Files         : {", ".join(written + [report_file])}']

    text = report.format(header)
    with open(report_file, 'w') as fh:
        fh.write(text)
    print(text, end='')

    return 0


if __name__ == '__main__':  # pragma: no cover
    sys.exit(ensemble_main())
