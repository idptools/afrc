"""
ensemble_report.py

Checks that a generated 3D ensemble reproduces the statistics of the chain model
it was drawn from.

Every model we generate ensembles for has exactly known mean-squared
inter-residue distances :math:`D_{ij} = \\langle r_{ij}^2 \\rangle`, and these
fix the root-mean-square radius of gyration and end-to-end distance too. For a
Gaussian chain (the AFRC) every distance also follows a Maxwell distribution set
by :math:`D_{ij}`, so the mean distances, the Kirkwood-Riseman hydrodynamic
radius and the full distance distributions can be checked as well. Models with
rigid geometry (fixed bond lengths or angles, or a maximum extension) add exact
structural checks.

``compare_ensemble_to_model()`` measures each of these in the ensemble, with a
standard error estimated from the ensemble itself, and compares it with the
model's exact expectation, described by a ``ModelExpectations``. A correct
ensemble agrees with every check to within sampling error.

Copyright Alex Holehouse 2018-2026 (holehouselab.com).

"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray
from scipy import stats
from scipy.spatial.distance import pdist

from .exceptions import AFRCException

# a check fails if the ensemble is more than this many standard errors from the
# model. The per-conformation quantities behind the checks (Rg^2, 1/r and so on)
# are skewed, so their means have somewhat heavier tails than a normal
# distribution. Calibrated on ~2000 simulated correct AFRC ensembles (2-154
# residues, 100-5000 conformations), the largest |z| seen was 4.0. For a
# 50-residue chain a 2% error in the size of the chain is always flagged with
# 1000 conformations, and a 1% error with 5000
Z_THRESHOLD: float = 5.0

# the distribution-shape (Kolmogorov-Smirnov) check fails below this p-value.
# Three pairs are tested, so a correct ensemble fails with probability ~3e-5
KS_P_THRESHOLD: float = 1e-5

# below this many conformations the standard errors are too rough to judge by
MIN_CONFORMATIONS: int = 100

# upper edges of the sequence-separation bands used by the pair-based checks
SEPARATION_BAND_EDGES: tuple[int, ...] = (1, 3, 10, 30, 100, 300, 1000)

# cap on the number of pairwise distances evaluated for the pair-based checks
# (conformations x residue pairs), so very long chains stay fast
MAX_PAIR_EVALUATIONS: int = 100_000_000

# relative tolerance for the exact (non-statistical) structural checks
STRUCTURAL_TOLERANCE: float = 1e-9


# .....................................................................................
#
@dataclass(frozen=True)
class Check:
    """
    One comparison between the ensemble and the model.

    Attributes
    ----------
    name : str
        What is being compared.

    passed : bool or None
        True if the ensemble agrees with the model, False if it does not, and
        None for a context-only comparison (or one that could not be assessed).

    summary : str
        Human-readable result (may span several lines).

    values : dict of str to float
        The numbers behind the summary, for programmatic use.

    structural : bool
        True for an exact geometric check (e.g. a fixed bond length), which
        does not depend on sampling and so is assessed for any number of
        conformations.

    """

    name: str
    passed: bool | None
    summary: str
    values: dict[str, float] = field(default_factory=dict)
    structural: bool = False


# .....................................................................................
#
@dataclass(frozen=True)
class ModelExpectations:
    """
    What a chain model says an ensemble drawn from it should look like.

    Attributes
    ----------
    model_name : str
        Name of the model, for display.

    mean_squared_distances : np.ndarray
        The model's exact [N x N] matrix of mean-squared inter-residue distances
        (Angstroms squared).

    gaussian : bool
        True if every inter-residue distance is exactly Gaussian (Maxwell
        distributed), which enables the mean-distance, hydrodynamic-radius and
        distribution checks. Default is False.

    bond_length : float, optional
        Exact length of every bond between consecutive beads, if the model fixes
        it (Angstroms).

    bond_angle_cosine : float, optional
        Exact cosine of the angle between consecutive bond vectors, if the model
        fixes it.

    contour_length_per_residue : float, optional
        The contour length per residue, if the model is finitely extensible: no
        two beads k residues apart can be further apart than k times this.

    reference_rg : float, optional
        A radius of gyration the model reports itself, shown as context.

    reference_rg_label : str
        What ``reference_rg`` is, for display.

    reference_re : float, optional
        A whole-chain end-to-end distance the model reports itself, shown as
        context.

    reference_re_label : str
        What ``reference_re`` is, for display.

    end_to_end_distribution : tuple of np.ndarray, optional
        The model's analytical end-to-end distribution ``(distances,
        probabilities)`` for the first-to-last pair (N - 1 residues apart), shown
        as context.

    end_to_end_label : str
        What that distribution is, for display.

    notes : tuple of str
        Model-specific caveats printed with the report.

    """

    model_name: str
    mean_squared_distances: NDArray[np.float64]
    gaussian: bool = False
    bond_length: float | None = None
    bond_angle_cosine: float | None = None
    contour_length_per_residue: float | None = None
    reference_rg: float | None = None
    reference_rg_label: str = ''
    reference_re: float | None = None
    reference_re_label: str = ''
    end_to_end_distribution: tuple[NDArray[np.float64], NDArray[np.float64]] | None = None
    end_to_end_label: str = ''
    notes: tuple[str, ...] = ()


# .....................................................................................
#
@dataclass(frozen=True)
class EnsembleReport:
    """
    The result of comparing an ensemble with its model.

    Attributes
    ----------
    model_name : str
        Name of the model the ensemble was drawn from.

    n_residues : int
        Number of residues (beads).

    n_conformations : int
        Number of conformations in the ensemble.

    n_pair_conformations : int
        Number of conformations used for the pair-based checks (all of them,
        unless that would exceed ``MAX_PAIR_EVALUATIONS`` distances).

    checks : list of Check
        Comparisons the ensemble must pass to be consistent with the model.

    context : list of Check
        Comparisons that are informative but not expected to match exactly.

    notes : list of str
        Model-specific caveats.

    """

    model_name: str
    n_residues: int
    n_conformations: int
    n_pair_conformations: int
    checks: list[Check]
    context: list[Check]
    notes: list[str] = field(default_factory=list)

    @property
    def assessed(self) -> bool:
        """Whether the statistical checks could be assessed (enough conformations and residues)."""
        return any(check.passed is not None for check in self.checks if not check.structural)

    @property
    def failed(self) -> list[str]:
        """Names of the checks that failed."""
        return [check.name for check in self.checks if check.passed is False]

    @property
    def passed(self) -> bool:
        """Whether the statistical checks were assessed and every check passed."""
        return self.assessed and not self.failed

    def format(self, header: list[str] | None = None) -> str:
        """
        Format the report as plain text.

        Parameters
        ----------
        header : list of str, optional
            Extra lines (e.g. the sequence and output files) shown under the title.

        Returns
        -------
        str
            The report, ready to print or save.

        """

        title = f'{self.model_name} ensemble report'
        lines = [title, '=' * len(title)]
        lines.extend(header or [])
        lines.append(f'Residues      : {self.n_residues}')
        lines.append(f'Conformations : {self.n_conformations}')
        if self.n_pair_conformations < self.n_conformations:
            lines.append(f'                (pair-based checks use the first {self.n_pair_conformations}, to keep the run fast)')
        lines.append('')

        lines.append(f'Checks against the model - a model-like ensemble agrees with each to within sampling error ({Z_THRESHOLD:g} standard errors)')
        lines.append('-' * 100)
        for check in self.checks:
            lines.extend(_format_check(check))

        if self.context:
            lines.append('')
            lines.append('Context - informative comparisons that are not expected to match exactly')
            lines.append('-' * 100)
            for check in self.context:
                lines.extend(_format_check(check))

        if self.notes:
            lines.append('')
            lines.append('Notes')
            lines.append('-' * 100)
            for note in self.notes:
                lines.append(f'* {note}')

        lines.append('')
        if self.failed:
            lines.append(f'Verdict: WARNING - {len(self.failed)} check(s) outside the expected range: {"; ".join(self.failed)}.')
        elif not self.assessed:
            lines.append(f'Verdict: NOT ASSESSED - use at least {MIN_CONFORMATIONS} conformations (and two residues) to check the ensemble.')
        else:
            lines.append('Verdict: PASS - the ensemble reproduces the model statistics checked above to within sampling error.')

        return '\n'.join(lines) + '\n'


# .....................................................................................
#
def _format_check(check: Check) -> list[str]:
    """
    Format one check: a status tag and name, then its (possibly multi-line) summary.

    Parameters
    ----------
    check : Check
        The check to format.

    Returns
    -------
    list of str
        The formatted lines.

    """
    tag = {True: '[OK]', False: '[FAIL]', None: '[info]'}[check.passed]
    return [f'{tag:<8}{check.name}'] + [f'{"":<8}{line}' for line in check.summary.split('\n')]


# .....................................................................................
#
def _separation_bands(n_res: int) -> list[tuple[int, int]]:
    """
    Split sequence separations 1 to n_res - 1 into roughly logarithmic bands.

    Parameters
    ----------
    n_res : int
        Number of residues.

    Returns
    -------
    list of tuple of int
        ``(lowest, highest)`` separation in each band, covering every separation
        once.

    """
    bands = []
    lower = 1
    for upper in SEPARATION_BAND_EDGES + (n_res - 1,):
        upper = min(upper, n_res - 1)
        if upper >= lower:
            bands.append((lower, upper))
            lower = upper + 1
    return bands


# .....................................................................................
#
def _z_scores(mean: NDArray[np.float64], expected: NDArray[np.float64], se: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    Compute z-scores, allowing for quantities that do not fluctuate at all.

    Some quantities are fixed by a model's geometry - the distance between
    neighbouring beads of a freely jointed or freely rotating chain, say - so
    across conformations they differ only by floating-point rounding. Their
    standard error is then ~1e-17 rather than exactly zero, and dividing a
    rounding-level difference by it gives a meaningless, huge z-score. So the
    standard error is floored at ``STRUCTURAL_TOLERANCE`` times the expected
    value: genuine sampling errors are always far larger than this, while a
    real mismatch in a fixed quantity (anything above ~1e-9 relative) still
    gives a very large z-score.

    Parameters
    ----------
    mean : np.ndarray
        Ensemble means.

    expected : np.ndarray
        Model expectations.

    se : np.ndarray
        Standard errors of the means (NaN if they could not be estimated).

    Returns
    -------
    np.ndarray
        The z-scores (NaN where the standard error could not be estimated).

    """
    floor = STRUCTURAL_TOLERANCE * np.abs(expected)
    effective_se = np.where(np.isnan(se), np.nan, np.maximum(se, floor))
    with np.errstate(divide='ignore', invalid='ignore'):
        z: NDArray[np.float64] = (mean - expected) / effective_se
    return z


# .....................................................................................
#
def _banded_ratio_check(name: str, per_conformation: NDArray[np.float64], bands: list[tuple[int, int]], what: str) -> Check:
    """
    Check that the ensemble-to-model ratio of a pair quantity is 1 in every band.

    Each column of ``per_conformation`` is, for one band of sequence
    separations, the average over that band's residue pairs of (ensemble value /
    model expectation) in a single conformation. Its mean over conformations
    should be 1, and because each conformation contributes one number per band,
    the standard error from the spread across conformations accounts for every
    correlation between pairs.

    Parameters
    ----------
    name : str
        What is being compared.

    per_conformation : np.ndarray
        Array of shape [n_conformations x n_bands].

    bands : list of tuple of int
        The separation range of each band, for display.

    what : str
        Short description of the ratio, for display.

    Returns
    -------
    Check
        The comparison; it passes if every band is within ``Z_THRESHOLD``
        standard errors of 1.

    """

    n = per_conformation.shape[0]
    mean = np.mean(per_conformation, axis=0)
    se = np.std(per_conformation, axis=0, ddof=1) / np.sqrt(n) if n > 1 else np.full(len(bands), np.nan)
    z = _z_scores(mean, np.ones_like(mean), se)

    lines = []
    values: dict[str, float] = {}
    for (low, high), m, s, zz in zip(bands, mean, se, z):
        label = f'|i - j| = {low}' if low == high else f'|i - j| = {low}-{high}'
        lines.append(f'{label:<20} {what} {m:.4f} ± {s:.4f} (z = {zz:+.2f})')
        values[f'ratio_{low}_{high}'] = float(m)
        values[f'z_{low}_{high}'] = float(zz)

    passed = bool(np.all(np.abs(z) < Z_THRESHOLD)) if n >= MIN_CONFORMATIONS else None
    return Check(name, passed, '\n'.join(lines), values)


# .....................................................................................
#
def _scalar_check(name: str, samples: NDArray[np.float64], model: float, unit: str = 'A', transform: str = 'none') -> Check:
    """
    Compare the mean of per-conformation values with the model's expectation.

    Parameters
    ----------
    name : str
        What is being compared.

    samples : np.ndarray
        One value per conformation.

    model : float
        The model's expected mean of those values.

    unit : str
        Unit for display.

    transform : str
        ``'none'`` to report the mean itself, or ``'sqrt'`` to report the square
        root of the mean (e.g. an RMS value from mean-squared samples). The
        z-score is always computed on the mean of ``samples``.

    Returns
    -------
    Check
        The comparison.

    """

    n = len(samples)
    mean = float(np.mean(samples))
    se = float(np.std(samples, ddof=1) / np.sqrt(n)) if n > 1 else float('nan')
    z = float(_z_scores(np.array([mean]), np.array([model]), np.array([se]))[0])

    if transform == 'sqrt':
        shown, shown_model, shown_se = np.sqrt(mean), np.sqrt(model), se / (2 * np.sqrt(mean))
    else:
        shown, shown_model, shown_se = mean, model, se

    passed = bool(abs(z) < Z_THRESHOLD) if n >= MIN_CONFORMATIONS and not np.isnan(z) else None
    summary = (f'ensemble {shown:.3f} ± {shown_se:.3f} {unit} | model {shown_model:.3f} {unit} | '
               f'difference {100 * (shown / shown_model - 1):+.2f}% (z = {z:+.2f})')
    return Check(name, passed, summary, {'ensemble': float(shown), 'ensemble_se': float(shown_se), 'model': float(shown_model), 'z': z})


# .....................................................................................
#
def _structural_check(name: str, worst: float, description: str) -> Check:
    """
    An exact geometric check, which passes if the worst deviation is negligible.

    Parameters
    ----------
    name : str
        What is being checked.

    worst : float
        The largest relative deviation from the model's geometry.

    description : str
        What the model requires, for display.

    Returns
    -------
    Check
        The check (always assessed, whatever the number of conformations).

    """
    passed = bool(worst <= STRUCTURAL_TOLERANCE)
    return Check(name, passed, f'{description}; largest relative deviation {worst:.1e}', {'max_relative_deviation': worst}, structural=True)


# .....................................................................................
#
def compare_ensemble_to_model(conformations: NDArray[np.float64],
                              expectations: ModelExpectations,
                              max_pair_evaluations: int = MAX_PAIR_EVALUATIONS) -> EnsembleReport:
    """
    Check that an ensemble reproduces the statistics of the model it came from.

    For every model the ensemble is checked against the exact mean-squared
    inter-residue distances :math:`D_{ij}`:

    * the root-mean-square radius of gyration,
      :math:`\\langle R_g^2 \\rangle = \\frac{1}{N^2}\\sum_{i<j} D_{ij}`;
    * the first-to-last bead distance;
    * the mean-squared distance of every residue pair, in bands of sequence
      separation.

    For a Gaussian model it is also checked against the Maxwell distribution
    each distance follows: the mean distances (again in bands), the
    Kirkwood-Riseman hydrodynamic radius
    (:math:`\\langle 1/r_{ij} \\rangle = \\sqrt{6/(\\pi D_{ij})}`), and the full
    distribution of three representative distances (Kolmogorov-Smirnov test).
    Models with fixed bond lengths, fixed bond angles or a maximum extension get
    exact structural checks too.

    Parameters
    ----------
    conformations : np.ndarray
        Array of shape [n_conformations x N x 3] with the bead coordinates, in
        Angstroms.

    expectations : ModelExpectations
        What the model says the ensemble should look like.

    max_pair_evaluations : int
        Upper limit on conformations x residue pairs for the pair-based checks.
        If the ensemble exceeds it, those checks use the first conformations only.

    Returns
    -------
    EnsembleReport
        The checks, context comparisons and notes.

    Raises
    ------
    AFRCException
        If the conformations and distance matrix do not have matching shapes.

    """

    xyz = np.asarray(conformations, dtype=float)
    D = np.asarray(expectations.mean_squared_distances, dtype=float)
    if xyz.ndim != 3 or xyz.shape[2] != 3 or xyz.shape[0] == 0:
        raise AFRCException(f'Conformations must have shape [n_conformations x N x 3] (got {xyz.shape})')
    n_frames, n_res, _ = xyz.shape
    if D.shape != (n_res, n_res):
        raise AFRCException(f'The mean-squared distance matrix must be {n_res} x {n_res} to match the conformations (got {D.shape})')

    name = expectations.model_name
    checks: list[Check] = []
    context: list[Check] = []
    notes = list(expectations.notes)

    # a single bead has no internal structure to check
    if n_res < 2:
        return EnsembleReport(name, n_res, n_frames, n_frames, checks, context, notes)

    gaussian = expectations.gaussian

    # ------------------------------------------------------------------
    # per-conformation quantities, from every conformation
    centred = xyz - xyz.mean(axis=1, keepdims=True)
    rg_squared = np.mean(np.sum(centred * centred, axis=-1), axis=1)
    first_to_last = np.linalg.norm(xyz[:, -1] - xyz[:, 0], axis=-1)

    upper_i, upper_j = np.triu_indices(n_res, 1)
    D_pairs = D[upper_i, upper_j]

    checks.append(_scalar_check('Root-mean-square radius of gyration', rg_squared,
                                float(np.sum(D_pairs) / n_res**2), transform='sqrt'))
    if gaussian:
        checks.append(_scalar_check('Mean first-to-last bead distance', first_to_last,
                                    float(np.sqrt(8 * D[0, -1] / (3 * np.pi)))))
    else:
        checks.append(_scalar_check('Root-mean-square first-to-last bead distance', first_to_last**2,
                                    float(D[0, -1]), transform='sqrt'))

    # ------------------------------------------------------------------
    # pair-based quantities, from as many conformations as the budget allows
    n_pairs = len(D_pairs)
    n_pair_frames = int(min(n_frames, max(1, max_pair_evaluations // n_pairs)))

    separation = upper_j - upper_i
    pairs_at_separation = np.bincount(separation, minlength=n_res)[1:]

    # each pair's band of sequence separation, for the banded checks
    bands = _separation_bands(n_res)
    band_of_pair = np.zeros(n_pairs, dtype=int)
    for band, (low, high) in enumerate(bands):
        band_of_pair[(separation >= low) & (separation <= high)] = band
    pairs_in_band = np.bincount(band_of_pair, minlength=len(bands))

    model_mean_distance = np.sqrt(8.0 * D_pairs / (3.0 * np.pi))
    contour = expectations.contour_length_per_residue
    max_pair_length = separation * contour if contour is not None else None

    # with a fixed bond length, neighbouring beads sit exactly at their contour
    # length (and are covered by the bond-length check), so the extension worth
    # reporting is that of pairs two or more residues apart
    extension_pairs = np.ones(n_pairs, dtype=bool)
    if expectations.bond_length is not None and n_res > 2:
        extension_pairs = separation >= 2

    # accumulate, conformation by conformation, only what the checks need - one
    # number per band (rather than per pair) keeps memory small for long chains
    mean_inverse = np.empty(n_pair_frames)
    squared_ratio = np.empty((n_pair_frames, len(bands)))
    distance_ratio = np.empty((n_pair_frames, len(bands)))
    by_separation = np.empty((n_pair_frames, n_res - 1))
    worst_extension = 0.0
    for frame in range(n_pair_frames):
        # pdist returns the pairs in the same (row-major, upper-triangle) order as triu_indices
        distances = pdist(xyz[frame])
        squared_ratio[frame] = np.bincount(band_of_pair, weights=distances * distances / D_pairs, minlength=len(bands)) / pairs_in_band
        by_separation[frame] = np.bincount(separation, weights=distances, minlength=n_res)[1:] / pairs_at_separation
        if gaussian:
            with np.errstate(divide='ignore'):
                mean_inverse[frame] = np.mean(1.0 / distances)
            distance_ratio[frame] = np.bincount(band_of_pair, weights=distances / model_mean_distance, minlength=len(bands)) / pairs_in_band
        if max_pair_length is not None:
            worst_extension = max(worst_extension, float(np.max(distances[extension_pairs] / max_pair_length[extension_pairs])))

    if gaussian:
        # Kirkwood-Riseman: the ensemble Rh is 1/<1/r>, so check <1/r> (with its
        # standard error) and report it as Rh
        model_mean_inverse = float(np.mean(np.sqrt(6.0 / (np.pi * D_pairs))))
        inverse_check = _scalar_check('Kirkwood-Riseman hydrodynamic radius', mean_inverse, model_mean_inverse)
        ens, se, z = inverse_check.values['ensemble'], inverse_check.values['ensemble_se'], inverse_check.values['z']
        rh, rh_model, rh_se = 1.0 / ens, 1.0 / model_mean_inverse, se / ens**2
        checks.append(Check(inverse_check.name, inverse_check.passed,
                            f'ensemble {rh:.3f} ± {rh_se:.3f} A | model {rh_model:.3f} A | '
                            f'difference {100 * (rh / rh_model - 1):+.2f}% (z = {-z:+.2f})',
                            {'ensemble': rh, 'ensemble_se': rh_se, 'model': rh_model, 'z': -z}))

    checks.append(_banded_ratio_check(f'Mean-squared inter-residue distances (all {n_pairs} pairs, by sequence separation)',
                                      squared_ratio, bands, '<r^2> ensemble/model'))

    if gaussian:
        checks.append(_banded_ratio_check('Mean inter-residue distances (internal scaling, by sequence separation)',
                                          distance_ratio, bands, '<r> ensemble/model  '))

        # distribution shape: a local pair, a middle pair and the chain ends
        ks_pairs = sorted({(n_res // 2 - 1, n_res // 2) if n_res > 2 else (0, 1), (n_res // 4, (3 * n_res) // 4), (0, n_res - 1)})
        ks_pairs = [(i, j) for (i, j) in ks_pairs if 0 <= i < j < n_res]
        p_values = []
        for i, j in ks_pairs:
            distances = np.linalg.norm(xyz[:, j] - xyz[:, i], axis=-1)
            p_values.append(float(stats.kstest(distances, stats.maxwell(scale=np.sqrt(D[i, j] / 3.0)).cdf).pvalue))
        passed = bool(min(p_values) > KS_P_THRESHOLD) if n_frames >= MIN_CONFORMATIONS else None
        checks.append(Check('Inter-residue distance distributions (Kolmogorov-Smirnov vs the exact distribution)', passed,
                            'pairs ' + ', '.join(f'({i}, {j})' for i, j in ks_pairs) + ': p = ' + ', '.join(f'{p:.3f}' for p in p_values)
                            + f' (fails below {KS_P_THRESHOLD:g})',
                            {f'p_{i}_{j}': p for (i, j), p in zip(ks_pairs, p_values)}))

    # ------------------------------------------------------------------
    # exact structural checks
    bonds = np.diff(xyz, axis=1)
    bond_lengths = np.linalg.norm(bonds, axis=-1)
    if expectations.bond_length is not None:
        b = expectations.bond_length
        checks.append(_structural_check('Bond length', float(np.max(np.abs(bond_lengths / b - 1))),
                                        f'every bond is {b:.3f} A'))
    if expectations.bond_angle_cosine is not None and n_res >= 3:
        a = expectations.bond_angle_cosine
        units = bonds / bond_lengths[..., None]
        cosines = np.sum(units[:, 1:] * units[:, :-1], axis=-1)
        checks.append(_structural_check('Bond angle', float(np.max(np.abs(cosines - a))),
                                        f'every bond angle is {180.0 - np.degrees(np.arccos(a)):.2f} degrees '
                                        f'(consecutive bond vectors have cosine {a:+.4f}; deviation measured in that cosine)'))
    if contour is not None:
        worst_extension = max(worst_extension, float(np.max(first_to_last / ((n_res - 1) * contour))))
        which = 'pair (2+ residues apart)' if expectations.bond_length is not None and n_res > 2 else 'pair'
        checks.append(_structural_check('Finite extensibility', max(0.0, worst_extension - 1.0),
                                        f'no two beads are further apart than their contour length ({contour:.3f} A per residue); '
                                        f'the most extended {which} reaches {100 * worst_extension:.2f}% of it'))

    # ------------------------------------------------------------------
    # context
    rg = np.sqrt(rg_squared)
    if expectations.reference_rg is not None:
        context.append(Check(f'Radius of gyration vs the {name} {expectations.reference_rg_label}', None,
                             f'ensemble RMS {np.sqrt(np.mean(rg_squared)):.3f} A, mean {np.mean(rg):.3f} A | '
                             f'{name} {expectations.reference_rg:.3f} A',
                             {'ensemble_rms': float(np.sqrt(np.mean(rg_squared))), 'ensemble_mean': float(np.mean(rg)),
                              'model': float(expectations.reference_rg)}))
    if expectations.reference_re is not None:
        context.append(Check(f'First-to-last distance vs the {name} {expectations.reference_re_label}', None,
                             f'ensemble RMS {np.sqrt(np.mean(first_to_last**2)):.3f} A, mean {np.mean(first_to_last):.3f} A | '
                             f'{name} {expectations.reference_re:.3f} A (the whole-chain value treats the chain as N residues, '
                             f'the first-to-last distance as N-1)',
                             {'ensemble_rms': float(np.sqrt(np.mean(first_to_last**2))), 'ensemble_mean': float(np.mean(first_to_last)),
                              'model': float(expectations.reference_re)}))
    if expectations.end_to_end_distribution is not None:
        r, p = (np.asarray(a, dtype=float) for a in expectations.end_to_end_distribution)
        model_cdf = np.cumsum(p)
        ensemble_cdf = np.searchsorted(np.sort(first_to_last), r, side='right') / n_frames
        ks_distance = float(np.max(np.abs(ensemble_cdf - model_cdf)))
        context.append(Check(f'First-to-last distance distribution vs {expectations.end_to_end_label}', None,
                             f'mean {np.mean(first_to_last):.3f} A vs {np.sum(r * p):.3f} A | '
                             f'RMS {np.sqrt(np.mean(first_to_last**2)):.3f} A vs {np.sqrt(np.sum(r * r * p)):.3f} A | '
                             f'largest difference between the cumulative distributions {ks_distance:.3f}',
                             {'ensemble_mean': float(np.mean(first_to_last)), 'model_mean': float(np.sum(r * p)),
                              'cdf_difference': ks_distance}))

    context.append(Check('Adjacent bead spacing', None,
                         f'mean {np.mean(bond_lengths):.2f} A, 5-95% range {np.percentile(bond_lengths, 5):.2f}-{np.percentile(bond_lengths, 95):.2f} A',
                         {'mean': float(np.mean(bond_lengths))}))

    if n_res >= 4:
        k = np.arange(1, n_res)
        nu_ensemble = float(np.polyfit(np.log(k), np.log(np.mean(by_separation, axis=0)), 1)[0])
        rms_by_separation = np.sqrt(np.bincount(separation, weights=D_pairs, minlength=n_res)[1:] / pairs_at_separation)
        nu_model = float(np.polyfit(np.log(k), np.log(rms_by_separation), 1)[0])
        context.append(Check('Apparent scaling exponent (log-log fit of the internal scaling profile)', None,
                             f'ensemble {nu_ensemble:.3f} (from mean distances) | model {nu_model:.3f} (from RMS distances)',
                             {'ensemble': nu_ensemble, 'model': nu_model}))

    return EnsembleReport(name, n_res, n_frames, n_pair_frames, checks, context, notes)


# .....................................................................................
#
def compare_ensemble_to_gaussian_model(conformations: NDArray[np.float64],
                                       mean_squared_distances: NDArray[np.float64],
                                       model_name: str = 'AFRC',
                                       reference_mean_rg: float | None = None,
                                       reference_mean_re: float | None = None,
                                       max_pair_evaluations: int = MAX_PAIR_EVALUATIONS) -> EnsembleReport:
    """
    Check that an ensemble reproduces the statistics of a Gaussian chain model.

    A convenience wrapper around ``compare_ensemble_to_model()`` for a model in
    which every inter-residue distance is Gaussian, fully specified by its
    mean-squared inter-residue distances (as for the AFRC).

    Parameters
    ----------
    conformations : np.ndarray
        Array of shape [n_conformations x N x 3] with the bead coordinates, in
        Angstroms.

    mean_squared_distances : np.ndarray
        The model's symmetric [N x N] matrix of mean-squared inter-residue
        distances, in Angstroms squared.

    model_name : str
        Name of the model, for display. Default is 'AFRC'.

    reference_mean_rg : float, optional
        The model's own mean radius of gyration, if it defines one separately
        (the AFRC's is calibrated independently of its distances). Reported as
        context.

    reference_mean_re : float, optional
        The model's whole-chain mean end-to-end distance (the AFRC's uses N
        rather than N-1 residues). Reported as context.

    max_pair_evaluations : int
        Upper limit on conformations x residue pairs for the pair-based checks.

    Returns
    -------
    EnsembleReport
        The checks and context comparisons.

    Raises
    ------
    AFRCException
        If the conformations and distance matrix do not have matching shapes.

    """
    expectations = ModelExpectations(model_name, np.asarray(mean_squared_distances, dtype=float), gaussian=True,
                                     reference_rg=reference_mean_rg,
                                     reference_rg_label='<Rg> (calibrated separately from its distances; expect up to ~2.5% difference)',
                                     reference_re=reference_mean_re, reference_re_label='whole-chain <Re>')
    return compare_ensemble_to_model(conformations, expectations, max_pair_evaluations)
