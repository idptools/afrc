## Changelog

* **0.4.3** (unreleased):
  * **Behavior-changing scientific fixes:**
    * **Zhou worm-like chain (`WormLikeChain`):** the end-to-end distribution assigned real probability to distances *beyond the contour length*, which is physically impossible. The Zhou (2004) expression is a series expansion in `Lp/Lc` and `r/Lc` and is simply not meaningful for `r > Lc`, but for short and/or stiff chains it returns substantial positive values there - around 30% of the weight for a 2-residue chain at the default `lp = 3 A`, and over half for a 5-residue chain at `lp = 10 A`. That region is now explicitly zeroed. Chains of more than ~20 residues at the default persistence length are unaffected (the change is below 1e-6 relative), and the distribution reproduces the exact worm-like chain `<R^2> = 2*Lp*Lc - 2*Lp^2*(1 - exp(-Lc/Lp))` to machine precision for long chains; short chains now return smaller (and physical) values. A chain whose contour length is comparable to or shorter than the persistence length has no valid region at all and previously produced a NaN-filled distribution; it now raises a `WLCException`.
    * **Zhou worm-like chain grid:** the `r` grid was fixed at `21*sqrt(N)`, which is comfortably wide for `lp = 3 A` but cut into the tail of `P(r)` for stiffer chains (it reached only 2.4 times the size scale at `lp = 10 A`). The grid now always extends to at least four times the ideal-chain size scale `sqrt(2*Lp*Lc)`.
    * **Kirkwood-Riseman hydrodynamic radius:** `get_mean_hydrodynamic_radius('kirkwood-riseman')` computed `1/<r_ij>` - the inverse of the *mean* inter-residue distance, taken from the distance map - whereas the Kirkwood-Riseman equation requires the mean of the *inverse* distance, `<1/r_ij>`. For a Gaussian chain the two differ by exactly `4/pi`, so the returned Rh was systematically 27% too large (and `Rg/Rh` came out near 1.05 rather than the ~1.3-1.5 expected for an ideal chain). Each inter-residue `PolymerObject` now exposes `get_mean_inverse_end_to_end_distance()`, which returns the exact closed form `<1/r> = sqrt(6/(pi <r^2>))` for the Gaussian distribution, and the hydrodynamic radius is `1/<1/r_ij>` averaged over every residue pair. This is the same form of the equation used by Nygaard et al. (2017), Pesce et al. (2023) and SOURSOP's `mode='kr'`. For the 154-residue test sequence Rh drops from 29.35 A to 23.00 A. The `'nygaard'` mode is unchanged. A single-residue chain now raises an `AFRCException` in this mode rather than returning a divide-by-zero.

  * **Bug fixes (code):**
    * `get_contact_fraction()` raised a bare `IndexError` when the threshold was smaller than one grid step (e.g. `threshold = 0.03` at the default 0.05 A resolution), because it tried to take a grid spacing from a single-element slice. It now returns the (zero) weight in the `r = 0` bin. It also validates its residue indices before the `R1 == R2` short-circuit, so an invalid pair such as `(-1, -1)` raises an `AFRCException` instead of returning 1.0.
    * `get_pre_profile()` now routes `label_position` through the shared residue-index validator, so a non-integer label (e.g. a string) raises an `AFRCException` rather than a `TypeError`.
    * Replaced the per-element Python loops in the two worm-like chain distribution calculations with vectorized NumPy (the O'Brien model additionally zeroes the undefined `r >= Lc` region explicitly rather than relying on NaN detection).

  * **Documentation:**
    * The Zhou WLC theory page and class docstring now state the range of validity of the expansion (contour length much larger than the persistence length), that no probability is assigned beyond the contour length, and that too-short chains raise.
    * `get_internal_scaling()` described its return value as a `[2 x n]` matrix; it is `[n-1 x 2]` (one row per sequence separation), and the docstring now says so.
    * The AFRC theory page and the `get_mean_hydrodynamic_radius()` docstring now give the Kirkwood-Riseman expression actually evaluated, and note the `4/pi` distinction between `<1/r>` and `1/<r>`.

  * **Tests:**
    * Added regression tests for each of the above: no WLC probability beyond the contour length, agreement with the exact WLC second moment for long chains (which also pins the sign convention of the correction series), tail coverage for stiff chains, the too-short-chain exception, the O'Brien model's `r >= Lc` region, sub-resolution contact thresholds, contact-fraction index validation, non-integer PRE label positions, the internal-scaling return shape, and the Kirkwood-Riseman fix (agreement with an explicit `<1/r_ij>` sum, the `4/pi` ratio to the old estimator, a theta-like `Rg/Rh`, and the closed-form `<1/r>` against numerical integration of `P(r)/r`).


* **0.4.2** (August 2026):
  * **Bug fixes (code):**
    * Made the version-resolution logic in `afrc/__init__.py` robust: it previously keyed off the `READTHEDOCS` environment variable and a path relative to the *working directory* (`../afrc/_version.py`), so it could mis-resolve depending on where Python was launched from. It now tries three sources in turn: the generated `_version.py` (authoritative for an installed package, though in a git checkout it can lag the repository until the next build), the installed-package metadata via `importlib.metadata`, and finally `versioningit` anchored to the package's own location for raw source checkouts. A corrupt or truncated `_version.py` (e.g. from an interrupted build) no longer makes the package unimportable, and if every source fails the import fails loudly rather than silently reporting a placeholder version.
    * Removed `RE_RG_CORRECTION_FACTORS` from `config.py`. This constant encoded the ~0.949 correction for the old `<Re>/sqrt(6)` Rg scaling law and became dead code once that mode switched to the calibrated `RG_R0` prefactor; nothing in the package or tests referenced it.

  * **Documentation:**
    * Cleaned up remaining typos and stale prose across the docstrings and docs (a garbled paragraph in the docs overview, a docstring that incorrectly described `get_mean_end_to_end_distance()` as possibly returning an RMS value, a `Returns` block that referenced a nonexistent `mode` argument, corrected parameter ordering in the `validate_keyword()` docstring, and assorted spelling fixes).


* **0.4.1** (July 2026):
  * **Behavior-changing scientific fixes:**
    * **O'Brien worm-like chain (`WormLikeChain2`):** fixed a sign error in the closed-form radius of gyration. The `1/(4*C2^2)` term (which is `-Lp^2` in the Benoit-Doty form) was being added rather than subtracted, so `get_mean_radius_of_gyration()` over-estimated Rg - by ~5% for a 50-residue chain, and badly for short chains, where it left a spurious constant offset instead of the correct rigid-rod limit of `Lc^2/12`. The corrected expression now reproduces Benoit-Doty exactly and respects the rigid-rod bound (`Rg^2 <= Lc^2/12`) at every chain length.
    * **nu-dependent SAW (`NuDepSAW`):** fixed the `A1` and `A2` normalization prefactors. The gamma functions were being evaluated at `5 + g/delta` and `3 + g/delta` instead of `(5+g)/delta` and `(3+g)/delta`, and in `A1` the `(3+g)/2` and `(5+g)/2` terms were applied as multipliers rather than exponents. Because `A2` sets the width of the distribution, this meant the root-mean-square end-to-end distance did not match the requested size scale - by 33% at `nu = 0.4`, 12% at `nu = 0.5` and -5% at `nu = 0.598`. The unexplained factor of `pi` in the size scale (flagged in a source comment as a suspected symptom of an error elsewhere) was compensating for this and has been removed. With the corrected prefactors `sqrt(<Re^2>) = prefactor * N^nu` holds exactly for any `nu`, and setting `nu = 0.598` now reproduces the fixed-exponent `SAW` model. The corrected constants agree with those hard-coded in `saw.py`, which were derived independently.
    * **Mean Rg via the scaling law:** `get_mean_radius_of_gyration('scaling law')` returned `<Re>/sqrt(6)`. That ideal-chain relation holds between *root-mean-square* radii, so applying it to the mean end-to-end distance under-estimated `<Rg>` by ~5% and left the `'scaling law'` and `'distribution'` modes disagreeing. It now uses the calibrated composition-weighted `RG_R0` prefactor (which was already in `config.py` but unused, and is back-calculated from the same Lhuillier fits that generate the distribution). The two modes now agree to better than 0.05%. This also affects `get_mean_interresidue_radius_of_gyration()`, which defaults to `'scaling law'`.
    * **Self-avoiding walk (`SAW`) radius of gyration:** the universal `<Rg^2>/<Re^2>` ratio was evaluated with a hard-coded `nu = 0.589` while the model's end-to-end distance scales as `N^0.598`, so the Rg and Re returned by the same object described chains with different scaling exponents. Both now read a single `self.nu = 0.598` attribute (alongside `self.gamma`), set once in the constructor, so the two cannot drift apart again. `get_mean_radius_of_gyration()` returns values 0.65% smaller than before (`Rg/sqrt(<Re^2>)` goes from 0.4008 to 0.3982), and now agrees with `NuDepSAW` evaluated at `nu = 0.598` to seven significant figures. Note that `theta = 0.3` and `delta = 2.5` are deliberately left as they are: together with `a = 3.67853` and `b = 1.23152` they form the matched set of constants tabulated by O'Brien et al., and are only a rounded description of the same exponent.
    * **SAW / SAW-nu distribution grid:** the `r` grid was fixed at `21*sqrt(N)` while these models' size scale grows as `N^nu` with `nu > 0.5`, so for long chains the grid cut into the tail of `P(r)` and biased the mean and RMS low (by 0.5% at N = 1000 and 1.8% at N = 3000). The grid now always extends to at least four times the model's own size scale.
    * **`WormLikeChain2` short-chain check:** the constructor compared the number of residues against the persistence length, mixing a residue count with a length in Angstroms. It now compares the contour length (`N * aa_size`) with the persistence length, as the error message always claimed.

  * **Bug fixes (code):**
    * `sample_inter_residue_distance_distribution()` did not validate its residue indices. Out-of-range indices raised a bare `IndexError`, and - worse - negative indices silently wrapped around and sampled the wrong pair of residues. Both now raise an `AFRCException`. The redundant index-swapping logic was removed (the inter-residue matrix is symmetric).
    * Zero-length `PolymerObject`s crashed rather than behaving as the code's own comments intended: `get_end_to_end_distribution()` raised a divide-by-zero and `get_radius_of_gyration_distribution()` raised a `ValueError` from `np.vectorize`. Both now return a degenerate distribution with all weight at zero, as `FreelyJointedChain` already did. The same handling was added to `WormLikeChain`, `SAW` and `NuDepSAW`.
    * `PolymerObject.get_mean_end_to_end_distance()` and `get_mean_radius_of_gyration()` silently returned `None` when handed an unrecognized `calculation_mode`; they now raise an `AFRCException`.
    * `compute_apparent_rms_bond_length()` divided by zero for a single-residue chain (a chain with no bonds); it now raises an `AFRCException`.
    * `NuDepSAW` raised a bare `ZeroDivisionError` for `nu = 1` and produced nonsense for `nu <= 0`; `nu` is now validated to lie strictly between 0 and 1.
    * `__validate_residue_index()` caught only `ValueError`, so a non-castable, non-string residue index (e.g. `None`) escaped as a `TypeError` rather than an `AFRCException`.
    * `AFRCException` moved into its own `afrc/exceptions.py` module, removing the circular-import workaround in `iofunctions.py`. The existing `from afrc.afrc import AFRCException` import path is unchanged.
    * Removed the unused `other_models` dictionary from the `AnalyticalFRC` constructor, and replaced the remaining per-element Python loops in the SAW and SAW-nu distribution code with vectorized NumPy.

  * **Documentation:**
    * Updated the O'Brien WLC and SAW-nu theory pages for the corrected expressions (the latter loses the note rationalising the factor of `pi`), the SAW theory page for the corrected scaling exponent in the Rg ratio, and the AFRC theory page for the corrected Rg scaling law.
    * Fixed a duplicated bullet on the documentation front page (the Rg distribution was listed as the end-to-end distribution) and noted the freely jointed and freely rotating chain models there.
    * Documented the previously-undocumented `calculation_mode`, `sample_size`, `R1` and `R2` arguments, and fixed a docstring that broke Sphinx's `Returns` block.
    * Re-executed every notebook in `demo/` so the stored outputs reflect the corrected models, and updated `compare_end_to_end.ipynb`, which still called the pre-0.4.0 method names (`get_re_distribution`, `get_mean_re`, `get_mean_rg`) and could no longer be run.

  * **Tests:**
    * Added regression tests pinning each of the fixes above: Benoit-Doty and rigid-rod-bound checks for the O'Brien Rg, `sqrt(<Re^2>) = prefactor * N^nu` for both SAW models, SAW/SAW-nu agreement (in both Re and Rg) at `nu = 0.598`, the SAW's Rg ratio and size scaling sharing one exponent, long-chain grid coverage, scaling-law/distribution agreement for Rg, index validation, degenerate zero-length distributions, and invalid `calculation_mode`/`nu` handling. 105 tests, all passing.


* **0.4.0** (June 2026):
  * **New models:**
    * Added a `FreelyJointedChain` reference model (`afrc.polymer_models.fjc`) - a composition-independent FJC of `N` segments of length `b` (default 3.8 Å). It uses the non-Gaussian Kuhn-Grün end-to-end distribution (with a Cohen Padé inverse-Langevin), so it respects the finite extensibility of the chain (`r < Nb`) while reducing to the Gaussian AFRC at small extensions. It exposes the same `get_end_to_end_distribution()`, `get_mean_end_to_end_distance()`, `get_root_mean_squared_end_to_end_distance()`, `get_mean_radius_of_gyration()`, and `sample_end_to_end_distribution()` interface as the other auxiliary models.
    * Added a `FreelyRotatingChain` reference model (`afrc.polymer_models.frc`) - a composition-independent ideal chain with a tunable characteristic ratio `c_inf` (stiffness) and bond length `b`. It uses the exact finite-N freely-rotating-chain mean-squared end-to-end distance with a Gaussian distribution (`c_inf = 1` recovers the freely jointed chain; `c_inf = 2` is the tetrahedral value). It exposes the same interface as the other auxiliary models.

  * **Bug fixes (code):**
    * Fixed `PolymerObject.compute_apparent_rms_bond_length()`, which referenced a non-existent attribute (`RMS_Re`) and raised `AttributeError` on every call; it now uses `RMS_Re_scaling`.
    * Fixed `validate_keyword()` in `iofunctions.py`, which raised `NameError` (rather than the intended `AFRCException`) on invalid keyword input because `AFRCException` was never imported.
    * The auxiliary-model exception classes (`WLCException`, `SAWException`, `NuDepSAWException`) now subclass `Exception`, so raising them works instead of throwing `TypeError`.
    * `WormLikeChain2` (O'Brien WLC) no longer references an undefined `WLCException` name in its input-validation checks; it now raises its own `WLC2Exception`.
    * Fixed the lower-bound check in `AnalyticalFRC.get_pre_profile()`, which previously allowed `label_position == len(seq)` (an out-of-range index).
    * Replaced the deprecated `np.trapz` with `np.trapezoid` (with a fallback for older NumPy) in `get_contact_fraction()`.
    * Corrected `afrc/__init__.py`'s Read the Docs version branch, which checked a stale path (`../goose/_version.py`) copied from another package.
    * Removed dead/unused code (unused imports, dead variables and commented-out blocks) and silenced an expected divide-by-zero warning in `get_pre_profile()`.

  * **Behaviour-changing scientific fixes:**
    * **Zhou worm-like chain (`WormLikeChain`):** fixed an operator-precedence error in the `zeta` correction term (`5*Lp/4*Lc` was evaluated as `(5*Lp/4)*Lc` instead of `5*Lp/(4*Lc)`). The Zhou-WLC end-to-end distribution and derived means now reflect the corrected expression. Spurious negative probabilities in the far tail of the series expansion are now clamped to zero.
    * **SAW / SAW-ν models:** `get_mean_end_to_end_distance()` now returns the true mean (`Σ r·P(r)`) rather than the root-mean-square value, matching the convention used by the AFRC and WLC models. A dedicated `get_root_mean_squared_end_to_end_distance()` method was added, and `get_mean_radius_of_gyration()` now uses the RMS value internally (as its analytical ratio requires), so Rg values are unchanged.

  * **Documentation / docstrings:**
    * Corrected numerous inaccurate docstrings (e.g. distribution getters return a `(distances, probabilities)` tuple, not a 2D array), filled in empty docstrings, fixed a malformed parameter header and incorrect default-mode descriptions, and fixed many typos throughout the source.
    * Added API documentation for the auxiliary polymer models (`WormLikeChain`, `WormLikeChain2`, `SAW`, `NuDepSAW`, `FreelyJointedChain`) and refreshed the quickstart examples.
    * Added a new "Polymer Models (Theory)" documentation section with a subpage per model, each describing the mathematical formalism, the free parameters (with physical intuition and typical values for a protein), and the primary citations. The AFRC theory page was expanded with the model's origin (rejection-free rotational-isomeric-state simulations, absence of end effects) and behaviour, citing Alston et al. (2023).
    * Added a parallel "Polymer Models (Application)" documentation section with a subpage per model giving quick-start usage examples plus the full autogenerated code reference for each class. The previous combined API reference page (`modules/afrc`) was retired in favour of these per-model pages.
    * Added per-class demo notebooks to the `demo/` directory.

  * **Tests:**
    * Replaced the near-empty test module (which used a faulty `a - b < tol` assertion) with a systematic suite covering construction/validation, distribution sanity, analytic relationships (scaling exponent, Rg = Re/√6), distance/contact maps, hydrodynamic radius, PRE profiles, the auxiliary models, and golden-value regression snapshots.

  * **Packaging / project:**
    * Overhauled the README (badges, model overview table, installation/quick-start, citation).
    * Enriched `pyproject.toml` metadata (keywords, Trove classifiers, and `[project.urls]`).
    * Added the previously-declared-but-missing `afrc/py.typed` marker file.
    * Updated supported Python versions: `requires-python` is now `>=3.10` (dropping end-of-life 3.7-3.9) and support is declared and verified through Python 3.14. Bumped the Read the Docs build to Python 3.12 and added `scipy` to the docs requirements.
    * Removed defunct CI/service configuration: `.travis.yml`, `.lgtm.yml` (LGTM was retired by GitHub), and the Travis helper under `devtools/`.


* **0.3.7** (May 2025):
  * Moved to `pyproject.toml` and `versioningit` for build and versioning to enable compatiblity with later Python versions.
