<p align="center">
  <img src="https://raw.githubusercontent.com/idptools/afrc/main/afrc_logo.png" width="420" alt="AFRC logo"/>
</p>

<h1 align="center">afrc</h1>

<p align="center">
  <strong>The Analytical Flory Random Coil &mdash; a sequence-specific reference model for unfolded and intrinsically disordered proteins.</strong>
</p>

<p align="center">
  <a href="https://pypi.org/project/afrc/"><img src="https://img.shields.io/pypi/v/afrc.svg" alt="PyPI version"/></a>
  <a href="https://pypi.org/project/afrc/"><img src="https://img.shields.io/pypi/pyversions/afrc.svg" alt="Supported Python versions"/></a>
  <a href="https://afrc.readthedocs.io/en/latest/?badge=latest"><img src="https://readthedocs.org/projects/afrc/badge/?version=latest" alt="Documentation status"/></a>
  <a href="https://www.gnu.org/licenses/lgpl-3.0"><img src="https://img.shields.io/badge/License-LGPL%20v3-blue.svg" alt="License: LGPL v3"/></a>
</p>

---

## Overview

`afrc` implements the **Analytical Flory Random Coil (AFRC)**: a closed-form, sequence-specific
polymer model that reports the dimensions a polypeptide would adopt if it behaved as an ideal
chain in a theta solvent (apparent scaling exponent of 0.5, no finite-size effects). Given only
an amino acid sequence, it returns a wide range of polymeric properties &mdash; instantly and
without running any simulations.

The AFRC is intended as a **reference (null) model**, not a predictor of real dimensions.
Real unfolded-state dimensions depend on sequence-encoded chain&ndash;chain and chain&ndash;solvent
interactions that the AFRC deliberately omits. Its value is as a fixed, sequence-matched
touchstone: deviations of a simulation or experiment *from* the AFRC directly report
sequence-specific intramolecular interactions, and normalising to the AFRC lets chains of
different length and composition be compared on a common footing.

Alongside the AFRC, the package ships a family of classic analytical polymer models that share
a common interface, so the same sequence can be compared against several reference frames.

## Installation

`afrc` requires **Python 3.10 or later** (tested through Python 3.14) and depends only on
NumPy and SciPy.

```bash
pip install afrc
```

## Quick start

```python
from afrc import AnalyticalFRC

# create an AnalyticalFRC object from a sequence
P = AnalyticalFRC('MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDI')

# ensemble-average dimensions (Angstroms)
mean_rg  = P.get_mean_radius_of_gyration()
mean_e2e = P.get_mean_end_to_end_distance()
mean_rh  = P.get_mean_hydrodynamic_radius()

# full probability distributions, returned as (distances, probabilities)
rg_r,  rg_p  = P.get_radius_of_gyration_distribution()
re_r,  re_p  = P.get_end_to_end_distribution()

# distribution of the distance between two specific residues
d_r, d_p = P.get_interresidue_distance_distribution(4, 20)

# whole-chain maps
distance_map = P.get_distance_map()
contact_map  = P.get_contact_map(15.0)   # contact fractions at a 15 A threshold
```

## Polymer models included

Every model takes an amino acid sequence and exposes a common interface
(`get_end_to_end_distribution`, `get_mean_end_to_end_distance`, ...), so they can be swapped in
and compared directly.

| Model | Import | Description |
| --- | --- | --- |
| Analytical Flory Random Coil | `from afrc import AnalyticalFRC` | Sequence-specific ideal (theta-state) chain; the reference null model. |
| Freely jointed chain | `from afrc.polymer_models.fjc import FreelyJointedChain` | Ideal chain with finite extensibility (non-Gaussian Kuhn&ndash;Grun). |
| Freely rotating chain | `from afrc.polymer_models.frc import FreelyRotatingChain` | Ideal chain with a tunable characteristic ratio (stiffness). |
| Worm-like chain (Zhou) | `from afrc.polymer_models.wlc import WormLikeChain` | Semiflexible chain parameterised by a persistence length. |
| Worm-like chain (O'Brien) | `from afrc.polymer_models.wlc2 import WormLikeChain2` | Semiflexible chain; better large-chain stability, also gives Rg. |
| Self-avoiding walk | `from afrc.polymer_models.saw import SAW` | Good-solvent (excluded-volume) chain at a fixed scaling exponent. |
| nu-dependent SAW | `from afrc.polymer_models.nudep_saw import NuDepSAW` | Excluded-volume chain with a tunable Flory scaling exponent. |

The full mathematical formalism, parameters, and usage examples for each model are in the
[documentation](https://afrc.readthedocs.io/).

## Google Colab

[Use the AFRC (and the other polymer models) directly in a Google Colab notebook.](https://colab.research.google.com/drive/1WHw8ous7IgcKd2LKYuJLeBTlkdEYoRAk?usp=sharing)

## Documentation

Full documentation &mdash; including a theory section (formalism, parameters, references) and an
application section (usage examples and code reference) for every model &mdash; is hosted at
[afrc.readthedocs.io](https://afrc.readthedocs.io/). Worked, plotted examples for each model are
in the [`demo/`](https://github.com/idptools/afrc/tree/main/demo) directory.

## Citation

If you use the AFRC in your work, please cite:

> Alston, J. J., Ginell, G. M., Soranno, A., & Holehouse, A. S. (2023). The Analytical Flory
> Random Coil is a simple-to-use reference model for unfolded and disordered proteins.
> *The Journal of Physical Chemistry B*, 127(21), 4746&ndash;4760.
> https://doi.org/10.1021/acs.jpcb.3c01619

## Help and contributing

If you find a bug or have a feature request, please open
[an issue on GitHub](https://github.com/idptools/afrc/issues).

## Authors

The AFRC was developed by Garrett Ginell, Jhullian Alston, and Alex Holehouse in the
[Holehouse Lab](https://www.holehouselab.com/).

## License

Distributed by the [Holehouse Lab](https://www.holehouselab.com/) under the GNU Lesser General
Public License (LGPL v3).

## References

- Alston, J. J., Ginell, G. M., Soranno, A., & Holehouse, A. S. (2023). The Analytical Flory
  Random Coil is a simple-to-use reference model for unfolded and disordered proteins.
  *J. Phys. Chem. B*, 127(21), 4746&ndash;4760.
- Flory, P. J. (1969). *Statistical Mechanics of Chain Molecules*. Wiley-Interscience.
- Mao, A. H., Lyle, N., & Pappu, R. V. (2013). Describing sequence&ndash;ensemble relationships
  for intrinsically disordered proteins. *Biochemical Journal*, 449(2), 307&ndash;318.
