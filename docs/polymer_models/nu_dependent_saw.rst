nu-dependent self-avoiding walk
=========================================================

The :math:`\nu`-dependent SAW, exposed through
:class:`~afrc.polymer_models.nudep_saw.NuDepSAW`, generalises the :doc:`self-avoiding walk
<self_avoiding_walk>` so that the Flory scaling exponent :math:`\nu` becomes a free parameter.
This lets a single model span the full range from a collapsed globule to a fully solvated coil.
It uses the form derived by Zheng et al. (2018) as written by Soranno (2020).

Mathematical formalism
---------------------------------------------------------

The end-to-end distribution is

.. math::

   P(r) = \frac{4\pi A_1}{R_{ee}}
          \left( \frac{r}{R_{ee}} \right)^{2+g}
          \exp\!\left[ -A_2 \left( \frac{r}{R_{ee}} \right)^{\delta} \right],

with the exponents

.. math::

   g = \frac{\gamma - 1}{\nu}, \qquad \delta = \frac{1}{1 - \nu},
   \qquad \gamma = 1.1615,

and normalisation prefactors :math:`A_1`, :math:`A_2` expressed through gamma functions of
:math:`g` and :math:`\delta` (Soranno 2020, Eq. 9b). The size scale is

.. math::

   R_{ee} = \texttt{prefactor} \cdot N^{\nu} \cdot \pi.

.. note::

   The factor of :math:`\pi` in :math:`R_{ee}` is an empirical correction (noted in the source)
   that brings the model into quantitative agreement with the :doc:`AFRC <afrc>` at
   :math:`\nu = 0.5` and with the :doc:`SAW <self_avoiding_walk>` at :math:`\nu \approx 0.598`.

The radius of gyration uses the same universal ratio as the SAW, evaluated at the chosen
:math:`\nu`.

Parameters
---------------------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Parameter
     - Default
     - Meaning and typical values
   * - ``nu``
     - 0.5
     - Flory scaling exponent. Physically meaningful values run from :math:`\approx 1/3`
       (collapsed globule, poor solvent), through :math:`0.5` (ideal / theta solvent), to
       :math:`\approx 0.588` (good solvent, fully expanded). Sweeping :math:`\nu` lets you
       place a measured chain on the collapse-to-expansion axis.
   * - ``prefactor``
     - 5.5 Å
     - Sets the absolute per-monomer length scale (as for the SAW). Around 5-6 Å is typical.

**What to expect for a protein.** At :math:`\nu = 0.5` the model reproduces theta-state
(AFRC-like) dimensions; increasing :math:`\nu` toward 0.588 swells the chain to good-solvent
dimensions, while decreasing toward 1/3 collapses it. Most aqueous IDRs sit somewhere between
:math:`\nu \approx 0.5` and :math:`0.6`.

Citations
---------------------------------------------------------

1. Zheng, W., Zerze, G. H., Borgia, A., Mittal, J., Schuler, B., & Best, R. B. (2018).
   Inferring properties of disordered chains from FRET transfer efficiencies.
   *The Journal of Chemical Physics*, 148(12), 123329.
2. Soranno, A. (2020). Physical basis of the disorder-order transition. *Archives of
   Biochemistry and Biophysics*, 685, 108305.
3. Le Guillou, J. C., & Zinn-Justin, J. (1977). Critical exponents for the n-vector model in
   three dimensions from field theory. *Physical Review Letters*, 39(2), 95-98.
