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

and the normalisation prefactors (Soranno 2020, Eq. 9b)

.. math::

   A_1 = \frac{\delta}{4\pi}\,
         \frac{\Gamma\!\left[\frac{5+g}{\delta}\right]^{(3+g)/2}}
              {\Gamma\!\left[\frac{3+g}{\delta}\right]^{(5+g)/2}},
   \qquad
   A_2 = \left(
            \frac{\Gamma\!\left[\frac{5+g}{\delta}\right]}
                 {\Gamma\!\left[\frac{3+g}{\delta}\right]}
         \right)^{\delta/2}.

:math:`A_1` normalises :math:`P(r)` to unity and :math:`A_2` is fixed by requiring that the
root-mean-square end-to-end distance of the distribution is exactly the size scale

.. math::

   R_{ee} = \texttt{prefactor} \cdot N^{\nu}.

.. note::

   Because :math:`A_2` is defined by that requirement,
   :math:`\sqrt{\langle R_e^2 \rangle} = \texttt{prefactor} \cdot N^{\nu}` holds exactly,
   for any :math:`\nu`. Setting :math:`\nu = 0.598` recovers the fixed-exponent
   :doc:`SAW <self_avoiding_walk>` to within the small differences between the two models'
   rounded exponents.

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
     - Flory scaling exponent. Must lie strictly between 0 and 1 (a ``NuDepSAWException``
       is raised otherwise). Physically meaningful values run from :math:`\approx 1/3`
       (collapsed globule, poor solvent), through :math:`0.5` (ideal / theta solvent), to
       :math:`\approx 0.588` (good solvent, fully expanded). Sweeping :math:`\nu` lets you
       place a measured chain on the collapse-to-expansion axis.
   * - ``prefactor``
     - 5.5 Å
     - Sets the absolute per-monomer length scale (as for the SAW): it *is* the
       root-mean-square end-to-end distance of a one-residue chain, and
       :math:`\sqrt{\langle R_e^2 \rangle} = \texttt{prefactor}\,N^{\nu}`. Around 5-6 Å is
       typical.

**What to expect for a protein.** At :math:`\nu = 0.5` the model gives theta-state
(ideal-chain) scaling; increasing :math:`\nu` toward 0.588 swells the chain to good-solvent
dimensions, while decreasing toward 1/3 collapses it. Most aqueous IDRs sit somewhere between
:math:`\nu \approx 0.5` and :math:`0.6`.

Note that the absolute dimensions are set by ``prefactor``, not by :math:`\nu` alone. With the
default 5.5 Å a :math:`\nu = 0.5` chain comes out roughly 12% more compact than the
:doc:`AFRC <afrc>`, because the AFRC's own composition-weighted prefactor
(:math:`R_0^{\mathrm{rms}} \approx 6.3` Å for a typical sequence) is larger. To use this model
as a :math:`\nu`-tunable version of the AFRC, set ``prefactor`` to the AFRC's
:math:`R_0^{\mathrm{rms}}` for your sequence; the two then agree closely at :math:`\nu = 0.5`.

Citations
---------------------------------------------------------

1. Zheng, W., Zerze, G. H., Borgia, A., Mittal, J., Schuler, B., & Best, R. B. (2018).
   Inferring properties of disordered chains from FRET transfer efficiencies.
   *The Journal of Chemical Physics*, 148(12), 123329.
2. Soranno, A. (2020). Physical basis of the disorder-order transition. *Archives of
   Biochemistry and Biophysics*, 685, 108305.
3. Le Guillou, J. C., & Zinn-Justin, J. (1977). Critical exponents for the n-vector model in
   three dimensions from field theory. *Physical Review Letters*, 39(2), 95-98.
