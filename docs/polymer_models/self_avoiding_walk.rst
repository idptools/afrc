Self-avoiding walk
=========================================================

The self-avoiding walk (SAW), exposed through :class:`~afrc.polymer_models.saw.SAW`, models a
chain in a good solvent, where excluded-volume interactions swell the chain relative to an
ideal coil. It uses the universal end-to-end distribution form of des Cloizeaux as
implemented by O'Brien et al. (2009), at the fixed good-solvent scaling exponent
:math:`\nu = 0.598`. It is composition-independent.

Mathematical formalism
---------------------------------------------------------

The end-to-end distribution is the scaling form

.. math::

   P(r) = \frac{a}{R_{ee}}
          \left( \frac{r}{R_{ee}} \right)^{2+\theta}
          \exp\!\left[ -b \left( \frac{r}{R_{ee}} \right)^{\delta} \right],

with the small-:math:`r` exponent :math:`\theta = 0.3`, the large-:math:`r` exponent
:math:`\delta = 2.5`, normalisation constants :math:`a = 3.679`, :math:`b = 1.232`, and a
chain-length-dependent scale

.. math::

   R_{ee} = \texttt{prefactor} \cdot N^{0.598}.

The radius of gyration is obtained from the end-to-end size via the universal ratio

.. math::

   \frac{\langle R_g^2 \rangle}{\langle R_e^2 \rangle}
   = \frac{\gamma(\gamma + 1)}{2(\gamma + 2\nu)(\gamma + 2\nu + 1)},
   \qquad \gamma = 1.1615,\ \nu = 0.598.

The exponent :math:`\nu` here is the same one that sets :math:`R_{ee}` above - both read
the object's ``nu`` attribute - so the :math:`R_g` and :math:`R_e` returned by the model
always describe the same chain.

Parameters
---------------------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Parameter
     - Default
     - Meaning and typical values
   * - ``prefactor``
     - 5.5 Å
     - Sets the absolute per-monomer length scale in :math:`R_{ee} = \texttt{prefactor}\,N^{0.598}`.
       The normalisation constants :math:`a` and :math:`b` are chosen so that
       :math:`R_{ee}` is exactly the root-mean-square end-to-end distance of the
       distribution. Values around 5-6 Å are reasonable, but the prefactor should be tuned
       to match explicit excluded-volume simulations for quantitative work.

The scaling exponent is fixed at the good-solvent value (:math:`\nu = 0.598`, held on the
object as ``nu``). To vary :math:`\nu` continuously, use the
:doc:`nu-dependent SAW <nu_dependent_saw>` instead - evaluated at
:math:`\nu = 0.598` with the same ``prefactor`` it reproduces this model.

**What to expect for a protein.** Because :math:`\nu \approx 0.6 > 0.5`, the SAW is *more
expanded* than the AFRC and is an appropriate reference for a strongly solvated, expanded IDR.
The prefactor sets where the absolute dimensions land; with a value near 5.5 Å the SAW gives
end-to-end and :math:`R_g` values noticeably larger than the theta-state AFRC.

Citations
---------------------------------------------------------

1. O'Brien, E. P., Morrison, G., Brooks, B. R., & Thirumalai, D. (2009). How accurate are
   polymer models in the analysis of Förster resonance energy transfer experiments on
   proteins? *The Journal of Chemical Physics*, 130(12), 124903.
2. des Cloizeaux, J. (1974). Lagrangian theory for a self-avoiding random chain.
   *Physical Review A*, 10(5), 1665-1669.
