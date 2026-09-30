Freely rotating chain
=========================================================

The freely rotating chain, exposed through
:class:`~afrc.polymer_models.frc.FreelyRotatingChain`, models the chain as :math:`N` bonds of
length :math:`b` joined at a fixed bond angle but with unrestricted (free) torsion angles. It
is an ideal chain - Gaussian end-to-end statistics with scaling exponent :math:`\nu = 0.5` -
whose absolute size is set by a single stiffness parameter, the characteristic ratio
:math:`C_\infty`. (The module is called ``frc``, but note this is not the Flory Random Coil
that the AFRC is built on.)

Mathematical formalism
---------------------------------------------------------

The mean-squared end-to-end distance uses the exact finite-:math:`N` freely-rotating-chain
result

.. math::

   \langle R^2 \rangle = C_\infty N b^2
        - 2 b^2\, \frac{\alpha\,(1 - \alpha^{N})}{(1 - \alpha)^2},
   \qquad \alpha = \frac{C_\infty - 1}{C_\infty + 1},

where :math:`\alpha` is the cosine of the angle between successive bonds and
:math:`C_\infty = (1+\alpha)/(1-\alpha)`. The first term is the long-chain limit and the
second is the finite-size correction (which vanishes when :math:`\alpha = 0`). The end-to-end
distribution is then the standard Gaussian chain form

.. math::

   P(r) = 4\pi r^2 \left( \frac{3}{2\pi \langle R^2 \rangle} \right)^{3/2}
          \exp\!\left( -\frac{3 r^2}{2 \langle R^2 \rangle} \right),

and the radius of gyration uses the ideal-chain relation
:math:`\sqrt{\langle R_g^2 \rangle} = \sqrt{\langle R^2 \rangle}/\sqrt{6}`.

Parameters
---------------------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Parameter
     - Default
     - Meaning and typical values
   * - ``b``
     - 3.8 Å
     - Bond (segment) length. The default is the Cα-Cα distance, i.e. one virtual bond per
       residue, so that the contour length is :math:`N b`.
   * - ``c_inf``
     - 2.0
     - Characteristic ratio :math:`C_\infty` - a dimensionless stiffness. ``c_inf = 1``
       (:math:`\alpha = 0`) recovers the freely jointed chain; ``c_inf = 2`` corresponds to a
       tetrahedral bond angle. Larger values give a stiffer, more extended ideal chain. Must
       be > 0.

.. note::

   A *genuine* freely rotating chain (free torsion) cannot reach the large characteristic
   ratio of a real polypeptide (:math:`C_\infty \approx 9`); that value arises from hindered
   rotation between backbone dihedrals. The :doc:`Analytical Flory Random Coil <afrc>`
   captures those local restrictions directly, so for a sequence-specific theta-state
   reference use the AFRC. The freely rotating chain is best thought of as a tunable,
   composition-independent ideal-chain reference.

**What to expect for a protein.** With one virtual bond per residue
(:math:`b = 3.8` Å), :math:`R_e \approx \sqrt{C_\infty}\, b\sqrt{N}` and
:math:`R_g \approx R_e/\sqrt{6}`. The default :math:`C_\infty = 2` gives dimensions in the
ballpark of the theta-state AFRC; raising :math:`C_\infty` swells the chain while preserving
ideal (:math:`\nu = 0.5`) scaling.

3D ensembles
---------------------------------------------------------

``sample_conformations()`` builds freely rotating chains directly: bonds of exactly :math:`b`, every bond angle fixed (at :math:`180^\circ - \arccos\alpha`, e.g. 109.47° for ``c_inf = 2``) and every torsion uniformly random. Its mean-squared distances are exactly the finite-:math:`N` result above for every pair of residues, but for short separations its distance distributions are not Gaussian, unlike the analytical :math:`P(r)`. Note that ``c_inf = 1`` gives a fixed 90° bond angle - the same mean-squared size as the freely jointed chain, but a different chain. ``check_ensemble()`` and ``afrc-ensemble -m frc`` (see :doc:`../cli`) check an ensemble against the exact mean-squared distances, the fixed bond length and angle, and finite extensibility.

Citations
---------------------------------------------------------

1. Flory, P. J. (1969). *Statistical Mechanics of Chain Molecules*. Wiley-Interscience.
2. Rubinstein, M., & Colby, R. H. (2003). *Polymer Physics*. Oxford University Press.
