Analytical Flory Random Coil
=========================================================

The Analytical Flory Random Coil (AFRC) is the central model of the package, exposed through
the :class:`~afrc.AnalyticalFRC` object. It reproduces the dimensions of a polypeptide
behaving as an ideal (Gaussian) chain - one in which the apparent Flory scaling exponent is
:math:`\nu = 0.5`, analogous to a chain in a :math:`\theta`-solvent - and, unlike a finite
self-avoiding chain, it carries no finite-size ("dangling-end") effects. It is a
*parameter-free reference model*: it is fully determined by the amino acid sequence and has
nothing to fit. The model and its parameterisation are described in Alston et al. (2023).

Origin: numerical Flory Random Coil simulations
---------------------------------------------------------

The AFRC is a closed-form fit to numerical *Flory Random Coil* (FRC) ensembles. Those
ensembles are generated with all-atom Monte Carlo using Flory's rotational isomeric state
(RIS) approximation: at each step a residue is chosen at random and its backbone dihedrals
(:math:`\phi, \psi`) are reassigned to one of a precomputed set of residue-specific allowed
states, drawn from all-atom Ramachandran maps. The moves are **rejection-free** - only
sterically allowed local dihedrals are proposed and the resulting global conformation is
accepted unconditionally, with no through-space (chain-chain, chain-solvent, or chain-self)
interactions of any kind.

Two consequences follow, and together they are what make an analytical description possible:

* **The chain is ideal.** Because every monomer is "agnostic" to its surroundings, both
  global and internal dimensions scale with an apparent exponent of :math:`\nu = 0.5`, exactly
  as expected for a Gaussian chain in a :math:`\theta`-solvent.
* **There are no end effects.** Terminal residues sample the same conformational space as
  internal ones, so the "dangling-end" finite-size deviations seen in finite self-avoiding
  chains are absent. The internal scaling profiles for chains of every length superimpose,
  which means a single set of closed-form expressions can describe the chain over all length
  scales.

From simulations to a closed-form model
---------------------------------------------------------

Residue-specific prefactors are fit once, against homopolymer FRC simulations, and then used
analytically:

* The root-mean-square inter-residue distance follows
  :math:`\sqrt{\langle r_{ij}^2 \rangle} = A_0\,|i-j|^{\nu}`; fitting the FRC internal scaling
  profiles yields a per-residue prefactor :math:`A_0` (the :math:`R_0` / :math:`R_0^{\mathrm{rms}}`
  constants in this package).
* The radius-of-gyration prefactor :math:`X_0` is fit so that the analytical Lhuillier
  distribution matches the numerically generated :math:`P(R_g)`.

Because the RIS construction treats each residue independently, **heteropolymers** are handled
by taking a composition-weighted average of these homopolymer prefactors. This generalisation
was validated against FRC simulations of hundreds of heteropolymeric sequences (10-500
residues), reproducing both end-to-end and :math:`R_g` distributions with sub-angstrom
accuracy.

Mathematical formalism
---------------------------------------------------------

**End-to-end distance.** The end-to-end distribution is the standard Gaussian chain result

.. math::

   P(r) = 4\pi r^2 \left( \frac{3}{2\pi \langle R_e^2 \rangle} \right)^{3/2}
          \exp\!\left( -\frac{3 r^2}{2 \langle R_e^2 \rangle} \right),

where the root-mean-square size follows the ideal-chain scaling law

.. math::

   \sqrt{\langle R_e^2 \rangle} = R_0^{\mathrm{rms}}\, N^{1/2}.

The prefactor :math:`R_0^{\mathrm{rms}}` is the composition-weighted average of the per-residue
:math:`A_0` constants described above. An analogous prefactor :math:`R_0` gives the mean
end-to-end distance, :math:`\langle R_e \rangle = R_0\, N^{1/2}`.

**Radius of gyration.** The :math:`R_g` distribution uses the analytical fractal-polymer
form of Lhuillier (1988):

.. math::

   P(R_g) \propto N^{-\nu d}\, \frac{\rho}{N^{\nu}}
   \exp\!\left[ -\left(\frac{N^{\nu}}{\rho}\right)^{\alpha d}
               -\left(\frac{\rho}{N^{\nu}}\right)^{\delta} \right],

with :math:`\rho = X_0\, R_g`, dimensionality :math:`d = 3`, :math:`\nu = 1/2`,
:math:`\alpha = 1/(\nu d - 1) = 2`, and :math:`\delta = 1/(1-\nu) = 2`. The
composition-weighted prefactor :math:`X_0` again comes from the calibrated per-residue table.
The mean radius of gyration can be taken either as the expectation of this distribution or
from the ideal-chain relation :math:`R_g = \langle R_e \rangle / \sqrt{6}`.

**Hydrodynamic radius.** :math:`R_h` is available either from the Kirkwood-Riseman relation
applied to the full inter-residue distance map, or from the empirical :math:`R_g \to R_h`
conversion of Nygaard et al. (2017).

Because the model also exposes every inter-residue distance, it additionally provides distance
maps, contact-fraction maps, and per-residue PRE profiles for the same theta-state reference.

Behaviour and relationship to other models
---------------------------------------------------------

By construction the AFRC behaves like a :doc:`nu-dependent SAW <nu_dependent_saw>` evaluated
at :math:`\nu = 0.5`; the two distributions sit essentially on top of one another. Relative to
the other reference models in this package, the AFRC is *slightly more expanded* than the
:doc:`worm-like chain <worm_like_chain_zhou>` (at a persistence length of 3 Å) and
*substantially more compact* than the good-solvent :doc:`self-avoiding walk
<self_avoiding_walk>` (:math:`\nu \approx 0.588`). It therefore occupies the theta-point
between the collapsed and fully solvated extremes.

Intended use
---------------------------------------------------------

.. note::

   The AFRC is a **reference (null) model, not a predictor** of unfolded-protein dimensions.
   Real dimensions depend on sequence-encoded chain-chain and chain-solvent interactions that
   the AFRC deliberately omits. Its value is as a fixed, sequence-matched touchstone:
   deviations of a simulation or experiment *from* the AFRC are a direct readout of
   sequence-specific intramolecular interactions, and normalising to the AFRC lets chains of
   different lengths and compositions be compared on a common footing.

Parameters
---------------------------------------------------------

The AFRC is deliberately **parameter-free**: the per-residue calibration constants
(:math:`R_0`, :math:`R_0^{\mathrm{rms}}`, :math:`X_0`) are fixed and the only sequence input
is composition and length. There is consequently nothing to tune.

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Argument
     - Default
     - Meaning and typical values
   * - ``adaptable_P_res``
     - ``False``
     - Numerical only. If ``True`` the distribution grid spacing is set to
       :math:`d_{max}/500` (with :math:`d_{max} = 3.7N`) rather than the fixed 0.05 Å.
       Does not change the model, only its discretisation.

**What to expect for a protein.** The apparent scaling exponent is :math:`\nu^{app} = 0.5`
by construction. With :math:`R_0 \approx 6` Å, a disordered region of :math:`N` residues has
:math:`R_e \approx 6\sqrt{N}` Å and :math:`R_g \approx R_e/\sqrt{6} \approx 2.5\sqrt{N}` Å.
Real intrinsically disordered regions *scatter around* these values: in the original study the
ratio of simulated/measured to AFRC dimensions ranged from roughly 0.7 (more compact) to 1.4
(more expanded), so the AFRC is best read as a theta-point touchstone rather than a strict
bound.

Citations
---------------------------------------------------------

1. Alston, J. J., Ginell, G. M., Soranno, A., & Holehouse, A. S. (2023). The Analytical Flory
   Random Coil is a simple-to-use reference model for unfolded and disordered proteins.
   *The Journal of Physical Chemistry B*, 127(21), 4746-4760.
   https://doi.org/10.1021/acs.jpcb.3c01619
2. Flory, P. J. (1969). *Statistical Mechanics of Chain Molecules*. Wiley-Interscience.
3. Mao, A. H., Lyle, N., & Pappu, R. V. (2013). Describing sequence-ensemble relationships
   for intrinsically disordered proteins. *Biochemical Journal*, 449(2), 307-318.
4. Lhuillier, D. (1988). A simple model for polymeric fractals in a good solvent and an
   improved version of the Flory approximation. *Journal de Physique*, 49(5), 705-710.
5. Rubinstein, M., & Colby, R. H. (2003). *Polymer Physics*. Oxford University Press.
6. Nygaard, M., Kragelund, B. B., Papaleo, E., & Lindorff-Larsen, K. (2017). An efficient
   method for estimating the hydrodynamic radius of disordered protein conformations.
   *Biophysical Journal*, 113(3), 550-557.
7. Kirkwood, J. G., & Riseman, J. (1948). The intrinsic viscosities and diffusion constants
   of flexible macromolecules in solution. *The Journal of Chemical Physics*, 16(6), 565-573.
