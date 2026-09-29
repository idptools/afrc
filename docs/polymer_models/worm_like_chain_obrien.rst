Worm-like chain (O'Brien)
=========================================================

A second worm-like chain, :class:`~afrc.polymer_models.wlc2.WormLikeChain2`, using the closed form of O'Brien et al. (2009). It describes the same physics as the :doc:`Zhou model <worm_like_chain_zhou>`, but enforces finite extensibility exactly, stays accurate for long and stiff chains, and also provides a closed-form radius of gyration.

Mathematical formalism
---------------------------------------------------------

With contour length :math:`L_c = N b` and :math:`\alpha = 3 L_c / (4 L_p)`, the end-to-end distribution is

.. math::

   P(r) = \frac{4\pi C_1\, r^2}{L_c^3 \left(1 - (r/L_c)^2\right)^{9/2}}
          \exp\!\left( -\frac{3 L_c}{4 L_p \left(1 - (r/L_c)^2\right)} \right),

with normalization constant

.. math::

   C_1 = \left[ \pi^{3/2} e^{-\alpha} \alpha^{-3/2}
        \left( 1 + 3\alpha^{-1} + \tfrac{15}{4}\alpha^{-2} \right) \right]^{-1}.

The :math:`\left(1 - (r/L_c)^2\right)` factors enforce :math:`r < L_c`. Because :math:`C_1 \propto e^{\alpha}` overflows for long chains, the code evaluates :math:`\log P(r)` and normalizes numerically; :math:`C_1` never enters the calculation.

The radius of gyration is given in closed form (with :math:`C_2 = 1/(2 L_p)`):

.. math::

   \langle R_g^2 \rangle = \frac{L_c}{6 C_2} - \frac{1}{4 C_2^2}
        + \frac{1}{4 C_2^3 L_c}
        - \frac{1 - e^{-L_c/L_p}}{8 C_2^4 L_c^2}.

This is the Benoit-Doty worm-like chain result, :math:`L_c L_p/3 - L_p^2 + 2 L_p^3/L_c - 2 L_p^4/L_c^2 \left(1 - e^{-L_c/L_p}\right)`, rewritten in terms of :math:`C_2`. It reduces to :math:`L_c L_p/3` when :math:`L_c \gg L_p` and to the rigid-rod value :math:`L_c^2/12` when :math:`L_c \ll L_p`. :meth:`~afrc.polymer_models.wlc2.WormLikeChain2.get_mean_radius_of_gyration` returns :math:`\sqrt{\langle R_g^2 \rangle}`.

Parameters
---------------------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Parameter
     - Default
     - Meaning and typical values
   * - ``lp``
     - 3.0 Å
     - Persistence length (chain stiffness). As for the Zhou model, 3-5 Å is typical for unfolded polypeptides; larger values give a stiffer, more extended chain.
   * - ``aa_size``
     - 3.8 Å
     - Contour length per residue, :math:`b` (the Cα-Cα distance); sets :math:`L_c = N b`.

The contour length must be at least one persistence length (:math:`N b \ge L_p`), otherwise the constructor raises a ``WLC2Exception``.

**What to expect for a protein.** Results closely track the Zhou model for typical disordered-protein parameters, and match the exact worm-like chain :math:`\langle R^2 \rangle` to within a fraction of a percent across chain lengths and stiffnesses.

Citations
---------------------------------------------------------

1. O'Brien, E. P., Morrison, G., Brooks, B. R., & Thirumalai, D. (2009). How accurate are
   polymer models in the analysis of Förster resonance energy transfer experiments on
   proteins? *The Journal of Chemical Physics*, 130(12), 124903.
2. Benoit, H., & Doty, P. (1953). Light scattering from non-Gaussian chains. *The Journal
   of Physical Chemistry*, 57(9), 958-963.
3. Rubinstein, M., & Colby, R. H. (2003). *Polymer Physics*. Oxford University Press.
