Worm-like chain (O'Brien)
=========================================================

A second worm-like chain implementation,
:class:`~afrc.polymer_models.wlc2.WormLikeChain2`, using the closed form of
O'Brien et al. (2009). It is conceptually equivalent to the :doc:`Zhou model
<worm_like_chain_zhou>` but is numerically better behaved at large contour lengths and
additionally provides a closed-form radius of gyration.

Mathematical formalism
---------------------------------------------------------

With contour length :math:`L_c = N b` and :math:`\alpha = 3 L_c / (4 L_p)`, the end-to-end
distribution is

.. math::

   P(r) = \frac{4\pi C_1\, r^2}{L_c \left(1 - (r/L_c)^2\right)^{9/2}}
          \exp\!\left( -\frac{3 L_c}{4 L_p \left(1 - (r/L_c)^2\right)} \right),

where the normalisation constant is

.. math::

   C_1 = \left[ \pi^{3/2} e^{-\alpha} \alpha^{-3/2}
        \left( 1 + 3\alpha^{-1} + \tfrac{15}{4}\alpha^{-2} \right) \right]^{-1}.

The :math:`\left(1 - (r/L_c)^2\right)` factors enforce finite extensibility
(:math:`r < L_c`). The radius of gyration is given in closed form (with
:math:`C_2 = 1/(2 L_p)`):

.. math::

   \langle R_g^2 \rangle = \frac{L_c}{6 C_2} + \frac{1}{4 C_2^2}
        + \frac{1}{4 C_2^3 L_c}
        - \frac{1 - e^{-L_c/L_p}}{8 C_2^4 L_c^2},

and :meth:`get_mean_radius_of_gyration` returns :math:`\sqrt{\langle R_g^2 \rangle}`.

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
     - Persistence length (chain stiffness). As for the Zhou model, 3-5 Å is typical for
       unfolded polypeptides; larger values give a stiffer, more extended chain.
   * - ``aa_size``
     - 3.8 Å
     - Segment length :math:`b` (Cα-Cα distance); sets the contour length
       :math:`L_c = N b`.

The constructor additionally requires the sequence to be at least as long as the persistence
length (:math:`N \ge L_p`); otherwise a ``WLC2Exception`` is raised.

**What to expect for a protein.** Results closely track the Zhou model for typical
disordered-protein parameters, with the practical advantages of numerical stability for long
chains and a directly available :math:`R_g`.

Citations
---------------------------------------------------------

1. O'Brien, E. P., Morrison, G., Brooks, B. R., & Thirumalai, D. (2009). How accurate are
   polymer models in the analysis of Förster resonance energy transfer experiments on
   proteins? *The Journal of Chemical Physics*, 130(12), 124903.
2. Rubinstein, M., & Colby, R. H. (2003). *Polymer Physics*. Oxford University Press.
