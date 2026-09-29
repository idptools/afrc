Worm-like chain (Zhou)
=========================================================

The worm-like chain (WLC) describes a semiflexible polymer whose stiffness is set by a persistence length :math:`L_p`. This implementation, :class:`~afrc.polymer_models.wlc.WormLikeChain`, uses the closed-form approximation of Zhou (2004). It is composition-independent: the sequence only sets the number of residues.

Mathematical formalism
---------------------------------------------------------

With contour length :math:`L_c = N b`, the end-to-end distribution is

.. math::

   P(r) = 4\pi A\, r^2
          \exp\!\left( -\frac{3 r^2}{4 L_p L_c} \right) \zeta(r),
   \qquad A = \left( \frac{3}{4\pi L_p L_c} \right)^{3/2},

where the Gaussian core is the flexible-limit result :math:`\langle R^2 \rangle = 2 L_p L_c` and :math:`\zeta(r)` is Zhou's polynomial correction series (Zhou 2004, Eq. 5) in powers of :math:`r/L_c` and :math:`L_p/L_c`. Mean and root-mean-square end-to-end distances are computed numerically from :math:`P(r)`.

.. note::

   Because :math:`\zeta(r)` is an expansion, it is only accurate when :math:`L_c \gg L_p` - with the default parameters, chains of more than roughly 10-20 residues. In that regime the distribution reproduces the exact worm-like chain :math:`\langle R^2 \rangle = 2 L_p L_c - 2 L_p^2 (1 - e^{-L_c/L_p})` to high precision. No probability is assigned beyond :math:`L_c`, and negative values of the series are set to zero. If :math:`L_c` is comparable to or shorter than :math:`L_p` the expansion has no valid region, and requesting the distribution raises a ``WLCException``.

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
     - Persistence length - the distance over which the chain "forgets" its direction. Larger values give a stiffer, more extended chain. Values of 3-5 Å are common for unfolded polypeptides (4 Å is also frequent in the literature).
   * - ``aa_size``
     - 3.8 Å
     - Contour length per residue, :math:`b`, so :math:`L_c = N b`. 3.8 Å is the Cα-Cα distance.

**What to expect for a protein.** Because :math:`\langle R^2 \rangle \approx 2 L_p L_c`, the chain grows as :math:`\sqrt{L_p}` and shows Gaussian-coil scaling when :math:`L_c \gg L_p`. With :math:`L_p` in the 3-4 Å range the WLC is roughly 10-25% more compact than the :doc:`AFRC <afrc>`; matching the AFRC's dimensions needs :math:`L_p \approx 5` Å.

Citations
---------------------------------------------------------

1. Zhou, H.-X. (2004). Polymer models of protein stability, folding, and interactions.
   *Biochemistry*, 43(8), 2141-2154.
2. Rubinstein, M., & Colby, R. H. (2003). *Polymer Physics*. Oxford University Press.
