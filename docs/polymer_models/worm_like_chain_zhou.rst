Worm-like chain (Zhou)
=========================================================

The worm-like chain (WLC) describes a semiflexible polymer whose stiffness is set by a
persistence length :math:`L_p`. This implementation,
:class:`~afrc.polymer_models.wlc.WormLikeChain`, uses the closed-form approximation of
Zhou (2004). It is composition-independent: the sequence enters only through the number of
residues.

Mathematical formalism
---------------------------------------------------------

With contour length :math:`L_c = N b`, the end-to-end distribution is

.. math::

   P(r) = 4\pi A\, r^2
          \exp\!\left( -\frac{3 r^2}{4 L_p L_c} \right) \zeta(r),
   \qquad A = \left( \frac{3}{4\pi L_p L_c} \right)^{3/2},

where the Gaussian core corresponds to the flexible-limit WLC result
:math:`\langle R^2 \rangle = 2 L_p L_c`, and :math:`\zeta(r)` is the polynomial correction
series of Zhou (2004, Eq. 5) in powers of :math:`r/L_c`, :math:`L_p/L_c`. The series is
accurate for :math:`r` up to about the contour length; spurious negative values in the far
tail are clamped to zero so the result remains a valid probability distribution. Mean and
root-mean-square end-to-end distances are obtained by numerical integration over
:math:`P(r)`.

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
     - Persistence length - the length scale over which the chain "forgets" its direction.
       Larger :math:`L_p` means a stiffer, more extended chain. For unfolded/disordered
       polypeptides values of roughly 3-5 Å are commonly used (4 Å is also frequent in the
       literature).
   * - ``aa_size``
     - 3.8 Å
     - Segment length :math:`b` (the per-residue contribution to the contour length,
       :math:`L_c = N b`). 3.8 Å is the Cα-Cα distance.

**What to expect for a protein.** With :math:`L_p` in the 3-4 Å range the WLC produces
dimensions comparable to a slightly expanded coil. Because :math:`\langle R^2 \rangle = 2 L_p L_c`,
the chain becomes more extended as :math:`L_p` increases, and recovers Gaussian-coil scaling
when :math:`L_c \gg L_p`.

Citations
---------------------------------------------------------

1. Zhou, H.-X. (2004). Polymer models of protein stability, folding, and interactions.
   *Biochemistry*, 43(8), 2141-2154.
2. Rubinstein, M., & Colby, R. H. (2003). *Polymer Physics*. Oxford University Press.
