Freely jointed chain
=========================================================

The freely jointed chain (FJC), exposed through
:class:`~afrc.polymer_models.fjc.FreelyJointedChain`, models the chain as :math:`N` rigid
segments of length :math:`b` connected by perfectly flexible joints. Unlike the Gaussian
AFRC it uses the *non-Gaussian* Kuhn-Grün distribution, which respects the chain's finite
extensibility: the end-to-end distance can never exceed the contour length
:math:`L = N b`.

Mathematical formalism
---------------------------------------------------------

The radial end-to-end probability is

.. math::

   P(r) \propto 4\pi r^2
   \exp\!\left[ -N \left( x\,\beta + \ln\frac{\beta}{\sinh\beta} \right) \right],
   \qquad x = \frac{r}{N b},

where :math:`x \in [0, 1)` is the fractional extension and :math:`\beta = \mathcal{L}^{-1}(x)`
is the inverse Langevin function. The inverse Langevin is evaluated with the Cohen Padé
approximant

.. math::

   \beta = \mathcal{L}^{-1}(x) \approx \frac{x\,(3 - x^2)}{1 - x^2},

which diverges as :math:`x \to 1`, correctly suppressing all probability beyond the contour
length. At small extension the exponent reduces to :math:`\tfrac{3}{2} N x^2`, recovering the
Gaussian chain with :math:`\langle R_e^2 \rangle = N b^2`; the root-mean-square size therefore
approaches :math:`b\sqrt{N}` from below. The radius of gyration uses the ideal-chain relation
:math:`R_g = \sqrt{\langle R_e^2 \rangle}/\sqrt{6}`.

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
     - Segment (Kuhn) length. The default corresponds to the Cα-Cα distance, i.e. one
       residue per segment. To represent a stiffer effective chain one can use a larger
       Kuhn length (a polypeptide Kuhn length is often quoted around 7-10 Å); :math:`N`
       and :math:`b` together fix the contour length :math:`L = Nb`.

**What to expect for a protein.** With one segment per residue (:math:`b = 3.8` Å) the FJC is
an ideal chain: :math:`R_e \approx b\sqrt{N}` and :math:`R_g \approx R_e/\sqrt{6}`, essentially
matching the AFRC through the bulk of the distribution. The differences appear in the far tail
(the FJC has a hard cutoff at :math:`L = Nb`) and at short chain lengths, where finite
extensibility pulls the mean and RMS slightly below the Gaussian values.

Citations
---------------------------------------------------------

1. Kuhn, W., & Grün, F. (1942). Beziehungen zwischen elastischen Konstanten und
   Dehnungsdoppelbrechung hochelastischer Stoffe. *Kolloid-Zeitschrift*, 101(3), 248-271.
2. Cohen, A. (1991). A Padé approximant to the inverse Langevin function.
   *Rheologica Acta*, 30(3), 270-273.
3. Rubinstein, M., & Colby, R. H. (2003). *Polymer Physics*. Oxford University Press.
