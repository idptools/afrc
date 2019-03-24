
AnalyticalFRC
=========================================================

The Analytical Flory Random Coil (AFRC) provides a reference state that reproduces the dimensions of a polypeptide in a :math:`\theta`-solvent, and provides a touchstone for comparison with experimental and computational results alike. 

A Flory Random Coil reflects a polymer behaving as a true Gaussian chain - one for which the scaling exponent (:math:`\nu`) is 0.5 [Holehouse2017]_. Even for polypeptides in a (:math:`\theta`)-solvent local dihedral preferences give rise to small sequence-specific differences regarding local and global dimensions. To get around this challenge we had previously developed a numerical approach to generate FRC ensembles via explicit chain simulations that take advantage of the rotational isomeric state approximation developed by Flory and Volkenstein.

In short, the rotational isomeric state approximation states that for a polymer that consists of *n* monomeric units, where each unit can exist in one of *i* possible isomeric states, it is possible to generate an ensemble of states by evolving the system through a series of Monte Carlo moves that 1) randomly select a monomer and 2) and randomly select an isomeric state associated with that monomer. In this way, one can numerically reproduce the conformational ensemble expected for a Gaussian chain - one in which neither repulsive nor attractive interactions are experienced by the chain, but local sterics (which determine the possible isomeric states) coupled with the intrinsic bond lengthscales define the experienced dimensions. This is an attractive approach as it allows you to bake local steric behaviour directly into the chain while ensuring a *bona fide* Gaussian chain ensemble is generated.

A numerical implementation of the FRC has been used extensively as a reference state (see [Mao2013]_, [Das2013]_, [Holehouse2015]_). However, obtaining this ensemble requires numerical simulation at all-atom resolution.

The analytical FRC takes advantage of the fact that for the FRC the finite-size effects that invalidate standard analytical polymer models when describing finite chains no longer apply. Specifically, each monomer is treated as an entirely independent unit, such that monomers at the ends of the chain sample example the same conformational space as thsoe in the center, mittigate finite chain effects. Consequently it is possible the fully parameterize an analytical model to reproduce the numerical simulations over all length-scales. This analytical FRC (AFRC) is implemented using the standard Gaussian chain model to obtain end-to-end distances [Rubinstein2003]_ and an analytical model for the distribution of the radius of gyration to obtain Rg values [Lhuillier1988]_

.. automodule:: afrc
.. autoclass:: AnalyticalFRC
   :members: 

      .. automethod:: __init__



.. rubric:: References
.. [Holehouse2017] Holehouse, A.S., and Pappu, R.V. (2018). Collapse Transitions of Proteins and the Interplay Among Backbone, Sidechain, and Solvent Interactions. Annu. Rev. Biophys. 47, 19-39.
.. [Mao2013] Mao, A.H., Lyle, N., and Pappu, R.V. (2013). Describing sequence-ensemble relationships for intrinsically disordered proteins. Biochem. J 449, 307-318.
.. [Das2013] Das, R.K., and Pappu, R.V. (2013). Conformations of intrinsically disordered proteins are influenced by linear sequence distributions of oppositely charged residues. Proc. Natl. Acad. Sci. U. S. A. 110, 13392-13397.
.. [Holehouse2015] Holehouse, A.S., Garai, K., Lyle, N., Vitalis, A., and Pappu, R.V. (2015). Quantitative assessments of the distinct contributions of polypeptide backbone amides versus side chain groups to chain expansion via chemical denaturation. J. Am. Chem. Soc. 137, 2984-2995.
.. [Rubinstein2003] Rubinstein, M., and Colby, R.H. (2003). Polymer Physics (New York: Oxford University Press).
.. [Lhuillier1988] Lhuillier, D. (1988). A Simple-Model for Polymeric Fractals in a Good Solvent and an Improved Version of the Flory Approximation. Journal De Physique 49, 705-710.
.. [Zhou2004] Zhou, H.-X. (2004). Polymer models of protein stability, folding, and interactions. Biochemistry 43, 2141-2154.

