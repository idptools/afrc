
Overview
=========================================================

The Analytical Flory Random Coil (AFRC) reproduces the dimensions of a polypeptide in a :math:`\theta`-solvent, giving a fixed reference against which experimental and computational results can be compared.

A Flory Random Coil (FRC) is a polymer that behaves as a true Gaussian chain, with a scaling exponent :math:`\nu` of 0.5 [Holehouse2018]_. Even in a :math:`\theta`-solvent, local dihedral preferences give polypeptides small sequence-specific differences in local and global dimensions. FRC ensembles capture these by explicit chain simulation using the rotational isomeric state (RIS) approximation of Flory and Volkenstein ([Flory1969]_ and [Volkenstein1977]_).

In the RIS approximation each monomer can occupy one of a set of isomeric states. Monte Carlo moves repeatedly pick a monomer at random and assign it a random allowed state. The chain experiences no attractive or repulsive interactions - only the local sterics that define the allowed states, together with its bond lengths - so the result is a *bona fide* Gaussian chain ensemble with local steric behaviour built in.

The numerical FRC has been used widely as a reference state (see [Mao2013]_, [Das2013]_, [Holehouse2015]_), but generating it requires all-atom simulations. The AFRC removes that step. Because every FRC monomer is independent, residues at the chain ends sample the same conformations as those in the middle, and the finite-size effects that normally break analytical polymer models are absent. An analytical model can therefore be parameterized to match the simulations at every length scale. The AFRC uses the standard Gaussian chain model for end-to-end distances [Rubinstein2003]_ and an analytical radius of gyration distribution for :math:`R_g` [Lhuillier1988]_.


.. rubric:: References
.. [Holehouse2018] Holehouse, A.S., and Pappu, R.V. (2018). Collapse Transitions of Proteins and the Interplay Among Backbone, Sidechain, and Solvent Interactions. Annu. Rev. Biophys. 47, 19-39.
.. [Mao2013] Mao, A.H., Lyle, N., and Pappu, R.V. (2013). Describing sequence-ensemble relationships for intrinsically disordered proteins. Biochem. J 449, 307-318.
.. [Das2013] Das, R.K., and Pappu, R.V. (2013). Conformations of intrinsically disordered proteins are influenced by linear sequence distributions of oppositely charged residues. Proc. Natl. Acad. Sci. U. S. A. 110, 13392-13397.
.. [Holehouse2015] Holehouse, A.S., Garai, K., Lyle, N., Vitalis, A., and Pappu, R.V. (2015). Quantitative assessments of the distinct contributions of polypeptide backbone amides versus side chain groups to chain expansion via chemical denaturation. J. Am. Chem. Soc. 137, 2984-2995.
.. [Rubinstein2003] Rubinstein, M., and Colby, R.H. (2003). Polymer Physics (New York: Oxford University Press).
.. [Lhuillier1988] Lhuillier, D. (1988). A Simple-Model for Polymeric Fractals in a Good Solvent and an Improved Version of the Flory Approximation. Journal De Physique 49, 705-710.
.. [Zhou2004] Zhou, H.-X. (2004). Polymer models of protein stability, folding, and interactions. Biochemistry 43, 2141-2154.
.. [Flory1969] Flory, P. J. (1969). Statistical Mechanics of Chain Molecules. Oxford University Press.
.. [Volkenstein1977] Volkenstein, M. V. (1977). Molecular Biophysics. Academic Press, New York.
