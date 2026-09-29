afrc – the Analytical Flory Random Coil
=========================================================
afrc is a Python package that computes the polymer properties of an unfolded polypeptide using an analytical implementation of the Flory Random Coil (the AFRC).

The AFRC reproduces the dimensions of a polypeptide in a theta solvent: an ideal chain with a scaling exponent of exactly 0.5 and no finite-size effects. It is a reference model, not a predictor. It tells you how an IDR of a given sequence would behave if chain-chain and chain-solvent interactions exactly cancelled out, which makes it a useful null model for simulations and experiments alike.

The model is pre-parameterized against numerical simulations, so the only input is an amino acid sequence. From that, the AFRC instantly returns:

1. The mean end-to-end distance and its distribution.
2. The mean radius of gyration and its distribution.
3. The mean hydrodynamic radius.
4. Every inter-residue mean distance and distance distribution (distance maps and internal scaling profiles).
5. Inter-residue contact fractions (contact maps).
6. Expected paramagnetic relaxation enhancement (PRE) profiles.

The package also implements several other analytical polymer models: the worm-like chain [Zhou2004]_ [Brien2009]_, the self-avoiding walk [Brien2009]_, a SAW with a tunable scaling exponent [Zheng2018]_, the freely jointed chain and the freely rotating chain. They share a common interface with the AFRC, so one sequence can be compared against several reference models. See :doc:`polymer_models/index` for the theory behind each and :doc:`polymer_models_application/index` for usage.


.. toctree::
   :maxdepth: 3
   :caption: Contents:

   overview
   installation
   quickstart
   polymer_models/index
   polymer_models_application/index



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. rubric:: References
.. [Brien2009] O’Brien, E. P., Morrison, G., Brooks, B. R., & Thirumalai, D. (2009). How accurate are polymer models in the analysis of Forster resonance energy transfer experiments on proteins? The Journal of Chemical Physics, 130(12), 124903.
.. [Zheng2018] Zheng, W., Zerze, G. H., Borgia, A., Mittal, J., Schuler, B., & Best, R. B. (2018). Inferring properties of disordered chains from FRET transfer efficiencies. The Journal of Chemical Physics, 148(12), 123329.
