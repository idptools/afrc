.. afrc documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

afrc – the Analytical Flory Random Coil
=========================================================
afrc is a Python-based package for computing polymeric properties for unfolded polypeptides using an analytical implementation of the so-called Flory Random Coil (the AFRC). 

Briefly, the AFRC is a pre-parameterized polymer model that recapitulates the dimensions of a polypeptide in a theta solvent. Technically speaking, this means both the second and third virial coefficients are set to zero, such that the AFRC enjoys fractal scaling with a true scaling exponent of 0.5 and no finite-size effects. This makes it well-suited as a reference mode for developing intuition, providing a comparison against experimental data, or offering normalization factors for simulations or experiments alike.

We developed the AFRC as a convenient tool for contextualizing simulations and experiments of disordered proteins. The AFRC is *not* a predictor of the dimensions of intrinsically disordered proteins or protein regions, but it does offer a 'null model' for how one might expect an IDR of a given sequence to behave if chain-chain and chain-solvent interactions were perfectly counterbalanced. 

The AFRC is parameterized against numerical simulations that recapitulate *bona fide* theta solvent behavior. As such, you only need to provide an amino acid sequence as input, and the AFRC can provide a variety of information back instantaneously, including:

1. The ensemble-average end-to-end distance.
2. The end-to-end distance distribution.
3. The ensemble-average radius of gyration.
4. The end-to-end distance distribution.
5. The ensemble-average hydrodynamic radius
6. All inter-residue average distances and distance distributions
7. Inter-residue contact fractions


Finally, the afrc package also implements several additional polymer models, including the Worm-like chain (WLC) model [Brien2009]_, the self-avoid walk (SAW) model [Brien2009]_, and a scaling-exponent SAW model (SAW-ν) [Zheng2018]_.


.. toctree::
   :maxdepth: 3
   :caption: Contents:

   overview   
   installation
   quickstart   
   modules/afrc
   


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. rubric:: References
.. [Brien2009] O’Brien, E. P., Morrison, G., Brooks, B. R., & Thirumalai, D. (2009). How accurate are polymer models in the analysis of Forster resonance energy transfer experiments on proteins? The Journal of Chemical Physics, 130(12), 124903.
.. [Zheng2018] Zheng, W., Zerze, G. H., Borgia, A., Mittal, J., Schuler, B., & Best, R. B. (2018). Inferring properties of disordered chains from FRET transfer efficiencies. The Journal of Chemical Physics, 148(12), 123329.
