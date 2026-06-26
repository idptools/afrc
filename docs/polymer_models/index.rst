Polymer Models (Theory)
=========================================================

Alongside the headline Analytical Flory Random Coil, the ``afrc`` package implements
several additional analytical polymer models. Each one takes an amino acid sequence and
returns an end-to-end distance distribution together with associated mean values, exposing
a common interface (``get_end_to_end_distribution``, ``get_mean_end_to_end_distance``, ...).

The pages below describe, for each model: (1) the mathematical formalism that is actually
implemented, (2) the free parameters, what they mean physically, and sensible values for a
polypeptide, and (3) the primary references. For runnable usage examples and the full code
reference for each class, see :doc:`Polymer Models (Application)
<../polymer_models_application/index>`.

.. list-table::
   :header-rows: 1
   :widths: 30 25 45

   * - Model
     - Class
     - In one line
   * - :doc:`Analytical Flory Random Coil <afrc>`
     - ``AnalyticalFRC``
     - Sequence-specific ideal (theta-state) chain; the reference null model.
   * - :doc:`Freely jointed chain <freely_jointed_chain>`
     - ``FreelyJointedChain``
     - Ideal chain with finite extensibility (non-Gaussian Kuhn-Grün).
   * - :doc:`Freely rotating chain <freely_rotating_chain>`
     - ``FreelyRotatingChain``
     - Ideal chain with a tunable characteristic ratio (stiffness).
   * - :doc:`Worm-like chain (Zhou) <worm_like_chain_zhou>`
     - ``WormLikeChain``
     - Semiflexible chain parameterised by a persistence length.
   * - :doc:`Worm-like chain (O'Brien) <worm_like_chain_obrien>`
     - ``WormLikeChain2``
     - Semiflexible chain; better large-chain stability, also gives Rg.
   * - :doc:`Self-avoiding walk <self_avoiding_walk>`
     - ``SAW``
     - Good-solvent (excluded-volume) chain at fixed scaling exponent.
   * - :doc:`nu-dependent SAW <nu_dependent_saw>`
     - ``NuDepSAW``
     - Excluded-volume chain with a tunable Flory scaling exponent.

.. toctree::
   :maxdepth: 1
   :caption: Models

   afrc
   freely_jointed_chain
   freely_rotating_chain
   worm_like_chain_zhou
   worm_like_chain_obrien
   self_avoiding_walk
   nu_dependent_saw

A note on conventions
---------------------------------------------------------

Throughout, :math:`N` is the number of residues in the sequence, :math:`r` is the
end-to-end distance, :math:`R_e` the mean end-to-end distance, and :math:`R_g` the radius
of gyration. All distances are in Angstroms. Distributions are returned as discrete,
normalised probability mass functions ``(distances, probabilities)`` evaluated on a grid
whose spacing is set by ``p_of_r_resolution`` (0.05 Å by default); this is a numerical
discretisation parameter, not a model parameter.
