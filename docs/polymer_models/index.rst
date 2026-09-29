Polymer Models (Theory)
=========================================================

Alongside the Analytical Flory Random Coil, the ``afrc`` package implements several other analytical polymer models. Each takes an amino acid sequence and returns an end-to-end distance distribution plus mean values through a common interface (``get_end_to_end_distribution``, ``get_mean_end_to_end_distance``, ...).

For each model, the pages below give (1) the formalism that is actually implemented, (2) the free parameters, what they mean and sensible values for a polypeptide, and (3) the primary references. For usage examples and the code reference, see :doc:`Polymer Models (Application) <../polymer_models_application/index>`.

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
     - Semiflexible chain parameterized by a persistence length.
   * - :doc:`Worm-like chain (O'Brien) <worm_like_chain_obrien>`
     - ``WormLikeChain2``
     - Semiflexible chain; exact finite extensibility, also gives Rg.
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

Throughout, :math:`N` is the number of residues, :math:`r` is the end-to-end distance, :math:`R_e` the end-to-end distance and :math:`R_g` the radius of gyration. All distances are in Angstroms.

Distributions are returned as discrete, normalized probability mass functions ``(distances, probabilities)`` on a grid whose spacing is set by ``p_of_r_resolution`` (0.05 Å by default). This is a numerical discretization, not a model parameter.

``get_mean_end_to_end_distance`` always returns the mean, :math:`\langle R_e \rangle`, and ``get_root_mean_squared_end_to_end_distance`` returns :math:`\sqrt{\langle R_e^2 \rangle}`. For the AFRC, ``get_mean_radius_of_gyration`` returns the true mean :math:`\langle R_g \rangle`. The other models have no :math:`R_g` distribution, so their ``get_mean_radius_of_gyration`` returns the root-mean-square value :math:`\sqrt{\langle R_g^2 \rangle}` from a closed-form ratio; for an ideal chain this is about 3% larger than :math:`\langle R_g \rangle`.
