Command-line tools
=========================================================

``afrc-ensemble``
---------------------------------------------------------

Installing afrc also installs ``afrc-ensemble``, which generates a 3D conformational ensemble (one bead per residue) from any of the package's polymer models, writes it to disk, and reports how well it reproduces the model:

.. code-block:: bash

   afrc-ensemble -s MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDI -n 5000 -o p53

This writes ``p53.pdb`` (topology and first conformation), ``p53.xtc`` (every conformation) and ``p53_report.txt``. The PDB/XTC pair loads directly into mdtraj or SOURSOP; writing the XTC needs mdtraj (``pip install "afrc[ensemble]"``).

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Option
     - Meaning
   * - ``-s``, ``--sequence``
     - Amino acid sequence (required).
   * - ``-n``, ``--number-of-conformers``
     - Number of conformations (default 1000).
   * - ``-m``, ``--model``
     - Polymer model to draw from (default ``afrc``; see below).
   * - ``-o``, ``--out``
     - Output prefix; writes ``PREFIX.pdb``, ``PREFIX.xtc`` and ``PREFIX_report.txt`` (default ``out``).
   * - ``--pdb-only``
     - Write a single multi-model PDB instead of a PDB/XTC pair (no mdtraj needed; at most 9999 conformations).
   * - ``--seed``
     - Random seed. Without one a fresh seed is drawn and recorded in the report, so any run can be reproduced.
   * - ``--segment-length``
     - Length per residue in Å: the bond length ``b`` (``fjc``, ``frc``) or ``aa_size`` (``wlc``, ``wlc2``). Default 3.8.
   * - ``--c-inf``
     - Characteristic ratio (``frc``). Default 2.0.
   * - ``--lp``
     - Persistence length in Å (``wlc``, ``wlc2``). Default 3.0.
   * - ``--prefactor``
     - Size prefactor in Å (``saw``, ``saw-nu``). Default 5.5.
   * - ``--nu``
     - Flory scaling exponent, between 0 and 1 (``saw-nu``). Default 0.5.

A model parameter given for a model that does not use it (``--lp`` with ``-m afrc``, say) is an error.

Models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 12 30 58

   * - ``-m``
     - Model
     - How the ensemble is generated
   * - ``afrc``
     - :doc:`Analytical Flory Random Coil <polymer_models/afrc>`
     - Exactly: every inter-residue distance is Gaussian, so bead coordinates are drawn from the multivariate Gaussian those distances define.
   * - ``fjc``
     - :doc:`Freely jointed chain <polymer_models/freely_jointed_chain>`
     - Exactly: bonds of length ``b`` in independent, uniformly random directions.
   * - ``frc``
     - :doc:`Freely rotating chain <polymer_models/freely_rotating_chain>`
     - Exactly: bonds of length ``b`` at a fixed bond angle (set by ``c_inf``) with uniformly random torsions.
   * - ``wlc``, ``wlc2``
     - Worm-like chain (:doc:`Zhou <polymer_models/worm_like_chain_zhou>`, :doc:`O'Brien <polymer_models/worm_like_chain_obrien>`)
     - Each residue is split into short straight sub-segments that bend like a worm-like chain, finely enough that the mean-squared distances match the continuous chain to within 0.02% (chains that are flexible on the scale of a residue use a correlation between sub-segments that gets there with far fewer of them). Both options generate the same ensemble; they differ only in which analytical P(r) the report compares it with.
   * - ``saw``, ``saw-nu``
     - :doc:`Self-avoiding walk <polymer_models/self_avoiding_walk>`, :doc:`nu-dependent SAW <polymer_models/nu_dependent_saw>`
     - A Gaussian approximation: the mean-squared distances, :math:`(\texttt{prefactor}\,|i-j|^{\nu})^2`, are exact, but the distance distributions are Gaussian rather than the SAW's and there is no excluded volume. Genuine self-avoiding conformations would need an explicit simulation.

The report
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The report (printed and saved) checks the ensemble against everything the model fixes exactly, each with a standard error estimated from the ensemble itself:

* for every model, the mean-squared distance of every residue pair (in bands of sequence separation) and the root-mean-square :math:`R_g` and first-to-last distance that follow from them;
* for the AFRC, whose distances are Gaussian, also the mean distances, the Kirkwood-Riseman :math:`R_h` and Kolmogorov-Smirnov tests of individual distance distributions;
* for models with rigid geometry, exact structural checks: bond lengths (``fjc``, ``frc``), bond angles (``frc``) and finite extensibility (``fjc``, ``frc``, ``wlc``, ``wlc2``).

A model-like ensemble passes every check (it fails beyond 5 standard errors); across thousands of simulated correct ensembles of every model, none was flagged, while for a 50-residue chain a 2% error in the chain's size is always caught with 1000 conformations and a 1% error with 5000. The report also lists context that is not expected to match exactly - the model's own radius of gyration and whole-chain end-to-end distance, and a comparison with its analytical end-to-end distribution, which for most models is an approximation - and model-specific notes, such as that the SAW ensembles are a Gaussian approximation. Ensembles of fewer than 100 conformations are reported but not assessed.

From Python, every model offers the same pieces: ``sample_conformations()``, ``save_ensemble()``, ``check_ensemble()`` (which returns the report) and ``get_mean_squared_distance_map()``.
