Analytical Flory Random Coil
=========================================================

Usage examples and full code reference for :class:`~afrc.AnalyticalFRC`. For the underlying
theory, see :doc:`../polymer_models/afrc`.

Quick start
---------------------------------------------------------

Build an object from a sequence and read off ensemble-average dimensions:

.. code-block:: python

   from afrc import AnalyticalFRC

   P = AnalyticalFRC('MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYG')

   P.get_mean_end_to_end_distance()                 # mean Re (A)
   P.get_mean_radius_of_gyration()                  # mean Rg (A)
   P.get_mean_hydrodynamic_radius()                 # mean Rh (A), Kirkwood-Riseman
   P.get_mean_hydrodynamic_radius('nygaard')        # mean Rh (A), Nygaard et al. (2017)

Pull out full probability distributions (each returns ``(distances, probabilities)``):

.. code-block:: python

   re_r, re_p = P.get_end_to_end_distribution()
   rg_r, rg_p = P.get_radius_of_gyration_distribution()

   # distribution between two specific residues
   d_r, d_p = P.get_interresidue_distance_distribution(4, 40)

Inter-residue quantities and whole-chain maps (residue indices start at 0):

.. code-block:: python

   P.get_mean_interresidue_distance(4, 40)  # mean distance between two residues (A)
   P.get_contact_fraction(4, 40, 10.0)      # fraction of time they are within 10 A
   P.get_internal_scaling()                 # [|i-j|, mean distance] profile
   P.get_distance_map()                     # n x n mean inter-residue distances
   P.get_contact_map(15.0)                  # contact fractions at a 15 A threshold
   P.get_pre_profile(0)                     # expected PRE profile for a spin label at residue 0

Draw a size-matched sample (e.g. to compare against a simulation trajectory):

.. code-block:: python

   samples = P.sample_end_to_end_distribution(n=5000)

See also the ``demo/demo_AnalyticalFRC.ipynb`` notebook for a worked, plotted example.

Generating 3D ensembles
---------------------------------------------------------

The AFRC can also generate explicit 3D conformations, one bead per residue, drawn exactly from the model (see :doc:`../polymer_models/afrc` for how). Save them as a PDB/XTC pair for analysis or visualization, or work with the coordinates directly:

.. code-block:: python

   # writes p53.pdb (topology + first conformation) and p53.xtc (all 5000 conformations)
   xyz = P.save_ensemble('p53', n=5000, seed=1)

   # or just the coordinates: an [n x N x 3] array in Angstroms
   xyz = P.sample_conformations(n=5000, seed=1)

Writing the XTC file needs mdtraj (``pip install "afrc[ensemble]"``). The pair loads directly into mdtraj or SOURSOP:

.. code-block:: python

   import mdtraj as md
   traj = md.load('p53.xtc', top='p53.pdb')

   from soursop.sstrajectory import SSTrajectory
   protein = SSTrajectory('p53.xtc', 'p53.pdb').proteinTrajectoryList[0]

Note that adjacent beads are not at a fixed 3.8 Å spacing and can overlap: the ensemble follows the AFRC's own inter-residue distance distributions, including for neighbouring residues. Pass ``pdb_only=True`` to ``save_ensemble()`` to write every conformation to a single multi-model PDB file instead (no mdtraj needed, at most 9999 conformations).

Checking an ensemble
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``check_ensemble()`` reports how well an ensemble reproduces the AFRC - every quantity the model fixes exactly, each with a standard error estimated from the ensemble:

.. code-block:: python

   report = P.check_ensemble(xyz)
   report.passed          # True for a model-like ensemble
   print(report.format())

The same thing is available from the command line, for this and every other model, with ``afrc-ensemble`` - see :doc:`../cli`.

Code reference
---------------------------------------------------------

.. autoclass:: afrc.AnalyticalFRC
   :members:

   .. automethod:: __init__

Ensemble utilities
---------------------------------------------------------

The functions behind ``sample_conformations()``, ``save_ensemble()`` and ``afrc-ensemble``, which can also be used to write, or check, your own one-bead-per-residue coordinates.

.. automodule:: afrc.ensemble
   :members: write_ensemble, write_pdb, write_multimodel_pdb, write_xtc, save_conformations, gaussian_chain_factor, sample_gaussian_chain, sample_freely_jointed_chain, sample_freely_rotating_chain, sample_worm_like_chain, freely_rotating_chain_msd, worm_like_chain_msd, discrete_worm_like_chain_msd, worm_like_chain_discretization, mean_squared_distance_map

.. automodule:: afrc.ensemble_report
   :members: compare_ensemble_to_model, compare_ensemble_to_gaussian_model, ModelExpectations, EnsembleReport, Check
