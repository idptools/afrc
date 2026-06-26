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

   P.get_mean_end_to_end_distance()        # mean Re (A)
   P.get_mean_radius_of_gyration()         # mean Rg (A)
   P.get_mean_hydrodynamic_radius()        # mean Rh (A), Kirkwood-Riseman

Pull out full probability distributions (each returns ``(distances, probabilities)``):

.. code-block:: python

   re_r, re_p = P.get_end_to_end_distribution()
   rg_r, rg_p = P.get_radius_of_gyration_distribution()

   # distribution between two specific residues
   d_r, d_p = P.get_interresidue_distance_distribution(4, 40)

Inter-residue and whole-chain maps:

.. code-block:: python

   P.get_internal_scaling()                # [|i-j|, mean distance] profile
   P.get_distance_map()                    # n x n mean inter-residue distances
   P.get_contact_map(15.0)                 # contact fractions at a 15 A threshold
   P.get_pre_profile(0)                    # hypothetical PRE profile for a label at residue 0

Draw a size-matched sample (e.g. to compare against a simulation trajectory):

.. code-block:: python

   samples = P.sample_end_to_end_distribution(n=5000)

See also the ``demo/demo_AnalyticalFRC.ipynb`` notebook for a worked, plotted example.

Code reference
---------------------------------------------------------

.. autoclass:: afrc.AnalyticalFRC
   :members:

   .. automethod:: __init__
