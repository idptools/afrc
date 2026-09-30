Quickstart
=========================================================

Everything in the AFRC is accessed through an ``AnalyticalFRC`` object, built from an amino acid sequence. More examples are in the `demo directory <https://github.com/idptools/afrc/tree/main/demo>`_.

.. code-block:: python

   from afrc import AnalyticalFRC

   protein = AnalyticalFRC('MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYG')

   # mean dimensions, all in Angstroms
   mean_re = protein.get_mean_end_to_end_distance()
   mean_rg = protein.get_mean_radius_of_gyration()
   mean_rh = protein.get_mean_hydrodynamic_radius()

   # distributions are returned as (distances, probabilities)
   re_distances, re_probabilities = protein.get_end_to_end_distribution()
   rg_distances, rg_probabilities = protein.get_radius_of_gyration_distribution()

   # internal scaling profile: [|i-j|, mean distance] for every separation
   internal_scaling = protein.get_internal_scaling()

   # contact fractions at a 15 A threshold for every pair of residues
   contact_map = protein.get_contact_map(15.0)

   # 1000 3D conformations (one bead per residue), written as ens.pdb + ens.xtc
   # (the XTC needs mdtraj), and a report on how well they reproduce the model
   xyz = protein.save_ensemble('ens', n=1000, seed=1)
   print(protein.check_ensemble(xyz).format())

Residue indices start at 0. See :doc:`polymer_models_application/afrc` for the full set of methods, and :doc:`cli` for generating ensembles from the command line with ``afrc-ensemble``.
