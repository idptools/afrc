Freely jointed chain
=========================================================

Usage examples and full code reference for
:class:`~afrc.polymer_models.fjc.FreelyJointedChain`. For the underlying theory, see
:doc:`../polymer_models/freely_jointed_chain`.

Quick start
---------------------------------------------------------

.. code-block:: python

   from afrc.polymer_models.fjc import FreelyJointedChain

   # default segment length b = 3.8 A (one Cα-Cα virtual bond per residue)
   model = FreelyJointedChain('MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYG')

   model.get_mean_end_to_end_distance()
   model.get_root_mean_squared_end_to_end_distance()
   model.get_mean_radius_of_gyration()

   # full end-to-end distribution (distances, probabilities)
   r, p = model.get_end_to_end_distribution()

   # draw a size-matched sample
   samples = model.sample_end_to_end_distribution(n=1000)

Change the segment length to make the chain stiffer or more flexible:

.. code-block:: python

   stiff = FreelyJointedChain('MASNDYTQQATQSYG', b=4.5)
   stiff.get_root_mean_squared_end_to_end_distance()

See also the ``demo/demo_FreelyJointedChain.ipynb`` notebook for a worked, plotted example.

Generating 3D ensembles
---------------------------------------------------------

Generate explicit conformations (one bead per residue), write them to disk, and check how well they reproduce the model:

.. code-block:: python

   model = FreelyJointedChain('MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYG')

   xyz = model.sample_conformations(n=5000, seed=1)                 # [n x N x 3] array, Angstroms
   model.save_ensemble('ens', n=5000, seed=1)                       # ens.pdb + ens.xtc (needs mdtraj)
   print(model.check_ensemble(xyz).format())                  # the model-check report

Or from the command line: ``afrc-ensemble -s <sequence> -m fjc`` (see :doc:`../cli`).

Code reference
---------------------------------------------------------

.. autoclass:: afrc.polymer_models.fjc.FreelyJointedChain
   :members:

   .. automethod:: __init__
