Worm-like chain (O'Brien)
=========================================================

Usage examples and full code reference for
:class:`~afrc.polymer_models.wlc2.WormLikeChain2`. For the underlying theory, see
:doc:`../polymer_models/worm_like_chain_obrien`.

Quick start
---------------------------------------------------------

.. code-block:: python

   from afrc.polymer_models.wlc2 import WormLikeChain2

   # defaults: persistence length lp = 3.0 A, segment size aa_size = 3.8 A
   model = WormLikeChain2('MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYG')

   model.get_mean_end_to_end_distance()
   model.get_root_mean_squared_end_to_end_distance()
   model.get_mean_radius_of_gyration()       # closed-form root-mean-square Rg

   r, p = model.get_end_to_end_distribution()

Vary the persistence length:

.. code-block:: python

   for lp in (2.0, 3.0, 4.0):
       wlc = WormLikeChain2('MASNDYTQQATQSYG', lp=lp)
       print(lp, wlc.get_mean_radius_of_gyration())

.. note::

   ``WormLikeChain2`` requires the chain's contour length (``len(seq) * aa_size``) to be at
   least one persistence length, otherwise a ``WLC2Exception`` is raised.

See also the ``demo/demo_WormLikeChain2.ipynb`` notebook for a worked, plotted example.

Generating 3D ensembles
---------------------------------------------------------

Generate explicit conformations (one bead per residue), write them to disk, and check how well they reproduce the model:

.. code-block:: python

   model = WormLikeChain2('MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYG', lp=3.0)

   xyz = model.sample_conformations(n=5000, seed=1)                 # [n x N x 3] array, Angstroms
   model.save_ensemble('ens', n=5000, seed=1)                       # ens.pdb + ens.xtc (needs mdtraj)
   print(model.check_ensemble(xyz).format())                  # the model-check report

Or from the command line: ``afrc-ensemble -s <sequence> -m wlc2`` (see :doc:`../cli`).

Code reference
---------------------------------------------------------

.. autoclass:: afrc.polymer_models.wlc2.WormLikeChain2
   :members:

   .. automethod:: __init__
