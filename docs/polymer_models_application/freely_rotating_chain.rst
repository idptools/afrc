Freely rotating chain
=========================================================

Usage examples and full code reference for
:class:`~afrc.polymer_models.frc.FreelyRotatingChain`. For the underlying theory, see
:doc:`../polymer_models/freely_rotating_chain`.

Quick start
---------------------------------------------------------

.. code-block:: python

   from afrc.polymer_models.frc import FreelyRotatingChain

   # defaults: bond length b = 3.8 A, characteristic ratio c_inf = 2.0
   model = FreelyRotatingChain('MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYG')

   model.get_mean_end_to_end_distance()
   model.get_root_mean_squared_end_to_end_distance()
   model.get_mean_radius_of_gyration()

   r, p = model.get_end_to_end_distribution()
   samples = model.sample_end_to_end_distribution(n=1000)

Tune the stiffness via the characteristic ratio (``c_inf = 1`` recovers the freely jointed
chain; larger values give a stiffer, more extended ideal chain):

.. code-block:: python

   flexible = FreelyRotatingChain('MASNDYTQQATQSYG', c_inf=1.0)
   stiff    = FreelyRotatingChain('MASNDYTQQATQSYG', c_inf=4.0)

   flexible.get_root_mean_squared_end_to_end_distance()
   stiff.get_root_mean_squared_end_to_end_distance()

See also the ``demo/demo_FreelyRotatingChain.ipynb`` notebook for a worked, plotted example.

Generating 3D ensembles
---------------------------------------------------------

Generate explicit conformations (one bead per residue), write them to disk, and check how well they reproduce the model:

.. code-block:: python

   model = FreelyRotatingChain('MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYG', c_inf=2.0)

   xyz = model.sample_conformations(n=5000, seed=1)                 # [n x N x 3] array, Angstroms
   model.save_ensemble('ens', n=5000, seed=1)                       # ens.pdb + ens.xtc (needs mdtraj)
   print(model.check_ensemble(xyz).format())                  # the model-check report

Or from the command line: ``afrc-ensemble -s <sequence> -m frc`` (see :doc:`../cli`).

Code reference
---------------------------------------------------------

.. autoclass:: afrc.polymer_models.frc.FreelyRotatingChain
   :members:

   .. automethod:: __init__
