Worm-like chain (Zhou)
=========================================================

Usage examples and full code reference for
:class:`~afrc.polymer_models.wlc.WormLikeChain`. For the underlying theory, see
:doc:`../polymer_models/worm_like_chain_zhou`.

Quick start
---------------------------------------------------------

.. code-block:: python

   from afrc.polymer_models.wlc import WormLikeChain

   # defaults: persistence length lp = 3.0 A, segment size aa_size = 3.8 A
   model = WormLikeChain('MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYG')

   model.get_mean_end_to_end_distance()
   model.get_root_mean_squared_end_to_end_distance()

   # full end-to-end distribution (distances, probabilities)
   r, p = model.get_end_to_end_distribution()

Explore the effect of chain stiffness by varying the persistence length:

.. code-block:: python

   for lp in (2.0, 3.0, 4.0):
       wlc = WormLikeChain('MASNDYTQQATQSYG', lp=lp)
       print(lp, wlc.get_mean_end_to_end_distance())

See also the ``demo/demo_WormLikeChain.ipynb`` notebook for a worked, plotted example.

Code reference
---------------------------------------------------------

.. autoclass:: afrc.polymer_models.wlc.WormLikeChain
   :members:

   .. automethod:: __init__
