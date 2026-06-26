Self-avoiding walk
=========================================================

Usage examples and full code reference for :class:`~afrc.polymer_models.saw.SAW`. For the
underlying theory, see :doc:`../polymer_models/self_avoiding_walk`.

Quick start
---------------------------------------------------------

.. code-block:: python

   from afrc.polymer_models.saw import SAW

   model = SAW('MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYG')

   model.get_mean_end_to_end_distance()
   model.get_root_mean_squared_end_to_end_distance()
   model.get_mean_radius_of_gyration()

   # full end-to-end distribution (distances, probabilities)
   r, p = model.get_end_to_end_distribution()

The dimensions are tuned with the ``prefactor`` argument (default 5.5 A), which can be passed
to any of the methods:

.. code-block:: python

   for pref in (4.5, 5.5, 6.5):
       print(pref, model.get_mean_end_to_end_distance(prefactor=pref))

See also the ``demo/demo_SAW.ipynb`` notebook for a worked, plotted example.

Code reference
---------------------------------------------------------

.. autoclass:: afrc.polymer_models.saw.SAW
   :members:

   .. automethod:: __init__
