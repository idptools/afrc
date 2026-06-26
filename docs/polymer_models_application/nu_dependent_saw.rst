nu-dependent self-avoiding walk
=========================================================

Usage examples and full code reference for
:class:`~afrc.polymer_models.nudep_saw.NuDepSAW`. For the underlying theory, see
:doc:`../polymer_models/nu_dependent_saw`.

Quick start
---------------------------------------------------------

.. code-block:: python

   from afrc.polymer_models.nudep_saw import NuDepSAW

   model = NuDepSAW('MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYG')

   # the scaling exponent nu is supplied per call (default nu = 0.5)
   model.get_mean_end_to_end_distance(nu=0.5)
   model.get_root_mean_squared_end_to_end_distance(nu=0.5)
   model.get_mean_radius_of_gyration(nu=0.5)

   # full end-to-end distribution (distances, probabilities)
   r, p = model.get_end_to_end_distribution(nu=0.5)

   # draw a size-matched sample
   samples = model.sample_end_to_end_distribution(n=1000, nu=0.5)

Sweep the scaling exponent from collapsed (nu ~ 1/3) to fully expanded (nu ~ 0.588):

.. code-block:: python

   for nu in (0.40, 0.50, 0.59):
       print(nu, model.get_mean_end_to_end_distance(nu=nu))

See also the ``demo/demo_NuDepSAW.ipynb`` notebook for a worked, plotted example.

Code reference
---------------------------------------------------------

.. autoclass:: afrc.polymer_models.nudep_saw.NuDepSAW
   :members:

   .. automethod:: __init__
