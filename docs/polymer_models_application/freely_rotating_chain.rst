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

Code reference
---------------------------------------------------------

.. autoclass:: afrc.polymer_models.frc.FreelyRotatingChain
   :members:

   .. automethod:: __init__
