AFRC Quickstart
=========================================================

Installation
************************
To install the Analytical Flory Random Coil (AFRC) package download the zip file, unzip, and then using pip install from the wheels file:

    ``pip install dist/afrc-0.0.0+1.g537e9d4-py3-none-any.whl``

This will install AFRC in a system-wide manner using the version of ``pip`` specified. We recommend using ``pip`` embedded within a ``conda`` environment to avoid any possible dependency issues. If you are unfamiliar with ``conda`` or ``pip`` we highly recommend reading this `this introductory material <http://geohackweek.github.io/Introductory/01-conda-tutorial//>`_.

AFRC requires ``numpy``, although this is dealt with automatically through the ``pip`` installation.

AFRC was developed for Python 3, but should be compatiable with Python 2. 


Usage
************************

AFRC gives you a way to obtain a variety of inter-residue distances from an analytical version of the Flory Random Coil. The input required is an amino acid sequence, and from this distribution and mean values are available through what is essentially an API.

As an example


.. code-block:: python

   import afrc

   # create an AFRC object
   protein = afrc.AnalyticalFRC('MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQG')

   # build the internal scaling profile 
   internal_scaling = protein.get_internal_scaling()

   # build the end-to-end distribution
   end_to_end_distribution = protein.get_re_distribution()

   # build the rg distribution
   rg_distribution = protein.get_rg_distribution()

   # compute the average rg
   mean_rg = protein.get_mean_rg()

   # compute the mean end-to-end distance (Re)
   mean_re = protein.get_mean_re()




