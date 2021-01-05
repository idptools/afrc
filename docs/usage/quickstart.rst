AFRC Quickstart
=========================================================

Installation
************************
To install the Analytical Flory Random Coil (afrc) package download the zip file from GitHub.

Then using pip install from the .zip file:

    ``pip install afrc-main.zip``

afrc requires ``numpy``, although this is dealt with automatically through the ``pip`` installation.


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




