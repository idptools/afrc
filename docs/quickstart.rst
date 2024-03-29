AFRC Quickstart
=========================================================
The code below provides some initial examples of how to use the AFRC. For more examples see the demo directory at https://github.com/idptools/afrc/

Usage
************************

AFRC gives you a way to obtain a variety of inter-residue distances from an analytical version of the Flory Random Coil. The input required is an amino acid sequence, and from this distribution and mean values are available through what is essentially an API.

As an example


.. code-block:: python

   from afrc import AnalyticalFRC

   # create an AFRC object
   protein = AnalyticalFRC('MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQG')

   # build the internal scaling profile 
   internal_scaling = protein.get_internal_scaling()

   # build the end-to-end distribution
   end_to_end_distribution = protein.get_end_to_end_distribution()

   # build the rg distribution
   rg_distribution = protein.get_radius_of_gyration_distribution()

   # compute the average rg
   mean_rg = protein.get_mean_radius_of_gyration()

   # compute the mean end-to-end distance (Re)
   mean_re = protein.get_mean_end_to_end_distance()




