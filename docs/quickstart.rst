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

   # build the end-to-end distribution. The distribution functions return a
   # 2-pair tuple (distances, probabilities), so they can be unpacked directly:
   distances, probabilities = protein.get_end_to_end_distribution()

   # build the rg distribution (also returned as (distances, probabilities))
   rg_distances, rg_probabilities = protein.get_radius_of_gyration_distribution()

   # compute the average rg
   mean_rg = protein.get_mean_radius_of_gyration()

   # compute the mean end-to-end distance (Re)
   mean_re = protein.get_mean_end_to_end_distance()

   # compute the mean hydrodynamic radius (Rh)
   mean_rh = protein.get_mean_hydrodynamic_radius()

   # build a contact map using a 15 angstrom contact threshold
   contact_map = protein.get_contact_map(15.0)




