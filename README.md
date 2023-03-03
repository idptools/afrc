afrc
==============================
[//]: # (Badges)
[![Documentation Status](https://readthedocs.org/projects/afrc/badge/?version=latest)](https://afrc.readthedocs.io/en/latest/?badge=latest)

`afrc` is a Python package that implements an analytical version of the Flory Random Coil (i.e. the AFRC) for polypeptides. 

This analytical solution is based on the rotational isomeric state approximation of Flory and Volkenstein and parameterized using numerical simulations of residue-specific Flory Random Coil. It provides an interface into sequence-specific polymeric properties (i.e. intra-molecular distances) expected for a given sequence and behaves like a polymer in a true theta solvent. In this way it provides a convenient reference state though which real simulations or experiments can be normalized against.

## Installation
To install `afrc`:

	pip install afrc 


## Quickstart
There is a single user-facing object that is built from the `afrc` package which is the AnalyticalFRC object. This object gives access to a bunch of additional object functions. As an example

	from afrc import AnalyticalFRC
	
	A = AnalyticalFRC('APPAPAPAPPAPAPAPPAPPAPPAPAPPA')
	
	# prints the expected radius of gyration if the associated sequence 
	# behaved like a bona fide Flory Random Coil
	print(A.get_mean_radius_of_gyration())  
	
	# prints the expected end-to-end distance if the associated sequence 
	# behaved like a bona fide Flory Random Coil
	print(A.get_mean_end_to_end_distance())
	

## Documentation
For full documentation, [see here](https://afrc.readthedocs.io/)

## Help
If you find a bug or have any feature requests please submit [an issue here on GitHub](https://github.com/idptools/afrc/issues). Also feel free to [shoot Alex an email]()

## License
The afrc package is distributed by the [Holehouse Lab](https://www.holehouselab.com/) under the GNU LESSER GENERAL PUBLIC LICENSE.


