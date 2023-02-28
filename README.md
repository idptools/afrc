afrc
==============================
[//]: # (Badges)
[![Documentation Status](https://readthedocs.org/projects/afrc/badge/?version=latest)](https://afrc.readthedocs.io/en/latest/?badge=latest)

`afrc` is a Python package that implements an analytical version of the Flory Random Coil (FRC) for polypeptides. This analytical solution is based on the rotational isomeric state approximation of Flory and Volkenstein and parameterized on the excluded volumed dihedral backbone maps. It provides an interface into sequence-specific polymeric properties (i.e. intra-molecular distances) expected for a given sequence behaves like a polymer in a true theta solvent. In this way it provides a convenient reference state though which real simulations or experiments can be normalized against.


### Installation
To install `afrc` 

1. Download from github as a `.zip` file
2. Run

		pip install afrc-main.zip

That should do it! To install a development version (i.e. so you can edit the code and see changes appear globally in real time) download and unpack and install from within the main source directory using.

	pip install -e .
	

### Usage
There is a single user-facing object that is built from the `afrc` package which is the AnalyticalFRC object. This object gives access to a bunch of additional object functions. As an example

	from afrc import AnalyticalFRC
	
	A = AnalyticalFRC('APPAPAPAPPAPAPAPPAPPAPPAPAPPA')
	
	# prints the expected radius of gyration if the associated sequence 
	# behaved like a bona fide Flory Random Coil
	print(A.get_mean_rg())  
	
	# prints the expected end-to-end distance if the associated sequence 
	# behaved like a bona fide Flory Random Coil
	print(A.get_mean_re())
	

### Documentation
For full documentation, build the docs in `/docs`.

### Copyright

Copyright (c) 2019-2023 Holehouse Lab


#### Acknowledgements
