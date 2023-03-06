afrc
==============================
[//]: # (Badges)
[![Documentation Status](https://readthedocs.org/projects/afrc/badge/?version=latest)](https://afrc.readthedocs.io/en/latest/?badge=latest)

![AFRC logo](afrc_logo.png)

## About
#### What?
`afrc` is a Python package that implements an analytical version of the Flory Random Coil (i.e. the AFRC) for polypeptides. By way of an example, if you have a protein sequence, one can calculate a variety of polymeric properties simply by passing in the sequence:


	from afrc import AnalyticalFRC
	
	# creat an AnalyticalFRC object
	P = AnalyticalFRC('MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDI')
	
	## from this object you can calculate various polymeric properties
	## directly without any additional information
	
	# get the ensemble-average radius of gyration
	mean_rg = P.get_mean_radius_of_gyration()
	
	# get the ensemble-average end-to-end distance
	mean_e2e = P.get_mean_end_to_end_distance()
	
	# get the full distribution of the radius of gyration
	[bins, p_rg] = P.get_radius_of_gyration_distribution()
	
	# get the full distribution of the distances between residue 4 and 20
	[bins, p_r] = P.get_interresidue_distance_distribution(4,20)

#### Why?
When studying disordered or unfolded polypeptides we often lack a relevant "reference" frame to calibrate our expectations or results. The AFRC provides a pre-parameterized model that recapitulates a polypeptide if it behaved as an ideal chain (which approximates the behavior in a theta solvent). A theta solvent is a solvent where chain-chain and chain-solvent are equally attractive, canceling out excluded volume effects, giving rise to a chain with a polymer scaling exponent of 0.5. Unlike a real chain in a theta solvent, the AFRC also shows ideal chain behavior, such that intra-residue distance distributions can be analytically calculated without concern for finite size effects.

#### How?
Read the preprint! But the TL/DR is we performed numerical simulations to generate Flory Random Coil ensembles (see Mao et al. 2013), and then parameterized closed-form models to reproduce the end-to-end distance distribution and the radius of gyration distance distribution. From these, a variety of additional parameters can be calculated, all without needing to run any simulations. 

#### Who?
The AFRC was developed by Garrett Ginell, Jhullian Alston, and Alex Holehouse in the [Holehouse lab](https://www.holehouselab.com/). For any questions, please contact Alex.


## Implementation details
This analytical solution is based on the rotational isomeric state approximation of Flory and Volkenstein and parameterized using numerical simulations of residue-specific Flory Random Coil. It provides an interface into sequence-specific polymeric properties (i.e., intra-molecular distances) expected for a given sequence and behaves like a polymer in a true theta solvent. In this way, it provides a convenient reference state through which real simulations or experiments can be normalized against.

## Installation
To install `afrc`:

	pip install afrc 


## Google colab notebook
[Click here to use AFRC (and other polymer models) via the Google colab notebook](https://colab.research.google.com/drive/1WHw8ous7IgcKd2LKYuJLeBTlkdEYoRAk?usp=sharing)

## Quickstart
There is a single user-facing object that is built from the `afrc` package, which is the AnalyticalFRC object. This object gives access to a bunch of additional object functions. As an example

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
If you find a bug or have any feature requests, please submit [an issue here on GitHub](https://github.com/idptools/afrc/issues). Also, feel free to [shoot Alex an email]()

## License
The afrc package is distributed by the [Holehouse Lab](https://www.holehouselab.com/) under the GNU LESSER GENERAL PUBLIC LICENSE.

## References
Flory, P. J. (1969). Statistical Mechanics of Chain Molecules. Oxford University Press.

Volkenstein, M. V. (1958). The configurational statistics of polymeric chains. Journal of Polymer Science, 29(120), 441–454.

Mao, A. H., Lyle, N., & Pappu, R. V. (2013). Describing sequence–ensemble relationships for intrinsically disordered proteins. Biochemical Journal, 449(2), 307–318.


