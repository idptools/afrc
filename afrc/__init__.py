"""
AFRC
An analytical version of the Flory Random Coil for polypeptides, implemented using the rotational isomeric state 
approximation of Flory and Volkenstein and parameterized on the excluded volumed dihedral backbone maps.

Copyright Alex Holehouse 2018-2022 (holehouselab.com).

For any questions please contact Alex.

"""

# Add imports here
from .afrc import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
