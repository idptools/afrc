"""
AFRC
An analytical version of the Flory Random Coil for polypeptides, implemented using the rotational isomeric state 
approximation of Flory and Volkenstein and parameterized on the excluded volumed dihedral backbone maps.

Copyright Alex Holehouse 2018-2022 (holehouselab.com).

For any questions please contact Alex.

"""

# Add imports here
from afrc.afrc import *  # noqa: F401,F403
import os

# Generate _version.py if missing and in the Read the Docs environment
if os.getenv("READTHEDOCS") == "True" and not os.path.isfile('../afrc/_version.py'):
    import versioningit
    __version__ = versioningit.get_version('../')
else:
    from ._version import __version__  # noqa: F401
