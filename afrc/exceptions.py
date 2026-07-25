"""
exceptions.py

Exception types used across the afrc package.

These live in their own module (rather than in afrc.py) so that any module in the
package can raise them without creating a circular import - afrc.py imports both
polymer.py and iofunctions.py at load time, so neither of those can import from
afrc.py directly.

Copyright Alex Holehouse 2018-2026 (holehouselab.com).

"""


class AFRCException(Exception):
    """
    Exception class specific for the Analytical FRC.

    """
    pass
