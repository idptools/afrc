"""
iofunctions.py

A collection of utilities for working with input/output data in the afrc module.

Copyright Alex Holehouse 2018-2026 (holehouselab.com).

For any questions please contact Alex.

"""

from .exceptions import AFRCException


# .....................................................................................
#
def validate_keyword(viable_keywords, input_keyword, keyword_name):
    """
    Validate a user-supplied keyword against a fixed set of options.

    Keywords are case insensitive, so the input is lower-cased before it is
    checked.

    Parameters
    ----------
    viable_keywords : list of str
        The complete set of (lower-case) options that are accepted. This is
        hard-coded by the calling function.

    input_keyword : str
        The value provided by the user.

    keyword_name : str
        The name of the argument being validated, used in the error message.

    Returns
    -------
    str
        The lower-cased keyword.

    Raises
    ------
    AFRCException
        If ``input_keyword`` is not a string or is not one of ``viable_keywords``.

    """

    # build a custom error message
    error_message = f"{keyword_name} must be set to one of {viable_keywords} (was set to {input_keyword})"

    # first see if you can cast the keyword to lower (if this fails the input is probably not
    # even a string, but we use the same error message
    try:
        input_keyword = input_keyword.lower()
    except AttributeError:
        raise AFRCException(error_message)

    # next check if the input keyword was one of the allowed words, and, if not, we raise an exception
    if input_keyword not in viable_keywords:
        raise AFRCException(error_message)

    # if everything was ok, just return the lower() version of the input keyword
    return input_keyword
