"""
iofunctions.py

A collection of utilities for working with input/output data in the afrc module.

Copyright Alex Holehouse 2018-2022 (holehouselab.com).

For any questions please contact Alex.

"""


# .....................................................................................
#        
def validate_keyword(viable_keywords, input_keyword, keyword_name):
    """
    General function that protects against poorly defined user input. This ensures
    that input keywords are case insensitive (which in general we want).

    Parameters
    ..........
    
    input_keyword : str
         This defines the input value provided by the user. All keywords are 
         lowercase so this function also automatically casts the keyword to 
         lowercase to ensure keywords are actually case insensitive.
        
    viable_keywords : list
         viable_options is a list of strings which are the complete set of 
         options that will be considered good. This is provided by the code 
         (i.e. is hard-coded in for a given function).

    keyword_name : string
         This is the name of the keyword being defined, and again is 
         hard-coded by the function.

    Returns
    ........
        Returns the lower-case cast keyword if valid, else raises an 
        exception.

    """

    # build a custom error message
    error_message = f"{keyword_name} must be set to one of {viable_keywords} (was set to {input_keyword})" 

    # first see if you can cast the keyword to lower (if this fails the input is probably not
    # even a string, but we use the same error message
    try:
        input_keyword = input_keyword.lower()
    except:
        raise AFRCException(error_message)

    # next check if the input keyword was one of te allowed words, and, if not, we raise an exception
    if input_keyword not in viable_keywords:
        raise AFRCException(error_message)

    # if everything was ok, just return the lower() version of the input keyword
    return input_keyword
