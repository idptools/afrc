# .....................................................................................
#        
def validate_keyword(viable_keywords, input_keyword, keyword_name):
    """
    General function that protects against poorly defined user input.

    Parameters
    ..........
    
    input_keyword : string
         This defines the input value provided by the user. All keywords are lowercase so this
         function also automatically casts the keyword to lowercase to ensure keywords are actually 
         case insensitive.
        
    viable_keywords : list
         viable_options is a list of strings which are the complete set of options that will be
         considered good. This is provided by the code (i.e. is hard-coded in for a given
         function).

    keyword_name : string
         This is the name of the keyword being defined, and again is hard-coded by the function

    Returns
    ........
    Returns the lower-case cast keyword if valid, else raises an exception

    """

    error_message = "%s must be set to one of %s (was set to %s)" % (keyword_name, viable_keywords, input_keyword)
        
    try:
        input_keyword = input_keyword.lower()
    except:
        raise AFRCException(error_message)

    if input_keyword not in viable_keywords:
        raise AFRCException(error_message)

    return input_keyword
