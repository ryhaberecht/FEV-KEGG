def getPercent(x, baseValue):
    """
    Calculate percentage.
    
    Parameters
    ----------
    x : float or int
    baseValue : float or int
    
    Returns
    -------
    float
        Percentage of `x` out of `baseValue`, normalised to 100.
    """
    return x/baseValue*100

def getPercentString(x, baseValue):
    """
    Calculate percentage and return as string.
    
    Parameters
    ----------
    x : float or int
    baseValue : float or int
    
    Returns
    -------
    str
        String of result of :func:`getPercent`.
    """
    return str( getPercent(x, baseValue) )

def getPercentStringShort(x, baseValue, decimalPlaces = 1):
    """
    Calculate percentage and return as shortened string.
    
    Parameters
    ----------
    x : float or int
    baseValue : float or int
    decimalPlaces : int, optional
        The number of decimal places to the right to conserve.
    
    Returns
    -------
    str
        String of result of :func:`getPercent`, shortened to `decimalPlaces` decimal places.
    """
    return ("%2." + str( decimalPlaces ) + "f") % getPercent(x, baseValue)

def getPercentSentence(x, baseValue, decimalPlaces = 1):
    """
    Calculate percentage and return shortened string within a fancy sentence.
    
    Parameters
    ----------
    x : float or int
    baseValue : float or int
    decimalPlaces : int, optional
        The number of decimal places to the right to conserve.
    
    Returns
    -------
    str
        String of result of :func:`getPercent`, shortened to `decimalPlaces` decimal places, and put into a complete sentence of the form "`x`/`baseValue` -> result%"", i.e. ``23/42 -> 54.76%``.
    """
    return str(x) + "/" + str(baseValue) + " -> " + getPercentStringShort(x, baseValue, decimalPlaces) + '%'

def floatToPercentString(x, decimalPlaces = 1):
    """
    Format a float as a truncated percent string.
    
    Parameters
    ----------
    x : float
        Must be between 0 and 1 to represent 0% to 100%.
    decimalPlaces : int, optional
        The number of decimal places to the right to conserve.
    
    Returns
    -------
    str
        String of `x`, shortened to `decimalPlaces` decimal places, concatenated with the '%' sign.
    """
    return ("%2." + str( decimalPlaces ) + "f") % (x * 100) + '%'