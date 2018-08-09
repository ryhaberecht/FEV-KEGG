from builtins import int
import math
from FEV_KEGG import settings


def getExpectationValue(bitScore: float, searchSequenceLength: int, foundSequenceLength: int, numberOfSequencesInDatabase: int) -> float:
    """
    Returns the E-value of a single query in a sequence database.
    
    Parameters
    ----------
    bitScore : float
        Comparison score, normalised to base 2 (bits).
    searchSequenceLength : int
        Length of the sequence that was the input of the query.
    foundSequenceLength : int
        Length of the sequence that was the result of the query. If you have several results, run this function for each of them individually.
    numberOfSequencesInDatabase : int
        Count of all sequences that could have potentially be found.
        For example, a search in all genes of eco would mean the count of all sequenced genes for eco -> numberOfSequencesInDatabase = 4,498.
        A search in all organisms, however, would mean the count of all sequenced genes in KEGG -> numberOfSequencesInDatabase = 25,632,969.
        
    Returns
    -------
    float
        Statistical E-value (expectation value) for the occurence of a match of the same confidence with a totally unrelated, e.g. random, sequence.
    """
    return numberOfSequencesInDatabase * searchSequenceLength * foundSequenceLength * math.pow(2, -bitScore)


def isMatchSignificant(bitScore: float, searchSequenceLength: int, foundSequenceLength: int, numberOfSequencesInDatabase: int, significanceThreshold: float = settings.defaultEvalue) -> bool:
    """
    Check if a sequence match is significant.
    
    Calculates E-value, using :func:`getExpectationValue`. If E-value is smaller than `significanceThreshold`, match is significant, return True.
    
    Parameters
    ----------
    bitScore : float
        Comparison score, normalised to base 2 (bits).
    searchSequenceLength : int
        Length of the sequence that was the input of the query.
    foundSequenceLength : int
        Length of the sequence that was the result of the query. If you have several results, run this function for each of them individually.
    numberOfSequencesInDatabase : int
        Count of all sequences that could have potentially be found.
        For example, a search in all genes of eco would mean the count of all sequenced genes for eco -> numberOfSequencesInDatabase = 4,498.
        A search in all organisms, however, would mean the count of all sequenced genes in KEGG -> numberOfSequencesInDatabase = 25,632,969.
    significanceThreshold : float, optional
        Threshold of E-value below which to consider a match significant.
    
    Returns
    -------
    bool
        Whether a sequence match is significant.
    """
    if getExpectationValue(bitScore, searchSequenceLength, foundSequenceLength, numberOfSequencesInDatabase) <= significanceThreshold:
        return True
    else:
        return False