from typing import List, Tuple, Dict
def deduplicateList(anyList: List, preserveOrder = False):
    """
    Deduplicates a list.
    
    Does not preserve order by default, because it is faster.
    
    Parameters
    ----------
    anyList : list
        The list to be deduplicated.
    preserveOrder : bool, optional
        If *True*, preserves order of elements in the list.
    
    Returns
    -------
    List
        A new list, containing all elements of `anyList`, but only once. 
    """
    if preserveOrder == False:
        return _deduplicateListNonPreserving(anyList)
    else:
        return _deduplicateListPreserving(anyList)
    

def _deduplicateListNonPreserving(anyList: List):
    keys = {}
    for e in anyList:
        keys[e] = 1
    return list( keys.keys() )


def _deduplicateListPreserving(seq: List, idfun=None):
    # order preserving
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        # in old Python versions:
        # if seen.has_key(marker)
        # but in new ones:
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    return result


def chunks(iterable, chunk_size):
    """
    Chops Iterable into chunks.
    
    Parameters
    ----------
    iterable : Iterable
        The Iterable object to be chunked.
    chunk_size : int
        The size if each chunk. Except for the last, of course.
    
    Yields
    -------
    Iterable
        Iterable of the same type as `iterable`, but with length `chunk_size`. Except for the last, of course.
    """
    iterable = list(iterable)
    for i in range(0, len(iterable), chunk_size):
        yield iterable[i:i + chunk_size]

def prettyPrintDict(dictionary, byValueFirst = False) -> List[Tuple]:
    
    sortedValues = dict()
    for key, value in dictionary.items():
        sortedValues[key] = sorted(value)
    
    if byValueFirst:
        anonymousFunctionA = lambda item: item[0]
        anonymousFunctionB = lambda item: item[1]
                
    else:
        anonymousFunctionA = lambda item: item[1]
        anonymousFunctionB = lambda item: item[0]
    
    preSortedDict = sorted(sortedValues.items(), key=anonymousFunctionA)
    sortedDict = sorted(preSortedDict, key=anonymousFunctionB)
    
    return sortedDict

def inverseDictKeepingAllKeys(dictionary) -> Dict:
    
    inversed = {}
    for key, value in dictionary.items():
        currentSet = inversed.get(value, None)
        
        if currentSet is None:
            currentSet = set()
            inversed[value] = currentSet
        
        currentSet.add( key )
    
    return inversed
    