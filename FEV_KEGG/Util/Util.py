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


def prettySortDict(dictionary, byValueFirst = False) -> List[Tuple]:
    
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
    sortedList = sorted(preSortedDict, key=anonymousFunctionB)
    
    return sortedList


def inverseDictKeepingAllKeys(dictionary) -> Dict:
    
    inversed = {}
    for key, value in dictionary.items():
        currentSet = inversed.get(value, None)
        
        if currentSet is None:
            currentSet = set()
            inversed[value] = currentSet
        
        currentSet.add( key )
    
    return inversed


def updateDictUpdatingValue(dictA, dictB):
    """
    Update `dictA` using `dictB`. However, if key already exists in dictA, does not overwrite dictA[key] with dictB[key], as the defautl update() function does, instead does an update: dictA[key].update( dictB[key] ).
    
    Warnings
    --------
    Only works if dictA's values have an update() function, e.g. are sets!
    """
    for key, value in dictB.items():
        
        if dictA.get(key, None) is None: # does not exist in A, yet. Copy!
            dictA[key] = value
            
        else: # already exists in A. Update!
            dictA[key].update(value)
    
    return dictA


def dictToHtml(dictionary, byValueFirst = False, addEcDescriptions = False, headingDescriptionForHeading = None) -> str:
    
    if addEcDescriptions is not False:        
        from FEV_KEGG.Graph.Elements import EcNumber
        EcNumber.addEcDescriptions(addEcDescriptions)
    
    sortedList = prettySortDict(dictionary, byValueFirst)
    
    from yattag import Doc
    doc, tag, _ = Doc().tagtext()
    
    with tag('html'):
        with tag('body'):
            
            for heading, rows in sortedList:
                
                with tag('p'):
                    with tag('table'):
                        doc.asis(heading.toHtml())
                        
                    if headingDescriptionForHeading is not None:
                        headingDescription = headingDescriptionForHeading.get(heading, None)
                        
                        if headingDescription is not None:
                            with tag('table'):
                                with tag('tr'):
                                    with tag('td'):
                                        doc.asis('{')
                                
                                for line in headingDescription:
                                    with tag('tr'):
                                        doc.asis(line.toHtml(short = True))
                                        
                                with tag('tr'):
                                    with tag('td'):
                                        doc.asis('}')
                    
                    with tag('table'):
                        
                        for row in rows:
                            with tag('tr'):
                                doc.asis(row.toHtml(short = True))
                                
                with tag('br'):
                    pass
                            
    return doc.getvalue()
    
    
def dictToHtmlFile(dictionary, file, byValueFirst = False, inCacheFolder = False, addEcDescriptions = False, headingDescriptionForHeading = None):
    """
    Parameters
    ----------
    file : str
        Path and name of the exported file. See `inCacheFolder`.
    inCacheFolder : bool, optional
        If *True*, interpret `file` relative to the cache folder. See :attr:`FEV_KEGG.settings.cachePath`.
        If *False*, interpret `file` relative to the current working directory.
    """
    import os
    from FEV_KEGG import settings
    from FEV_KEGG.KEGG import File
    
    htmlString = dictToHtml(dictionary, byValueFirst, addEcDescriptions, headingDescriptionForHeading)
    
    if not file.endswith('.html'):
        file += '.html'
    
    if inCacheFolder is True:
        file = os.path.join(settings.cachePath, file)
    
    dirName = os.path.dirname(file)
    if not os.path.isdir(dirName) and dirName != '':
        os.makedirs(os.path.dirname(file))
    
    File.writeToFile(htmlString, file)
        