from builtins import int
import urllib

from FEV_KEGG.lib.Biopython.KEGG import REST
from bs4 import BeautifulSoup
from retrying import retry
import tqdm
from typing import List, Tuple

from FEV_KEGG.Util import Util, Parallelism
from FEV_KEGG.KEGG import SSDB
import FEV_KEGG.settings as settings
import concurrent.futures
import re


def is_not_404(exception):
    """
    Checks if `exception` is **not** an HTTP 404 error.
    
    This is useful, because for most types of downloads 404 means 'does not exist' and is in itself a valid reply. Therefore, retrying a download because it raised a 404 error is usually unnecessary.
    
    Parameters
    ----------
    exception : Exception
        Any error class instance. Must be of class :class:`urllib.error.HTTPError` to be recognised as a 404 error.
    
    Returns
    -------
    bool
        *False*, if `exception` is an HTTP 404 error. *True*, if it is anythign else.
    """
    return not ( isinstance(exception, urllib.error.HTTPError) and exception.code == 404 )


@retry(wait_exponential_multiplier=settings.retryDownloadBackoffFactor, wait_exponential_max=settings.retryDownloadBackoffMax, stop_max_delay=settings.retryDownloadMax, retry_on_exception=is_not_404)
def downloadPathwayList(organismString: 'eco') -> str:
    """
    Downloads list of all pathways for a given organism from KEGG.
    
    Tries several times before giving up, see :attr:`FEV_KEGG.settings.retryDownloadBackoffFactor`.
    
    Parameters
    ----------
    organismString : str
        Abbreviation of the organism, e.g. 'eco'.
    
    Returns
    -------
    str
        List of pathways, delimited by '\\\\n'.
    
    Raises
    ------
    HTTPError
        If pathway list does not exist.
    URLError
        If connection to KEGG fails.
    """
    return REST.kegg_list('pathway', organismString, timeout=settings.downloadTimeoutSocket).read()


@retry(wait_exponential_multiplier=settings.retryDownloadBackoffFactor, wait_exponential_max=settings.retryDownloadBackoffMax, stop_max_delay=settings.retryDownloadMax, retry_on_exception=is_not_404) # do not retry on HTTP error 404, raise immediately instead
def downloadPathway(organismString: 'eco', pathwayName: '00260') -> str:
    """
    Downloads pathway as KGML for a given organism from KEGG.
    
    Tries several times before giving up, see :attr:`FEV_KEGG.settings.retryDownloadBackoffFactor`.
    
    Parameters
    ----------
    organismString : str
        Abbreviation of the organism, e.g. 'eco'.
    pathwayName : str
        Name of the pathway, e.g. '00260'. Will be automatically concatenated with `organismString` to form the pathway ID, e.g. 'eco:00260'.
    
    Returns
    -------
    str
        Pathway in KGML format.
    
    Raises
    ------
    HTTPError
        If pathway does not exist.
    URLError
        If connection to KEGG fails.
    """
    return REST.kegg_get(organismString + pathwayName, 'kgml', timeout=settings.downloadTimeoutSocket).read()
    

@retry(wait_exponential_multiplier=settings.retryDownloadBackoffFactor, wait_exponential_max=settings.retryDownloadBackoffMax, stop_max_delay=settings.retryDownloadMax)
def downloadGene(geneID: 'eco:b0004') -> str:
    """
    Downloads gene description for a given gene ID (includes organism) from KEGG.
    
    Tries several times before giving up, see :attr:`FEV_KEGG.settings.retryDownloadBackoffFactor`.
    
    Parameters
    ----------
    geneID : str
        ID of the gene, including organism abbreviation, e.g. 'eco:b0004'.
    
    Returns
    -------
    str
        Gene in KEGG GENE format.
    
    Raises
    ------
    HTTPError
        If gene does not exist.
    URLError
        If connection to KEGG fails.
    """
    result = REST.kegg_get(geneID, timeout=settings.downloadTimeoutSocket).read()
    if len( result ) < 3:
        raise urllib.error.HTTPError( "Download too small:\n" + result)
    else:
        return result

def downloadGeneBulk(geneIDs: '[eco:b0004, eco:b0015,...]') -> str:
    """
    Downloads gene descriptions for a given list of gene IDs (includes organism) from KEGG
    
    Tries several times before giving up, see :attr:`FEV_KEGG.settings.retryDownloadBackoffFactor`.
    
    Parameters
    ----------
    geneIDs : Iterable
        IDs of the genes, including organism abbreviation, e.g. '[eco:b0004, eco:b0015,...]'.
    
    Returns
    -------
    str
        Genes in KEGG GENE format, delimited by a line of '///'. You will have to split them! Order is arbitrary.
    
    Raises
    ------
    IOError
        If result is too small. Possibly because none of the genes exist.
    URLError
        If connection to KEGG fails.
    """    
    max_query_count = 10 # hard limit imposed by KEGG server
    
    # split list of GeneIDs into chunks of size max_query_count
    geneIDs_chunks = Util.chunks(geneIDs, max_query_count)
    
    # form sub-queries
    query_parts = []
    for chunk in geneIDs_chunks:
        query_parts.append( '+'.join(chunk) )
    
    tqdmPosition = Parallelism.getTqdmPosition()
    threadPool = concurrent.futures.ThreadPoolExecutor(Parallelism.getNumberOfThreadsDownload())
    futures = []
    iterator = None
    
    try:
        # query KEGG in parallel
        
        for query_part in query_parts:
            futures.append( threadPool.submit(_downloadGeneBulk, query_part) )
        
        iterator = concurrent.futures.as_completed(futures)
        
        if settings.verbosity >= 1:
            if settings.verbosity >= 2:
                print( 'Downloading ' + str(len(geneIDs)) + ' genes, max. ' + str(max_query_count) + ' per chunk...' )
            iterator = tqdm.tqdm(iterator, total = len(query_parts), unit = ' *10 genes', position = tqdmPosition)
        
        result = ''
        for future in iterator:
            try:
                result_part = future.result()
            except concurrent.futures.CancelledError:
                continue
            result += result_part
            
        threadPool.shutdown(wait = False)
    
    except KeyboardInterrupt: # only raised in main thread (once in each process!)
        
        Parallelism.keyboardInterruptHandler(threadPool=threadPool, threadPoolFutures=futures, terminateProcess=True)
        raise
    
    except BaseException:
        
        if Parallelism.isMainThread():
            Parallelism.keyboardInterruptHandler(threadPool=threadPool, threadPoolFutures=futures, silent=True)
        raise
        
    finally:
        
        if iterator is not None: iterator.close()
    
    return result

@retry(wait_exponential_multiplier=settings.retryDownloadBackoffFactor, wait_exponential_max=settings.retryDownloadBackoffMax, stop_max_delay=settings.retryDownloadMax)
def _downloadGeneBulk(query_part):
    if Parallelism.getShallCancelThreads() is True:
        raise concurrent.futures.CancelledError()
    else:
        result = REST.kegg_get(query_part, timeout=settings.downloadTimeoutSocket).read()
        if len( result ) < 3:
            raise IOError( "Download too small:\n" + result)
        else:
            return result


@retry(wait_exponential_multiplier=settings.retryDownloadBackoffFactor, wait_exponential_max=settings.retryDownloadBackoffMax, stop_max_delay=settings.retryDownloadMax)
def downloadOrganismList() -> str:
    """
    Download the list of all organisms known to KEGG.
    
    Tries several times before giving up, see :attr:`FEV_KEGG.settings.retryDownloadBackoffFactor`.
    
    Returns
    -------
    str
        List of organism descriptions known to KEGG, delimited by '\\\\n'.
    
    Raises
    ------
    URLError
        If connection to KEGG fails.
    """
    return REST.kegg_list('organism', timeout=settings.downloadTimeoutSocket).read()


@retry(wait_exponential_multiplier=settings.retryDownloadBackoffFactor, wait_exponential_max=settings.retryDownloadBackoffMax, stop_max_delay=settings.retryDownloadMax, retry_on_exception=is_not_404) # do not retry on HTTP error 404, raise immediately instead
def downloadEnzymeEcNumbers(enzymeAbbreviation) -> str:
    """
    Download the list of all EC numbers for a given enzyme, identified by its abbreviation, from KEGG.
    
    Also works for everything else in the description of an enzyme, not just the abbreviation.
    Tries several times before giving up, see :attr:`FEV_KEGG.settings.retryDownloadBackoffFactor`.
    
    Parameters
    ----------
    enzymeAbbreviation : str
        Common abbreviation of the desired enzyme, as it appears in its description, e.g. 'MiA'. Also works for everything else in the description of an enzyme, not just the abbreviation.
    
    Returns
    -------
    str
        EC numbers, delimited by '\\\\n'.
    
    Raises
    ------
    URLError
        If connection to KEGG fails.
    """
    ecNumbers = []
    
    # look up enzyme EC numbers
    searchResult = REST.kegg_find('enzyme', enzymeAbbreviation, timeout=settings.downloadTimeoutSocket).read().split('\n')
    for line in searchResult:
        
        if len( line ) < 10:
            continue
        
        ecNumber = line.split('\t')[0].split(':')[1]
        ecNumbers.append(ecNumber)
    
    return '\n'.join(ecNumbers)


@retry(wait_exponential_multiplier=settings.retryDownloadBackoffFactor, wait_exponential_max=settings.retryDownloadBackoffMax, stop_max_delay=settings.retryDownloadMax, retry_on_exception=is_not_404) # do not retry on HTTP error 404, raise immediately instead
def downloadOrganismInfo(organismAbbreviation) -> str:
    """
    Downloads the info file of an organism.
    
    Parameters
    ----------
    organismAbbreviation : str
        Abbreviation of the organism to check, e.g. 'eco'.
    
    Returns
    -------
    str
        Raw organism info. *None*, if download was empty (400 Bad Request), because this organism does not exist.
        
    Raises
    ------
    URLError
        If connection to KEGG fails.
    """
    try:
        return REST.kegg_info(organismAbbreviation, timeout=settings.downloadTimeoutSocket).read()
    except urllib.error.HTTPError as e:
        if isinstance(e, urllib.error.HTTPError) and e.code == 400:
            return None
        else:
            raise




def downloadOrthologs(geneID: 'GeneID', comparisonOrganismString: 'eco') -> Tuple[int, List[SSDB.PreMatch]]:
    """
    Download orthologs of gene `geneID` found in organism `comparisonOrganismString`.
    
    Parameters
    ----------
    geneID : GeneID
        GeneID object of the gene to be compared against, i.e. against its amino acid sequence.
    comparisonOrganismString : str
        Abbreviation of the organism to search for orthologs of `geneID`.
    
    Returns
    -------
    List[SSDB.PreMatch]
        List of Pre-Matches, containing gene IDs of orthologs, and other data necessary for sequence matching. Will be empty, if nothing is found.
    
    Raises
    ------
    URLError
        If connection to KEGG fails.
    """
    # download list of orthologs
    data = _downloadHomologs(geneID.geneIDString, comparisonOrganismString)
    
    # parse HTML
    orthologData = _parseSsdbOrthologView(data)
    if orthologData is None:
        print( 'Failed to get SSDB data for ' + str(geneID) + ' in ' + comparisonOrganismString)
        return None
    
    foundGenes = orthologData
    
    return foundGenes
    

def downloadParalogs(geneID: 'GeneID') -> Tuple[int, List[SSDB.PreMatch]]:
    """
    Download paralogs of gene `geneID`.
    
    Parameters
    ----------
    geneID : GeneID
        GeneID object of the gene to be compared against, i.e. against its amino acid sequence.
    
    Returns
    -------
    List[SSDB.PreMatch]
        List of Pre-Matches, containing gene IDs of paralogs, and other data necessary for sequence matching. Will be empty, if nothing is found.
    
    Raises
    ------
    URLError
        If connection to KEGG fails.
    """
    # download list of paralogs
    data = _downloadHomologs(geneID.geneIDString, geneID.organismAbbreviation)
    
    # parse HTML
    orthologData = _parseSsdbOrthologView(data)
    if orthologData is None:
        print( 'Failed to get SSDB data for ' + str(geneID))
        return None
    
    searchedSequenceLength, foundGenes = orthologData 
    
    # remove the gene that was searched for
    for preMatch in foundGenes:
        if preMatch.foundGeneIdString == geneID.geneIDString:
            foundGenes.remove(preMatch)
            break
    
    return (searchedSequenceLength, foundGenes)

@retry(wait_exponential_multiplier=settings.retryDownloadBackoffFactor, wait_exponential_max=settings.retryDownloadBackoffMax, stop_max_delay=settings.retryDownloadMax, retry_on_exception=is_not_404) # do not retry on HTTP error 404, raise immediately instead
def _downloadHomologs(geneIdString, organismAbbreviationString):
    return str(urllib.request.urlopen('https://www.kegg.jp/ssdb-bin/ssdb_ortholog_view?org_gene=' + geneIdString + '&org=' + organismAbbreviationString, timeout=settings.downloadTimeoutSocket).read()).replace('\\n', '')

AA_SEQ_LENGTH_REGEX_PATTERN = re.compile('\(([0-9]+) a\.a\.\)')
NT_SEQ_LENGTH_REGEX_PATTERN = re.compile('\(([0-9]+) n\.t\.\)') # length in AA == length in NT / 3 - 1

def _parseSsdbOrthologView(htmlString) -> Tuple[int, List[SSDB.PreMatch]]:
    
    try:
        html = BeautifulSoup(htmlString, 'html.parser')
        searchedSequenceLengthMatch = AA_SEQ_LENGTH_REGEX_PATTERN.search(html.body.a.next_sibling)
        if searchedSequenceLengthMatch is None: # length in amino acids not found, maybe it is given in nucleic acids
            searchedSequenceLengthMatch = NT_SEQ_LENGTH_REGEX_PATTERN.search(html.body.a.next_sibling)
        searchedSequenceLength = int(searchedSequenceLengthMatch.group(1))
        
        matches = []
        
        for index, tr in enumerate( html.table.children ):
            
            # ignore head of table
            if index == 0:
                continue
            
            for index2, td in enumerate( tr.children ):
                 
                if index2 == 0: # read gene ID
                    foundGeneIdString = td.text
                     
                elif index2 == 1: # read Smith-Waterman score
                    swScore = int(td.text)
                     
                elif index2 == 2: # read bit score
                    bitScore = float(td.text)
                     
                elif index2 == 3: # read identity
                    identity = float(td.text)
                     
                elif index2 == 4: # read overlap
                    overlap = int(td.text)
            
            matches.append( SSDB.PreMatch(foundGeneIdString, swScore, bitScore, identity, overlap) )
        
        return (searchedSequenceLength, matches)
    
    except BaseException as _:
        return None




def downloadOrthologOverview(geneID: 'GeneID') -> Tuple[int, List[SSDB.Match]]:
    """
    Download overview of orthologs of gene `geneID` found in any organism.
    
    Parameters
    ----------
    geneID : GeneID
        GeneID object of the gene to be compared against, i.e. against its amino acid sequence.
    
    Returns
    -------
    List[SSDB.Match]
        List of Matches, containing gene IDs of orthologs, and other data necessary for sequence matching. Will be empty, if nothing is found.
    
    Raises
    ------
    URLError
        If connection to KEGG fails.
    """
    # download list of orthologs
    data = _downloadOrthologOverview(geneID.geneIDString)
    
    # parse HTML
    try:
        foundGenes = _parseSsdbBestView(data)
    except (IndexError, ValueError, AttributeError):
        print("\nError when parsing ortholog overview of gene: " + str(geneID))
        return None
    
    return foundGenes

@retry(wait_exponential_multiplier=settings.retryDownloadBackoffFactor, wait_exponential_max=settings.retryDownloadBackoffMax, stop_max_delay=settings.retryDownloadMax, retry_on_exception=is_not_404) # do not retry on HTTP error 404, raise immediately instead
def _downloadOrthologOverview(geneIdString):
    return str(urllib.request.urlopen('https://www.kegg.jp/ssdb-bin/ssdb_best_best?threshold=400&org_gene=' + geneIdString, timeout=settings.downloadTimeoutSocket).read()).replace('\\n', '')

SSDB_OVERVIEW_REGEX = re.compile("\)\s*|\s*[\(]{0,1}\s*")

def _parseSsdbBestView(htmlString) -> Tuple[int, List[SSDB.Match]]:
    
    html = BeautifulSoup(htmlString.replace('&#', ''), 'html.parser')
    searchedSequenceLengthMatch = AA_SEQ_LENGTH_REGEX_PATTERN.search(html.table.tr.text)
    if searchedSequenceLengthMatch is None: # length in amino acids not found, maybe it is given in nucleic acids
        searchedSequenceLengthMatch = NT_SEQ_LENGTH_REGEX_PATTERN.search(html.table.tr.text)
    searchedSequenceLength = int(searchedSequenceLengthMatch.group(1))
    
    matches = []
    
    content = str(html.pre).split('<input')
    content.pop(0) # ignore first line, containing table heading
    
    for line in content:
        line = '<input' + line
        #print(line)
        lineHtml = BeautifulSoup(line, 'html.parser')
        #print(lineHtml)
        foundGeneIdString = lineHtml.input['value']
        
#         textFields = lineHtml.text.split()#.replace('&', '').replace('#', '')
        textFields = SSDB_OVERVIEW_REGEX.split(lineHtml.text)
#         print(textFields)
        
        try:
            length = int(textFields[-8])
            swScore = int(textFields[-7])
            bitScore = float(textFields[-5])
            identity = float(textFields[-4])
            overlap= int(textFields[-3])
        except (IndexError, ValueError):
            print(line)
            raise
        
        try:
            match = SSDB.Match(foundGeneIdString, swScore, bitScore, identity, overlap, length)
            matches.append( match )
        except ValueError: # ignore GeneID which is not correctly formatted. Can happen when it is an incomplete "Addendum" organism, or a virus, etc.
            continue        
    
    return (searchedSequenceLength, matches)
    






@retry(wait_exponential_multiplier=settings.retryDownloadBackoffFactor, wait_exponential_max=settings.retryDownloadBackoffMax, stop_max_delay=settings.retryDownloadMax)
def downloadTaxonomyNCBI() -> str:
    """
    Download NCBI taxonomy from KEGG BRITE.
    
    Returns
    -------
    str
        NCBI taxonomy in special text format.
    
    Raises
    ------
    URLError
        If connection to KEGG fails.
    """
    return REST.kegg_get('br:br08610', timeout=settings.downloadTimeoutSocket).read()


@retry(wait_exponential_multiplier=settings.retryDownloadBackoffFactor, wait_exponential_max=settings.retryDownloadBackoffMax, stop_max_delay=settings.retryDownloadMax)
def downloadTaxonomyKEGG():
    """
    Download KEGG taxonomy from KEGG BRITE.
    
    Returns
    -------
    str
        KEGG taxonomy in special text format.
    
    Raises
    ------
    URLError
        If connection to KEGG fails.
    """
    return REST.kegg_get('br:br08601', timeout=settings.downloadTimeoutSocket).read()



@retry(wait_exponential_multiplier=settings.retryDownloadBackoffFactor, wait_exponential_max=settings.retryDownloadBackoffMax, stop_max_delay=settings.retryDownloadMax, retry_on_exception=is_not_404) # do not retry on HTTP error 404, raise immediately instead
def downloadSubstance(substanceID):
    """
    Download a substance description file from KEGG, compound or glycan.
    
    Parameters
    ----------
    substanceID : str
        The ID string of the substance to download, e.g. 'C00084'.
    
    Returns
    -------
    str
        Content of the substance's description.
    
    Raises
    ------
    URLError
        If connection to KEGG fails.
    """
    return REST.kegg_get(substanceID, timeout=settings.downloadTimeoutSocket).read()


@retry(wait_exponential_multiplier=settings.retryDownloadBackoffFactor, wait_exponential_max=settings.retryDownloadBackoffMax, stop_max_delay=settings.retryDownloadMax, retry_on_exception=is_not_404) # do not retry on HTTP error 404, raise immediately instead
def downloadEcEnzyme(ecNumberID):
    """
    Download an enzyme description file from KEGG, defined by its EC number.
    
    Parameters
    ----------
    ecNumber : str
        The EC number string of the enzyme to download, e.g. '4.1.2.48'.
    
    Returns
    -------
    str
        Content of the EcEnzymes's description.
    
    Raises
    ------
    URLError
        If connection to KEGG fails.
    """
    return REST.kegg_get('ec:' + ecNumberID, timeout=settings.downloadTimeoutSocket).read()

    