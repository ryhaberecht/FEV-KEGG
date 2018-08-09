from FEV_KEGG.Graph.Elements import Enzyme, GeneID
from typing import Set, Iterable, Dict
from FEV_KEGG.KEGG.Database import getOrthologsBulk
from FEV_KEGG import settings

def getOrthologsWithinGeneIDs(geneIDs: Iterable[GeneID], eValue : float = settings.defaultEvalue) -> Dict[GeneID, Set[GeneID]]:
    """
    Get orthologs within `geneIDs`.
    
    Parameters
    ----------
    geneIDs : Iterable[GeneID]
        Gene IDs among which to search for orthology.
    eValue : float, optional
        Statistical expectation value (E-value), below which a sequence alignment is considered significant.
    
    Returns
    -------
    Dict[GeneID, Set[GeneID]]
        Dictionary of every gene in `geneIDs` which has at least one ortholog in `geneIDs`, pointing to a set of all its orthologs.
    
    Raises
    ------
    ValueError
        If any organism does not exist.
    URLError
        If connection to KEGG fails.
    
    Warnings
    --------
    This operation has a worst-case complexity class of O(n(n+1)/2)! To make matters worse, each step involves at least one, if not several, network operations with KEGG SSDB/GENE, which are inherently slow.
    The best-case complexity class is O(n-1).
    """
    # sort enzymes by their organism
    genesByOrganism = dict()
    for geneID in geneIDs:
        organismAbbreviation = geneID.organismAbbreviation
        organismSet = genesByOrganism.get(organismAbbreviation)
        if organismSet is None:
            genesByOrganism[organismAbbreviation] = set()
        genesByOrganism[organismAbbreviation].add(geneID)
    
    orthologousGenes = dict()
    
    # search orthologs, for each organism's genes in every other organism
    partiallySearchedOrganisms = list(genesByOrganism.keys())
    for organismAbbreviation in genesByOrganism.keys():
        # get all other organisms
        partiallySearchedOrganisms.remove(organismAbbreviation)
        if len(partiallySearchedOrganisms) == 0:
            break
        
        # for each other organism            
        # get orthologs for all of current organism's genes
        currentOrganism_Genes = genesByOrganism[organismAbbreviation]
        matchingsDict = getOrthologsBulk(currentOrganism_Genes, partiallySearchedOrganisms, eValue)
        
        # link found orthologous genes to searched genes
        for searchedGeneID, matchingList in matchingsDict.items():
            for matching in matchingList:
                matchedGeneIDs = {x.foundGeneID for x in matching.matches}
                matchedGeneIDs.intersection_update(geneIDs) # leave only genes we actually search for
                
                if len(matchedGeneIDs) > 0:
                    orthologousGenesSet = orthologousGenes.get(searchedGeneID)
                    if orthologousGenesSet is None:
                        orthologousGenes[searchedGeneID] = set()
                    orthologousGenes[searchedGeneID].update(matchedGeneIDs)
        
    return orthologousGenes

def getOrthologsWithinEnzymes(enzymes: Iterable[Enzyme], eValue : float = settings.defaultEvalue) -> Dict[Enzyme, Set[Enzyme]]:
    """
    Get orthologs within `enzymes`.
    
    Parameters
    ----------
    enzymes : Iterable[Enzyme]
        Enzymes among which to search for orthology.
    eValue : float, optional
        Statistical expectation value (E-value), below which a sequence alignment is considered significant.
    
    Returns
    -------
    Dict[Enzyme, Set[Enzyme]]
        Dictionary of every enzyme in `enzymes` which has at least one ortholog in `enzymes`, pointing to a set of all its orthologs.
    
    Raises
    ------
    ImpossiblyOrthologousError
        If any gene ID in `geneIDs` is from `comparisonOrganism`.
    ValueError
        If any organism does not exist.
    URLError
        If connection to KEGG fails.
    
    Warnings
    --------
    This operation has a worst-case complexity class of O(n(n+1)/2)! To make matters worse, each step involves at least one, if not several, network operations with KEGG SSDB/GENE, which are inherently slow.
    The best-case complexity class is still O(n-1).
    """
    # create reverse mapping for GeneID -> Enzyme
    gene2Enzyme = dict()
    for enzyme in enzymes:
        gene2Enzyme[enzyme.geneID] = enzyme
    
    orthologousGenesDict = getOrthologsWithinGeneIDs([x.geneID for x in enzymes], eValue)
    
    # reverse-map found GeneID -> Enzyme
    orthologousEnzymes = dict()
    for geneID, orthologousGeneIDs in orthologousGenesDict.items():
        orthologousEnzymes[gene2Enzyme[geneID]] = {gene2Enzyme[geneID] for geneID in orthologousGeneIDs}
    
    return orthologousEnzymes
