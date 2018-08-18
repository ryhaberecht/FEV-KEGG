from typing import Set, Dict, Tuple

from FEV_KEGG.Graph.Elements import Enzyme, EcNumber
from FEV_KEGG.Graph.SubstanceGraphs import SubstanceEnzymeGraph, SubstanceEcGraph
from FEV_KEGG.KEGG import Database
from builtins import set
from FEV_KEGG.KEGG.Organism import Group
from _collections_abc import Iterable
from FEV_KEGG import settings
import itertools
from FEV_KEGG.Drawing import Export

defaultEValue = settings.defaultEvalue
"""
Default threshold for the statistical expectation value (E-value), below which a sequence alignment is considered significant.
This can be overridden in each relevant method's `eValue` parameter in this module.
"""


class GeneFunctionConservation(object):
    """
    Evolutionary event of conserving a gene function (EC number) between a pair of arbitrary ancestor and descendant.
    
    The conditions for a gene function conservation are:
        - The EC number has been conserved, along the way from an older group of organisms to a newer one.
    """
    @staticmethod
    def getECs(ancestorEcGraph: SubstanceEcGraph, descendantEcGraph: SubstanceEcGraph) -> Set[EcNumber]:
        """
        Get EC numbers which have been conserved between ancestor and descendant, existing in both.
        
        Parameters
        ----------
        ancestorEcGraph : SubstanceEcGraph
        descendantEcGraph : SubstanceEcGraph
        
        Returns
        -------
        Set[EcNumber]
            Set of EC numbers which occur in the ancestor's EC graph and in the decendants's, i.e. EC numbers which are conserved in the descendant.
        """
        conservedECs = ancestorEcGraph.getECs()
        conservedECs.intersection_update(descendantEcGraph.getECs())
        return conservedECs
    
    @staticmethod
    def getGraph(ancestorEcGraph: SubstanceEcGraph, descendantEcGraph: SubstanceEcGraph) -> SubstanceEcGraph:
        """
        Get graph containing EC numbers which have been conserved between ancestor and descendant, existing in both.
        
        Parameters
        ----------
        ancestorEcGraph : SubstanceEcGraph
        descendantEcGraph : SubstanceEcGraph
        
        Returns
        -------
        SubstanceEcGraph
            Graph of EC numbers which occur in the ancestor's EC graph, and in the decendants's, i.e. EC numbers which are conserved in the descendant.
            Substance-EC-product edges are only included if both graphs, ancestor and descendant, have both nodes, substrate and product.
        """
        conservedGraph = ancestorEcGraph.intersection(descendantEcGraph, addCount=False)
        conservedGraph.removeIsolatedNodes()
        return conservedGraph
    
    
class GeneFunctionAddition(object):
    """
    Evolutionary event of adding a gene function (EC number) between a pair of arbitrary ancestor and descendant.
    
    The conditions for a gene function addition are:
        - The EC number has been added, from an unknown origin, along the way from an older group of organisms to a newer one.
    """
    @staticmethod
    def getECs(ancestorEcGraph: SubstanceEcGraph, descendantEcGraph: SubstanceEcGraph) -> Set[EcNumber]:
        """
        Get EC numbers which have been added between ancestor and descendant, existing only in the descendant.
        
        Parameters
        ----------
        ancestorEcGraph : SubstanceEcGraph
        descendantEcGraph : SubstanceEcGraph
        
        Returns
        -------
        Set[EcNumber]
            Set of EC numbers which occur in the descendant's EC graph, but not in the ancestor's, i.e. EC numbers which are new to the descendant.
        """
        addedECs = descendantEcGraph.getECs()
        addedECs.difference_update(ancestorEcGraph.getECs())
        return addedECs
    
    @staticmethod
    def getGraph(ancestorEcGraph: SubstanceEcGraph, descendantEcGraph: SubstanceEcGraph) -> SubstanceEcGraph:
        """
        Get graph containing EC numbers which have been added between ancestor and descendant, existing only in the descendant.
        
        Parameters
        ----------
        ancestorEcGraph : SubstanceEcGraph
        descendantEcGraph : SubstanceEcGraph
        
        Returns
        -------
        SubstanceEcGraph
            Graph of EC numbers which occur in the descendant's EC graph, but not in the ancestor's, i.e. EC numbers which are new in the descendant.
            Substance-EC-product edges are only included if both graphs, ancestor and descendant, have both nodes, substrate and product.
        """
        addedGraph = descendantEcGraph.difference(ancestorEcGraph, subtractNodes=False)
        addedGraph.removeIsolatedNodes()
        return addedGraph


class GeneFunctionLoss(object):
    """
    Evolutionary event of losing a gene function (EC number) between a pair of arbitrary ancestor and descendant.
    
    The conditions for a gene function loss are:
        - The EC number has been lost, along the way from an older group of organisms to a newer one.
    """
    @staticmethod
    def getECs(ancestorEcGraph: SubstanceEcGraph, descendantEcGraph: SubstanceEcGraph) -> Set[EcNumber]:
        """
        Get EC numbers which have been lost between ancestor and descendant, existing only in the ancestor.
        
        Parameters
        ----------
        ancestorEcGraph : SubstanceEcGraph
        descendantEcGraph : SubstanceEcGraph
        
        Returns
        -------
        Set[EcNumber]
            Set of EC numbers which occur in the ancestor's EC graph, but not in the decendants's, i.e. EC numbers which are lost to the descendant.
        """
        lostECs = ancestorEcGraph.getECs()
        lostECs.difference_update(descendantEcGraph.getECs())
        return lostECs
    
    @staticmethod
    def getGraph(ancestorEcGraph: SubstanceEcGraph, descendantEcGraph: SubstanceEcGraph) -> SubstanceEcGraph:
        """
        Get graph containing EC numbers which have been lost between ancestor and descendant, existing only in the ancestor.
        
        Parameters
        ----------
        ancestorEcGraph : SubstanceEcGraph
        descendantEcGraph : SubstanceEcGraph
        
        Returns
        -------
        SubstanceEcGraph
            Graph of EC numbers which occur in the ancestor's EC graph, but not in the decendants's, i.e. EC numbers which are lost in the descendant.
            Substance-EC-product edges are only included if both graphs, ancestor and descendant, have both nodes, substrate and product.
        """
        lostGraph = ancestorEcGraph.difference(descendantEcGraph, subtractNodes=False)
        lostGraph.removeIsolatedNodes()
        return lostGraph


class GeneFunctionDivergence(object):
    """
    Evolutionary event of diverging (adding or losing) a gene function (EC number) between a pair of arbitrary ancestor and descendant.
    
    The conditions for a gene function divergence are:
        - The EC number exists in an older group of organisms, but not in a newer one, or the other way around.
    """
    @staticmethod
    def getECs(ancestorEcGraph: SubstanceEcGraph, descendantEcGraph: SubstanceEcGraph) -> Set[EcNumber]:
        """
        Get EC numbers which have diverged between ancestor and descendant, existing only in either one of them.
        
        Obviously, `ancestorEcGraph` and `descendantEcGraph` can be swapped here without changing the result.
        
        Parameters
        ----------
        ancestorEcGraph : SubstanceEcGraph
        descendantEcGraph : SubstanceEcGraph
        
        Returns
        -------
        Set[EcNumber]
            Set of EC numbers which occur in the ancestor's EC graph, but not in the decendants's and vice versa, i.e. EC numbers which only exist in either one of the organism groups.
        """
        divergedECs = ancestorEcGraph.getECs()
        divergedECs.symmetric_difference(descendantEcGraph.getECs())
        return divergedECs
    
    @staticmethod
    def getGraph(ancestorEcGraph: SubstanceEcGraph, descendantEcGraph: SubstanceEcGraph) -> SubstanceEcGraph:
        """
        Get graph containing EC numbers which have diverged between ancestor and descendant, existing only in either one of them.
        
        Parameters
        ----------
        ancestorEcGraph : SubstanceEcGraph
        descendantEcGraph : SubstanceEcGraph
        
        Returns
        -------
        SubstanceEcGraph
            Graph of EC numbers which occur in the ancestor's EC graph, but not in the decendants's and vice versa, i.e. EC numbers which only exist in either one of the organism groups.
            Substance-EC-product edges are only included if both graphs, ancestor and descendant, have both nodes, substrate and product.
        """
        addedGraph = GeneFunctionAddition.getGraph(ancestorEcGraph, descendantEcGraph)
        lostGraph = GeneFunctionLoss.getGraph(ancestorEcGraph, descendantEcGraph)
        divergedGraph = addedGraph.union(lostGraph, addCount=False)
        return divergedGraph






class GeneDuplication(object): 
    """
    Abstract class for any type of gene duplication.
    """

class SimpleGroupGeneDuplication(GeneDuplication):
    """
    Evolutionary event of duplicating a gene, regardless of ancestoral bonds in the comparison group.
    
    The conditions for a 'simple group' gene duplication are:
        - The gene has at least one homolog within the set of organisms its organism belongs to.
    """
    def __init__(self, sameGroupOrganisms: 'Iterable[Organism] or KEGG.Organism.Group'):
        """
        Simple group gene duplication extends simple gene duplication by expanding the term 'paralog' to every organism in the set of organisms the gene's organism blongs to.
        In contrast to :class:`SimpleGeneDuplication`, this class has to be instantiated, using the aforementioned set of organisms belonging to each other.
        This would usually be a :class:`FEV_KEGG.KEGG.Organism.Group` of the same :class:`FEV_KEGG.Evolution.Clade.Clade`. 
        
        Parameters
        ----------
        sameGroupOrganisms : Iterable[Organism] or Organism.Group
            Organisms which will be searched for the occurence of homologs, i.e. are considered "semi-paralogously" related.
        
        Attributes
        ----------
        self.sameGroupOrganisms : Iterable[Organism]
        
        Raises
        ------
        ValueError
            If `sameGroupOrganisms` is of wrong type.
        
        Warnings
        --------
        This takes much longer than :class:`SimpleGeneDuplication`, because each sought gene is compared between all organisms of the same group, not only within its own organism.
        """
        if isinstance(sameGroupOrganisms, Group):
            self.sameGroupOrganisms = sameGroupOrganisms.organisms
        elif isinstance(sameGroupOrganisms, Iterable):
            self.sameGroupOrganisms = sameGroupOrganisms
        else:
            raise ValueError("'sameGroupOrganisms' must be of type Iterable or KEGG.Organism.Group")
    
    def getEnzymes(self, enzymes: 'Set[Enzyme] or SubstanceEnzymeGraph', eValue = defaultEValue, returnMatches = False, ignoreDuplicatesOutsideSet = None):
        """
        Get gene-duplicated enzymes.
        
        Parameters
        ----------
        enzymes : Set[Enzyme] or SubstanceEnzymeGraph
            Set of enzymes to be checked for gene duplication, or a graph.
        eValue : float, optional
            Threshold for the statistical expectation value (E-value), below which a sequence alignment is considered significant.
        returnMatches : bool, optional
            If *True*, return not only `enzymes` that have homologs, but also which homologs they have. Useful for filtering for relevant homologs afterwards.
        ignoreDuplicatesOutsideSet : Set[GeneID], optional
            If *None*, report all found duplicates.
            If not *None*, count only such enzymes as gene duplicated, which have at least one of their duplicates inside this set.
            This can, for example, serve to exclude duplicates in secondary metabolism.
        
        Returns
        -------
        Set[Enzyme] or Dict[Element, Set[GeneID]]
            If returnMatches == *False*, all enzymes in `enzymes` which fulfil the conditions of this gene duplication definition.
            If returnMatches == *True*, all enzymes in `enzymes` which fulfil the conditions of this gene duplication definition, pointing to a set of gene IDs of the found homologs.
        
        Raises
        ------
        ValueError
            If any organism does not exist.
        HTTPError
            If any gene does not exist.
        URLError
            If connection to KEGG fails.
        """
        if isinstance(enzymes, SubstanceEnzymeGraph):
            enzymes = enzymes.getEnzymes()
        
        # find gene duplicated enzymes according to simple gene duplication
        shallReturnMatches = returnMatches or (ignoreDuplicatesOutsideSet is not None)
        simpleGeneDuplicates = SimpleGeneDuplication.getEnzymes(enzymes, eValue, returnMatches = shallReturnMatches, ignoreDuplicatesOutsideSet = ignoreDuplicatesOutsideSet)
        
        if returnMatches or ignoreDuplicatesOutsideSet is not None:
            simpleGeneDuplicatesSet = simpleGeneDuplicates.keys()
        else:
            simpleGeneDuplicatesSet = simpleGeneDuplicates
        
        duplicatedEnzymes = dict()
        
        # those enzymes already duplicated according to simple gene duplication do not have to be tested again
        soughtEnzymes = enzymes.copy()
        soughtEnzymes.difference_update( simpleGeneDuplicatesSet )
            
        # get orthologs
        geneIDs = [enzyme.geneID for enzyme in soughtEnzymes]
        
        if returnMatches or ignoreDuplicatesOutsideSet is not None: # need to get all orthologs
            
            matchingsDict = Database.getOrthologsBulk(geneIDs, self.sameGroupOrganisms, eValue) # GeneID -> List[ Matching ]
            
            for enzyme in enzymes: # for all enzymes, even the ones already found as paralogs
                
                # add paralogs
                paralogousGeneIDs = simpleGeneDuplicates.get(enzyme.geneID)
                if paralogousGeneIDs is not None:
                    
                    if len(paralogousGeneIDs) > 0:
                        
                        duplicatedEnzymes[enzyme] = paralogousGeneIDs
                    
                    continue # can not be in orthologs if it was already paralogous
                
                # add orthologs
                orthologousMatchings = matchingsDict.get(enzyme.geneID)
                if orthologousMatchings is not None:
                    
                    if len(orthologousMatchings) > 0:
                        
                        matches = []
                        for matching in orthologousMatchings:
                            matches.extend(matching.matches)
                        
                        matchedGeneIDs = {match.foundGeneID for match in matches}
                        
                        if len(matchedGeneIDs) > 0:
                            currentDuplicates = duplicatedEnzymes.get(enzyme)
                            
                            if currentDuplicates is None:
                                currentDuplicates = set()
                                duplicatedEnzymes[enzyme] = currentDuplicates
                                
                            currentDuplicates.update(matchedGeneIDs)
            
        else: # only interesting IF there are orthologs
        
            orthologousOrganismsDict = Database.hasOrthologsBulk(geneIDs, self.sameGroupOrganisms, eValue) # GeneID -> List[ organismAbbreviation ]
            
            # regarding each enzyme sought
            for enzyme in soughtEnzymes:
                # count orthologous GeneIDs
                orthologousOrganisms = orthologousOrganismsDict.get(enzyme.geneID)
                if orthologousOrganisms is None:
                    continue
                if len(orthologousOrganisms) > 0: # has at least one ortholog in any other organism
                    duplicatedEnzymes[enzyme] = None #orthologousOrganisms
        
        if ignoreDuplicatesOutsideSet is not None:
            filteredDuplicatedEnzymes = dict()
            
            for enzyme, matchGeneIDs in duplicatedEnzymes.items():
                
                matchesInSearchedSet = ignoreDuplicatesOutsideSet.intersection( matchGeneIDs )
                
                if len(matchesInSearchedSet) > 0: # some of the matches are in the set of enzymes to be checked for duplicates 
                    filteredDuplicatedEnzymes[enzyme] = matchesInSearchedSet
            
            duplicatedEnzymes = filteredDuplicatedEnzymes        
        
        if returnMatches:
            return duplicatedEnzymes
        else:
            return set(duplicatedEnzymes.keys())
    
    def getEnzymePairs(self, enzymes: 'Set[Enzyme] or SubstanceEnzymeGraph', eValue = defaultEValue, ignoreDuplicatesOutsideSet = None, geneIdToEnzyme = None) -> Set[Tuple[Enzyme, Enzyme]]:
        """
        Get gene-duplicated enzymes, in pairs of duplicates.
        
        If enzymeA is a duplicate of enzymeB and vice versa, this returns symmetric duplicates of the form (enzymeA, enzymeB) and (enzymeB, enzymeA).
        
        Parameters
        ----------
        enzymes : Set[Enzyme] or SubstanceEnzymeGraph
            Set of enzymes to be checked for gene duplication, or a graph.
        eValue : float, optional
            Threshold for the statistical expectation value (E-value), below which a sequence alignment is considered significant.
        ignoreDuplicatesOutsideSet : Set[GeneID], optional
            If *None*, report all found duplicates.
            If not *None*, count only such enzymes as gene duplicated, which have at least one of their duplicates inside this set.
            This can, for example, serve to exclude duplicates in secondary metabolism.
        geneIdToEnzyme : Dict[GeneID, Enzyme], optional
            Dictionary for mapping each gene ID of every found duplicate to an enzyme object.
            If *None*, gets the enzyme from the database. This avoids the KeyError, but can cause a lot of network load.
        
        Returns
        -------
        Set[Tuple[Enzyme, Enzyme]]
            Set of pairs of gene-duplicated enzymes, realised as tuples. The order is arbitrary and there will almost certainly be 100% duplicates.
        
        Raises
        ------
        KeyError
            If `geneIdToEnzyme` is passed, but does not contain the gene ID of every duplicate.
        """
        duplicatedEnzymeMatches = self.getEnzymes(enzymes, eValue, returnMatches = True, ignoreDuplicatesOutsideSet = ignoreDuplicatesOutsideSet)
        
        # expand matches of homologous gene IDs to pairs of duplicated enzymes. Which we can do (without wasting further resources) only here, because only here we have the geneID -> enzyme dict!
        duplicatedEnzymePairs = set()
        
        
        if geneIdToEnzyme is None: # need to get enzyme objects from database
            allGeneIDs = set()
            for geneIDs in duplicatedEnzymeMatches.values():
                allGeneIDs.update( geneIDs )
            
            geneIdToEnzyme = dict()
            
            for geneID, gene in Database.getGeneBulk(allGeneIDs).items():
                enzyme = Enzyme.fromGene(gene)
                geneIdToEnzyme[geneID] = enzyme
            
        
        for enzymeA, geneIDs in duplicatedEnzymeMatches.items():
            
            for geneID in geneIDs:
                
                enzymeB = geneIdToEnzyme[geneID]
                
                duplicatedEnzymePairs.add( (enzymeA, enzymeB) )
        
        return duplicatedEnzymePairs
    
    def filterEnzymes(self, substanceEnzymeGraph: SubstanceEnzymeGraph, eValue = defaultEValue) -> SubstanceEnzymeGraph:
        """
        Remove all enzymes from a graph which have not been gene-duplicated.
        
        Parameters
        ----------
        substanceEnzymeGraph : SubstanceEnzymeGraph
            Graph of enzymes to be checked for gene duplication.
        eValue : float, optional
            Threshold for the statistical expectation value (E-value), below which a sequence alignment is considered significant.
        
        Returns
        -------
        SubstanceEnzymeGraph
            A copy of the `substanceEnzymeGraph` containing only enzymes which fulfil the conditions of this gene duplication definition.
        
        Raises
        ------
        ValueError
            If any organism does not exist.
        HTTPError
            If any gene does not exist.
        URLError
            If connection to KEGG fails.
        """
        graph = substanceEnzymeGraph.copy()
        possibleGeneDuplicates = self.getEnzymes(substanceEnzymeGraph, eValue)
        graph.removeAllEnzymesExcept( possibleGeneDuplicates )
        return graph


class SimpleGeneDuplication(GeneDuplication):
    """
    Evolutionary event of duplicating a gene, regardless of ancestoral bonds.
    
    The conditions for a 'simple' gene duplication are:
        - The gene has at least one paralog.
    """    
    @staticmethod
    def getEnzymes(enzymes: 'Set[Enzyme] or SubstanceEnzymeGraph', eValue = defaultEValue, returnMatches = False, ignoreDuplicatesOutsideSet = None, preCalculatedEnzymes = None):
        """
        Get gene-duplicated enzymes.
        
        Parameters
        ----------
        enzymes : Set[Enzyme] or SubstanceEnzymeGraph
            Set of enzymes to be checked for gene duplication, or a graph.
        eValue : float, optional
            Threshold for the statistical expectation value (E-value), below which a sequence alignment is considered significant.
        returnMatches : bool, optional
            If *True*, return not only `enzymes` that have homologs, but also which homologs they have. Useful for filtering for relevant homologs afterwards.
        ignoreDuplicatesOutsideSet : Set[GeneID] or *True*, optional
            If *None*, report all found duplicates.
            If *True*, automatically restrict to all enzymes in `substanceEnzymeGraph`.
            If a set, count only such enzymes as gene duplicated, which have at least one of their duplicates inside this set. Beware, the set has to contain the enzymes' gene ID!
            This can, for example, serve to exclude duplicates in secondary metabolism.
        
        Returns
        -------
        Set[Enzyme] or Dict[Element, Set[GeneID]]
            If returnMatches == *False*, all enzymes in `enzymes` which fulfil the conditions of this gene duplication definition.
            If returnMatches == *True*, all enzymes in `enzymes` which fulfil the conditions of this gene duplication definition, pointing to a set of gene IDs of the found homologs.
        
        Raises
        ------
        ValueError
            If any organism does not exist.
        HTTPError
            If any gene does not exist.
        URLError
            If connection to KEGG fails.
        """
        if isinstance(enzymes, SubstanceEnzymeGraph):
            enzymes = enzymes.getEnzymes()
        
        
        if preCalculatedEnzymes is not None:
            
            if returnMatches is True:
                result = dict()
            else:
                result = set()
            
            # filter preCalculatedEnzymes by enzymes
            for enzyme, duplicatesSet in preCalculatedEnzymes.items():
                
                if enzyme in enzymes:
                        
                    validDuplicates = set()
                    for duplicate in duplicatesSet:
                        
                        if ignoreDuplicatesOutsideSet is None or duplicate.geneID in ignoreDuplicatesOutsideSet:
                            validDuplicates.add( duplicate.geneID )
                    
                    if len(validDuplicates) > 0:
                        
                        if returnMatches is True:
                            result[enzyme] = validDuplicates
                        
                        else:
                            result.add( enzyme )
            
            return result
        
        
        possibleGeneDuplicates = dict()
        
        geneIDs = {enzyme.geneID for enzyme in enzymes}
        
        matchingsDict = Database.getParalogsBulk(geneIDs, eValue)
        
        for enzyme in enzymes:
            matching = matchingsDict.get(enzyme.geneID, None)
            
            if matching is None:
                print( 'WARNING: data for GeneID ' + str(enzyme.geneID) + ' could not be downloaded. Maybe you want to exclude this erroneous organism completely, before it skews statistics? See quirks.py')
                continue
            
            paralogs = matching.matches
            if len(paralogs) > 0:
                possibleGeneDuplicates[enzyme] = [match.foundGeneID for match in paralogs]
        
        if ignoreDuplicatesOutsideSet is True:
            ignoreDuplicatesOutsideSet = geneIDs
        
        if ignoreDuplicatesOutsideSet is not None and len(ignoreDuplicatesOutsideSet) > 0:
            filteredPossibleGeneDuplicates = dict()
            
            for enzyme, matchGeneIDs in possibleGeneDuplicates.items():
                
                # which genes have been found AND are in the set of relevant duplicates?
                matchesInSearchedSet = ignoreDuplicatesOutsideSet.intersection( matchGeneIDs )
                
                if len(matchesInSearchedSet) > 0: # some of the matches are in the set of relevant duplicates
                    filteredPossibleGeneDuplicates[enzyme] = matchesInSearchedSet
            
            possibleGeneDuplicates = filteredPossibleGeneDuplicates
        
        if returnMatches:
            return possibleGeneDuplicates
        else:
            return set(possibleGeneDuplicates.keys())
    
    @classmethod
    def getEnzymePairs(cls, enzymes: 'Set[Enzyme] or SubstanceEnzymeGraph', eValue = defaultEValue, ignoreDuplicatesOutsideSet = None, geneIdToEnzyme = None, preCalculatedEnzymes = None) -> Set[Tuple[Enzyme, Enzyme]]:
        """
        Get gene-duplicated enzymes, in pairs of duplicates.
        
        If enzyme A is a duplicate of enzyme B and vice versa, this does not return duplicates, but returns only one pair, with the "smaller" enzyme as the first value. An enzyme is "smaller" if its gene ID string is "smaller".
        
        Parameters
        ----------
        enzymes : Set[Enzyme] or SubstanceEnzymeGraph
            Set of enzymes to be checked for gene duplication, or a graph.
        eValue : float, optional
            Threshold for the statistical expectation value (E-value), below which a sequence alignment is considered significant.
        ignoreDuplicatesOutsideSet : Set[GeneID] or *True*, optional
            If *None*, report all found duplicates.
            If *True*, automatically restrict to all enzymes in `substanceEnzymeGraph`.
            If a set, count only such enzymes as gene duplicated, which have at least one of their duplicates inside this set. Beware, the set has to contain the enzymes' gene ID!
            This can, for example, serve to exclude duplicates in secondary metabolism.
        geneIdToEnzyme : Dict[GeneID, Enzyme], optional
            Dictionary for mapping each gene ID of every found duplicate to an enzyme object.
            If *None*, gets the enzyme from the database. This avoids the KeyError, but can cause a lot of network load.
        
        Returns
        -------
        Set[Tuple[Enzyme, Enzyme]]
            Set of pairs of gene-duplicated enzymes, realised as tuples. The order is arbitrary.
        
        Raises
        ------
        KeyError
            If `geneIdToEnzyme` is passed, but does not contain the gene ID of every duplicate.
        """
        duplicatedEnzymeMatches = cls.getEnzymes(enzymes, eValue, returnMatches = True, ignoreDuplicatesOutsideSet = ignoreDuplicatesOutsideSet, preCalculatedEnzymes = preCalculatedEnzymes)
        
        # expand matches of homologous gene IDs to pairs of duplicated enzymes. Which we can do (without wasting further resources) only here, because only here we have the geneID -> enzyme dict!
        duplicatedEnzymePairs = set()
        
        
        if geneIdToEnzyme is None: # need to get enzyme objects from database
            allGeneIDs = set()
            for geneIDs in duplicatedEnzymeMatches.values():
                allGeneIDs.update( geneIDs )
            
            geneIdToEnzyme = dict()
            
            for geneID, gene in Database.getGeneBulk(allGeneIDs).items():
                enzyme = Enzyme.fromGene(gene)
                geneIdToEnzyme[geneID] = enzyme
            
        
        for enzymeA, geneIDs in duplicatedEnzymeMatches.items():
            
            for geneID in geneIDs:
                
                enzymeB = geneIdToEnzyme[geneID]
                
                duplicatedEnzymePairs.add( (enzymeA, enzymeB) )
        
        # filter symmetric duplicates
        deduplicatedEnzymePairs = set()
        for enzymeA, enzymeB in duplicatedEnzymePairs:
            if enzymeA <= enzymeB:
                deduplicatedEnzymePairs.add( (enzymeA, enzymeB) )
            else:
                deduplicatedEnzymePairs.add( (enzymeB, enzymeA) )
        
        return deduplicatedEnzymePairs
        
    
    @classmethod
    def filterEnzymes(cls, substanceEnzymeGraph: SubstanceEnzymeGraph, eValue = defaultEValue, ignoreDuplicatesOutsideSet = None, preCalculatedEnzymes = None) -> SubstanceEnzymeGraph:
        """
        Remove all enzymes from a graph which have not been gene-duplicated.
        
        Parameters
        ----------
        substanceEnzymeGraph : SubstanceEnzymeGraph
            Graph of enzymes to be checked for gene duplication.
        eValue : float, optional
            Threshold for the statistical expectation value (E-value), below which a sequence alignment is considered significant.
        ignoreDuplicatesOutsideSet : Set[GeneID] or *True*, optional
            If *None*, report all found duplicates.
            If *True*, automatically restrict to all enzymes in `substanceEnzymeGraph`.
            If a set, count only such enzymes as gene duplicated, which have at least one of their duplicates inside this set. Beware, the set has to contain the enzymes' gene ID!
            This can, for example, serve to exclude duplicates in secondary metabolism.
        
        Returns
        -------
        SubstanceEnzymeGraph
            A copy of the `substanceEnzymeGraph` containing only enzymes which fulfil the conditions of this gene duplication definition.
        
        Raises
        ------
        ValueError
            If any organism does not exist.
        HTTPError
            If any gene does not exist.
        URLError
            If connection to KEGG fails.
        """
        graph = substanceEnzymeGraph.copy()
        possibleGeneDuplicates = cls.getEnzymes(substanceEnzymeGraph, eValue, returnMatches = False, ignoreDuplicatesOutsideSet = ignoreDuplicatesOutsideSet, preCalculatedEnzymes = preCalculatedEnzymes)
        graph.removeAllEnzymesExcept( possibleGeneDuplicates )
        return graph
        
    

class ChevronGeneDuplication(GeneDuplication):
    """
    Evolutionary event of duplicating a gene, in dependence of a certain ancestoral bond.
    
    The conditions for a 'chevron' gene duplication are:
        - The gene has at least one paralog.
        - The gene has at least one ortholog in a pre-defined set of organisms.
    """
    def __init__(self, possiblyOrthologousOrganisms: 'Iterable[Organism] or KEGG.Organism.Group'):
        """
        Chevron gene duplication extends simple gene duplication by limiting the possibly duplicated genes via a set of possibly orthologous organisms.
        In contrast to :class:`SimpleGeneDuplication`, this class has to be instantiated, using the aforementioned set of possibly orthologous organisms. 
        
        Parameters
        ----------
        possiblyOrthologousOrganisms : Iterable[Organism] or Organism.Group
            Organisms which will be searched for the occurence of orthologs, i.e. are considered ancestoral.
        
        Attributes
        ----------
        self.possiblyOrthologousOrganisms : Iterable[Organism]
        
        Raises
        ------
        ValueError
            If `possiblyOrthologousOrganisms` is of wrong type.
        
        Warnings
        --------
        This takes much longer than :class:`SimpleGeneDuplication`, because additionally, each found paralog is searched for an ortholog in all organisms of the other group.
        However, if you set `returnMatches` == *False* and `ignoreDuplicatesOutsideSelf` == *False*, the search is aborted with the very first ortholog, which is much faster than getting all orthologs.
        Because in this model even a single orthologous match is enough to prove gene duplication, we do not necessarily have to fully search all organisms.
        """
        if isinstance(possiblyOrthologousOrganisms, Group):
            self.possiblyOrthologousOrganisms = possiblyOrthologousOrganisms.organisms
        elif isinstance(possiblyOrthologousOrganisms, Iterable):
            self.possiblyOrthologousOrganisms = possiblyOrthologousOrganisms
        else:
            raise ValueError("'possiblyOrthologusOrganisms' must be of type Iterable or KEGG.Organism.Group")
            
    
    def getEnzymes(self, enzymes: 'Set[Enzyme] or SubstanceEnzymeGraph', eValue = defaultEValue, returnMatches = False, ignoreDuplicatesOutsideSelf = False):
        """
        Get gene-duplicated enzymes.
        
        Parameters
        ----------
        enzymes : Set[Enzyme] or SubstanceEnzymeGraph
            Set of enzymes to be checked for gene duplication, or a graph.
        eValue : float, optional
            Threshold for the statistical expectation value (E-value), below which a sequence alignment is considered significant.
        returnMatches : bool, optional
            If *True*, return not only `enzymes` that have homologs, but also which homologs they have. Useful for filtering for relevant homologs afterwards.
        ignoreDuplicatesOutsideSelf : bool, optional
            If *True*, count only such enzymes as gene duplicated, which have at least one of their duplicates inside the set of enzymes searched for duplicates.
            This can, for example, serve to exclude duplicates in secondary metabolism.
        
        Returns
        -------
        Set[Enzyme] or Dict[Element, Set[GeneID]]
            If returnMatches == *False*, all enzymes in `enzymes` which fulfil the conditions of this gene duplication definition.
            If returnMatches == *True*, all enzymes in `enzymes` which fulfil the conditions of this gene duplication definition, pointing to a set of gene IDs of the found homologs.
        
        Raises
        ------
        ValueError
            If any organism does not exist.
        HTTPError
            If any gene does not exist.
        URLError
            If connection to KEGG fails.
        """
        # graph -> set
        if isinstance(enzymes, SubstanceEnzymeGraph):
            enzymes = enzymes.getEnzymes()
        
        # get paralogs first
        possibleGeneDuplicates = SimpleGeneDuplication.getEnzymes(enzymes, eValue, returnMatches = returnMatches or ignoreDuplicatesOutsideSelf, ignoreDuplicatesOutsideSet = {enzyme.geneID for enzyme in enzymes}) # always need returnMatches if ignoreDuplicatesOutsideSelf is True! 
        
        if returnMatches or ignoreDuplicatesOutsideSelf:
            possibleGeneDuplicatesSet = possibleGeneDuplicates.keys()
        else:
            possibleGeneDuplicatesSet = possibleGeneDuplicates
        
        if len( possibleGeneDuplicatesSet ) == 0: # nothing to do, because there are no paralogs
            return possibleGeneDuplicates
        
        duplicatedEnzymes = dict()
        
        # get orthologs
        geneIDs = {enzyme.geneID for enzyme in possibleGeneDuplicatesSet}
        
        if returnMatches or ignoreDuplicatesOutsideSelf: # need to get all orthologs
            
            matchingsDict = Database.getOrthologsBulk(geneIDs, self.possiblyOrthologousOrganisms, eValue) # GeneID -> List[ Matching ]
            
            for enzyme in possibleGeneDuplicatesSet:
                
                # add orthologs
                orthologousMatchings = matchingsDict.get(enzyme.geneID)
                if orthologousMatchings is not None:

                    if len(orthologousMatchings) > 0:
                        
                        matches = []
                        for matching in orthologousMatchings:
                            matches.extend(matching.matches)
                            
                        duplicatedEnzymes[enzyme] = {match.foundGeneID for match in matches}
                
                # add paralogs
                paralogousGeneIDs = possibleGeneDuplicates.get(enzyme.geneID)
                if paralogousGeneIDs is not None:
                    
                    if len(paralogousGeneIDs) > 0:
                        
                        currentDuplicates = duplicatedEnzymes.get(enzyme)
                        if currentDuplicates is None:
                            currentDuplicates = set()
                            duplicatedEnzymes[enzyme] = currentDuplicates
                        
                        currentDuplicates.update(paralogousGeneIDs)
            
        else: # only interesting IF there are orthologs
            
            orthologousOrganismsDict = Database.hasOrthologsBulk(geneIDs, self.possiblyOrthologousOrganisms, eValue) # GeneID -> List[ organismAbbreviation ]
            
            for enzyme in possibleGeneDuplicatesSet:
            
                orthologousOrganisms = orthologousOrganismsDict.get(enzyme.geneID)
                if orthologousOrganisms is None:
                    continue
                if len(orthologousOrganisms) > 0:
                    duplicatedEnzymes[enzyme] = None #orthologousOrganisms
        
        if ignoreDuplicatesOutsideSelf:
            filteredDuplicatedEnzymes = dict()
            
            for enzyme, matchGeneIDs in duplicatedEnzymes.items():
                
                matchesInSearchedSet = geneIDs.intersection( matchGeneIDs )
                
                if len(matchesInSearchedSet) > 0: # some of the matches are in the set of enzymes to be checked for duplicates 
                    filteredDuplicatedEnzymes[enzyme] = matchesInSearchedSet
            
            duplicatedEnzymes = filteredDuplicatedEnzymes        
        
        if returnMatches:
            return duplicatedEnzymes
        else:
            return set(duplicatedEnzymes.keys())
    
    def getEnzymePairs(self, enzymes: 'Set[Enzyme] or SubstanceEnzymeGraph', eValue = defaultEValue, ignoreDuplicatesOutsideSelf = False, geneIdToEnzyme = None) -> Set[Tuple[Enzyme, Enzyme]]:
        """
        Get gene-duplicated enzymes, in pairs of duplicates.
        
        If enzyme A is a duplicate of enzyme B and vice versa, this does not return duplicates, but returns only one pair, with the "smaller" enzyme as the first value. An enzyme is "smaller" if its gene ID string is "smaller".
        
        Parameters
        ----------
        enzymes : Set[Enzyme] or SubstanceEnzymeGraph
            Set of enzymes to be checked for gene duplication, or a graph.
        eValue : float, optional
            Threshold for the statistical expectation value (E-value), below which a sequence alignment is considered significant.
        ignoreDuplicatesOutsideSelf : bool, optional
            If *True*, count only such enzymes as gene duplicated, which have at least one of their duplicates inside the set of enzymes searched for duplicates.
            This can, for example, serve to exclude duplicates in secondary metabolism.
        geneIdToEnzyme : Dict[GeneID, Enzyme], optional
            Dictionary for mapping each gene ID of every found duplicate to an enzyme object.
            If *None*, gets the enzyme from the database. This avoids the KeyError, but can cause a lot of network load.
        
        Returns
        -------
        Set[Tuple[Enzyme, Enzyme]]
            Set of pairs of gene-duplicated enzymes, realised as tuples. The order is arbitrary.
        
        Raises
        ------
        KeyError
            If `geneIdToEnzyme` is passed, but does not contain the gene ID of every duplicate.
        """
        duplicatedEnzymeMatches = self.getEnzymes(enzymes, eValue, returnMatches = True, ignoreDuplicatesOutsideSelf = ignoreDuplicatesOutsideSelf)
        
        # expand matches of homologous gene IDs to pairs of duplicated enzymes. Which we can do (without wasting further resources) only here, because only here we have the geneID -> enzyme dict!
        duplicatedEnzymePairs = set()
        
        
        if geneIdToEnzyme is None: # need to get enzyme objects from database
            allGeneIDs = set()
            for geneIDs in duplicatedEnzymeMatches.values():
                allGeneIDs.update( geneIDs )
            
            geneIdToEnzyme = dict()
            
            for geneID, gene in Database.getGeneBulk(allGeneIDs).items():
                enzyme = Enzyme.fromGene(gene)
                geneIdToEnzyme[geneID] = enzyme
            
        
        for enzymeA, geneIDs in duplicatedEnzymeMatches.items():
            
            for geneID in geneIDs:
                
                enzymeB = geneIdToEnzyme[geneID]
                
                duplicatedEnzymePairs.add( (enzymeA, enzymeB) )
        
        # filter symmetric duplicates
        deduplicatedEnzymePairs = set()
        for enzymeA, enzymeB in duplicatedEnzymePairs:
            if enzymeA <= enzymeB:
                deduplicatedEnzymePairs.add( (enzymeA, enzymeB) )
            else:
                deduplicatedEnzymePairs.add( (enzymeB, enzymeA) )
        
        return deduplicatedEnzymePairs
    
    def filterEnzymes(self, substanceEnzymeGraph: SubstanceEnzymeGraph, eValue = defaultEValue, ignoreDuplicatesOutsideSelf = False) -> SubstanceEnzymeGraph:
        """
        Remove all enzymes from a graph which have not been gene-duplicated.
        
        Parameters
        ----------
        substanceEnzymeGraph : SubstanceEnzymeGraph
            Graph of enzymes to be checked for gene duplication.
        eValue : float, optional
            Threshold for the statistical expectation value (E-value), below which a sequence alignment is considered significant.
        ignoreDuplicatesOutsideSelf : bool, optional
            If *True*, count only such enzymes as gene duplicated, which have at least one of their duplicates inside the set of enzymes searched for duplicates.
            This can, for example, serve to exclude duplicates in secondary metabolism.
        
        Returns
        -------
        SubstanceEnzymeGraph
            A copy of the `substanceEnzymeGraph` containing only enzymes which fulfil the conditions of this gene duplication definition.
        
        Raises
        ------
        ValueError
            If any organism does not exist.
        HTTPError
            If any gene does not exist.
        URLError
            If connection to KEGG fails.
        """
        graph = substanceEnzymeGraph.copy()
        possibleGeneDuplicates = self.getEnzymes(substanceEnzymeGraph, eValue, returnMatches = False, ignoreDuplicatesOutsideSelf = ignoreDuplicatesOutsideSelf)
        graph.removeAllEnzymesExcept( possibleGeneDuplicates )
        return graph








class Neofunctionalisation():
    
    def __init__(self, enzymeA: Enzyme, enzymeB: Enzyme):
        """
        Evolutionary event of Neofunctionalisation between a pair of enzymes.
        
        The conditions for a neofunctionalisation are:
            - The enzyme's gene has been duplicated, according to a certain class of GeneDuplication.
            - The duplicated enzyme is associated with a different EC number than its duplicate.
        
        The order of the two enzymes has no meaning, it has been arbitrarily chosen to reflect the lexicographic order of their associated EC numbers. The enzyme posessing the "smallest" EC number comes first.
        This absolute ordering prevents duplicate events, because without an order there would have always been a second event with the exact same enzymes, but in swapped positions, because neofunctionalisation has no direction here and is, thus, symmetric.
        
        Parameters
        ----------
        enzymeA : Enzyme
            An enzyme, which is a gene duplicate of `enzymeB`. The order is arbitrary.
        enzymeB : Enzyme
            An enzyme, which is a gene duplicate of `enzymeA`. The order is arbitrary.
        
        Attributes
        ----------
        self.enzymePair : Tuple[Enzyme, Enzyme]
            Tuple of the two enzymes, sorted by the lexicographic order of their "smallest" EC number.
        
        Raises
        ------
        ValueError
            If the enzymes are equal, have the same set of EC numbers, or one has no EC number.
        """
        if enzymeA == enzymeB:
            raise ValueError('The enzymes must be unequal!')
        
        ecNumbersAset = enzymeA.ecNumbers
        ecNumbersBset = enzymeB.ecNumbers
        
        if ecNumbersAset == ecNumbersBset:
            raise ValueError('The enzymes must have differing EC numbers to be NEOfunctionalised!')
        
        if len(ecNumbersAset) == 0 or len(ecNumbersBset) == 0:
            raise ValueError('The enzymes have to be associated with  at least one EC number!')
        
        ecNumbersA = list(ecNumbersAset)
        ecNumbersA.sort()
        
        ecNumbersB = list(ecNumbersBset)
        ecNumbersB.sort()
        
        if len(ecNumbersA) <= len(ecNumbersB):
            smallerSet = ecNumbersA
            biggerSet = ecNumbersB
            isAsmaller = True
            
        else:
            smallerSet = ecNumbersB
            biggerSet = ecNumbersA
            isAsmaller = False
        
        for index, ec1 in enumerate(smallerSet):
            ec2 = biggerSet[index]
            
            if ec1 == ec2:
                if index == len(smallerSet) - 1: # end of smaller set reached
                    if isAsmaller:
                        self.enzymePair = (enzymeA, enzymeB)
                    else:
                        self.enzymePair = (enzymeB, enzymeA)
                
                else: # not yet at the end, compare next indizes
                    continue
            
            elif ec1 < ec2:
                if isAsmaller:
                    self.enzymePair = (enzymeA, enzymeB)
                else:
                    self.enzymePair = (enzymeB, enzymeA)
            
            else:
                if isAsmaller:
                    self.enzymePair = (enzymeB, enzymeA)
                else:
                    self.enzymePair = (enzymeA, enzymeB)
    
    def getEnzymes(self) -> Tuple[Enzyme, Enzyme]:
        """
        Get the pair of enzymes.
        
        Returns
        -------
        Tuple[Enzyme, Enzyme]
        """
        return self.enzymePair
    
    def getEcNumbers(self) -> Tuple[Set[EcNumber], Set[EcNumber]]:
        """
        Get the enzymes' EC numbers.
        
        Returns
        -------
        Tuple[Set[EcNumber], Set[EcNumber]]
            Same order as in :func:`getEnzymes`.
            Because an enzyme could have multiple EC numbers, they are given as sets.
        """
        return (self.enzymePair[0].ecNumbers, self.enzymePair[1].ecNumbers)
    
    def getDifferingEcLevels(self) -> int:
        """
        Get the maximum number of EC levels in which the enzymes' EC numbers differ.
        
        Returns
        -------
        int
            Number of differing EC levels between the two enzymes' EC numbers, starting with the substrate-level.
            If an enzyme has multiple EC numbers, returns the biggest difference.
            For example 1.2.3.4 and 1.2.3.7 returns 1, while 1.2.3.4 and 1.8.9.10 returns 3.
            However, wildcards do **not** match numbers: 1.2.3.4 and 1.2.3.- returns 1!
        """
        biggestDifference = 0
        
        for ecA in self.getEcNumbers()[0]:
            
            for ecB in self.getEcNumbers()[1]:
                
                difference = 4 - ecA.matchingLevels(ecB, wildcardMatchesNumber = False)
                
                if difference > biggestDifference:
                    biggestDifference = difference
        
        return biggestDifference
    
    def isSameEcReaction(self) -> bool:
        """
        Whether the enzymes' EC numbers describe a different reaction, or merely a different substrate
        
        Returns
        -------
        bool
            *True*, if the biggest difference in EC levels of the two enzymes is still on the fourth (substrate) level.
        """
        if self.getDifferingEcLevels() > 1:
            return False
        else:
            return True
    
    def toHtml(self, short = False):
        """
        Get the string representation as an HTML line.
        """
        enzymePair = self.enzymePair
        ecPair = self.getEcNumbers()
        return '<td>' + enzymePair[0].toHtml() + '<td>' + ',&nbsp;'.join([ec.toHtml(short) for ec in sorted(ecPair[0])]) + '</td></td><td><-></td><td>'+ enzymePair[1].toHtml() + '<td>' + ',&nbsp;'.join([ec.toHtml(short) for ec in sorted(ecPair[1])]) + '</td></td>'
    
    def __str__(self):
        enzymePair = self.enzymePair
        ecPair = self.getEcNumbers()
        return '(' + str(enzymePair[0]) + ' [' + ', '.join([str(ec) for ec in sorted(ecPair[0])]) + '],\t'+ str(enzymePair[1]) + ' [' + ', '.join([str(ec) for ec in sorted(ecPair[1])]) + '])'
    
    def __repr__(self):
        return self.__str__()
        
    def __eq__(self, other):
        if isinstance(self, other.__class__):
            return self.enzymePair == other.enzymePair
        return False
        
    def __ne__(self, other):
        return not self == other
    
    def __hash__(self):
        return self.enzymePair.__hash__()
    
    def __lt__(self, other):
        
        # sort by EC number first
        selfEnzyme1 = self.enzymePair[0]
        selfEnzyme2 = self.enzymePair[1]
        
        selfEnzyme1EcList = list(selfEnzyme1.ecNumbers)
        selfEnzyme2EcList = list(selfEnzyme2.ecNumbers)
        
        otherEnzyme1 = other.enzymePair[0]
        otherEnzyme2 = other.enzymePair[1]
        
        otherEnzyme1EcList = list(otherEnzyme1.ecNumbers)
        otherEnzyme2EcList = list(otherEnzyme2.ecNumbers)

        if selfEnzyme1EcList == otherEnzyme1EcList:
            
            if selfEnzyme2EcList == otherEnzyme2EcList:
                
                # then by gene ID
                if selfEnzyme1.uniqueID == otherEnzyme1.uniqueID:
                    return selfEnzyme2.uniqueID < otherEnzyme2.uniqueID
                
                else:
                    return selfEnzyme1.uniqueID < otherEnzyme1.uniqueID
                
            else:
                return selfEnzyme2EcList < otherEnzyme2EcList

        else:
            return selfEnzyme1EcList < otherEnzyme1EcList
    
    def __gt__(self, other):
        
        # sort by EC number first
        selfEnzyme1 = self.enzymePair[0]
        selfEnzyme2 = self.enzymePair[1]
        
        selfEnzyme1EcList = list(selfEnzyme1.ecNumbers)
        selfEnzyme2EcList = list(selfEnzyme2.ecNumbers)
        
        otherEnzyme1 = other.enzymePair[0]
        otherEnzyme2 = other.enzymePair[1]
        
        otherEnzyme1EcList = list(otherEnzyme1.ecNumbers)
        otherEnzyme2EcList = list(otherEnzyme2.ecNumbers)

        if selfEnzyme1EcList == otherEnzyme1EcList:
            
            if selfEnzyme2EcList == otherEnzyme2EcList:
                
                # then by gene ID
                if selfEnzyme1.uniqueID == otherEnzyme1.uniqueID:
                    return selfEnzyme2.uniqueID > otherEnzyme2.uniqueID
                
                else:
                    return selfEnzyme1.uniqueID > otherEnzyme1.uniqueID
                
            else:
                return selfEnzyme2EcList > otherEnzyme2EcList

        else:
            return selfEnzyme1EcList > otherEnzyme1EcList
    
    def __le__(self, other):
        
        # sort by EC number first
        selfEnzyme1 = self.enzymePair[0]
        selfEnzyme2 = self.enzymePair[1]
        
        selfEnzyme1EcList = list(selfEnzyme1.ecNumbers)
        selfEnzyme2EcList = list(selfEnzyme2.ecNumbers)
        
        otherEnzyme1 = other.enzymePair[0]
        otherEnzyme2 = other.enzymePair[1]
        
        otherEnzyme1EcList = list(otherEnzyme1.ecNumbers)
        otherEnzyme2EcList = list(otherEnzyme2.ecNumbers)

        if selfEnzyme1EcList == otherEnzyme1EcList:
            
            if selfEnzyme2EcList == otherEnzyme2EcList:
                
                # then by gene ID
                if selfEnzyme1.uniqueID == otherEnzyme1.uniqueID:
                    return selfEnzyme2.uniqueID <= otherEnzyme2.uniqueID
                
                else:
                    return selfEnzyme1.uniqueID <= otherEnzyme1.uniqueID
                
            else:
                return selfEnzyme2EcList <= otherEnzyme2EcList

        else:
            return selfEnzyme1EcList <= otherEnzyme1EcList
    
    def __ge__(self, other):
        
        # sort by EC number first
        selfEnzyme1 = self.enzymePair[0]
        selfEnzyme2 = self.enzymePair[1]
        
        selfEnzyme1EcList = list(selfEnzyme1.ecNumbers)
        selfEnzyme2EcList = list(selfEnzyme2.ecNumbers)
        
        otherEnzyme1 = other.enzymePair[0]
        otherEnzyme2 = other.enzymePair[1]
        
        otherEnzyme1EcList = list(otherEnzyme1.ecNumbers)
        otherEnzyme2EcList = list(otherEnzyme2.ecNumbers)

        if selfEnzyme1EcList == otherEnzyme1EcList:
            
            if selfEnzyme2EcList == otherEnzyme2EcList:
                
                # then by gene ID
                if selfEnzyme1.uniqueID == otherEnzyme1.uniqueID:
                    return selfEnzyme2.uniqueID >= otherEnzyme2.uniqueID
                
                else:
                    return selfEnzyme1.uniqueID >= otherEnzyme1.uniqueID
                
            else:
                return selfEnzyme2EcList >= otherEnzyme2EcList

        else:
            return selfEnzyme1EcList >= otherEnzyme1EcList




class NeofunctionalisedEnzymes():
    
    def __init__(self, enzymes: Set[Enzyme], geneDuplicationModel: GeneDuplication, eValue = defaultEValue, ignoreDuplicatesOutsideSet: bool = True):
        """
        Neofunctionalisation events among certain `enzymes`.
        
        Parameters
        ----------
        enzymes : Set[Enzyme]
            Enzymes among which to test for neofunctionalisation. Neofunctionalisations involving enzymes outside this set are **not** reported.
        geneDuplicationModel : GeneDuplication
            The model of gene duplication to use.
        eValue : float, optional
            Threshold for the statistical expectation value (E-value), below which a sequence alignment is considered significant.
        ignoreDuplicatesOutsideSet : bool, optional
            If *True*, any neofunctionalisation involving an enzyme outside the `enzymes` set is not reported.
            This helps to exclude secondary metabolism when examining core metabolism.
        
        Raises
        ------
        ValueError
            If a gene duplication model is used which requires instantiation, but only its class was given.
        """
        if geneDuplicationModel == ChevronGeneDuplication:
            raise ValueError("Chevron gene duplication model requires you to instantiate an object, parametrised with the set of possibly orthologous organisms.")
        
        elif geneDuplicationModel == SimpleGroupGeneDuplication:
            raise ValueError("Simple group gene duplication model requires you to instantiate an object, parametrised with the set of organisms belonging to the same group.")
        
        self._geneDuplicationModel = geneDuplicationModel
        
        self._neofunctionalisations = set()
        
        # get possibly gene-duplicated enzymes, pointing to their duplicates
        geneIdToEnzyme = dict()
        for enzyme in enzymes:
            geneIdToEnzyme[enzyme.geneID] = enzyme
        
        if ignoreDuplicatesOutsideSet is True: # restrict allowed gene duplicates to the set of passed enzymes. This helps to avoid parts of the metabolism outside the core metabolism
            
            if isinstance(geneDuplicationModel, ChevronGeneDuplication):
                duplicatedEnzymePairs = geneDuplicationModel.getEnzymePairs(enzymes, eValue, ignoreDuplicatesOutsideSelf = ignoreDuplicatesOutsideSet, geneIdToEnzyme = geneIdToEnzyme)
            else:
                duplicatedEnzymePairs = geneDuplicationModel.getEnzymePairs(enzymes, eValue, ignoreDuplicatesOutsideSet = set(geneIdToEnzyme.keys()), geneIdToEnzyme = geneIdToEnzyme)
        
        else: # you want to check for neofunctionalisation outside the passed enzymes (usually core metabolism)
            duplicatedEnzymePairs = geneDuplicationModel.getEnzymePairs(enzymes, eValue)
        
        # in all pairs of duplicated enzymes, find neofunctionalised ones
        for enzymePair in duplicatedEnzymePairs:
            try:
                neofunctionalisation = Neofunctionalisation(enzymePair[0], enzymePair[1])
                self._neofunctionalisations.add( neofunctionalisation )
            
            except ValueError: # obviously not neofunctionalised
                pass # ignore
    
    
    def getNeofunctionalisations(self, minimumEcDifference: int = None) -> Set[Neofunctionalisation]:
        """
        Get neofunctionalisation events between two enzymes each.
        
        Parameters
        ----------
        minimumEcDifference : int, optional
            May only be one of [1, 2, 3, 4].
            If *None* or *1*, all neofunctionalisations are returned.
            If > *1*, return only neofunctionalisations in which the EC numbers differ in more than the `minimumEcDifference` lowest levels.
            They then describe a different reaction, instead of only a different substrate.
            For example, `minimumEcDifference` == *2* means that 1.2.3.4/1.2.3.5 is not reported, while 1.2.3.4/1.2.5.6 is.
        
        Returns
        -------
        Set[Neofunctionalisation]
            Set of possible neofunctionalisation events.
        """
        if minimumEcDifference is not None and minimumEcDifference > 1:
            filteredNeofunctionalisations = set()
            
            for neofunctionalisation in self._neofunctionalisations:
                if neofunctionalisation.getDifferingEcLevels() >= minimumEcDifference:
                    filteredNeofunctionalisations.add( neofunctionalisation )
            
            return filteredNeofunctionalisations
            
        else:
            return self._neofunctionalisations.copy()            
        
    def getEnzymes(self, minimumEcDifference: int = None) -> Set[Enzyme]:
        """
        Get all neofunctionalised enzymes.
        
        Parameters
        ----------
        minimumEcDifference : int, optional
            May only be one of [1, 2, 3, 4].
            If *None* or *1*, all neofunctionalisations are returned.
            If > *1*, return only neofunctionalisations in which the EC numbers differ in more than the `minimumEcDifference` lowest levels.
            They then describe a different reaction, instead of only a different substrate.
            For example, `minimumEcDifference` == *2* means that 1.2.3.4/1.2.3.5 is not reported, while 1.2.3.4/1.2.5.6 is.
        
        Returns
        -------
        Set[Enzyme]
            Set of all possibly neofunctionalised enzymes, regardless of the real direction of neofunctionalisation, which we can not determine here.
        """
        neofunctionalisedEnzymes = set()
        
        for neofunctionalisation in self.getNeofunctionalisations(minimumEcDifference):
            
            neofunctionalisedEnzymes.update( neofunctionalisation.getEnzymes() )
        
        return neofunctionalisedEnzymes
    
    def filterGraph(self, enzymeGraph: SubstanceEnzymeGraph, minimumEcDifference: int = None) -> SubstanceEnzymeGraph:
        """
        Filter enzyme graph to only contain neofunctionalised enzymes.
        
        Parameters
        ----------
        enzymeGraph : SubstanceEnzymeGraph
            The enzyme graph to filter.
        minimumEcDifference : int, optional
            May only be one of [1, 2, 3, 4].
            If *None* or *1*, all neofunctionalisations are returned.
            If > *1*, return only neofunctionalisations in which the EC numbers differ in more than the `minimumEcDifference` lowest levels.
            They then describe a different reaction, instead of only a different substrate.
            For example, `minimumEcDifference` == *2* means that 1.2.3.4/1.2.3.5 is not reported, while 1.2.3.4/1.2.5.6 is.
        
        Returns
        -------
        SubstanceEnzymeGraph
            A copy of `enzymeGraph`, leaving only edges with a neofunctionalised enzyme as key.
        """
        graph = enzymeGraph.copy()
        neofunctionalisedEnzymes = self.getEnzymes(minimumEcDifference)
        graph.removeAllEnzymesExcept( neofunctionalisedEnzymes )
        return graph
    
    def colourGraph(self, enzymeGraph: SubstanceEnzymeGraph, colour: Export.Colour = Export.Colour.GREEN, minimumEcDifference: int = None) -> SubstanceEnzymeGraph:
        """
        Colour enzyme graph's neofunctionalised enzyme edges.
        
        Parameters
        ----------
        enzymeGraph : SubstanceEnzymeGraph
            The enzyme graph to colour.
        colour : Export.Colour, optional
            The colour to use for edges with neofunctionalised enzymes as key.
        minimumEcDifference : int, optional
            May only be one of [1, 2, 3, 4].
            If *None* or *1*, all neofunctionalisations are returned.
            If > *1*, return only neofunctionalisations in which the EC numbers differ in more than the `minimumEcDifference` lowest levels.
            They then describe a different reaction, instead of only a different substrate.
            For example, `minimumEcDifference` == *2* means that 1.2.3.4/1.2.3.5 is not reported, while 1.2.3.4/1.2.5.6 is.
        
        Returns
        -------
        SubstanceEnzymeGraph
            A copy of `enzymeGraph` in which edges with neofunctionalised enzymes as key have an additional colour attribute, see :func:`FEV_KEGG.Drawing.Export.addColourAttribute`.
        """
        graph = enzymeGraph.copy()
        neofunctionalisedEnzymes = self.getEnzymes(minimumEcDifference)
        Export.addColourAttribute(graph, colour, nodes = False, edges = neofunctionalisedEnzymes)
        return graph




class FunctionChange():
    
    def __init__(self, ecA: EcNumber, ecB: EcNumber):
        """
        Possible evolutionary change of enzymatic function, from one EC number to another.
        
        The direction of change, or if it really happened, can not be determined here!
        The order of EC numbers in this object is arbitrarily chosen to reflect their lexicographic order.
        
        A function change resembles the possibility that the first EC number has evolutionarily changed into the second one. Or the other way around, since the direction of evolution can not be determined here.
        A function change can never have the same EC number twice, nor can the first be lexicographically "bigger" than the second, see the examples section.
        
        Parameters
        ----------
        ecA: EcNumber
        ecB: EcNumber
            Must be lexicographically "bigger" than `ecA`. Can not be equal to `ecA`.
        
        Raises
        ------
        ValueError
            If `ecA` is lexicographically "bigger" than, or equal to, `ecB`.
        """
        if ecA >= ecB:
            raise ValueError("EC number 1 is bigger than or equal to EC number 2.")
        
        self.ecA = ecA
        self.ecB = ecB
        self.ecPair = (ecA, ecB)
    
    @classmethod
    def fromNeofunctionalisation(cls, neofunctionalisation: Neofunctionalisation) -> Set['FunctionChange']:
        """
        Create combinations of function changes from a `neofunctionalisation`.
        
        Parameters
        ----------
        neofunctionalisation : Neofunctionalisation
        
        Returns
        -------
        Set[FunctionChange]
            Set of function changes which might have been caused by the `neofunctionalisation`.
            
            Since an enzyme of a neofunctionalisation can have multiple EC numbers, all combinations of the two enzymes' EC numbers are formed and treated as separate possible function changes.
        
        Examples
        --------
        A: 1
        B: 2
        = (1, 2)
        
        A: 1
        B: 1, 2
        = (1, 2)
        
        A: 1, 2
        B: 1, 3
        = (1, 3) (1, 2) (2, 3)
        
        A: 1, 2
        B: 3, 4
        = (1, 3) (1, 4) (2, 3) (2, 4)
        
        A: 1, 2
        B: 1, 2, 3
        = (1, 3) (2, 3)
        """
        ecNumbersA, ecNumbersB = neofunctionalisation.getEcNumbers()
            
        # cross product
        questionableEcPairs = itertools.product(ecNumbersA, ecNumbersB)
        
        # filter illegal products
        functionChanges = set()
        for pair in questionableEcPairs:
            try:
                functionChanges.add(cls(pair[0], pair[1]))
            except ValueError:
                pass
        
        return functionChanges
    
    def getDifferingEcLevels(self) -> int:
        """
        Get the number of EC levels in which the EC numbers differ.
        
        Returns
        -------
        int
            Number of differing EC levels between the two EC numbers, starting with the substrate-level.
            For example 1.2.3.4 and 1.2.3.7 returns 1, while 1.2.3.4 and 1.8.9.10 returns 3.
            However, wildcards do **not** match numbers: 1.2.3.4 and 1.2.3.- returns 1!
        """
        return 4 - self.ecA.matchingLevels(self.ecB, wildcardMatchesNumber = False)
    
    def toHtml(self, short = False):
        """
        Get the string representation as an HTML line.
        """
        return '<td>' + self.ecPair[0].toHtml(short) + '</td><td><-></td><td>' + self.ecPair[1].toHtml(short) + '</td>'
    
    def __str__(self):
        return self.ecPair.__str__()
    
    def __repr__(self):
        return self.__str__()
        
    def __eq__(self, other):
        if isinstance(self, other.__class__):
            return self.ecPair == other.ecPair
        return False
        
    def __ne__(self, other):
        return not self == other
    
    def __hash__(self):
        return self.ecPair.__hash__()
    
    def __lt__(self, other):
        return self.ecPair < other.ecPair
    
    def __gt__(self, other):
        return self.ecPair > other.ecPair
    
    def __le__(self, other):
        return self.ecPair <= other.ecPair
    
    def __ge__(self, other):
        return self.ecPair >= other.ecPair




class NeofunctionalisedECs():
    
    def __init__(self, neofunctionalisedEnzymes: NeofunctionalisedEnzymes):
        """
        EC numbers which are affected by neofunctionalisation events.
        
        Parameters
        ----------
        neofunctionalisedEnzymes : NeofunctionalisedEnzymes
            Neofunctionalisation events among certain enzymes.
        """
        self._neofunctionalisedEnzymes = neofunctionalisedEnzymes
    
    
    def getNeofunctionalisationsForFunctionChange(self, minimumEcDifference: int = None, minimumOrganismsCount: int = None) -> Dict[FunctionChange, Set[Neofunctionalisation]]:
        """
        Get neofunctionalsation events, keyed by a change of function between the two enzymes.        
        
        Parameters
        ----------
        minimumEcDifference : int, optional
            May only be one of [1, 2, 3, 4].
            If *None* or *1*, all neofunctionalisations are returned.
            If > *1*, return only neofunctionalisations in which the EC numbers differ in more than the `minimumEcDifference` lowest levels.
            They then describe a different reaction, instead of only a different substrate.
            For example, `minimumEcDifference` == *2* means that 1.2.3.4/1.2.3.5 is not reported, while 1.2.3.4/1.2.5.6 is.
        minimumOrganismsCount : int, optional
            Minimum number of organisms which have to be involved in the neofunctionalisations of each function change.
            If *None*, there is no filtering due to organism involvement.
            For example, the function change 1->2 is associated with two neofunctionalisations 'eco:12345'->'eco:69875' and 'obc:76535'->'abc:41356', this involves three organisms in total (eco, obc, abc), finally, if `minimumOrganismsCount` <= 3, the function change 1->2 is returned.
        
        Returns
        -------
        Dict[FunctionChange, Set[Neofunctionalisation]]
            Dictionary of function changes, pointing to a set of neofunctionalisations which might have caused them.
            
            Since an enzyme of a neofunctionalisation can have multiple EC numbers, all combinations of the two enzymes' EC numbers are formed and treated as separate possible function changes.
            The neofunctionalisation is then saved again for each function change, which obviously leads to duplicated neofunctionalisation objects.
        
        Examples
        --------
        A: 1
        B: 2
        = (1, 2)
        
        A: 1
        B: 1, 2
        = (1, 2)
        
        A: 1, 2
        B: 1, 3
        = (1, 3) (1, 2) (2, 3)
        
        A: 1, 2
        B: 3, 4
        = (1, 3) (1, 4) (2, 3) (2, 4)
        
        A: 1, 2
        B: 1, 2, 3
        = (1, 3) (2, 3)
        """
        neofunctionalisationForFunctionChange = dict()
        
        # get neofunctionalisations
        for neofunctionalisation in self._neofunctionalisedEnzymes.getNeofunctionalisations(minimumEcDifference):
            
            # split sets of EC numbers into pair-wise combinations
            functionChanges = FunctionChange.fromNeofunctionalisation(neofunctionalisation)
            
            for functionChange in functionChanges:
                currentSet = neofunctionalisationForFunctionChange.get(functionChange, None)
                
                if currentSet is None:
                    currentSet = set()
                    neofunctionalisationForFunctionChange[functionChange] = currentSet
                
                currentSet.add(neofunctionalisation)
        
        # filter function changes with neofunctionalised enzymes which stem from too few organisms
        if minimumOrganismsCount is not None:
            keysToDelete = []
            
            for functionChange, neofunctionalisations in neofunctionalisationForFunctionChange.items():
                
                enzymes = set()
                
                for neofunctionalisation in neofunctionalisations:    
                    enzymes.update( neofunctionalisation.getEnzymes() )
                
                organisms = set()
                
                for enzyme in enzymes:
                    organisms.add( enzyme.organismAbbreviation )
                
                # enough occuring organisms?
                if not len(organisms) >= minimumOrganismsCount:
                    keysToDelete.append( functionChange )
            
            # delete keys which failed the test
            for key in keysToDelete:
                del neofunctionalisationForFunctionChange[key]
        
        return neofunctionalisationForFunctionChange
    
    def getFunctionChanges(self, minimumEcDifference: int = None, minimumOrganismsCount: int = None) -> Set[FunctionChange]:
        """
        Get all possible changes of function between the two enzymes of every neofunctionalisation.
        
        Parameters
        ----------
        minimumEcDifference : int, optional
            May only be one of [1, 2, 3, 4].
            If *None* or *1*, all neofunctionalisations are returned.
            If > *1*, return only neofunctionalisations in which the EC numbers differ in more than the `minimumEcDifference` lowest levels.
            They then describe a different reaction, instead of only a different substrate.
            For example, `minimumEcDifference` == *2* means that 1.2.3.4/1.2.3.5 is not reported, while 1.2.3.4/1.2.5.6 is.
        minimumOrganismsCount : int, optional
            Minimum number of organisms which have to be involved in the neofunctionalisations of each function change.
            If *None*, there is no filtering due to organism involvement.
        
        Returns
        -------
        Set[FunctionChange]
            Set of all function changes, which meet the criteria.
        """
        return set( self.getNeofunctionalisationsForFunctionChange(minimumEcDifference, minimumOrganismsCount).keys() )
    
    def getEnzymesForFunctionChange(self, minimumEcDifference: int = None, minimumOrganismsCount: int = None) -> Dict[FunctionChange, Set[Enzyme]]:
        """
        Get enzymes of neofunctionalisations, keyed by a possible change of function.
        
        Parameters
        ----------
        minimumEcDifference : int, optional
            May only be one of [1, 2, 3, 4].
            If *None* or *1*, all neofunctionalisations are returned.
            If > *1*, return only neofunctionalisations in which the EC numbers differ in more than the `minimumEcDifference` lowest levels.
            They then describe a different reaction, instead of only a different substrate.
            For example, `minimumEcDifference` == *2* means that 1.2.3.4/1.2.3.5 is not reported, while 1.2.3.4/1.2.5.6 is.
        minimumOrganismsCount : int, optional
            Minimum number of organisms which have to be involved in the neofunctionalisations of each function change.
            If *None*, there is no filtering due to organism involvement.
            For example, the function change 1->2 is associated with two neofunctionalisations 'eco:12345'->'eco:69875' and 'obc:76535'->'abc:41356', this involves three organisms in total (eco, obc, abc), finally, if `minimumOrganismsCount` <= 3, the function change 1->2 is returned.
        
        Returns
        -------
        Dict[FunctionChange, Set[Enzyme]]
            Dictionary of function changes, pointing to a set of enzymes involved in the neofunctionalisations which might have caused the function change.
            This can lead to many duplicated enzymes.
        """
        enzymesForFunctionChange = dict()
        
        for functionChange, neofunctionalisations in self.getNeofunctionalisationsForFunctionChange(minimumEcDifference, minimumOrganismsCount).items():
            currentSet = enzymesForFunctionChange.get(functionChange, None)
                
            if currentSet is None:
                currentSet = set()
                enzymesForFunctionChange[functionChange] = currentSet
            
            for neofunctionalisation in neofunctionalisations:
                currentSet.update( neofunctionalisation.getEnzymes() )
            
        return enzymesForFunctionChange
    
    
    
    def getNeofunctionalisationsForEC(self, minimumEcDifference: int = None, minimumOrganismsCount: int = None) -> Dict[EcNumber, Set[Neofunctionalisation]]:
        """
        Get neofunctionalisation events, keyed by an EC number participating in the change of function between the two enzymes.
        
        Parameters
        ----------
        minimumEcDifference : int, optional
            May only be one of [1, 2, 3, 4].
            If *None* or *1*, all neofunctionalisations are returned.
            If > *1*, return only neofunctionalisations in which the EC numbers differ in more than the `minimumEcDifference` lowest levels.
            They then describe a different reaction, instead of only a different substrate.
            For example, `minimumEcDifference` == *2* means that 1.2.3.4/1.2.3.5 is not reported, while 1.2.3.4/1.2.5.6 is.
        minimumOrganismsCount : int, optional
            Minimum number of organisms which have to be involved in the neofunctionalisations of each EC number.
            If *None*, there is no filtering due to organism involvement.
            This sums the occurences of organisms across function changes, for each EC number the function changes overlap with. Hence, it is much less likely that a neofunctionalisation is filtered, compared to filtering per function change.
            For example, the function change 1->2 is associated with two neofunctionalisations 'eco:12345'->'eco:69875' and 'obc:76535'->'abc:41356', this involves three organisms in total (eco, obc, abc).
            Also, the function change 1->3 involves two organisms ('eco:53235'->'iuf:34587'). If `minimumOrganismsCount` == 4, neither 1->2, nor 1->3 are reported.
            However, if we look at single EC numbers, 1 is involved in function changes affecting four organisms (eco, obc, abc, iuf). Thus, 1 would be reported here, but neither 2 nor 3.
        
        Returns
        -------
        Dict[EcNumber, Set[Neofunctionalisation]]
            Dictionary of EC numbers which are part of function changes, pointing to a set of neofunctionalisations which might have caused them.
            Very likely has duplicated neofunctionalisations, because there are always at least two EC numbers involved in a neofunctionalisation.
        """
        neofunctionalisationForEC = dict()
        
        # get neofunctionalisations
        for neofunctionalisation in self._neofunctionalisedEnzymes.getNeofunctionalisations(minimumEcDifference):
            
            # split sets of EC numbers into pair-wise combinations
            functionChanges = FunctionChange.fromNeofunctionalisation(neofunctionalisation)
            
            for functionChange in functionChanges:
                
                # for each EC number of a function change, save neofunctionalisation
                for ec in functionChange.ecPair:
                    currentSet = neofunctionalisationForEC.get(ec, None)
                    
                    if currentSet is None:
                        currentSet = set()
                        neofunctionalisationForEC[ec] = currentSet
                    
                    currentSet.add(neofunctionalisation)
        
        # filter ECs with neofunctionalised enzymes which stem from too few organisms
        if minimumOrganismsCount is not None:
            keysToDelete = []
            
            for ec, neofunctionalisations in neofunctionalisationForEC.items():
                
                enzymes = set()
                
                for neofunctionalisation in neofunctionalisations:    
                    enzymes.update( neofunctionalisation.getEnzymes() )
                
                organisms = set()
                
                for enzyme in enzymes:
                    organisms.add( enzyme.organismAbbreviation )
                
                # enough occuring organisms?
                if not len(organisms) >= minimumOrganismsCount:
                    keysToDelete.append( ec )
            
            # delete keys which failed the test
            for key in keysToDelete:
                del neofunctionalisationForEC[key]
        
        return neofunctionalisationForEC
    
    def getECs(self, minimumEcDifference: int = None, minimumOrganismsCount: int = None) -> Set[EcNumber]:
        """
        Get EC numbers participating in the change of function due to neofunctionalisations.
        
        They could also be called "neofunctionalised" EC numbers.
        
        Parameters
        ----------
        minimumEcDifference : int, optional
            May only be one of [1, 2, 3, 4].
            If *None* or *1*, all neofunctionalisations are returned.
            If > *1*, return only neofunctionalisations in which the EC numbers differ in more than the `minimumEcDifference` lowest levels.
            They then describe a different reaction, instead of only a different substrate.
            For example, `minimumEcDifference` == *2* means that 1.2.3.4/1.2.3.5 is not reported, while 1.2.3.4/1.2.5.6 is.
        minimumOrganismsCount : int, optional
            Minimum number of organisms which have to be involved in the neofunctionalisations of each EC number.
            If *None*, there is no filtering due to organism involvement.
            This sums the occurences of organisms across function changes, for each EC number the function changes overlap with. Hence, it is much less likely that a neofunctionalisation is filtered, compared to filtering per function change.
            For example, the function change 1->2 is associated with two neofunctionalisations 'eco:12345'->'eco:69875' and 'obc:76535'->'abc:41356', this involves three organisms in total (eco, obc, abc).
            Also, the function change 1->3 involves two organisms ('eco:53235'->'iuf:34587'). If `minimumOrganismsCount` == 4, neither 1->2, nor 1->3 are reported.
            However, if we look at single EC numbers, 1 is involved in function changes affecting four organisms (eco, obc, abc, iuf). Thus, 1 would be reported here, but neither 2 nor 3.
        
        Returns
        -------
        Set[EcNumber]
            Set of EC numbers which are part of function changes which possibly happened due to neofunctionalisations.
        """
        return set( self.getNeofunctionalisationsForEC(minimumEcDifference, minimumOrganismsCount).keys() )
    
    def getEnzymesForEC(self, minimumEcDifference: int = None, minimumOrganismsCount: int = None) -> Dict[EcNumber, Set[Enzyme]]:
        """
        Get enzymes of neofunctionalisations, keyed by an EC number of a possible function change.
        
        Parameters
        ----------
        minimumEcDifference : int, optional
            May only be one of [1, 2, 3, 4].
            If *None* or *1*, all neofunctionalisations are returned.
            If > *1*, return only neofunctionalisations in which the EC numbers differ in more than the `minimumEcDifference` lowest levels.
            They then describe a different reaction, instead of only a different substrate.
            For example, `minimumEcDifference` == *2* means that 1.2.3.4/1.2.3.5 is not reported, while 1.2.3.4/1.2.5.6 is.
        minimumOrganismsCount : int, optional
            Minimum number of organisms which have to be involved in the neofunctionalisations of each EC number.
            If *None*, there is no filtering due to organism involvement.
            This sums the occurences of organisms across function changes, for each EC number the function changes overlap with. Hence, it is much less likely that a neofunctionalisation is filtered, compared to filtering per function change.
            For example, the function change 1->2 is associated with two neofunctionalisations 'eco:12345'->'eco:69875' and 'obc:76535'->'abc:41356', this involves three organisms in total (eco, obc, abc).
            Also, the function change 1->3 involves two organisms ('eco:53235'->'iuf:34587'). If `minimumOrganismsCount` == 4, neither 1->2, nor 1->3 are reported.
            However, if we look at single EC numbers, 1 is involved in function changes affecting four organisms (eco, obc, abc, iuf). Thus, 1 would be reported here, but neither 2 nor 3.
        
        Returns
        -------
        Dict[EcNumber, Set[Enzyme]]
            Dictionary of EC numbers, pointing to a set of enzymes involved in the neofunctionalisations which might have caused the function changes the EC number is part of.
            This can lead to many duplicated enzymes.
        """
        enzymesForEC = dict()
        
        for ec, neofunctionalisations in self.getNeofunctionalisationsForEC(minimumEcDifference, minimumOrganismsCount).items():
            currentSet = enzymesForEC.get(ec, None)
                
            if currentSet is None:
                currentSet = set()
                enzymesForEC[ec] = currentSet
            
            for neofunctionalisation in neofunctionalisations:
                currentSet.update( neofunctionalisation.getEnzymes() )
            
        return enzymesForEC
    
    
    
    def filterGraph(self, ecGraph: SubstanceEcGraph, minimumEcDifference: int = None, minimumOrganismsCount: int = None) -> SubstanceEcGraph:
        """
        Filter EC graph to only contain "neofunctionalised" EC numbers.
        
        Parameters
        ----------
        minimumEcDifference : int, optional
            May only be one of [1, 2, 3, 4].
            If *None* or *1*, all neofunctionalisations are returned.
            If > *1*, return only neofunctionalisations in which the EC numbers differ in more than the `minimumEcDifference` lowest levels.
            They then describe a different reaction, instead of only a different substrate.
            For example, `minimumEcDifference` == *2* means that 1.2.3.4/1.2.3.5 is not reported, while 1.2.3.4/1.2.5.6 is.
        minimumOrganismsCount : int, optional
            Minimum number of organisms which have to be involved in the neofunctionalisations of each EC number.
            If *None*, there is no filtering due to organism involvement.
            This sums the occurences of organisms across function changes, for each EC number the function changes overlap with. Hence, it is much less likely that a neofunctionalisation is filtered, compared to filtering per function change.
            For example, the function change 1->2 is associated with two neofunctionalisations 'eco:12345'->'eco:69875' and 'obc:76535'->'abc:41356', this involves three organisms in total (eco, obc, abc).
            Also, the function change 1->3 involves two organisms ('eco:53235'->'iuf:34587'). If `minimumOrganismsCount` == 4, neither 1->2, nor 1->3 are reported.
            However, if we look at single EC numbers, 1 is involved in function changes affecting four organisms (eco, obc, abc, iuf). Thus, 1 would be reported here, but neither 2 nor 3.
        
        Returns
        -------
        SubstanceEcGraph
            A copy of `ecGraph`, leaving only edges with a "neofunctionalised" EC as key.
        """
        graph = ecGraph.copy()
        neofunctionalisedECs = self.getECs(minimumEcDifference, minimumOrganismsCount)
        graph.removeAllECsExcept( neofunctionalisedECs )
        return graph
    
    def colourGraph(self, ecGraph: SubstanceEcGraph, colour: Export.Colour = Export.Colour.GREEN, minimumEcDifference: int = None, minimumOrganismsCount: int = None) -> SubstanceEcGraph:
        """
        Colour EC graph's "neofunctionalised" EC number edges.
        
        Parameters
        ----------
        minimumEcDifference : int, optional
            May only be one of [1, 2, 3, 4].
            If *None* or *1*, all neofunctionalisations are returned.
            If > *1*, return only neofunctionalisations in which the EC numbers differ in more than the `minimumEcDifference` lowest levels.
            They then describe a different reaction, instead of only a different substrate.
            For example, `minimumEcDifference` == *2* means that 1.2.3.4/1.2.3.5 is not reported, while 1.2.3.4/1.2.5.6 is.
        minimumOrganismsCount : int, optional
            Minimum number of organisms which have to be involved in the neofunctionalisations of each EC number.
            If *None*, there is no filtering due to organism involvement.
            This sums the occurences of organisms across function changes, for each EC number the function changes overlap with. Hence, it is much less likely that a neofunctionalisation is filtered, compared to filtering per function change.
            For example, the function change 1->2 is associated with two neofunctionalisations 'eco:12345'->'eco:69875' and 'obc:76535'->'abc:41356', this involves three organisms in total (eco, obc, abc).
            Also, the function change 1->3 involves two organisms ('eco:53235'->'iuf:34587'). If `minimumOrganismsCount` == 4, neither 1->2, nor 1->3 are reported.
            However, if we look at single EC numbers, 1 is involved in function changes affecting four organisms (eco, obc, abc, iuf). Thus, 1 would be reported here, but neither 2 nor 3.
        
        Returns
        -------
        SubstanceEcGraph
            A copy of `ecGraph` in which edges with "neofunctionalised" ECs as key have an additional colour attribute, see :func:`FEV_KEGG.Drawing.Export.addColourAttribute`.
        """
        graph = ecGraph.copy()
        neofunctionalisedECs = self.getECs(minimumEcDifference, minimumOrganismsCount)
        Export.addColourAttribute(graph, colour, nodes = False, edges = neofunctionalisedECs)
        return graph
    