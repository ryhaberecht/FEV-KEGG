from FEV_KEGG.Graph.SubstanceGraphs import SubstanceEcGraph, SubstanceEnzymeGraph
from FEV_KEGG.Evolution.Taxonomy import NCBI, Taxonomy
from FEV_KEGG.KEGG.Organism import Group 
from FEV_KEGG.Evolution.Events import GeneFunctionAddition, GeneFunctionLoss, GeneFunctionDivergence, GeneFunctionConservation, SimpleGeneDuplication,\
    NeofunctionalisedECs, NeofunctionalisedEnzymes, Neofunctionalisation, FunctionChange
from FEV_KEGG import settings
from builtins import str
from FEV_KEGG.Drawing import Export
import math
from typing import Dict, Set, Tuple
from FEV_KEGG.Graph.Elements import Enzyme, GeneID, EcNumber

defaultExcludeUnclassified = True
"""
If *True*, ignore taxons with a path containing the string 'unclassified'.
This can be overridden in each relevant method's `excludeUnclassified` parameter in this module.
"""

defaultExcludeMultifunctionalEnzymes = settings.defaultNoMultifunctional
"""
If *True*, ignore enzymes with more than one EC number.
This can be overridden in each relevant method's `excludeMultifunctionalEnzymes` parameter in this module.
"""

defaultMajorityPercentageCoreMetabolism = 80
"""
Default percentage of organisms in the clade, which have to possess an EC number, for it to be included in the core metabolism of the clade.
See :func:`FEV_KEGG.KEGG.Organism.Group.majorityEcGraph`.
This can be overridden in each relevant method's `majorityPercentageCoreMetabolism` parameter in this module.
"""

defaultMajorityPercentageNeofunctionalisation = 0
"""
Default percentage of organisms in the clade, which have to possess the same "neofunctionalised" EC number, for it to be included in the set of "neofunctionalised" EC numbers of the clade.
See :class:`FEV_KEGG.KEGG.Evolution.Events.NeofunctionalisedECs`.
This can be overridden in each relevant method's `majorityPercentageNeofunctionalisation` parameter in this module.
"""

defaultEValue = settings.defaultEvalue
"""
Default threshold for the statistical expectation value (E-value), below which a sequence alignment is considered significant.
"""

defaultOneOrganismPerSpecies = settings.defaultOneOrganismPerSpecies
"""
Default descision whether to use only the first organism for each species in NCBI taxonomy.
"""

class Clade(object):
    
    def __init__(self, ncbiNames: 'e.g. Enterobacter or Proteobacteria/Gammaproteobacteria. Allows list of names, e.g. ["Gammaproteobacteria", "/Archaea"]', excludeUnclassified = defaultExcludeUnclassified, oneOrganismPerSpecies = defaultOneOrganismPerSpecies):
        """
        A clade in NCBI taxonomy, containing all leaf taxon's KEGG organisms.
        
        Parameters
        ----------
        ncbiNames : str or Iterable[str]
            String(s) a taxon's path must contain to be included in this clade.
        excludeUnclassified : bool, optional
            If *True*, ignore taxons with a path containing the string 'unclassified'.
        oneOrganismPerSpecies : bool, optional
            If *True*, use only the first organism of each species.
        
        Attributes
        ----------
        self.ncbiNames : Iterable[str]
            Part of the path of each leaf taxon to be included in this clade. A single string is wrapped in a list.
        self.group
            The :class:`FEV_KEGG.KEGG.Organism.Group` of KEGG organisms created from the found leaf taxons.
        
        Raises
        ------
        ValueError
            If no clade with `ncbiNames` in its path could be found.
        
        Warnings
        --------
        It is possible to include organisms of several clades in the same Clade object!
        For example, if you were to search for `ncbiNames` == 'Donaldus Duckus', you would get every organism within '/Bacteria/Donaldus Duckus' **and** '/Archaea/Order/Donaldus Duckus'.
        Use the slash (/) notation to make sure you only get the taxon you want, e.g. 'Proteobacteria/Gammaproteobacteria' or '/Archaea'.
        """
        taxonomy = NCBI.getTaxonomy()
        
        if isinstance(ncbiNames, str):
            ncbiNames = [ncbiNames]
            
        self.ncbiNames = ncbiNames
        
        allOrganisms = set()
        for ncbiName in ncbiNames:
            organisms = taxonomy.getOrganismAbbreviationsByPath(ncbiName, exceptPaths=('unclassified' if excludeUnclassified else None), oneOrganismPerSpecies=oneOrganismPerSpecies)
            if organisms is None or len(organisms) == 0:
                raise ValueError("No clade of this path found: " + ncbiName)
            allOrganisms.update(organisms)
        
        self.group = Group( allOrganisms )
        
        self._lastNeofunctionalisedEnzymesCache = None
        self._lastGeneDuplicatedEnzymesMatches = None
    
    
    def collectiveMetabolism(self, excludeMultifunctionalEnzymes = defaultExcludeMultifunctionalEnzymes, addEcDescriptions = False) -> SubstanceEcGraph:
        """
        The Substance-EC graph representing the collective metabolic network, occuring in any organism of the clade.
        
        This includes each and every EC number which occurs in any organism of this clade.
        
        Parameters
        ----------
        excludeMultifunctionalEnzymes : bool, optional
            If *True*, ignore enzymes with more than one EC number.
        
        Returns
        -------
        SubstanceEcGraph
            Collective metabolic network of EC numbers, including counts of occurence in each of the clade's organisms.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        graph = self.group.collectiveEcGraph(noMultifunctional = excludeMultifunctionalEnzymes, addCount = True, keepOnHeap = True, addEcDescriptions = addEcDescriptions)
        graph.name = 'Collective metabolism ECs ' + ' '.join(self.ncbiNames)
        return graph
    
    def collectiveMetabolismEnzymes(self, excludeMultifunctionalEnzymes = defaultExcludeMultifunctionalEnzymes) -> SubstanceEnzymeGraph:
        """
        The Substance-Enzyme graph representing the collective metabolic network, occuring in any organism of the clade.
        
        This includes each and every enzyme of every organism of this clade.
        
        Parameters
        ----------
        excludeMultifunctionalEnzymes : bool, optional
            If *True*, ignore enzymes with more than one EC number.
        
        Returns
        -------
        SubstanceEnzymeGraph
            Collective metabolic network of enzymes.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        graph = self.group.collectiveEnzymeGraph(noMultifunctional = excludeMultifunctionalEnzymes, keepOnHeap = True)
        graph.name = 'Collective metabolism enzymes ' + ' '.join(self.ncbiNames)
        return graph
    
    def coreMetabolism(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, excludeMultifunctionalEnzymes = defaultExcludeMultifunctionalEnzymes) -> SubstanceEcGraph:
        """
        The Substance-EC graph representing the common metabolic network, shared among all organisms of the clade.
        
        This includes only EC numbers which occur in at least `majorityPercentageCoreMetabolism` % of all organisms of this clade.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            A path (substance -> EC -> product) has to occur in `majorityPercentageCoreMetabolism` % of the clade's organisms to be included.
        excludeMultifunctionalEnzymes : bool, optional
            If *True*, ignore enzymes with more than one EC number.
        
        Returns
        -------
        SubstanceEcGraph
            Core metabolic network of EC numbers.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        graph = self.group.majorityEcGraph(majorityPercentage = majorityPercentageCoreMetabolism, noMultifunctional = excludeMultifunctionalEnzymes, keepOnHeap = True)
        graph.name = 'Core metabolism ECs ' + ' '.join(self.ncbiNames)
        return graph
    
    def coreMetabolismEnzymes(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, excludeMultifunctionalEnzymes = defaultExcludeMultifunctionalEnzymes) -> SubstanceEnzymeGraph:
        """
        The Substance-Enzyme graph representing the common metabolic network, shared among all organisms of the clade.
        
        This includes every Enzyme associated with an EC number occuring in core metabolism (see :func:`substanceEcGraph`), no matter from which organism it stems.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            A path (substance -> EC -> product) has to occur in `majorityPercentageCoreMetabolism` % of the clade's organisms to be included.
        excludeMultifunctionalEnzymes : bool, optional
            If *True*, ignore enzymes with more than one EC number.
        
        Returns
        -------
        SubstanceEnzymeGraph
            Core metabolic network of enzymes.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        graph = self.group.collectiveEnzymeGraphByEcMajority(majorityPercentage = majorityPercentageCoreMetabolism, majorityTotal = None, noMultifunctional = excludeMultifunctionalEnzymes)
        graph.name = 'Core metabolism Enzymes ' + ' '.join(self.ncbiNames)
        return graph
    
    @property
    def organismsCount(self) -> int:
        """
        The number of organisms (leaf taxons) in this clade.
        
        Returns
        -------
        int
            The number of organisms (leaf taxons) in this clade.
        """
        return self.group.organismsCount
    
    
    
    
    # gene duplication
    
    def geneDuplicatedEnzymes(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, colour = False) -> SubstanceEnzymeGraph:
        """
        The substance-Enzyme graph of all gene duplicated enzymes of the core metabolism.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            Every substance-EC-product edge has to occur in `majorityPercentageCoreMetabolism` % of organisms constituting the clade, to be included in the core metabolism.
        colour : bool, optional
            If *True*, colours the gene-duplicated enzyme edges in green. The colouring is realised by adding a 'colour' attribute to each edge. Nodes are not coloured.
            Alternatively, you can specify a :class:`Export.Colour`.
        
        Returns
        -------
        SubstanceEnzymeGraph
            Substance-Enzyme graph containing all gene-duplicated enzymes, and nothing else.
            If `colour` == *True*, returns the full core metabolism enzyme graph, colouring gene-duplicated enzymes green.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        
                
        
        enzymeGraph = self.coreMetabolismEnzymes(majorityPercentageCoreMetabolism)
        
        geneDuplicationModel = SimpleGeneDuplication
#         geneDuplicationModel = SimpleGroupGeneDuplication(sameGroupOrganisms = self.group)
        
        # filter core metabolism enzyme graph    
        geneDuplicatedEnzymes = geneDuplicationModel.filterEnzymes(enzymeGraph, eValue = defaultEValue, ignoreDuplicatesOutsideSet = True, preCalculatedEnzymes = None)
        
        # colour core metabolism
        if colour is not False:
            
            if colour is True:
                colourToUse = Export.Colour.GREEN
            else:
                colourToUse = colour
            
            geneDuplicatedEnzymesOnly = geneDuplicatedEnzymes
            geneDuplicatedEnzymes = enzymeGraph
            Export.addColourAttribute(geneDuplicatedEnzymes, colourToUse, nodes = False, edges = geneDuplicatedEnzymesOnly.getEdges())
        
        geneDuplicatedEnzymes.name = 'Gene-duplicated core metabolism enzymes ' + ' '.join(self.ncbiNames)
        
        return geneDuplicatedEnzymes
    
    
    def geneDuplicatedEnzymesDict(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism) -> Dict[Enzyme, Set[GeneID]]:
        """
        All gene duplicated enzymes of the core metabolism, pointing to all their duplicates.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            Every substance-EC-product edge has to occur in `majorityPercentageCoreMetabolism` % of organisms constituting the clade, to be included in the core metabolism.
        
        Returns
        -------
        Dict[Enzyme, Set[GeneID]]
            Each gene ID on the right usually has an entry of its own, as an enzyme object, on the left, because they are each others homologs.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        
        
        enzymeGraph = self.coreMetabolismEnzymes(majorityPercentageCoreMetabolism)
        geneDuplicationModel = SimpleGeneDuplication
        
        geneIDsForEnzyme = geneDuplicationModel.getEnzymes(enzymeGraph, returnMatches = True, ignoreDuplicatesOutsideSet = True, preCalculatedEnzymes = None)
        
#         if keepOnHeap is True:
#             self._geneDuplicatedEnzymesObject = geneIDsForEnzyme
        
        return geneIDsForEnzyme
    
    
    def geneDuplicatedEnzymePairs(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism) -> Set[Tuple[Enzyme, Enzyme]]:
        """
        All gene duplicated enzymes of the core metabolism, paired with each of their duplicates.
        
        If enzyme A is a duplicate of enzyme B and vice versa, this does not return duplicates, but returns only one pair, with the "smaller" enzyme as the first value. An enzyme is "smaller" if its gene ID string is "smaller".
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            Every substance-EC-product edge has to occur in `majorityPercentageCoreMetabolism` % of organisms constituting the clade, to be included in the core metabolism.
        
        Returns
        -------
        Set[Tuple[Enzyme, Enzyme]]
            Set of gene-duplicated enzymes, broken down into pairs of enzymes.
            Can obviously create many duplicates left and right.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        
        
        enzymes = self.coreMetabolismEnzymes(majorityPercentageCoreMetabolism).getEnzymes()
        geneDuplicationModel = SimpleGeneDuplication
        
        geneIdToEnzyme = dict()
        for enzyme in enzymes:
            geneIdToEnzyme[enzyme.geneID] = enzyme
        
        enzymePairs = geneDuplicationModel.getEnzymePairs(enzymes, ignoreDuplicatesOutsideSet = True, geneIdToEnzyme = geneIdToEnzyme, preCalculatedEnzymes = None)
        
        return enzymePairs
        
            
        
        
    
    
    
    # neofunctionalisation
    
    def _neofunctionalisedEnzymes(self, majorityPercentageCoreMetabolism, eValue = defaultEValue, considerOnlyECs = None):
        
        # check if the last calculation can be returned
        if hasattr(self, '_lastNeofunctionalisedEnzymesCache') and self._lastNeofunctionalisedEnzymesCache is not None and considerOnlyECs is None:
            
            lastMajorityPercentage, lastNeofunctionalisedEnzymes = self._lastNeofunctionalisedEnzymesCache
            
            if lastMajorityPercentage == majorityPercentageCoreMetabolism:
                
                return lastNeofunctionalisedEnzymes
            
            else:
                self._lastNeofunctionalisedEnzymesCache = None
        
        # calculate
        enzymes = self.coreMetabolismEnzymes(majorityPercentageCoreMetabolism)
        
        if considerOnlyECs is not None:
            
            enzymes.keepEnzymesByEC(considerOnlyECs)

        enzymes = enzymes.getEnzymes()
            
            
        geneDuplicationModel = SimpleGeneDuplication
#             geneDuplicationModel = SimpleGroupGeneDuplication(sameGroupOrganisms = self.group)
        neofunctionalisedEnzymes = NeofunctionalisedEnzymes(enzymes, geneDuplicationModel, eValue = eValue)
        
        # Cache calculation
        if considerOnlyECs is None:
            self._lastNeofunctionalisedEnzymesCache = (majorityPercentageCoreMetabolism, neofunctionalisedEnzymes)
        
        return neofunctionalisedEnzymes
        
    
    def neofunctionalisedEnzymes(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, colour = False, eValue = defaultEValue, considerOnlyECs = None) -> SubstanceEnzymeGraph:
        """
        The substance-Enzyme graph of all neofunctionalised enzymes of the core metabolism.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            Every substance-EC-product edge has to occur in `majorityPercentageCoreMetabolism` % of organisms constituting the clade, to be included in the core metabolism.
        colour : bool, optional
            If *True*, colours the neofunctionalised enzyme edges in green. The colouring is realised by adding a 'colour' attribute to each edge. Nodes are not coloured.
            Alternatively, you can specify a :class:`Export.Colour`.
        eValue : float, optional
            Threshold for the statistical expectation value (E-value), below which a sequence alignment is considered significant.
        considerOnlyECs : Iterable[EcNumber], optional
            If given, only enzymes with an EC number in `considerOnlyECs` are tested for neofunctionalisation.
        
        Returns
        -------
        SubstanceEnzymeGraph
            Substance-Enzyme graph containing all neofunctionalised enzymes, and nothing else.
            If `colour` == *True*, returns the full core metabolism enzyme graph, colouring neofunctionalised enzymes green.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        # get neofunctionalisations        
        neofunctionalisedEnzymes = self._neofunctionalisedEnzymes(majorityPercentageCoreMetabolism, eValue, considerOnlyECs)
        
        # filter core metabolism enzyme graph
        enzymeGraph = self.coreMetabolismEnzymes(majorityPercentageCoreMetabolism)        
        neofunctionalisedMetabolism = neofunctionalisedEnzymes.filterGraph(enzymeGraph, minimumEcDifference = None)
        
        # colour core metabolism            
        if colour is not False:
            
            if colour is True:
                colourToUse = Export.Colour.GREEN
            else:
                colourToUse = colour
            
            neofunctionalisedMetabolismOnly = neofunctionalisedMetabolism
            neofunctionalisedMetabolism = enzymeGraph
            Export.addColourAttribute(neofunctionalisedMetabolism, colourToUse, nodes = False, edges = neofunctionalisedMetabolismOnly.getEdges())
        
        neofunctionalisedMetabolism.name = 'Neofunctionalised core metabolism enzymes ' + ' '.join(self.ncbiNames)
        
        return neofunctionalisedMetabolism
    
    
    def neofunctionalisedECs(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation = defaultMajorityPercentageNeofunctionalisation, colour = False, eValue = defaultEValue, considerOnlyECs = None) -> SubstanceEcGraph:
        """
        The substance-EC graph of EC numbers belonging to function changes of neofunctionalised enzymes.
        
        Only EC numbers which could have actually taken part in a function change are reported. This is because enzymes can have multiple EC numbers, while only some might be eligible for a function change.
        For example, consider enzyme A (1.2.3.4, 6.5.4.3) and enzyme B (1.2.3.4, 4.5.6.7). 1.2.3.4 can never change its function to itself, which leaves 1.2.3.4 <-> 6.5.4.3, 1.2.3.4 <-> 4.5.6.7, and 4.5.6.7 <-> 6.5.4.3 as possible function changes.
        This obviously requires a function to change to a single other function, without splitting or merging, which might be biologically inacurate. However, this should happen rarely and you can exclude all enzymes with multiple functions from the core metabolism in the first place.
        
        The maximum expectation value (e-value) necessary for a sequence alignment to constitute a "similar sequence" can be changed via :attr:`defaultEValue`.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            Every substance-EC-product edge has to occur in `majorityPercentageCoreMetabolism` % of organisms constituting the clade, to be included in the core metabolism. 
        majorityPercentageNeofunctionalisation : int, optional
            Every EC number considered for neofunctionalisation has to be associated with a function change of neofunctionalisations whose enzymes involve at least `majorityPercentageNeofunctionalisation` % of of the clade's organisms.
            A high `majorityPercentageNeofunctionalisation` disallows us to detect neofunctionalisations which happened a long time ago, with their genes having diverged significantly; 
            or only recently, with not all organisms of the child clade having picked up the new function, yet.
        colour : bool, optional
            If *True*, colours the neofunctionalised EC edges in green. The colouring is realised by adding a 'colour' attribute to each edge. Nodes are not coloured.
            Alternatively, you can specify a :class:`Export.Colour`.
        eValue : float, optional
            Threshold for the statistical expectation value (E-value), below which a sequence alignment is considered significant.
        considerOnlyECs : Iterable[EcNumber], optional
            If given, only enzymes with an EC number in `considerOnlyECs` are tested for neofunctionalisation.
        
        Returns
        -------
        SubstanceEcGraph
            The substance-EC graph representing the metabolic network which was probably affected due to neofunctionalisations of the core metabolism of the clade.
            If `colour` == *True*, returns the full union of parent and child, colouring neofunctionalised ECs green.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        # get neofunctionalisations        
        neofunctionalisedECs = NeofunctionalisedECs(self._neofunctionalisedEnzymes(majorityPercentageCoreMetabolism, eValue, considerOnlyECs))
        
        # filter core metabolism EC graph
        coreMetabolism = self.coreMetabolism(majorityPercentageCoreMetabolism)
        minimumOrganismsCount = math.ceil(self.organismsCount * (majorityPercentageNeofunctionalisation / 100))
        
        neofunctionalisedMetabolism = neofunctionalisedECs.filterGraph(coreMetabolism, minimumEcDifference = None, minimumOrganismsCount = minimumOrganismsCount)
        
        # colour core metabolism
        if colour is not False:
            
            if colour is True:
                colourToUse = Export.Colour.GREEN
            else:
                colourToUse = colour
                
            neofunctionalisedMetabolismOnly = neofunctionalisedMetabolism
            neofunctionalisedMetabolism = coreMetabolism
            Export.addColourAttribute(neofunctionalisedMetabolism, colourToUse, nodes = False, edges = neofunctionalisedMetabolismOnly.getEdges())
        
        neofunctionalisedMetabolism.name = 'Neofunctionalised core metabolism ' + ' '.join(self.ncbiNames)
        
        return neofunctionalisedMetabolism
    
    
    def neofunctionalisations(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, eValue = defaultEValue, considerOnlyECs = None) -> Set[Neofunctionalisation]:
        """
        Get neofunctionalisation events of all enzymes in the core metabolism.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            Every substance-EC-product edge has to occur in `majorityPercentageCoreMetabolism` % of organisms constituting the clade, to be included in the core metabolism.
        eValue : float, optional
            Threshold for the statistical expectation value (E-value), below which a sequence alignment is considered significant.
        considerOnlyECs : Iterable[EcNumber], optional
            If given, only enzymes with an EC number in `considerOnlyECs` are tested for neofunctionalisation.
        
        Returns
        -------
        Set[Neofunctionalisation]
            Set of possible neofunctionalisation events.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        # get neofunctionalisations        
        return self._neofunctionalisedEnzymes(majorityPercentageCoreMetabolism, eValue, considerOnlyECs).getNeofunctionalisations()
    
    
    def neofunctionalisationsForFunctionChange(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation = defaultMajorityPercentageNeofunctionalisation, eValue = defaultEValue, considerOnlyECs = None) -> Dict[FunctionChange, Set[Neofunctionalisation]]:
        """
        Get neofunctionalisation events of all enzymes in the core metabolism, grouped by each possible function change event.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            Every substance-EC-product edge has to occur in `majorityPercentageCoreMetabolism` % of organisms constituting the clade, to be included in the core metabolism. 
        majorityPercentageNeofunctionalisation : int, optional
            Every EC number considered for neofunctionalisation has to be associated with a function change of neofunctionalisations whose enzymes involve at least `majorityPercentageNeofunctionalisation` % of of the clade's organisms.
            A high `majorityPercentageNeofunctionalisation` disallows us to detect neofunctionalisations which happened a long time ago, with their genes having diverged significantly; 
            or only recently, with not all organisms of the child clade having picked up the new function, yet.
        eValue : float, optional
            Threshold for the statistical expectation value (E-value), below which a sequence alignment is considered significant.
        considerOnlyECs : Iterable[EcNumber], optional
            If given, only enzymes with an EC number in `considerOnlyECs` are tested for neofunctionalisation.
            
        Returns
        -------
        Dict[FunctionChange, Set[Neofunctionalisation]]
            Dictionary of function changes, pointing to a set of neofunctionalisations which might have caused them.
            
            Since an enzyme of a neofunctionalisation can have multiple EC numbers, all combinations of the two enzymes' EC numbers are formed and treated as separate possible function changes.
            The neofunctionalisation is then saved again for each function change, which obviously leads to duplicated neofunctionalisation objects.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        # get neofunctionalisations
        minimumOrganismsCount = math.ceil(self.organismsCount * (majorityPercentageNeofunctionalisation / 100))        
        return NeofunctionalisedECs(self._neofunctionalisedEnzymes(majorityPercentageCoreMetabolism, eValue, considerOnlyECs)).getNeofunctionalisationsForFunctionChange(minimumOrganismsCount = minimumOrganismsCount)

    
    
    
    
    # redundancy of neofunctionalisation
    
    def redundantECsForContributingNeofunctionalisation(self, 
                                                        majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, 
                                                        majorityPercentageNeofunctionalisation = defaultMajorityPercentageNeofunctionalisation, 
                                                        eValue = defaultEValue, 
                                                        redundancyType: 'RedundancyType' = None,
                                                        considerOnlyECs = None) -> Dict[Neofunctionalisation, Set[EcNumber]]:
        """
        Get neofunctionalisation events of all enzymes in the core metabolism, which contribute to redundancy, pointing to the EC numbers their function changes' EC numbers provides redundancy for.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            Every substance-EC-product edge has to occur in `majorityPercentageCoreMetabolism` % of organisms constituting the clade, to be included in the core metabolism. 
        majorityPercentageNeofunctionalisation : int, optional
            Every EC number considered for neofunctionalisation has to be associated with a function change of neofunctionalisations whose enzymes involve at least `majorityPercentageNeofunctionalisation` % of of the clade's organisms.
            A high `majorityPercentageNeofunctionalisation` disallows us to detect neofunctionalisations which happened a long time ago, with their genes having diverged significantly; 
            or only recently, with not all organisms of the child clade having picked up the new function, yet.
        eValue : float, optional
            Threshold for the statistical expectation value (E-value), below which a sequence alignment is considered significant.
        redundancyType : RedundancyType
            Definition of redundancy for which to check the neofunctionalisation's contribution. Default to `RedundancyType.default`.
        considerOnlyECs : Iterable[EcNumber], optional
            If given, only enzymes with an EC number in `considerOnlyECs` are tested for neofunctionalisation.
            
        Returns
        -------
        Dict[FunctionChange, Set[Neofunctionalisation]]
            Dictionary of function changes, pointing to a set of neofunctionalisations which might have caused them.
            
            Since an enzyme of a neofunctionalisation can have multiple EC numbers, all combinations of the two enzymes' EC numbers are formed and treated as separate possible function changes.
            The neofunctionalisation is then saved again for each function change, which obviously leads to duplicated neofunctionalisation objects.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        from FEV_KEGG.Robustness.Topology.Redundancy import Redundancy, RedundancyContribution, RedundancyType
        
        if redundancyType is None:
            redundancyType = RedundancyType.default
        
        #- calculate "neofunctionalised" ECs
        neofunctionalisedMetabolismSet = self.neofunctionalisedECs(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation, eValue, considerOnlyECs).getECs()
        neofunctionalisationsForFunctionChange = self.neofunctionalisationsForFunctionChange(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation, eValue, considerOnlyECs)
        
        #- calculate redundancy
        redundancy = Redundancy( self.coreMetabolism(majorityPercentageCoreMetabolism) )
        redundancyContribution = RedundancyContribution(redundancy, neofunctionalisedMetabolismSet)
        
        contributedECsForContributingNeofunctionalisedEC = redundancyContribution.getContributedKeysForSpecial(redundancyType)
        contributingNeofunctionalisedECs = set(contributedECsForContributingNeofunctionalisedEC.keys())
        
        #- REPEAT for each function change consisting of "neofunctionalised" ECs, which also contribute to redundancy
        contributingNeofunctionalisations = dict()
        
        for functionChange, neofunctionalisations in neofunctionalisationsForFunctionChange.items():
            #-     report enzyme pairs of neofunctionalisations, which caused the EC to be considered "neofunctionalised", and are in return contributing to redundancy        
            
            if functionChange.ecA in contributingNeofunctionalisedECs or functionChange.ecB in contributingNeofunctionalisedECs: # function change contributes to redundancy
                
                for neofunctionalisation in neofunctionalisations:
                    currentSetOfContributedECs = contributingNeofunctionalisations.get(neofunctionalisation, None)
                    
                    if currentSetOfContributedECs is None:
                        currentSetOfContributedECs = set()
                        contributingNeofunctionalisations[neofunctionalisation] = currentSetOfContributedECs
                    
                    for ec in functionChange.ecPair:
                        contributedECs = contributedECsForContributingNeofunctionalisedEC.get(ec, None)
                        if contributedECs is not None:
                            currentSetOfContributedECs.update(contributedECs)
        
        return contributingNeofunctionalisations
        
    
    







class CladePair(object):
    
    def __init__(self, parent, child, excludeUnclassified = defaultExcludeUnclassified, oneOrganismPerSpecies = defaultOneOrganismPerSpecies):
        """
        Two clades in NCBI taxonomy, 'child' is assumed younger than 'parent'.
        
        Does not check if the child taxon is actually a child of the parent taxon.
        Therefore, it would be possible to pass a list of NCBI names to the underlying :class:`Clade` objects by instantiating `parent` = List[str] and/or `child` = List[str].
        This is useful when comparing groups of organisms which are, according to NCBI, not related.
        
        Parameters
        ----------
        parent : str or List[str] or Clade
            Path(s) of the parent clade's taxon, as defined by NCBI taxonomy, e.g. 'Proteobacteria/Gammaproteobacteria'. Or a ready :class:`Clade` object.
        child : str or List[str] or Clade
            Path(s) of the child clade's taxon, as defined by NCBI taxonomy, e.g. 'Enterobacter'. Or a ready :class:`Clade` object.
        excludeUnclassified : bool, optional
            If *True*, ignore taxons with a path containing the string 'unclassified'. Only used if one of `parent` and/or `child` is not already a :class:`Clade`.
        oneOrganismPerSpecies : bool, optional
            If *True*, use only the first organism of each species.
        
        Attributes
        ----------
        self.childClade : :class:`Clade`
        self.parentClade : :class:`Clade`
        """
        # read NCBI names from Clade object, if necessary
        if isinstance(parent, Clade):
            self.parentClade = parent
        else:
            self.parentClade = Clade(parent, excludeUnclassified, oneOrganismPerSpecies=oneOrganismPerSpecies)
        
        if isinstance(child, Clade):
            self.childClade = child
        else:
            self.childClade = Clade(child, excludeUnclassified, oneOrganismPerSpecies=oneOrganismPerSpecies)
    
    
    @property
    def parentNCBInames(self):
        """
        All names/paths in NCBI taxonomy used to create the parent clade.
        """
        return self.parentClade.ncbiNames
    
    @property
    def childNCBInames(self):
        """
        All names/paths in NCBI taxonomy used to create the child clade.
        """
        return self.childClade.ncbiNames
    
    
    
    
    
    # set-operations on core metabolism
    ## for EC graphs
    def conservedMetabolism(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism) -> SubstanceEcGraph:
        """
        Substance-EC graph of the conserved core metabolism.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            Every substance-EC-product edge has to occur in `majorityPercentageCoreMetabolism` % of organisms constituting the clade, to be included in the core metabolism. This is individually true for both parent clade and child clade.
            The parent clade fully includes the child clade, therefore, the occurence of a substance-EC-product edge in the child clade's core metabolism counts towards the percentage for the parent clade's core metabolism.
            Meaning: if an EC number does not occur in the child clade's core metabolism, it is unlikely that it will occur in the parent clade's core metabolism, unless `majorityPercentageCoreMetabolism` is consecutively lowered towards 0.
        
        Returns
        -------
        SubstanceEcGraph
            The substance-EC graph representing the metabolic network which stayed the same between the core metabolism of the parent (assumed older) and the core metabolism of the child (assumed younger).
            
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        parentCoreMetabolism = self.parentClade.coreMetabolism(majorityPercentageCoreMetabolism)
        childCoreMetabolism = self.childClade.coreMetabolism(majorityPercentageCoreMetabolism)
        graph = GeneFunctionConservation.getGraph(parentCoreMetabolism, childCoreMetabolism)
        graph.name = 'Conserved metabolism ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
        return graph
    
    
    def addedMetabolism(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism) -> SubstanceEcGraph:
        """
        Substance-EC graph of the added core metabolism.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            Every substance-EC-product edge has to occur in `majorityPercentageCoreMetabolism` % of organisms constituting the clade, to be included in the core metabolism. This is individually true for both parent clade and child clade.
            The parent clade fully includes the child clade, therefore, the occurence of a substance-EC-product edge in the child clade's core metabolism counts towards the percentage for the parent clade's core metabolism.
            Meaning: if an EC number does not occur in the child clade's core metabolism, it is unlikely that it will occur in the parent clade's core metabolism, unless `majorityPercentageCoreMetabolism` is consecutively lowered towards 0.
        
        Returns
        -------
        SubstanceEcGraph
            The substance-EC graph representing the metabolic network which was added to the core metabolism of the parent (assumed older) on the way to the core metabolism of the child (assumed younger).
            
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        parentCoreMetabolism = self.parentClade.coreMetabolism(majorityPercentageCoreMetabolism)
        childCoreMetabolism = self.childClade.coreMetabolism(majorityPercentageCoreMetabolism)
        graph = GeneFunctionAddition.getGraph(parentCoreMetabolism, childCoreMetabolism)
        graph.name = 'Added metabolism ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
        return graph
    
    
    def lostMetabolism(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism) -> SubstanceEcGraph:
        """
        Substance-EC graph of the lost core metabolism.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            Every substance-EC-product edge has to occur in `majorityPercentageCoreMetabolism` % of organisms constituting the clade, to be included in the core metabolism. This is individually true for both parent clade and child clade.
            The parent clade fully includes the child clade, therefore, the occurence of a substance-EC-product edge in the child clade's core metabolism counts towards the percentage for the parent clade's core metabolism.
            Meaning: if an EC number does not occur in the child clade's core metabolism, it is unlikely that it will occur in the parent clade's core metabolism, unless `majorityPercentageCoreMetabolism` is consecutively lowered towards 0.
        
        Returns
        -------
        SubstanceEcGraph
            The substance-EC graph representing the metabolic network which got lost from the core metabolism of the parent (assumed older) on the way to the core metabolism of the child (assumed younger).
            
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        parentCoreMetabolism = self.parentClade.coreMetabolism(majorityPercentageCoreMetabolism)
        childCoreMetabolism = self.childClade.coreMetabolism(majorityPercentageCoreMetabolism)        
        graph = GeneFunctionLoss.getGraph(parentCoreMetabolism, childCoreMetabolism)
        graph.name = 'Lost metabolism ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
        return graph
    
    
    def divergedMetabolism(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, colour = False) -> SubstanceEcGraph:
        """
        Substance-EC graph of the diverged core metabolism.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            Every substance-EC-product edge has to occur in `majorityPercentageCoreMetabolism` % of organisms constituting the clade, to be included in the core metabolism. This is individually true for both parent clade and child clade.
            The parent clade fully includes the child clade, therefore, the occurence of a substance-EC-product edge in the child clade's core metabolism counts towards the percentage for the parent clade's core metabolism.
            Meaning: if an EC number does not occur in the child clade's core metabolism, it is unlikely that it will occur in the parent clade's core metabolism, unless `majorityPercentageCoreMetabolism` is consecutively lowered towards 0.
        colour : bool, optional
            If *True*, colours the lost EC edges in blue, and the added EC edges in red. The colouring is realised by adding a 'colour' attribute to each edge. Nodes are not coloured.
        
        Returns
        -------
        SubstanceEcGraph
            The substance-EC graph representing the metabolic network which changed between the core metabolism of the parent (assumed older) and the core metabolism of the child (assumed younger).
            
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        parentCoreMetabolism = self.parentClade.coreMetabolism(majorityPercentageCoreMetabolism)
        childCoreMetabolism = self.childClade.coreMetabolism(majorityPercentageCoreMetabolism)
        
        if colour is True:
            lostGraph = GeneFunctionLoss.getGraph(parentCoreMetabolism, childCoreMetabolism)
            lostEdges = lostGraph.getEdges()
            
            addedGraph = GeneFunctionAddition.getGraph(parentCoreMetabolism, childCoreMetabolism)
            addedEdges = addedGraph.getEdges()
            
            graph = lostGraph.union(addedGraph, addCount = False, updateName = False) 
            
            Export.addColourAttribute(graph, colour = Export.Colour.BLUE, nodes = False, edges = lostEdges)
            Export.addColourAttribute(graph, colour = Export.Colour.RED, nodes = False, edges = addedEdges)
            
        else:       
            graph = GeneFunctionDivergence.getGraph(parentCoreMetabolism, childCoreMetabolism)
        
        graph.name = 'Diverged metabolism ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
            
        return graph
    
    
    def unifiedMetabolism(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, colour = False) -> SubstanceEcGraph:
        """
        Substance-EC graph of the unified core metabolisms.
        
        The lost metabolism of the parent is coloured in blue, the conserved metabolism of both in red, and the added metabolism of the child in pink.
        The colouring is realised by adding a 'colour' attribute to each edge. Nodes are not coloured.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`conservedMetabolism`.
        colour : bool, optional
            If *True*, colours the parent's EC edges in blue, the child's EC edges in red, and the shared EC edges in pink. The colouring is realised by adding a 'colour' attribute to each edge. Nodes are not coloured.
        
        Returns
        -------
        SubstanceEcGraph
            The substance-EC graph representing the combined metabolic networks of both, child and parent. If `colour` == *True*, coloured differently for the lost, conserved, and added edges. Nodes are not coloured.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        
        See Also
        --------
        :mod:`FEV_KEGG.Drawing.Export` : Export the graph into a file, e.g. for visualisation in Cytoscape.
        """
        parentCoreMetabolism = self.parentClade.coreMetabolism(majorityPercentageCoreMetabolism)
        childCoreMetabolism = self.childClade.coreMetabolism(majorityPercentageCoreMetabolism)
        
        graph = parentCoreMetabolism.union(childCoreMetabolism, addCount = False, updateName = False)
        graph.name = 'Unified metabolism ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
        
        if colour is True:
            lostGraph = GeneFunctionLoss.getGraph(parentCoreMetabolism, childCoreMetabolism)
            lostEdges = lostGraph.getEdges()
            
            addedGraph = GeneFunctionAddition.getGraph(parentCoreMetabolism, childCoreMetabolism)
            addedEdges = addedGraph.getEdges()
            
            conservedGraph = GeneFunctionConservation.getGraph(parentCoreMetabolism, childCoreMetabolism)
            conservedEdges = conservedGraph.getEdges()            
            
            Export.addColourAttribute(graph, colour = Export.Colour.BLUE, nodes = False, edges = lostEdges)
            Export.addColourAttribute(graph, colour = Export.Colour.RED, nodes = False, edges = addedEdges)
            Export.addColourAttribute(graph, colour = Export.Colour.PINK, nodes = False, edges = conservedEdges)
            
        return graph
    
    
    
    
    
    
    ## for enzyme graphs
    def conservedMetabolismEnzymes(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, colour = False):
        """
        Two Substance-Enzyme graphs derived from the conserved core metabolism, see :func:`conservedMetabolism`.
        
        First, the conserved core metabolism is calculated. Then, the enzymes associated with the conserved EC numbers are extracted from the collective parent's and child's metabolism individually.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`conservedMetabolism`.
        colour : bool, optional
            If *True*, colours the enzyme edges from the parent in blue, and from the child in red. When doing so, a single :class:`SubstanceEnzymeGraph` is returned, not a :class:`Tuple`. The colouring is realised by adding a 'colour' attribute to each edge. Nodes are not coloured.
        
        Returns
        -------
        Tuple[SubstanceEnzymeGraph, SubstanceEnzymeGraph] or SubstanceEnzymeGraph
            Tuple of two Substance-Enzyme graphs calculated using the conserved EC numbers found by :func:`conservedMetabolism`. The first graph is from the parent clade, the second graph from the child clade.
            If `colour` == *True*, returns a single Substance-Enzyme graph, coloured blue for parent and red for child.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        parentCoreMetabolism = self.parentClade.coreMetabolism(majorityPercentageCoreMetabolism)
        childCoreMetabolism = self.childClade.coreMetabolism(majorityPercentageCoreMetabolism)
        conservedECs = GeneFunctionConservation.getECs(parentCoreMetabolism, childCoreMetabolism)
        
        parentGraph = self.parentClade.collectiveMetabolismEnzymes().keepEnzymesByEC(conservedECs)        
        childGraph = self.childClade.collectiveMetabolismEnzymes().keepEnzymesByEC(conservedECs)    
    
        if colour is True:
            parentEdges = parentGraph.getEdges()
            childEdges = childGraph.getEdges()
            
            graph = parentGraph.union(childGraph, addCount = False, updateName = False)
            
            Export.addColourAttribute(graph, colour = Export.Colour.BLUE, nodes = False, edges = parentEdges)
            Export.addColourAttribute(graph, colour = Export.Colour.RED, nodes = False, edges = childEdges)
            
            graph.name = 'Conserved metabolism enzymes ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
            
            return graph
        else:
            parentGraph.name = 'Conserved metabolism enzymes *' + ' '.join(self.parentNCBInames) + '* -> ' + ' '.join(self.childNCBInames)
            childGraph.name = 'Conserved metabolism enzymes ' + ' '.join(self.parentNCBInames) + ' -> *' + ' '.join(self.childNCBInames) + '*'
        
            return (parentGraph, childGraph)
    
    
    def addedMetabolismEnzymes(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism) -> SubstanceEnzymeGraph:
        """
        Substance-Enzyme graph derived from the added core metabolism, see :func:`addedMetabolism`.
        
        First, the added core metabolism is calculated. Then, the enzymes associated with the added EC numbers are extracted from the child's enzyme metabolism.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`addedMetabolism`.
        
        Returns
        -------
        SubstanceEnzymeGraph
            Substance-Enzyme graph of enzymes from the child clade. Calculated using the added EC numbers found by :func:`addedMetabolism`.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        parentCoreMetabolism = self.parentClade.coreMetabolism(majorityPercentageCoreMetabolism)
        childCoreMetabolism = self.childClade.coreMetabolism(majorityPercentageCoreMetabolism)
        addedECs = GeneFunctionAddition.getECs(parentCoreMetabolism, childCoreMetabolism)
        
        childGraph = self.childClade.collectiveMetabolismEnzymes().keepEnzymesByEC(addedECs)
        childGraph.name = 'Added metabolism enzymes ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
        
        return childGraph
    
    
    def lostMetabolismEnzymes(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism) -> SubstanceEnzymeGraph:
        """
        Substance-Enzyme graph derived from the lost core metabolism, see :func:`lostMetabolism`.
        
        First, the lost core metabolism is calculated. Then, the enzymes associated with the added EC numbers are extracted from the parent's enzyme metabolism.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`lostMetabolism`.
        
        Returns
        -------
        SubstanceEnzymeGraph
            Substance-Enzyme graph of enzymes from the parent clade. Calculated using the lost EC numbers found by :func:`lostMetabolism`.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        parentCoreMetabolism = self.parentClade.coreMetabolism(majorityPercentageCoreMetabolism)
        childCoreMetabolism = self.childClade.coreMetabolism(majorityPercentageCoreMetabolism)
        lostECs = GeneFunctionLoss.getECs(parentCoreMetabolism, childCoreMetabolism)
        
        parentGraph = self.parentClade.collectiveMetabolismEnzymes().keepEnzymesByEC(lostECs)
        parentGraph.name = 'Lost metabolism enzymes ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
        
        return parentGraph
    
    
    def divergedMetabolismEnzymes(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, colour = False):
        """
        Two Substance-Enzyme graphs derived from the diverged core metabolism, see :func:`divergedMetabolism`.
        
        First, the diverged core metabolism is calculated. Then, the enzymes associated with the added EC numbers are extracted from the collective parent's and child's metabolism individually.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`divergedMetabolism`.
        colour : bool, optional
            If *True*, colours the lost enzyme edges in blue, and the added enzyme edges in red. When doing so, a single :class:`SubstanceEnzymeGraph` is returned, not a :class:`Tuple`. The colouring is realised by adding a 'colour' attribute to each edge. Nodes are not coloured.
        
        Returns
        -------
        Tuple[SubstanceEnzymeGraph, SubstanceEnzymeGraph] or SubstanceEnzymeGraph
            Tuple of two Substance-Enzyme graphs calculated using the diverged EC numbers found by :func:`divergedMetabolism`. The first graph is from the parent clade, the second graph from the child clade.
            If `colour` == *True*, returns a single Substance-Enzyme graph, coloured blue for parent and red for child.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        parentCoreMetabolism = self.parentClade.coreMetabolism(majorityPercentageCoreMetabolism)
        childCoreMetabolism = self.childClade.coreMetabolism(majorityPercentageCoreMetabolism)
        divergedECs = GeneFunctionDivergence.getECs(parentCoreMetabolism, childCoreMetabolism)
        
        parentGraph = self.parentClade.collectiveMetabolismEnzymes().keepEnzymesByEC(divergedECs)        
        childGraph = self.childClade.collectiveMetabolismEnzymes().keepEnzymesByEC(divergedECs)
        
        if colour is True:
            parentEdges = parentGraph.getEdges()
            childEdges = childGraph.getEdges()
            
            graph = parentGraph.union(childGraph, addCount = False, updateName = False)
            
            Export.addColourAttribute(graph, colour = Export.Colour.BLUE, nodes = False, edges = parentEdges)
            Export.addColourAttribute(graph, colour = Export.Colour.RED, nodes = False, edges = childEdges)
            
            graph.name = 'Diverged metabolism enzymes ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
            
            return graph
        else:
            parentGraph.name = 'Diverged metabolism enzymes *' + ' '.join(self.parentNCBInames) + '* -> ' + ' '.join(self.childNCBInames)
            childGraph.name = 'Diverged metabolism enzymes ' + ' '.join(self.parentNCBInames) + ' -> *' + ' '.join(self.childNCBInames) + '*'
        
            return (parentGraph, childGraph)
    
    
    def unifiedMetabolismEnzymes(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, colour = False) -> SubstanceEnzymeGraph:
        """
        Substance-Enzyme graph derived from the unified core metabolisms.
        
        The lost metabolism of the parent is coloured in blue, the conserved metabolism of both in red, and the added metabolism of the child in pink.
        The colouring is realised by adding a 'colour' attribute to each edge. Nodes are not coloured.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`conservedMetabolism`.
        colour : bool, optional
            If *True*, colours the parent's enzyme edges in blue, and the child's enzyme edges in red. The colouring is realised by adding a 'colour' attribute to each edge. Nodes are not coloured.
        
        Returns
        -------
        SubstanceEnzymeGraph
            The substance-Enzyme graph representing the combined metabolic networks of both, child and parent. If `colour` == *True*, coloured differently for the lost, conserved, and added edges. Nodes are not coloured.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        parentGraph = self.parentClade.coreMetabolismEnzymes(majorityPercentageCoreMetabolism)
        childGraph = self.childClade.coreMetabolismEnzymes(majorityPercentageCoreMetabolism)
        
        graph = parentGraph.union(childGraph, addCount = False, updateName = False)
        graph.name = 'Unified metabolism enzymes ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
        
        if colour is True:
            parentEdges = parentGraph.getEdges()
            childEdges = childGraph.getEdges()
            
            Export.addColourAttribute(graph, colour = Export.Colour.BLUE, nodes = False, edges = parentEdges)
            Export.addColourAttribute(graph, colour = Export.Colour.RED, nodes = False, edges = childEdges)
                
        return graph
    
    
    
    
    
    
    
    # set-operations on gene-duplicated core metabolism
    ## for enzymes
    ### for enzyme graphs
    def conservedMetabolismGeneDuplicatedEnzymes(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, colour = False):
        """
        Two Substance-Enzyme graphs of gene-duplicated enzymes, derived from the conserved core metabolism.
        
        First, the conserved core metabolism is calculated. Then, the enzymes associated with the conserved EC numbers are extracted from the collective parent's and child's metabolism individually.
        Then, for parent and child, the gene-duplicated enzymes are calculated. Finally, the gene-duplicated enzymes of the conserved core metabolism enzymes are reported.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`conservedMetabolism`.
        colour : bool, optional
            If *True*, colours the enzyme edges from the parent in blue, and from the child in red. Gene-duplicated enzyme edges of the parent are coloured in green, the ones of the child in yellow.
            When doing so, a single :class:`SubstanceEnzymeGraph` is returned, not a :class:`Tuple`. The colouring is realised by adding a 'colour' attribute to each edge. Nodes are not coloured.
        
        Returns
        -------
        Tuple[SubstanceEnzymeGraph, SubstanceEnzymeGraph] or SubstanceEnzymeGraph
            Tuple of two Substance-Enzyme graphs calculated using the conserved EC numbers found by :func:`conservedMetabolism`. The first graph is from the parent clade, the second graph from the child clade.
            If `colour` == *True*, returns a single Substance-Enzyme graph.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        conservedMetabolismEnzymes = self.conservedMetabolismEnzymes(majorityPercentageCoreMetabolism, colour = colour)
        
        parentGeneDuplicated = self.parentClade.geneDuplicatedEnzymes(majorityPercentageCoreMetabolism, colour = False)
        childGeneDuplicated = self.childClade.geneDuplicatedEnzymes(majorityPercentageCoreMetabolism, colour = False)
        
        if colour is True:
            parentEdges = parentGeneDuplicated.getEdges()
            childEdges = childGeneDuplicated.getEdges()
            
            graph = conservedMetabolismEnzymes
            
            Export.addColourAttribute(graph, colour = Export.Colour.GREEN, nodes = False, edges = parentEdges)
            Export.addColourAttribute(graph, colour = Export.Colour.YELLOW, nodes = False, edges = childEdges)
            
            graph.name = 'Conserved metabolism gene-duplicated enzymes ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
            
            return graph
        else:
            parentGraph = conservedMetabolismEnzymes[0].removeAllEnzymesExcept(parentGeneDuplicated.getEnzymes())
            childGraph = conservedMetabolismEnzymes[1].removeAllEnzymesExcept(childGeneDuplicated.getEnzymes())
            
            parentGraph.name = 'Conserved metabolism gene-duplicated enzymes *' + ' '.join(self.parentNCBInames) + '* -> ' + ' '.join(self.childNCBInames)
            childGraph.name = 'Conserved metabolism gene-duplicated enzymes ' + ' '.join(self.parentNCBInames) + ' -> *' + ' '.join(self.childNCBInames) + '*'
        
            return (parentGraph, childGraph)
    
    
    def addedMetabolismGeneDuplicatedEnzymes(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism) -> SubstanceEnzymeGraph:
        """
        Substance-Enzyme graph of gene-duplicated enzymes, derived from the added core metabolism.
        
        First, the added core metabolism is calculated. Then, the enzymes associated with the added EC numbers are extracted from the child's enzyme metabolism.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`addedMetabolism`.
        
        Returns
        -------
        SubstanceEnzymeGraph
            Substance-Enzyme graph of enzymes from the child clade. Calculated using the added EC numbers found by :func:`addedMetabolism`.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        parentCoreMetabolism = self.parentClade.coreMetabolism(majorityPercentageCoreMetabolism)
        childCoreMetabolism = self.childClade.coreMetabolism(majorityPercentageCoreMetabolism)
        addedECs = GeneFunctionAddition.getECs(parentCoreMetabolism, childCoreMetabolism)
        
        childGraph = self.childClade.geneDuplicatedEnzymes(majorityPercentageCoreMetabolism, colour = False).keepEnzymesByEC(addedECs)
        childGraph.name = 'Added metabolism gene-duplicated enzymes ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
        
        return childGraph
    
    
    def lostMetabolismGeneDuplicatedEnzymes(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism) -> SubstanceEnzymeGraph:
        """
        Substance-Enzyme graph of gene-duplicated enzymes, derived from the lost core metabolism.
        
        First, the lost core metabolism is calculated. Then, the enzymes associated with the added EC numbers are extracted from the parent's enzyme metabolism.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`lostMetabolism`.
        
        Returns
        -------
        SubstanceEnzymeGraph
            Substance-Enzyme graph of enzymes from the parent clade. Calculated using the lost EC numbers found by :func:`lostMetabolism`.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        parentCoreMetabolism = self.parentClade.coreMetabolism(majorityPercentageCoreMetabolism)
        childCoreMetabolism = self.childClade.coreMetabolism(majorityPercentageCoreMetabolism)
        lostECs = GeneFunctionLoss.getECs(parentCoreMetabolism, childCoreMetabolism)
        
        parentGraph = self.parentClade.geneDuplicatedEnzymes(majorityPercentageCoreMetabolism, colour = False).keepEnzymesByEC(lostECs)
        parentGraph.name = 'Lost metabolism gene-duplicated enzymes ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
        
        return parentGraph
    
    
    def divergedMetabolismGeneDuplicatedEnzymes(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, colour = False):
        """
        Two Substance-Enzyme graphs of gene-duplicated enzymes, derived from the diverged core metabolism.
        
        First, the diverged core metabolism is calculated. Then, the enzymes associated with the added EC numbers are extracted from the collective parent's and child's metabolism individually.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`divergedMetabolism`.
        colour : bool, optional
            If *True*, colours the lost enzyme edges in blue, and the added enzyme edges in red. Gene-duplicated enzyme edges of the parent are coloured in green, the ones of the child in yellow.
            When doing so, a single :class:`SubstanceEnzymeGraph` is returned, not a :class:`Tuple`. The colouring is realised by adding a 'colour' attribute to each edge. Nodes are not coloured.
        
        Returns
        -------
        Tuple[SubstanceEnzymeGraph, SubstanceEnzymeGraph] or SubstanceEnzymeGraph
            Tuple of two Substance-Enzyme graphs calculated using the diverged EC numbers found by :func:`divergedMetabolism`. The first graph is from the parent clade, the second graph from the child clade.
            If `colour` == *True*, returns a single Substance-Enzyme graph, coloured blue for parent and red for child.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        divergedMetabolismEnzymes = self.divergedMetabolismEnzymes(majorityPercentageCoreMetabolism, colour = colour)
        
        parentGeneDuplicated = self.parentClade.geneDuplicatedEnzymes(majorityPercentageCoreMetabolism, colour = False)
        childGeneDuplicated = self.childClade.geneDuplicatedEnzymes(majorityPercentageCoreMetabolism, colour = False)
        
        if colour is True:
            parentEdges = parentGeneDuplicated.getEdges()
            childEdges = childGeneDuplicated.getEdges()
            
            graph = divergedMetabolismEnzymes
            
            Export.addColourAttribute(graph, colour = Export.Colour.GREEN, nodes = False, edges = parentEdges)
            Export.addColourAttribute(graph, colour = Export.Colour.YELLOW, nodes = False, edges = childEdges)
            
            graph.name = 'Diverged metabolism gene-duplicated enzymes ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
            
            return graph
        else:
            parentGraph = divergedMetabolismEnzymes[0].removeAllEnzymesExcept(parentGeneDuplicated.getEnzymes())
            childGraph = divergedMetabolismEnzymes[1].removeAllEnzymesExcept(childGeneDuplicated.getEnzymes())
            
            parentGraph.name = 'Diverged metabolism gene-duplicated enzymes *' + ' '.join(self.parentNCBInames) + '* -> ' + ' '.join(self.childNCBInames)
            childGraph.name = 'Diverged metabolism gene-duplicated enzymes ' + ' '.join(self.parentNCBInames) + ' -> *' + ' '.join(self.childNCBInames) + '*'
        
            return (parentGraph, childGraph)
    
    
    def unifiedMetabolismGeneDuplicatedEnzymes(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, colour = False) -> SubstanceEnzymeGraph:
        """
        Substance-Enzyme graph of gene-duplicated enzymes, derived from the unified core metabolisms.
        
        The lost metabolism of the parent is coloured in blue, the conserved metabolism of both in red, and the added metabolism of the child in pink.
        The colouring is realised by adding a 'colour' attribute to each edge. Nodes are not coloured.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`conservedMetabolism`.
        colour : bool, optional
            If *True*, colours the parent's enzyme edges in blue, and the child's enzyme edges in red. Gene-duplicated enzyme edges of the parent are coloured in green, the ones of the child in yellow.
            The colouring is realised by adding a 'colour' attribute to each edge. Nodes are not coloured.
        
        Returns
        -------
        SubstanceEcGraph
            The substance-Enzyme graph representing the combined metabolic networks of both, child and parent. If `colour` == *True*, coloured differently for the lost, conserved, and added edges. Nodes are not coloured.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """        
        parentGeneDuplicated = self.parentClade.geneDuplicatedEnzymes(majorityPercentageCoreMetabolism, colour = False)
        childGeneDuplicated = self.childClade.geneDuplicatedEnzymes(majorityPercentageCoreMetabolism, colour = False)
        
        if colour is False:
            graph = parentGeneDuplicated.union(childGeneDuplicated, addCount = False, updateName = False)
        
        else:
            unifiedMetabolismEnzymes = self.unifiedMetabolismEnzymes(majorityPercentageCoreMetabolism, colour = True)
            
            parentEdges = parentGeneDuplicated.getEdges()
            childEdges = childGeneDuplicated.getEdges()
            
            graph = unifiedMetabolismEnzymes
            
            Export.addColourAttribute(graph, colour = Export.Colour.GREEN, nodes = False, edges = parentEdges)
            Export.addColourAttribute(graph, colour = Export.Colour.YELLOW, nodes = False, edges = childEdges)
                
        return graph
    
    
    
    ### for enzyme pairs
    def conservedMetabolismGeneDuplicatedEnzymePairs(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism) -> Tuple[Set[Tuple[Enzyme, Enzyme]]]:
        """
        Pairs of gene-duplicated enzymes, derived from the conserved core metabolism.
        
        First, the conserved core metabolism is calculated. Then, the enzymes associated with the conserved EC numbers are extracted from the collective parent's and child's metabolism individually.
        Then, for parent and child, the gene-duplicated enzyme pairs are calculated. Finally, the gene-duplicated enzymes where both enzymes are in the conserved core metabolism are reported.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`conservedMetabolism`.
        
        Returns
        -------
        Tuple[Set[Tuple[Enzyme, Enzyme]]]
            Tuple of two sets of tuples of gene-duplicated enzyme pairs calculated using the conserved EC numbers found by :func:`conservedMetabolism`. The first set is from the parent clade, the second set from the child clade.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        # get conserved metabolism
        conservedMetabolismEnzymes = self.conservedMetabolismEnzymes(majorityPercentageCoreMetabolism).getEnzymes()
        
        # get gene-duplicate enzyme pairs
        parentGeneDuplicated = self.parentClade.geneDuplicatedEnzymePairs(majorityPercentageCoreMetabolism)
        childGeneDuplicated = self.childClade.geneDuplicatedEnzymePairs(majorityPercentageCoreMetabolism)
        
        # filter gene-duplicated enzyme pairs for the ones with both enzymes in the conserved metabolism
        parentGeneDuplicatedConserved = set()
        childGeneDuplicatedConserved = set()
        
        for enzymeTuple in parentGeneDuplicated:
            if enzymeTuple[0] in conservedMetabolismEnzymes and enzymeTuple[1] in conservedMetabolismEnzymes:
                parentGeneDuplicatedConserved.add(enzymeTuple)
        
        for enzymeTuple in childGeneDuplicated:
            if enzymeTuple[0] in conservedMetabolismEnzymes and enzymeTuple[1] in conservedMetabolismEnzymes:
                childGeneDuplicatedConserved.add(enzymeTuple)
        
        return (parentGeneDuplicatedConserved, childGeneDuplicatedConserved)
    
    
    def addedMetabolismGeneDuplicatedEnzymePairs(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism) -> Set[Tuple[Enzyme, Enzyme]]:
        """
        Pairs of gene-duplicated enzymes, derived from the added core metabolism.
        
        First, the added core metabolism is calculated. Then, the enzymes associated with the added EC numbers are extracted from the child's enzyme metabolism.
        Then the gene-duplicated enzymes are calculated. Finally, the gene-duplicated enzyme pairs of the conserved core metabolism enzymes are reported.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`addedMetabolism`.
        
        Returns
        -------
        Set[Tuple[Enzyme, Enzyme]]
            Pairs of enzymes from the child clade. Calculated using the added EC numbers found by :func:`addedMetabolism`.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """        
        # get added metabolism
        addedMetabolismEnzymes = self.addedMetabolismEnzymes(majorityPercentageCoreMetabolism).getEnzymes()
        
        # get gene-duplicated enzyme pairs
        geneDuplicated = self.childClade.geneDuplicatedEnzymePairs(majorityPercentageCoreMetabolism)
        
        # filter gene-duplicated enzyme pairs for the ones with both enzymes in the added metabolism
        geneDuplicatedAdded = set()
        
        for enzymeTuple in geneDuplicated:
            if enzymeTuple[0] in addedMetabolismEnzymes and enzymeTuple[1] in addedMetabolismEnzymes:
                geneDuplicatedAdded.add(enzymeTuple)
        
        return geneDuplicatedAdded
    
    
    def lostMetabolismGeneDuplicatedEnzymePairs(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism) -> Set[Tuple[Enzyme, Enzyme]]:
        """
        Pairs of gene-duplicated enzymes, derived from the lost core metabolism.
        
        First, the lost core metabolism is calculated. Then, the enzymes associated with the added EC numbers are extracted from the parent's enzyme metabolism.
        Then the gene-duplicated enzymes are calculated. Finally, the gene-duplicated enzyme pairs of the conserved core metabolism enzymes are reported.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`lostMetabolism`.
        
        Returns
        -------
        Set[Tuple[Enzyme, Enzyme]]
            Pairs of enzymes from the parent clade. Calculated using the lost EC numbers found by :func:`lostMetabolism`.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        # get added metabolism
        lostMetabolismEnzymes = self.lostMetabolismEnzymes(majorityPercentageCoreMetabolism).getEnzymes()
        
        # get gene-duplicated enzyme pairs
        geneDuplicated = self.childClade.geneDuplicatedEnzymePairs(majorityPercentageCoreMetabolism)
        
        # filter gene-duplicated enzyme pairs for the ones with both enzymes in the lost metabolism
        geneDuplicatedLost = set()
        
        for enzymeTuple in geneDuplicated:
            if enzymeTuple[0] in lostMetabolismEnzymes and enzymeTuple[1] in lostMetabolismEnzymes:
                geneDuplicatedLost.add(enzymeTuple)
        
        return geneDuplicatedLost
    
    
    def divergedMetabolismGeneDuplicatedEnzymePairs(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism) -> Set[Tuple[Enzyme, Enzyme]]:
        """
        Pairs of gene-duplicated enzymes, derived from the diverged core metabolism.
        
        First, the diverged core metabolism is calculated. Then, the enzymes associated with the added EC numbers are extracted from the collective parent's and child's metabolism individually.
        Then, for parent and child, the gene-duplicated enzyme pairs are calculated. Finally, the gene-duplicated enzymes where both enzymes are in the conserved core metabolism are reported.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`divergedMetabolism`.
        colour : bool, optional
            If *True*, colours the lost enzyme edges in blue, and the added enzyme edges in red. Gene-duplicated enzyme edges of the parent are coloured in green, the ones of the child in yellow.
            When doing so, a single :class:`SubstanceEnzymeGraph` is returned, not a :class:`Tuple`. The colouring is realised by adding a 'colour' attribute to each edge. Nodes are not coloured.
        
        Returns
        -------
        Set[Tuple[Enzyme, Enzyme]
            Sets of tuples of gene-duplicated enzyme pairs calculated using the diverged EC numbers found by :func:`divergedMetabolism`.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        # get diverged metabolism
        divergedMetabolismEnzymes = self.divergedMetabolismEnzymes(majorityPercentageCoreMetabolism).getEnzymes()
        
        # get gene-duplicate enzyme pairs
        parentGeneDuplicated = self.parentClade.geneDuplicatedEnzymePairs(majorityPercentageCoreMetabolism)
        childGeneDuplicated = self.childClade.geneDuplicatedEnzymePairs(majorityPercentageCoreMetabolism)
        
        # filter gene-duplicated enzyme pairs for the ones with both enzymes in the diverged metabolism
        parentGeneDuplicatedDiverged = set()
        childGeneDuplicatedDiverged = set()
        
        for enzymeTuple in parentGeneDuplicated:
            if enzymeTuple[0] in divergedMetabolismEnzymes and enzymeTuple[1] in divergedMetabolismEnzymes:
                parentGeneDuplicatedDiverged.add(enzymeTuple)
        
        for enzymeTuple in childGeneDuplicated:
            if enzymeTuple[0] in divergedMetabolismEnzymes and enzymeTuple[1] in divergedMetabolismEnzymes:
                childGeneDuplicatedDiverged.add(enzymeTuple)
        
        return parentGeneDuplicatedDiverged.union(childGeneDuplicatedDiverged)
    
    
    def unifiedMetabolismGeneDuplicatedEnzymePairs(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism) -> Set[Tuple[Enzyme, Enzyme]]:
        """
        Pairs of gene-duplicated enzymes, derived from the unified core metabolisms.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`conservedMetabolism`.
        
        Returns
        -------
        Set[Tuple[Enzyme, Enzyme]
            Set of enzyme pairs representing the gene-duplicated enzymes of the combined metabolic networks of both, child and parent.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """        
        parentGeneDuplicated = self.parentClade.geneDuplicatedEnzymePairs(majorityPercentageCoreMetabolism)
        childGeneDuplicated = self.childClade.geneDuplicatedEnzymePairs(majorityPercentageCoreMetabolism)
        
        return parentGeneDuplicated.union(childGeneDuplicated)
    
    
    
    
    
    # set-operations on neofunctionalised core metabolism
    ## for enzyme graphs
    def conservedMetabolismNeofunctionalisedEnzymes(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, colour = False):
        """
        Two Substance-Enzyme graphs of neofunctionalised enzymes, derived from the conserved core metabolism.
        
        First, the conserved core metabolism is calculated. Then, the enzymes associated with the conserved EC numbers are extracted from the collective parent's and child's metabolism individually.
        Then, for parent and child, the gene-duplicated enzymes are calculated. Then, the gene-duplicated enzymes of the conserved core metabolism enzymes are identified.
        Finally, the pairs of enzymes in which EC numbers differ are reported.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`conservedMetabolism`.
        colour : bool, optional
            If *True*, colours the enzyme edges from the parent in blue, and from the child in red. Neofunctionalised enzyme edges of the parent are coloured in green, the ones of the child in yellow.
            When doing so, a single :class:`SubstanceEnzymeGraph` is returned, not a :class:`Tuple`. The colouring is realised by adding a 'colour' attribute to each edge. Nodes are not coloured.
        
        Returns
        -------
        Tuple[SubstanceEnzymeGraph, SubstanceEnzymeGraph] or SubstanceEnzymeGraph
            Tuple of two Substance-Enzyme graphs calculated using the conserved EC numbers found by :func:`conservedMetabolism`. The first graph is from the parent clade, the second graph from the child clade.
            If `colour` == *True*, returns a single Substance-Enzyme graph.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        conservedMetabolismEnzymes = self.conservedMetabolismEnzymes(majorityPercentageCoreMetabolism, colour = colour)
        
        parentNeofunctionalised= self.parentClade.neofunctionalisedEnzymes(majorityPercentageCoreMetabolism, colour = False)
        childNeofunctionalised = self.childClade.neofunctionalisedEnzymes(majorityPercentageCoreMetabolism, colour = False)
        
        if colour is True:
            parentEdges = parentNeofunctionalised.getEdges()
            childEdges = childNeofunctionalised.getEdges()
            
            graph = conservedMetabolismEnzymes
            
            Export.addColourAttribute(graph, colour = Export.Colour.GREEN, nodes = False, edges = parentEdges)
            Export.addColourAttribute(graph, colour = Export.Colour.YELLOW, nodes = False, edges = childEdges)
            
            graph.name = 'Conserved metabolism neofunctionalised enzymes ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
            
            return graph
        else:
            parentGraph = conservedMetabolismEnzymes[0].removeAllEnzymesExcept(parentNeofunctionalised.getEnzymes())
            childGraph = conservedMetabolismEnzymes[1].removeAllEnzymesExcept(childNeofunctionalised.getEnzymes())
            
            parentGraph.name = 'Conserved metabolism neofunctionalised enzymes *' + ' '.join(self.parentNCBInames) + '* -> ' + ' '.join(self.childNCBInames)
            childGraph.name = 'Conserved metabolism neofunctionalised enzymes ' + ' '.join(self.parentNCBInames) + ' -> *' + ' '.join(self.childNCBInames) + '*'
        
            return (parentGraph, childGraph)
    
    
    def addedMetabolismNeofunctionalisedEnzymes(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism) -> SubstanceEnzymeGraph:
        """
        Substance-Enzyme graph of neofunctionalised enzymes, derived from the added core metabolism.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`addedMetabolism`.
        
        Returns
        -------
        SubstanceEnzymeGraph
            Substance-Enzyme graph of enzymes from the child clade. Calculated using the added EC numbers found by :func:`addedMetabolism`.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        parentCoreMetabolism = self.parentClade.coreMetabolism(majorityPercentageCoreMetabolism)
        childCoreMetabolism = self.childClade.coreMetabolism(majorityPercentageCoreMetabolism)
        addedECs = GeneFunctionAddition.getECs(parentCoreMetabolism, childCoreMetabolism)
        
        childGraph = self.childClade.neofunctionalisedEnzymes(majorityPercentageCoreMetabolism, colour = False).keepEnzymesByEC(addedECs)
        childGraph.name = 'Added metabolism neofunctionalised enzymes ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
        
        return childGraph
    
    
    def lostMetabolismNeofunctionalisedEnzymes(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism) -> SubstanceEnzymeGraph:
        """
        Substance-Enzyme graph of neofunctionalised enzymes, derived from the lost core metabolism.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`lostMetabolism`.
        
        Returns
        -------
        SubstanceEnzymeGraph
            Substance-Enzyme graph of enzymes from the parent clade. Calculated using the lost EC numbers found by :func:`lostMetabolism`.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        parentCoreMetabolism = self.parentClade.coreMetabolism(majorityPercentageCoreMetabolism)
        childCoreMetabolism = self.childClade.coreMetabolism(majorityPercentageCoreMetabolism)
        lostECs = GeneFunctionLoss.getECs(parentCoreMetabolism, childCoreMetabolism)
        
        parentGraph = self.parentClade.neofunctionalisedEnzymes(majorityPercentageCoreMetabolism, colour = False).keepEnzymesByEC(lostECs)
        parentGraph.name = 'Lost metabolism neofunctionalised enzymes ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
        
        return parentGraph
    
    
    def divergedMetabolismNeofunctionalisedEnzymes(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, colour = False):
        """
        Two Substance-Enzyme graphs of neofunctionalised enzymes, derived from the diverged core metabolism.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`divergedMetabolism`.
        colour : bool, optional
            If *True*, colours the lost enzyme edges in blue, and the added enzyme edges in red. Neofunctionalised enzyme edges of the parent are coloured in green, the ones of the child in yellow.
            When doing so, a single :class:`SubstanceEnzymeGraph` is returned, not a :class:`Tuple`. The colouring is realised by adding a 'colour' attribute to each edge. Nodes are not coloured.
        
        Returns
        -------
        Tuple[SubstanceEnzymeGraph, SubstanceEnzymeGraph] or SubstanceEnzymeGraph
            Tuple of two Substance-Enzyme graphs calculated using the diverged EC numbers found by :func:`divergedMetabolism`. The first graph is from the parent clade, the second graph from the child clade.
            If `colour` == *True*, returns a single Substance-Enzyme graph, coloured blue for parent and red for child.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        divergedMetabolismEnzymes = self.divergedMetabolismEnzymes(majorityPercentageCoreMetabolism, colour = colour)
        
        parentNeofunctionalised = self.parentClade.neofunctionalisedEnzymes(majorityPercentageCoreMetabolism, colour = False)
        childNeofunctionalised = self.childClade.neofunctionalisedEnzymes(majorityPercentageCoreMetabolism, colour = False)
        
        if colour is True:
            parentEdges = parentNeofunctionalised.getEdges()
            childEdges = childNeofunctionalised.getEdges()
            
            graph = divergedMetabolismEnzymes
            
            Export.addColourAttribute(graph, colour = Export.Colour.GREEN, nodes = False, edges = parentEdges)
            Export.addColourAttribute(graph, colour = Export.Colour.YELLOW, nodes = False, edges = childEdges)
            
            graph.name = 'Diverged metabolism neofunctionalised enzymes ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
            
            return graph
        else:
            parentGraph = divergedMetabolismEnzymes[0].removeAllEnzymesExcept(parentNeofunctionalised.getEnzymes())
            childGraph = divergedMetabolismEnzymes[1].removeAllEnzymesExcept(childNeofunctionalised.getEnzymes())
            
            parentGraph.name = 'Diverged metabolism neofunctionalised enzymes *' + ' '.join(self.parentNCBInames) + '* -> ' + ' '.join(self.childNCBInames)
            childGraph.name = 'Diverged metabolism neofunctionalised enzymes ' + ' '.join(self.parentNCBInames) + ' -> *' + ' '.join(self.childNCBInames) + '*'
        
            return (parentGraph, childGraph)
    
    
    def unifiedMetabolismNeofunctionalisedEnzymes(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, colour = False) -> SubstanceEnzymeGraph:
        """
        Substance-Enzyme graph of neofunctionalised enzymes, derived from the unified core metabolisms.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`conservedMetabolism`.
        colour : bool, optional
            If *True*, colours the parent's enzyme edges in blue, and the child's enzyme edges in red. Neofunctionalised enzyme edges of the parent are coloured in green, the ones of the child in yellow.
            The colouring is realised by adding a 'colour' attribute to each edge. Nodes are not coloured.
        
        Returns
        -------
        SubstanceEcGraph
            The substance-Enzyme graph representing the combined metabolic networks of both, child and parent. If `colour` == *True*, coloured differently for the lost, conserved, and added edges. Nodes are not coloured.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """        
        parentNeofunctionalised = self.parentClade.neofunctionalisedEnzymes(majorityPercentageCoreMetabolism, colour = False)
        childNeofunctionalised = self.childClade.neofunctionalisedEnzymes(majorityPercentageCoreMetabolism, colour = False)
        
        if colour is False:
            graph = parentNeofunctionalised.union(childNeofunctionalised, addCount = False, updateName = False)
        
        else:
            unifiedMetabolismEnzymes = self.unifiedMetabolismEnzymes(majorityPercentageCoreMetabolism, colour = True)
            
            parentEdges = parentNeofunctionalised.getEdges()
            childEdges = childNeofunctionalised.getEdges()
            
            graph = unifiedMetabolismEnzymes
            
            Export.addColourAttribute(graph, colour = Export.Colour.GREEN, nodes = False, edges = parentEdges)
            Export.addColourAttribute(graph, colour = Export.Colour.YELLOW, nodes = False, edges = childEdges)
            
            graph.name = 'Diverged metabolism neofunctionalised enzymes ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
                
        return graph
    
    
    
    
    
    ## for EC graphs
    def conservedMetabolismNeofunctionalisedECs(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation = defaultMajorityPercentageNeofunctionalisation, colour = False):
        """
        Two Substance-EC graphs of "neofunctionalised" EC numbers, derived from the conserved core metabolism.
        
        First, the conserved core metabolism is calculated. Then, the enzymes associated with the conserved EC numbers are extracted from the collective parent's and child's metabolism individually.
        Then, for parent and child, the gene-duplicated enzymes are calculated. Then, the gene-duplicated enzymes of the conserved core metabolism enzymes are identified.
        Then, the pairs of enzymes in which EC numbers differ are identified. Finally, the EC numbers which are part of these function changes are reported.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`conservedMetabolism`.
        colour : bool, optional
            If *True*, colours the EC edges from the parent in blue, and from the child in red. "Neofunctionalised" EC edges of the parent are coloured in green, the ones of the child in yellow.
            When doing so, a single :class:`SubstanceEcGraph` is returned, not a :class:`Tuple`. The colouring is realised by adding a 'colour' attribute to each edge. Nodes are not coloured.
        
        Returns
        -------
        Tuple[SubstanceEcGraph, SubstanceEcGraph] or SubstanceEcGraph
            Tuple of two Substance-EC graphs calculated using the conserved EC numbers found by :func:`conservedMetabolism`. The first graph is from the parent clade, the second graph from the child clade.
            If `colour` == *True*, returns a single Substance-EC graph.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        conservedMetabolism = self.conservedMetabolism(majorityPercentageCoreMetabolism)
        
        parentNeofunctionalised= self.parentClade.neofunctionalisedECs(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation, colour = False)
        childNeofunctionalised = self.childClade.neofunctionalisedECs(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation, colour = False)
        
        if colour is True:
            parentEdges = parentNeofunctionalised.getEdges()
            childEdges = childNeofunctionalised.getEdges()
            
            graph = conservedMetabolism
            
            Export.addColourAttribute(graph, colour = Export.Colour.GREEN, nodes = False, edges = parentEdges)
            Export.addColourAttribute(graph, colour = Export.Colour.YELLOW, nodes = False, edges = childEdges)
            
            graph.name = 'Conserved metabolism neofunctionalised ECs ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
            
            return graph
        else:
            parentGraph = conservedMetabolism[0].removeAllECsExcept(parentNeofunctionalised.getECs())
            childGraph = conservedMetabolism[1].removeAllECsExcept(childNeofunctionalised.getECs())
            
            parentGraph.name = 'Conserved metabolism neofunctionalised ECs *' + ' '.join(self.parentNCBInames) + '* -> ' + ' '.join(self.childNCBInames)
            childGraph.name = 'Conserved metabolism neofunctionalised ECs ' + ' '.join(self.parentNCBInames) + ' -> *' + ' '.join(self.childNCBInames) + '*'
        
            return (parentGraph, childGraph)
    
    
    def addedMetabolismNeofunctionalisedECs(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation = defaultMajorityPercentageNeofunctionalisation) -> SubstanceEcGraph:
        """
        Substance-EC graph of "neofunctionalised" EC numbers, derived from the added core metabolism.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`addedMetabolism`.
        
        Returns
        -------
        SubstanceEcGraph
            Substance-EC graph of ECs from the child clade. Calculated using the added EC numbers found by :func:`addedMetabolism`.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        parentCoreMetabolism = self.parentClade.coreMetabolism(majorityPercentageCoreMetabolism)
        childCoreMetabolism = self.childClade.coreMetabolism(majorityPercentageCoreMetabolism)
        addedECs = GeneFunctionAddition.getECs(parentCoreMetabolism, childCoreMetabolism)
        
        childGraph = self.childClade.neofunctionalisedECs(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation, colour = False).removeAllECsExcept(addedECs)
        childGraph.name = 'Added metabolism neofunctionalised ECs ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
        
        return childGraph
    
    
    def lostMetabolismNeofunctionalisedECs(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation = defaultMajorityPercentageNeofunctionalisation) -> SubstanceEcGraph:
        """
        Substance-EC graph of "neofunctionalised" EC numbers, derived from the lost core metabolism.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`lostMetabolism`.
        
        Returns
        -------
        SubstanceEcGraph
            Substance-EC graph of ECs from the parent clade. Calculated using the lost EC numbers found by :func:`lostMetabolism`.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        parentCoreMetabolism = self.parentClade.coreMetabolism(majorityPercentageCoreMetabolism)
        childCoreMetabolism = self.childClade.coreMetabolism(majorityPercentageCoreMetabolism)
        lostECs = GeneFunctionLoss.getECs(parentCoreMetabolism, childCoreMetabolism)
        
        parentGraph = self.parentClade.neofunctionalisedECs(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation, colour = False).removeAllECsExcept(lostECs)
        parentGraph.name = 'Lost metabolism neofunctionalised ECs ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
        
        return parentGraph
    
    
    def divergedMetabolismNeofunctionalisedECs(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation = defaultMajorityPercentageNeofunctionalisation, colour = False):
        """
        Two Substance-EC graphs of "neofunctionalised" EC numbers, derived from the diverged core metabolism.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`divergedMetabolism`.
        colour : bool, optional
            If *True*, colours the lost EC edges in blue, and the added EC edges in red. "Neofunctionalised" EC edges of the parent are coloured in green, the ones of the child in yellow.
            When doing so, a single :class:`SubstanceEcGraph` is returned, not a :class:`Tuple`. The colouring is realised by adding a 'colour' attribute to each edge. Nodes are not coloured.
        
        Returns
        -------
        Tuple[SubstanceEcGraph, SubstanceEcGraph] or SubstanceEcGraph
            Tuple of two Substance-EC graphs calculated using the diverged EC numbers found by :func:`divergedMetabolism`. The first graph is from the parent clade, the second graph from the child clade.
            If `colour` == *True*, returns a single Substance-EC graph, coloured blue for parent and red for child.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """
        divergedMetabolism = self.divergedMetabolism(majorityPercentageCoreMetabolism, colour = colour)
        
        parentNeofunctionalised = self.parentClade.neofunctionalisedECs(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation, colour = False)
        childNeofunctionalised = self.childClade.neofunctionalisedECs(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation, colour = False)
        
        if colour is True:
            parentEdges = parentNeofunctionalised.getEdges()
            childEdges = childNeofunctionalised.getEdges()
            
            graph = divergedMetabolism
            
            Export.addColourAttribute(graph, colour = Export.Colour.GREEN, nodes = False, edges = parentEdges)
            Export.addColourAttribute(graph, colour = Export.Colour.YELLOW, nodes = False, edges = childEdges)
            
            graph.name = 'Diverged metabolism neofunctionalised ECs ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
            
            return graph
        else:
            parentGraph = divergedMetabolism[0].removeAllECsExcept(parentNeofunctionalised.getECs())
            childGraph = divergedMetabolism[1].removeAllECsExcept(childNeofunctionalised.getECs())
            
            parentGraph.name = 'Diverged metabolism neofunctionalised ECs *' + ' '.join(self.parentNCBInames) + '* -> ' + ' '.join(self.childNCBInames)
            childGraph.name = 'Diverged metabolism neofunctionalised ECs ' + ' '.join(self.parentNCBInames) + ' -> *' + ' '.join(self.childNCBInames) + '*'
        
            return (parentGraph, childGraph)
    
    
    def unifiedMetabolismNeofunctionalisedECs(self, majorityPercentageCoreMetabolism = defaultMajorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation = defaultMajorityPercentageNeofunctionalisation, colour = False) -> SubstanceEcGraph:
        """
        Substance-EC graph of "neofunctionalised" EC numbers, derived from the unified core metabolisms.
        
        Parameters
        ----------
        majorityPercentageCoreMetabolism : int, optional
            See :func:`conservedMetabolism`.
        colour : bool, optional
            If *True*, colours the parent's EC edges in blue, and the child's EC edges in red. "Neofunctionalised" EC edges of the parent are coloured in green, the ones of the child in yellow.
            The colouring is realised by adding a 'colour' attribute to each edge. Nodes are not coloured.
        
        Returns
        -------
        SubstanceEcGraph
            The substance-EC graph representing the combined metabolic networks of both, child and parent. If `colour` == *True*, coloured differently for the lost, conserved, and added edges. Nodes are not coloured.
        
        Raises
        ------
        TypeError
            If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`FEV_KEGG.KEGG.Organism.Group._getGraphsParallelly`.
        HTTPError
            If fetching any of the underlying graphs fails.
        URLError
            If connection to KEGG fails.
        """        
        parentNeofunctionalised = self.parentClade.neofunctionalisedECs(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation, colour = False)
        childNeofunctionalised = self.childClade.neofunctionalisedECs(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation, colour = False)
        
        if colour is False:
            graph = parentNeofunctionalised.union(childNeofunctionalised, addCount = False, updateName = False)
        
        else:
            unifiedMetabolism = self.unifiedMetabolism(majorityPercentageCoreMetabolism, colour = True)
            
            parentEdges = parentNeofunctionalised.getEdges()
            childEdges = childNeofunctionalised.getEdges()
            
            graph = unifiedMetabolism
            
            Export.addColourAttribute(graph, colour = Export.Colour.GREEN, nodes = False, edges = parentEdges)
            Export.addColourAttribute(graph, colour = Export.Colour.YELLOW, nodes = False, edges = childEdges)
            
            graph.name = 'Diverged metabolism neofunctionalised ECs ' + ' '.join(self.parentNCBInames) + ' -> ' + ' '.join(self.childNCBInames)
                
        return graph
    
    
    
    
    
    



    
    


class NestedCladePair(CladePair):
    
    def __init__(self, parent, child, excludeUnclassified = defaultExcludeUnclassified):
        """
        Two clades in NCBI taxonomy, 'child' is assumed younger and must be nested somewhere inside 'parent'.
        
        This only checks nestedness for the first node found in taxonomy, by the first parent's/child's NCBI name, respectively. The latter being relevant if you pass a :class:`Clade`, which has a list of NCBI names, or a list of NCBI names itself.
        
        Parameters
        ----------
        parent : str or List[str] or Clade
            Path(s) of the parent clade's taxon, as defined by NCBI taxonomy, e.g. 'Proteobacteria/Gammaproteobacteria'. Or a ready :class:`Clade` object.
        child : str or List[str] or Clade
            Path(s) of the child clade's taxon, as defined by NCBI taxonomy, e.g. 'Enterobacter'. Or a ready :class:`Clade` object.
        excludeUnclassified : bool, optional
            If *True*, ignore taxons with a path containing the string 'unclassified'.
        
        Attributes
        ----------
        self.childClade : :class:`Clade`
        self.parentClade : :class:`Clade`
        
        Raises
        ------
        ValueError
            If parent or child are unknown taxons. Or if the child taxon is not actually a child of the parent taxon.
        """
        # read first NCBI name from Clade object, if necessary
        if isinstance(parent, Clade):
            parentNCBIname = parent.ncbiNames[0]
        elif not isinstance(parent, str):
            # must be iterable, else fail
            parentNCBIname = parent[0]
        
        if isinstance(child, Clade):
            childNCBIname = child.ncbiNames[0]
        elif not isinstance(child, str):
            # must be iterable, else fail
            childNCBIname = child[0]
            
        # check if child is really a child of parent
        taxonomy = NCBI.getTaxonomy()
        parentNode = taxonomy.searchNodesByPath(parentNCBIname, exceptPaths=('unclassified' if excludeUnclassified else None))
        if parentNode is None or len(parentNode) == 0:
            raise ValueError("No clade of this path found: " + parentNCBIname)
        else: # only consider first element
            parentNode = parentNode[0]
        
        childNode = taxonomy.searchNodesByPath(childNCBIname, exceptPaths=('unclassified' if excludeUnclassified else None))
        if childNode is None or len(childNode) == 0:
            raise ValueError("No clade of this path found: " + childNCBIname)
        else: # only consider first element
            childNode = childNode[0]
        
        foundParent = False
        for ancestor in childNode.ancestors:
            if Taxonomy.nodePath2String(ancestor) == Taxonomy.nodePath2String(parentNode):
                foundParent = True
                break
        
        if foundParent == False:
            raise ValueError("Child is not a descendant of parent.")
        
        super().__init__(parent, child, excludeUnclassified)
