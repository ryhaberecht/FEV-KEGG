from FEV_KEGG.Graph.Elements import EcNumber
from FEV_KEGG.Graph.SubstanceGraphs import SubstanceEcGraph
from FEV_KEGG.KEGG.File import cache
from FEV_KEGG.KEGG.NUKA import NUKA
from FEV_KEGG.settings import verbosity as init_verbosity
from enum import Enum
import FEV_KEGG.Evolution.Clade
from typing import Set


class GoldmanLUCA(object):

    def __init__(self):
        """
        Last Universal Common Ancestor by Goldman et al.
        
        This is the Last Universal Common Ancestor, as described in [1]_.
        The original work on LUCA, however, does not specify enzyme function, but merely COGs [3]_.
        This class already contains the list of LUCA's enzymes from the above paper, as depicted in the first table of said paper [2]_.
        As the most plausible minimal set of enzymatic functions, the authors chose the intersection of EC numbers found in universal sequence + structure, combined with the ones found in universal sequence + structure + function. See the original source for details.
        
        This list is parsed and converted into a SubstanceEcGraph.
        Conversion is done by using the graph of a hypothetical 'complete' organism - NUKA - which possesses all EC numbers known to all metabolic KEGG pathways, see :mod:`FEV_KEGG.KEGG.NUKA`
        All EC numbers not present in LUCA are filtered out.
        Keep in mind, though, that LUCA's EC numbers only contain three levels, to more adequately model the likely patchwork evolution in ancient times.
        Therefore, all EC numbers starting with the sub-class remain, regardless of substrate specificity.
        
        Conversion into another type of graph is not supported, because LUCA is a strictly hypothetical organism without any exactly known genes.
        
        Attributes
        ----------
        self.nameAbbreviation : str
        
        References
        ----------
        .. [1] Goldman et al. (2012), "The Enzymatic and Metabolic Capabilities of Early Life", `<https://doi.org/10.1371/journal.pone.0039912>`_
        .. [2] Goldman et al. (2012), Table 1, `<https://doi.org/10.1371/journal.pone.0039912.t001>`_
        .. [3] Mirkin et al. (2003), "Algorithms for computing parsimonious evolutionary scenarios for genome evolution, the last universal common ancestor and dominance of horizontal gene transfer in the evolution of prokaryotes", `<https://doi.org/10.1186/1471-2148-3-2>`_
        """
        self.nameAbbreviation = 'GoldmanLUCA'        
    
    @property
    def ecNumbers(self) -> Set[EcNumber]:
        """
        GoldmanLUCA's EC numbers.
        
        Generalised to the first three levels. The last level is always a wildcard.
        
        Returns
        -------
        Set[EcNumber]
            Set of the EC numbers predicted by Goldman et al. to belong to LUCA.
        """
        lucaEcNumberStrings = [
            '1.3.1.-',
            '2.4.1.-',
            '2.7.1.-',
            '2.7.7.-',
            '3.1.2.-',
            '3.1.4.-',
            '3.2.1.-',
            '3.5.1.-',
            '4.1.2.-',
            '6.3.2.-'
            ]
        
        lucaEcNumbers = set()
        for string in lucaEcNumberStrings:
            lucaEcNumbers.add(EcNumber(string))
        
        return lucaEcNumbers
    
    @property
    @cache(folder_path = 'GoldmanLUCA/graph', file_name = 'SubstanceEcGraph')
    def substanceEcGraph(self) -> SubstanceEcGraph:
        """
        GoldmanLUCA's substance-EC graph.
        
        Returns
        -------
        SubstanceEcGraph
            Contains all substrates/products and all EC numbers in :mod:`FEV_KEGG.KEGG.NUKA` filtered by the EC numbers predicted by Goldman et al. for LUCA.
        
        Raises
        ------
        HTTPError
            If any underlying organism, pathway, or gene of NUKA does not exist.
        URLError
            If connection to KEGG fails.
        """
        lucaEcNumbers = self.ecNumbers
        
        # copy SubstanceEcGraph from NUKA
        graph = NUKA().substanceEcGraph.copy()
        
        # delete all edges with EC numbers that are not contained in the EC numbers of LUCA
        for edge in list(graph.getEdges()): # use a copy, to avoid modification in-place during iteration
            substrate, product, nukaEC = edge
            
            nukaEcContainedInLuca = False
            
            for lucaEC in lucaEcNumbers:
                if lucaEC.contains(nukaEC):
                    nukaEcContainedInLuca = True
                    break
            
            if nukaEcContainedInLuca == False: # NUKA's EC number is not contained in any of LUCA'S EC numbers
                graph.removeEcEdge(substrate, product, nukaEC, bothDirections = False) # remove this edge, both directions will be removed eventually
        
        graph.name = 'Substance-EC GoldmanLUCA'
        
        # remove isolated nodes
        graph.removeIsolatedNodes()
        
        if init_verbosity > 0:
            print('calculated ' + graph.name)
        
        return graph
    
    
    


class CoreLUCA(object):

    class CladeType(Enum):
        """
        Possible types of CoreLUCA.
        
        Each accordings to a single, or a combination of, top-clades of NCBI.
        Only the 'universal' clade gives you the "true" LUCA.
        """
        universal = "/"
        archaea = "/Archaea"
        bacteria = "/Bacteria"
        eukaryota = "/Eukaryota"
        archaeaBacteria = [archaea, bacteria]
        archaeaEukaryota = [archaea, eukaryota]
        bacteriaEukaryota = [bacteria, eukaryota]
    
    def __init__(self, clade: 'CoreLUCA.CladeType'):
        """
        Last Universal Common Ancestor by intersection of many or all organisms in KEGG.
        
        This is the Last Universal Common Ancestor, as defined by a common "core metabolism" shared among all organisms known to KEGG within a certain NCBI top-clade.
        This would include Bacteria, Arachaea, and Eukaryota; which is a very big data set!
        Alternatively, you can specify which isolated top-clade to use, using `clade`, e.g. yielding the Bacteria-LUCA, or Archaea-LUCA.
        For each species only the first organism is considered, to prevent statistical overrepresentation.
        
        Conversion into another type of graph is not supported, because LUCA is a strictly hypothetical organism without any exactly known genes.
        
        Parameters
        ----------
        clade : CoreLUCA.CladeType
            Which clade to use for defining a LUCA. Using 'archae' obviously only gives an Archae-LUCA, not the "true" LUCA, etc.
        
        Attributes
        ----------
        self.nameAbbreviation : str
        self.clade : :class:`FEV_KEGG.Evolution.Clade.Clade`
        self.cladeType : :class:`CoreLUCA.CladeType`
        
        Raises
        ------
        HTTPError
            If any underlying organism, pathway, or gene does not exist.
        URLError
            If connection to KEGG fails.
        
        Warnings
        --------
        This function takes hours to days to complete, and requires several gigabytes of memory, disk space, and network traffic!
        """
        if not isinstance(clade, self.CladeType):
            raise ValueError("No valid top-clade of type CoreLUCA.CladeType specified")
        self.cladeType = clade 
        
        if clade == self.CladeType.universal:
            self.nameAbbreviation = 'Universal-CoreLUCA'
            self.clade = self._getUniversalClade()
        
        elif clade == self.CladeType.archaea:
            self.nameAbbreviation = 'Archaea-CoreLUCA'
            self.clade = self._getArchaeaClade()
        
        elif clade == self.CladeType.bacteria:
            self.nameAbbreviation = 'Bacteria-CoreLUCA'
            self.clade = self._getBacteriaClade()
            
        elif clade == self.CladeType.eukaryota:
            self.nameAbbreviation = 'Eukaryota-CoreLUCA'
            self.clade = self._getEukaryotaClade()
            
        elif clade == self.CladeType.archaeaBacteria:
            self.nameAbbreviation = 'Archaea-Bacteria-CoreLUCA'
            self.clade = self._getArchaeaBacteriaClade()
            
        elif clade == self.CladeType.archaeaEukaryota:
            self.nameAbbreviation = 'Archaea-Eukaryota-CoreLUCA'
            self.clade = self._getArchaeaEukaryotaClade()
            
        elif clade == self.CladeType.bacteriaEukaryota:
            self.nameAbbreviation = 'Bacteria-Eukaryota-CoreLUCA'
            self.clade = self._getBacteriaEukaryotaClade()
            
        else:
            raise ValueError("Unknown CoreLUCA.CladeType")
    
    @cache('CoreLUCA/graph', 'universal_clade')
    def _getUniversalClade(self):
        return self._getClade(self.CladeType.universal.value)
    
    @cache('CoreLUCA/graph', 'archaea_clade')
    def _getArchaeaClade(self):
        return self._getClade(self.CladeType.archaea.value)
    
    @cache('CoreLUCA/graph', 'bacteria_clade')
    def _getBacteriaClade(self):
        return self._getClade(self.CladeType.bacteria.value)
    
    @cache('CoreLUCA/graph', 'eukaryota_clade')
    def _getEukaryotaClade(self):
        return self._getClade(self.CladeType.eukaryota.value)
    
    @cache('CoreLUCA/graph', 'archaea_bacteria_clade')
    def _getArchaeaBacteriaClade(self):
        return self._getClade(self.CladeType.archaeaBacteria.value)
    
    @cache('CoreLUCA/graph', 'archaea_eukaryota_clade')
    def _getArchaeaEukaryotaClade(self):
        return self._getClade(self.CladeType.archaeaEukaryota.value)
    
    @cache('CoreLUCA/graph', 'bacteria_eukaryota_clade')
    def _getBacteriaEukaryotaClade(self):
        return self._getClade(self.CladeType.bacteriaEukaryota.value)
    
    def _getClade(self, ncbiName):
        clade = FEV_KEGG.Evolution.Clade.Clade(ncbiName, excludeUnclassified=True, oneOrganismPerSpecies=True)
        # pre-populate object memory with the full graph
        clade.collectiveMetabolism(excludeMultifunctionalEnzymes=False)
        return clade
    
    
    def coreMetabolism(self, majorityPercentage) -> SubstanceEcGraph:
        """
        CoreLUCA's core metabolism.
        
        Parameters
        ----------
        majorityPercentage : float
            Percentage for determining how many organisms have to possess an EC edge, for it to be included in this 'core metabolism'.
        
        Returns
        -------
        SubstanceEcGraph
            Contains all substrates/products and all EC numbers in the "core metabolism" of the top-clade you chose.
        """
        graph = self.clade.coreMetabolism(majorityPercentageCoreMetabolism=majorityPercentage, excludeMultifunctionalEnzymes=False)
        graph.removePartialEcNumbers()
        graph.removeIsolatedNodes()
        graph.name = 'Substance-EC ' + self.nameAbbreviation
        return graph
    
    def collectiveMetabolism(self) -> SubstanceEcGraph:
        """
        CoreLUCA's collective metabolism, i.e. core metabolism with the lowest possible majorityPercentage value.
        
        Returns
        -------
        SubstanceEcGraph
            Contains all substrates/products and all EC numbers of any organism in the top-clade you chose.
        """
        return self.clade.collectiveMetabolism(excludeMultifunctionalEnzymes=False)
        