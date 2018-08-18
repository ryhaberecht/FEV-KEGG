from builtins import str
import re

from FEV_KEGG.lib.Biopython.KEGG.KGML import KGML_pathway
from FEV_KEGG.KEGG.DataTypes import Gene
from typing import Set, List, Iterable, Generator

from FEV_KEGG.Graph.SubstanceGraphs import SubstanceEnzymeGraph, SubstanceEcGraph, SubstanceGeneGraph, SubstanceReactionGraph, Conversion
from FEV_KEGG.KEGG import Database
from FEV_KEGG.KEGG.Database import NoKnownPathwaysError
from FEV_KEGG.KEGG.File import cache, cacheEntry
import tqdm
import FEV_KEGG.settings as settings

import concurrent.futures
from FEV_KEGG.Util import Parallelism
import gc
import FEV_KEGG.quirks as quirks
from math import ceil



class Organism(object):
	GLOBAL_PATHWAY_PATTERN = re.compile('01[12][0-9]{2}')
	"""
	Pattern defining a global/overview pathway.
	
	These pathways should contain nothing more and nothing less than what is included in all non-global pathways. At least concerning metabolic pathways.
	Alas, they do not. Global/Overview pathways are often inconsistent with the union of all pathways. Also, they discard information about edge direction, containing only undirected edges.
	"""
	
	_DIGITS_PATTERN = re.compile('\d+')
	
	def __init__(self, nameAbbreviation: 'eco', skipExistsCheck = False):
		"""
		An Organism as listed in KEGG, e.g. Escherichia coli K-12 MG1655 "eco".
		
		Checks whether the organism actually exists before creating the object.
		
		Parameters
		----------
		nameAbbreviation : str
			The abbreviation of an organism as used in KEGG.
		skipExistsCheck : bool, optional
			If *True*, skips the check for existence of an organism identified by `nameAbbreviation`. Any subsequent method access may raise an error if the organism does not exist in KEGG!
		
		Attributes
		----------
		self.nameAbbreviation : str
		
		Raises
		------
		ValueError
			If the organism does not exist at all in KEGG.
		URLError
			If the connection to the KEGG server fails and the requested organism has not already been cached.
		
		Note
		----
		All operations are cached via :class:`FEV_KEGG.KEGG.Database`. This includes any downloads and the calculations available directly via this class.
		"""
		self.nameAbbreviation = nameAbbreviation
		
		if skipExistsCheck is False and Database.doesOrganismExist(nameAbbreviation) is False:
			raise ValueError('Organism abbreviation does not exist: ' + nameAbbreviation)
	
	@classmethod
	def __initBulk__(cls, nameAbbreviations: List[str]) -> List['Organism']:
		"""
		Creates many :class:`Organism` objects at once.
		
		Checking for existence is much faster when done this way, because it is parallelised.
		Pathologically non-existent organisms are automatically filtered, see :attr:`FEV_KEGG.quirks.NON_EXISTING_ORGANISMS`
		
		Parameters
		----------
		nameAbbreviations : list[str]
			List of abbreviation strings which will be used to create an :class:`Organism`.
		
		Returns
		-------
		List[Organism]
			List of :class:`Organism` objects. Or *None* if none of the `nameAbbreviations` existed.
			
		Raises
		------
		ValueError
			If any organism does not exist at all in KEGG.
		URLError
			If the connection to the KEGG server fails and the requested organism has not already been cached.
		"""
		existingOrganisms = Database.doesOrganismExistBulk(nameAbbreviations)
		if len(existingOrganisms) == 0:
			return None
		else:
			
			if len(nameAbbreviations) != len(existingOrganisms):
				# some organism do not exist
				nonExistingOrganisms = set(nameAbbreviations).difference_update(set(existingOrganisms))
				raise ValueError('Organism abbreviations do not exist: ' + ', '.join(nonExistingOrganisms))
			
			organisms = []
			
			# filter quirky organisms, without any pathway
			for nonExisting in quirks.NON_EXISTING_ORGANISMS:
				while True:
					try:
						existingOrganisms.remove(nonExisting)
					except ValueError:
						break
			
			for abbreviation in existingOrganisms:
				organism = cls(abbreviation, skipExistsCheck = True)
				organisms.append(organism)
			return organisms
	
	def __str__(self):
		return 'Organism(' + self.nameAbbreviation + ')'
	
	def __eq__(self, other):
		if isinstance(self, other.__class__):
			return self.nameAbbreviation == other.nameAbbreviation
		return False
	
	def __ne__(self, other):
		return not self == other
	
	def __hash__(self):
		return self.nameAbbreviation.__hash__()
	
	def __lt__(self, other):
		return self.nameAbbreviation.lower() < other.nameAbbreviation.lower()
	
	def __gt__(self, other):
		return self.nameAbbreviation.lower() > other.nameAbbreviation.lower()
	
	def __le__(self, other):
		return self.nameAbbreviation.lower() <= other.nameAbbreviation.lower()
	
	def __ge__(self, other):
		return self.nameAbbreviation.lower() >= other.nameAbbreviation.lower()
	
	def getPathway(self, pathwayName: '00260') -> KGML_pathway.Pathway:
		"""
		Gets a certain pathway for this organism. 
		
		Parameters
		----------
		pathwayName : str
			The name/number of the pathway, e.g. "00260". This will be expanded with `self.nameAbbreviation` to, e.g. "eco:00260".
		
		Returns
		-------
		Pathway
			The pathway object, or *None* if such a pathway does not exist.
		
		Raises
		------
		HTTPError
			If pathway does not exist.
		URLError
			If connection to KEGG fails.
		"""
		return Database.getPathway(self.nameAbbreviation, pathwayName)
	
	
	def getPathways(self, includeOverviewMaps = False) -> Set[KGML_pathway.Pathway]:
		"""
		Gets a set of all pathway objects for this organism.
		
		Parameters
		----------
		includeOverviewMaps : bool, optional
			Whether to include global/overview maps.
		
		Returns
		-------
		Set[Pathway]
			Set of pathways objects for all known pathways.
		
		Raises
		------
		NoKnownPathwaysError
			If the organism has no known pathways.
		HTTPError
			If pathway does not exist.
		URLError
			If connection to KEGG fails.
		"""
		return self.getPathwaysFromNames(self.getPathwayNames(self.getPathwayIDs(self.getPathwayDescriptions(includeOverviewMaps))))
	
	
	def getMetabolicPathways(self, includeOverviewMaps = False) -> Set[KGML_pathway.Pathway]:
		"""
		Gets a set of all metabolic pathway objects for this organism.
		
		Parameters
		----------
		includeOverviewMaps : bool, optional
			Whether to include global/overview maps.
		
		Returns
		-------
		Set[Pathway]
			Set of pathways objects for all known metabolic pathways.
		
		Raises
		------
		HTTPError
			If pathway does not exist.
		URLError
			If connection to KEGG fails.
		"""
		return self.getPathwaysFromNames(self.getPathwayNames(self.getPathwayIDs(self.getMetabolicPathwayDescriptions(includeOverviewMaps))))
	
	
	@property
	def metabolicPathways(self) -> Set[KGML_pathway.Pathway]:
		"""
		Pathways of metabolism, without overview or global maps.
		
		Returns
		-------
		Set[Pathway]
			Set of pathways objects for all known metabolic pathways, excluding global/overview pathways.
		
		Raises
		------
		HTTPError
			If pathway does not exist.
		URLError
			If connection to KEGG fails.
		"""
		return self.getMetabolicPathways(includeOverviewMaps = False)
	
	
	def getPathwaysFromNames(self, pathwayNameSet: Set[str]) -> Set[KGML_pathway.Pathway]:
		"""
		Gets a set of pathway objects for this organism, based on a set of pathway names, eg. {'00260', '01100'}.
		
		Parameters
		----------
		pathwayNameSet : Set[str]
			The names/numbers of the pathways, e.g. '00260'. This will be expanded with `self.nameAbbreviation` to, e.g. 'eco00260'.
		
		Returns
		-------
		Set[Pathway]
			Set of pathways objects for all pathways from `pathwayNameSet`. If pathway exists, but has no KGML format, the entry for this pathway is *None*.
			
		Raises
		------
		HTTPError
			If pathway does not exist.
		URLError
			If connection to KEGG fails.
		"""
		return Database.getPathwayBulk(self.nameAbbreviation, pathwayNameSet).values()
	
	
	def getPathwayDescriptions(self, includeOverviewMaps = False) -> Set[str]:
		"""
		Get pathway descriptions for this organism.
		
		Parameters
		----------
		includeOverviewMaps : bool, optional
			Whether to include global/overview maps.
		
		Returns
		-------
		Set[str]
			Set of pathway descriptions for all known pathways.
		
		Raises
		------
		NoKnownPathwaysError
			If the organism has no known pathways.
		HTTPError
			If any other HTTP error occurs.
		URLError
			If connection to KEGG fails.
		"""
		descriptions = Database.getPathwayDescriptions(self.nameAbbreviation)
		if includeOverviewMaps == True:
			return descriptions
		else:
			return self.__class__._filterGlobalAndOverview(descriptions)
	
	
	@classmethod
	def _filterGlobalAndOverview(cls, pathwayDescriptions: Set[str]) -> Set[str]:
		"""
		Removes pathway descriptions of pathways belonging to global or overview maps.
		
		Parameters
		----------
		pathwayDescriptions : Set[str]
			The pathway descriptions to filter.
			
		Returns
		-------
		Set[str]
			Pathway descriptions, leaving only the ones **not** from a global/overview pathway.
		"""
		newSet = set()
		for pathwayString in pathwayDescriptions:
			if cls.GLOBAL_PATHWAY_PATTERN.search(pathwayString) is None: # not a global/overview map
				newSet.add(pathwayString)
		return newSet
	
	
	def getMetabolicPathwayDescriptions(self, includeOverviewMaps = False) -> Set[str]:
		"""
		Get descriptions of pathways that are part of metabolism.
		
		Parameters
		----------
		includeOverviewMaps : bool, optional
			Whether to include global/overview maps.
		
		Returns
		-------
		Set[str]
			Set of pathway descriptions for all known metabolic pathways.
		
		Raises
		------
		NoKnownPathwaysError
			If the organism has no known pathways.
		HTTPError
			If pathway description list should not exist. Which would be odd, if we are certain and tested that the organism itself exists.
		URLError
			If connection to KEGG fails.
		"""
		descriptions = self.getPathwayDescriptions(includeOverviewMaps)
		return self._filterNonMetabolic(descriptions)
	
	
	def _filterNonMetabolic(self, pathwayDescriptions: Set[str]) -> Set[str]:
		"""
		Removes pathway descriptions of pathways not belonging to metabolism.
		
		Parameters
		----------
		pathwayDescriptions : Set[str]
			The pathway descriptions to filter.
			
		Returns
		-------
		Set[str]
			Pathway descriptions, leaving only the ones from pathways belonging to metabolism.
		
		See Also
		--------
		FEV_KEGG.quirks.METABOLIC_PATHWAYS : List of names for all pathways belonging to metabolism.
		"""
		newSet = set()
		for pathwayString in pathwayDescriptions:
			
			tempSet = set()
			tempSet.add(pathwayString)
			
			if self.getPathwayNames(self.getPathwayIDs(tempSet)).pop() in quirks.METABOLIC_PATHWAYS: # is a metabolism pathway
				newSet.add(pathwayString)
		
		return newSet
	
	
	def getPathwayIDs(self, pathwayDescriptions: Set[str]) -> Set[str]:
		"""
		Get pathway IDs from a set of descriptions.
		
		A pathway ID is the occurence of a specific pathway in a specific organism, e.g. pathway '00260' in 'eco' -> 'eco00260'. 
		
		Parameters
		----------
		pathwayDescriptions : Set[str]
			The pathway descriptions to be searched, e.g. 'path:eco00260	Glycine, serine and threonine metabolism - Escherichia coli K-12 MG1655'.
		
		Returns
		-------
		Set[str]
			Pathway IDs, e.g. 'eco00260'.
		"""
		pathwayIDSet = set()
		for pathway in pathwayDescriptions:
			pathwayID = pathway.split('\t')[0].replace('path:','')
			pathwayIDSet.add(pathwayID)
		return pathwayIDSet
	
	
	def getPathwayNames(self, pathwayIDs: Set[str]) -> Set[str]:
		"""
		Get pathway names from a set of IDs.
		
		A pathway name is a specific pathway, independent from any organism.
		
		Parameters
		----------
		pathwayDescriptions : Set[str]
			The pathway IDs to be searched, e.g. 'eco00260'.
		
		Returns
		-------
		Set[str]
			Pathway name, e.g. '00260'.
		"""
		pathwayNameSet = set()
		for pathwayID in pathwayIDs:
			pathwayNameSet.add(pathwayID.replace(self.nameAbbreviation, ''))
		return pathwayNameSet
	
	
	
	# Gene ====================================================================================
	
	
	
	def getGene(self, gene: 'eco:b0004 or b0004') -> Gene:
		"""
		Get a certain gene object for this organism.
		
		Automatically recognises format.
		
		Parameters
		----------
		gene : str
			Gene ID or name, e.g. 'eco:b0004' or 'b0004'.
		
		Returns
		-------
		Gene
			Gene object.
			
		Raises
		------
		HTTPError
			If gene does not exist.
		URLError
			If connection to KEGG fails.
		"""	
		if ':' in gene:
			return self.getGeneByID(gene)
		else:
			return self.getGeneByName(gene)
	
	
	def getGeneByName(self, geneName: 'b0004') -> Gene:
		"""
		Get a certain gene object for this organism.
		
		Automatically prepends organism, eg. 'eco:'+geneName.
		
		Parameters
		----------
		geneName : str
			Gene name, e.g. 'b0004'.
		
		Returns
		-------
		Gene
			Gene object.
			
		Raises
		------
		HTTPError
			If gene does not exist.
		URLError
			If connection to KEGG fails.
		"""
		gene = Database.getGene(self.nameAbbreviation + ':' + geneName)
		return gene
	
	
	def getGeneByID(self, geneID: 'eco:b0004') -> Gene:
		"""
		Get a certain gene object for this organism.
		
		Does not check if the prefix matches this organism!
		
		Parameters
		----------
		geneID : str
			Gene name, e.g. 'eco:b0004'.
		
		Returns
		-------
		Gene
			Gene object.
			
		Raises
		------
		HTTPError
			If gene does not exist.
		URLError
			If connection to KEGG fails.
		"""
		gene = Database.getGene(geneID)
		return gene
	
	
	def getGeneIDs(self, pathway: 'KGML_pathway.Pathway or 00260') -> Set[str]:
		"""
		Get the set of all gene IDs of this organism in a certain pathway.
		
		Automatically chooses :func:`getGeneIDsByName` or :func:`getGeneIDsByPathway`, depending on the type of `pathway`. Deduplicates original list.
		
		Parameters
		----------
		pathway : Pathway or str
			The pathway to search, either as :class:`FEV_KEGG.lib.Biopython.KEGG.KGML.KGML_pathway.Pathway` object or its name as a string, e.g. '00260'.
		
		Returns
		-------
		Set[str]
			List of gene IDs in `pathway`, e.g. ['eco:b0632', 'eco:b0839', 'eco:b2010'].
			
		Raises
		------
		HTTPError
			If name passed and pathway does not exist.
		URLError
			If name passed and connection fails.
		"""
		if pathway.__class__ == KGML_pathway.Pathway:
			return self.getGeneIDsByPathway(pathway)
		else:
			return self.getGeneIDsByName(pathway)
	
	
	def getGeneIDsByPathway(self, pathway: KGML_pathway.Pathway) -> Set[str]:
		"""
		Get the set of all gene IDs of this organism for a :class:`FEV_KEGG.lib.Biopython.KEGG.KGML.KGML_pathway.Pathway` object.
		
		Deduplicates original list.
		
		Parameters
		----------
		pathway : Pathway
			The pathway to search.
			
		Returns
		-------
		Set[str]
			List of gene IDs in `pathway`, e.g. ['eco:b0632', 'eco:b0839', 'eco:b2010'].
		"""
		pathwayNumber = self.__class__._DIGITS_PATTERN.findall(pathway.name)[0]
		geneNameList = Database.getPathwayGeneIDs(self.nameAbbreviation, pathwayNumber) # try to get list from disk
		
		# if not on disk, calculate the list
		if geneNameList == None:
			geneNameList = self._calculateGeneIDs(pathway, pathwayNumber)
			
		return geneNameList
		
	
	def getGeneIDsByName(self, pathwayName: '00260') -> Set[str]:
		"""
		Get the set of all gene IDs of this organism for a pathway name.
		
		Deduplicates original list.
		
		Parameters
		----------
		pathwayName : str
			The pathway name to search, e.g. '00260'.
		
		Returns
		-------
		Set[str]
			List of gene IDs in `pathwayName`, e.g. ['eco:b0632', 'eco:b0839', 'eco:b2010'].
		
		Raises
		------
		HTTPError
			If pathway does not exist.
		URLError
			If connection to KEGG fails.
		"""
		geneNameList = Database.getPathwayGeneIDs(self.nameAbbreviation, pathwayName) # try to get list from disk
		
		# if not on disk, calculate the list
		if geneNameList == None:
			geneNameList = self._calculateGeneIDs(self.getPathway(pathwayName), pathwayName)
			
		return geneNameList
	
	
	def _calculateGeneIDs(self, pathway: KGML_pathway.Pathway, pathwayName: '00260') -> Set[str]:
		geneNameSet = set()
		geneDescriptionList = pathway.genes
		for geneDescription in geneDescriptionList:
			namePossiblyList = geneDescription.name
			
			nameList = namePossiblyList.split(' ')
			
			for name in nameList:
				geneNameSet.add(name)
		
		Database.setPathwayGeneIDs(self.nameAbbreviation, pathwayName, geneNameSet) # cache the deduplicated list to disk
		
		return geneNameSet
			
	
	def getGenes(self, pathway: 'KGML_pathway.Pathway or 00260') -> Set[Gene]:
		"""		
		Get the set of all genes of this organism in a certain pathway.
		
		Automatically chooses :func:`getGenesByName` or :func:`getGenesByPathway`, depending on the type of `pathway`. Deduplicates original list.
		
		Parameters
		----------
		pathway : Pathway or str
			The pathway to search, either as :class:`FEV_KEGG.lib.Biopython.KEGG.KGML.KGML_pathway.Pathway` object or its name as a string, e.g. '00260'.
		
		Returns
		-------
		Set[Gene]
			List of gene IDs in `pathway`, e.g. ['eco:b0632', 'eco:b0839', 'eco:b2010'].
			
		Raises
		------
		HTTPError
			If name passed and pathway does not exist.
		URLError
			If name passed and connection fails.
		"""
		if pathway.__class__ == KGML_pathway.Pathway:
			return self.getGenesByPathway(pathway)
		else:
			return self.getGenesByName(pathway)
	
	
	def getGenesByName(self, pathwayName: '00260') -> Set[Gene]:
		"""		
		Get the set of all genes of this organism for a pathway name.
		
		Deduplicates original list.
		
		Parameters
		----------
		pathwayName : str
			The pathway name to search, e.g. '00260'.
		
		Returns
		-------
		Set[Gene]
			List of gene IDs in `pathwayName`, e.g. ['eco:b0632', 'eco:b0839', 'eco:b2010'].
		
		Raises
		------
		HTTPError
			If pathway does not exist.
		URLError
			If connection to KEGG fails.
		"""
		geneList = set()
		
		geneIDList = self.getGeneIDsByName(pathwayName)
		
		for geneID in geneIDList:
			geneList.add(self.getGeneByID(geneID))
			
		return geneList
	
	
	def getGenesByPathway(self, pathway: KGML_pathway.Pathway) -> Set[Gene]:
		"""		
		Get the set of all genes of this organism for a :class:`FEV_KEGG.lib.Biopython.KEGG.KGML.KGML_pathway.Pathway` object.
		
		Deduplicates original list.
		
		Parameters
		----------
		pathway : Pathway
			The pathway to search.
			
		Returns
		-------
		Set[Gene]
			List of gene IDs in `pathway`, e.g. ['eco:b0632', 'eco:b0839', 'eco:b2010'].
		"""
		geneList = set()
		
		geneIDList = self.getGeneIDsByPathway(pathway)
		
		for geneID in geneIDList:
			geneList.add(self.getGeneByID(geneID))
			
		return geneList
	
	
	def getNumberOfGenes(self):
		"""
		Get the number of known genes within this genome.
		
		Returns
		-------
		int
			Count of genes in this organism's genome. These do not necessarily have to be mentioned in any pathway.
		"""
		organismInfo = Database.getOrganismInfo(self.nameAbbreviation, checkExpiration = False)
		return Database._extractGeneEntries(organismInfo)
	
	
	# Graph ====================================================================================
	
	
	
	def substanceEcGraph(self, noMultifunctional = settings.defaultNoMultifunctional, returnCacheEntry = False) -> SubstanceEcGraph:
		"""
		Substance-EC graph of this organism.
		
		Parameters
		----------
		noMultifunctional : bool, optional
			If *True*, ignore enzymes with multiple EC numbers.
		returnCacheEntry : bool, optional
			If *True*, do not return the graph, but instead a :class:`FEV_KEGG.KEGG.File.CacheEntry`. This cache entry can be useful for parallel computation. 
		
		Returns
		-------
		SubstanceEcGraph
			The graph has substrates/products as its nodes and EC numbers as the connecting edges. Edges have a direction.
			
		Raises
		------
		HTTPError
			If any gene or pathway does not exist.
		URLError
			If connection to KEGG fails.
		"""
		file_name = 'SubstanceEcGraph'
		if noMultifunctional is True:
			file_name += '_noMultifunctional'
		
		folder_path = 'organism/' + self.nameAbbreviation + '/graph'
		
		if returnCacheEntry is False: # shall return result
			decorator = cache(folder_path = folder_path, file_name = file_name)
		else: # shall return CacheEntry object
			decorator = cacheEntry(folder_path = folder_path, file_name = file_name)
			
		function = lambda: Conversion.SubstanceEnzymeGraph2SubstanceEcGraph(self.substanceEnzymeGraph(noMultifunctional))
		return decorator(function)()
	
	def substanceEnzymeGraph(self, noMultifunctional = settings.defaultNoMultifunctional, returnCacheEntry = False) -> SubstanceEnzymeGraph:
		"""
		Substance-Enzyme graph of this organism.
		
		Parameters
		----------
		noMultifunctional : bool, optional
			If *True*, ignore enzymes with multiple EC numbers.
		returnCacheEntry : bool, optional
			If *True*, do not return the graph, but instead a :class:`FEV_KEGG.KEGG.File.CacheEntry`. This cache entry can be useful for parallel computation. 
		
		Returns
		-------
		SubstanceEnzymeGraph
			The graph has substrates/products as its nodes and enzymes as the connecting edges. Edges have a direction.
			
		Raises
		------
		HTTPError
			If any gene or pathway does not exist.
		URLError
			If connection to KEGG fails.
		"""
		file_name = 'SubstanceEnzymeGraph'
		if noMultifunctional is True:
			file_name += '_noMultifunctional'
		
		folder_path = 'organism/' + self.nameAbbreviation + '/graph'
		
		if returnCacheEntry is False: # shall return result
			decorator = cache(folder_path = folder_path, file_name = file_name)
		else: # shall return CacheEntry object
			decorator = cacheEntry(folder_path = folder_path, file_name = file_name)
		
		function = lambda: Conversion.SubstanceGeneGraph2SubstanceEnzymeGraph(self.substanceGeneGraph(), noMultifunctional)
		return decorator(function)()
	
	def substanceGeneGraph(self, returnCacheEntry = False) -> SubstanceGeneGraph:
		"""
		Substance-Gene graph of this organism.
		
		Parameters
		----------
		returnCacheEntry : bool, optional
			If *True*, do not return the graph, but instead a :class:`FEV_KEGG.KEGG.File.CacheEntry`. This cache entry can be useful for parallel computation. 
		
		Returns
		-------
		SubstanceGeneGraph
			The graph has substrates/products as its nodes and genes as the connecting edges. Edges have a direction.
			
		Raises
		------
		HTTPError
			If any gene or pathway does not exist.
		URLError
			If connection to KEGG fails.
		"""
		file_name = 'SubstanceGeneGraph'
		folder_path = 'organism/' + self.nameAbbreviation + '/graph'
		
		if returnCacheEntry is False: # shall return result
			decorator = cache(folder_path = folder_path, file_name = file_name)
		else: # shall return CacheEntry object
			decorator = cacheEntry(folder_path = folder_path, file_name = file_name)
		
		function = lambda: Conversion.SubstanceReactionGraph2SubstanceGeneGraph(self.substanceReactionGraph())
		return decorator(function)()
	
	def substanceReactionGraph(self, returnCacheEntry = False) -> SubstanceReactionGraph:
		"""
		Substance-Reaction graph of this organism.
		
		Parameters
		----------
		returnCacheEntry : bool, optional
			If *True*, do not return the graph, but instead a :class:`FEV_KEGG.KEGG.File.CacheEntry`. This cache entry can be useful for parallel computation. 
		
		Returns
		-------
		SubstanceReactionGraph
			The graph has substrates/products as its nodes and reactions as the connecting edges. Edges have a direction.
			
		Raises
		------
		HTTPError
			If any gene or pathway does not exist.
		URLError
			If connection to KEGG fails.
		"""
		file_name = 'SubstanceReactionGraph'
		folder_path = 'organism/' + self.nameAbbreviation + '/graph'
		
		if returnCacheEntry is False: # shall return result
			decorator = cache(folder_path = folder_path, file_name = file_name)
		else: # shall return CacheEntry object
			decorator = cacheEntry(folder_path = folder_path, file_name = file_name)
		
		function = lambda: Conversion.KeggPathwaySet2SubstanceReactionGraph(self.getMetabolicPathways(), name = self.nameAbbreviation)
		return decorator(function)()
	






	
	
class Group(object):
	
	def __init__(self, organismAbbreviations: Iterable[str] = None, searchString: 'any part of organism description' = None, name = None, minimalSize = None):
		"""
		A Group of :class:`Organism`.
		
		If both parameters, `organismAbbreviations` and `searchString`, are specified, both lists will be appended, forming this group's list of organisms.
		If none of the parameters is specified, this Group has an empty organism list.
		
		Parameters
		----------
		organismAbbreviations : Iterable[str], optional
			Abbreviations of the desired organisms. If != *None*, tries to find one :class:`Organism` for each abbreviation in the list.
		searchString : str, optional
			Any part of the desired organisms' description. If != *None*, searches the list of all organisms known to KEGG for the passed string.
			An example entry of the KEGG list of organisms looks as follows: *"T00338	eci	Escherichia coli O18:K1:H7 UTI89 (UPEC)	Prokaryotes;Bacteria;Gammaproteobacteria - Enterobacteria;Escherichia"*.
			Any list entry matching the search string creates one :class:`Organism`, aggregated into this group.
		name : str, optional
			Custom name of this group.
		minimalSize : int, optional
			If not *None*, incorporate only organisms with EC graphs with at least `minimalSize` edges. Can be useful to filter incompletely annotated organisms.
		
		Attributes
		----------
		self.searchString : str
		self.name : str
		
		Raises
		------
		ValueError
			If any organism does not exist at all in KEGG.
		URLError
			If connection to KEGG fails.
		"""
		self._collectiveEcGraph = None
		self._collectiveEcGraph_noMultifunctional = None
		self._collectiveEnzymeGraph = None
		self._collectiveEnzymeGraph_noMultifunctional = None
		
		self.searchString = searchString
		self.__organisms = set()
		
		if searchString is not None:
			organismList = Database.getOrganismList()
			
			matchList = []
			
			for entry in organismList:
				if searchString in entry:
					entrySplit = entry.split('\t')
					matchList.append(entrySplit[1])
			
			organisms = Organism.__initBulk__(matchList)
			self.__organisms.update(organisms)
		
		if organismAbbreviations is not None:
			listObject = []
			listObject.extend(organismAbbreviations)
			organisms = Organism.__initBulk__(listObject)
			self.__organisms.update(organisms)
		
		self.name = name
		
		if minimalSize is not None: # count sizes of each organisms graph
			
			organismsToRemove = set()
			for organism in organisms:
				ecGraph = organism.substanceEcGraph(noMultifunctional = False)
				if len(ecGraph.getEdges()) < minimalSize:
					organismsToRemove.add(organism)
			
			self.__organisms.difference_update( organismsToRemove )
	
	def __str__(self):
		return 'Group(' + self.organisms + ')'

	def __eq__(self, other):
		if isinstance(self, other.__class__):
			return self.organisms == other.organisms
		return False

	def __ne__(self, other):
		return not self == other

	def __hash__(self):
		return self.organisms.__hash__()
	
	def freeHeap(self):
		"""
		Free heap of memory-cached data.
		
		Removes objects kept on heap, i.e. which have a pointer kept in this object, because some function was called with `keepOnHeap` == *True*.
		Also calls garbage collector to break reference cycles in previously uncollected objects (generation 0).
		
		Note
		----
		Instead of using this, you will most likely want to never use any of the group methods with `keepOnHeap` == *True*. Then, calling this method would be unnecessary, as it would have no effect.
		"""
		self._collectiveEcGraph = None
		self._collectiveEcGraph_noMultifunctional = None
		self._collectiveEnzymeGraph = None
		self._collectiveEnzymeGraph_noMultifunctional = None
		gc.collect(0)
	
	
	@property
	def organisms(self) -> Set[Organism]:
		"""
		Organisms of this group.
		
		Returns
		-------
		Set[Organism]
			The set of organisms which are part of this group. Order is arbitrary.
		"""
		return self.__organisms
	
	@property
	def organismsCount(self) -> int:
		"""
		Number of organisms of this group.
		
		Returns
		-------
		int
			The number of organisms in the set of organisms of this group.
		"""
		return len( self.__organisms )
	
	
	
	
	def enzymeGraphs(self, noMultifunctional = settings.defaultNoMultifunctional) -> List[SubstanceEnzymeGraph]:
		"""
		All substance-enzyme graphs of this group.
		
		Parameters
		----------
		noMultifunctional : bool, optional
			If *True*, ignore enzymes with multiple EC numbers.
		
		Returns
		-------
		List[SubstanceEnzymeGraph]
			Substance-enzyme graphs of all organisms in this group. Order is arbitrary.
		
		Raises
		------
		TypeError
			If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`_getGraphsParallelly`.
		HTTPError
			If fetching any of the underlying graphs fails.
		URLError
			If connection to KEGG fails.
		"""
		return self._getGraphsParallelly(self._enzymeGraphsWorker, self.organisms, noMultifunctional, 'enzyme graphs')
	
	def _enzymeGraphsWorker(self, organism: Organism, noMultifunctional, returnCacheEntry) -> SubstanceEnzymeGraph:
		return organism.substanceEnzymeGraph(noMultifunctional, returnCacheEntry)
	
	def ecGraphs(self, noMultifunctional = settings.defaultNoMultifunctional, minimalSize = None) -> List[SubstanceEcGraph]:
		"""
		All substance-EC graphs of this group.
		
		Parameters
		----------
		noMultifunctional : bool, optional
			If *True*, ignore enzymes with multiple EC numbers.
		minimalSize : int, optional
			If not *None*, return only EC graphs with at least `minimalSize` edges. Can be useful to filter incompletely annotated organisms.
		
		Returns
		-------
		List[SubstanceEcGraph]
			Substance-EC graphs of all organisms in this group. Order is arbitrary.
		
		Raises
		------
		TypeError
			If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`_getGraphsParallelly`.
		HTTPError
			If fetching any of the underlying graphs fails.
		URLError
			If connection to KEGG fails.
		"""
		return self._getGraphsParallelly(self._ecGraphsWorker, self.organisms, noMultifunctional, 'EC graphs', minimalSize = minimalSize)
	
	def _ecGraphsWorker(self, organism: Organism, noMultifunctional, returnCacheEntry) -> SubstanceEcGraph:
		return organism.substanceEcGraph(noMultifunctional, returnCacheEntry)
	
	def _getGraphsParallelly(self, worker, organisms, noMultifunctional, debugText, minimalSize = None):
		"""
		Does the actual fetching and computing of the graphs in parallel.
		
		Raises
		------
		TypeError
			If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`.
		
		Warnings
		--------
		If you did not enable parallel computation by either enabling :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or providing :attr:`FEV_KEGG.Util.Parallelism.processPool`, this will fail with a TypeError.
		
		Note
		----
		If you enabled :attr:`FEV_KEGG.settings.automaticallyStartProcessPool`, this will run parallelly in multiple processes with multiple threads in each, depending on your settings and computer.
		"""
		threadPool = concurrent.futures.ThreadPoolExecutor(Parallelism.getNumberOfThreadsFile())
		futures = []
		futuresIO = []
		futuresGenerator = None
		resultFutures = None
		
		try:
			
			# submit work to process pool
			for organism in organisms:
				if Parallelism.processPool is None:
					raise TypeError("Process pool does not exist. Did you forget to FEV_KEGG.startProcessPool()?")
				futures.append( Parallelism.processPool.submit( worker, organism, noMultifunctional, True ) )
				
			futuresGenerator = concurrent.futures.as_completed( futures )
			
			# add progress bar
			if settings.verbosity >= 1:
				if settings.verbosity >= 2:
					print( 'Fetching ' + debugText + ' from ' + str(len(organisms)) + ' organisms...' )
				futuresGenerator = tqdm.tqdm(futuresGenerator, total = len(organisms), unit = ' organisms', position = 0)
			
			# when any work item in process pool finishes
			for future in futuresGenerator:
				
				try:
					cacheEntry = future.result()
				except KeyboardInterrupt:
					raise
				except concurrent.futures.CancelledError:
					Parallelism.printBelowProgress( "Future cancelled. Continuing anyway..." )
					continue
				except concurrent.futures.TimeoutError:
					Parallelism.printBelowProgress( "Future timed out. Continuing anyway..." )
					continue
				except NoKnownPathwaysError as e: # organism has no known pathways, ignore it
					Parallelism.printBelowProgress( "Future raised error: " + str(e) + " Ignoring this organism." )
					continue
				except Exception: # any non-exiting error
					Parallelism.printBelowProgress( "Future raised error, see stack trace above. Halting by KeyboardInterrupt..." )
					raise KeyboardInterrupt()
				
				futuresIO.append( threadPool.submit(cacheEntry.getResult) )
			
			resultFutures = concurrent.futures.as_completed( futuresIO )
			
			if settings.verbosity >= 2:
				if settings.verbosity >= 2:
					print( 'Doing I/O for ' + str(len(organisms)) + ' organisms...' )
				resultFutures = tqdm.tqdm(resultFutures, total = len(organisms), unit = ' organism I/Os', position = 0)
			
			graphs = []
			for future in resultFutures:
				graph = future.result()
				if minimalSize is None or len(graph.getEdges()) >= minimalSize:
					graphs.append( graph )
			
			return graphs
		
		except KeyboardInterrupt: # only raised in main thread (once in each process!)
			
			Parallelism.keyboardInterruptHandler(processPoolFutures=futures, threadPool=threadPool, threadPoolFutures=futuresIO, terminateProcess=True)
			raise

		except BaseException:
			
			if Parallelism.isMainThread():
				Parallelism.keyboardInterruptHandler(processPoolFutures=futures, threadPool=threadPool, threadPoolFutures=futuresIO, silent=True)
			raise
		
		finally:
			
			if threadPool is not None: threadPool.shutdown(wait = False)
			if futuresGenerator is not None: futuresGenerator.close()
			if resultFutures is not None: resultFutures.close()
			
			Parallelism.printBelowProgress(None)

	
	
	
	def collectiveEcGraph(self, noMultifunctional = settings.defaultNoMultifunctional, addCount = False, keepOnHeap = True, addEcDescriptions = False) -> SubstanceEcGraph:
		"""
		The collective of all SubstanceEcGraphs, by union operation, from all organisms in this group.
		
		Nodes of the same Substance are merged, all edges of differing ECs with a unique pair of nodes survive.
		
		Parameters
		----------
		noMultifunctional : bool, optional
			If *True*, ignore enzymes with multiple EC numbers.
		addCount : bool, optional
			If *True*, the returned graph contains extra dicts.
			1 ``graph.nodeCounts[node]`` = number of organisms which contained this node.
			2 ``graph.edgeCounts[(node, node, element)]`` = number of organisms which contained this edge.
			3 ``graph.edgeElementCounts[element]`` = number of organisms which contained this element.
			Attention! These counter dictionaries are **NOT** updated if your add or remove a node/edge/element!
		keepOnHeap : bool, optional
			Keeps a common graph in memory to speed up subsequent calls of this or other methods.
			This can take up a lot of memory! Once this object is garbage collected, the common graph will be, too.
		
		Returns
		-------
		SubstanceEcGraph
			The substance-EC graph composed of all this group's organism's substance-EC graphs.
		
		Raises
		------
		TypeError
			If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`_getGraphsParallelly`.
		HTTPError
			If fetching any of the underlying graphs fails.
		URLError
			If connection to KEGG fails.
		
		See Also
		--------
		FEV_KEGG.Graph.Models.CommonGraphApi.union : Union operator.
		"""
		if keepOnHeap is True: # shall keep the result on heap
			# check first if it already is on heap
			if noMultifunctional is True:
				collectiveEcGraph = self._collectiveEcGraph_noMultifunctional
			else:
				collectiveEcGraph = self._collectiveEcGraph
			
			if collectiveEcGraph is not None:
				return collectiveEcGraph.copy()
		
		if settings.verbosity >= 1:
			print('calculating collective EC graph...')
			
		allSubstanceEcGraphs = self.ecGraphs(noMultifunctional)
		if isinstance(allSubstanceEcGraphs, Generator):
			lastGraph = allSubstanceEcGraphs.next()
		else:
			lastGraph = allSubstanceEcGraphs.pop()
		
		collectiveGraph = lastGraph.union(allSubstanceEcGraphs, addCount = addCount)
		collectiveGraph.removeIsolatedNodes()
		
		if self.name is None:
			collectiveGraph.name = 'Collective EC graph'
		else:
			collectiveGraph.name = 'Collective EC graph ' + self.name
			
		if addEcDescriptions is True:
			collectiveGraph.addEcDescriptions()
		
		if keepOnHeap is True: # shall keep the result on heap
			# save the calculated result on heap
			if noMultifunctional is True:
				self._collectiveEcGraph_noMultifunctional = collectiveGraph
			else:
				self._collectiveEcGraph = collectiveGraph
		
		return collectiveGraph.copy()
			
		
	
	def consensusEcGraph(self, noMultifunctional = settings.defaultNoMultifunctional, keepOnHeap = True) -> SubstanceEcGraph:
		"""
		The consensus of all SubstanceEcGraphs, by intersection operation, from all organisms in this group.
		
		Afterwards, removes isolated nodes.
		
		Parameters
		----------
		noMultifunctional : bool, optional
			If *True*, ignore enzymes with multiple EC numbers.
		keepOnHeap : bool, optional
			Keeps a common graph in memory to speed up subsequent calls of this or other methods.
			This can take up a lot of memory! Once this object is garbage collected, the common graph will be, too.
		
		Returns
		-------
		SubstanceEcGraph
			The substance-EC graph intersected of all this group's organism's substance-EC graphs.
		
		Raises
		------
		TypeError
			If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`_getGraphsParallelly`.
		HTTPError
			If fetching any of the underlying graphs fails.
		URLError
			If connection to KEGG fails.
		
		See Also
		--------
		FEV_KEGG.Graph.Models.CommonGraphApi.intersection : Intersection operator.
		"""		
		if keepOnHeap is True:
			
			collectiveEcGraph = self.collectiveEcGraph(noMultifunctional, addCount = True, keepOnHeap = True) 
			edgeCounts = collectiveEcGraph.edgeCounts
			
			edgesToBeRemoved = []
			for edge, count in edgeCounts.items():
				if count < self.organismsCount: # has not occured in enough organisms
					edgesToBeRemoved.append(edge)
			
			collectiveEcGraph.removeEcEdges(edgesToBeRemoved)
			collectiveEcGraph.removeIsolatedNodes()
			
			if self.name is None:
				collectiveEcGraph.name = 'Consensus EC graph'
			else:
				collectiveEcGraph.name = 'Consensus EC graph ' + self.name
			
			return collectiveEcGraph
			
		else:
			
			allSubstanceEcGraphs = self.ecGraphs(noMultifunctional)
			if isinstance(allSubstanceEcGraphs, Generator):
				lastGraph = allSubstanceEcGraphs.next()
			else:
				lastGraph = allSubstanceEcGraphs.pop()
			consensusGraph = lastGraph.intersection(allSubstanceEcGraphs)
			consensusGraph.removeIsolatedNodes()
			
			if self.name is None:
				consensusGraph.name = 'Consensus EC graph'
			else:
				consensusGraph.name = 'Consensus EC graph ' + self.name
			
			return consensusGraph
		
	def majorityEcGraph(self, majorityPercentage = 90, majorityTotal = None, noMultifunctional = settings.defaultNoMultifunctional, keepOnHeap = True) -> SubstanceEcGraph:
		"""
		The majority-consensus of all SubstanceEcGraphs, by majority-intersection operation, from all organisms in this group.
		
		If the majority of organisms contains an edge, it is added to the majority-consensus. Afterwards, removes isolated nodes.
		
		Parameters
		----------
		majorityPercentage : float, optional
			Majority percentage means 'at least x%' and is rounded up. For example 90% of 11 organisms would be ceiling(9,9) = 10 organisms.
		majorityTotal : int, optional
			If given (not *None*), `majorityPercentage` is ignored and the percentage of organisms for a majority is calculated from `majorityTotal` alone.
		noMultifunctional : bool, optional
			If *True*, ignore enzymes with multiple EC numbers.
		keepOnHeap : bool, optional
			Keeps a common graph in memory to speed up subsequent calls of this or other methods.
			This can take up a lot of memory! Once this object is garbage collected, the common graph will be, too.
		
		Returns
		-------
		SubstanceEcGraph
			The substance-EC graph majority-intersected of all this group's organism's substance-EC graphs.
		
		Raises
		------
		TypeError
			If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`_getGraphsParallelly`.
		HTTPError
			If fetching any of the underlying graphs fails.
		URLError
			If connection to KEGG fails.
		
		See Also
		--------
		FEV_KEGG.Graph.Models.CommonGraphApi.majorityIntersection : Majority intersection operator.
		"""		
		if majorityTotal is not None:
			percentage = majorityTotal / self.organismsCount * 100
		else:
			percentage = majorityPercentage
			
		if keepOnHeap is True:
			
			# check if majorityPercentage is sane
			if majorityPercentage <= 0 or majorityPercentage > 100:
				raise ValueError('Majority percentage is not a sane value (0 < percentage <= 100): ' + str(majorityPercentage))
			majorityTotal = ceil((majorityPercentage/100) * self.organismsCount)
			
			collectiveEcGraph = self.collectiveEcGraph(noMultifunctional, addCount = True, keepOnHeap = True)
			edgeCounts = collectiveEcGraph.edgeCounts
			
			edgesToBeRemoved = []
			for edge, count in edgeCounts.items():
				if count < majorityTotal: # has not occured in enough organisms
					edgesToBeRemoved.append(edge)
			
			collectiveEcGraph.removeEcEdges(edgesToBeRemoved)
			collectiveEcGraph.removeIsolatedNodes()
			
			if self.name is None:
				collectiveEcGraph.name = 'Majority EC graph'
			else:
				collectiveEcGraph.name = 'Majority EC graph ' + self.name
			
			return collectiveEcGraph
		
		else:
			
			allSubstanceEcGraphs = self.ecGraphs(noMultifunctional)
			if isinstance(allSubstanceEcGraphs, Generator):
				lastGraph = allSubstanceEcGraphs.next()
			else:
				lastGraph = allSubstanceEcGraphs.pop()
			majorityGraph = lastGraph.majorityIntersection(allSubstanceEcGraphs, percentage)
			majorityGraph.removeIsolatedNodes()
			
			if self.name is None:
				majorityGraph.name = 'Majority EC graph'
			else:
				majorityGraph.name = 'Majority EC graph ' + self.name
			
			return majorityGraph
	
	
	
	
	
	def collectiveEnzymeGraph(self, noMultifunctional = settings.defaultNoMultifunctional, keepOnHeap = True) -> SubstanceEnzymeGraph:
		"""
		The collective of all SubstanceEnzymeGraphs, by union operation, from all organisms in this group.
		
		Nodes of the same Substance are merged, all edges of differing Enzymes with a unique pair of nodes survive. Enzymes are compared by their GeneID and should, thus, all be different!
		
		Parameters
		----------
		noMultifunctional : bool, optional
			If *True*, ignore enzymes with multiple EC numbers.
		keepOnHeap : bool, optional
			Keeps a common graph in memory to speed up subsequent calls of this or other methods.
			This can take up a lot of memory! Once this object is garbage collected, the common graph will be, too.
		
		Returns
		-------
		SubstanceEnzymeGraph
			The substance-enzyme graph composed of all this group's organism's substance-enzyme graphs.
		
		Raises
		------
		TypeError
			If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`_getGraphsParallelly`.
		HTTPError
			If fetching any of the underlying graphs fails.
		URLError
			If connection to KEGG fails.
		
		See Also
		--------
		FEV_KEGG.Graph.Models.CommonGraphApi.union : Union operator.
		"""
		if keepOnHeap is True: # shall keep the result on heap
			# check first if it already is on heap
			if noMultifunctional is True:
				collectiveEnzymeGraph = self._collectiveEnzymeGraph_noMultifunctional
			else:
				collectiveEnzymeGraph = self._collectiveEnzymeGraph
			
			if collectiveEnzymeGraph is not None:
				return collectiveEnzymeGraph.copy()
		
		
		if settings.verbosity >= 1:
			print('calculating collective enzyme graph...')
		
		allSubstanceEnzymeGraphs = self.enzymeGraphs(noMultifunctional)
		if isinstance(allSubstanceEnzymeGraphs, Generator):
			lastGraph = allSubstanceEnzymeGraphs.next()
		else:
			lastGraph = allSubstanceEnzymeGraphs.pop()
		collectiveGraph = lastGraph.union(allSubstanceEnzymeGraphs)
		
		if self.name is None:
			collectiveGraph.name = 'Collective Enzyme graph'
		else:
			collectiveGraph.name = 'Collective Enzyme graph ' + self.name
			
		
		if keepOnHeap is True: # shall keep the result on heap
			# save the calculated result on heap
			if noMultifunctional is True:
				self._collectiveEnzymeGraph_noMultifunctional = collectiveGraph
			else:
				self._collectiveEnzymeGraph = collectiveGraph
		
		return collectiveGraph.copy()
	

	def collectiveEnzymeGraphByEcConsensus(self, noMultifunctional = settings.defaultNoMultifunctional) -> SubstanceEnzymeGraph:
		"""
		The collective SubstanceEnzymGraph, but containing only Enzymes whose EC numbers occur in the consensus of all SubstanceEcGraphs.
		
		Parameters
		----------
		noMultifunctional : bool, optional
			If *True*, ignore enzymes with multiple EC numbers.
		
		Raises
		------
		TypeError
			If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`_getGraphsParallelly`.
		HTTPError
			If fetching any of the underlying graphs fails.
		URLError
			If connection to KEGG fails.
		"""
		return self.collectiveEnzymeGraphByEcMajority(majorityPercentage = 100, majorityTotal = None, noMultifunctional = noMultifunctional)
		
	
	def collectiveEnzymeGraphByEcMajority(self, majorityPercentage = 80, majorityTotal = None, noMultifunctional = settings.defaultNoMultifunctional) -> SubstanceEnzymeGraph:
		"""
		The collective SubstanceEnzymGraph, but containing only Enzymes whose EC numbers occur in the majority of all SubstanceEcGraphs.
		
		Parameters
		----------
		majorityPercentage : float, optional
			Majority percentage means 'at least x%' and is rounded up. For example 90% of 11 organisms would be ceiling(9,9) = 10 organisms.
		majorityTotal : int, optional
			If given (not *None*), `majorityPercentage` is ignored and the percentage of organisms for a majority is calculated from `majorityTotal` alone.
		noMultifunctional : bool, optional
			If *True*, ignore enzymes with multiple EC numbers.
		
		Raises
		------
		TypeError
			If you failed to enable :attr:`FEV_KEGG.settings.automaticallyStartProcessPool` or to provide a :attr:`FEV_KEGG.Util.Parallelism.processPool`. See :func:`_getGraphsParallelly`.
		HTTPError
			If fetching any of the underlying graphs fails.
		URLError
			If connection to KEGG fails.
		"""
		if majorityPercentage >= 100:
			ecNumbers = self.consensusEcGraph(noMultifunctional).getECs()
			
		else:
			ecNumbers = self.majorityEcGraph(majorityPercentage, majorityTotal, noMultifunctional).getECs()
		
		collectiveEnzymeGraph = self.collectiveEnzymeGraph(noMultifunctional)
		
		collectiveEnzymeGraph.keepEnzymesByEC(ecNumbers)

		collectiveEnzymeGraph.name += ' by EC majority'
		
		return collectiveEnzymeGraph
		
	