from enum import Enum
import re

import anytree
from anytree.node.node import Node
from typing import Iterable, List

from FEV_KEGG.KEGG import Database
from FEV_KEGG.KEGG.File import cache
from FEV_KEGG import settings


class TaxonType(Enum):
    """
    Type of a taxon.
    """
    ROOT = 0
    """
    Root taxon, i.e. '/'.
    """
    ORGANISM = 1
    """
    Organism taxon, i.e. a leaf with a unique sequenced genome.
    """
    SPECIES = 2
    """
    Species taxon, e.g. 'Escherichia Coli'.
    """
    OTHER = 3
    """
    Other taxon, i.e. any other taxonomic rank in between.
    """

class Taxonomy(object):
    
    def __init__(self, rawLines, isNCBI):
        """        
        Generic taxonomy of organisms in KEGG.
        
        Parameters
        ----------
        rawLines : List[str]
            List of lines making up the raw data of a known taxonomy, either NCBI or KEGG.
        isNCBI : bool
            If *True*, `rawLines` is parsed as NCBI taxonomy. If *False*, `rawLines` is parsed as KEGG taxonomy.
        
        Attributes
        ----------
        self.indexOnAbbreviation : Dict[str, :class:`anytree.node.node.Node`]
            Index to find a :class:`anytree.node.node.Node` for an organism abbreviation, with `.type` == :attr:`TaxonType.ORGANISM`. 
        self.tree : :class:`anytree.node.node.Node`
            The root node of the taxonomy, with `.type` == :attr:`TaxonType.ROOT`.
        """
        self.indexOnAbbreviation = dict()
        self.tree = self._parse(rawLines, isNCBI)
    
    def getOrganismNodeByAbbreviation(self, abbreviation: 'eco') -> Node:
        """
        Get node for an organism by its abbreviation.
        
        Parameters
        ----------
        abbreviation : str
            Abbreviation of the organism in KEGG.
        
        Returns
        -------
        :class:`anytree.node.node.Node`
            Node of `.type` == :attr:`TaxonType.ORGANISM` with `.abbreviation` == `abbreviation`. *None* if none can be found.
        """
        return self.indexOnAbbreviation.get(abbreviation, None)
    
    def getOrganismNodesByName(self, name: 'Escherichia', oneOrganismPerSpecies = settings.defaultOneOrganismPerSpecies) -> List[Node]:
        """
        Get nodes for organisms by a part of their name.
        
        Parameters
        ----------
        name : str
            Part of the name of the desired organism taxons. This does **not** search parts of the path! The name may be abbreviated, i.e. 'Escherichia' will match 'Escherichia Coli K-12 MG1655'.
        oneOrganismPerSpecies : bool, optional
            If *True*, return only the first organism node of each species node.
        
        Returns
        -------
        List[:class:`anytree.node.node.Node`]
            List of organism nodes containing `name` in their name attribute. *None* if none found.
        """
        if oneOrganismPerSpecies is True:
            organismNodes = []
            
            speciesNodes = self.searchNodesByName(name, TaxonType.SPECIES)
            for speciesNode in speciesNodes:
                organismNode = speciesNode.descendants[0]
                organismNodes.append(organismNode)
            
            return organismNodes
        
        else:
            return self.searchNodesByName(name, TaxonType.ORGANISM)
    
    def getOrganismNodesByPath(self, path: 'Gammaproteobacteria/Enterobacterales', exceptPaths: List['Gammaproteobacteria/unclassified'] = None, oneOrganismPerSpecies = settings.defaultOneOrganismPerSpecies) -> List[Node]:
        """
        Get nodes for organisms by a part of their `path`.
        
        Parameters
        ----------
        path : str
            Part of the path of the desired organism taxons. The parts of the path specified here have to match the wording of the path nodes exactly, i.e. 'Enterobac' will **not** match 'Enterobacterales'.
        exceptPaths : Iterable[str] or str
            Paths which match any of these will not be returned. Accepts iterables of exceptions or a single string exception.
        oneOrganismPerSpecies : bool, optional
            If *True*, return only the first organism node of each species node.
        
        Returns
        -------
        List[:class:`anytree.node.node.Node`]
            List of organism nodes containing `path` in their path. *None* if none found.
        """
        if oneOrganismPerSpecies is True:
            organismNodes = []
            
            speciesNodes = self.searchNodesByPath(path, TaxonType.SPECIES, exceptPaths)
            for speciesNode in speciesNodes:
                try:
                    organismNode = speciesNode.descendants[0]
                except IndexError:
                    continue
                organismNodes.append(organismNode)
            
            return organismNodes
        
        else:
            return self.searchNodesByPath(path, TaxonType.ORGANISM, exceptPaths)
    
    def getOrganismAbbreviations(self, nodes: Iterable[Node]) -> List[str]:
        """
        Get abbreviations of organisms for organism taxon `nodes`.
        
        Parameters
        ----------
        nodes : List[:class:`anytree.node.node.Node`]
            List of organism taxon nodes. These nodes are **not** traversed to find child nodes!
        
        Returns
        -------
        List[str]
            List of organism abbreviations from the `nodes` passed. *None* if no :attr:`TaxonType.ORGANISM` node was passed.
        """
        if nodes is None:
            return None
        
        abbreviations = []
        for node in nodes:
            
            if node.type == TaxonType.ORGANISM:
                abbreviations.append(node.abbreviation)
        
        if len(abbreviations) == 0:
            abbreviations = None
        
        return abbreviations
    
    def getOrganismAbbreviationsByPath(self, path: 'Gammaproteobacteria/Enterobacterales', exceptPaths: List['Gammaproteobacteria/unclassified'] = None, oneOrganismPerSpecies = settings.defaultOneOrganismPerSpecies) -> List[str]:
        """
        Get abbreviations of organisms by a part of their `path`.
        
        Parameters
        ----------
        path : str
            Part of the path of the desired organism taxons. The parts of the path specified here have to match the wording of the path nodes exactly, i.e. 'Enterobac' will **not** match 'Enterobacterales'.
        exceptPaths : Iterable[str] or str
            Paths which match any of these will not be returned. Accepts iterables of exceptions or a single string exception.
        oneOrganismPerSpecies : bool, optional
            If *True*, return only the first organism node of each species node.
            
        Returns
        -------
        List[str]
            List of organism abbreviations from the organism taxon nodes found at the end of `path`. *None* if no `path` leading to an :attr:`TaxonType.ORGANISM` node was passed.
        """
        return self.getOrganismAbbreviations( self.getOrganismNodesByPath(path, exceptPaths, oneOrganismPerSpecies=oneOrganismPerSpecies) )
    
    def getOrganismAbbreviationsByName(self, name: 'Escherichia', oneOrganismPerSpecies = settings.defaultOneOrganismPerSpecies) -> List[str]:
        """
        Get abbreviations of organisms by a part of their `name`.
        
        Parameters
        ----------
        name : str
            Part of the name of the desired organism taxons. The name may be abbreviated, i.e. 'Escherichia' will match 'Escherichia Coli K-12 MG1655'.
        oneOrganismPerSpecies : bool, optional
            If *True*, return only the first organism node of each species node.
            
        Returns
        -------
        List[str]
            List of organism abbreviations containing `name` in their name attribute. *None* if none found.
        """
        return self.getOrganismAbbreviations( self.getOrganismNodesByName(name, oneOrganismPerSpecies=oneOrganismPerSpecies) )
    
    def searchNodesByName(self, name: 'Escherichia', taxonType: TaxonType = None) -> List[Node]:
        """
        Search taxons of a certain type by their name.
        
        Parameters
        ----------
        name : str
            Name of the taxon to be found. The name may be abbreviated, i.e. 'Escherichia' will match 'Escherichia Coli K-12 MG1655'.
        taxonType : TaxonType, optional
            Type of the taxons to be searched. Taxons of any other type are ignored. If *None*, all taxon types are searched.
        
        Returns
        -------
        List[:class:`anytree.node.node.Node`]
            All Nodes containing `name` in their name attribute. *None* if none can be found.
            Only taxons of :class:`TaxonType` `taxonType` are returned. If *None*, all taxon types are considered.
        """
        resultsTuple = anytree.search.findall(self.tree, filter_ = lambda node: (taxonType is node.type or taxonType is None) and name in node.name)
        if len( resultsTuple ) == 0:
            return None
        else:
            return list(resultsTuple)
        
    def searchNodesByPath(self, path: 'Gammaproteobacteria/Enterobacterales', taxonType: TaxonType = None, exceptPaths: 'list of "Gammaproteobacteria/unclassified Bacteria" etc.' = None) -> List[Node]:
        """
        Search taxons of a certain type by their path, allowing exceptions.
        
        Parameters
        ----------
        path : str
            Part of the path of the desired organism taxons. The parts of the path specified here have to match the wording of the path nodes exactly, i.e. 'Enterobac' will **not** match 'Enterobacterales'.
        taxonType : TaxonType, optional
            Type of the taxons to be searched. Taxons of any other type are ignored. If *None*, all taxon types are searched.
        exceptPaths : Iterable[str] or str
            Paths which match any of these will not be returned. Accepts iterables of exceptions or a single string exception.
        
        Returns
        -------
        List[:class:`anytree.node.node.Node`]
            All nodes containing `path` in their path. *None* if none can be found.
            Each path element has to be delimited by a slash ('/').
            Each path element has to match the name of the intermediate taxon exactly, i.e. 'Enterobac' will **not** match 'Enterobacterales'.
        """
        pathElements = path.split('/')
        
        if path.startswith("/"):
            pathElements[0] = "root"
        
        lastPathElement = pathElements.pop()
        pathElements.reverse()
        
        nodesFound = anytree.search.findall(self.tree, filter_ = lambda node: node.name == lastPathElement)
        
        matchingNodes = []
        
        # find nodes down to the supplied path
        for node in nodesFound:
            
            parent = node.parent
            nodeMatches = True
            
            for lastPathElement in pathElements:
            
                if parent.name != lastPathElement:
                    nodeMatches = False
                    break
                else:
                    parent = parent.parent
            
            if nodeMatches is True:
                matchingNodes.append(node)
        
        # find all children of surviving nodes, filter by TaxonType
        validNodes = []
        for node in matchingNodes:
            descendants = node.descendants
            
            if taxonType is not None:
                for descendant in descendants:
                    if descendant.type == taxonType:
                        validNodes.append(descendant)
            else:
                validNodes = matchingNodes
                validNodes.extend(descendants)
        
        if exceptPaths is None:
            return validNodes
        else:
            if not isinstance(exceptPaths, Iterable) or isinstance(exceptPaths, str):
                exceptPaths = [exceptPaths]
                
            # filter excepted paths
            filteredNodes = []
            for node in validNodes:
                filterNode = False
                for exception in exceptPaths:
                    nodePath = Taxonomy.nodePath2String(node)
                    if exception in nodePath:
                        filterNode = True
                        break
                if filterNode is False:
                    filteredNodes.append(node)
             
            return filteredNodes
    
    @staticmethod
    def nodePath2String(node: Node) -> str:
        """
        Parameters
        ----------
        node : :class:`anytree.node.node.Node`
            Node which' path to be expressed as a string.
        
        Returns
        -------
        str
            Full path of `node`, expressed as string. Each taxon level is delimited by a slash ('/').
        """
        return '/'.join([''] + [str(x.name) for x in node.path])
    
    def _parse(self, raw, isNCBI) -> Node:
        
        speciesRegex = re.compile(' \[TAX:[\d]+\]$')
        organismRegex = re.compile('^([a-z]{3,4})  ')
        
        root = Node('root', type = TaxonType.ROOT)
        
        previousNode = root
        previousLevel = 0
        
        for line in raw:
            
            # filter empty line
            if len(line) == 0:
                continue
            
            levelCharacter = line[0]
            
            # filter comments etc. Eveything but line starting with an uppercase letter.
            if not levelCharacter.isupper():
                continue
            
            levelNumber = ord(levelCharacter) - 64
            
            # check if the level is [A-Z]
            if levelNumber < 1 or levelNumber > 26:
                continue
            
            entry = line[1:].strip()
            
            
            # has level changed?
            levelChange = levelNumber - previousLevel
            
            
            # find parent Node
            if levelChange == 0: # same level
                parent = previousNode.parent 
            elif levelChange < 0: # went up in tree
                parent = previousNode.parent
                for _ in range(-levelChange): # for each level we jumped up the tree
                    parent = parent.parent # trace back parents
            else: # went down in tree
                parent = previousNode
            
            
            # is this a species?
            isSpecies = False
            if parent.type is TaxonType.OTHER:
                speciesSplit = speciesRegex.split(entry)
                
                if len(speciesSplit) > 1: # is species
                    isSpecies = True
                    species = speciesSplit[0]
                
                
            # is this an organism?
            isOrganism = False
            if isNCBI is False or parent.type == TaxonType.SPECIES:
                organismSplit = organismRegex.split(entry)
                
                if len(organismSplit) > 1: # is organism
                    isOrganism = True
                    abbreviation = organismSplit[1]
                    name = organismSplit[2]
                elif isNCBI is True: # for NCBI only, this is an incomplete organism (no abbreviation, only Taxon number)
                    continue
                
            
            
            # create new Node
            if isOrganism is True:
                newNode = Node(name, parent, type = TaxonType.ORGANISM, abbreviation = abbreviation)
                self.indexOnAbbreviation[abbreviation] = newNode
            elif isSpecies is True:
                newNode = Node(species, parent, type = TaxonType.SPECIES)
            else:
                newNode = Node(entry, parent, type = TaxonType.OTHER)
            
            # save variables for next round
            previousNode = newNode
            previousLevel = levelNumber
            
        return root



class NCBI(Taxonomy):
    """
    The taxonomy of organisms in KEGG, following the NCBI scheme: `<http://www.kegg.jp/kegg-bin/get_htext?br08610.keg>`_
    """
    
    @staticmethod
    @cache(folder_path = 'taxonomy/', file_name = 'NCBI_parsed')
    def getTaxonomy() -> Taxonomy:
        """
        Downloads and parses raw taxonomy from KEGG into an anytree object.
        
        Returns
        -------
        Taxonomy
            NCBI taxonomy.
        
        Raises
        ------
        URLError
            If connection to KEGG fails.
        """
        raw = Database.getTaxonomyNCBI()
        return Taxonomy(raw, isNCBI = True)


class KEGG(Taxonomy):
    """
    The taxonomy of organisms in KEGG, following KEGG's own scheme: `<http://www.kegg.jp/kegg-bin/get_htext?br08601.keg>`_
    """
    
    @staticmethod
    @cache(folder_path = 'taxonomy/', file_name = 'KEGG_parsed')
    def getTaxonomy() -> Taxonomy:
        """
        Downloads and parses raw taxonomy from KEGG into an anytree object.
        
        Returns
        -------
        Taxonomy
            KEGG taxonomy.
        
        Raises
        ------
        URLError
            If connection to KEGG fails.
        """
        raw = Database.getTaxonomyKEGG()
        return Taxonomy(raw, isNCBI = False)
