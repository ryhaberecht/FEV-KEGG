"""
In this module robustness is reduced to the topological point of view. When doing so robustness equals redundancy.
The term 'robustness' is treated as a special case of 'flexibility', meaning both are concepts of redundancy.

Definitions
-----------
If a key element is deleted from a graph, e.g. all enzymes realising a certain EC number are deleted from an organism's genome:

'Flexibility' exists if the former substrates and/or products can still be (at least partially) metabolised via alternative paths, in their original respective direction.
This does **not** mean the alternative paths have to include both, the orginal substrate and product, at the same time.

'Robustness' exists if the former substrates and products are (at least partially) still connected, in their original direction, albeit via alternative paths.
This means the alternative paths have to include both, the orginal substrate and product, at the same time.
"""

from FEV_KEGG.Graph.Models import DirectedMultiGraph, Path, MarkedPath
import FEV_KEGG.KEGG.Organism as Organism
from FEV_KEGG.Evolution.Clade import Clade, CladePair
from typing import Dict, Set, Tuple
from FEV_KEGG.Graph.Elements import Element
from enum import Enum
from FEV_KEGG import settings
import tqdm
from FEV_KEGG.Util.Util import updateDictUpdatingValue


class Robustness():
    
    def __init__(self, graph: DirectedMultiGraph, onlyLargestComponent = False):
        """
        Robustness metrics for a `graph`.
        
        If a key element is deleted from a `graph`, e.g. all enzymes (edges from substrate [edge source] to product [edge target]) realising a certain EC number (key) are deleted from an organism's genome:
        'Robustness' exists if the former edges' sources and targets are (at least partially) still connected, in their original direction, albeit via alternative paths.
        This means the alternative paths have to include both, the orginal source and target, at the same time.
        
        Parameters
        ----------
        graph : DirectedMultiGraph
            The graph to be measured for its robustness.
        onlyLargestComponent : bool, optional
            If *True*, reduce `graph` to its largest component before measuring robustness.
        
        Attributes
        ----------
        self.redundantPathsForEdgeForKey : Dict[Element, Dict[Tuple[Element, Element], Set[Path]]]
            Each key element in `graph` pointing to a psuedo-set (dictionary.keys()) of its edges, represented by a tuple of both participating nodes.
            Each edge points to a set of paths, which would provide redundancy for this edge if it (or the whole key element) were to be removed.
        
        self.sumKeys : int
            Sum of all individual key elements in `graph`.
        self.sumBreakingKeys : int
            Sum of keys which cause the graph to break when removed, i.e. which have not a single edge with redundant paths between the same source and target nodes.
        self.sumPartiallyRedundantKeys : int
            Sum of keys which, if removed, have redundant edges, but not all of them are redundant.
        self.sumRedundantKeys : int
            Sum of keys which, if removed, have only redundant edges.
        
        self.sumEdges = 0 : int
            Sum of all edges in `graph`.
        self.sumBreakingEdges : int
            Sum of edges which, if removed, cause the graph to break, because they are not redundant.
        self.sumRedundantEdges : int
            Sum of edges which, if removed, still have other redundant edges between the same source and target nodes.
        self.sumPaths : int
            Sum of all paths providing redundancy for edges.
            
        self.partiallyRedundantKeyPathCounts : Dict[Element, int]
            Edge key element pointing to the number of redundant paths it can be replaced with, although only partially. This is possible if the edge key occurs in multiple edges connecting more than two different nodes.
        self.redundantKeyPathCounts : Dict[Element, int]
            Edge key element pointing to the number of redundant paths it can be replaced with.
            
        self.redundantEdgePathCounts : Dict[Tuple[Element, Element, Element], int]
            Full edge tuples including nodes AND the key element, pointing to the number of redundant paths it can be replaced with.
        
        self.nonRedundantKeys : Set[Element]
            Set of key elements which cause the `graph` to break when removed, i.e. which have not a single edge with redundant paths between the same source and target nodes.
        self.partiallyRedundantKeys : Set[Element]
            Set of key elements which, if removed, have redundant edges, but not all of them are redundant.
        self.redundantKeys : Set[Element]
            Set of key elements which, if removed, have only redundant edges.
        
        self.nonRedundantEdges : Set[Tuple[Element, Element, Element]]
            Set of edge tuples including nodes AND the key element, which cause the graph to break, because they have no redundant path.
        self.redundantEdges : Set[Tuple[Element, Element, Element]]
            Set of edge tuples including nodes AND the key element, which have redundant paths.
        
        self.paths : Set[Path]
            Set of all paths which act as redundancy for edges when a key has been removed.
            
        self.partiallyRedundantKeyPaths : Dict[Element, Set[Path]]
            Edge key element pointing to the set of redundant paths it can be replaced with, although only partially. This is possible if the edge key occurs in multiple edges connecting more than two different nodes.
        self.redundantKeyPaths : Dict[Element, Set[Path]]
            Edge key element pointing to the set of redundant paths it can be replaced with.
        
        self.redundantEdgesRatio : float
            Ratio of the sum of redundant edges to the sum of all edges.
        self.nonRedundantEdgesRatio : float
            Ratio of the sum of nonRedundant edges to the sum of all edges.
        
        self.redundantKeysRatio : float
            Ratio of the sum of redundant key elements to the sum of all key elements.
        self.partiallyRedundantKeysRatio : float
            Ratio of the sum of partially redundant key elements to the sum of all key elements.
        self.nonRedundantKeysRatio : float
            Ratio of the sum of non-redundant key elements to the sum of all key elements.
        """
        self.redundantPathsForEdgeForKey = dict()
        
        if onlyLargestComponent:
            # get largest component of graph, which is a copy
            graph = graph.getLargestComponent()
        else:
            # still use a only a copy of the original graph
            graph = graph.copy()        
        
        # remove each edge key independently, but all its edges at once
        edgesForKey = graph.getEdgesForKey()
        
        iterator = edgesForKey.items()
        if settings.verbosity >= 1:
            if settings.verbosity >= 2:
                print( 'Calculating robust paths for ' + str(len(edgesForKey)) + ' edge keys...' )
            iterator = tqdm.tqdm(iterator, total = len(edgesForKey), unit = ' edge keys')
        
        for element, edges in iterator:
            
            # pre-fill key => edges => ...
            self.redundantPathsForEdgeForKey[element] = dict.fromkeys([(source, target) for source, target, _ in edges], None)
            
            # remove all edges sharing the same key
            graph.removeEdges(edges)
            
            # for each edge containing element as its key
            for edge in edges:
                source, target, _ = edge
                
                # determine new shortest paths between the old nodes. Fill key => edges => paths
                shortestPaths = graph.getShortestPaths(source, target)
                self.redundantPathsForEdgeForKey[element][(source, target)] = shortestPaths
            
            # add the edges sharing the same key again, faster than copying the whole graph
            graph.addEdges(edges)            
        
        self._calculateMetrics()
                    
    @classmethod
    def fromOrganismGroup(cls, group: Organism.Group, majorityPercentage = None):
        """
        Robustness metrics for a `group` core metabolism.
        
        Parameters
        ----------
        group : Organism.Group
            The group from which to extract the graph.
        majorityPercentage : float, optional
            If *None*, use collective EC graph. If not *None*, use majority EC graph with `majorityPercentage`% majority.
        
        Returns
        -------
        Robustness
        """
        if majorityPercentage is None:
            return cls(group.collectiveEcGraph())
        else:
            return cls(group.majorityEcGraph(majorityPercentage))
        
    @classmethod
    def fromClade(cls, clade: Clade, majorityPercentage = None):
        """
        Robustness metrics for a `clade` core metabolism.
        
        Parameters
        ----------
        clade : Clade
            The clade from which to extract the graph.
        majorityPercentage : float, optional
            If *None*, use collective EC graph. If not *None*, use majority EC graph with `majorityPercentage`% majority.
        
        Returns
        -------
        Robustness
        """
        return cls.fromOrganismGroup(clade.group, majorityPercentage)
    
    def _calculateMetrics(self):
        
        # basic metrics
        ## sums
        ### keys
        self.sumKeys = 0
        self.sumBreakingKeys = 0
        self.sumPartiallyRedundantKeys = 0
        self.sumRedundantKeys = 0
        
        ### edges
        self.sumEdges = 0
        self.sumBreakingEdges = 0
        self.sumRedundantEdges = 0
        
        
        ## counts
        ### key -> paths
        self.partiallyRedundantKeyPathCounts = dict() # Dict[Element, int] edge key pointing to the number of redundant paths it can be replaced with, although only partially. This is possible if the edge key occurs in multiple edges connecting more than two different nodes.
        self.redundantKeyPathCounts = dict() # Dict[Element, int] edge key pointing to the number of redundant paths it can be replaced with.
        
        ### key -> edges
        
        ### edge -> paths
        self.redundantEdgePathCounts = dict() # Dict[Tuple[Element, Element, Element], int] full edge tuples including nodes AND the key element, pointing to the number of redundant paths it can be replaced with.
        
        
        ## listings
        ### keys
        self.nonRedundantKeys = set()
        #property self.partiallyRedundantKeys = set()
        #property self.redundantKeys = set()
        
        ### edges
        self.nonRedundantEdges = set() # Set[Tuple[Element, Element, Element]] full edge tuples!
        #property self.redundantEdges = set()
        
        ### paths
        self.paths = set()
        
        
        ## assignments
        ### key -> paths
        #property self.partiallyRedundantKeyPaths = dict()
        #property self.redundantKeyPaths = dict()
        
        
        # calculation
        iterator = self.redundantPathsForEdgeForKey.items()
        if settings.verbosity >= 1:
            if settings.verbosity >= 2:
                print( 'Gathering robust path statistics for ' + str(len(self.redundantPathsForEdgeForKey)) + ' edge keys...' )
            iterator = tqdm.tqdm(iterator, total = len(self.redundantPathsForEdgeForKey), unit = ' edge keys')
        
        for key, edgesDict in iterator:
            
            self.sumKeys += 1
            
            tmpRedundantEdges = 0
            tmpRedundantPaths = 0
            
            for edge, paths in edgesDict.items():
                
                source, target = edge
                self.sumEdges += 1
                
                pathsLength = len(paths)
                
                if pathsLength > 0: # edge has redundant paths
                    self.sumRedundantEdges += 1
                    self.redundantEdgePathCounts[(source, target, key)] = pathsLength
                    tmpRedundantEdges += 1
                    
                else: # edge has no redundant paths
                    self.sumBreakingEdges += 1
                    self.nonRedundantEdges.add((source, target, key))
                
                tmpRedundantPaths += pathsLength
                self.paths.update(paths)
            
            if tmpRedundantEdges == len(edgesDict.keys()): # all edges have redundant paths
                self.sumRedundantKeys += 1
                self.redundantKeyPathCounts[key] = tmpRedundantPaths
            
            elif tmpRedundantEdges > 0: # only part of the edges have redundant paths
                self.sumPartiallyRedundantKeys += 1
                self.partiallyRedundantKeyPathCounts[key] = tmpRedundantPaths
            
            else: # none of the edges have redundant paths
                self.sumBreakingKeys += 1
                self.nonRedundantKeys.add(key)
                
        
        # derived metrics
        ## sums
        self.sumPaths = len(self.paths) # sum of all redundant paths. Counting each path only once, because there are likely duplicates.
        
        ## ratios
        ### edges
        self.redundantEdgesRatio = self.sumRedundantEdges/self.sumEdges
        self.nonRedundantEdgesRatio = self.sumBreakingEdges/self.sumEdges
        
        ### keys
        self.redundantKeysRatio = self.sumRedundantKeys/self.sumKeys
        self.partiallyRedundantKeysRatio = self.sumPartiallyRedundantKeys/self.sumKeys
        self.nonRedundantKeysRatio = self.sumBreakingKeys/self.sumKeys 
    
    @property
    def partiallyRedundantKeys(self) -> Set[Element]:
        return set(self.partiallyRedundantKeyPathCounts.keys())
    
    @property
    def redundantKeys(self) -> Set[Element]:
        return set(self.redundantKeyPathCounts.keys())
    
    @property
    def redundantEdges(self) -> Set[Element]:
        return set(self.redundantEdgesCounts.keys())
    
    @property
    def partiallyRedundantKeyPaths(self) -> Dict[Element, Set[Path]]:
        return self._keyPaths(self.partiallyRedundantKeys)
    
    @property
    def redundantKeyPaths(self) -> Dict[Element, Set[Path]]:
        return self._keyPaths(self.redundantKeys)
        
    def _keyPaths(self, keys: Set[Element]) -> Dict[Element, Set[Path]]:
        pathsDict = dict()

        for key in keys:
            
            allPaths = pathsDict.get(key)
            if allPaths is None:
                allPaths = set()
                pathsDict[key] = allPaths
                
            edgesDict = self.redundantPathsForEdgeForKey[key]
            for paths in edgesDict.values():
                if len(paths) > 0:
                    allPaths.update(paths)
        
        return pathsDict

    
class RobustnessContribution():
    
    def __init__(self, robustness: Robustness, specialKeys: Set[Element]):
        """
        Contribution to robustness accountable to edges with `specialKeys`.
        
        Allows to answer the question how much certain key elements contribute to the robustness of a graph.
        
        Parameters
        ----------
        robustness : Robustness
        specialKeys : Set[Element]
            Set of key elements viewed to be somehow special. One type of special could be 'neofunctionalised'.
        
        Attributes
        ----------
        self.robustness : Robustness
        
        self.sumSpecialKeys : int
            Sum of special keys passed.
        self.sumRedundantKeysWithSpecialKeyOnPaths : int
            Sum of redundant keys which have a special key on an alternative path.
        self.sumPartiallyRedundantKeysWithSpecialKeyOnPaths : int
            Sum of partially redundant keys which have a special key on an alternative path.
        self.sumPathsWithSpecialKeys : int
            Sum of alternative paths with a special key on it.
        
        self.pathsWithSpecialKeys : Set[MarkedPath]
            Set of paths with special keys on them.
        
        self.redundantKeySpecialKeysOnPaths : Dict[Element, Set[Element]]
            Sets of special keys which are on an alternative path of a redundant key, keyed by the key.
        self.partiallyRedundantKeySpecialKeysOnPaths : Dict[Element, Set[Element]
            Sets of special keys which are on an alternative path of a partially redundant key, keyed by the key.
        
        self.specialKeyOnRedundantKeysPaths : Dict[Element, Set[Element]]
            Sets of redundant keys which have an alternative path with a special key on it, keyed by the special key.
        self.specialKeyOnPartiallyRedundantKeysPaths : Dict[Element, Set[Element]]
            Sets of partially redundant keys which have an alternative path with a special key on it, keyed by the special key.
        
        self.pathsWithSpecialKeyRatio : float
            Ratio of the sum of paths with a special key on it to the sum of all paths.
        self.redundantKeysWithSpecialKeyOnPathsRatio : float
            Ratio of the sum of redundant keys with a special key on its paths to the sum of all redundant keys.
        self.partiallyRedundantKeysWithSpecialKeyOnPathsRatio : float
            Ratio of the sum of partially redundant keys with a special key on its paths to the sum of all partially redundant keys.
        
        self.redundantKeyPathsWithSpecialKey : Dict[Element, Set[MarkedPath]]
            Redundant key element pointing to the set of marked redundant paths it can be replaced with, which contain a special key.
        self.partiallyRedundantKeyPathsWithSpecialKey : Dict[Element, Set[MarkedPath]]
            Partially redundant key element pointing to the set of marked redundant paths it can be replaced with, which contain a special key.
        """
        self.robustness = robustness
        
        # basic metrics
        ## sums
        ### special keys
        self.sumSpecialKeys = len(specialKeys)
        #self.sumSpecialKeysOnPaths = 0
        
        ### keys
        self.sumRedundantKeysWithSpecialKeyOnPaths = 0
        self.sumPartiallyRedundantKeysWithSpecialKeyOnPaths = 0
        
        ### edges
        #self.sumEdgesWithSpecialKeysOnPaths = 0
        
        
        ## counts
        ### special key -> paths
        #self.specialKeyOnPathsCounts = dict() # on how many redundant paths is a certain special key?
        
        ### special key -> edges
        #self.specialKeyOnEdgePathsCounts = dict() # for how many edge is a certain special key on a redundant path?
        
        ### special key -> keys
        #self.specialKeyOnRedundantKeyPathsCounts = dict() # for how many redundant keys is a certain special key on a redundant path?
        #self.specialKeyOnPartiallyRedundantKeyPathsCounts = dict() # for how many partially redundant keys is a certain special key on a redundant path?
        
        ### key -> special keys
        #self.redundantKeyWithSpecialKeyOnPathsCounts = dict()
        #self.partiallyRedundantKeyWithSpecialKeyOnPathsCounts = dict()
        
        ### edge -> special keys
        #self.edgeWithSpecialKeyOnPathsCounts = dict()
        
        ### path -> special keys
        #self.pathWithSpecialKeyCounts = dict()
        
        
        ## listings
        ### special keys
        #self.specialKeysOnPaths = set() # which are actually in the graph
        #self.specialKeysOnRedundantKeyPaths = set()
        #self.specialKeysOnPartiallyRedundantKeyPaths = set()
        
        ### keys
        #self.redundantKeysWithSpecialKeyOnPaths = set()
        #self.partiallyRedundantKeysWithSpecialKeyOnPaths = set()
        
        ### edges
        #self.edgesWithSpecialKeyOnPaths = set()
        
        ### paths
        self.pathsWithSpecialKeys = set()
        
        
        ## assignments
        ### key -> special keys
        self.redundantKeySpecialKeysOnPaths = dict()
        self.partiallyRedundantKeySpecialKeysOnPaths = dict()
        
        ### key -> paths
        self.redundantKeyPathsWithSpecialKey = dict()
        self.partiallyRedundantKeyPathsWithSpecialKey = dict()


        # calculation
        
        iterator = self.robustness.redundantPathsForEdgeForKey.items()
        if settings.verbosity >= 1:
            if settings.verbosity >= 2:
                print( 'Calculating robust path contributions for ' + str(len(self.robustness.redundantPathsForEdgeForKey)) + ' edge keys...' )
            iterator = tqdm.tqdm(iterator, total = len(self.robustness.redundantPathsForEdgeForKey), unit = ' edge keys')
        
        for key, edgesDict in iterator:
            
            hasKeySpecialKeyOnPath = False
            tmpRedundantEdges = 0
            tmpSpecialKeysOnPaths = set()
            tmpMarkedPaths = set()
            
            for _, paths in edgesDict.items():
                
                #source, target = edge               
                pathsLength = len(paths)
                
                if pathsLength > 0: # edge has redundant paths
                    
                    tmpRedundantEdges += 1
                    
                    for path in paths:
                        
                        markedPath = MarkedPath(path, specialKeys)
                        specialKeysOnPath = markedPath.specialKeys
                        specialKeysOnPathLength = len(specialKeysOnPath)
                        
                        if specialKeysOnPathLength > 0: # path has special keys
                            
                            tmpSpecialKeysOnPaths.update(specialKeysOnPath)
                            self.pathsWithSpecialKeys.add(markedPath)
                            tmpMarkedPaths.add(markedPath)
                            hasKeySpecialKeyOnPath = True
            
            if hasKeySpecialKeyOnPath is True:
                
                if tmpRedundantEdges == len(edgesDict.keys()): # all edges have redundant paths
                    self.sumRedundantKeysWithSpecialKeyOnPaths += 1
                    self.redundantKeySpecialKeysOnPaths[key] = tmpSpecialKeysOnPaths
                    self.redundantKeyPathsWithSpecialKey[key] = tmpMarkedPaths
                
                elif tmpRedundantEdges > 0: # only part of the edges have redundant paths
                    self.sumPartiallyRedundantKeysWithSpecialKeyOnPaths += 1
                    self.partiallyRedundantKeySpecialKeysOnPaths[key] = tmpSpecialKeysOnPaths
                    self.partiallyRedundantKeyPathsWithSpecialKey[key] = tmpMarkedPaths
                
                else: # none of the edges have redundant paths
                    pass
        
        
        
        # derived metrics
        ## sum
        self.sumPathsWithSpecialKeys = len(self.pathsWithSpecialKeys) # sum of all paths with at least one special key on it. Counting each path only once, because there are likely dupliactes.
        
        ## ratios
        self.pathsWithSpecialKeyRatio = 0 if self.robustness.sumPaths == 0 else self.sumPathsWithSpecialKeys/self.robustness.sumPaths
        self.redundantKeysWithSpecialKeyOnPathsRatio = 0 if self.robustness.sumRedundantKeys == 0 else self.sumRedundantKeysWithSpecialKeyOnPaths/self.robustness.sumRedundantKeys
        self.partiallyRedundantKeysWithSpecialKeyOnPathsRatio = 0 if self.robustness.sumPartiallyRedundantKeys == 0 else self.sumPartiallyRedundantKeysWithSpecialKeyOnPaths/self.robustness.sumPartiallyRedundantKeys
        
        ## assignments
        ### special key -> keys
        self.specialKeyOnRedundantKeysPaths = dict()
        self.specialKeyOnPartiallyRedundantKeysPaths = dict()
        
        for key, specialKeys in self.redundantKeySpecialKeysOnPaths.items():
            for specialKey in specialKeys:
                currentSet = self.specialKeyOnRedundantKeysPaths.get(specialKey)
                if currentSet is None:
                    currentSet = set()
                    self.specialKeyOnRedundantKeysPaths[specialKey] = currentSet
                currentSet.add(key)
        
        for key, specialKeys in self.partiallyRedundantKeySpecialKeysOnPaths.items():
            for specialKey in specialKeys:
                currentSet = self.specialKeyOnPartiallyRedundantKeysPaths.get(specialKey)
                if currentSet is None:
                    currentSet = set()
                    self.specialKeyOnPartiallyRedundantKeysPaths[specialKey] = currentSet
                currentSet.add(key)
        
 
        
                    



class Flexibility():
    #@profile
    def __init__(self, graph: DirectedMultiGraph, onlyLargestComponent = False):
        """
        Flexibility metrics for a `graph`.
        
        Parameters
        ----------
        graph : DirectedMultiGraph
            The graph to be measured for its flexibility.
        onlyLargestComponent : bool, optional
            If *True*, reduce `graph` to its largest component before measuring robustness.
        
        Attributes
        ----------
        self.redundantPathsTupleForEdgeForKey : Dict[Element, Dict[Tuple[Element, Element], Tuple[Set[Path], Set[Path]]]]
            Each key element in `graph` pointing to a pseudo-set (dictionary.keys()) of its edges, represented by a tuple of both participating nodes.
            Each edge points to a set of paths, which would provide redundancy for this edge if it (or the whole key element) were to be removed.
        
        self.sumKeys : int
            Sum of all individual key elements in `graph`.
        self.sumEdges = 0 : int
            Sum of all edges in `graph`.
        self.sumPaths : int
            Sum of all paths providing redundancy for edges.
        
        
        self.sumNonRedundantKeys : int
            Sum of keys which, if removed, have not a single redundant edge. An edge is only redundant if it has redundant paths for both its source and its target node.
        self.sumPartiallyRedundantKeys : int
            Sum of keys which, if removed, have redundant edges, but not all of them are redundant. An edge is only redundant if it has redundant paths for both its source and its target node.
        self.sumRedundantKeys : int
            Sum of keys which, if removed, have only redundant edges. An edge is only redundant if it has redundant paths for both its source and its target node.
        
        self.sumNonTargetRedundantKeys : int
            Sum of keys which, if removed, have not a single target-redundant edge. An edge is target-redundant if it has redundant paths for its target node.
        self.sumPartiallyTargetRedundantKeys : int
            Sum of keys which, if removed, have target-redundant edges, but not all of them are target-redundant. An edge is target-redundant if it has redundant paths for its target node.
        self.sumTargetRedundantKeys : int
            Sum of keys which, if removed, have only target-redundant edges. An edge is target-redundant if it has redundant paths for its target node.
        
        self.sumNonSourceRedundantKeys : int
            Sum of keys which, if removed, have not a single source-redundant edge. An edge is source-redundant if it has redundant paths for its source node.
        self.sumPartiallySourceRedundantKeys : int
            Sum of keys which, if removed, have source-redundant edges, but not all of them are source-redundant. An edge is source-redundant if it has redundant paths for its source node.
        self.sumSourceRedundantKeys : int
            Sum of keys which, if removed, have only source-redundant edges. An edge is source-redundant if it has redundant paths for its source node.
        
        
        self.nonRedundantKeys : Set[Element]
            Set of key elements which, if removed, have not a single redundant edge. An edge is only redundant if it has redundant paths for both its source and its target node.
        self.partiallyRedundantKeys : Set[Element]
            Set of key elements which, if removed, have redundant edges, but not all of them are redundant. An edge is only redundant if it has redundant paths for both its source and its target node.
        self.redundantKeys : Set[Element]
            Set of key elements which, if removed, have only redundant edges. An edge is only redundant if it has redundant paths for both its source and its target node.
        
        self.nonTargetRedundantKeys : Set[Element]
            Set of key elements which, if removed, have not a single target-redundant edge. An edge is target-redundant if it has redundant paths for its target node.
        self.partiallyTargetRedundantKeys : Set[Element]
            Set of key elements which, if removed, have target-redundant edges, but not all of them are target-redundant. An edge is target-redundant if it has redundant paths for its target node.
        self.targetRedundantKeys : Set[Element]
            Set of key elements which, if removed, have only target-redundant edges. An edge is target-redundant if it has redundant paths for its target node.
        
        self.nonSourceRedundantKeys : Set[Element]
            Set of key elements which, if removed, have not a single source-redundant edge. An edge is source-redundant if it has redundant paths for its source node.
        self.partiallySourceRedundantKeys : Set[Element]
            Set of key elements which, if removed, have source-redundant edges, but not all of them are source-redundant. An edge is source-redundant if it has redundant paths for its source node.
        self.sourceRedundantKeys : Set[Element]
            Set of key elements which, if removed, have only source-redundant edges. An edge is source-redundant if it has redundant paths for its source node.
        
        
        self.paths : Set[Path]
            Set of all paths which act as redundancy for edges when a key has been removed.
        self.targetPaths : Set[Path]
            Set of all paths which act as redundancy for the target node of edges when a key has been removed.
        self.sourcePaths : Set[Path]
            Set of all paths which act as redundancy for the source node of edges when a key has been removed.
        
        
        self.redundantKeysRatio : float
            Ratio of the sum of redundant key elements to the sum of all key elements. An edge is only redundant if it has redundant paths for both its source and its target node.
        self.partiallyRedundantKeysRatio : float
            Ratio of the sum of partially redundant key elements to the sum of all key elements. An edge is only redundant if it has redundant paths for both its source and its target node.
        self.nonRedundantKeysRatio : float
            Ratio of the sum of non-redundant key elements to the sum of all key elements. An edge is only redundant if it has redundant paths for both its source and its target node.
        
        self.targetRedundantKeysRatio : float
            Ratio of the sum of target-redundant key elements to the sum of all key elements. An edge is target-redundant if it has redundant paths for its target node.
        self.partiallyTargetRedundantKeysRatio : float
            Ratio of the sum of partially target-redundant key elements to the sum of all key elements. An edge is target-redundant if it has redundant paths for its target node.
        self.nonTargetRedundantKeysRatio : float
            Ratio of the sum of non-target-redundant key elements to the sum of all key elements. An edge is target-redundant if it has redundant paths for its target node.
        
        self.sourceRedundantKeysRatio : float
            Ratio of the sum of target-redundant key elements to the sum of all key elements. An edge is source-redundant if it has redundant paths for its source node.
        self.partiallySourceRedundantKeysRatio : float
            Ratio of the sum of partially target-redundant key elements to the sum of all key elements. An edge is source-redundant if it has redundant paths for its source node.
        self.nonSourceRedundantKeysRatio : float
            Ratio of the sum of non-target-redundant key elements to the sum of all key elements. An edge is source-redundant if it has redundant paths for its source node.
        """
        self.redundantPathsTupleForEdgeForKey = dict()
        
        if onlyLargestComponent:
            # get largest component of graph, which is a copy
            graph = graph.getLargestComponent()
        else:
            # still use a only a copy of the original graph
            graph = graph.copy()        
        
        # remove each edge key independently, but all its edges at once
        edgesForKey = graph.getEdgesForKey()
        
        iterator = edgesForKey.items()
        if settings.verbosity >= 1:
            if settings.verbosity >= 2:
                print( 'Calculating flexible paths for ' + str(len(edgesForKey)) + ' edge keys...' )
            iterator = tqdm.tqdm(iterator, total = len(edgesForKey), unit = ' edge keys')
        
        for element, edges in iterator:
            
            self.redundantPathsTupleForEdgeForKey[element] = dict()
            
            lookasideBuffer = dict()
            
            # remove all edges sharing the same key
            graph.removeEdges(edges)
            
            # for each edge containing element as its key
            for edge in edges:
                source, target, _ = edge
                
                # does the source still have other edges going away from it?
                outgoingPaths = lookasideBuffer.get((source, None)) # try the buffer first
                if outgoingPaths is None: # not in buffer
                    outgoingPaths = graph.getShortestPaths(source, None) # calculate
                    lookasideBuffer[(source, None)] = outgoingPaths
                
                # does the target still have other edges leading towards it?
                incomingPaths = lookasideBuffer.get((None, target)) # try the buffer first
                if incomingPaths is None: # not in buffer
                    incomingPaths = graph.getShortestPaths(None, target) # calculate
                    lookasideBuffer[(None, target)] = incomingPaths
                
                self.redundantPathsTupleForEdgeForKey[element][(source, target)] = (outgoingPaths, incomingPaths) # save both in result
            
            # add the edges sharing the same key again, faster than copying the whole graph
            graph.addEdges(edges)
        
        self._calculateMetrics()
    
    #@profile
    def _calculateMetrics(self):
        
        # basic metrics
        ## sums
        ### keys
        self.sumKeys = 0 # sum of all keys. Each key can be represented by multiple edges with different or overlapping pairs of source+target.
        
        #### both redundancies
        self.sumNonRedundantKeys = 0 # sum of keys for which none of their edges are redundant.
        self.sumPartiallyRedundantKeys = 0 # sum of keys for which only part of their edges are redundant. This counts all possible forms of non-redundancy, including non-source-redundancy where the edge has not other edge sharing the same source.        
        self.sumRedundantKeys = 0 # sum of keys for which all their edges are redundant, for both source and target.
        
        #### only target-redundancy
        self.sumNonTargetRedundantKeys = 0 # sum of keys for which none of their edges are target-redundant. In contrast to self.sumNonRedundantKeys, this ignores non-source-redundancy, counting only missing target-redundancy as non-redundancy.
        self.sumPartiallyTargetRedundantKeys = 0 # sum of keys for which only part of their edges are target-redundant. In contrast to self.sumPartiallyRedundantKeys, this ignores non-source-redundancy, counting only missing target-redundancy as non-redundancy.
        self.sumTargetRedundantKeys = 0 # sum of keys for which all their edges are target-redundant. In contrast to self.sumRedundantKeys, this ignores non-source-redundancy, counting only missing target-redundancy as non-redundancy.
        
        #### only source-redundancy
        self.sumNonSourceRedundantKeys = 0 # sum of keys for which none of their edges are source-redundant. In contrast to self.sumNonRedundantKeys, this ignores non-target-redundancy, counting only missing source-redundancy as non-redundancy.
        self.sumPartiallySourceRedundantKeys = 0 # sum of keys for which only part of their edges are source-redundant. In contrast to self.sumPartiallyRedundantKeys, this ignores non-target-redundancy, counting only missing source-redundancy as non-redundancy.
        self.sumSourceRedundantKeys = 0 # sum of keys for which all their edges are source-redundant. In contrast to self.sumRedundantKeys, this ignores non-target-redundancy, counting only missing source-redundancy as non-redundancy.
        
        
        ### edges
        self.sumEdges = 0 # sum of all edges
        
        #### both redundancies
        #self.sumNonRedundantEdges = 0 # sum of edges with no redundant edge sharing neither source nor target (not necessarily in the same edge).
        #self.sumPartiallyRedundantEdges = 0 # sum of edges with redundant edges sharing either soure or target.
        #self.sumRedundantEdges = 0 # sum of edges with redundant edges sharing both its source and target (not necessarily in the same edge).
        
        #### only target-redundancy
        #self.sumTargetRedundantEdges = 0 # sum of edges with redundant edges sharing its target.
        #self.sumNonTargetRedundantEdges = 0 # sum of edges with no redundant edge sharing its target.
        
        #### only source-redundancy
        #self.sumSourceRedundantEdges = 0 # sum of edges with redundant edges sharing its source.
        #self.sumNonSourceRedundantEdges = 0 # sum of edges with no redundant edge sharing its source.
        
    
        ### paths
        #self.sumSourcePaths = 0 # sum of all redundant paths coming from a source.
        #self.sumTargetPaths = 0 # sum of all redundant paths leading to a target.
        
        
        ## counts
        ### key -> paths
        #self.partiallyRedundantKeyPathCounts = dict() # Dict[Element, int] edge key pointing to the number of redundant paths it can be replaced with, although only partially. This is possible if the edge key occurs in multiple edges connectign more than two different nodes.
        #self.redundantKeyPathCounts = dict() # Dict[Element, int] edge key pointing to the number of redundant paths it can be replaced with.
        
        ### key -> edges
        
        ### edge -> paths
        #self.redundantEdgePathCounts = dict() # Dict[Tuple[Element, Element, Element], int] full edge tuples including nodes AND the key element, pointing to the number of redundant paths it can be replaced with.
        
        
        ## listings
        ### keys
        #### both redundancies
        self.nonRedundantKeys = set()
        self.partiallyRedundantKeys = set()
        self.redundantKeys = set()
        
        #### only target-redundancy
        self.nonTargetRedundantKeys = set()
        self.partiallyTargetRedundantKeys = set()
        self.targetRedundantKeys = set()
        
        #### only source-redundancy
        self.nonSourceRedundantKeys = set()
        self.partiallySourceRedundantKeys = set()
        self.sourceRedundantKeys = set()
        
        
        ### edges
        #self.nonRedundantEdges = set() # Set[Tuple[Element, Element, Element]] full edge tuples!
        #self.redundantEdges = set()
        
        ### paths
        #property self.paths = set()
        self.targetPaths = set()
        self.sourcePaths = set()
        
        
        ## assignments
        ### key -> paths
        #### both redundancies
        #property self.partiallyRedundantKeyPaths = dict()
        #property self.redundantKeyPaths = dict()
        
        #### only target-redundancy
        #property self.partiallyTargetRedundantKeyPaths = dict()
        #property self.targetRedundantKeyPaths = dict()
        
        #### only source-redundancy
        #property self.partiallySourceRedundantKeyPaths = dict()
        #property self.sourceRedundantKeyPaths = dict()

        
        # calculation
        iterator = self.redundantPathsTupleForEdgeForKey.items()
        if settings.verbosity >= 1:
            if settings.verbosity >= 2:
                print( 'Gathering flexible path statistics for ' + str(len(self.redundantPathsTupleForEdgeForKey)) + ' edge keys...' )
            iterator = tqdm.tqdm(iterator, total = len(self.redundantPathsTupleForEdgeForKey), unit = ' edge keys')
            
        for key, edgesDict in iterator: # for all keys
            
            self.sumKeys += 1
            
            tmpRedundantEdges = 0
            tmpTargetRedundantEdges = 0            
            tmpSourceRedundantEdges = 0
            
            #tmpRedundantPaths = 0
            #tmpTargetRedundantPaths = 0
            #tmpSourceRedundantPaths = 0
            
            for _, pathsTuple in edgesDict.items(): # for all edges
                
                #source, target = edge
                self.sumEdges += 1
                
                isSourceRedundant = False
                isTargetRedundant = False
                
                for index, paths in enumerate(pathsTuple): # for source/target
                    
                    pathsLength = len(paths)
                    
#                     if index == 0: # source-redundancy
#                         #tmpSumSourceWildcardEdges += 1
#                         pass
#                         
#                     elif index == 1: # target-redundancy
#                         #tmpSumTargetWildcardEdges += 1
#                         pass
#                     
#                     else:
#                         raise RuntimeError
                
                    if pathsLength > 0: # wildcardEdge has redundant paths
                        
                        if index == 0: # source-redundancy
                            #tmpSourceRedundantPaths += pathsLength
                            isSourceRedundant = True
                            self.sourcePaths.update(paths)
                            
                        elif index == 1: # target-redundancy
                            #tmpTargetRedundantPaths += pathsLength
                            isTargetRedundant = True
                            self.targetPaths.update(paths)
                        
                        else:
                            raise RuntimeError
                    
                    else: # wildcardEdge has no redundant paths
                        pass
                
                # both source and target redundant? -> edge redundant
                if isSourceRedundant and isTargetRedundant:
                    tmpRedundantEdges += 1
                
                # target redundant? -> edge target-redundant
                if isTargetRedundant:
                    tmpTargetRedundantEdges += 1
                
                # source redundant? -> edge source-redundant
                if isSourceRedundant:
                    tmpSourceRedundantEdges += 1
            
            # both redundancies
            if tmpRedundantEdges == len(edgesDict.keys()): # all edges have redundant paths
                self.sumRedundantKeys += 1
                self.redundantKeys.add(key)
            
            elif tmpRedundantEdges > 0: # only part of the edges have redundant paths
                self.sumPartiallyRedundantKeys += 1
                self.partiallyRedundantKeys.add(key)
            
            else: # none of the wildcard edges have redundant paths
                self.sumNonRedundantKeys += 1
                self.nonRedundantKeys.add(key)
            
            # target-redundancy
            if tmpTargetRedundantEdges == len(edgesDict.keys()): # all edges have target-redundant paths
                self.sumTargetRedundantKeys += 1
                self.targetRedundantKeys.add(key)
            
            elif tmpTargetRedundantEdges > 0: # only part of the edges have target-redundant paths
                self.sumPartiallyTargetRedundantKeys += 1
                self.partiallyTargetRedundantKeys.add(key)
            
            else: # none of the edges have target-redundant paths
                self.sumNonTargetRedundantKeys += 1
                self.nonTargetRedundantKeys.add(key)
            
            # source-redundancy
            if tmpSourceRedundantEdges == len(edgesDict.keys()): # all edges have source-redundant paths
                self.sumSourceRedundantKeys += 1
                self.sourceRedundantKeys.add(key)
            
            elif tmpSourceRedundantEdges > 0: # only part of the edges have source-redundant paths
                self.sumPartiallySourceRedundantKeys += 1
                self.partiallySourceRedundantKeys.add(key)
            
            else: # none of the edges have source-redundant paths
                self.sumNonSourceRedundantKeys += 1
                self.nonSourceRedundantKeys.add(key)


        
        # derived metrics
        ## sums
        #property self.sumPaths = len(self.paths) # sum of all redundant paths for both source and target. Counting each path only once, because there are likely duplicates.
        
        ## ratios
        ### keys
        #### both redundancies
        self.redundantKeysRatio = self.sumRedundantKeys/self.sumKeys
        self.partiallyRedundantKeysRatio = self.sumPartiallyRedundantKeys/self.sumKeys
        self.nonRedundantKeysRatio = self.sumNonRedundantKeys/self.sumKeys
        
        ### only target-redundancy
        self.targetRedundantKeysRatio = self.sumTargetRedundantKeys/self.sumKeys
        self.partiallyTargetRedundantKeysRatio = self.sumPartiallyTargetRedundantKeys/self.sumKeys
        self.nonTargetRedundantKeysRatio = self.sumNonTargetRedundantKeys/self.sumKeys
        
        ### only source-redundancy
        self.sourceRedundantKeysRatio = self.sumSourceRedundantKeys/self.sumKeys
        self.partiallySourceRedundantKeysRatio = self.sumPartiallySourceRedundantKeys/self.sumKeys
        self.nonSourceRedundantKeysRatio = self.sumNonSourceRedundantKeys/self.sumKeys
        
        
        ### edges
        #### both redundancies
        #self.redundantEdgesRatio = self.sumRedundantEdges/self.sumEdges
        #self.partiallyRedundantEdgesRatio = self.sumPartiallyRedundantEdges/self.sumEdges
        #self.nonRedundantEdgesRatio = self.sumNonRedundantEdges/self.sumEdges
        
        ### only target-redundancy
        
        ### only source-redundancy
    
    @property
    def paths(self):
        return self.sourcePaths.union(self.targetPaths)
    
    @property
    def sumPaths(self):
        return len(self.paths)
    
    @property
    def partiallyRedundantKeyPaths(self) -> Dict[Element, Set[Path]]:
        return self._keyPaths(self.partiallyRedundantKeys, 0)
    
    @property
    def redundantKeyPaths(self) -> Dict[Element, Set[Path]]:
        return self._keyPaths(self.redundantKeys, 0)
    
    @property
    def partiallyTargetRedundantKeyPaths(self) -> Dict[Element, Set[Path]]:
        return self._keyPaths(self.partiallyTargetRedundantKeys, 1)
    
    @property
    def targetRedundantKeyPaths(self) -> Dict[Element, Set[Path]]:
        return self._keyPaths(self.targetRedundantKeys, 1)

    @property
    def partiallySourceRedundantKeyPaths(self) -> Dict[Element, Set[Path]]:
        return self._keyPaths(self.partiallySourceRedundantKeys, -1)
    
    @property
    def sourceRedundantKeyPaths(self) -> Dict[Element, Set[Path]]:
        return self._keyPaths(self.sourceRedundantKeys, -1)
        
    def _keyPaths(self, keys: Set[Element], targetOrSource: int) -> Dict[Element, Set[Path]]:
        pathsDict = dict()

        for key in keys:
            
            allPaths = pathsDict.get(key)
            if allPaths is None:
                allPaths = set()
                pathsDict[key] = allPaths
                
            edgesDict = self.redundantPathsTupleForEdgeForKey[key]
            for pathsTuple in edgesDict.values():
                for index, paths in enumerate(pathsTuple): # for source/target
                    if len(paths) > 0:
                        
                        if index == 0 and targetOrSource <= 0: # source-redundancy
                            allPaths.update(paths)
                            
                        elif index == 1 and targetOrSource >= 0: # target-redundancy
                            allPaths.update(paths)
        
        return pathsDict
        
        
        
class FlexibilityContribution():
    
    def __init__(self, flexibility: Flexibility, specialKeys: Dict[str, Set[Element]]):
        """
        Contribution to flexibility accountable to edges with `specialKeys`.
        
        Allows to answer the question how much certain key elements contribute to the flexibility of a graph.
        
        Parameters
        ----------
        flexibility : Flexibility
        specialKeys : Set[Element]
            Set of key elements viewed to be somehow special. One type of special could be 'neofunctionalised'.
        
        Attributes
        ----------
        self.flexibility : Flexibility
        
        self.sumSpecialKeys : int
            Sum of special keys passed.
        self.sumPathsWithSpecialKeys : int
            Sum of alternative paths with a special key on it.
        
        self.sumRedundantKeysWithSpecialKeyOnPaths : int
            Sum of redundant keys which have a special key on an alternative path.
        self.sumPartiallyRedundantKeysWithSpecialKeyOnPaths : int
            Sum of partially redundant keys which have a special key on an alternative path.
        
        self.sumTargetRedundantKeysWithSpecialKeyOnPaths : int
            Sum of target-redundant keys which have a special key on an alternative path.
        self.sumPartiallyTargetRedundantKeysWithSpecialKeyOnPaths : int
            Sum of partially target-redundant keys which have a special key on an alternative path.
        
        self.sumSourceRedundantKeysWithSpecialKeyOnPaths : int
            Sum of source-redundant keys which have a special key on an alternative path.
        self.sumPartiallySourceRedundantKeysWithSpecialKeyOnPaths : int
            Sum of partially source-redundant keys which have a special key on an alternative path.
        
        
        self.pathsWithSpecialKeys : Set[MarkedPath]
            Set of redundant paths with special keys on them.
        self.targetPathsWithSpecialKeys : Set[MarkedPath]
            Set of target-redundant paths with special keys on them.
        self.sourcePathsWithSpecialKeys : Set[MarkedPath]
            Set of source-redundant paths with special keys on them.
        
        
        self.redundantKeySpecialKeysOnPaths : Dict[Element, Set[Element]]
            Sets of special keys which are on an alternative path of a redundant key, keyed by the key.
        self.partiallyRedundantKeySpecialKeysOnPaths : Dict[Element, Set[Element]
            Sets of special keys which are on an alternative path of a partially redundant key, keyed by the key.
        
        self.targetRedundantKeySpecialKeysOnPaths : Dict[Element, Set[Element]]
            Sets of special keys which are on an alternative path of a target-redundant key, keyed by the key.
        self.partiallyTargetRedundantKeySpecialKeysOnPaths : Dict[Element, Set[Element]
            Sets of special keys which are on an alternative path of a partially target-redundant key, keyed by the key.
        
        self.sourceRedundantKeySpecialKeysOnPaths : Dict[Element, Set[Element]]
            Sets of special keys which are on an alternative path of a source-redundant key, keyed by the key.
        self.partiallySourceRedundantKeySpecialKeysOnPaths : Dict[Element, Set[Element]
            Sets of special keys which are on an alternative path of a partially source-redundant key, keyed by the key.
        
        
        self.specialKeyOnRedundantKeysPaths : Dict[Element, Set[Element]]
            Sets of redundant keys which have an alternative path with a special key on it, keyed by the special key.
        self.specialKeyOnPartiallyRedundantKeysPaths : Dict[Element, Set[Element]]
            Sets of partially redundant keys which have an alternative path with a special key on it, keyed by the special key.
            
        self.specialKeyOnTargetRedundantKeysPaths : Dict[Element, Set[Element]]
            Sets of target-redundant keys which have an alternative path with a special key on it, keyed by the special key.
        self.specialKeyOnPartiallyTargetRedundantKeysPaths : Dict[Element, Set[Element]]
            Sets of partially target-redundant keys which have an alternative path with a special key on it, keyed by the special key.
        
        self.specialKeyOnSourceRedundantKeysPaths : Dict[Element, Set[Element]]
            Sets of source-redundant keys which have an alternative path with a special key on it, keyed by the special key.
        self.specialKeyOnPartiallySourceRedundantKeysPaths : Dict[Element, Set[Element]]
            Sets of partially source-redundant keys which have an alternative path with a special key on it, keyed by the special key.
        
        
        self.pathsWithSpecialKeyRatio : float
            Ratio of the sum of paths with a special key on it to the sum of all paths.
        
        
        self.redundantKeysWithSpecialKeyOnPathsRatio : float
            Ratio of the sum of redundant keys with a special key on its paths to the sum of all redundant keys.
        self.partiallyRedundantKeysWithSpecialKeyOnPathsRatio : float
            Ratio of the sum of partially redundant keys with a special key on its paths to the sum of all partially redundant keys.
        
        self.targetRedundantKeysWithSpecialKeyOnPathsRatio : float
            Ratio of the sum of target-redundant keys with a special key on its paths to the sum of all target-redundant keys.
        self.partiallyTargetRedundantKeysWithSpecialKeyOnPathsRatio : float
            Ratio of the sum of partially target-redundant keys with a special key on its paths to the sum of all partially target-redundant keys.
        
        self.sourceRedundantKeysWithSpecialKeyOnPathsRatio : float
            Ratio of the sum of source-redundant keys with a special key on its paths to the sum of all source-redundant keys.
        self.partiallySourceRedundantKeysWithSpecialKeyOnPathsRatio : float
            Ratio of the sum of partially source-redundant keys with a special key on its paths to the sum of all partially source-redundant keys.
        """
        self.flexibility = flexibility
        
        # basic metrics
        ## sums
        ### special keys
        self.sumSpecialKeys = 0
        #self.sumSpecialKeysOnPaths = 0
        
        ### keys
        #### both redundancies
        self.sumRedundantKeysWithSpecialKeyOnPaths = 0
        self.sumPartiallyRedundantKeysWithSpecialKeyOnPaths = 0
        
        #### only target-redundancy
        self.sumTargetRedundantKeysWithSpecialKeyOnPaths = 0
        self.sumPartiallyTargetRedundantKeysWithSpecialKeyOnPaths = 0
        
        #### only source-redundancy
        self.sumSourceRedundantKeysWithSpecialKeyOnPaths = 0
        self.sumPartiallySourceRedundantKeysWithSpecialKeyOnPaths = 0
        
        ### edges
        #self.sumEdgesWithSpecialKeysOnPaths = 0
        
        
        ## counts
        ### special key -> paths
        #self.specialKeyOnPathsCounts = dict() # on how many redundant paths is a certain special key?
        
        ### special key -> edges
        #self.specialKeyOnEdgePathsCounts = dict() # for how many edge is a certain special key on a redundant path?
        
        ### special key -> keys
        #self.specialKeyOnRedundantKeyPathsCounts = dict() # for how many redundant keys is a certain special key on a redundant path?
        #self.specialKeyOnPartiallyRedundantKeyPathsCounts = dict() # for how many partially redundant keys is a certain special key on a redundant path?
        
        ### key -> special keys
        #self.redundantKeyWithSpecialKeyOnPathsCounts = dict()
        #self.partiallyRedundantKeyWithSpecialKeyOnPathsCounts = dict()
        
        ### edge -> special keys
        #self.edgeWithSpecialKeyOnPathsCounts = dict()tmpRedundantEdges
        
        ### path -> special keys
        #self.pathWithSpecialKeyCounts = dict()
        
        
        ## listings
        ### special keys
        #self.specialKeysOnPaths = set() # which are actually in the graph
        #self.specialKeysOnRedundantKeyPaths = set()
        #self.specialKeysOnPartiallyRedundantKeyPaths = set()
        
        ### keys
        #self.redundantKeysWithSpecialKeyOnPaths = set()
        #self.partiallyRedundantKeysWithSpecialKeyOnPaths = set()
        
        ### edges
        #self.edgesWithSpecialKeyOnPaths = set()
        
        ### paths
        self.pathsWithSpecialKeys = set()
        self.targetPathsWithSpecialKeys = set()
        self.sourcePathsWithSpecialKeys = set()
        
        
        ## assignments
        ### key -> special keys
        #### both redundancies
        self.redundantKeySpecialKeysOnPaths = dict()
        self.partiallyRedundantKeySpecialKeysOnPaths = dict()
        
        #### only target-redundancy
        self.targetRedundantKeySpecialKeysOnPaths = dict()
        self.partiallyTargetRedundantKeySpecialKeysOnPaths = dict()
        
        #### only source-redundancy
        self.sourceRedundantKeySpecialKeysOnPaths = dict()
        self.partiallySourceRedundantKeySpecialKeysOnPaths = dict()
        
        ### key -> paths
        #### both redundancies
        self.redundantKeyPathsWithSpecialKey = dict()
        self.partiallyRedundantKeyPathsWithSpecialKey = dict()
        
        #### only target-redundancy
        self.targetRedundantKeyPathsWithSpecialKey = dict()
        self.partiallyTargetRedundantKeyPathsWithSpecialKey = dict()
        
        #### only source-redundancy
        self.sourceRedundantKeyPathsWithSpecialKey = dict()
        self.partiallySourceRedundantKeyPathsWithSpecialKey = dict()
        
        


        # calculation
        
        iterator = self.flexibility.redundantPathsTupleForEdgeForKey.items()
        if settings.verbosity >= 1:
            if settings.verbosity >= 2:
                print( 'Calculating flexible path contributions for ' + str(len(self.flexibility.redundantPathsTupleForEdgeForKey)) + ' edge keys...' )
            iterator = tqdm.tqdm(iterator, total = len(self.flexibility.redundantPathsTupleForEdgeForKey), unit = ' edge keys')
        
        for key, edgesDict in iterator:
            
            hasKeySpecialKeyOnPath = False
            tmpRedundantEdges = 0
            tmpTargetRedundantEdges = 0
            tmpSourceRedundantEdges = 0
            
            tmpSpecialKeysOnPaths = set()
            tmpSpecialKeysOnTargetPaths = set()
            tmpSpecialKeysOnSourcePaths = set()
            
            tmpMarkedPaths = set()
            tmpMarkedTargetPaths = set()
            tmpMarkedSourcePaths = set()
            
            for _, pathsTuple in edgesDict.items():
                
                isSourceRedundant = False
                isTargetRedundant = False
                
                for index, paths in enumerate(pathsTuple): # for source/target
                    
                    pathsLength = len(paths)
                    
                    if pathsLength > 0: # wildcardEdge has redundant paths
                        
                        if index == 0: # source-redundancy
                            isSourceRedundant = True
                            
                        elif index == 1: # target-redundancy
                            isTargetRedundant = True
                        
                        else:
                            raise RuntimeError
                        
                        for path in paths:
                            
                            markedPath = MarkedPath(path, specialKeys)
                            specialKeysOnPath = markedPath.specialKeys
                            specialKeysOnPathLength = len(specialKeysOnPath)
                            
                            if specialKeysOnPathLength > 0: # path has special keys
                                
                                tmpSpecialKeysOnPaths.update(specialKeysOnPath)
                                self.pathsWithSpecialKeys.add(markedPath)
                                
                                if index == 0: # source-redundancy
                                    self.sourcePathsWithSpecialKeys.add(markedPath)
                                    tmpSpecialKeysOnSourcePaths.update(specialKeysOnPath)
                                    tmpMarkedSourcePaths.add(markedPath)                                    
                                    
                                elif index == 1: # target-redundancy
                                    self.targetPathsWithSpecialKeys.add(markedPath)
                                    tmpSpecialKeysOnTargetPaths.update(specialKeysOnPath)
                                    tmpMarkedTargetPaths.add(markedPath)
                                
                                tmpMarkedPaths.add(markedPath)
                                
                                hasKeySpecialKeyOnPath = True
                    
                    else: # wildcardEdge has no redundant paths
                        pass
                                
                # both source and target redundant? -> edge redundant
                if isSourceRedundant and isTargetRedundant:
                    tmpRedundantEdges += 1
                
                # target redundant? -> edge target-redundant
                if isTargetRedundant:
                    tmpTargetRedundantEdges += 1
                
                # source redundant? -> edge source-redundant
                if isSourceRedundant:
                    tmpSourceRedundantEdges += 1
            
            if hasKeySpecialKeyOnPath is True:
                
                # both redundancies
                if tmpRedundantEdges == len(edgesDict.keys()): # all edges have redundant paths
                    self.sumRedundantKeysWithSpecialKeyOnPaths += 1
                    self.redundantKeySpecialKeysOnPaths[key] = tmpSpecialKeysOnPaths
                    self.redundantKeyPathsWithSpecialKey[key] = tmpMarkedPaths
                
                elif tmpRedundantEdges > 0: # only part of the edges have redundant paths
                    self.sumPartiallyRedundantKeysWithSpecialKeyOnPaths += 1
                    self.partiallyRedundantKeySpecialKeysOnPaths[key] = tmpSpecialKeysOnPaths
                    self.partiallyRedundantKeyPathsWithSpecialKey[key] = tmpMarkedPaths
                
                else: # none of the edges have redundant paths
                    pass
                
                # target-redundancy
                if tmpTargetRedundantEdges == len(edgesDict.keys()): # all edges have target-redundant paths
                    self.sumTargetRedundantKeysWithSpecialKeyOnPaths += 1
                    self.targetRedundantKeySpecialKeysOnPaths[key] = tmpSpecialKeysOnTargetPaths
                    self.targetRedundantKeyPathsWithSpecialKey[key] = tmpMarkedTargetPaths
                
                elif tmpTargetRedundantEdges > 0: # only part of the edges have target-redundant paths
                    self.sumPartiallyTargetRedundantKeysWithSpecialKeyOnPaths += 1
                    self.partiallyTargetRedundantKeySpecialKeysOnPaths[key] = tmpSpecialKeysOnTargetPaths
                    self.partiallyTargetRedundantKeyPathsWithSpecialKey[key] = tmpMarkedTargetPaths
                
                else: # none of the edges have target-redundant paths
                    pass
                
                # source-redundancy
                if tmpSourceRedundantEdges == len(edgesDict.keys()): # all edges have source-redundant paths
                    self.sumSourceRedundantKeysWithSpecialKeyOnPaths += 1
                    self.sourceRedundantKeySpecialKeysOnPaths[key] = tmpSpecialKeysOnSourcePaths
                    self.sourceRedundantKeyPathsWithSpecialKey[key] = tmpMarkedSourcePaths
                
                elif tmpSourceRedundantEdges > 0: # only part of the edges have source-redundant paths
                    self.sumPartiallySourceRedundantKeysWithSpecialKeyOnPaths += 1
                    self.partiallySourceRedundantKeySpecialKeysOnPaths[key] = tmpSpecialKeysOnSourcePaths
                    self.partiallySourceRedundantKeyPathsWithSpecialKey[key] = tmpMarkedSourcePaths
                
                else: # none of the edges have source-redundant paths
                    pass
        

        # derived metrics
        ## sum
        self.sumPathsWithSpecialKeys = len(self.pathsWithSpecialKeys) # sum of all paths with at least one special key on it. Counting each path only once, because there are likely dupliactes.
        
        ## ratios
        self.pathsWithSpecialKeyRatio = 0 if self.flexibility.sumPaths == 0 else self.sumPathsWithSpecialKeys/self.flexibility.sumPaths
        
        ### both redundancies
        self.redundantKeysWithSpecialKeyOnPathsRatio = 0 if self.flexibility.sumRedundantKeys == 0 else self.sumRedundantKeysWithSpecialKeyOnPaths/self.flexibility.sumRedundantKeys
        self.partiallyRedundantKeysWithSpecialKeyOnPathsRatio = 0 if self.flexibility.sumPartiallyRedundantKeys == 0 else self.sumPartiallyRedundantKeysWithSpecialKeyOnPaths/self.flexibility.sumPartiallyRedundantKeys
        
        ### only target-redundancy
        self.targetRedundantKeysWithSpecialKeyOnPathsRatio = 0 if self.flexibility.sumTargetRedundantKeys == 0 else self.sumTargetRedundantKeysWithSpecialKeyOnPaths/self.flexibility.sumTargetRedundantKeys
        self.partiallyTargetRedundantKeysWithSpecialKeyOnPathsRatio = 0 if self.flexibility.sumPartiallyTargetRedundantKeys == 0 else self.sumPartiallyTargetRedundantKeysWithSpecialKeyOnPaths/self.flexibility.sumPartiallyTargetRedundantKeys
        
        ### only source-redundancy
        self.sourceRedundantKeysWithSpecialKeyOnPathsRatio = 0 if self.flexibility.sumSourceRedundantKeys == 0 else self.sumSourceRedundantKeysWithSpecialKeyOnPaths/self.flexibility.sumSourceRedundantKeys
        self.partiallySourceRedundantKeysWithSpecialKeyOnPathsRatio = 0 if self.flexibility.sumPartiallySourceRedundantKeys == 0 else self.sumPartiallySourceRedundantKeysWithSpecialKeyOnPaths/self.flexibility.sumPartiallySourceRedundantKeys
        
        
        ## assignments
        ### special key -> keys
        #### both redundancies
        self.specialKeyOnRedundantKeysPaths = dict()
        self.specialKeyOnPartiallyRedundantKeysPaths = dict()
        
        #### only target-redundancy
        self.specialKeyOnTargetRedundantKeysPaths = dict()
        self.specialKeyOnPartiallyTargetRedundantKeysPaths = dict()
        
        #### only source-redundancy
        self.specialKeyOnSourceRedundantKeysPaths = dict()
        self.specialKeyOnPartiallySourceRedundantKeysPaths = dict()
        
        dictionaryPairs = [(self.redundantKeySpecialKeysOnPaths, self.specialKeyOnRedundantKeysPaths), 
                           (self.partiallyRedundantKeySpecialKeysOnPaths, self.specialKeyOnPartiallyRedundantKeysPaths), 
                           (self.targetRedundantKeySpecialKeysOnPaths, self.specialKeyOnTargetRedundantKeysPaths), 
                           (self.partiallyTargetRedundantKeySpecialKeysOnPaths, self.specialKeyOnPartiallyTargetRedundantKeysPaths), 
                           (self.sourceRedundantKeySpecialKeysOnPaths, self.specialKeyOnSourceRedundantKeysPaths), 
                           (self.partiallySourceRedundantKeySpecialKeysOnPaths, self.specialKeyOnPartiallySourceRedundantKeysPaths)]
        for dictA, dictB in dictionaryPairs:
            for key, specialKeys in dictA.items():
                for specialKey in specialKeys:
                    currentSet = dictB.get(specialKey)
                    if currentSet is None:
                        currentSet = set()
                        dictB[specialKey] = currentSet
                    currentSet.add(key)
        

        








class RedundancyType(Enum):
    """
    Type of a redundancy.
    
    Redundancy comes in many forms, some are defined here as constants pointing to their realising classes.
    """
    
    ROBUSTNESS = [Robustness, 0]
    """
    If a key element is deleted from a graph, e.g. all enzymes realising a certain EC number are deleted from an organism's genome:    
    Robustness exists if the former substrates and products are (at least partially) still connected, in their original direction, albeit via alternative paths.
    This means the alternative paths have to include both, the orginal substrate and product, at the same time.
    
    This means robustness is a sub-type of :attr:`FLEXIBILITY`. Only some flexible edges are also robust. Robustness and flexibility have an inheritance relation.
    """
    
    ROBUSTNESS_PARTIAL = [Robustness, 1]
    """
    This is the same as :attr:`ROBUSTNESS`, but only one of the edges of a key has to be redundant, not all its edges at once. This does *not* include the results of :attr:`ROBUSTNESS`!
    """
    
    ROBUSTNESS_BOTH = [Robustness, 2]
    """
    This combines the results of both :attr:`ROBUSTNESS` and :attr:`ROBUSTNESS_PARTIAL`.
    """
    
    
    FLEXIBILITY = [Flexibility, 0]
    """
    If a key element is deleted from a graph, e.g. all enzymes realising a certain EC number are deleted from an organism's genome:
    Flexibility exists if the former substrates and/or products can still be (at least partially) metabolised via alternative paths, in their original respective direction.
    This does **not** mean the alternative paths have to include both, the orginal substrate and product, at the same time.
    """
    
    FLEXIBILITY_PARTIAL = [Flexibility, 1]
    """
    This is the same as :attr:`FLEXIBILITY`, but only one of the edges of a key has to be redundant, not all its edges at once. This does *not* include the results of :attr:`FLEXIBILITY`!
    """
    
    FLEXIBILITY_BOTH = [Flexibility, 2]
    """
    This combines the results of both :attr:`FLEXIBILITY` and :attr:`FLEXIBILITY_PARTIAL`.
    """
    
    
    TARGET_FLEXIBILITY = [Flexibility, 3]
    """
    This is a super-type of :attr:`FLEXIBILITY`, where only the target node of each edge has to have redundant paths (leading to it), for the whole edge to be counted as redundant. It does not matter whether the source node also has redundant paths.
    Only some target-flexible edges are also flexible, in fact exactly the ones which are also source-flexible. This means that combining target-flexibility with source-flexibility yields flexibility, they have a composition relation.
    """
    
    TARGET_FLEXIBILITY_PARTIAL = [Flexibility, 4]
    """
    This is the same as :attr:`TARGET_FLEXIBILITY`, but only one of the edges of a key has to be redundant, not all its edges at once. This does *not* include the results of :attr:`TARGET_FLEXIBILITY`!
    """
    
    TARGET_FLEXIBILITY_BOTH = [Flexibility, 5]
    """
    This combines the results of both :attr:`TARGET_FLEXIBILITY` and :attr:`TARGET_FLEXIBILITY_PARTIAL`.
    """
    
    
    SOURCE_FLEXIBILITY = [Flexibility, 6]
    """
    This is a super-type of :attr:`FLEXIBILITY`, where only the source node of each edge has to have redundant paths (leaving from it), for the whole edge to be counted as redundant. It does not matter whether the target node also has redundant paths.
    Only some source-flexible edges are also flexible, in fact exactly the ones which are also target-flexible. This means that combining target-flexibility with source-flexibility yields flexibility, they have a composition relation.
    """
    
    SOURCE_FLEXIBILITY_PARTIAL = [Flexibility, 7]
    """
    This is the same as :attr:`SOURCE_FLEXIBILITY`, but only one of the edges of a key has to be redundant, not all its edges at once. This does *not* include the results of :attr:`SOURCE_FLEXIBILITY`!
    """
    
    SOURCE_FLEXIBILITY_BOTH = [Flexibility, 8]
    """
    This combines the results of both :attr:`SOURCE_FLEXIBILITY` and :attr:`SOURCE_FLEXIBILITY_PARTIAL`.
    """
    
    
    default = ROBUSTNESS
    """
    Defaults to robustness, which is (currently) the most picky measure of redundancy. If you absolutely have to over-simplify the question of redundancy, use this default redundancy type.
    Robustness was chosen as default, because it seems wise not to break the graph when dealing with incomplete datasets.
    Which can only be prevented (to the extent of our knowledge) by using the more narrow definition of 'robustness', not the broader 'flexibility',
    because the difference of 'flexibility' minus 'robustness' leaves the cases where the graph breaks, but source and/or target are still redundant.
    """
    
class Redundancy():
    
    def __init__(self, graph: DirectedMultiGraph, onlyLargestComponent = False, onlyType: RedundancyType = None):
        """
        Redundancy metrics, consisting of :class:`Flexibility` and class:`Robustness`.
        
        The most important metrics are realised as methods.
        
        Parameters
        ----------
        graph : DirectedMultiGraph
            The graph to calculate flexibility and robustness metrics for.
        onlyLargestComponent : bool, optional
            If *True*, reduce `graph` to its largest component before measuring redundancy.
        onlyType : RedundancyType, optional
            If give, only metrics for this type of redundancy are actually calculated. Requests for metrics of another type of redundancy will raise an error!
        
        Attributes
        ----------
        self.flexibility : Flexibility
        self.robustness : Robustness
        
        Warnings
        --------
        The underlying algorithms have a rather high memory-complexity. This is fine for small graphs, i.e. substance-EC graphs of the core metabolism.
        But for bigger graphs, i.e. substance-enzyme graphs of the core metabolism, memory consumption can easily exceed 16 GiB. Be sure to have swap space available!
        """
        self.onlyType = onlyType
        
        if onlyType is None or onlyType.value[0] is Flexibility:
            self.flexibility = Flexibility(graph, onlyLargestComponent)
        
        if onlyType is None or onlyType.value[0] is Robustness:
            self.robustness = Robustness(graph, onlyLargestComponent)
    
    def getRedundancyRatio(self, redundancyType: RedundancyType = RedundancyType.default) -> float:
        """
        Get ratio of redundant keys to all keys.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation.
        
        Returns
        -------
        float
        
        Raises
        ------
        ValueError
            If `onlyType` was given in the contructor, but metrics of another type of redundancy are to be returned here.
        """
        try:
            if redundancyType is RedundancyType.ROBUSTNESS:
                return self.robustness.redundantKeysRatio
            
            elif redundancyType is RedundancyType.ROBUSTNESS_PARTIAL:
                return self.robustness.partiallyRedundantKeysRatio
            
            elif redundancyType is RedundancyType.ROBUSTNESS_BOTH:
                return self.robustness.redundantKeysRatio + self.robustness.partiallyRedundantKeysRatio
            
            
            elif redundancyType is RedundancyType.FLEXIBILITY:
                return self.flexibility.redundantKeysRatio
            
            elif redundancyType is RedundancyType.FLEXIBILITY_PARTIAL:
                return self.flexibility.partiallyRedundantKeysRatio
            
            elif redundancyType is RedundancyType.FLEXIBILITY_BOTH:
                return self.flexibility.redundantKeysRatio + self.flexibility.partiallyRedundantKeysRatio
            
            
            elif redundancyType is RedundancyType.TARGET_FLEXIBILITY:
                return self.flexibility.targetRedundantKeysRatio
            
            elif redundancyType is RedundancyType.TARGET_FLEXIBILITY_PARTIAL:
                return self.flexibility.partiallyTargetRedundantKeysRatio
            
            elif redundancyType is RedundancyType.TARGET_FLEXIBILITY_BOTH:
                return self.flexibility.targetRedundantKeysRatio + self.flexibility.partiallyTargetRedundantKeysRatio
            
            
            elif redundancyType is RedundancyType.SOURCE_FLEXIBILITY:
                return self.flexibility.sourceRedundantKeysRatio
            
            elif redundancyType is RedundancyType.SOURCE_FLEXIBILITY_PARTIAL:
                return self.flexibility.partiallySourceRedundantKeysRatio
            
            elif redundancyType is RedundancyType.SOURCE_FLEXIBILITY_BOTH:
                return self.flexibility.sourceRedundantKeysRatio + self.flexibility.partiallySourceRedundantKeysRatio
            
            
            else:
                raise ValueError("This type of redundancy is unknown: " + str(redundancyType))
            
        except AttributeError:
            raise ValueError("When constructing the redundancy object, you excluded this type of redundancy!")
    
    def getRedundantKeys(self, redundancyType: RedundancyType = RedundancyType.default) -> Set[Element]:
        """
        Get redundant key elements.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation.
        
        Returns
        -------
        Set[Element]
        
        Raises
        ------
        ValueError
            If `onlyType` was given in the contructor, but metrics of another type of redundancy are to be returned here.
        """
        try:
            if redundancyType is RedundancyType.ROBUSTNESS:
                return self.robustness.redundantKeys
            
            elif redundancyType is RedundancyType.ROBUSTNESS_PARTIAL:
                return self.robustness.partiallyRedundantKeys            
            
            elif redundancyType is RedundancyType.ROBUSTNESS_BOTH:
                return self.robustness.redundantKeys.union( self.robustness.partiallyRedundantKeys )
            
            
            elif redundancyType is RedundancyType.FLEXIBILITY:
                return self.flexibility.redundantKeys
            
            elif redundancyType is RedundancyType.FLEXIBILITY_PARTIAL:
                return self.flexibility.partiallyRedundantKeys
            
            elif redundancyType is RedundancyType.FLEXIBILITY_BOTH:
                return self.flexibility.redundantKeys.union( self.flexibility.partiallyRedundantKeys )
            
            
            elif redundancyType is RedundancyType.TARGET_FLEXIBILITY:
                return self.flexibility.targetRedundantKeys
            
            elif redundancyType is RedundancyType.TARGET_FLEXIBILITY_PARTIAL:
                return self.flexibility.partiallyTargetRedundantKeys
            
            elif redundancyType is RedundancyType.TARGET_FLEXIBILITY_BOTH:
                return self.flexibility.targetRedundantKeys.union( self.flexibility.partiallyTargetRedundantKeys )
            
            
            elif redundancyType is RedundancyType.SOURCE_FLEXIBILITY:
                return self.flexibility.sourceRedundantKeys
            
            elif redundancyType is RedundancyType.SOURCE_FLEXIBILITY_PARTIAL:
                return self.flexibility.partiallySourceRedundantKeys
            
            elif redundancyType is RedundancyType.SOURCE_FLEXIBILITY_BOTH:
                return self.flexibility.sourceRedundantKeys.union( self.flexibility.partiallySourceRedundantKeys )
            
            
            else:
                raise ValueError("This type of redundancy is unknown: " + str(redundancyType))
        
        except AttributeError:
            raise ValueError("When constructing the redundancy object, you excluded this type of redundancy!")
        
    def getRedundancyPaths(self, redundancyType: RedundancyType = RedundancyType.default) -> Set[Path]:
        """
        Get all paths providing redundancy.
        
        This always includes partial redundancy, use :func:`getRedundancyPathsForKey` if you want to differentiate.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        Set[Path]
        
        Raises
        ------
        ValueError
            If `onlyType` was given in the contructor, but metrics of another type of redundancy are to be returned here.
        """
        try:
            if redundancyType in [RedundancyType.ROBUSTNESS, RedundancyType.ROBUSTNESS_PARTIAL, RedundancyType.ROBUSTNESS_BOTH]:
                return self.robustness.paths
            
            elif redundancyType in [RedundancyType.FLEXIBILITY, RedundancyType.FLEXIBILITY_PARTIAL, RedundancyType.FLEXIBILITY_BOTH]:
                return self.flexibility.paths
            
            elif redundancyType in [RedundancyType.TARGET_FLEXIBILITY, RedundancyType.TARGET_FLEXIBILITY_PARTIAL, RedundancyType.TARGET_FLEXIBILITY_BOTH]:
                return self.flexibility.targetPaths
            
            elif redundancyType in [RedundancyType.SOURCE_FLEXIBILITY, RedundancyType.SOURCE_FLEXIBILITY_PARTIAL, RedundancyType.SOURCE_FLEXIBILITY_BOTH]:
                return self.flexibility.sourcePaths
            
            else:
                raise ValueError("This type of redundancy is unknown: " + str(redundancyType))
        
        except AttributeError:
            raise ValueError("When constructing the redundancy object, you excluded this type of redundancy!")
    
    def getRedundancyPathsForKey(self, redundancyType: RedundancyType = RedundancyType.default) -> Dict[Element, Set[Path]]:
        """
        Get paths providing redundancy, keyed by the key element they provide redundancy for.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        Dict[Element, Set[Path]]
        
        Raises
        ------
        ValueError
            If `onlyType` was given in the contructor, but metrics of another type of redundancy are to be returned here.
        """
        try:
            if redundancyType is RedundancyType.ROBUSTNESS:
                return self.robustness.redundantKeyPaths
            
            elif redundancyType is RedundancyType.ROBUSTNESS_PARTIAL:
                return self.robustness.partiallyRedundantKeyPaths
            
            elif redundancyType is RedundancyType.ROBUSTNESS_BOTH:
                return self.robustness.redundantKeyPaths.update( self.robustness.partiallyRedundantKeyPaths )
            
            
            elif redundancyType is RedundancyType.FLEXIBILITY:
                return self.flexibility.redundantKeyPaths
            
            elif redundancyType is RedundancyType.FLEXIBILITY_PARTIAL:
                return self.flexibility.partiallyRedundantKeyPaths
            
            elif redundancyType is RedundancyType.FLEXIBILITY_BOTH:
                return self.flexibility.redundantKeyPaths.update( self.flexibility.partiallyRedundantKeyPaths )
            
            
            elif redundancyType is RedundancyType.TARGET_FLEXIBILITY:
                return self.flexibility.targetRedundantKeyPaths
            
            elif redundancyType is RedundancyType.TARGET_FLEXIBILITY_PARTIAL:
                return self.flexibility.partiallyTargetRedundantKeyPaths
            
            elif redundancyType is RedundancyType.TARGET_FLEXIBILITY_BOTH:
                return self.flexibility.targetRedundantKeyPaths.update( self.flexibility.partiallyTargetRedundantKeyPaths )
            
            
            elif redundancyType is RedundancyType.SOURCE_FLEXIBILITY:
                return self.flexibility.sourceRedundantKeyPaths
            
            elif redundancyType is RedundancyType.SOURCE_FLEXIBILITY_PARTIAL:
                return self.flexibility.partiallySourceRedundantKeyPaths
            
            elif redundancyType is RedundancyType.SOURCE_FLEXIBILITY_BOTH:
                return self.flexibility.sourceRedundantKeyPaths.update( self.flexibility.partiallySourceRedundantKeyPaths )
            
            
            else:
                raise ValueError("This type of redundancy is unknown: " + str(redundancyType))
        
        except AttributeError:
            raise ValueError("When constructing the redundancy object, you excluded this type of redundancy!")




class RedundancyContribution():

    def __init__(self, redundancy: Redundancy, specialKeys: Set[Element]):
        """
        Contribution to redundancy, consisting of :class:`Flexibility` and class:`Robustness`, accountable to edges with `specialKeys`.
        
        Allows to answer the question how much certain key elements contribute to the flexibility/robustness of a graph.
        
        Parameters
        ----------
        redundancy : Redundancy
        specialKeys : Set[Element]
            Set of key elements viewed to be somehow special. One type of special could be 'neofunctionalised'.
        
        Attributes
        ----------
        self.flexibilityContribution : FlexibilityContribution
        self.robustnessContribution : RobustnessContribution
        """
        self.redundancy = redundancy
        onlyType = redundancy.onlyType
        
        if onlyType is None or onlyType.value[0] is Flexibility:
            self.flexibilityContribution = FlexibilityContribution(redundancy.flexibility, specialKeys)
        
        if onlyType is None or onlyType.value[0] is Robustness:
            self.robustnessContribution = RobustnessContribution(redundancy.robustness, specialKeys)
    
    @classmethod
    def fromGraph(cls, graph: DirectedMultiGraph, specialKeys: Set[Element]):
        """
        Create RedundancyContribution object from `graph`.
        
        Parameters
        ----------
        graph : DirectedMultiGraph
        specialKeys : Set[Element]
        
        Returns
        -------
        RedundancyContribution
        """
        return cls(Redundancy(graph), specialKeys)
    
    
    
    def getKeyContributionRatio(self, redundancyType: RedundancyType = RedundancyType.default) -> float:
        """
        Get ratio of contribution to keys' redundancy, to all keys.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation.
        
        Returns
        -------
        float
            Ratio of redundant keys which have a special key on at least one of their alternative paths.
        
        Raises
        ------
        ValueError
            If `onlyType` was given in the contructor of the underlying redundancy object, but metrics of another type of redundancy are to be returned here.
        """
        try:
            if redundancyType is RedundancyType.ROBUSTNESS:
                return self.robustnessContribution.redundantKeysWithSpecialKeyOnPathsRatio            
            
            elif redundancyType is RedundancyType.ROBUSTNESS_PARTIAL:
                return self.robustnessContribution.partiallyRedundantKeysWithSpecialKeyOnPathsRatio
            
            elif redundancyType is RedundancyType.ROBUSTNESS_BOTH:
                return self.robustnessContribution.redundantKeysWithSpecialKeyOnPathsRatio + self.robustnessContribution.partiallyRedundantKeysWithSpecialKeyOnPathsRatio
            
            
            elif redundancyType is RedundancyType.FLEXIBILITY:
                return self.flexibilityContribution.redundantKeysWithSpecialKeyOnPathsRatio
            
            elif redundancyType is RedundancyType.FLEXIBILITY_PARTIAL:
                return self.flexibilityContribution.partiallyRedundantKeysWithSpecialKeyOnPathsRatio
            
            elif redundancyType is RedundancyType.FLEXIBILITY_BOTH:
                return self.flexibilityContribution.redundantKeysWithSpecialKeyOnPathsRatio + self.flexibilityContribution.partiallyRedundantKeysWithSpecialKeyOnPathsRatio
            
            
            elif redundancyType is RedundancyType.TARGET_FLEXIBILITY:
                return self.flexibilityContribution.targetRedundantKeysWithSpecialKeyOnPathsRatio
            
            elif redundancyType is RedundancyType.TARGET_FLEXIBILITY_PARTIAL:
                return self.flexibilityContribution.partiallyTargetRedundantKeysWithSpecialKeyOnPathsRatio
            
            elif redundancyType is RedundancyType.TARGET_FLEXIBILITY_BOTH:
                return self.flexibilityContribution.targetRedundantKeysWithSpecialKeyOnPathsRatio + self.flexibilityContribution.partiallyTargetRedundantKeysWithSpecialKeyOnPathsRatio
            
            
            elif redundancyType is RedundancyType.SOURCE_FLEXIBILITY:
                return self.flexibilityContribution.sourceRedundantKeysWithSpecialKeyOnPathsRatio
            
            elif redundancyType is RedundancyType.SOURCE_FLEXIBILITY_PARTIAL:
                return self.flexibilityContribution.partiallySourceRedundantKeysWithSpecialKeyOnPathsRatio
            
            elif redundancyType is RedundancyType.SOURCE_FLEXIBILITY_BOTH:
                return self.flexibilityContribution.sourceRedundantKeysWithSpecialKeyOnPathsRatio + self.flexibilityContribution.partiallySourceRedundantKeysWithSpecialKeyOnPathsRatio
            
            
            else:
                raise ValueError("This type of redundancy is unknown: " + str(redundancyType))
        
        except AttributeError:
            raise ValueError("When constructing the redundancy object, you excluded this type of redundancy!")
    
    def getContributedKeysForSpecial(self, redundancyType: RedundancyType = RedundancyType.default) -> Dict[Element, Set[Element]]:
        """
        Get keys for which a special key contributes to redundancy.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation.
        
        Returns
        -------
        Dict[Element, Set[Element]]
            Sets of key elements which have a redundant path that contains a special key, keyed by the special key.
        
        Raises
        ------
        ValueError
            If `onlyType` was given in the contructor of the underlying redundancy object, but metrics of another type of redundancy are to be returned here.
        """
        try:
            if redundancyType is RedundancyType.ROBUSTNESS:
                return self.robustnessContribution.specialKeyOnRedundantKeysPaths            
            
            elif redundancyType is RedundancyType.ROBUSTNESS_PARTIAL:
                return self.robustnessContribution.specialKeyOnPartiallyRedundantKeysPaths
            
            elif redundancyType is RedundancyType.ROBUSTNESS_BOTH:
                dictA = self.robustnessContribution.specialKeyOnRedundantKeysPaths  
                dictB = self.robustnessContribution.specialKeyOnPartiallyRedundantKeysPaths
                return updateDictUpdatingValue(dictA, dictB)
            
            
            elif redundancyType is RedundancyType.FLEXIBILITY:
                return self.flexibilityContribution.specialKeyOnRedundantKeysPaths
            
            elif redundancyType is RedundancyType.FLEXIBILITY_PARTIAL:
                return self.flexibilityContribution.specialKeyOnPartiallyRedundantKeysPaths
            
            elif redundancyType is RedundancyType.FLEXIBILITY_BOTH:
                dictA = self.flexibilityContribution.specialKeyOnRedundantKeysPaths  
                dictB = self.flexibilityContribution.specialKeyOnPartiallyRedundantKeysPaths
                return updateDictUpdatingValue(dictA, dictB)
            
            
            elif redundancyType is RedundancyType.TARGET_FLEXIBILITY:
                return self.flexibilityContribution.specialKeyOnTargetRedundantKeysPaths            
            
            elif redundancyType is RedundancyType.TARGET_FLEXIBILITY_PARTIAL:
                return self.flexibilityContribution.specialKeyOnPartiallyTargetRedundantKeysPaths            
            
            elif redundancyType is RedundancyType.TARGET_FLEXIBILITY_BOTH:
                dictA = self.flexibilityContribution.specialKeyOnTargetRedundantKeysPaths  
                dictB = self.flexibilityContribution.specialKeyOnPartiallyTargetRedundantKeysPaths
                return updateDictUpdatingValue(dictA, dictB)
            
            
            elif redundancyType is RedundancyType.SOURCE_FLEXIBILITY:
                return self.flexibilityContribution.specialKeyOnSourceRedundantKeysPaths
            
            elif redundancyType is RedundancyType.SOURCE_FLEXIBILITY_PARTIAL:
                return self.flexibilityContribution.specialKeyOnPartiallySourceRedundantKeysPaths
            
            elif redundancyType is RedundancyType.SOURCE_FLEXIBILITY_BOTH:
                dictA = self.flexibilityContribution.specialKeyOnSourceRedundantKeysPaths  
                dictB = self.flexibilityContribution.specialKeyOnPartiallySourceRedundantKeysPaths
                return updateDictUpdatingValue(dictA, dictB)
            
            
            else:
                raise ValueError("This type of redundancy is unknown: " + str(redundancyType))
        
        except AttributeError:
            raise ValueError("When constructing the redundancy object, you excluded this type of redundancy!")
    
    def getContributingSpecialForKey(self, redundancyType: RedundancyType = RedundancyType.default) -> Dict[Element, Set[Element]]:
        """
        Get special keys which contribute to redundancy of a key.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation.
        
        Returns
        -------
        Dict[Element, Set[Element]]
            Sets of special key elements which contribute to a redundant path of a key, keyed by that key.
        
        Raises
        ------
        ValueError
            If `onlyType` was given in the contructor of the underlying redundancy object, but metrics of another type of redundancy are to be returned here.
        """
        try:
            if redundancyType is RedundancyType.ROBUSTNESS:
                return self.robustnessContribution.redundantKeySpecialKeysOnPaths            
            
            elif redundancyType is RedundancyType.ROBUSTNESS_PARTIAL:
                return self.robustnessContribution.partiallyRedundantKeySpecialKeysOnPaths
            
            elif redundancyType is RedundancyType.ROBUSTNESS_BOTH:
                return self.robustnessContribution.redundantKeySpecialKeysOnPaths.update( self.robustnessContribution.partiallyRedundantKeySpecialKeysOnPaths )
            
            
            elif redundancyType is RedundancyType.FLEXIBILITY:
                return self.flexibilityContribution.redundantKeySpecialKeysOnPaths
            
            elif redundancyType is RedundancyType.FLEXIBILITY_PARTIAL:
                return self.flexibilityContribution.partiallyRedundantKeySpecialKeysOnPaths
            
            elif redundancyType is RedundancyType.FLEXIBILITY_BOTH:
                return self.flexibilityContribution.redundantKeySpecialKeysOnPaths.update( self.flexibilityContribution.partiallyRedundantKeySpecialKeysOnPaths )
            
            
            elif redundancyType is RedundancyType.TARGET_FLEXIBILITY:
                return self.flexibilityContribution.targetRedundantKeySpecialKeysOnPaths
            
            elif redundancyType is RedundancyType.TARGET_FLEXIBILITY_PARTIAL:
                return self.flexibilityContribution.partiallyTargetRedundantKeySpecialKeysOnPaths
            
            elif redundancyType is RedundancyType.TARGET_FLEXIBILITY_BOTH:
                return self.flexibilityContribution.targetRedundantKeySpecialKeysOnPaths.update( self.flexibilityContribution.partiallyTargetRedundantKeySpecialKeysOnPaths )
            
            
            elif redundancyType is RedundancyType.SOURCE_FLEXIBILITY:
                return self.flexibilityContribution.sourceRedundantKeySpecialKeysOnPaths
            
            elif redundancyType is RedundancyType.SOURCE_FLEXIBILITY_PARTIAL:
                return self.flexibilityContribution.partiallySourceRedundantKeySpecialKeysOnPaths
            
            elif redundancyType is RedundancyType.SOURCE_FLEXIBILITY_BOTH:
                return self.flexibilityContribution.sourceRedundantKeySpecialKeysOnPaths.update( self.flexibilityContribution.partiallySourceRedundantKeySpecialKeysOnPaths )
            
            
            else:
                raise ValueError("This type of redundancy is unknown: " + str(redundancyType))
        
        except AttributeError:
            raise ValueError("When constructing the redundancy object, you excluded this type of redundancy!")
        
    def getContributedPaths(self, redundancyType: RedundancyType = RedundancyType.default) -> Set[MarkedPath]:
        """
        Get all paths on which any special key contributes to redundancy.
        
        This always includes partial redundancy, use :func:`getContributedPathsForKey` if you want to differentiate.

        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        Set[MarkedPath]
            Set of marked redundant paths that contain a special key.
        
        Raises
        ------
        ValueError
            If `onlyType` was given in the contructor of the underlying redundancy object, but metrics of another type of redundancy are to be returned here.
        """
        try:
            if redundancyType in [RedundancyType.ROBUSTNESS, RedundancyType.ROBUSTNESS_PARTIAL, RedundancyType.ROBUSTNESS_BOTH]:
                return self.robustnessContribution.pathsWithSpecialKeys
            
            elif redundancyType in [RedundancyType.FLEXIBILITY, RedundancyType.FLEXIBILITY_PARTIAL, RedundancyType.FLEXIBILITY_BOTH]:
                return self.flexibilityContribution.pathsWithSpecialKeys
            
            elif redundancyType in [RedundancyType.TARGET_FLEXIBILITY, RedundancyType.TARGET_FLEXIBILITY_PARTIAL, RedundancyType.TARGET_FLEXIBILITY_BOTH]:
                return self.flexibilityContribution.targetPathsWithSpecialKeys
            
            elif redundancyType in [RedundancyType.SOURCE_FLEXIBILITY, RedundancyType.SOURCE_FLEXIBILITY_PARTIAL, RedundancyType.SOURCE_FLEXIBILITY_BOTH]:
                return self.flexibilityContribution.sourcePathsWithSpecialKeys
            
            else:
                raise ValueError("This type of redundancy is unknown: " + str(redundancyType))
        
        except AttributeError:
            raise ValueError("When constructing the redundancy object, you excluded this type of redundancy!")
    
    def getContributedPathsForKey(self, redundancyType: RedundancyType = RedundancyType.default) -> Dict[Element, Set[MarkedPath]]:
        """
        Get paths on which a special key contributes to redundancy, keyed by the key element they provide redundancy for.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        Dict[Element, Set[MarkedPath]]
            Set of marked redundant paths that contain a special key, keyed by the key element they provide redundancy for.
        
        Raises
        ------
        ValueError
            If `onlyType` was given in the contructor of the underlying redundancy object, but metrics of another type of redundancy are to be returned here.
        """
        try:
            if redundancyType is RedundancyType.ROBUSTNESS:
                return self.robustnessContribution.redundantKeyPathsWithSpecialKey
            
            elif redundancyType is RedundancyType.ROBUSTNESS_PARTIAL:
                return self.robustnessContribution.partiallyRedundantKeyPathsWithSpecialKey
            
            elif redundancyType is RedundancyType.ROBUSTNESS_BOTH:
                return self.robustnessContribution.redundantKeyPathsWithSpecialKey.update( self.robustnessContribution.partiallyRedundantKeyPathsWithSpecialKey )
            
            
            elif redundancyType is RedundancyType.FLEXIBILITY:
                return self.flexibilityContribution.redundantKeyPathsWithSpecialKey
            
            elif redundancyType is RedundancyType.FLEXIBILITY_PARTIAL:
                return self.flexibilityContribution.partiallyRedundantKeyPathsWithSpecialKey
            
            elif redundancyType is RedundancyType.FLEXIBILITY_BOTH:
                return self.flexibilityContribution.redundantKeyPathsWithSpecialKey.update( self.flexibilityContribution.partiallyRedundantKeyPathsWithSpecialKey )
            
            
            elif redundancyType is RedundancyType.TARGET_FLEXIBILITY:
                return self.flexibilityContribution.targetRedundantKeyPathsWithSpecialKey
            
            elif redundancyType is RedundancyType.TARGET_FLEXIBILITY_PARTIAL:
                return self.flexibilityContribution.partiallyTargetRedundantKeyPathsWithSpecialKey
            
            elif redundancyType is RedundancyType.TARGET_FLEXIBILITY_BOTH:
                return self.flexibilityContribution.targetRedundantKeyPathsWithSpecialKey.update( self.flexibilityContribution.partiallyTargetRedundantKeyPathsWithSpecialKey )
            
            
            elif redundancyType is RedundancyType.SOURCE_FLEXIBILITY:
                return self.flexibilityContribution.sourceRedundantKeyPathsWithSpecialKey
            
            elif redundancyType is RedundancyType.SOURCE_FLEXIBILITY_PARTIAL:
                return self.flexibilityContribution.partiallySourceRedundantKeyPathsWithSpecialKey
            
            elif redundancyType is RedundancyType.SOURCE_FLEXIBILITY_BOTH:
                return self.flexibilityContribution.sourceRedundantKeyPathsWithSpecialKey.update( self.flexibilityContribution.partiallySourceRedundantKeyPathsWithSpecialKey )
            
            
            else:
                raise ValueError("This type of redundancy is unknown: " + str(redundancyType))
        
        except AttributeError:
            raise ValueError("When constructing the redundancy object, you excluded this type of redundancy!")
    
    


class Comparison():
    
    def __init__(self, graphA: DirectedMultiGraph, graphB: DirectedMultiGraph):
        """
        Compare redundancy between two graphs.
        
        Parameters
        ----------
        graphA : DirectedMultiGraph
            First graph to measure redundancy metrics for.
        graphB : DirectedMultiGraph
            Second graph to measure redundancy metrics for.
        
        Attributes
        ----------
        self.graphA : DirectedMultiGraph
        self.redundancyA : Redundancy
        
        self.graphB : DirectedMultiGraph
        self.redundancyB : Redundancy
        """
        self.graphA = graphA
        self.redundancyA = Redundancy(graphA)
        self.graphB = graphB
        self.redundancyB = Redundancy(graphB)
    
    @classmethod
    def fromOrganismGroups(cls, groupA: Organism.Group, groupB: Organism.Group, majorityPercentage = None):
        """
        Compare redundancy between two groups' core metabolisms.
        
        Parameters
        ----------
        groupA : Organism.Group
            First group from which to extract the graph.
        groupB : Organism.Group
            Second group from which to extract the graph.
        majorityPercentage : float, optional
            If *None*, use collective EC graph. If not *None*, use majority EC graph with `majorityPercentage`% majority.
        
        Returns
        -------
        Comparison
        """
        if majorityPercentage is None:
            return cls(groupA.collectiveEcGraph(), groupB.collectiveEcGraph())
        else:
            return cls(groupA.majorityEcGraph(majorityPercentage), groupB.majorityEcGraph(majorityPercentage))
    
    @classmethod
    def fromCladePair(cls, cladePair: CladePair, majorityPercentage = None):
        """
        Compare redundancy between two clades' core metabolisms.
        
        Parameters
        ----------
        cladePair : CladePair
            The pair of clades from which to extract the graph.
        majorityPercentage : float, optional
            If *None*, use collective EC graph. If not *None*, use majority EC graph with `majorityPercentage`% majority.
        
        Returns
        -------
        Comparison
        """
        return cls.fromOrganismGroups(cladePair.parentClade.group, cladePair.childClade.group, majorityPercentage)
    
    
    
    def getLostRedundancyRatio(self, redundancyType: RedundancyType = RedundancyType.default) -> float:
        """
        Get ratio of keys, which have lost redundancy from graph A to graph B, to all keys.
        
        This only counts keys which exist in both graphs, for both 'redundant keys' and 'all keys'.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        float
        """
        return self._getRedundancyRatio(redundancyType, -1)
        
    def getConservedRedundancyRatio(self, redundancyType: RedundancyType = RedundancyType.default) -> float:
        """
        Get ratio of keys, which have conserved redundancy from graph A to graph B, to all keys.
        
        This only counts keys which exist in both graphs, for both 'redundant keys' and 'all keys'.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        float
        """
        return self._getRedundancyRatio(redundancyType, 0)
        
    def getAddedRedundancyRatio(self, redundancyType: RedundancyType = RedundancyType.default) -> float:
        """
        Get ratio of keys, which have become redundant from graph A to graph B, to all keys.
        
        This only counts keys which exist in both graphs, for both 'redundant keys' and 'all keys'.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        float
        """
        return self._getRedundancyRatio(redundancyType, 1)
    
    def _getRedundancyRatio(self, redundancyType: RedundancyType,  direction: int) -> float:
        keysA = self.graphA.getEdgeKeys()
        keysB = self.graphB.getEdgeKeys()        
        keysBoth = keysA.intersection(keysB)
        
        redundantKeysA = self.redundancyA.getRedundantKeys(redundancyType)
        redundantKeysB = self.redundancyB.getRedundantKeys(redundancyType)
        
        # find relevant redundant keys
        relevantKeys = set()
        for key in keysBoth:
            if direction == -1: # lost
                if key in redundantKeysA and not key in redundantKeysB:
                    relevantKeys.add(key)
                
            elif direction == 0: # conserved
                if key in redundantKeysA and redundantKeysB:
                    relevantKeys.add(key)
                    
            elif direction == 1: # added
                if key in redundantKeysB and not key in redundantKeysA:
                    relevantKeys.add(key)
        
        return len(relevantKeys)/len(keysBoth)
    
    
    
    
    def getLostRedundancyKeys(self, redundancyType: RedundancyType = RedundancyType.default) -> Set[Element]:
        """
        Get keys which have lost redundancy from graph A to graph B.
        
        This only counts keys which exist in both graphs.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        Set[Element]
        """
        self._getRedundancyKeys(redundancyType, -1)
    
    def getConservedRedundancyKeys(self, redundancyType: RedundancyType = RedundancyType.default) -> Set[Element]:
        """
        Get keys which have conserved redundancy from graph A to graph B.
        
        This only counts keys which exist in both graphs.
        Beware, this only means there is some redundancy in A and B, not the exact same paths providing redundancy in A and B!
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        Set[Element]
        """
        self._getRedundancyKeys(redundancyType, 0)
    
    def getAddedRedundancyKeys(self, redundancyType: RedundancyType = RedundancyType.default) -> Set[Element]:
        """
        Get keys which have become redundant from graph A to graph B.
        
        This only counts keys which exist in both graphs.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        Set[Element]
        """
        self._getRedundancyKeys(redundancyType, 1)
        
    def _getRedundancyKeys(self, redundancyType: RedundancyType, direction: int) -> Set[Element]:
        keysA = self.graphA.getEdgeKeys()
        keysB = self.graphB.getEdgeKeys()        
        keysBoth = keysA.intersection(keysB)
        
        redundantKeysA = self.redundancyA.getRedundantKeys(redundancyType)
        redundantKeysB = self.redundancyB.getRedundantKeys(redundancyType)
        
        # find relevant redundant keys
        relevantKeys = set()
        for key in keysBoth:
            if direction == -1: # lost
                if key in redundantKeysA and not key in redundantKeysB:
                    relevantKeys.add(key)
                
            elif direction == 0: # conserved
                if key in redundantKeysA and redundantKeysB:
                    relevantKeys.add(key)
                    
            elif direction == 1: # added
                if key in redundantKeysB and not key in redundantKeysA:
                    relevantKeys.add(key)
        
        return relevantKeys
    
    
    
    def getLostRedundancyPathsForKey(self, redundancyType: RedundancyType = RedundancyType.default) -> Dict[Element, Set[Path]]:
        """
        Get alternative paths of keys which have lost redundancy from graph A to graph B.
        
        This only counts keys which exist in both graphs.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        Dict[Element, Set[Path]]
            Set of alternative paths belonging to a key which lost redundancy, keyed by this key.
        """
        return self._getRedundancyPathsForKey(redundancyType, -1)
    
    def getConservedRedundancyKeysPathsForKey(self, redundancyType: RedundancyType = RedundancyType.default) -> Dict[Element, Tuple[Set[Path], Set[Path], Set[Path]]]:
        """
        Get tuple of alternative paths of keys which have conserved redundancy from graph A to graph B.
        
        This only counts keys which exist in both graphs.
        However, even if a key is redundant in both graphs, it does not have to be redundant due to the exact same paths.
        This is why the tuple has three sets, the first for paths which only exist only in A, the second for paths which exist in both, and the third for paths which exist only in B.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        Dict[Element, Tuple[Set[Path], Set[Path], Set[Path]]]
            Tuple of sets of alternative paths belonging to a key which conserved redundancy, keyed by this key.
            The first set of the tuple contains paths which only exists in graph A.
            The second set of the tuple contains paths which exist in both graphs.
            The third set of the tuple contains paths which only exists in graph B.
        """
        return self._getRedundancyPathsForKey(redundancyType, 0)
    
    def getAddedRedundancyKeysPathsForKey(self, redundancyType: RedundancyType = RedundancyType.default) -> Dict[Element, Set[Path]]:
        """
        Get alternative paths of keys which have become redundant from graph A to graph B.
        
        This only counts keys which exist in both graphs.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        Dict[Element, Set[Path]]
            Set of alternative paths belonging to a key which has become redundant, keyed by this key.
        """
        return self._getRedundancyPathsForKey(redundancyType, 1)
    
    def _getRedundancyPathsForKey(self, redundancyType: RedundancyType, direction: int):
        redundantKeys = self._getRedundancyKeys(redundancyType, direction)
        
        if direction == -1: # lost
            paths = self.redundancyA.getRedundancyPathsForKey(redundancyType)
            
        elif direction == 0: # conserved
            pathsA = self.redundancyA.getRedundancyPathsForKey(redundancyType)
            pathsB = self.redundancyB.getRedundancyPathsForKey(redundancyType)
            
        elif direction == 1: # added
            paths = self.redundancyB.getRedundancyPathsForKey(redundancyType)
        
        resultPaths = dict()
        for key in redundantKeys:
            
            if direction == 0: # conserved
                currentPathsA = pathsA[key]
                currentPathsB = pathsB[key]
                
                currentPathsBoth = currentPathsA.intersection(currentPathsB)
                currentPathsA.difference_update(currentPathsBoth)
                currentPathsB.difference_update(currentPathsBoth)
                
                resultPaths[key] = (currentPathsA, currentPathsBoth, currentPathsB)
                
            else: # lost or added
                resultPaths[key] = paths[key]
        
        return resultPaths
    
    
    
    def getLostRedundancyPaths(self, redundancyType: RedundancyType = RedundancyType.default) -> Set[Path]:
        """
        Get all alternative paths which have been lost from graph A to graph B.
        
        This counts all paths, no matter which key they provide redundancy for!
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        Set[Path]
            Set of alternative paths which have been lost in graph B.
        """
        return self._getRedundancyPaths(redundancyType, -1)
    
    def getConservedRedundancyPaths(self, redundancyType: RedundancyType = RedundancyType.default) -> Set[Path]:
        """
        Get all alternative paths which have been conserved between graph A and graph B.
        
        This counts all paths, no matter which key they provide redundancy for!
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        Set[Path]
            Set of alternative paths which are both in graph A and in graph B.
        """
        return self._getRedundancyPaths(redundancyType, 0)
    
    def getAddedRedundancyPaths(self, redundancyType: RedundancyType = RedundancyType.default) -> Set[Path]:
        """
        Get all alternative paths which have been added from graph A to graph B.
        
        This counts all paths, no matter which key they provide redundancy for!
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        Set[Path]
            Set of alternative paths which have been added in graph B.
        """
        return self._getRedundancyPaths(redundancyType, 1)
    
    def _getRedundancyPaths(self, redundancyType: RedundancyType, direction: int) -> Set[Path]:
        pathsA = self.redundancyA.getRedundancyPaths(redundancyType)
        pathsB = self.redundancyB.getRedundancyPaths(redundancyType)
        
        if direction == -1: # lost
            return pathsA.difference(pathsB)
            
        elif direction == 0: # conserved
            return pathsA.intersection(pathsB)
                
        elif direction == 1: # added
            return pathsB.difference(pathsA)
        
    
        
        
        


class ContributionComparison():
    
    def __init__(self, comparison: Comparison, specialKeysA: Set[Element], specialKeysB: Set[Element]):
        """
        Compare contribution to redundancy between two graphs.
        
        Allows to answer the question how much certain key elements contribute to certain comparing aspects of the flexibility/robustness of two graphs.
        Special keys will usually differ between the two graphs, however, they are also allowed to be the same set.
        
        Parameters
        ----------
        comparison : Comparison
            You first have to calculate a comparison object using your two graphs.
        specialKeysA : Set[Element]
            Set of key elements of graph A in `comparison` viewed to be somehow special. One type of special could be 'neofunctionalised'.
        specialKeysB : Set[Element]
            Set of key elements of graph B in `comparison` viewed to be somehow special. One type of special could be 'neofunctionalised'.
        
        Attributes
        ----------
        self.comparison : Comparison
        self.redundancyContributionA : RedundancyContribution
        self.redundancyContributionB : RedundancyContribution
        """
        self.comparison = comparison
        self.redundancyContributionA = RedundancyContribution(comparison.redundancyA, specialKeysA)
        self.redundancyContributionB = RedundancyContribution(comparison.redundancyB, specialKeysB)
    
    
    def getLostRedundancyKeyContributionRatio(self, redundancyType: RedundancyType = RedundancyType.default) -> float:
        """
        Get ratio of contribution to rendundancy of keys, which have lost redundancy.
        
        This only counts keys which exist in both graphs. The ratio is to the number of redundant keys, not to all keys.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation.
        
        Returns
        -------
        float
            Ratio of contribution to redundant keys, which have lost redundancy from graph A to graph B, and for which a special key contributes to redundancy.
        """
        return self._getRedundancyKeyContributionRatio(redundancyType, -1)
    
    def getConservedRedundancyKeyContributionRatio(self, redundancyType: RedundancyType = RedundancyType.default) -> float:
        """
        Get ratio of contribution to rendundancy of keys, which have conserved redundancy.
        
        A *or* B may have a special key on an alternative path of a redundant key, for the key to be counted here.
        It is **not** necessary, that both A *and* B have such a special key.
        This only counts keys which exist in both graphs. The ratio is to the number of redundant keys, not to all keys.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation.
        
        Returns
        -------
        float
            Ratio of contribution to redundant keys, which have conserved redundancy from graph A to graph B, and for which a special key contributes to redundancy.
        """
        return self._getRedundancyKeyContributionRatio(redundancyType, 0)
    
    def getAddedRedundancyKeyContributionRatio(self, redundancyType: RedundancyType = RedundancyType.default) -> float:
        """
        Get ratio of contribution to rendundancy of keys, which have become redundant.
        
        This only counts keys which exist in both graphs. The ratio is to the number of redundant keys, not to all keys.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation.
        
        Returns
        -------
        float
            Ratio of contribution to redundant keys, which have become redundant from graph A to graph B, and for which a special key contributes to redundancy.
        """
        return self._getRedundancyKeyContributionRatio(redundancyType, 1)
    
    def _getRedundancyKeyContributionRatio(self, redundancyType: RedundancyType,  direction: int) -> float:
        redundantKeys = self.comparison._getRedundancyKeys(redundancyType, direction)
        
        if direction == -1: # lost
            specialForRedundantKeyA = self.redundancyContributionA.getContributingSpecialForKey(redundancyType)
            
        elif direction == 0: # conserved
            specialForRedundantKeyA = self.redundancyContributionA.getContributingSpecialForKey(redundancyType)
            specialForRedundantKeyB = self.redundancyContributionB.getContributingSpecialForKey(redundancyType)
                
        elif direction == 1: # added
            specialForRedundantKeyB = self.redundancyContributionB.getContributingSpecialForKey(redundancyType)
        
        contributedKeys = set() # which have a redundancy due to a special key
        for key in redundantKeys:
            if direction == -1: # lost
                specials = specialForRedundantKeyA.get(key, None)
                
            elif direction == 0: # conserved
                specialsA = specialForRedundantKeyA.get(key, None)
                specialsB = specialForRedundantKeyB.get(key, None)
                
                if specialsA is not None:
                    if specialsB is not None:
                        specialsA.update( specialsB ) # A or B can have a special key on alternative path of a redundant key, for the key to be reported here. It is **not** necessary, that both A and B have to have such a special key. 
                        specials = specialsA
                    else:
                        specials = specialsA
                else:
                    specials = specialsB
                    
            elif direction == 1: # added
                specials = specialForRedundantKeyB.get(key, None)
            
            if specials is not None: # has special key in redundant paths
                contributedKeys.add( key )
        
        return len(contributedKeys)/len(redundantKeys)
    
    
    
    
    def getContributedLostRedundancyKeysForSpecial(self, redundancyType: RedundancyType = RedundancyType.default) -> Dict[Element, Set[Element]]:
        """
        Get keys, which have lost redundancy, for which a special key contributes to redundancy.
        
        This only counts keys which exist in both graphs.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation.
        
        Returns
        -------
        Dict[Element, Set[Element]]
            Sets of key elements, which have lost redundancy from graph A to graph B, and which have a redundant path that contains a special key, keyed by the special key.
            All of these paths (except for maybe one) exist only in A.
        """
        self._getContributedRedundancyKeyForSpecial(redundancyType, -1)
    
    def getContributedConservedRedundancyKeysForSpecial(self, redundancyType: RedundancyType = RedundancyType.default) -> Dict[Element, Set[Element]]:
        """
        Get keys, which have conserved redundancy, for which a special key contributes to redundancy.
        
        A *or* B may have a special key on an alternative path of a redundant key, for the key to be reported here.
        It is **not** necessary, that both A *and* B have such a special key.
        This only counts keys which exist in both graphs.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation.
        
        Returns
        -------
        Dict[Element, Set[Element]]
            Sets of key elements, which have conserved redundancy from graph A to graph B, and which have a redundant path that contains a special key, keyed by the special key.
            All of these paths may exist in only A, only B, or in both. One occurence is enough to be reported here.
        """
        self._getContributedRedundancyKeyForSpecial(redundancyType, 0)
    
    def getContributedAddedRedundancyKeysForSpecial(self, redundancyType: RedundancyType = RedundancyType.default) -> Dict[Element, Set[Element]]:
        """
        Get keys, which have become redundant, for which a special key contributes to redundancy.
        
        This only counts keys which exist in both graphs.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation.
        
        Returns
        -------
        Dict[Element, Set[Element]]
            Sets of key elements, which have become redundant from graph A to graph B, and which have a redundant path that contains a special key, keyed by the special key.
            All of these paths (except for maybe one) exist only in B.
        """
        self._getContributedRedundancyKeyForSpecial(redundancyType, 1)
    
    def _getContributedRedundancyKeyForSpecial(self, redundancyType: RedundancyType,  direction: int) -> Dict[Element, Set[Element]]:
        redundantKeys = self.comparison._getRedundancyKeys(redundancyType, direction)
        contributedKeysForSpecial = dict() # which have a redundancy due to a special key
        
        if direction == -1: # lost
            redundantKeyForSpecialA = self.redundancyContributionA.getContributedKeysForSpecial(redundancyType)
            
            for special, keys in redundantKeyForSpecialA.items():
                overlappingKeys = keys.intersection(redundantKeys)
                
                if len(overlappingKeys) > 0:
                    contributedKeysForSpecial[special] = overlappingKeys
            
        elif direction == 0: # conserved
            redundantKeyForSpecialA = self.redundancyContributionA.getContributedKeysForSpecial(redundancyType)
            redundantKeyForSpecialB = self.redundancyContributionB.getContributedKeysForSpecial(redundancyType)
            
            for special, keys in redundantKeyForSpecialA.items():
                overlappingKeys = keys.intersection(redundantKeys)
                
                if len(overlappingKeys) > 0:
                    contributedKeysForSpecial[special] = overlappingKeys
            
            for special, keys in redundantKeyForSpecialB.items():
                overlappingKeys = keys.intersection(redundantKeys)
                
                if len(overlappingKeys) > 0:
                    currentSet = contributedKeysForSpecial.get(special, None)
                    
                    if currentSet is None:
                        currentSet = set()
                        contributedKeysForSpecial[special] = currentSet
                        
                    currentSet.update(overlappingKeys)
                
        elif direction == 1: # added
            redundantKeyForSpecialB = self.redundancyContributionB.getContributedKeysForSpecial(redundancyType)
            
            for special, keys in redundantKeyForSpecialB.items():
                overlappingKeys = keys.intersection(redundantKeys)
                
                if len(overlappingKeys) > 0:
                    contributedKeysForSpecial[special] = overlappingKeys
        
        return contributedKeysForSpecial
        
    
    
    
    def getContributedLostRedundancyPaths(self, redundancyType: RedundancyType = RedundancyType.default) -> Set[MarkedPath]:
        """
        Get all alternative paths which have been lost, and have a special key on them.
        
        This counts all paths, no matter which key they provide redundancy for!
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        Set[MarkedPath]
            Set of alternative paths which have been lost from graph A to graph B.
        """
        return self._getContributedRedundancyPaths(redundancyType, -1)
    
    def getContributedConservedRedundancyPaths(self, redundancyType: RedundancyType = RedundancyType.default) -> Set[MarkedPath]:
        """
        Get all alternative paths which have been conserved, and have a special key on them.
        
        This counts all paths, no matter which key they provide redundancy for!
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        Set[MarkedPath]
            Set of alternative paths which have been conserved from graph A to graph B.
        """
        return self._getContributedRedundancyPaths(redundancyType, 0)
    
    def getContributedAddedRedundancyPaths(self, redundancyType: RedundancyType = RedundancyType.default) -> Set[MarkedPath]:
        """
        Get all alternative paths which have become redundant, and have a special key on them.
        
        This counts all paths, no matter which key they provide redundancy for!
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        Set[MarkedPath]
            Set of alternative paths which have become redundant from graph A to graph B.
        """
        return self._getContributedRedundancyPaths(redundancyType, 1)
    
    def _getContributedRedundancyPaths(self, redundancyType: RedundancyType,  direction: int) -> Set[MarkedPath]:
        pathsA = self.redundancyContributionA.getContributedPaths(redundancyType)
        pathsB = self.redundancyContributionB.getContributedPaths(redundancyType)
        
        if direction == -1: # lost
            return pathsA.difference(pathsB)
            
        elif direction == 0: # conserved
            return pathsA.intersection(pathsB)
                
        elif direction == 1: # added
            return pathsB.difference(pathsA)
        
        
        
        
    
    
    def getContributedLostRedundancyPathsForKey(self, redundancyType: RedundancyType = RedundancyType.default) -> Dict[Element, Set[MarkedPath]]:
        """
        Get alternative paths of keys which have lost redundancy, and have a special key on them.
        
        This only counts keys which exist in both graphs.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        Dict[Element, Set[MarkedPath]]
            Set of alternative marked paths belonging to a key which lost redundancy from graph A to graph B, keyed by this key.
        """
        return self._getContributedRedundancyPathsForKey(redundancyType, -1)
    
    def getContributedConservedRedundancyPathsForKey(self, redundancyType: RedundancyType = RedundancyType.default) -> Dict[Element, Tuple[Set[MarkedPath], Set[MarkedPath], Set[MarkedPath]]]:
        """
        Get alternative paths of keys which have conserved redundancy, and have a special key on them.
        
        This only counts keys which exist in both graphs.
        However, even if a key is redundant in both graphs, it does not have to be redundant due to the exact same paths.
        This is why the tuple has three sets, the first for paths which only exist only in A, the second for paths which exist in both, and the third for paths which exist only in B.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        Dict[Element, Tuple[Set[MarkedPath], Set[MarkedPath], Set[MarkedPath]]]
            Tuple of sets of alternative marked paths belonging to a key which conserved redundancy from graph A to graph B, keyed by this key.
            The first set of the tuple contains paths which only exists in graph A.
            The second set of the tuple contains paths which exist in both graphs.
            The third set of the tuple contains paths which only exists in graph B.
        """
        return self._getContributedRedundancyPathsForKey(redundancyType, 0)
    
    def getContributedAddedRedundancyPathsForKey(self, redundancyType: RedundancyType = RedundancyType.default) -> Dict[Element, Set[MarkedPath]]:
        """
        Get alternative paths of keys which have become redundant, and have a special key on them.
        
        This only counts keys which exist in both graphs.
        
        Parameters
        ----------
        redundancyType : RedundancyType
            Type of redundancy to use for computation. For target-/source-flexibility, only the paths for target/source nodes are reported, paths of the respective other node are ignored.
        
        Returns
        -------
        Dict[Element, Set[MarkedPath]]
            Set of alternative marked paths belonging to a key which has become redundant from graph A to graph B, keyed by this key.
        """
        return self._getContributedRedundancyPathsForKey(redundancyType, 1)
    
    def _getContributedRedundancyPathsForKey(self, redundancyType: RedundancyType,  direction: int):
        redundantKeys = self.comparison._getRedundancyKeys(redundancyType, direction)
        
        if direction == -1: # lost
            pathsForKey = self.redundancyContributionA.getContributedPathsForKey(redundancyType)
        
        elif direction == 0: # conserved
            pathsForKeyA = self.redundancyContributionA.getContributedPathsForKey(redundancyType)
            pathsForKeyB = self.redundancyContributionB.getContributedPathsForKey(redundancyType)
            
        elif direction == 1: # added
            pathsForKey = self.redundancyContributionB.getContributedPathsForKey(redundancyType)
        
        markedPathsForKey = dict()
        for key in redundantKeys:
            
            if direction == 0: # conserved
                currentPathsA = pathsForKeyA[key]
                currentPathsB = pathsForKeyB[key]
                
                currentPathsBoth = currentPathsA.intersection(currentPathsB)
                currentPathsA.difference_update(currentPathsBoth)
                currentPathsB.difference_update(currentPathsBoth)
                
                markedPathsForKey[key] = (currentPathsA, currentPathsBoth, currentPathsB)
            
            else: # lost or added
                markedPathsForKey[key] = pathsForKey[key]
        
        return markedPathsForKey
    