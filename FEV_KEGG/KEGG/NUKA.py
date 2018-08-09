from FEV_KEGG.Graph import SubstanceGraphs
from FEV_KEGG.Graph.Elements import ReactionID, EcNumber
from FEV_KEGG.Graph.SubstanceGraphs import SubstanceEcGraph, SubstanceReactionGraph
from FEV_KEGG.KEGG.File import cache
from FEV_KEGG.KEGG.Organism import Organism
from FEV_KEGG.settings import verbosity as init_verbosity


class NUKA(object):
    
    def __init__(self):
        """
        This is a hypothetical 'complete' organism - NUKA - which possesses all EC numbers known to all metabolic KEGG pathways.
        
        Conversions to other graph types are not possible, because as a hypothetical organism, NUKA has no genes.
        
        Attributes
        ----------
        self.nameAbbreviation : str
        """
        self.nameAbbreviation = 'NUKA'
        
    
    @property
    @cache(folder_path = 'NUKA/graph', file_name = 'SubstanceReactionGraph')
    def substanceReactionGraph(self) -> SubstanceReactionGraph:
        """
        NUKA's substance-reaction graph.
        
        Returns
        -------
        SubstanceReactionGraph
            Contains all substrates/products and all reactions known to KEGG's metabolic pathways.
        
        Raises
        ------
        HTTPError
            If any underlying organism, pathway, or gene does not exist.
        URLError
            If connection to KEGG fails.
        
        Note
        ----
        This SubstanceReactionGraph can **NOT** be converted into a SubstanceGeneGraph, as the pathways do not contain gene information!
        """
        mockOrganism = Organism('ec') # 'ec' is not an organism abbreviation, but merely desribes that pathways shall contain EC numbers as edges. This returns the full pathways not specific to any species.
        pathwaysSet = mockOrganism.getMetabolicPathways(includeOverviewMaps = False)
        substanceReactionGraph = SubstanceGraphs.Conversion.KeggPathwaySet2SubstanceReactionGraph(pathwaysSet, localVerbosity = 0)
        substanceReactionGraph.name = 'Substance-Reaction NUKA'
        
        if init_verbosity > 0:
            print('calculated ' + substanceReactionGraph.name)
        
        return substanceReactionGraph
    
    
    @property
    @cache(folder_path = 'NUKA/graph', file_name = 'SubstanceEcGraph')
    def substanceEcGraph(self) -> SubstanceEcGraph:
        """
        NUKA's substance-EC graph.
        
        Returns
        -------
        SubstanceEcGraph
            Contains all substrates/products and all EC numbers known to KEGG's metabolic pathways.
            
        Raises
        ------
        HTTPError
            If any underlying organism, pathway, or gene does not exist.
        URLError
            If connection to KEGG fails.
        """
        return self._SubstanceReactionGraph2SubstanceEcGraph(self.substanceReactionGraph)
        
        
    def _SubstanceReactionGraph2SubstanceEcGraph(self, speciesSubstanceReactionGraph: SubstanceReactionGraph) -> SubstanceEcGraph:
        """
        Converts NUKA's substance-reaction graph into a substance-EC graph. Uses pathway information embedded into the graph object.
        
        Parameters
        ----------
        speciesSubstanceReactionGraph : SubstanceReactionGraph
            NUKA's substance-reaction graph.
        
        Returns
        -------
        SubstanceEcGraph
            NUKA's substance-EC graph.
        
        Warnings
        --------
        This function is special to NUKA and **MUST NOT** be used anywhere else!
        """
        # shallow-copy old graph to new graph
        graph = SubstanceEcGraph(speciesSubstanceReactionGraph.underlyingRawGraph)
        graph.name = 'Substance-EC NUKA'
        
        # create dict of replacements: reaction -> {EC numbers}
        replacementDict = dict()
        
        # for each embedded pathway, get list of 'enzyme' entries
        for pathway in speciesSubstanceReactionGraph.pathwaySet:
            ecEntryList = [e for e in pathway.entries.values() if e.type == 'enzyme']
            
            # for each EC number, get reactions in which it is involved
            for ecEntry in ecEntryList:
                reactionIDList = ecEntry.reaction.split()
                if len(reactionIDList) > 0: # filter EC numbers not associated with any reaction
                    ecNumberList = ecEntry.name.split()
                    
                    # replace each reaction with its associated EC number
                    for reactionID in reactionIDList:
                        reactionName = reactionID.split(':', 1)[1]
                        reaction = ReactionID(reactionName)
                        
                        # save associated EC numbers in a set
                        ecNumberSet = set()
                        for ecNumberString in ecNumberList:
                            ecNumber = EcNumber(ecNumberString.replace('ec:', ''))
                            ecNumberSet.add(ecNumber)
                        
                        # update the replacement dict for the current reaction, adding the newly created EC number set
                        replacementSet = replacementDict.get(reaction, None)
                        if replacementSet == None or replacementSet.__class__ != set:
                            replacementSet = set()
                        replacementSet.update(ecNumberSet)
                        replacementDict[reaction] = replacementSet
        
        # get list of all reaction edges. Copy edge list to prevent changes in-place, which would NOT work
        edgeList = list(graph.getEdges())
            
        # replace reaction edges with EC number edges, using replacement dict
        for edge in edgeList:
            substrate, product, reaction = edge
            
            # delete old edge
            graph.removeEdge(substrate, product, reaction, False)
            
            # add new edges, according to replacement dict
            replacementSet = replacementDict[reaction]
            for ecNumber in replacementSet:
                graph.addEC(substrate, product, ecNumber, False)
        
        if init_verbosity > 0:
            print('calculated ' + graph.name)
        
        return graph
    