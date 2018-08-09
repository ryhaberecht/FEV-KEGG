"""
Question
--------
Which percentage of EC numbers is shared between all subspecies of Escherichia coli K-12?

Method
------
- Get all metabolic pathways of all E. coli species from KEGG.
- For each species, combine all pathways to the metabolic network, by UNION operation.
- Convert this metabolic network into a substance-ecNumber graph.
- Combine all species' networks to a single unified E. coli network, by UNION operation.
- Combine all species' networks to a consensus network, by INTERSECT operation, leaving only substances and EC numbers that occur in all species.
- Print number of EC numbers in the unified network -> numberUnified.
- Print number of EC numbers in the consensus network -> numberConsensus.
- Print percentage numberConsensus/numberUnified.

Result
------

::

    unified: 763
    consensus: 699
    91.61205766710354%

Conclusion
----------
About 91% of all EC numbers present in all known Escherichia coli K-12 subspecies occur in all subspecies. Thus, 9% of EC numbers only occur in some of the subspecies.
This indicates that subspecies evolve by acquiring new enzymatic functionalities. Whether they stem from neofunctionalisation, horizontal gene transfer, or genesis remains to be investigated.
"""
from FEV_KEGG.Graph.SubstanceGraphs import SubstanceReactionGraph, SubstanceGeneGraph, SubstanceEcGraph
import FEV_KEGG.KEGG.Organism


if __name__ == '__main__':
    
    #- Get all metabolic pathways of all E. coli species from KEGG.
    eColiSpecies = FEV_KEGG.KEGG.Organism.Group('Escherichia coli K-12').getOrganisms()
    
    #- For each species, combine all pathways to the metabolic network, by UNION operation.
    speciesEcGraphs = []
    for species in eColiSpecies:
        speciesPathways = species.getMetabolicPathways()
        speciesSubstanceReactionGraph = SubstanceReactionGraph.fromPathway(speciesPathways)
    
        #- Convert this metabolic network into a substance-ecNumber graph.
        speciesSubstanceGeneGraph = SubstanceGeneGraph.fromSubstanceReactionGraph(speciesSubstanceReactionGraph)
        speciesSubstanceEcGraph = SubstanceEcGraph.fromSubstanceGeneGraph(speciesSubstanceGeneGraph)
        
        speciesEcGraphs.append(speciesSubstanceEcGraph)
    
    firstGraph = speciesEcGraphs.pop(0)
    
    #- Combine all species' networks to a single unified E. coli network, by UNION operation.
    unifiedEcGraph = firstGraph
    unifiedEcGraph = unifiedEcGraph.union(speciesEcGraphs)
    
    #- Combine all species' networks to a consensus network, by INTERSECT operation, leaving only substances and EC numbers that occur in all species.
    intersectedEcGraph = firstGraph
    intersectedEcGraph = intersectedEcGraph.intersection(speciesEcGraphs)
    
    
    #- Print number of EC numbers in the unified network -> numberUnified.
    output = []
    numberUnified = len( unifiedEcGraph.getECs() )
    output.append('unified: ' + str(numberUnified))
    
    #- Print number of EC numbers in the consensus network -> numberConsensus.
    numberConsensus = len( intersectedEcGraph.getECs() )
    output.append('consensus: ' + str(numberConsensus))
    
    #- Print percentage numberConsensus/numberUnified.
    output.append(str(numberConsensus/numberUnified * 100) + '%')
    
    for line in output:
        print(line)