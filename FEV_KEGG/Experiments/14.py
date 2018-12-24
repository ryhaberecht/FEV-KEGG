"""
Question
--------
Which percentage of EC numbers is shared between all organisms of Escherichia coli K-12?

Method
------
- Get all metabolic pathways of all E. coli organisms from KEGG.
- For each organisms, combine all pathways to the metabolic network, by UNION operation.
- Convert this metabolic network into a substance-ecNumber graph.
- Combine all organisms' networks to a single unified E. coli network, by UNION operation.
- Combine all organisms' networks to a consensus network, by INTERSECT operation, leaving only substances and EC numbers that occur in all organisms.
- Print number of EC numbers in the unified network -> numberUnified.
- Print number of EC numbers in the consensus network -> numberConsensus.
- Print percentage numberConsensus/numberUnified.

Result
------

::
    
    unified: 617
    consensus: 560
    90.76175040518638%

Conclusion
----------
About 91% of all EC numbers present in all known Escherichia coli K-12 organisms occur in all organisms. Thus, 9% of EC numbers only occur in some of the organisms.
This indicates that organisms evolve by acquiring new enzymatic functionalities. Whether they stem from neofunctionalisation, horizontal gene transfer, or genesis remains to be investigated.
"""
from FEV_KEGG.Graph.SubstanceGraphs import SubstanceReactionGraph, SubstanceGeneGraph, SubstanceEcGraph
import FEV_KEGG.KEGG.Organism


if __name__ == '__main__':
    
    #- Get all metabolic pathways of all E. coli organisms from KEGG.
    eColiOrganisms = FEV_KEGG.KEGG.Organism.Group(searchString = 'Escherichia coli K-12').organisms
    
    #- For each organism, combine all pathways to the metabolic network, by UNION operation.
    organismEcGraphs = []
    for organism in eColiOrganisms:
        organismPathways = organism.getMetabolicPathways()
        organismSubstanceReactionGraph = SubstanceReactionGraph.fromPathway(organismPathways)
    
        #- Convert this metabolic network into a substance-ecNumber graph.
        organismSubstanceGeneGraph = SubstanceGeneGraph.fromSubstanceReactionGraph(organismSubstanceReactionGraph)
        organismSubstanceEcGraph = SubstanceEcGraph.fromSubstanceGeneGraph(organismSubstanceGeneGraph)
        
        organismEcGraphs.append(organismSubstanceEcGraph)
    
    firstGraph = organismEcGraphs.pop(0)
    
    #- Combine all organisms' networks to a single unified E. coli network, by UNION operation.
    unifiedEcGraph = firstGraph
    unifiedEcGraph = unifiedEcGraph.union(organismEcGraphs)
    
    #- Combine all organisms' networks to a consensus network, by INTERSECT operation, leaving only substances and EC numbers that occur in all organisms.
    intersectedEcGraph = firstGraph
    intersectedEcGraph = intersectedEcGraph.intersection(organismEcGraphs)
    
    
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