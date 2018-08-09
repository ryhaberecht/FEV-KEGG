"""
Question
--------
Which EC numbers are present in all Escherichia coli K-12 subspecies?

Method
------
- Get all metabolic pathways of all E. coli species from KEGG.
- For each species, combine all pathways to the metabolic network, by UNION operation.
- Convert this metabolic network into a substance-ecNumber graph.
- Combine all species' networks to a consensus network, by INTERSECT operation, leaving only substances and EC numbers that occur in all species.
- Print all EC numbers that occur in all species.
- Print number of EC numbers that occur in each species combined, by UNION operation.

Result
------

::

    699 results
    ...shortened

Conclusion
----------
Alls subspecies of E. coli K-12 share a high number of EC numbers.
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
    
    #- Combine all species' networks to a consensus network, by INTERSECT operation, leaving only substances and EC numbers that occur in all species.
    intersectedEcGraph = firstGraph
    intersectedEcGraph = intersectedEcGraph.intersection(speciesEcGraphs)
    
    #- Print all EC numbers that occur in all species.
    output = []
    for ecNumber in intersectedEcGraph.getECs():
        output.append(ecNumber.__str__())
    
    output.sort()
    print(str(len(output)) + ' results')
    for line in output:
        print(line)