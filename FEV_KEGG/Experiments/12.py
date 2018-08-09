"""
Question
--------
Which EC numbers are not present in all Escherichia coli K-12 subspecies?

Method
------
- Get all metabolic pathways of all E. coli species from KEGG.
- For each species, combine all pathways to the metabolic network, by UNION operation.
- Convert this metabolic network into a substance-ecNumber graph.
- Combine all species' networks to a single unified E. coli network, by UNION operation.
- Combine all species' networks to a consensus network, by INTERSECT operation, leaving only substances and EC numbers that occur in all species.
- Subtract the consensus network from the E. coli network, by DIFFERENCE operation, leaving only substances and EC numbers that do not occur in all species.
- Print all EC numbers that do not occur in all species.

Result
------

::

    69 results
    1.1.1.-
    1.1.1.133
    1.1.1.157
    1.1.1.251
    1.1.1.271
    1.1.1.28
    1.1.1.381
    1.1.1.57
    1.1.1.58
    1.1.1.65
    1.1.1.85
    1.1.3.15
    1.14.11.17
    1.14.13.149
    1.16.3.1
    1.17.1.9
    1.2.1.16
    1.2.1.2
    1.2.1.20
    1.2.1.39
    1.2.1.79
    1.2.1.91
    1.2.7.-
    1.2.7.1
    1.4.3.21
    1.5.1.34
    2.1.1.10
    2.1.2.10
    2.3.1.174
    2.3.1.223
    2.3.3.13
    2.3.3.5
    2.5.1.7
    2.6.1.16
    2.7.1.16
    2.7.1.200
    2.7.1.6
    2.7.7.13
    2.7.8.7
    2.8.3.-
    2.9.1.1
    3.1.2.-
    3.1.3.7
    3.2.1.14
    3.2.1.17
    3.2.1.23
    3.2.1.37
    3.3.2.12
    3.5.1.25
    3.5.4.1
    4.1.3.30
    4.1.3.39
    4.2.1.33
    4.2.1.35
    4.2.1.47
    4.2.1.79
    4.2.1.80
    4.2.1.9
    4.4.1.15
    4.6.1.1
    5.1.3.13
    5.3.1.4
    5.3.2.6
    5.3.3.18
    5.4.2.8
    5.4.99.2
    5.4.99.9
    6.2.1.17
    6.2.1.30

Conclusion
----------
Some EC numbers are not shared between subspecies.
It could make sense to ignore incomplete EC numbers, as they may represent identical reactions on identical substances and could, thus, be counted twice.
For example, 1.2.7.- might merely represent incomplete data, while the associated enzyme actually performs 1.2.7.1., causing a duplicate in the result list.
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
    
    #- Subtract the consensus network from the E. coli network, by DIFFERENCE operation, leaving only substances and EC numbers that do not occur in all species.
    differenceEcGraph = unifiedEcGraph.difference(intersectedEcGraph)
    
    #- Print all EC numbers that do not occur in all species.
    output = []
    for ecNumber in differenceEcGraph.getECs():
        output.append(ecNumber.__str__())
    
    output.sort()
    print(str(len(output)) + ' results')
    for line in output:
        print(line)