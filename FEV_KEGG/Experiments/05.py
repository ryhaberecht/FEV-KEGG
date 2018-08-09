"""
Question
--------
Which EC numbers, independent from substrate/product, are present in eco00260, but not in eco01100?

Method
------
- Download pathway description as KGML.
- Convert to substance-reaction graph.
- Convert to substance-gene graph.
- Convert to substance-ec graph
- Get set of EC numbers for each graph.
- Calculate difference of EC number sets.
- Print genes (EC number).

Result
------

::

    3 results
    1.1.1.103
    2.3.1.29
    4.3.1.18


Conclusion
----------
Global map 01100 does not necessarily contain all EC numbers for an organism.
"""
from FEV_KEGG.Graph.SubstanceGraphs import SubstanceReactionGraph, SubstanceGeneGraph, SubstanceEcGraph
from FEV_KEGG.KEGG.Organism import Organism


if __name__ == '__main__':
    
    #- Download pathway description as KGML.
    eco = Organism('eco')
    
    eco00260 = eco.getPathway('00260')
    eco01100 = eco.getPathway('01100')
    
    #- Convert to substance-reaction graph.
    eco00260_reactionGraph = SubstanceReactionGraph.fromPathway(eco00260)
    eco01100_reactionGraph = SubstanceReactionGraph.fromPathway(eco01100)
    
    #- Convert to substance-gene graph
    eco00260_geneGraph = SubstanceGeneGraph.fromSubstanceReactionGraph(eco00260_reactionGraph)
    eco01100_geneGraph = SubstanceGeneGraph.fromSubstanceReactionGraph(eco01100_reactionGraph)
    
    #- Convert to substance-ec graph
    eco00260_ecGraph = SubstanceEcGraph.fromSubstanceGeneGraph(eco00260_geneGraph)
    eco01100_ecGraph = SubstanceEcGraph.fromSubstanceGeneGraph(eco01100_geneGraph)
    
    #- Get set of EC numbers for each graph.
    eco00260_ecNumbers = eco00260_ecGraph.getECs()
    eco01100_ecNumbers = eco01100_ecGraph.getECs()
    
    #- Calculate difference of EC number sets.
    difference_ecNumbers = eco00260_ecNumbers.difference(eco01100_ecNumbers)
    
    #- Print genes (EC number).
    output = []
    for ecNumber in difference_ecNumbers:
        output.append(ecNumber.__str__())
    output.sort()
    print(str(len(output)) + ' results')
    for line in output:
        print(line)