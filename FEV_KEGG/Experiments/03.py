"""
Question
--------
Which EC edges (substrate -> EC number -> product) are present in eco00260, but not in eco01100?

Method
------
- Download pathway description as KGML.
- Convert to substance-reaction graph.
- Convert to substance-gene graph.
- Convert to substance-ec graph
- Calculate difference (eco00260 - eco01100). Do not subtract nodes, to keep all surviving edges.
- Print EC number edges (substrate -> EC number -> product).

Result
------

::

    19 results
    C00022 -> 4.3.1.18 -> C00740
    C00037 -> 1.4.4.2 -> C00011
    C00037 -> 2.3.1.29 -> C03508
    C00065 -> 2.7.8.8 -> C02737
    C00101 -> 2.1.2.10 -> C00014
    C00101 -> 2.1.2.10 -> C02972
    C00168 -> 1.1.1.215 -> C00258
    C00188 -> 1.1.1.103 -> C03508
    C00258 -> 1.1.1.215 -> C00168
    C00740 -> 4.3.1.18 -> C00022
    C01242 -> 2.1.2.10 -> C00014
    C01242 -> 2.1.2.10 -> C02972
    C01888 -> 1.4.3.21 -> C00546
    C02051 -> 1.4.4.2 -> C00011
    C02051 -> 1.4.4.2 -> C01242
    C02972 -> 1.8.1.4 -> C02051
    C03508 -> 2.3.1.29 -> C00037
    C05519 -> 1.1.1.- -> C03508
    C05519 -> 1.1.1.381 -> C03508


Conclusion
----------
Global map 01100 does not necessarily contain all EC number edges for an organism.
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
    
    #- Calculate difference (eco00260 - eco01100). Do not subtract nodes, to keep all surviving edges.
    difference_ecGraph = eco00260_ecGraph.difference(eco01100_ecGraph, subtractNodes = False)
    
    #- Print EC number edges (substrate -> EC number -> product).
    output = []
    edges = difference_ecGraph.getEdges()
    for edge in edges:
        substrate, product, ecNumber = edge
        output.append(substrate.__str__() + ' -> ' + ecNumber.__str__() + ' -> ' + product.__str__())
    
    output.sort()
    print(str(len(output)) + ' results')
    for line in output:
        print(line)