"""
Question
--------
Which reaction edges (substrate -> reaction ID -> product) are present in eco00260, but not in eco01100?

Method
------
- Download pathway description as KGML.
- Convert to substance-reaction graph.
- Calculate difference (eco00260 - eco01100). Do not subtract nodes, to keep all surviving edges.
- Print reaction edges (substrate -> reaction ID -> product).

Result
------

::

    17 results
    C00022 -> R00221 -> C00740
    C00037 -> R00371 -> C03508
    C00037 -> R03425 -> C00011
    C00065 -> R01800 -> C02737
    C00101 -> R04125 -> C00014
    C00101 -> R04125 -> C02972
    C00114 -> R01025 -> C00576
    C00188 -> R01465 -> C03508
    C00740 -> R00221 -> C00022
    C01242 -> R04125 -> C00014
    C01242 -> R04125 -> C02972
    C01888 -> R02529 -> C00546
    C02051 -> R03425 -> C00011
    C02051 -> R03425 -> C01242
    C02972 -> R03815 -> C02051
    C03508 -> R00371 -> C00037
    C05519 -> R10851 -> C03508

Conclusion
----------
Global map 01100 does not necessarily contain all reaction edges for an organism.
"""
from FEV_KEGG.Graph.SubstanceGraphs import SubstanceReactionGraph
from FEV_KEGG.KEGG.Organism import Organism


if __name__ == '__main__':
    
    #- Download pathway description as KGML.
    eco = Organism('eco')
    
    eco00260 = eco.getPathway('00260')
    eco01100 = eco.getPathway('01100')
    
    #- Convert to substance-reaction graph.
    eco00260_reactionGraph = SubstanceReactionGraph.fromPathway(eco00260)
    eco01100_reactionGraph = SubstanceReactionGraph.fromPathway(eco01100)
    
    #- Calculate difference (eco00260 - eco01100). Do not subtract nodes, to keep all surviving edges.
    difference_reactionGraph = eco00260_reactionGraph.difference(eco01100_reactionGraph, subtractNodes = False)
    
    #- Print reactions (substance -> reaction ID -> product).
    output = []
    edges = difference_reactionGraph.getEdges()
    for edge in edges:
        substrate, product, reactionID = edge
        output.append(substrate.__str__() + ' -> ' + reactionID.__str__() + ' -> ' + product.__str__())
    
    output.sort()
    print(str(len(output)) + ' results')
    for line in output:
        print(line)