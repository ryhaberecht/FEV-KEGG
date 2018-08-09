"""
Question
--------
Does a simple UNION operation on node sets and edge sets suffice for combining several pathways to a global metabolic network?

Method
------
- Download pathway descriptions as KGML -> eco00260 and eco00270.
- Convert to substance-reaction graphs -> eco00260_graph and eco00270_graph.
- Combine the pathways via UNION operation -> union_graph.
- Subtract the individual parts of the resulting network. Firstly, without subtracting nodes, to leave any surviving edges -> subtracted_graph_keeping_nodes. Secondly, with subtracting nodes -> subtracted_graph_deleting_nodes.
- Print any edge that remains in subtracted_graph_keeping_nodes. Print any node that remains in subtracted_graph_deleting_nodes.

Result
------

::

    0 results

Conclusion
----------
UNION is indeed correctly implemented and fit for the task of combining several pathways to a global network.
"""
from FEV_KEGG.Graph.SubstanceGraphs import SubstanceReactionGraph
from FEV_KEGG.KEGG.Organism import Organism


if __name__ == '__main__':
    
    #- Download pathway description as KGML.
    eco = Organism('eco')
    
    eco00260 = eco.getPathway('00260')
    eco00270 = eco.getPathway('00270')
    
    #- Convert to substance-reaction graph.
    eco00260_graph = SubstanceReactionGraph.fromPathway(eco00260)
    eco00270_graph = SubstanceReactionGraph.fromPathway(eco00270)
    
    #- Combine the pathways via intersect operation -> union_graph.
    union_graph = SubstanceReactionGraph.fromPathway([eco00260, eco00270])

    #- Subtract the individual parts of the resulting network.
    # Firstly, without subtracting nodes, to leave any surviving edges -> subtracted_graph_keeping_nodes.
    subtracted_graph_keeping_nodes = union_graph.difference(eco00260_graph, subtractNodes = False).difference(eco00270_graph, subtractNodes = False)
    
    # Secondly, with subtracting nodes -> subtracted_graph_deleting_nodes.
    subtracted_graph_deleting_nodes = union_graph.difference(eco00260_graph, subtractNodes = True).difference(eco00270_graph, subtractNodes = True)
    
    #- Print any edge that remains in subtracted_graph_keeping_nodes.
    output = []
    for reaction in subtracted_graph_keeping_nodes.getReactions():
        output.append(reaction.__str__())
    
    #Print any node that remains in subtracted_graph_deleting_nodes.
    output = []
    for substance in subtracted_graph_deleting_nodes.getNodes():
        output.append(substance.__str__())
    
    output.sort()
    print(str(len(output)) + ' results')
    for line in output:
        print(line)