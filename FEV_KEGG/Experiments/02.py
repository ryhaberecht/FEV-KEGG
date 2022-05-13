"""
Question
--------
Which enzyme edges (substrate -> enzyme -> product) are present in eco00260, but not in eco01100?

Method
------
- Download pathway definition as KGML.
- Convert to substance-reaction graph.
- Convert to substance-gene graph.
- Calculate difference (eco00260 - eco01100). Do not subtract nodes, to keep all surviving edges.
- Convert difference graph to substance-enzyme graph. Not necessary earlier, because each unique gene object ist mapped to the same unique enzyme object.
- Print enzyme edges (substrate -> gene name, enzyme name, all associated EC numbers, enzyme long name -> product).

Result
------

::

    18 results
    C00022 -> b2366 dsdA [4.3.1.18] D-serine dehydratase -> C00740
    C00037 -> b2903 gcvP [1.4.4.2] glycine dehydrogenase -> C00011
    C00037 -> b3617 kbl [2.3.1.29] glycine C-acetyltransferase -> C03508
    C00065 -> b2585 pssA [2.7.8.8] CDP-diacylglycerol---serine O-phosphatidyltransferase -> C02737
    C00101 -> b2905 gcvT [2.1.2.10] aminomethyltransferase -> C00014
    C00101 -> b2905 gcvT [2.1.2.10] aminomethyltransferase -> C02972
    C00168 -> b3553 ghrB [1.1.1.215, 1.1.1.81, 1.1.1.79] glyoxylate/hydroxypyruvate/2-ketogluconate reductase -> C00258
    C00188 -> b3616 tdh [1.1.1.103] threonine 3-dehydrogenase -> C03508
    C00258 -> b3553 ghrB [1.1.1.215, 1.1.1.81, 1.1.1.79] glyoxylate/hydroxypyruvate/2-ketogluconate reductase -> C00168
    C00740 -> b2366 dsdA [4.3.1.18] D-serine dehydratase -> C00022
    C01242 -> b2905 gcvT [2.1.2.10] aminomethyltransferase -> C00014
    C01242 -> b2905 gcvT [2.1.2.10] aminomethyltransferase -> C02972
    C01888 -> b1386 tynA [1.4.3.21] primary-amine oxidase -> C00546
    C02051 -> b2903 gcvP [1.4.4.2] glycine dehydrogenase -> C00011
    C02051 -> b2903 gcvP [1.4.4.2] glycine dehydrogenase -> C01242
    C02972 -> b0116 lpd [1.8.1.4] dihydrolipoamide dehydrogenase -> C02051
    C03508 -> b3617 kbl [2.3.1.29] glycine C-acetyltransferase -> C00037
    C05519 -> b1539 ydfG [1.1.1.-, 1.1.1.381] 3-hydroxy acid dehydrogenase / malonic semialdehyde reductase -> C03508


Conclusion
----------
Global map 01100 does not necessarily contain all enzyme edges for an organism.
"""
from FEV_KEGG.Graph.SubstanceGraphs import SubstanceReactionGraph, SubstanceGeneGraph, SubstanceEnzymeGraph
from FEV_KEGG.KEGG.Organism import Organism


if __name__ == '__main__':
    
    #- Download pathway definition as KGML.
    eco = Organism('eco')
    
    eco00260 = eco.getPathway('00260')
    eco01100 = eco.getPathway('01100')
    
    #- Convert to substance-reaction graph.
    eco00260_reactionGraph = SubstanceReactionGraph.fromPathway(eco00260)
    eco01100_reactionGraph = SubstanceReactionGraph.fromPathway(eco01100)
    
    #- Convert to substance-gene graph
    eco00260_geneGraph = SubstanceGeneGraph.fromSubstanceReactionGraph(eco00260_reactionGraph)
    eco01100_geneGraph = SubstanceGeneGraph.fromSubstanceReactionGraph(eco01100_reactionGraph)
    
    #- Calculate difference (eco00260 - eco01100). Do not subtract nodes, to keep all surviving edges.
    difference_geneGraph = eco00260_geneGraph.difference(eco01100_geneGraph, subtractNodes = False)
    
    #- Convert difference graph to substance-enzyme graph.
    difference_enzymeGraph = SubstanceEnzymeGraph.fromSubstanceGeneGraph(difference_geneGraph)
    
    #- Print enzymes (substrate -> gene name, enzyme name, all associated EC numbers, enzyme long name -> product).
    output = []
    edges = difference_enzymeGraph.getEdges()
    for edge in edges:
        substrate, product, enzyme = edge
        output.append(substrate.__str__() + ' -> ' + enzyme.geneName + ' ' + enzyme.name + ' [' + enzyme.getEcNumbersString() + '] ' + enzyme.description + ' -> ' + product.__str__())
    
    output.sort()
    print(str(len(output)) + ' results')
    for line in output:
        print(line)