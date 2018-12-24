"""
Question
--------
Which Genes cause the occurence of EC numbers in either E. coli (eco) or S. solfataricus (sso) within the Serine, Glycine, Threonine pathway (00260)?

Method
------
- Get 00260 pathways of eco and sso from KEGG.
- Convert each pathway into a substance-ecNumber graph. Using a substance-enzyme graph as intermediate. This is the incomplete metabolic network.
- Combine both species' networks to a consensus network, by INTERSECTION operation.
- Create the networks of EC numbers that only occur in 1) eco, and 2) sso.
- Print the genes associated with those EC numbers, using the intermediate substance-enzyme graph.

Result
------

::

    EC numbers occuring only in eco (associated genes):
    18 results
    1.1.1.- (eco:b1539)
    1.1.1.103 (eco:b3616)
    1.1.1.215 (eco:b3553)
    1.1.1.3 (eco:b3940, eco:b0002)
    1.1.1.381 (eco:b1539)
    1.1.1.79 (eco:b1033, eco:b3553)
    1.1.1.81 (eco:b1033, eco:b3553)
    1.1.99.1 (eco:b0311)
    1.2.1.8 (eco:b0312)
    2.3.1.29 (eco:b3617)
    2.6.1.52 (eco:b0907)
    2.7.2.4 (eco:b3940, eco:b0002, eco:b4024)
    2.7.8.8 (eco:b2585)
    3.1.3.3 (eco:b4388)
    4.1.2.48 (eco:b0870)
    4.3.1.17 (eco:b2797, eco:b4471, eco:b1814)
    4.3.1.18 (eco:b2366)
    5.4.2.11 (eco:b0755)
    EC numbers occuring only in sso (associated genes):
    1 results
    1.5.3.1 (sso:SSO0186, sso:SSO0187)

Conclusion
----------
Many unique enzyme functions are realised by multiple genes. Therefore, multiple genes have to be considered in researching enzyme function evolution.
It is therefore best to usually only work on EC graphs, because redundant enzymes are not considered in the context of this work.
"""
from FEV_KEGG.Graph.SubstanceGraphs import SubstanceReactionGraph, SubstanceGeneGraph, SubstanceEcGraph, SubstanceEnzymeGraph
import FEV_KEGG.KEGG.Organism


if __name__ == '__main__':
    
    #- Get 00260 pathways of eco and sso from KEGG.
    eco = FEV_KEGG.KEGG.Organism.Organism('eco')
    sso = FEV_KEGG.KEGG.Organism.Organism('sso')
    
    eco00260 = eco.getPathway('00260')
    sso00260 = sso.getPathway('00260')
    
    #- Convert each pathway into a substance-ecNumber graph. This is the incomplete metabolic network.
    ecoReactionGraph = SubstanceReactionGraph.fromPathway(eco00260)
    ssoReactionGraph = SubstanceReactionGraph.fromPathway(sso00260)
    
    ecoGeneGraph = SubstanceGeneGraph.fromSubstanceReactionGraph(ecoReactionGraph)
    ssoGeneGraph = SubstanceGeneGraph.fromSubstanceReactionGraph(ssoReactionGraph)
    
    ecoEnzymeGraph = SubstanceEnzymeGraph.fromSubstanceGeneGraph(ecoGeneGraph)
    ssoEnzymeGraph = SubstanceEnzymeGraph.fromSubstanceGeneGraph(ssoGeneGraph)
    
    ecoEcGraph = SubstanceEcGraph.fromSubstanceEnzymeGraph(ecoEnzymeGraph)
    ssoEcGraph = SubstanceEcGraph.fromSubstanceEnzymeGraph(ssoEnzymeGraph)
    
    #- Combine both species' networks to a consensus network, by INTERSECTION operation.
    intersectionEcGraph = ecoEcGraph.intersection(ssoEcGraph)
    
    #- Create the networks of EC numbers that only occur in 1) eco, and 2) sso.
    onlyEcoEcGraph = ecoEcGraph.difference(intersectionEcGraph, subtractNodes = False)
    onlySsoEcGraph = ssoEcGraph.difference(intersectionEcGraph, subtractNodes = False)
    
    #- Print those EC numbers.
    print('EC numbers occuring only in eco (associated genes):')
    output = []
    for ecNumber in onlyEcoEcGraph.getECs():
        geneIDs = ecoEnzymeGraph.getGeneIDsForEcNumber(ecNumber)
        output.append(ecNumber.__str__() + ' (' + ', '.join(str(geneID) for geneID in geneIDs) + ')')
    
    output.sort()
    print(str(len(output)) + ' results')
    for line in output:
        print(line)
    
    print('EC numbers occuring only in sso (associated genes):')
    output = []
    for ecNumber in onlySsoEcGraph.getECs():
        geneIDs = ssoEnzymeGraph.getGeneIDsForEcNumber(ecNumber)
        output.append(ecNumber.__str__() + ' (' + ', '.join(str(geneID) for geneID in geneIDs) + ')')
    
    output.sort()
    print(str(len(output)) + ' results')
    for line in output:
        print(line)