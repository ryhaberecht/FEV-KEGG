"""
Question
--------
Which EC numbers occur in either E. coli (eco) or S. solfataricus (sso) within the Serine, Glycine, Threonine pathway (00260)?

Method
------
- Get 00260 pathways of eco and sso from KEGG.
- Convert each pathway into a substance-ecNumber graph. This is the incomplete metabolic network.
- Combine both species' networks to a consensus network, by INTERSECTION operation.
- Create the networks of EC numbers that only occur in 1) eco, and 2) sso.
- Print those EC numbers.

Result
------

::

    EC numbers occuring only in eco:
    18 results
    1.1.1.-
    1.1.1.103
    1.1.1.215
    1.1.1.3
    1.1.1.381
    1.1.1.79
    1.1.1.81
    1.1.99.1
    1.2.1.8
    2.3.1.29
    2.6.1.52
    2.7.2.4
    2.7.8.8
    3.1.3.3
    4.1.2.48
    4.3.1.17
    4.3.1.18
    5.4.2.11
    EC numbers occuring only in sso:
    1 results
    1.5.3.1

Conclusion
----------
The results are consistent with the pathway image to be found on KEGG.
"""
from FEV_KEGG.Graph.SubstanceGraphs import SubstanceReactionGraph, SubstanceGeneGraph, SubstanceEcGraph
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
    
    ecoEcGraph = SubstanceEcGraph.fromSubstanceGeneGraph(ecoGeneGraph)
    ssoEcGraph = SubstanceEcGraph.fromSubstanceGeneGraph(ssoGeneGraph)
    
    #- Combine both species' networks to a consensus network, by INTERSECTION operation.
    intersectionEcGraph = ecoEcGraph.intersection(ssoEcGraph)
    
    #- Create the networks of EC numbers that only occur in 1) eco, and 2) sso.
    onlyEcoEcGraph = ecoEcGraph.difference(intersectionEcGraph, subtractNodes = False)
    onlySsoEcGraph = ssoEcGraph.difference(intersectionEcGraph, subtractNodes = False)
    
    #- Print those EC numbers.
    print('EC numbers occuring only in eco:')
    output = []
    for ecNumber in onlyEcoEcGraph.getECs():
        output.append(ecNumber.__str__())
    
    output.sort()
    print(str(len(output)) + ' results')
    for line in output:
        print(line)
    
    print('EC numbers occuring only in sso:')
    output = []
    for ecNumber in onlySsoEcGraph.getECs():
        output.append(ecNumber.__str__())
    
    output.sort()
    print(str(len(output)) + ' results')
    for line in output:
        print(line)