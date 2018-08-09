"""
Question
--------
Which EC numbers, independent of substrate and product, are present in eco01100, but not in a combination of all non-overview metabolic pathways?

Method
------
- Download pathway description as KGML.
- Convert to substance-reaction graph.
- Convert to substance-gene graph.
- Convert to substance-ec graph.
- Get set of EC numbers for each graph.
- Calculate difference of EC number sets.
- Print EC numbers.

Result
------

::

    8 results
    1.14.13.-
    2.1.1.222
    2.1.1.64
    2.5.1.39
    2.5.1.75
    2.7.1.19
    2.7.8.20
    4.1.1.98

Conclusion
----------
The overview map eco01100 contains genes, leading to EC numbers, that are not present in any of the single metabolism pathways.

Further research shows the reasons for such strange results:

- Some single pathway files in KGML format do not contain a "reaction tag", even though they should, according to the "entry tags" present in the same file and the pathway image.
  Example: Reaction R05000 only exists in eco011000, but not in eco00130, where it should, because eco00130 contains all compound and gene entries necessary for this reaction.
  Eco00130 simply does not contain the "reaction tag", probably due to an error in KEGG's algorithms for creating KGML files.
    
- The overview map eco01100 contains more gene entries and thus enzymes than the single pathways combined. In the above cases, this was probably not an error in eco01100, but in the respective single pathway maps.
  Although it may occur in different organisms, that 01100 falsly lists a non-existing gene, genes are still organism-specific and the attempt to download a non-existing gene would simply throw a download error and would, thus, not remain undetected.
  In reverse conclusion, all the genes (leading to the EC numbers) listed above actually do exist in eco and the error lies within the single pathway maps, not within the overview map 01100.
  Example: 2.7.8.20 is the function of gene eco:b4359, which occurs in eco01100, but not in eco00561, where it probably should occur, too, because 2.7.8.20 is part of the general map 00561.
    
- Maybe more.
    
Neiter using only 01100, nor only a combination of single pathways is perfect.
01100 has the following shortcomings:
- No directed edges, irreversible reactions are handled as reversible.

- Erroneous additional reactions, which are linked to EC numbers not present in any gene of the organism. 
  This only affects the substance-reaction graph, since the erroneous reaction will not be translated into a gene for the substance-gene graph.

- Some EC number edges are missing. The EC number itself is present, but not between all substrate/product combinations.
- Some EC numbers themselves are missing.

The combination of all single pathways has the following shortcomings:

- Some EC numbers themselves are missing. Due to missing "reaction tags" (fixable without 01100) AND due to missing gene entries (not fixable without 01100).

Solution A: Use combined single pathways, ignore missing EC numbers.
Solution B: Use combined single pathways, repair by comparison with overview map 01100, adding missing "reaction tags" and missing gene entries.
"""
from FEV_KEGG.Graph.SubstanceGraphs import SubstanceReactionGraph, SubstanceGeneGraph, SubstanceEcGraph
from FEV_KEGG.KEGG.Organism import Organism


if __name__ == '__main__':
    
    #- Download pathway description as KGML.
    eco = Organism('eco')
    
    eco01100 = eco.getPathway('01100')
    allNonOverviewPathways = eco.getMetabolicPathways(includeOverviewMaps = False)
    
    #- Convert to substance-reaction graph.
    eco01100_reactionGraph = SubstanceReactionGraph.fromPathway(eco01100)
    allNonOverviewPathways_reactionGraph = SubstanceReactionGraph.fromPathway(allNonOverviewPathways)
    
    #- Convert to substance-gene graph.
    eco01100_geneGraph = SubstanceGeneGraph.fromSubstanceReactionGraph(eco01100_reactionGraph)
    allNonOverviewPathways_geneGraph = SubstanceGeneGraph.fromSubstanceReactionGraph(allNonOverviewPathways_reactionGraph)
    
    #- Convert to substance-ec graph.
    eco01100_ecGraph = SubstanceEcGraph.fromSubstanceGeneGraph(eco01100_geneGraph)
    allNonOverviewPathways_ecGraph = SubstanceEcGraph.fromSubstanceGeneGraph(allNonOverviewPathways_geneGraph)
    
    #- Get set of EC numbers for each graph.
    eco01100_ecNumbers = eco01100_ecGraph.getECs()
    allNonOverviewPathways_ecNumbers = allNonOverviewPathways_ecGraph.getECs()
    
    #- Calculate difference of EC number sets.
    difference_ecNumbers = eco01100_ecNumbers.difference(allNonOverviewPathways_ecNumbers)
    
    #- Print EC numbers.
    output = []
    for ecNumber in difference_ecNumbers:
        output.append(ecNumber.__str__())
    
    output.sort()
    print(str(len(output)) + ' results')
    for line in output:
        print(line)