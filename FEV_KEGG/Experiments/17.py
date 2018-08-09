"""
Question
--------
Assuming eco and sso had a common ancestor, and further assuming eco has evolved into a more versatile metabolism, while sso stayed closer to the common ancestor. Surrogate the common ancestor with sso.
Which functions - defined by a new EC number - occuring in eco, but not in sso, may have originated from a gene duplication with a following neofunctionalisation? Using triangle method. 

Method
------
Similar to :mod:`16`:
- Get all known pathways of eco and sso from KEGG.
- Convert each pathway into a substance-ecNumber graph. Using a substance-enzyme graph as intermediate. This is the (a little incomplete) metabolic network.
- Remove multifunctional enzymes, meaning enzymes associated with more than one EC number. Helps to reduce false gene duplications.
- Combine both species' networks to a consensus network, by INTERSECTION operation.
- Create the networks of EC numbers that only occur in eco.
New steps, triangle condition:
- For each eco-only EC number, gather any eco gene which encodes an enzyme associated with this EC number.
- For each such gene, find homologs within eco (paralogs, threshold=200).
- If there exists at least one paralog with a different function, find homologs in sso (orthologs, threshold=200) for all paralogs.

- If two paralogs share the same ortholog, there has probably been a gene duplication.
  The one ortholog and two paralogs can be visualised as a triangle, nodes being genes, edges being a connection of type ortholog/paralog. 

- If one of the paralogs shares the same EC number with their ortholog, but the other does not, there has probably been a gene duplication with neofunctionalisation.

- Report all such probable gene duplications, in the following format:
  [orthologous genes, all of them are possible ancestors]
  => paralog with the same function as the ancestor (EC number of ancestor and paralog) + new gene (new function)
  This may report several genes as the source (ancestor) of the gene duplication, because when there are multiple matching orthologs it is unknown which of the orthologs was duplicated. May also report other neofunctionalised genes as paralogs, because a single group of paralogs can have multiple neofunctionalisations.

- Report percentage gene duplication and neofunctionalisations of all new functions.

Result
------

::

    63 results
    [sso:SSO0004]
      => eco:b3752 (2.7.1.15) + eco:b1723 (2.7.1.11)
    
    [sso:SSO0090]
      => eco:b2890 (6.1.1.6) + eco:b0930 (6.1.1.22)
    
    [sso:SSO0090]
      => eco:b4129 (6.1.1.6) + eco:b0930 (6.1.1.22)
    
    [sso:SSO0093]
      => eco:b2400 (6.1.1.17) + eco:b0680 (6.1.1.18)
    
    [sso:SSO0166]
      => eco:b0759 (5.1.3.2) + eco:b3619 (5.1.3.20)
    
    [sso:SSO0173]
      => eco:b1866 (6.1.1.12) + eco:b0930 (6.1.1.22)
    
    [sso:SSO0182]
      => eco:b0154 (5.4.3.8) + eco:b0774 (2.6.1.62)
    
    [sso:SSO0182]
      => eco:b0154 (5.4.3.8) + eco:b1748 (2.6.1.81)
    
    [sso:SSO0182]
      => eco:b0154 (5.4.3.8) + eco:b3073 (2.6.1.82)
    
    [sso:SSO0299]
      => eco:b2465 (2.2.1.1) + eco:b0420 (2.2.1.7)
    
    [sso:SSO0299]
      => eco:b2935 (2.2.1.1) + eco:b0420 (2.2.1.7)
    
    [sso:SSO0306]
      => eco:b3281 (1.1.1.25) + eco:b1692 (1.1.1.282)
    
    [sso:SSO0366, sso:SSO2554, sso:SSO2440]
      => eco:b3870 (6.3.1.2) + eco:b1297 (6.3.1.11)
    
    [sso:SSO0369, sso:SSO3064]
      => eco:b1805 (6.2.1.3) + eco:b0586 (6.3.2.14)
    
    [sso:SSO0369, sso:SSO3064]
      => eco:b1805 (6.2.1.3) + eco:b2260 (6.2.1.26)
    
    [sso:SSO0579]
      => eco:b0077 (2.2.1.6) + eco:b0507 (4.1.1.47)
    
    [sso:SSO0579]
      => eco:b0077 (2.2.1.6) + eco:b2373 (4.1.1.8)
    
    [sso:SSO0579]
      => eco:b3671 (2.2.1.6) + eco:b0507 (4.1.1.47)
    
    [sso:SSO0579]
      => eco:b3671 (2.2.1.6) + eco:b2264 (2.2.1.9)
    
    [sso:SSO0579]
      => eco:b3671 (2.2.1.6) + eco:b2373 (4.1.1.8)
    
    [sso:SSO0810]
      => eco:b2028 (1.1.1.22) + eco:b3787 (1.1.1.336)
    
    [sso:SSO0893]
      => eco:b1264 (4.1.3.27) + eco:b0593 (5.4.4.2)
    
    [sso:SSO0893]
      => eco:b1264 (4.1.3.27) + eco:b1812 (2.6.1.85)
    
    [sso:SSO0893]
      => eco:b1264 (4.1.3.27) + eco:b2265 (5.4.4.2)
    
    [sso:SSO0997]
      => eco:b2574 (1.4.3.16) + eco:b4154 (1.3.5.4)
    
    [sso:SSO1077]
      => eco:b1611 (4.2.1.2) + eco:b4139 (4.3.1.1)
    
    [sso:SSO1123, sso:SSO1565]
      => eco:b0116 (1.8.1.4) + eco:b3365 (1.7.1.15)
    
    [sso:SSO1530, sso:SSO1529]
      => eco:b0115 (2.3.1.12) + eco:b0727 (2.3.1.61)
    
    [sso:SSO1565, sso:SSO1123, sso:SSO1524, sso:SSO2559]
      => eco:b0116 (1.8.1.4) + eco:b3500 (1.8.1.7)
    
    [sso:SSO1565, sso:SSO1123, sso:SSO1524, sso:SSO2559]
      => eco:b0116 (1.8.1.4) + eco:b3962 (1.6.1.1)
    
    [sso:SSO1781, sso:SSO0830]
      => eco:b2041 (4.2.1.46) + eco:b2052 (1.1.1.271)
    
    [sso:SSO1781, sso:SSO0830]
      => eco:b2041 (4.2.1.46) + eco:b2053 (4.2.1.47)
    
    [sso:SSO1781, sso:SSO0830]
      => eco:b3788 (4.2.1.46) + eco:b2052 (1.1.1.271)
    
    [sso:SSO1781, sso:SSO0830]
      => eco:b3788 (4.2.1.46) + eco:b2053 (4.2.1.47)
    
    [sso:SSO2133, sso:SSO1600]
      => eco:b3926 (2.7.1.30) + eco:b2803 (2.7.1.51)
    
    [sso:SSO2133, sso:SSO1600]
      => eco:b3926 (2.7.1.30) + eco:b3564 (2.7.1.17)
    
    [sso:SSO2133, sso:SSO1600]
      => eco:b3926 (2.7.1.30) + eco:b3580 (2.7.1.53)
    
    [sso:SSO2289, sso:SSO2276, sso:SSO3114, sso:SSO0975, sso:SSO2500]
      => eco:b1093 (1.1.1.100) + eco:b2842 (1.1.1.127)
    
    [sso:SSO2289, sso:SSO2276, sso:SSO3114, sso:SSO0975]
      => eco:b1093 (1.1.1.100) + eco:b2705 (1.1.1.140)
    
    [sso:SSO2289, sso:SSO2276, sso:SSO3114, sso:SSO2500]
      => eco:b1093 (1.1.1.100) + eco:b0596 (1.3.1.28)
    
    [sso:SSO2368]
      => eco:b3939 (2.5.1.48) + eco:b3008 (4.4.1.8)
    
    [sso:SSO2451]
      => eco:b2507 (6.3.5.2) + eco:b3360 (2.6.1.85)
    
    [sso:SSO2536, sso:SSO1220, sso:SSO1646, sso:SSO2494, sso:SSO2878]
      => eco:b1478 (1.1.1.1) + eco:b1580 (1.1.1.380)
    
    [sso:SSO2536, sso:SSO1220, sso:SSO1646, sso:SSO2494, sso:SSO2878]
      => eco:b1478 (1.1.1.1) + eco:b2091 (1.1.1.251)
    
    [sso:SSO2536, sso:SSO1220, sso:SSO1646, sso:SSO2494, sso:SSO2878]
      => eco:b1478 (1.1.1.1) + eco:b3616 (1.1.1.103)
    
    [sso:SSO2589]
      => eco:b0720 (2.3.3.1) + eco:b0333 (2.3.3.5)
    
    [sso:SSO2665]
      => eco:b4478 (4.2.1.6) + eco:b1581 (4.2.1.8)
    
    [sso:SSO2665]
      => eco:b4478 (4.2.1.6) + eco:b2247 (4.2.1.90)
    
    [sso:SSO2863, sso:SSO2059, sso:SSO2041, sso:SSO1903, sso:SSO2216, sso:SSO1342, sso:SSO2070, sso:SSO1340]
      => eco:b4069 (6.2.1.1) + eco:b0335 (6.2.1.17)
    
    [sso:SSO2863, sso:SSO2059, sso:SSO2041, sso:SSO1903, sso:SSO2216, sso:SSO2070]
      => eco:b4069 (6.2.1.1) + eco:b0586 (6.3.2.14)
    
    [sso:SSO2869]
      => eco:b1479 (1.1.1.38) + eco:b2463 (1.1.1.40)
    
    [sso:SSO2962, sso:SSO1012]
      => eco:b0871 (1.2.5.1) + eco:b0507 (4.1.1.47)
    
    [sso:SSO2962, sso:SSO1012]
      => eco:b0871 (1.2.5.1) + eco:b2373 (4.1.1.8)
    
    [sso:SSO3035, sso:SSO3072, sso:SSO2274]
      => eco:b0268 (4.3.3.7) + eco:b3225 (4.1.3.3)
    
    [sso:SSO3035, sso:SSO3072, sso:SSO2274]
      => eco:b2478 (4.3.3.7) + eco:b3225 (4.1.3.3)
    
    [sso:SSO3035, sso:SSO3072, sso:SSO2274]
      => eco:b4298 (4.3.3.7) + eco:b3225 (4.1.3.3)
    
    [sso:SSO3036]
      => eco:b1617 (3.2.1.31) + eco:b0344 (3.2.1.23)
    
    [sso:SSO3036]
      => eco:b1617 (3.2.1.31) + eco:b3076 (3.2.1.23)
    
    [sso:SSO3064]
      => eco:b1805 (6.2.1.3) + eco:b0335 (6.2.1.17)
    
    [sso:SSO3107]
      => eco:b3771 (4.2.1.9) + eco:b1851 (4.2.1.12)
    
    [sso:SSO3211, sso:SSO2727]
      => eco:b1302 (2.6.1.19) + eco:b0774 (2.6.1.62)
    
    [sso:SSO3211, sso:SSO2727]
      => eco:b1302 (2.6.1.19) + eco:b1748 (2.6.1.81)
    
    [sso:SSO3211, sso:SSO2727]
      => eco:b1302 (2.6.1.19) + eco:b3073 (2.6.1.82)
    
    42/413 -> 10.2 percent gene duplicated and neofunctionalised of all new functions.

Conclusion
----------
Many gene duplication events can be identified using the triangle method. However,
a) the triangle assumes that sso is the ancestor of eco, which it is not. Using groups of species could alleviate this problem: form consensus between Gammaproteobacteria, form consensus between all Proteobacteria, replace sso with Proteobacteria and eco with Gammaproteobacteria. This approach resembles a simple simulation of reconstructing a common ancestor, but without the need to reconstruct ancestral genes.
b) searching for paralogs with differing function only, might suffice for finding gene duplications. Because if there is a paralog with a new function, there surely has been a gene duplication somewhere in the tree of life. Nevertheless, a change in function for the second paralog does not tell us that the gene duplication happened in the examined part of the tree of life, it merely states that a new function was developed within this part of the tree. But the interesting event here is the neofunctionalisation, not the gene duplication, which renders this method valid.
"""
from FEV_KEGG.Graph.SubstanceGraphs import SubstanceReactionGraph, SubstanceGeneGraph, SubstanceEcGraph, SubstanceEnzymeGraph
from FEV_KEGG.KEGG import Database
import FEV_KEGG.KEGG.Organism


if __name__ == '__main__':
    
    #- Get all known pathways of eco and sso from KEGG.
    eco = FEV_KEGG.KEGG.Organism.Organism('eco')
    sso = FEV_KEGG.KEGG.Organism.Organism('sso')
    
    ecoPathways = eco.getMetabolicPathways()
    ssoPathways = sso.getMetabolicPathways()
    
    #- Convert each pathway into a substance-ecNumber graph. This is the incomplete metabolic network.
    ecoReactionGraph = SubstanceReactionGraph.fromPathway(ecoPathways)
    ssoReactionGraph = SubstanceReactionGraph.fromPathway(ssoPathways)
    
    ecoGeneGraph = SubstanceGeneGraph.fromSubstanceReactionGraph(ecoReactionGraph)
    ssoGeneGraph = SubstanceGeneGraph.fromSubstanceReactionGraph(ssoReactionGraph)
    
    ecoEnzymeGraph = SubstanceEnzymeGraph.fromSubstanceGeneGraph(ecoGeneGraph)
    ssoEnzymeGraph = SubstanceEnzymeGraph.fromSubstanceGeneGraph(ssoGeneGraph)
    
    #- Remove multifunctional enzymes, meaning enzymes associated with more than one EC number. Helps to reduce false gene duplications.
    ecoEnzymeGraph.removeMultifunctionalEnzymes()
    ssoEnzymeGraph.removeMultifunctionalEnzymes()
    
    ecoEcGraph = SubstanceEcGraph.fromSubstanceEnzymeGraph(ecoEnzymeGraph)
    ssoEcGraph = SubstanceEcGraph.fromSubstanceEnzymeGraph(ssoEnzymeGraph)
    
    #- Combine both species' networks to a consensus network, by INTERSECTION operation.
    intersectionEcGraph = ecoEcGraph.intersection(ssoEcGraph)
    
    #- Create the networks of EC numbers that only occur in eco.
    onlyEcoEcGraph = ecoEcGraph.difference(intersectionEcGraph, subtractNodes = False)
    
    
    
    # New steps:
    output = []
    
    #- For each eco-only EC number, gather any eco gene which encodes an enzyme associated with this EC number.
    new_ec_numbers = onlyEcoEcGraph.getECs()
    GD_ec_numbers = set()
    for ecNumber in new_ec_numbers:
        geneIDs = ecoEnzymeGraph.getGeneIDsForEcNumber(ecNumber)
        
        #- For each such gene, find homologs within eco (paralogs).
        for geneID in geneIDs:
            paralogs = Database.getParalogsOnlyGeneID(geneID)
                        
            #- If there exists at least one paralog with a different function, find homologs in sso (orthologs) for all such paralogs.
            paralogs_different_function = []
            
            for paralog in paralogs:
                paralog_enzyme = ecoEnzymeGraph.getEnzymeForGeneID(paralog)
                
                if paralog_enzyme is None: # can happen, when the paralogous gene encodes a non-enzymatic protein
                    continue
                elif ecNumber not in paralog_enzyme.ecNumbers:
                    paralogs_different_function.append([paralog, paralog_enzyme.ecNumbers])
                    
            if len( paralogs_different_function ) == 0:
                continue
            
            geneID_orthologs = Database.getOrthologsOnlyGeneID(geneID, sso)
            
            for entry in paralogs_different_function:
                paralog, paralog_ecNumbers = entry
                paralog_orthologs = Database.getOrthologsOnlyGeneID(paralog, sso)
                
                #- If two paralogs share the same ortholog, there has probably been a gene duplication.
                shared_orthologs = paralog_orthologs.intersection(geneID_orthologs)
                possible_GD_orthologs = []
                
                for shared_ortholog in shared_orthologs:
                    
                    #- If one of the paralogs shares the same EC number with their ortholog, but the other does not, there has probably been a gene duplication with neofunctionalisation.
                    shared_ortholog_enzyme = ssoEnzymeGraph.getEnzymeForGeneID(shared_ortholog)
                    if shared_ortholog_enzyme is None: # can happen, when the orthologous gene encodes a non-enzymatic protein
                        continue
                    shared_ortholog_ecNumbers = shared_ortholog_enzyme.ecNumbers
                    if shared_ortholog_ecNumbers == paralog_ecNumbers:
                        possible_GD_orthologs.append(shared_ortholog)
                
                if len( possible_GD_orthologs ) > 0:
                    output.append("[" + ", ".join([x.__str__() for x in possible_GD_orthologs]) + "]\n  => " + paralog.__str__() + " (" + ", ".join([x.__str__() for x in paralog_ecNumbers]) + ") + " + geneID.__str__() + " (" + ecNumber.__str__() + ")\n")
                    
                    #- Report percentage gene duplication and neofunctionalisations of all new functions.
                    GD_ec_numbers.add(ecNumber)
    
    output.sort()
    print(str(len(output)) + ' results')
    for line in output:
        print(line)
    
    print("\n" + str(len(GD_ec_numbers)) + "/" + str(len(new_ec_numbers)) + " -> %2.1f percent gene duplicated and neofunctionalised of all new functions." % (len(GD_ec_numbers)/len(new_ec_numbers)*100))
    