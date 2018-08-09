"""
Context
-------
Similar to :mod:`21`, where both group and super group were calculated with the same majority threshold. The super group seemed even more diverse than the group.
Maybe the diversity can be leveled out when using different majority thresholds for group and super group.

Question
--------
Taking a complete set of organisms, forming a group (Brachyspirales, Order) and a super group (Spirochaetes, Class) of organisms, 
how many functions occur in the majority network of the group, but not in the majority network of the super group?
For a fixed and high majority threshold for the smaller group, but varying majority thresholds for the bigger super group.

Method
------
Similar to :mod:`21`, only fixed majority threshold for group.
- Create a group of all organisms of Order Brachyspirales to be found in KEGG.
- Create a super group of all organisms of Class Spirochaetes to be found in KEGG, including the same organsims as the group of Order.
- Calculate majority substance-ec graph for group. Leaving only EC numbers which occur in all organisms of the group.
- REPEAT for varying majority-percentages:
-     Calculate majority substance-ec graph for super group. Leaving only EC numbers which occur in all organisms of the super group.
-     Calculate the difference of the two sets of majority EC numbers, leaving only the EC numbers which occur in Order majority, but not in Class majority.
-     Print the amount of EC numbers and their percentage of all EC numbers in Order, ie. how many of the EC numbers in Order do not exist in Class majority.


Result
------

::

    302/328 -> 92.1% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 100% majority for super group and 80% majority for group
    280/328 -> 85.4% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 90% majority for super group and 80% majority for group
    268/328 -> 81.7% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 80% majority for super group and 80% majority for group
    258/328 -> 78.7% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 70% majority for super group and 80% majority for group
    245/328 -> 74.7% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 60% majority for super group and 80% majority for group
    220/328 -> 67.1% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 50% majority for super group and 80% majority for group
    153/328 -> 46.6% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 40% majority for super group and 80% majority for group
    108/328 -> 32.9% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 30% majority for super group and 80% majority for group
     43/328 -> 13.1% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 20% majority for super group and 80% majority for group
      1/328 ->  0.3% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 10% majority for super group and 80% majority for group


Conclusion
----------
Percentages above 80%/80% grow a litte bit larger, below a bit smaller. Generally not much change, except for the 10% area.
This low area is probably heavily influneced by the fact that the Brachyspirales group is also included in the super group.
Maybe the group should be excluded from the super group to get less extreme results. The question answered would then be 'EC numbers new to group, compared to ALL OTHER groups in super group'.
But it could just boost the lower threshold areas.
"""
import FEV_KEGG.KEGG.Organism
from FEV_KEGG.Evolution.Taxonomy import NCBI
from FEV_KEGG.Statistics.Percent import getPercentSentence



groupMajority = 80

def groupEcGraph(taxonomy):
    #- Create a group of example organisms of Order.
    nodes = taxonomy.getOrganismNodesByPath('Bacteria/Spirochaetes/Brachyspirales')
    brachyspirales_organisms = FEV_KEGG.KEGG.Organism.Group( taxonomy.getOrganismAbbreviations( nodes ) )
    
    return brachyspirales_organisms.majorityEcGraph(groupMajority)

def supergroupEcGraph(taxonomy):
    #- Create a super group of example organisms of Class, including the same organsims as the group of Order.
    nodes = taxonomy.getOrganismNodesByPath('Bacteria/Spirochaetes')
    spirochaetes_organisms = FEV_KEGG.KEGG.Organism.Group( taxonomy.getOrganismAbbreviations( nodes ) )
    
    return spirochaetes_organisms.majorityEcGraph(majorityPercentage)
    

if __name__ == '__main__':
    
    taxonomy = NCBI.getTaxonomy()
    
    output = []
    
    #- Calculate majority substance-ec graph for group. Leaving only EC numbers which occur in all organisms of the group.
    brachyspirales_EC_graph = groupEcGraph(taxonomy)
    
    #- REPEAT for varying majority-percentages:
    for majorityPercentage in range(100, 0, -10):
    
        #-     Calculate majority substance-ec graph for super group. Leaving only EC numbers which occur in all organisms of the super group.
        spirochaetes_EC_graph = supergroupEcGraph(taxonomy)
        
        #-     Calculate the difference of the two sets of majority EC numbers, leaving only the EC numbers which occur in Order majority, but not in Class majority.
        brachyspirales_EC_set = brachyspirales_EC_graph.getECs()
        spirochaetes_EC_set = spirochaetes_EC_graph.getECs()
        
        only_brachyspirales_EC_set = brachyspirales_EC_set.difference(spirochaetes_EC_set)
        
        #-     Print the amount of EC numbers and their percentage of all EC numbers in Order, ie. how many of the EC numbers in Order do not exist in Class majority.
        output.append( getPercentSentence(len(only_brachyspirales_EC_set), len(brachyspirales_EC_set)) + ' of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ ' + str(majorityPercentage) + '% majority for super group and 80% majority for group' )
    
    
    for line in output:
        print(line)
    