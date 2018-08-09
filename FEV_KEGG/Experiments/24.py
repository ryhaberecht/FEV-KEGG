"""
Context
-------
Similar to :mod:`23`, where group and super group were calculated with a fixed and a variable majority threshold.
Maybe the group should be excluded from the super group to get less extreme results. The question answered would then be 'EC numbers new to group, compared to ALL OTHER groups in super group'.

Question
--------
Taking a complete set of organisms, forming a group (Brachyspirales, Order) and a super group (Spirochaetes, Class) of organisms, 
how many functions occur in the majority network of the group, but not in the majority network of the super group, excluding the organisms of the group?
For a fixed and high majority threshold for the smaller group, but varying majority thresholds for the bigger super group.

Method
------
Similar to :mod:`21`, only fixed majority threshold for group.
- Create a group of all organisms of Order Brachyspirales to be found in KEGG.
- Create a super group of all organisms of Class Spirochaetes to be found in KEGG, including the same organsims as the group of Order.
- Exclude organisms of group from super group.
- Calculate majority substance-ec graph for group. Leaving only EC numbers which occur in all organisms of the group.
- REPEAT for varying majority-percentages:
-     Calculate majority substance-ec graph for super group. Leaving only EC numbers which occur in all organisms of the super group.
-     Calculate the difference of the two sets of majority EC numbers, leaving only the EC numbers which occur in Order majority, but not in Class majority.
-     Print the amount of EC numbers and their percentage of all EC numbers in Order, ie. how many of the EC numbers in Order do not exist in Class majority.

Result
------

::

    300/328 -> 91.5% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 100% majority for super group and 80% majority for group
    281/328 -> 85.7% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 90% majority for super group and 80% majority for group
    271/328 -> 82.6% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 80% majority for super group and 80% majority for group
    264/328 -> 80.5% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 70% majority for super group and 80% majority for group
    254/328 -> 77.4% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 60% majority for super group and 80% majority for group
    224/328 -> 68.3% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 50% majority for super group and 80% majority for group
    203/328 -> 61.9% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 40% majority for super group and 80% majority for group
    142/328 -> 43.3% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 30% majority for super group and 80% majority for group
     85/328 -> 25.9% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 20% majority for super group and 80% majority for group
     35/328 -> 10.7% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 10% majority for super group and 80% majority for group


Conclusion
----------
Leaving out the group and comparing to the rest of the super group merely boosts the percentages for the lower thresholds. For higher thresholds, the group's redundancy does not have a noticeable effect.
If the group were bigger, though, the effect may well be noticeable and skew the higher thresholds, too. It should always be assured, that the super group is sufficiently larger than the group.
"""
from FEV_KEGG.Evolution.Taxonomy import NCBI
import FEV_KEGG.KEGG.Organism
from FEV_KEGG.Statistics.Percent import getPercentSentence


if __name__ == '__main__':
    
    output = []
    
    taxonomy = NCBI.getTaxonomy()
    groupMajority = 80
    
    #- Create a group of example organisms of Order.
    nodesGroup = set(taxonomy.getOrganismNodesByPath('Bacteria/Spirochaetes/Brachyspirales'))
    brachyspirales_organisms = FEV_KEGG.KEGG.Organism.Group( taxonomy.getOrganismAbbreviations( nodesGroup ) )
    
    #- Create a super group of all organisms of Class Spirochaetes to be found in KEGG, including the same organsims as the group of Order.
    nodesSupergroup = set(taxonomy.getOrganismNodesByPath('Bacteria/Spirochaetes'))
    #- Exclude organisms of group from super group.
    nodesSupergroup = nodesSupergroup.difference(nodesGroup)
    spirochaetes_organisms = FEV_KEGG.KEGG.Organism.Group( taxonomy.getOrganismAbbreviations( nodesSupergroup ) )
    
    #- Calculate majority substance-ec graph for group. Leaving only EC numbers which occur in all organisms of the group.
    brachyspirales_EC_graph = brachyspirales_organisms.majorityEcGraph(groupMajority)
    
    
    #- REPEAT for varying majority-percentages:
    for majorityPercentage in range(100, 0, -10):
    
        #-     Calculate majority substance-ec graph for super group. Leaving only EC numbers which occur in all organisms of the super group.
        spirochaetes_EC_graph = spirochaetes_organisms.majorityEcGraph(majorityPercentage)#FIXME nicht ganz richtig, laufne lassen udn mit oben vergleichen!!!
        
        #-     Calculate the difference of the two sets of majority EC numbers, leaving only the EC numbers which occur in Order majority, but not in Class majority.
        brachyspirales_EC_set = brachyspirales_EC_graph.getECs()
        spirochaetes_EC_set = spirochaetes_EC_graph.getECs()
        
        only_brachyspirales_EC_set = brachyspirales_EC_set.difference(spirochaetes_EC_set)
        
        #-     Print the amount of EC numbers and their percentage of all EC numbers in Order, ie. how many of the EC numbers in Order do not exist in Class majority.
        output.append( getPercentSentence(len(only_brachyspirales_EC_set), len(brachyspirales_EC_set)) + ' of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ ' + str(majorityPercentage) + '% majority for super group and 80% majority for group' )
    
    
    for line in output:
        print(line)
    