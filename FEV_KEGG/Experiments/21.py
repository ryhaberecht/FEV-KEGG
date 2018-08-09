"""
Context
-------
Similar to :mod:`20`, where results seemingly showed very diverse metabolism. Possible reasons might have been a poor choice of example organisms or not enough organisms per group.
These reasons can be mitigated by this experiment, up to the amount of real-world organisms available in KEGG.
A small common 'core metabolism' is expected to be found, except for symbionts.

Question
--------
Taking a complete set of organisms, forming a group (Brachyspirales, Order) and a super group (Spirochaetes, Class) of organisms, 
how many functions occur in the majority network of the group, but not in the majority network of the super group?
For varying majority-percentages.

Method
------
Similar to :mod:`20`, only different sets, complete sets, and 'how many?' instead of 'which?' functions.
- Create a group of all organisms of Order Brachyspirales to be found in KEGG.
- Create a super group of all organisms of Class Spirochaetes to be found in KEGG, including the same organsims as the group of Order.
- REPEAT for varying majority-percentages:
-     Calculate majority substance-ec graphs for both groups. Leaving only EC numbers which occur in all organisms of the group.
-     Calculate the difference of the two sets of majority EC numbers, leaving only the EC numbers which occur in Order majority, but not in Class majority.
-     Print the amount of EC numbers and their percentage of all EC numbers in Order, ie. how many of the EC numbers in Order do not exist in Class majority.


Result
------

::

    269/295 -> 91.2% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 100% majority for both
    249/295 -> 84.4% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 90% majority for both
    268/328 -> 81.7% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 80% majority for both
    272/342 -> 79.5% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 70% majority for both
    271/354 -> 76.6% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 60% majority for both
    269/380 -> 70.8% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 50% majority for both
    200/380 -> 52.6% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 40% majority for both
    146/389 -> 37.5% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 30% majority for both
     78/395 -> 19.7% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 20% majority for both
     36/425 ->  8.5% of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ 10% majority for both


Conclusion
----------
This group (Brachyspirales) is even more diverse than Enterobacteriales in :mod:`20`. It is, therefore, unlikely that a poor selection of example organisms caused the results in :mod:`20`.
It might still be the case that KEGG only contains poor representatives of the real group, which can never be ruled out.

Even though the threshold of majority decreases the same for group and super group, the size of the majority network is rather constant for Brachyspirales.
The size of the majority network for Spirochaetes, however, seems to increase quicker, yet linearly, with decreasing majority threshold. This indicates an even greater diversity in the super group Spirochaetes.
"""
from FEV_KEGG.Evolution.Taxonomy import NCBI
import FEV_KEGG.KEGG.Organism
from FEV_KEGG.Statistics.Percent import getPercentSentence


if __name__ == '__main__':
    FEV_KEGG.startProcessPool()
    taxonomy = NCBI.getTaxonomy()
    
    output = []
    
    #- Create a group of example organisms of Order.
    nodes = taxonomy.getOrganismNodesByPath('Bacteria/Spirochaetes/Brachyspirales')
    brachyspirales_organisms = FEV_KEGG.KEGG.Organism.Group( taxonomy.getOrganismAbbreviations( nodes ) )
    
    #- Create a super group of example organisms of Class, including the same organsims as the group of Order.
    nodes = taxonomy.getOrganismNodesByPath('Bacteria/Spirochaetes')
    spirochaetes_organisms = FEV_KEGG.KEGG.Organism.Group( taxonomy.getOrganismAbbreviations( nodes ) )
    
    #- REPEAT for varying majority-percentages:
    for majorityPercentage in range(100, 0, -10):
    
        #-     Calculate majority substance-ec graphs for both groups. Leaving only EC numbers which occur in all organisms of the group.
        brachyspirales_EC_graph = brachyspirales_organisms.majorityEcGraph(majorityPercentage)
        spirochaetes_EC_graph = spirochaetes_organisms.majorityEcGraph(majorityPercentage)
        
        #-     Calculate the difference of the two sets of majority EC numbers, leaving only the EC numbers which occur in Order majority, but not in Class majority.
        brachyspirales_EC_set = brachyspirales_EC_graph.getECs()
        spirochaetes_EC_set = spirochaetes_EC_graph.getECs()
        
        only_brachyspirales_EC_set = brachyspirales_EC_set.difference(spirochaetes_EC_set)
        
        #-     Print the amount of EC numbers and their percentage of all EC numbers in Order, ie. how many of the EC numbers in Order do not exist in Class majority.
        output.append( getPercentSentence(len(only_brachyspirales_EC_set), len(brachyspirales_EC_set)) + ' of EC numbers in Brachyspirales are new, compared to Spirochaetes majority @ ' + str(majorityPercentage) + '% majority for both' )
    
    
    for line in output:
        print(line)
    