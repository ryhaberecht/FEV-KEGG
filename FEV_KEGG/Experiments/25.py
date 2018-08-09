"""
Context
-------
Similar to :mod:`21`, only with a different Class and Order, to see if the distribution of new enzymes by majority threshold is comparable between very different Classes.

Question
--------
Taking a complete set of organisms, forming a group (Desulfobacterales, Order) and a super group (Deltaproteobacteria, Class) of organisms, 
how many functions occur in the majority network of the group, but not in the majority network of the super group?
For varying majority-percentages.

Method
------
Same as :mod:`21`, only different Class and Order.
- Create a group of all organisms of Order to be found in KEGG.
- Create a super group of all organisms of Class to be found in KEGG, including the same organsims as the group of Order.
- REPEAT for varying majority-percentages:
-     Calculate majority substance-ec graphs for both groups. Leaving only EC numbers which occur in all organisms of the group.
-     Calculate the difference of the two sets of majority EC numbers, leaving only the EC numbers which occur in Order majority, but not in Class majority.
-     Print the amount of EC numbers and their percentage of all EC numbers in Order, ie. how many of the EC numbers in Order do not exist in Class majority.




Result
------

::

    224/238 -> 94.1% of EC numbers in Desulfobacterales are new, compared to Deltaproteobacteria majority @ 100% majority for both
     63/238 -> 26.5% of EC numbers in Desulfobacterales are new, compared to Deltaproteobacteria majority @ 90% majority for both
     65/299 -> 21.7% of EC numbers in Desulfobacterales are new, compared to Deltaproteobacteria majority @ 80% majority for both
     49/333 -> 14.7% of EC numbers in Desulfobacterales are new, compared to Deltaproteobacteria majority @ 70% majority for both
     37/373 ->  9.9% of EC numbers in Desulfobacterales are new, compared to Deltaproteobacteria majority @ 60% majority for both
     56/418 -> 13.4% of EC numbers in Desulfobacterales are new, compared to Deltaproteobacteria majority @ 50% majority for both
     45/470 ->  9.6% of EC numbers in Desulfobacterales are new, compared to Deltaproteobacteria majority @ 40% majority for both
     49/526 ->  9.3% of EC numbers in Desulfobacterales are new, compared to Deltaproteobacteria majority @ 30% majority for both
     59/602 ->  9.8% of EC numbers in Desulfobacterales are new, compared to Deltaproteobacteria majority @ 20% majority for both
     93/736 -> 12.6% of EC numbers in Desulfobacterales are new, compared to Deltaproteobacteria majority @ 10% majority for both


Conclusion
----------
In contrast to Spirochaetes/Brachyspirales from :mod:`21`, ther percentage of new EC numbers quickly drops.
Since both groups are calculated with the same decreasing majority threshold, the percentage of new EC numbers even oscillates.
Deltaproteobacteria/Desulfobacterales seem much less diverse than Spirochaetes/Brachyspirales.
"""
from FEV_KEGG.Evolution.Taxonomy import NCBI
import FEV_KEGG.KEGG.Organism
from FEV_KEGG.Statistics.Percent import getPercentSentence


def groupEcGraph(taxonomy):
    #- Create a group of example organisms of Order.
    nodes = taxonomy.getOrganismNodesByPath('Bacteria/Proteobacteria/Deltaproteobacteria/Desulfobacterales')
    brachyspirales_organisms = FEV_KEGG.KEGG.Organism.Group( taxonomy.getOrganismAbbreviations( nodes ) )
    
    return brachyspirales_organisms.majorityEcGraph(majorityPercentage)

def supergroupEcGraph(taxonomy):
    #- Create a super group of example organisms of Class, including the same organsims as the group of Order.
    nodes = taxonomy.getOrganismNodesByPath('Bacteria/Proteobacteria/Deltaproteobacteria')
    spirochaetes_organisms = FEV_KEGG.KEGG.Organism.Group( taxonomy.getOrganismAbbreviations( nodes ) )
    
    return spirochaetes_organisms.majorityEcGraph(majorityPercentage)
    

if __name__ == '__main__':
    
    taxonomy = NCBI.getTaxonomy()
    
    output = []
    
    #- REPEAT for varying majority-percentages:
    for majorityPercentage in range(100, 0, -10):
    
        #-     Calculate majority substance-ec graphs for both groups. Leaving only EC numbers which occur in all organisms of the group.
        brachyspirales_EC_graph = groupEcGraph(taxonomy)
        spirochaetes_EC_graph = supergroupEcGraph(taxonomy)
        
        #-     Calculate the difference of the two sets of majority EC numbers, leaving only the EC numbers which occur in Order majority, but not in Class majority.
        brachyspirales_EC_set = brachyspirales_EC_graph.getECs()
        spirochaetes_EC_set = spirochaetes_EC_graph.getECs()
        
        only_brachyspirales_EC_set = brachyspirales_EC_set.difference(spirochaetes_EC_set)
        
        #-     Print the amount of EC numbers and their percentage of all EC numbers in Order, ie. how many of the EC numbers in Order do not exist in Class majority.
        output.append( getPercentSentence(len(only_brachyspirales_EC_set), len(brachyspirales_EC_set)) + ' of EC numbers in Desulfobacterales are new, compared to Deltaproteobacteria majority @ ' + str(majorityPercentage) + '% majority for both' )
    
    
    for line in output:
        print(line)
    