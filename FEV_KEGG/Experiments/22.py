"""
Context
-------
:mod:`19` found many EC numbers new to an example group of Enterobacteriales vs. the super group of Gammaproteobacteria.
"108/190 -> 56.8% of EC numbers in Enterobacteriales are new, compared to Gammaproteobacteria consensus"

Question
--------
Could this be due to incomplete data in KEGG? Because substance-ec graphs are intersected to form the consensus.
Does the same result occur when intersecting the set of EC numbers itself?

Method
------
Similar to :mod:`19`, only intersect sets instead of networks.
- Create a group of example organisms of Order Enterobacteriales.
- Create a group of example organisms of Class Gammaproteobacteria, including the same organsims as the group of Enterobacteriales.
- Get sets of EC numbers from each graph.
- Calculate consensus set for both groups (Order and Class). Leaving only EC numbers which occur in all organisms of the group.
- Calculate the difference of the two sets of consensus EC numbers, leaving only the EC numbers which occur in Enterobacteriales consensus, but not in Gammaproteobacteria consensus.
- Print these EC numbers and their percentage of all EC numbers in Enterobacteriales, ie. how many of the EC numbers in Enterobacteriales do not exist in Gammaproteobacteria consensus.

Result
------

::

    107/190 -> 56.3% of EC numbers in Enterobacteriales are new, compared to Gammaproteobacteria consensus

Conclusion
----------
Only one EC number (0.5%) was lost due to the way the consensus is calculated. This kind of error in KEGG data should not be able to cause much difference in results.
"""
import FEV_KEGG.KEGG.Organism
from FEV_KEGG.Statistics.Percent import getPercentSentence


def enterobacterialesEcSets():
    #- Create a group of example organisms of Order Enterobacteriales.
    enterobacteriales_organisms_abbreviations = ['eco', 'ses', 'sfl', 'ent', 'esa', 'kpn', 'cko', 'ype', 'spe', 'buc']
    
    enterobacteriales_organisms = FEV_KEGG.KEGG.Organism.Group(organismAbbreviations = enterobacteriales_organisms_abbreviations)
    
    ecNumberSets = []
    
    for graph in enterobacteriales_organisms.ecGraphs().values():
        ecNumberSets.append(graph.getECs())
    
    return ecNumberSets


def gammaproteobacteriaEcSets():
    #- Create a group of example organisms of Class Gammaproteobacteria, including the same organsims as the group of Enterobacteriales.
    enterobacteriales_organisms_abbreviations = ['eco', 'ses', 'sfl', 'ent', 'esa', 'kpn', 'cko', 'ype', 'spe', 'buc']
    gammaproteobacteria_organisms_abbreviations = ['hin', 'mht', 'xcc', 'vch', 'pae', 'acb', 'son', 'pha', 'amc', 'lpn', 'ftu', 'aha']
    gammaproteobacteria_organisms_abbreviations.extend(enterobacteriales_organisms_abbreviations) # extend with the sub-set, because they are also part of the set
    
    gammaproteobacteria_organisms = FEV_KEGG.KEGG.Organism.Group(organismAbbreviations = gammaproteobacteria_organisms_abbreviations)
    
    ecNumberSets = []
    
    for graph in gammaproteobacteria_organisms.ecGraphs().values():
        ecNumberSets.append(graph.getECs())
    
    return ecNumberSets
    

if __name__ == '__main__':
    
    #- Get sets of EC numbers from each graph.
    enterobacteriales_EC_setList = enterobacterialesEcSets()
    gammaproteobacteria_EC_setList = gammaproteobacteriaEcSets()
    
    #- Calculate consensus set for both groups (Order and Class). Leaving only EC numbers which occur in all organisms of the group.
    enterobacteriales_EC_set = enterobacteriales_EC_setList.pop()
    for ecNumberSet in enterobacteriales_EC_setList:
        enterobacteriales_EC_set.intersection_update(ecNumberSet)
        
    gammaproteobacteria_EC_set = gammaproteobacteria_EC_setList.pop()
    for ecNumberSet in gammaproteobacteria_EC_setList:
        gammaproteobacteria_EC_set.intersection_update(ecNumberSet)
    
    #- Calculate the difference of the two sets of consensus EC numbers, leaving only the EC numbers which occur in Enterobacteriales consensus, but not in Gammaproteobacteria consensus.
    only_enterobacteriales_EC_set = enterobacteriales_EC_set.difference(gammaproteobacteria_EC_set)
    
    #- Print these EC numbers and their percentage of all EC numbers in Enterobacteriales, ie. how many of the EC numbers in Enterobacteriales do not exist in Gammaproteobacteria consensus.
    output = []
    
    for ec in only_enterobacteriales_EC_set:
        output.append(ec.__str__())
    
    output.sort()
    print(str(len(output)) + ' results')
    for line in output:
        print(line)
        
    print( getPercentSentence(len(only_enterobacteriales_EC_set), len(enterobacteriales_EC_set)) + ' of EC numbers in Enterobacteriales are new, compared to Gammaproteobacteria consensus' )