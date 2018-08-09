"""
Question
--------
Taking an incomplete set of example organisms, forming a group (Enterobacteriales) and a super group (Gammaproteobacteria) of organisms, 
which functions occur in the consensus network of the group, but not in the consensus network of the super group?

Method
------
- Create a group of example organisms of Order Enterobacteriales.
- Create a group of example organisms of Class Gammaproteobacteria, including the same organsims as the group of Enterobacteriales.
- Calculate consensus substance-ec graphs for both groups. Leaving only EC numbers which occur in all organisms of the group.
- Calculate the difference of the two sets of consensus EC numbers, leaving only the EC numbers which occur in Enterobacteriales consensus, but not in Gammaproteobacteria consensus.
- Print these EC numbers and their percentage of all EC numbers in Enterobacteriales, ie. how many of the EC numbers in Enterobacteriales do not exist in Gammaproteobacteria consensus.

Result
------

::

    108 results
    1.1.1.100
    1.1.1.17
    1.1.1.193
    1.1.1.23
    1.1.1.267
    1.1.1.3
    1.1.1.343
    1.1.1.363
    1.1.1.44
    1.1.1.49
    1.1.1.85
    1.1.1.86
    1.17.1.8
    1.17.7.1
    1.17.7.3
    1.17.7.4
    1.2.1.38
    1.3.1.76
    1.5.1.20
    1.7.1.7
    1.8.1.2
    2.1.1.107
    2.1.1.14
    2.1.2.11
    2.1.2.3
    2.1.3.2
    2.1.3.3
    2.2.1.2
    2.2.1.6
    2.2.1.7
    2.3.1.1
    2.3.1.117
    2.3.1.157
    2.3.1.30
    2.3.1.41
    2.3.1.8
    2.3.3.13
    2.4.2.-
    2.4.2.1
    2.4.2.17
    2.4.2.18
    2.4.2.22
    2.4.2.8
    2.5.1.16
    2.5.1.19
    2.5.1.47
    2.5.1.61
    2.6.1.9
    2.7.1.-
    2.7.1.11
    2.7.1.148
    2.7.1.197
    2.7.1.199
    2.7.1.25
    2.7.1.26
    2.7.1.39
    2.7.1.40
    2.7.2.1
    2.7.2.4
    2.7.2.8
    2.7.4.16
    2.7.6.1
    2.7.7.18
    2.7.7.2
    2.7.7.23
    2.7.7.4
    2.7.7.60
    2.7.7.8
    2.7.8.7
    2.8.1.7
    3.1.1.31
    3.1.2.-
    3.1.3.15
    3.1.3.7
    3.2.2.9
    3.5.1.16
    3.5.1.18
    3.5.2.3
    3.5.4.10
    3.5.4.13
    3.5.4.19
    3.5.4.26
    3.6.1.23
    3.6.1.27
    3.6.1.31
    4.1.1.50
    4.1.3.-
    4.2.1.19
    4.2.1.33
    4.2.1.35
    4.2.1.9
    4.3.2.1
    4.3.3.7
    4.6.1.12
    4.99.1.4
    5.1.1.3
    5.1.1.7
    5.3.1.16
    5.3.1.28
    5.4.2.11
    5.4.2.7
    5.4.99.5
    6.1.1.17
    6.3.2.1
    6.3.2.3
    6.3.4.21
    6.3.4.5
    6.3.5.5
    108/190 -> 56.8% of EC numbers in Enterobacteriales are new, compared to Gammaproteobacteria consensus


Conclusion
----------
57% is a very high percentage. Maybe the example set of Gammaproteobacteria contained an organism with limited metabolic capabilities, eg. a parasite.
Next step is to limit 'consensus' to x percent of all organisms, ie. 1-x percent of organisms in a group may not contain an EC number, and it would still be added to the 'consensus', or rather majority.
"""
from FEV_KEGG.KEGG.File import cache
import FEV_KEGG.KEGG.Organism
from FEV_KEGG.Statistics.Percent import getPercentSentence


@cache(folder_path = 'experiments/19', file_name = 'enterobacteriales_SubstanceEcGraph')
def enterobacterialesEcGraph():
    #- Create a group of example organisms of Order Enterobacteriales.
    enterobacteriales_organisms_abbreviations = ['eco', 'ses', 'sfl', 'ent', 'esa', 'kpn', 'cko', 'ype', 'spe', 'buc']
    
    enterobacteriales_organisms = FEV_KEGG.KEGG.Organism.Group(organismAbbreviations = enterobacteriales_organisms_abbreviations)
    
    return enterobacteriales_organisms.consensusEcGraph()

@cache(folder_path = 'experiments/19', file_name = 'gammaproteobacteria_SubstanceEcGraph')
def gammaproteobacteriaEcGraph():
    #- Create a group of example organisms of Class Gammaproteobacteria, including the same organsims as the group of Enterobacteriales.
    enterobacteriales_organisms_abbreviations = ['eco', 'ses', 'sfl', 'ent', 'esa', 'kpn', 'cko', 'ype', 'spe', 'buc']
    gammaproteobacteria_organisms_abbreviations = ['hin', 'mht', 'xcc', 'vch', 'pae', 'acb', 'son', 'pha', 'amc', 'lpn', 'ftu', 'aha']
    gammaproteobacteria_organisms_abbreviations.extend(enterobacteriales_organisms_abbreviations) # extend with the sub-set, because they are also part of the set
    
    gammaproteobacteria_organisms = FEV_KEGG.KEGG.Organism.Group(organismAbbreviations = gammaproteobacteria_organisms_abbreviations)
    
    return gammaproteobacteria_organisms.consensusEcGraph()
    

if __name__ == '__main__':
    
    #- Calculate consensus substance-ec graphs for both groups. Leaving only EC numbers which occur in all organisms of the group.
    enterobacteriales_EC_graph = enterobacterialesEcGraph()
    gammaproteobacteria_EC_graph = gammaproteobacteriaEcGraph()
    
    #- Calculate the difference of the two sets of consensus EC numbers, leaving only the EC numbers which occur in Enterobacteriales consensus, but not in Gammaproteobacteria consensus.
    enterobacteriales_EC_set = enterobacteriales_EC_graph.getECs()
    gammaproteobacteria_EC_set = gammaproteobacteria_EC_graph.getECs()
    
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
    