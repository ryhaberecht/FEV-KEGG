"""
Question
--------
Same as :mod:`19`, only majority instead of consensus.
Taking an incomplete set of example organisms, forming a group (Enterobacteriales, Order) and a super group (Gammaproteobacteria, Class) of organisms, 
which functions occur in the majority network of the group, but not in the majority network of the super group?

Method
------
Same as :mod:`19`, only majority instead of consensus.
- Create a group of example organisms of Order Enterobacteriales.
- Create a group of example organisms of Class Gammaproteobacteria, including the same organsims as the group of Enterobacteriales.
- Calculate majority substance-ec graphs for both groups. Leaving only EC numbers which occur in all organisms of the group.
- Calculate the difference of the two sets of majority EC numbers, leaving only the EC numbers which occur in Enterobacteriales majority, but not in Gammaproteobacteria majority.
- Print these EC numbers and their percentage of all EC numbers in Enterobacteriales, ie. how many of the EC numbers in Enterobacteriales do not exist in Gammaproteobacteria majority.

Result
------

::

    222 results @ 80%
    1.1.-.-
    1.1.1.-
    1.1.1.103
    1.1.1.127
    1.1.1.169
    1.1.1.17
    1.1.1.215
    1.1.1.22
    1.1.1.262
    1.1.1.28
    1.1.1.284
    1.1.1.290
    1.1.1.336
    1.1.1.343
    1.1.1.381
    1.1.1.44
    1.1.1.57
    1.1.1.79
    1.1.1.81
    1.1.2.3
    1.1.5.12
    1.1.98.6
    1.11.1.21
    1.11.1.9
    1.14.11.17
    1.14.14.5
    1.16.1.3
    1.2.1.-
    1.2.1.10
    1.2.1.19
    1.2.1.71
    1.2.1.72
    1.2.5.1
    1.2.7.-
    1.2.7.1
    1.3.1.-
    1.3.1.10
    1.3.1.28
    1.3.5.3
    1.3.98.3
    1.4.1.13
    1.4.1.14
    1.4.1.4
    1.4.3.16
    1.4.4.2
    1.4.5.1
    1.5.1.38
    1.5.1.41
    1.6.1.1
    1.6.5.2
    1.7.1.15
    1.7.1.7
    1.7.5.1
    1.7.99.-
    1.7.99.1
    1.8.4.10
    1.8.4.14
    1.8.4.8
    2.1.1.37
    2.1.2.10
    2.2.1.9
    2.3.1.109
    2.3.1.183
    2.3.1.243
    2.3.1.29
    2.3.1.40
    2.3.1.46
    2.3.1.54
    2.3.1.57
    2.3.3.9
    2.4.1.1
    2.4.1.18
    2.4.2.1
    2.4.2.22
    2.4.2.3
    2.4.2.4
    2.5.1.-
    2.5.1.129
    2.5.1.16
    2.5.1.17
    2.5.1.48
    2.5.1.74
    2.6.1.19
    2.6.1.57
    2.6.1.81
    2.7.1.-
    2.7.1.11
    2.7.1.12
    2.7.1.15
    2.7.1.16
    2.7.1.167
    2.7.1.17
    2.7.1.191
    2.7.1.193
    2.7.1.194
    2.7.1.196
    2.7.1.197
    2.7.1.199
    2.7.1.2
    2.7.1.201
    2.7.1.202
    2.7.1.205
    2.7.1.208
    2.7.1.22
    2.7.1.25
    2.7.1.35
    2.7.1.45
    2.7.1.5
    2.7.1.50
    2.7.1.56
    2.7.1.59
    2.7.1.6
    2.7.1.73
    2.7.1.89
    2.7.7.1
    2.7.7.12
    2.7.7.13
    2.7.7.24
    2.7.7.27
    2.7.7.58
    2.7.7.70
    2.7.7.73
    2.7.7.75
    2.7.8.-
    2.7.8.37
    2.7.8.7
    2.7.9.3
    2.8.1.2
    2.8.1.4
    2.9.1.1
    3.1.1.11
    3.1.1.32
    3.1.1.4
    3.1.1.45
    3.1.1.5
    3.1.2.12
    3.1.2.28
    3.1.3.10
    3.1.3.102
    3.1.3.104
    3.1.3.2
    3.1.3.4
    3.1.3.6
    3.1.3.81
    3.1.3.82
    3.1.3.83
    3.1.3.89
    3.1.4.14
    3.1.4.16
    3.2.1.1
    3.2.1.196
    3.2.1.20
    3.2.1.21
    3.2.1.22
    3.2.1.23
    3.2.1.86
    3.2.1.93
    3.2.2.4
    3.3.2.1
    3.4.11.23
    3.4.13.-
    3.5.1.-
    3.5.1.16
    3.5.1.19
    3.5.1.96
    3.5.3.11
    3.5.3.23
    3.5.4.4
    3.5.99.6
    3.6.1.26
    3.6.1.45
    3.6.1.67
    3.6.1.7
    3.6.1.9
    4.1.1.11
    4.1.1.17
    4.1.1.18
    4.1.1.50
    4.1.2.4
    4.1.2.40
    4.1.2.48
    4.1.2.50
    4.1.3.1
    4.1.3.36
    4.1.99.17
    4.1.99.19
    4.2.1.113
    4.2.1.126
    4.2.1.7
    4.2.1.8
    4.2.3.12
    4.2.3.3
    4.2.99.20
    4.3.1.1
    4.4.1.15
    4.4.1.16
    4.4.1.21
    4.6.1.17
    5.1.3.14
    5.1.3.20
    5.1.3.3
    5.1.3.4
    5.1.3.9
    5.3.1.12
    5.3.1.14
    5.3.1.4
    5.3.1.5
    5.3.1.8
    5.3.3.2
    5.3.3.8
    5.4.2.11
    5.4.2.12
    5.4.2.7
    5.4.4.2
    6.2.1.1
    6.2.1.20
    6.2.1.26
    6.3.1.1
    6.3.1.20
    6.3.1.5
    6.3.2.14
    6.3.4.21
    222/585 -> 37.9% of EC numbers in Enterobacteriales are new, compared to Gammaproteobacteria majority


Conclusion
----------
Old results, with 100% majority, from :mod:`19`: 108/190 -> 56.8%
At 80% majority, many more EC numbers exist within both, group and super group. Now, 585 instead of 190 in Enterobacteriales. From this bigger set, however, only 38% are now considered 'new', in comparison to Gammaproteobacteria.

Both sets seem to be very diverse, but this could be a misinterpretation due to,
a) poor choice of example organisms, however, this problem should have already been mitigated by introducing a majority graph.
b) not enough organisms per group to really form a characteristic 'core metabolism'.

On the other hand, it might be that with growing set sizes, there will be no difference left, meaning there is no characteristic difference in 'core metabolisms' between Order and Class of bacteria. This might be because,
a) differences are only visible between groups farther apart in taxonomy, eg. Order and Phylum, or even Domain. At least LUCA should definitely show differences.
b) there are, in general, no significant differences in 'core metabolisms' in prokaryotes, maybe due to horizontal gene transfer.

But the general expectation is to find a common 'core metabolism', which might be rather small, and should usually contain, for example, amino acid synthesis. Expected exceptions are symbionts.
"""
from FEV_KEGG.KEGG.File import cache
import FEV_KEGG.KEGG.Organism
from FEV_KEGG.Statistics.Percent import getPercentSentence


majorityPercentage = 80

enterobacteriales_organisms_abbreviations = ['eco', 'ses', 'sfl', 'ent', 'esa', 'kpn', 'cko', 'ype', 'spe', 'buc']
gammaproteobacteria_organisms_abbreviations = ['hin', 'mht', 'xcc', 'vch', 'pae', 'acb', 'son', 'pha', 'amc', 'lpn', 'ftu', 'aha']
gammaproteobacteria_organisms_abbreviations.extend(enterobacteriales_organisms_abbreviations) # extend with the sub-set, because they are also part of the set

@cache(folder_path = 'experiments/20', file_name = 'enterobacteriales_SubstanceEcGraph')
def enterobacterialesEcGraph():
    #- Create a group of example organisms of Order Enterobacteriales.
    enterobacteriales_organisms = FEV_KEGG.KEGG.Organism.Group(organismAbbreviations = enterobacteriales_organisms_abbreviations)
    
    return enterobacteriales_organisms.majorityEcGraph(majorityPercentage)


@cache(folder_path = 'experiments/20', file_name = 'gammaproteobacteria_SubstanceEcGraph')
def gammaproteobacteriaEcGraph():
    #- Create a group of example organisms of Class Gammaproteobacteria, including the same organsims as the group of Enterobacteriales.
    gammaproteobacteria_organisms = FEV_KEGG.KEGG.Organism.Group(organismAbbreviations = gammaproteobacteria_organisms_abbreviations)
    
    return gammaproteobacteria_organisms.majorityEcGraph(majorityPercentage)
    

if __name__ == '__main__':
    
    #- Calculate majority substance-ec graphs for both groups. Leaving only EC numbers which occur in all organisms of the group.
    enterobacteriales_EC_graph = enterobacterialesEcGraph()
    gammaproteobacteria_EC_graph = gammaproteobacteriaEcGraph()
    
    #- Calculate the difference of the two sets of majority EC numbers, leaving only the EC numbers which occur in Enterobacteriales majority, but not in Gammaproteobacteria majority.
    enterobacteriales_EC_set = enterobacteriales_EC_graph.getECs()
    gammaproteobacteria_EC_set = gammaproteobacteria_EC_graph.getECs()
    
    only_enterobacteriales_EC_set = enterobacteriales_EC_set.difference(gammaproteobacteria_EC_set)
    
    #- Print these EC numbers and their percentage of all EC numbers in Enterobacteriales, ie. how many of the EC numbers in Enterobacteriales do not exist in Gammaproteobacteria majority.
    output = []
    
    for ec in only_enterobacteriales_EC_set:
        output.append(ec.__str__())
    
    output.sort()
    print(str(len(output)) + ' results @ ' + str(majorityPercentage) + '%')
    for line in output:
        print(line)
        
    print( getPercentSentence(len(only_enterobacteriales_EC_set), len(enterobacteriales_EC_set)) + ' of EC numbers in Enterobacteriales are new, compared to Gammaproteobacteria majority' )
    