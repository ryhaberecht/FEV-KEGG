"""
Context
-------
The approach of building a consensus/majority graph of enzymes/EC numbers to find a core metabolism shared among several organisms has to be validated against previous research.
One such previous research deals with seven representative genomes of Thaumarchaeota, but manually curates all genes, including their associated EC numbers.
Kerou et al. (2016) list EC numbers in the additional file "Dataset_S02.xlsx", from which we extracted the ones on the sheets "Cell surface & glycosyl" and "Metabolism".
As with other validations before, we have to filter EC numbers which are outdated, or somehow not represented by KEGG's standard pathways. This is done here by restricting any EC number to the ones in NUKA, whis is done for both, their and our set of EC numbers.


Question
--------
Does the consensus/majority graph approach to core metabolism yield a similar set of EC numbers as the approach of Kerou et al. (2016)?


Method
------
- extract EC numbers from Kerou et al. (2016) by hand
- sanitise them by leaving only the ones found in NUKA
- also remove the ones with wildcards
- REPEAT with different groups
-     1. get group of organisms by clade 'Thaumarchaeota'
-     2. get only the seven organisms used by Kerou et al. (2016)
-     REPEAT for varying majority-percentages:
-         calculate EC numbers occuring in group's core metabolism
-         sanitise them by leaving only the ones found in NUKA
-         also remove the ones with wildcards
-         overlap their set with ours and print amount of EC numbers inside the intersection and falling off either side


Result
------

::

    Maj. %  others both   ours
    Thaumarchaeota:
    100%:   102     65     80
    90%:    73     94    139
    80%:    66    101    154
    70%:    65    102    156
    60%:    62    105    162
    50%:    60    107    167
    40%:    58    109    174
    30%:    53    114    188
    20%:    51    116    203
    10%:    46    121    240
    1%:    38    129    334
    
    Representative organisms:
    100%:    74     93    142
    90%:    74     93    142
    80%:    65    102    155
    70%:    65    102    161
    60%:    65    102    161
    50%:    61    106    174
    40%:    53    114    191
    30%:    53    114    191
    20%:    50    117    209
    10%:    47    120    245
    1%:    47    120    245


Conclusion
----------
When comparing the core metabolisms of all Thaumarchaeota known today, with only the ones from the seven representative organisms, there is not much difference.
This shows that these seven organisms are indeed very well chosen representatives.

Considering the amount of EC numbers falling off to either side:
The number of ECs only in our set is larger than the overlap, thus, we again see that core metabolisms created with our approach tend to be bigger than manually curated ones. The latter is most likely due to the fact that ECs which occur in all genomes do not necessarily have to be essential, while Kerou et al. aimed at only including essential ECs.
The number of ECs only in their set is also very high, accounting to roughly 65% of the overlap, or 40% of their set, and 20% of the overall set.
These ECs only in their set can not stem from ECs not in KEGG pathways at all, since we pre-filtered them using NUKA.
The most likely explanations seems to be that Kerou et al. were able to annotate many more EC numbers manually than KEGG's GENE database has stored to this date. This, again, would mean that KEGG's data is incomplete, which is strongly implied by the fact that even the collective graph (1% majority) does not contain 47 of the representative's EC numbers, which can only happen if these EC numbers are nowhere to be found in any of today's seven organisms in KEGG.

In conclusion of the effectiveness of our approach of building a core metabolism, we are left to say that completeness and quality of EC number annoations vary greatly, both within literature and KEGG.
Therefore, to achieve the most exact model of an organisms metabolism, one needs to apply further steps beyond our approach. Such steps may involve flux balance analysis with a manually curated list of 'essential' metabolites.
Still, however, when reducing the set of EC numbers to the ones known to standard KEGG pathways (using NUKA), core metabolisms created via our approach can be used to roughly compare the metabolic capabilities of closely, or even remotely related organisms, groups of organisms, and whole clades. 
"""
from FEV_KEGG.Evolution.Taxonomy import NCBI
import FEV_KEGG.KEGG.Organism as Organism
from FEV_KEGG.Graph.Elements import EcNumber
from FEV_KEGG.KEGG.NUKA import NUKA


if __name__ == '__main__':
    
    output = ['Maj. %\tothers\tboth\tours']
    
    #- extract EC numbers 
    theirECnumberStrings = ['1.1.1.-', '1.1.1.1', '1.1.1.100', '1.1.1.157', '1.1.1.22', '1.1.1.23', '1.1.1.25', '1.1.1.26', '1.1.1.261', '1.1.1.298', '1.1.1.3', '1.1.1.302', '1.1.1.37', '1.1.1.38', '1.1.1.40', '1.1.1.41', '1.1.1.60', '1.1.1.81', '1.1.1.85', '1.1.1.86', '1.1.1.88', '1.1.1.95', '1.1.5.2', '1.1.98.2', '1.11.1.15', '1.14.14.3', '1.14.99.39', '1.15.1.1', '1.16.1.1', '1.2.1.11', '1.2.1.16', '1.2.1.38', '1.2.1.59', '1.2.1.75', '1.2.1.76', '1.2.3.3', '1.2.7.-', '1.20.4.1', '1.3.1.-', '1.3.1.12', '1.3.1.84', '1.3.99.1', '1.4.1.1', '1.4.1.3', '1.5.1.12', '1.5.99.11', '1.5.99.8', '1.6.99.5', '1.7.2.1', '1.8.1.9', '1.8.7.1', '1.9.3.1', '2.1.1.13', '2.1.1.137', '2.1.2.1', '2.1.3.3', '2.2.1.1', '2.2.1.2', '2.2.1.6', '2.3.1.-', '2.3.1.1', '2.3.1.129', '2.3.1.157', '2.3.1.16', '2.3.1.182', '2.3.1.9', '2.3.3.1', '2.3.3.10', '2.3.3.13', '2.4.1.117', '2.4.1.212', '2.4.1.217', '2.4.1.227', '2.4.1.83', '2.4.2.-', '2.4.2.17', '2.4.2.18', '2.4.99.18', '2.5.1.1', '2.5.1.19', '2.5.1.41', '2.5.1.42', '2.5.1.47', '2.5.1.54', '2.5.1.6', '2.5.1.78', '2.5.1.9', '2.6.1.1', '2.6.1.11', '2.6.1.16', '2.6.1.42', '2.6.1.51', '2.6.1.9', '2.7.1.-', '2.7.1.161', '2.7.1.36', '2.7.1.39', '2.7.1.40', '2.7.1.71', '2.7.13.3', '2.7.2.3', '2.7.2.4', '2.7.2.8', '2.7.6.1', '2.7.7.12', '2.7.7.13', '2.7.7.2', '2.7.7.68', '2.7.8.13', '2.7.8.28', '2.7.8.5', '2.7.9.1', '2.7.9.2', '2.8.1.7', '3.1.2.6', '3.1.3.11', '3.1.3.15', '3.1.3.18', '3.1.3.25', '3.1.3.3', '3.4.11.1', '3.4.11.18', '3.4.11.2', '3.5.1.-', '3.5.1.102', '3.5.1.5', '3.5.2.10', '3.5.3.1', '3.5.3.11', '3.5.3.3', '3.5.4.-', '3.5.4.19', '3.5.4.29', '3.6.1.1', '3.6.1.27', '3.6.1.31', '3.6.1.7', '3.6.3.14', '4.1.1.19', '4.1.1.20', '4.1.1.48', '4.1.1.49', '4.1.2.13', '4.1.3.1', '4.1.3.27', '4.1.3.34', '4.1.3.6', '4.1.99.12', '4.2.1.1', '4.2.1.16', '4.2.1.20', '4.2.1.10', '4.2.1.11', '4.2.1.19', '4.2.1.2', '4.2.1.20', '4.2.1.3', '4.2.1.33', '4.2.1.51', '4.2.1.55', '4.2.1.9', '4.2.3.1', '4.2.3.4', '4.2.3.5', '4.3.1.19', '4.3.2.1', '4.4.1.16', '5.1.1.3', '5.1.3.1', '5.1.3.14', '5.1.3.2', '5.1.99.1', '5.3.1.1', '5.3.1.16', '5.3.1.24', '5.3.1.6', '5.3.1.8', '5.3.1.9', '5.3.3.2', '5.4.2.-', '5.4.2.10', '5.4.2.11', '5.4.2.12', '5.4.99.2', '5.4.99.5', '5.5.1.4', '6.2.1.1', '6.2.1.2', '6.2.1.3', '6.2.1.36', '6.2.1.5', '6.3.1.2', '6.3.2.-', '6.3.4.15', '6.3.4.5', '6.3.5.4', '6.3.5.5']
    theirECnumbers = set()
    for string in theirECnumberStrings:
        theirECnumbers.add( EcNumber(string) )
    
    #- sanitise them by leaving only the ones found in NUKA
    nuka = NUKA()
    nukaECs = nuka.substanceEcGraph.getECs()
    theirECnumbers.intersection_update( nukaECs )
    
    #- also remove the ones with wildcards
    theirECnumbers = EcNumber.removeWildcards(theirECnumbers)
    
    taxonomy = NCBI.getTaxonomy()
    #- REPEAT with different groups
    for i in range(1, 3):
        
        if i == 1:
            #-     1. get group of organisms by clade 'Thaumarchaeota'
            organisms = taxonomy.getOrganismAbbreviationsByPath('/Archaea/Thaumarchaeota', oneOrganismPerSpecies=False)
            output.append( 'Thaumarchaeota:' )
            
        elif i == 2:
            
            #-     2. get only the seven organisms used by Kerou et al. (2016)
            organisms = ['nga', 'nvn', 'nev', 'nmr', 'nkr', 'nbv', 'ndv']
            output.append( '\nRepresentative organisms:' )
        
        group = Organism.Group( organisms )
        
        #-     REPEAT for varying majority-percentages:
        for percentage in [100, 90, 80, 70, 60, 50, 40, 30, 20, 10 , 1]:
        
            #-         calculate EC numbers occuring in group's core metabolism
            ourECnumbers = group.majorityEcGraph(majorityPercentage = percentage, noMultifunctional = False).getECs()
            
            #-         sanitise them by leaving only the ones found in NUKA
            ourECnumbers.intersection_update( nukaECs )
            
            #-         also remove the ones with wildcards
            ourECnumbers = EcNumber.removeWildcards(ourECnumbers)
            
            #-         overlap their set with ours and print amount of EC numbers inside the intersection and falling off either side
            onlyInTheirs = theirECnumbers.difference( ourECnumbers )
            inBoth = theirECnumbers.intersection( ourECnumbers )
            onlyInOurs = ourECnumbers.difference( theirECnumbers )
            
            output.append(str(percentage) + '%:\t' + str(len(onlyInTheirs)) + '\t' + str(len(inBoth)) + '\t' + str(len(onlyInOurs)) )
        
    for line in output:
        print(line)
