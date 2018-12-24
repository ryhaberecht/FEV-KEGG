"""
Context
-------
The approach of building a consensus/majority graph of enzymes/EC numbers to find a core metabolism shared among several organisms has to be validated against previous research.
One such previous research deals with Thermus thermophilus HB27, but not with its complete metabolism, but only the essential "core" parts. They enhance the information provided by the KEGG database with information from MetaCyc and, most importantly, from manual curation with organism specific knowledge.
Lee et al. (2014) list EC numbers in additional file 1 "12934_2014_968_MOESM1_ESM.xls".
As with other validations before, we have to filter EC numbers which are outdated, or somehow not represented by KEGG's standard pathways. This is done here by restricting any EC number to the ones in NUKA, whis is done for both, their and our set of EC numbers.



Question
--------
Does the consensus/majority graph approach to core metabolism yield a similar set of EC numbers as the approach of Lee et al. (2014)?


Method
------
- extract EC numbers from Lee et al. (2014) by hand
- sanitise them by leaving only the ones found in NUKA
- also remove the ones with wildcards
- REPEAT with different groups
-     1. get group of organisms 'Thermus thermophilus'
-     2. get only the organism 'Thermus thermophilus HB27'
-     REPEAT for varying majority-percentages:
-         calculate EC numbers occuring in group's core metabolism
-         sanitise them by leaving only the ones found in NUKA
-         also remove the ones with wildcards
-         overlap their set with ours and print amount of EC numbers inside the intersection and falling off either side


Result
------

::

    Maj. %  others  both   ours
    All Thermus thermophilus:
    100%:    118    306    106
    90%:    118    306    106
    80%:    118    306    106
    70%:    111    313    119
    60%:    111    313    119
    50%:    102    322    138
    40%:    102    322    138
    30%:    102    322    138
    20%:     95    329    160
    10%:     95    329    160
    1%:     95    329    160
    
    Thermus thermophilus HB27:
    100%:    98    326    119
    90%:    98    326    119
    80%:    98    326    119
    70%:    98    326    119
    60%:    98    326    119
    50%:    98    326    119
    40%:    98    326    119
    30%:    98    326    119
    20%:    98    326    119
    10%:    98    326    119
    1%:    98    326    119



Conclusion
----------
Comparing the core metabolism of all Thermus thermophilus, we see less overlap between both sets of EC numbers, which was to be expected.
Still, this shows that the variance between Thermus thermophilus subspecies is not huge, but certainly visible.
Three EC numbers only overlap when using the consensus metabolism of the whole group (329 vs. 326), indicating that Lee et al. manually added these EC numbers to their HB27 model, while they also exist in other subspecies known to KEGG.

For further analysis, we only regard the metabolism of the HB27 subspecies.

About 30% of overlapping EC numbers fall off either side. This shows a significant disrepancy between today's data and/or the data Lee et al. added manually.
If there were mainly EC numbers occuring only in their set, they would clearly stem from the manual addition.
However, even more EC numbers occur only in our set, which raises the question if that many EC numbers could have been added to these organisms in KEGG since 2014.
While this is possible, we sadly have no way to verify this.
All in all, the overlap is about 60% of the total sum of all EC numbers occuring in either set. This at least shows a fundamental consensus between both approaches and data sets.
"""
from FEV_KEGG.Evolution.Taxonomy import NCBI
import FEV_KEGG.KEGG.Organism as Organism
from FEV_KEGG.Graph.Elements import EcNumber
from FEV_KEGG.KEGG.NUKA import NUKA


if __name__ == '__main__':
    
    output = ['Maj. %\tothers\tboth\tours']
    
    #- extract EC numbers from Lee et al. (2014) by hand
    theirECnumberStrings = ['1.1.1.1', '1.1.1.100', '1.1.1.103', '1.1.1.125', '1.1.1.157', '1.1.1.158', '1.1.1.169', '1.1.1.193', '1.1.1.205', '1.1.1.23', '1.1.1.26', '1.1.1.267', '1.1.1.27', '1.1.1.282', '1.1.1.3', '1.1.1.31', '1.1.1.35', '1.1.1.37', '1.1.1.38', '1.1.1.40', '1.1.1.42', '1.1.1.69', '1.1.1.81', '1.1.1.85', '1.1.1.86', '1.1.1.87', '1.1.1.94', '1.1.1.95', '1.1.19.14', '1.1.3.15', '1.1.5.2', '1.10.2.2', '1.13.11.15', '1.13.12.16', '1.13.12.4', '1.14.11.-', '1.14.13.-', '1.14.13.149', '1.14.14.9', '1.16.8.1', '1.17.1.17', '1.17.1.2', '1.17.4.1', '1.17.7.1', '1.2.1.-', '1.2.1.10', '1.2.1.11', '1.2.1.12', '1.2.1.16', '1.2.1.2', '1.2.1.22', '1.2.1.3', '1.2.1.38', '1.2.1.39', '1.2.1.41', '1.2.1.43', '1.2.1.60', '1.2.3.13', '1.2.4.2', '1.2.4.4', '1.2.7.2', '1.2.7.3', '1.2.7.5', '1.3.1.13', '1.3.1.2', '1.3.1.76', '1.3.3.3', '1.3.3.4', '1.3.5.1', '1.3.5.2', '1.3.5.3', '1.3.5.5', '1.3.8.-', '1.3.8.1', '1.3.8.6', '1.3.8.7', '1.3.99.1', '1.3.99.22', '1.3.99.25', '1.3.99.28', '1.4.1.1', '1.4.1.13', '1.4.1.14', '1.4.1.3', '1.4.3.16', '1.4.4.2', '1.5.1.12', '1.5.1.2', '1.5.1.20', '1.5.1.3', '1.5.1.5', '1.5.7.1', '1.5.99.8', '1.6.1.1', '1.6.1.2', '1.6.5.-', '1.6.5.10', '1.6.5.3', '1.7.7.1', '1.8.1.4', '1.8.1.9', '1.8.2.1', '1.8.2.2', '1.8.4.8', '1.9.3.1', '2.1.1.107', '2.1.1.13', '2.1.1.130', '2.1.1.131', '2.1.1.132', '2.1.1.133', '2.1.1.148', '2.1.1.151', '2.1.1.163', '2.1.1.17', '2.1.1.201', '2.1.1.222', '2.1.1.45', '2.1.1.64', '2.1.1.71', '2.1.2.1', '2.1.2.10', '2.1.2.11', '2.1.2.2', '2.1.2.3', '2.1.3.2', '2.1.3.3', '2.2.1.1', '2.2.1.2', '2.2.1.6', '2.2.1.7', '2.3.1.-', '2.3.1.1', '2.3.1.15', '2.3.1.157', '2.3.1.16', '2.3.1.168', '2.3.1.174', '2.3.1.179', '2.3.1.180', '2.3.1.181', '2.3.1.29', '2.3.1.30', '2.3.1.31', '2.3.1.35', '2.3.1.51', '2.3.1.54', '2.3.1.61', '2.3.1.8', '2.3.1.9', '2.3.3.1', '2.3.3.13', '2.3.3.14', '2.3.3.9', '2.4.1.1', '2.4.1.21', '2.4.1.217', '2.4.1.227', '2.4.2.-', '2.4.2.1', '2.4.2.10', '2.4.2.11', '2.4.2.14', '2.4.2.17', '2.4.2.18', '2.4.2.19', '2.4.2.21', '2.4.2.22', '2.4.2.3', '2.4.2.4', '2.4.2.7', '2.4.2.8', '2.4.2.9', '2.5.1.-', '2.5.1.1', '2.5.1.10', '2.5.1.15', '2.5.1.16', '2.5.1.17', '2.5.1.19', '2.5.1.29', '2.5.1.3', '2.5.1.32', '2.5.1.39', '2.5.1.47', '2.5.1.49', '2.5.1.54', '2.5.1.6', '2.5.1.61', '2.5.1.68', '2.5.1.7', '2.5.1.72', '2.5.1.78', '2.5.1.9', '2.5.1.90', '2.6.1.-', '2.6.1.1', '2.6.1.11', '2.6.1.13', '2.6.1.16', '2.6.1.19', '2.6.1.39', '2.6.1.42', '2.6.1.44', '2.6.1.45', '2.6.1.51', '2.6.1.52', '2.6.1.66', '2.6.1.82', '2.6.1.85', '2.6.1.9', '2.7.1.107', '2.7.1.11', '2.7.1.12', '2.7.1.148', '2.7.1.15', '2.7.1.165', '2.7.1.2', '2.7.1.21', '2.7.1.23', '2.7.1.24', '2.7.1.25', '2.7.1.26', '2.7.1.33', '2.7.1.39', '2.7.1.4', '2.7.1.40', '2.7.1.45', '2.7.1.48', '2.7.1.49', '2.7.1.6', '2.7.1.71', '2.7.2.-', '2.7.2.1', '2.7.2.11', '2.7.2.2', '2.7.2.3', '2.7.2.4', '2.7.2.8', '2.7.4.14', '2.7.4.15', '2.7.4.16', '2.7.4.22', '2.7.4.3', '2.7.4.6', '2.7.4.7', '2.7.4.8', '2.7.4.9', '2.7.6.1', '2.7.6.3', '2.7.7.12', '2.7.7.13', '2.7.7.18', '2.7.7.2', '2.7.7.22', '2.7.7.23', '2.7.7.27', '2.7.7.3', '2.7.7.4', '2.7.7.41', '2.7.7.60', '2.7.7.62', '2.7.8.13', '2.7.8.26', '2.7.8.8', '2.7.9.1', '2.7.9.2', '2.8.1.-', '2.8.1.10', '2.8.1.6', '2.8.1.7', '2.8.1.8', '2.8.3.5', '3.1.1.17', '3.1.2.20', '3.1.2.4', '3.1.2.6', '3.1.3.1', '3.1.3.10', '3.1.3.11', '3.1.3.12', '3.1.3.15', '3.1.3.18', '3.1.3.2', '3.1.3.25', '3.1.3.27', '3.1.3.3', '3.1.3.5', '3.1.3.64', '3.1.3.70', '3.1.3.71', '3.1.5.1', '3.1.7.2', '3.2.1.20', '3.2.1.21', '3.2.1.22', '3.2.1.23', '3.2.1.28', '3.2.1.93', '3.2.2.14', '3.2.2.3', '3.2.2.5', '3.2.2.9', '3.5.1.-', '3.5.1.1', '3.5.1.10', '3.5.1.14', '3.5.1.19', '3.5.1.20', '3.5.1.25', '3.5.1.42', '3.5.1.54', '3.5.1.6', '3.5.2.2', '3.5.2.3', '3.5.3.1', '3.5.3.11', '3.5.4.10', '3.5.4.13', '3.5.4.16', '3.5.4.19', '3.5.4.25', '3.5.4.26', '3.5.4.5', '3.5.4.9', '3.5.99.6', '3.6.1.-', '3.6.1.1', '3.6.1.13', '3.6.1.15', '3.6.1.17', '3.6.1.19', '3.6.1.22', '3.6.1.27', '3.6.1.31', '3.6.1.40', '3.6.1.5', '3.6.1.6', '3.6.1.7', '3.6.3.14', '3.6.5.-', '3.7.1.16', '4.1.1.-', '4.1.1.11', '4.1.1.19', '4.1.1.21', '4.1.1.23', '4.1.1.31', '4.1.1.36', '4.1.1.37', '4.1.1.43', '4.1.1.48', '4.1.1.49', '4.1.1.50', '4.1.1.65', '4.1.1.68', '4.1.2.-', '4.1.2.13', '4.1.2.14', '4.1.2.17', '4.1.2.25', '4.1.2.4', '4.1.2.5', '4.1.3.-', '4.1.3.1', '4.1.3.16', '4.1.3.24', '4.1.3.27', '4.1.3.38', '4.1.3.4', '4.1.3.40', '4.1.99.12', '4.2.1.10', '4.2.1.11', '4.2.1.113', '4.2.1.114', '4.2.1.12', '4.2.1.17', '4.2.1.18', '4.2.1.19', '4.2.1.2', '4.2.1.20', '4.2.1.24', '4.2.1.3', '4.2.1.33', '4.2.1.51', '4.2.1.52', '4.2.1.59', '4.2.1.75', '4.2.1.9', '4.2.3.1', '4.2.3.3', '4.2.3.4', '4.2.3.5', '4.3.1.17', '4.3.1.19', '4.3.2.1', '4.3.2.2', '4.4.1.21', '4.4.1.5', '4.4.1.8', '4.6.1.12', '4.99.1.1', '4.99.1.4', '5.1.1.1', '5.1.1.3', '5.1.3.1', '5.1.3.14', '5.1.3.2', '5.1.99.1', '5.3.1.1', '5.3.1.16', '5.3.1.23', '5.3.1.24', '5.3.1.6', '5.3.1.8', '5.3.1.9', '5.3.3.-', '5.3.3.10', '5.3.3.18', '5.3.3.2', '5.3.99.10', '5.4.1.2', '5.4.2.1', '5.4.2.10', '5.4.2.2', '5.4.2.7', '5.4.2.8', '5.4.3.8', '5.4.99.18', '5.4.99.2', '5.4.99.5', '5.5.1.18', '5.5.1.4', '6.2.1.1', '6.2.1.14', '6.2.1.17', '6.2.1.3', '6.2.1.30', '6.2.1.5', '6.2.1.9', '6.3.1.10', '6.3.1.2', '6.3.1.5', '6.3.2.-', '6.3.2.1', '6.3.2.10', '6.3.2.13', '6.3.2.17', '6.3.2.4', '6.3.2.5', '6.3.2.6', '6.3.2.8', '6.3.2.9', '6.3.3.1', '6.3.3.2', '6.3.4.13', '6.3.4.15', '6.3.4.18', '6.3.4.2', '6.3.4.3', '6.3.4.4', '6.3.4.5', '6.3.4.6', '6.3.5.1', '6.3.5.10', '6.3.5.11', '6.3.5.2', '6.3.5.3', '6.3.5.4', '6.3.5.5', '6.3.5.9', '6.4.1.2', '6.4.1.3', '6.4.1.4', '6.6.1.1', '6.6.1.2']
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
            #-     1. get group of organisms 'Thermus thermophilus'
            organisms = taxonomy.getOrganismAbbreviationsByName('Thermus thermophilus')
            output.append( 'All Thermus thermophilus:' )
            
        elif i == 2:
            
            #-     2. get only the organism 'Thermus thermophilus HB27'
            organisms = ['tth']
            output.append( '\nThermus thermophilus HB27:' )
        
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
