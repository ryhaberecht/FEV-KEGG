"""
Context
-------
As concluded in experiment :mod:`30`, we first need to remove outdated EC numbers from Poot-Hernandez' set.
Then, we have to reduce both sets, theirs and ours, to the first three levels, and then compare them.
Now only the group of representative organisms and the Gammaproteobacteria group excluding 'unclassified' organisms is used.

Question
--------
Does the consensus/majority graph approach to core metabolism yield a similar set of EC numbers as the approach of Poot-Hernandez et al. (2015)?

Method
------
- extract EC numbers from Poot-Hernandez et al. (2015) by hand, any which are marked as blue (preserved)
- remove outdated EC numbers
- reduce set of EC numbers to first three levels
- REPEAT with different groups
-     1. get group of organisms deemed representative by Poot-Hernandez et al.
-     2. get group of organisms 'Gammaproteobacteria', excluding unclassified
-     REPEAT for varying majority-percentages:
-         calculate EC numbers occuring in group's core metabolism
-         reduce set of EC numbers to first three levels
-         overlap Poot-Hernandez' set with ours and print amount of EC numbers inside the intersection and falling off either side

Result
------

::

    Maj. %    others    both    ours
    Representative:
    100%:    43    7    1
     90%:    8    42    11
     80%:    1    49    21
     70%:    0    50    32
     60%:    0    50    39
     50%:    0    50    41
     40%:    0    50    49
     30%:    0    50    56
     20%:    0    50    67
     10%:    0    50    83
      1%:    0    50    108
    
    Gammaproteobacteria without unclassified:
    100%:    49    1    0
     90%:    2    48    18
     80%:    0    50    36
     70%:    0    50    40
     60%:    0    50    44
     50%:    0    50    51
     40%:    0    50    56
     30%:    0    50    66
     20%:    0    50    80
     10%:    0    50    88
      1%:    0    50    104

Conclusion
----------
Starting at 70% majority, there a no enzyme reactions in their core metabolism which do not also occur in ours.
While this method can not verify our approach, it at least rules out the most obvious path of falsification.

Regarding the full group of Gammaproteobacteria on our side, their set of EC numbers is fully covered by ours at a higher percentage, as to be expected.
Therefore, it seems to be a good idea to use the whole taxon, instead of only representative organisms, even though the taxon is more diverse than the chosen representatives, 
which shows in the higher count of EC numbers only in our set at any majority percentage.

On the other hand, our core metabolism is consistently larger, apart from the special case of a consensus (100%) core metabolism.
This could indicate that Poot-Hernandez et al. used a high percentage of occurence - between 100% and 90% - to define 'preserved'. Which percentage they used is, sadly, not documented.
But still, there is no point of clean overlap, showing that the two approaches yield fundamentally different results, with ours yielding a bigger core metabolism.
"""
from FEV_KEGG.Graph.Elements import EcNumber
from FEV_KEGG.Evolution.Taxonomy import NCBI
from FEV_KEGG.KEGG import Organism


if __name__ == '__main__':
    
    output = ['Maj. %\tothers\tboth\tours']
    
    #- extract EC numbers from Poot-Hernandez et al. (2015) by hand, any which are marked as blue (preserved)
    theirECnumberStrings = ['1.1.1.-', '1.1.1.100', '1.1.1.157', '1.1.1.158', '1.1.1.193', '1.1.1.205', '1.1.1.22', '1.1.1.267', '1.1.1.29', '1.1.1.30', '1.1.1.35', '1.1.1.36', '1.1.1.37', '1.1.1.41', '1.1.1.42', '1.1.1.43', '1.1.1.44', '1.1.1.49', '1.1.1.85', '1.1.1.86', '1.1.1.94', '1.17.1.2', '1.17.4.1', '1.17.4.2', '1.17.7.1', '1.2.1.-', '1.2.1.11', '1.2.1.12', '1.2.1.41', '1.2.1.70', '1.2.4.1', '1.2.4.2', '1.2.4.4', '1.3.1.12', '1.3.1.26', '1.3.3.1', '1.3.99.1', '1.4.3.16', '1.5.1.20', '1.5.1.3', '1.5.1.5', '1.8.1.4', '1.8.1.7', '1.8.1.9', '2.-.-.-', '2.1.1.45', '2.1.2.1', '2.1.2.2', '2.1.2.3', '2.1.2.9', '2.1.3.2', '2.2.1.6', '2.2.1.7', '2.3.1.-', '2.3.1.1', '2.3.1.109', '2.3.1.117', '2.3.1.12', '2.3.1.129', '2.3.1.15', '2.3.1.157', '2.3.1.16', '2.3.1.179', '2.3.1.180', '2.3.1.181', '2.3.1.39', '2.3.1.40', '2.3.1.41', '2.3.1.47', '2.3.1.51', '2.3.1.54', '2.3.1.61', '2.3.1.8', '2.3.1.89', '2.3.1.9', '2.3.3.1', '2.3.3.13', '2.4.1.182', '2.4.1.227', '2.4.2.1', '2.4.2.10', '2.4.2.11', '2.4.2.14', '2.4.2.17', '2.4.2.19', '2.4.2.2', '2.4.2.3', '2.4.2.4', '2.4.2.7', '2.5.1.-', '2.5.1.1', '2.5.1.10', '2.5.1.15', '2.5.1.3', '2.5.1.30', '2.5.1.31', '2.5.1.47', '2.5.1.48', '2.5.1.49', '2.5.1.55', '2.5.1.61', '2.5.1.7', '2.5.1.78', '2.5.1.9', '2.5.1.90', '2.6.1.-', '2.6.1.1', '2.6.1.11', '2.6.1.16', '2.6.1.17', '2.6.1.42', '2.6.1.62', '2.6.1.66', '2.6.1.85', '2.6.1.9', '2.7.1.130', '2.7.1.148', '2.7.1.23', '2.7.1.24', '2.7.1.26', '2.7.1.33', '2.7.1.40', '2.7.2.-', '2.7.2.11', '2.7.2.3', '2.7.2.4', '2.7.4.14', '2.7.4.16', '2.7.4.22', '2.7.4.3', '2.7.4.6', '2.7.4.7', '2.7.4.8', '2.7.4.9', '2.7.6.1', '2.7.6.3', '2.7.7.1', '2.7.7.18', '2.7.7.3', '2.7.7.38', '2.7.7.41', '2.7.7.6', '2.7.7.60', '2.7.7.7', '2.7.7.8', '2.7.8.-', '2.7.8.13', '2.7.8.24', '2.7.8.5', '2.7.8.7', '2.7.8.8', '2.8.1.6', '2.8.1.8', '3.1.3.27', '3.1.3.45', '3.1.3.5', '3.1.3.6', '3.5.1.-', '3.5.1.18', '3.5.1.47', '3.5.2.3', '3.5.4.19', '3.5.4.25', '3.6.1.-', '3.6.1.11', '3.6.1.13', '3.6.1.15', '3.6.1.17', '3.6.1.19', '3.6.1.27', '3.6.1.40', '3.6.1.41', '3.6.1.5', '4.1.1.11', '4.1.1.20', '4.1.1.23', '4.1.1.31', '4.1.1.36', '4.1.1.37', '4.1.1.49', '4.1.1.65', '4.1.3.38', '4.2.1.-', '4.2.1.11', '4.2.1.17', '4.2.1.18', '4.2.1.2', '4.2.1.20', '4.2.1.22', '4.2.1.24', '4.2.1.3', '4.2.1.35', '4.2.1.51', '4.2.1.52', '4.2.1.55', '4.2.1.60', '4.2.1.75', '4.2.1.9', '4.2.3.1', '4.3.1.17', '4.3.1.19', '4.3.2.1', '4.3.2.2', '4.6.1.12', '5.1.1.1', '5.1.1.3', '5.1.1.7', '5.3.1.1', '5.3.1.16', '5.3.1.9', '5.4.2.1', '5.4.2.10', '5.4.2.2', '5.4.2.7', '5.4.3.8', '5.4.99.18', '6.1.1.10', '6.1.1.11', '6.1.1.17', '6.1.1.4', '6.1.1.5', '6.1.1.9', '6.2.1.1', '6.2.1.5', '6.3.2.1', '6.3.2.10', '6.3.2.12', '6.3.2.13', '6.3.2.2', '6.3.2.3', '6.3.2.4', '6.3.2.6', '6.3.2.8', '6.3.2.9', '6.3.3.1', '6.3.3.3', '6.3.4.13', '6.3.4.14', '6.3.4.15', '6.3.4.18', '6.3.4.2', '6.3.4.4', '6.3.4.5', '6.3.5.2', '6.3.5.3', '6.3.5.4', '6.3.5.5', '6.4.1.2']
    theirECnumbers = set()
    
    for string in theirECnumberStrings:
        theirECnumbers.add( EcNumber(string) )
    
    #- remove outdated EC numbers
    outdatedEcNumberStrings = ['1.1.1.158', '1.17.1.2', '1.17.4.2', '1.3.1.26', '1.3.3.1', '1.3.99.1', '2.3.1.89', '2.4.2.11', '2.7.4.14', '3.5.1.47', '3.6.1.15', '3.6.1.19', '4.2.1.52', '4.2.1.60', '5.4.2.1']
    outdatedEcNumbers = set()
    for string in outdatedEcNumberStrings:
        outdatedEcNumbers.add( EcNumber(string) )
    theirECnumbers.difference_update( outdatedEcNumbers )
    
    #- reduce set of EC numbers to first three levels
    theirECnumbers = EcNumber.insertWildcards(theirECnumbers, keepLevels = 3, allowHigherWildcards = False)
    
    taxonomy = NCBI.getTaxonomy()
    #- REPEAT with different groups
    for i in range(1, 3):
        
        if i == 1:
            #-     1. get group of organisms deemed representative by Poot-Hernandez et al.
            representativeOrganisms = ['abo', 'aci', 'aeh', 'aha', 'bcc', 'bci', 'bfl', 'bpn', 'buc', 'cbu', 'cps', 'crp', 'csa', 'eco', 'ftu', 'hch', 'hdu', 'hha', 'ilo', 'lpn', 'maq', 'mca', 'msu', 'noc', 'pat', 'pcr', 'pfl', 'pha', 'pin', 'plu', 'ppr', 'rma', 'saz', 'sde', 'sdn', 'shm', 'tcx', 'vfi', 'vvu', 'xca']
            organisms = representativeOrganisms
            output.append( 'Representative:' )
            
        elif i == 2:
            
            #-     2. get group of organisms 'Gammaproteobacteria', excluding unclassified
            organisms = taxonomy.getOrganismAbbreviationsByPath('Gammaproteobacteria', exceptPaths='unclassified')
            output.append( '\nGammaproteobacteria without unclassified:' )
        
        group = Organism.Group( organisms )
        
        #-     REPEAT for varying majority-percentages:
        for percentage in [100, 90, 80, 70, 60, 50, 40, 30, 20, 10 , 1]:
        
            #-         calculate EC numbers occuring in group's core metabolism
            ourECnumbers = group.majorityEcGraph(majorityPercentage = percentage, noMultifunctional = False).getECs()
            
            #-         reduce set of EC numbers to first three levels
            ourECnumbers = EcNumber.insertWildcards(ourECnumbers, keepLevels = 3, allowHigherWildcards = False)
            
            #-         overlap Poot-Hernandez' set with ours and print amount of EC numbers inside the intersection and falling off either side
            onlyInTheirs = theirECnumbers.difference( ourECnumbers ) 
            inBoth = theirECnumbers.intersection( ourECnumbers )
            onlyInOurs = ourECnumbers.difference( theirECnumbers )
            
            output.append(str(percentage) + '%:\t' + str(len(onlyInTheirs)) + '\t' + str(len(inBoth)) + '\t' + str(len(onlyInOurs)) )
        
    for line in output:
        print(line)
