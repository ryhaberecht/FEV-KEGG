"""
Context
-------
The approach of building a consensus/majority graph of enzymes/EC numbers to find a core metabolism shared among several organisms has to be validated against previous research.
One such previous research deals with Gammaproteobacteria, but uses a slightly different approach, first calculating Enzymatic Step Sequences (ESS) and comparing these among species.
Poot-Hernandez et al. (2015) list EC numbers associated with enzymes from ESS shared among most of the 40 representative species of Gammaproteobacteria in File S2, marked with the hexadecimal code for blue.
It does not matter in which category the EC numbers occur, since "highly preserved" etc. only represents the percentage of EC numbers in the pathway, which are actually common in Gammaproteobacteria, and are thus marked as blue. 
These EC numbers are used to validate the approach of this library. EC numbers containing wildcards (e.g. 1.2.-.-) are excluded on both sides, to minimise statistical skew.
Multifunctinal enzymes are, however, included, because they were not excluded by Poot-Hernandez et al.
Our group consists of:
1) The exact same organisms deemed representative by Poot-Hernandez et al. (Table 2 lists them)
2) All organisms of Gammaproteobacteria
3) All organisms of Gammaproteobacteria, excluding 'unclassified' organisms


Question
--------
Does the consensus/majority graph approach to core metabolism yield a similar set of EC numbers as the approach of Poot-Hernandez et al. (2015)?


Method
------
- extract EC numbers from Poot-Hernandez et al. (2015) by hand, any which are marked as blue (preserved)
- REPEAT with different groups
-     1. get group of organisms deemed representative by Poot-Hernandez et al.
-     2. get group of organisms 'Gammaproteobacteria'
-     3. get group of organisms 'Gammaproteobacteria', excluding unclassified
-     REPEAT for varying majority-percentages:
-         calculate EC numbers occuring in group's core metabolism
-         overlap Poot-Hernandez' set with ours and print amount of EC numbers inside the intersection and falling off either side


Result
------

::

    Maj. %    others    both    ours
    Representative:
    100%:    228    9    14
    90%:    128    109    65
    80%:    71    166    113
    70%:    54    183    163
    60%:    44    193    209
    50%:    39    198    259
    40%:    33    204    315
    30%:    26    211    383
    20%:    24    213    466
    10%:    22    215    659
    1%:    15    222    1059
    
    Gammaproteobacteria:
    100%:    235    2    0
    90%:    85    152    98
    80%:    51    186    174
    70%:    43    194    219
    60%:    39    198    263
    50%:    32    205    336
    40%:    30    207    402
    30%:    25    212    497
    20%:    23    214    620
    10%:    21    216    760
    1%:    18    219    1070
    
    Gammaproteobacteria without unclassified:
    100%:    235    2    0
    90%:    85    152    98
    80%:    51    186    174
    70%:    43    194    219
    60%:    39    198    264
    50%:    32    205    335
    40%:    30    207    402
    30%:    25    212    500
    20%:    23    214    620
    10%:    21    216    760
    1%:    18    219    1070



Conclusion
----------
Regarding the result at 1% of the representative group:
There are 15 EC numbers which occur in Poot-Hernandez' but nowhere in our organisms.
Looking at the EC numbers in detail, this is easily explained by outdated used by Poot-Hernandez:
1.1.1.158   deleted 2013
1.17.1.2    deleted 2016
1.17.4.2    in none of the 40 organisms of today
1.3.1.26    deleted 2013
1.3.3.1     deleted 2011
1.3.99.1    deleted 2014
2.3.1.89    in none of the 40 organisms of today
2.4.2.11    deleted 2013
2.7.4.14    in none of the 40 organisms of today
3.5.1.47    in none of the 40 organisms of today
3.6.1.15    in none of the 40 organisms of today
3.6.1.19    deleted 2016
4.2.1.52    deleted 2012
4.2.1.60    deleted 2012
5.4.2.1     deleted 2013

The experiment should be re-run without the now obsolete EC numbers above, to avoid skewed results.


Regarding the difference between Gammaproteobacteria with and without 'unclassified' organisms:
In the upper (100/90%) and lower (20/10/1%) regions of majority, there is no difference. Only in the middle regions some results differ by one to three counts.
This comes as no surprise, because the group of Gammaproteobacteria consists of 1074 organisms, while only one of them can be excepted as 'unclassified'.
This does not leave much room for difference. Still, it might be best practice to always exclude 'unclassified' organisms, if only for their unknown position in taxonomy.


As to be expected, using all Gammaproteobacteria organisms, e.g. at 70% majority, yields a higher overlap between our core metabolism and theirs.
But the difference is small, implying that the 40 organisms were well-selected representatives. 
Still, however, there is a high number of EC numbers only found in their core metabolism. This might result from their methodology. 
After creating the Enzymatic Step Sequences (ESS), they reduced the EC numbers in the ESS to contain only three levels, i.e. abstracting from substrate specificity keeping only reaction types. 
Then, at some undocumented point, they translated the set of three-level EC numbers back to the original set of EC numbers. But this results in a list of EC numbers which reactions are preserved, 
not which substrate specificities are preserved! In order to accurately follow their definition of preservation, we have to reduce both sets of EC numbers, ours and theirs, to their first three levels, 
and then re-run the experiment.
"""
from FEV_KEGG.Evolution.Taxonomy import NCBI
import FEV_KEGG.KEGG.Organism as Organism
from FEV_KEGG.Graph.Elements import EcNumber


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
    
    theirECnumbers = EcNumber.removeWildcards(theirECnumbers)
    
    taxonomy = NCBI.getTaxonomy()
    #- REPEAT with different groups
    for i in range(1, 4):
        
        if i == 1:
            #-     1. get group of organisms deemed representative by Poot-Hernandez et al.
            representativeOrganisms = ['abo', 'aci', 'aeh', 'aha', 'bcc', 'bci', 'bfl', 'bpn', 'buc', 'cbu', 'cps', 'crp', 'csa', 'eco', 'ftu', 'hch', 'hdu', 'hha', 'ilo', 'lpn', 'maq', 'mca', 'msu', 'noc', 'pat', 'pcr', 'pfl', 'pha', 'pin', 'plu', 'ppr', 'rma', 'saz', 'sde', 'sdn', 'shm', 'tcx', 'vfi', 'vvu', 'xca']
            organisms = representativeOrganisms
            output.append( 'Representative:' )
            
        elif i == 2:
            
            #-     2. get group of organisms 'Gammaproteobacteria'
            organisms = taxonomy.getOrganismAbbreviationsByPath('Gammaproteobacteria', oneOrganismPerSpecies=False)
            output.append( '\nGammaproteobacteria:' )
            
        elif i == 3:
            
            #-     3. get group of organisms 'Gammaproteobacteria', excluding unclassified
            organisms = taxonomy.getOrganismAbbreviationsByPath('Gammaproteobacteria', exceptPaths='unclassified', oneOrganismPerSpecies=False)
            output.append( '\nGammaproteobacteria without unclassified:' )
        
        group = Organism.Group( organisms )
        
        #-     REPEAT for varying majority-percentages:
        for percentage in [100, 90, 80, 70, 60, 50, 40, 30, 20, 10 , 1]:
        
            #-         calculate EC numbers occuring in group's core metabolism
            ourECnumbers = group.majorityEcGraph(majorityPercentage = percentage, noMultifunctional = False).getECs()
            ourECnumbers = EcNumber.removeWildcards(ourECnumbers)
            
            #-         overlap Poot-Hernandez' set with ours and print amount of EC numbers inside the intersection and falling off either side
            onlyInTheirs = theirECnumbers.difference( ourECnumbers ) 
            inBoth = theirECnumbers.intersection( ourECnumbers )
            onlyInOurs = ourECnumbers.difference( theirECnumbers )
            
            output.append(str(percentage) + '%:\t' + str(len(onlyInTheirs)) + '\t' + str(len(inBoth)) + '\t' + str(len(onlyInOurs)) )
        
    for line in output:
        print(line)
