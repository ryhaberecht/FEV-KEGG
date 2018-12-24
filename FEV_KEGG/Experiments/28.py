"""
Context
-------
The approach of building a consensus/majority graph of enzymes/EC numbers to find a core metabolism shared among several organisms has to be validated against previous research.
One such previous research deals with E. coli, but uses a different approach, asking not what the core metabolism 'can do', but what it 'always does'.
Almaas et al. (2005) list core reactions calculated via flux analysis in table S1 (https://doi.org/10.1371/journal.pcbi.0010068.st001), some of which have annotated EC numbers.
These EC numbers are used to validate the approach of this library. Multifunctional enzymes and EC numbers containing wildcards (e.g. 1.2.-.-) are excluded on both sides, to minimise statistical skew.
This leaves 62 EC numbers in Almaas' approach.

Question
--------
Does the consensus/majority graph approach to core metabolism yield a similar set of EC numbers as the approach of Almaas et al. (2005)?

Method
------
- extract EC numbers from Almaas et al. (2005) by hand
- get group of organisms 'Escherichia coli'
- REPEAT for varying majority-percentages:
-    calculate EC numbers occuring in group's core metabolism
-    overlap Almaas' set with ours and print amount of EC numbers inside the intersection and falling off either side

Result
------

::

    Maj. %   others    both    ours
    100%:    28        34      381
     90%:    19        43      491
     80%:    19        43      499
     70%:    19        43      510
     60%:    19        43      518
     50%:    19        43      522
     40%:    19        43      531
     30%:    19        43      542
     20%:    19        43      550
     10%:    19        43      564
      1%:    19        43      602

Conclusion
----------
With a 90% majority and below, the number of overlapping ECs does not increase any more. This indicates that, at least for E. coli, a 90% majority is enough 
to create a stable core metabolism, diminishing the skew excerted by unusually specialised organisms.
In the case of E. coli these could be soil-based E. coli strains, which remains to be researched.

About 69% of the ECs in the reaction-based core metabolism, as postulated by Almaas et al., are also included in the majority-based core metabolism of 
our approach. Due to some ECs missing in Almaas' table S1, this percentage could have been even bigger.
This substantial overlap shows most essential reactions are also covered by a majority approach.
However, this goes along with two interesting observations:

1) 31% of essential reactions are not included in any majority, not even in a single organism from KEGG at 1% majority (effectively n=1).
    This could be because of a flaw in either approach, or because the data Almaas et al. use stems from the year 2000 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC25862/)
    and it might be that then Escherichia coli MG1655 was said to include different ECs than in today's KEGG database. This has to be investigated.

2) Only 8% of majority ECs (at 90%) are essential reactions. This indicates that while E. coli organisms share many ECs, most of them are only active at 
    certain times.
"""
from FEV_KEGG.Evolution.Taxonomy import NCBI
import FEV_KEGG.KEGG.Organism as Organism
from FEV_KEGG.Graph.Elements import EcNumber


if __name__ == '__main__':
    
    output = ['Maj. %\tothers\tboth\tours']
    
    #- extract EC numbers from Almaas et al. (2005) by hand
    theirECnumberStrings = ['1.1.1.158', '1.1.1.193', '1.1.1.25', '1.5.1.3', '1.6.4.5', '1.7.99.5', '2.2.1.2', '2.3.1.129', '2.3.1.39', '2.4.1.182', '2.4.1.21', '2.5.1.15', '2.5.1.19', '2.5.1.7', '2.5.1.9', '2.6.1.16', '2.7.1.107', '2.7.1.130', '2.7.1.23', '2.7.1.24', '2.7.1.26', '2.7.1.33', '2.7.2.3', '2.7.4.6', '2.7.4.8', '2.7.4.9', '2.7.6.3', '2.7.7.18', '2.7.7.2', '2.7.7.23', '2.7.7.27', '2.7.7.3', '2.7.7.38', '2.7.7.41', '2.7.8.5', '2.7.8.8', '3.1.3.45', '3.5.4.16', '3.5.4.25', '3.5.4.26', '3.6.1.1', '3.6.1.34', '3.6.1.45', '4.1.1.36', '4.1.1.65', '4.1.2.13', '4.1.2.16', '4.1.2.25', '4.2.1.10', '4.2.1.11', '4.6.1.3', '4.6.1.4', '5.1.1.3', '5.3.1.1', '5.3.1.13', '6.3.2.12', '6.3.2.13', '6.3.2.15', '6.3.2.4', '6.3.2.5', '6.3.2.8', '6.3.2.9']
    theirECnumbers = set()
    for string in theirECnumberStrings:
        theirECnumbers.add( EcNumber(string) )
    
    #- get group of organisms 'Escherichia coli'
    taxonomy = NCBI.getTaxonomy()
    group = Organism.Group( taxonomy.getOrganismAbbreviationsByPath('Escherichia coli', oneOrganismPerSpecies=False) )
    
    #- REPEAT for varying majority-percentages:
    for percentage in [100, 90, 80, 70, 60, 50, 40, 30, 20, 10 , 1]:
    
    #-    calculate EC numbers occuring in group's core metabolism
        ourECnumbersWithWildcard = group.majorityEcGraph(majorityPercentage = percentage, noMultifunctional = True).getECs()
        ourECnumbers = EcNumber.removeWildcards(ourECnumbersWithWildcard)
    
    #-    overlap Almaas' set with ours and print amount of EC numbers inside the intersection and falling off either side
        onlyInTheirs = theirECnumbers.difference( ourECnumbers ) 
        inBoth = theirECnumbers.intersection( ourECnumbers )
        onlyInOurs = ourECnumbers.difference( theirECnumbers )
        
        output.append(str(percentage) + '%:\t' + str(len(onlyInTheirs)) + '\t' + str(len(inBoth)) + '\t' + str(len(onlyInOurs)) )
        
    for line in output:
        print(line)