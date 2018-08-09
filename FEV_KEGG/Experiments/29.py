"""
Context
-------
See experiment :mod:`28` conclusion.

Question
--------
Do the missing 19 ECs even exist in today's Escherichia coli MG1655 (eco)?

Method
------
- extract EC numbers from Almaas et al. (2005) by hand
- get organism 'Escherichia coli MG1655' (eco) from current KEGG
- calculate EC numbers occuring in eco's core metabolism
- overlap Almaas' set with ours and print amount of EC numbers inside the intersection and falling off either side

Result
------

::

    others    both    ours
    19        43      530

Conclusion
----------
As in experiment :mod:`28`, there are still 19 ECs from 'essential reactions' missing from today's Escherichia coli MG1655, even though they were obviously 
present in the year 2000. The 19 ECs not contained in our approach of a core metabolism are, therefore, to be attributed to a different set of base data, 
and thus not necessarily to a flaw in either approach.

All in all, there is no contradiction between both approaches towards a core metabolism. To some extend, the findings could be interpreted as a partial 
validation of our approach.
"""
import FEV_KEGG.KEGG.Organism as Organism
from FEV_KEGG.Graph.Elements import EcNumber


if __name__ == '__main__':
    
    output = ['others\tboth\tours']
    
    #- extract EC numbers from Almaas et al. (2005) by hand
    theirECnumberStrings = ['1.1.1.158', '1.1.1.193', '1.1.1.25', '1.5.1.3', '1.6.4.5', '1.7.99.5', '2.2.1.2', '2.3.1.129', '2.3.1.39', '2.4.1.182', '2.4.1.21', '2.5.1.15', '2.5.1.19', '2.5.1.7', '2.5.1.9', '2.6.1.16', '2.7.1.107', '2.7.1.130', '2.7.1.23', '2.7.1.24', '2.7.1.26', '2.7.1.33', '2.7.2.3', '2.7.4.6', '2.7.4.8', '2.7.4.9', '2.7.6.3', '2.7.7.18', '2.7.7.2', '2.7.7.23', '2.7.7.27', '2.7.7.3', '2.7.7.38', '2.7.7.41', '2.7.8.5', '2.7.8.8', '3.1.3.45', '3.5.4.16', '3.5.4.25', '3.5.4.26', '3.6.1.1', '3.6.1.34', '3.6.1.45', '4.1.1.36', '4.1.1.65', '4.1.2.13', '4.1.2.16', '4.1.2.25', '4.2.1.10', '4.2.1.11', '4.6.1.3', '4.6.1.4', '5.1.1.3', '5.3.1.1', '5.3.1.13', '6.3.2.12', '6.3.2.13', '6.3.2.15', '6.3.2.4', '6.3.2.5', '6.3.2.8', '6.3.2.9']
    theirECnumbers = set()
    for string in theirECnumberStrings:
        theirECnumbers.add( EcNumber(string) )
    
    #- get group of organisms 'Escherichia coli'
    eco = Organism.Organism( 'eco' )
    
    
    #- calculate EC numbers occuring in eco's core metabolism
    ourECnumbersWithWildcard = eco.substanceEcGraph(noMultifunctional = True).getECs()
    ourECnumbers = EcNumber.removeWildcards(ourECnumbersWithWildcard)
    
    #- overlap Almaas' set with ours and print amount of EC numbers inside the intersection and falling off either side
    onlyInTheirs = theirECnumbers.difference( ourECnumbers ) 
    inBoth = theirECnumbers.intersection( ourECnumbers )
    onlyInOurs = ourECnumbers.difference( theirECnumbers )
        
    output.append(str(len(onlyInTheirs)) + '\t' + str(len(inBoth)) + '\t' + str(len(onlyInOurs)) )
        
    for line in output:
        print(line)