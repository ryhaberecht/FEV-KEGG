"""
Question
--------
Which EC numbers are linked with neofunctionalisations in Deltaproteobacteria and Spirochaetes?
How do the two groups compare?

Method
------
- get both clades
- calculate set of "neofunctionalised" ECs for Alphaproteobacteria
- calculate set of "neofunctionalised" ECs for Betaproteobacteria
- REPEAT for varying core metabolism majority-percentages:
-     overlap sets and print amount of EC numbers inside the intersection and falling off either side

Result
------

::
    
    core metabolism majority: 80%
    neofunctionalisation majority: 0% (this means that gene duplication within a single organism is enough)
    
    Deltaproteobacteria:
    
    core metabolism ECs: 228
    
    "neofunctionalised" ECs: 36 (16%)
    1.1.1.42
    1.1.1.85
    2.1.3.2
    2.1.3.3
    2.2.1.1
    2.2.1.7
    2.4.2.14
    2.5.1.47
    2.6.1.1
    2.6.1.16
    2.6.1.62
    2.6.1.9
    4.1.3.-
    4.2.1.46
    5.1.1.1
    5.1.3.2
    5.1.3.6
    5.3.1.16
    5.4.3.8
    5.4.99.18
    6.1.1.12
    6.1.1.16
    6.1.1.17
    6.1.1.18
    6.1.1.22
    6.1.1.4
    6.1.1.5
    6.1.1.6
    6.1.1.9
    6.2.1.1
    6.2.1.3
    6.3.2.10
    6.3.2.13
    6.3.2.8
    6.3.2.9
    6.3.4.13
    
    
    Spirochaetes:
    
    core metabolism ECs: 80
    
    "neofunctionalised" ECs: 10 (12%)
    6.1.1.10
    6.1.1.12
    6.1.1.20
    6.1.1.22
    6.1.1.4
    6.1.1.5
    6.1.1.6
    6.1.1.9
    6.3.2.13
    6.3.2.8
    
    
    Comparison:
    100%:   0       2       3
    80%:    28      8       2
    60%:    47      8       4
    40%:    72      20      8
    20%:    125     40      14

Conclusion
----------
Deltaproteobacteria consistently have more "neofunctionalised" EC numbers, which is not surprising, as the core metabolism is much bigger.
The same effect explains the constant increase in "neofunctionalised" EC numbers with growing core metabolism, due to lower majority-percentage.
This indicates that neofunctionalisation is widespread and not limited to the most important parts of metabolism. However, there might still be a statistical which remains to be investigated.

Occurence in these results above does not mean, that those ECs are actually "neofunctionalised" in all of the organisms! Because in this experiment, a single neofunctionalisation is enough to mark the associated ECs as "neofunctionalised".
It might be interesting to increase the neofunctionalisation majority percentage, to only include ECs which have neofunctionalisations in x% of the clade's organisms. 
"""

from FEV_KEGG.KEGG.File import cache
from FEV_KEGG.Evolution.Clade import Clade
from FEV_KEGG.Statistics import Percent

@cache(folder_path='experiments', file_name='deltaproteobacteria_clade')
def getCladeA():
    clade = Clade('Deltaproteobacteria')
    # pre-fetch collective metabolism into memory
    clade.collectiveMetabolism(excludeMultifunctionalEnzymes=True)
    # pre-fetch collective enzyme metabolism into memory
    clade.collectiveMetabolismEnzymes(excludeMultifunctionalEnzymes=True)
    return clade

@cache(folder_path='experiments', file_name='spirochaetes_clade')
def getCladeB():
    clade = Clade('Spirochaetes')
    # pre-fetch collective metabolism into memory
    clade.collectiveMetabolism(excludeMultifunctionalEnzymes=True)
    # pre-fetch collective enzyme metabolism into memory
    clade.collectiveMetabolismEnzymes(excludeMultifunctionalEnzymes=True)
    return clade

if __name__ == '__main__':

    output = []
    
    cladeA = getCladeA()
    cladeB = getCladeB()
    
    majorityPercentageCoreMetabolism = 80
    majorityPercentageNeofunctionalisation = 0    
    output.append( 'core metabolism majority: ' + str(majorityPercentageCoreMetabolism) + '%' )
    output.append( 'neofunctionalisation majority: ' + str(majorityPercentageNeofunctionalisation) + '% (this means that gene duplication within a single organism is enough)' )
    output.append('')
     
     
    # Clade A
    output.append(', '.join(cladeA.ncbiNames) + ':')
    output.append('')
    cladeAEcCount = len(cladeA.coreMetabolism(majorityPercentageCoreMetabolism).getECs())
    output.append( 'core metabolism ECs: ' + str(cladeAEcCount) )
    output.append('')
      
    cladeANeofunctionalisedMetabolismSet = cladeA.neofunctionalisedECs(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation).getECs()
        
    cladeASortedList = list(cladeANeofunctionalisedMetabolismSet)
    cladeASortedList.sort()
       
    output.append( '"neofunctionalised" ECs: ' + str(len(cladeASortedList)) + ' (' + str(Percent.getPercentStringShort(len(cladeASortedList), cladeAEcCount, 0)) + '%)' )
    for ecNumber in cladeASortedList:
        output.append( ecNumber )
     
     
     
     
     
    # Clade B
    output.append('')
    output.append('')
    output.append(', '.join(cladeB.ncbiNames) + ':')
    output.append('')
    cladeBEcCount = len(cladeB.coreMetabolism(majorityPercentageCoreMetabolism).getECs())
    output.append( 'core metabolism ECs: ' + str(cladeBEcCount) )
    output.append('')
     
    cladeBNeofunctionalisedMetabolismSet = cladeB.neofunctionalisedECs(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation).getECs()
       
    cladeBSortedList = list(cladeBNeofunctionalisedMetabolismSet)
    cladeBSortedList.sort()
      
    output.append( '"neofunctionalised" ECs: ' + str(len(cladeBSortedList)) + ' (' + str(Percent.getPercentStringShort(len(cladeBSortedList), cladeBEcCount, 0)) + '%)' )
    for ecNumber in cladeBSortedList:
        output.append( ecNumber )
    
    
    
    
    # Clade comparison
    output.append('')
    output.append('')
    output.append('Comparison:')
    
    for percentage in [100, 80, 60, 40, 20]:
        
        aECs = cladeA.neofunctionalisedECs(percentage).getECs()
        bECs = cladeB.neofunctionalisedECs(percentage).getECs()
        bothECs = aECs.intersection(bECs)
        onlyAECs = aECs.difference(bECs)
        onlyBECs = bECs.difference(aECs)
        
        output.append( str(percentage) + '%:\t' + str(len(onlyAECs)) + '\t' + str(len(bothECs)) + '\t' + str(len(onlyBECs)) )
    
    for line in output:
        print( line )
