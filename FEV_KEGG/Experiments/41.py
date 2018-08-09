"""
Question
--------
When comparing the core metabolism of Archaea and Bacteria, what differences and similarities occur?

Method
------
- build Bacteria clade
- build Archaea clade
- REPEAT for varying majority-percentages:
-     overlap core metabolisms and print amount of EC numbers inside the intersection and falling off either side
-     remove wildcard EC numbers
- END
- build clade pair
- export unified metabolism, coloured by only-Archaea/both/only-Bacteria

Result
------

::

    Maj. %    Bacteria    both    Archaea
    100%:    0    1    7
    90%:    50    40    36
    80%:    83    67    47
    70%:    103    103    59
    60%:    125    129    58
    50%:    153    163    72
    40%:    191    192    75
    30%:    235    229    90
    20%:    304    279    83
    10%:    400    386    87
    1%:    631    653    91

    See bacteria_vs_archaea.jpeg

Conclusion
----------
Bacteria and Archaea always share a significant amount of EC numbers, but never all of them. The much bigger group of Bacteria also has many more EC numbers which never occur in Archaea.
This might be because there are more known Bacteria organisms than known Archaea organisms, i.e. a statistical skew. Or it might be because Bacteria are, as a group, more versatile in habitat than Archaea.

The exported graph comparing Bacteria and Archaea directly (at 80% majority) shows several regions (more or less complete pathways) which only occur in either of the clades' core metabolisms.
Which does not mean they do not occur in any individual organism of the other clade!
For example:
Only in Bacteria: 00061 Fatty acid biosynthesis and 00550 Peptidoglycan biosynthesis
Only in Archaea: 00790 Folate biosynthesis and 00900 Terpenoid backbone biosynthesis

Apart from these regions standing out, both clades seem to have evolved different ways of providing redundancy to their common metabolism.
"""
from FEV_KEGG.Drawing import Export
from FEV_KEGG.KEGG.File import cache
from FEV_KEGG.Evolution.Clade import CladePair, Clade
from FEV_KEGG.Graph.Elements import EcNumber

@cache(folder_path='experiments/41', file_name='bacteria_clade')
def getBacteriaClade():
    bacteriaClade = Clade('/Bacteria')
    # pre-fetch collective metabolism into memory
    bacteriaClade.collectiveMetabolism(excludeMultifunctionalEnzymes=True)
    return bacteriaClade

@cache(folder_path='experiments/41', file_name='archaea_clade')
def getArchaeaClade():
    archaeaClade = Clade('/Archaea')
    # pre-fetch collective metabolism into memory
    archaeaClade.collectiveMetabolism(excludeMultifunctionalEnzymes=True)
    return archaeaClade

if __name__ == '__main__':
    
    output = ['Maj. %\tBacteria\tboth\tArchaea']
    
    #- build Bacteria clade
    bacteriaClade = getBacteriaClade()
    #- build Archaea clade
    archaeaClade = getArchaeaClade()
    
    #- REPEAT for varying majority-percentages:
    for percentage in [100, 90, 80, 70, 60, 50, 40, 30, 20, 10 , 1]:
        
        #-     overlap core metabolisms and print amount of EC numbers inside the intersection and falling off either side
        bacteriaECs = bacteriaClade.coreMetabolism(percentage).getECs()
        archaeaECs = archaeaClade.coreMetabolism(percentage).getECs()
        bothECs = bacteriaECs.intersection(archaeaECs)
        onlyBacteriaECs = bacteriaECs.difference(archaeaECs)
        onlyArchaeaECs = archaeaECs.difference(bacteriaECs)
        
        #-     remove wildcard EC numbers
        onlyBacteriaECs = EcNumber.removeWildcards(onlyBacteriaECs)
        bothECs = EcNumber.removeWildcards(bothECs)
        onlyArchaeaECs = EcNumber.removeWildcards(onlyArchaeaECs)
        
        output.append( str(percentage) + '%:\t' + str(len(onlyBacteriaECs)) + '\t' + str(len(bothECs)) + '\t' + str(len(onlyArchaeaECs)) )
    
    for line in output:
        print(line)
    
    #- build clade pair
    cladePair = CladePair(bacteriaClade, archaeaClade)
    #- export unified metabolism, coloured by only Archaea/both/only Bacteria
    unifiedEcGraph = cladePair.unifiedMetabolism(colour = True)
    Export.forCytoscape(unifiedEcGraph, file = 'experiments/41/bacteria_vs_archaea', inCacheFolder=True)
