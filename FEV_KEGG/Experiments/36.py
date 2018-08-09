"""
Context
-------
Our method of creating a LUCA was introduced in :mod:`35`. Now it needs to be compared with other approaches, e.g. with Goldman et al.
Beware, though, that this experiment is limited to Archaea and Bacteria, while Goldman is not!
We have to reduce our resulting set of EC numbers to three levels, to be able to compare with Goldman's set.

Question
--------
Which EC numbers from our LUCA approach match the approach of Goldman et al?

Method
------
- get Goldman-LUCA
- REPEAT for each clade in Archaea, Bacteria, Archaea+Bacteria
-     get collective metabolism of clade
-     REPEAT for varying majority percentages
-         calculate core metabolism
-         reduce EC numbers to first three levels
-         overlap Goldman's set with ours and print amount of EC numbers inside the intersection and falling off either side

Result
------

::

    Archaea-CoreLUCA
    100%:    9    1    3
    90%:    6    4    37
    80%:    6    4    45
    70%:    5    5    53
    60%:    4    6    57
    50%:    3    7    62
    40%:    1    9    69
    30%:    1    9    78
    20%:    1    9    88
    10%:    0    10    104
    
    
    Bacteria-CoreLUCA
    100%:    9    1    0
    90%:    4    6    30
    80%:    3    7    45
    70%:    2    8    55
    60%:    1    9    63
    50%:    1    9    73
    40%:    0    10    83
    30%:    0    10    96
    20%:    0    10    107
    10%:    0    10    125
    
    
    Archaea-Bacteria-CoreLUCA
    100%:    9    1    0
    90%:    6    4    26
    80%:    4    6    42
    70%:    2    8    56
    60%:    1    9    61
    50%:    1    9    71
    40%:    0    10    83
    30%:    0    10    96
    20%:    0    10    105
    10%:    0    10    123

Conclusion
----------
At the highest majority percentages, there is little overlap between Goldman's LUCA and ours. Around 80%, however, there are more shared ECs than mismatched.
Starting with 40% majority, almost all ECs found in Goldman's LUCA are also in ours.
Lastly, at  10% and below, all of Goldman's ECs are included in our LUCA.

However, our LUCA is obviously much bigger than Goldman's. This would probably not be much different when including Eukaryota.
It might be because Goldman's sources carefully compiled their list of 'essential' ECs by hand, implying that our list includes 'common' but not 'essential' ECs.
This is a question of biochemistry and definition outside the scope of this work.

All this seems to prove that our approach is not wrong, albeit only moderately consistent with Goldman's analysis, and generally more generous about the definition of 'essential'.
"""
from FEV_KEGG.Evolution.LUCA import CoreLUCA, GoldmanLUCA
from FEV_KEGG.Graph.Elements import EcNumber

if __name__ == '__main__':
    
    output = []
    
    #- get Goldman-LUCA
    goldmanLucaEcNumbers = GoldmanLUCA().ecNumbers
    
    #- REPEAT for each clade in Archaea, Bacteria, Archaea+Bacteria
    for clade in [CoreLUCA.CladeType.archaea, CoreLUCA.CladeType.bacteria, CoreLUCA.CladeType.archaeaBacteria]:
        
        #-     get collective metabolism of clade
        coreLuca = CoreLUCA(clade)
        output.append(coreLuca.nameAbbreviation)
        
        #-     REPEAT for varying majority percentages
        for percentage in range(100, 0, -10):
            
            #-         calculate core metabolism
            ourECnumbers = coreLuca.substanceEcGraph( percentage ).getECs()
            
            #-         reduce EC numbers to first three levels
            ourECnumbers = EcNumber.insertWildcards(ourECnumbers, keepLevels = 3, allowHigherWildcards = False, returnSet = True)
            
            onlyInTheirs = goldmanLucaEcNumbers.difference( ourECnumbers ) 
            inBoth = goldmanLucaEcNumbers.intersection( ourECnumbers )
            onlyInOurs = ourECnumbers.difference( goldmanLucaEcNumbers )
            
            output.append(str(percentage) + '%:\t' + str(len(onlyInTheirs)) + '\t' + str(len(inBoth)) + '\t' + str(len(onlyInOurs)) )
            
        output.append("\n")
    
    for line in output:
        print(line)
