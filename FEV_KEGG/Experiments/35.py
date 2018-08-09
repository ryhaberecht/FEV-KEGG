"""
Context
-------
LUCA has previously been derived from various methods. This is a method of creating a LUCA using only annotated genomes from KEGG. It is expected to be far less precise than methods carefully crafted by experts.
Furthermore, this experiment is limited to the common ancestor of Archaea, Bacteria, and Archaea+Bacteria; Eukaryota are not in the scope of this work and, thus, the truley universal LUCA is not.

Question
--------
How many EC numbers are shared among Archaea, Bacteria, and Archaea+Bacteria?
Including the ones from multifunctional enzymes.

Method
------
- REPEAT for each clade in Archaea, Bacteria, Archaea+Bacteria
-     get collective metabolism of clade
-     REPEAT for varying majority percentages
-         calculate core metabolism
-         print number of ECs by majority percentage
-         IF majority percentage is 100%, i.e. consensus, also print the ECs

Result
------

::

    Archaea-CoreLUCA
    100% ECs: 8 (6.1.1.2, 6.3.5.7, 6.1.1.20, 6.1.1.4, 2.7.7.6, 2.7.7.7, 6.1.1.21, 2.7.4.3) ->  1% of collective
    90% ECs: 101 ->  9% of collective
    80% ECs: 152 -> 14% of collective
    70% ECs: 207 -> 19% of collective
    60% ECs: 241 -> 23% of collective
    50% ECs: 293 -> 27% of collective
    40% ECs: 328 -> 31% of collective
    30% ECs: 400 -> 37% of collective
    20% ECs: 460 -> 43% of collective
    10% ECs: 603 -> 56% of collective
    
    Bacteria-CoreLUCA
    100% ECs: 1 (2.7.7.7) ->  0% of collective
    90% ECs: 86 ->  4% of collective
    80% ECs: 179 ->  9% of collective
    70% ECs: 262 -> 12% of collective
    60% ECs: 329 -> 16% of collective
    50% ECs: 394 -> 19% of collective
    40% ECs: 472 -> 22% of collective
    30% ECs: 582 -> 28% of collective
    20% ECs: 720 -> 34% of collective
    10% ECs: 969 -> 46% of collective
    
    Archaea-Bacteria-CoreLUCA
    100% ECs: 1 (2.7.7.7) ->  0% of collective
    90% ECs: 70 ->  3% of collective
    80% ECs: 164 ->  8% of collective
    70% ECs: 255 -> 12% of collective
    60% ECs: 321 -> 15% of collective
    50% ECs: 388 -> 18% of collective
    40% ECs: 466 -> 22% of collective
    30% ECs: 571 -> 27% of collective
    20% ECs: 708 -> 33% of collective
    10% ECs: 965 -> 45% of collective

Conclusion
----------
There are many EC numbers shared among all species of the Archaea, Bacteria, and even the combination of both clades.
Even when defining 'core metabolism' by consensus, 2.7.7.7 occurs in all organisms.

However, this approach suffers harshly from incomplete annotations, not only resulting in one less count of an EC number which should have occured, but also resulting in a decreased chance of occuring within a fixed majority threshold.
"""
from FEV_KEGG.Evolution.LUCA import CoreLUCA
from FEV_KEGG.Statistics import Percent

if __name__ == '__main__':
    
    output = []
    
    #- REPEAT for each clade in Archaea, Bacteria, Archaea+Bacteria
    for clade in [CoreLUCA.CladeType.archaea, CoreLUCA.CladeType.bacteria, CoreLUCA.CladeType.archaeaBacteria]:
        
        #-     get collective metabolism of clade
        luca = CoreLUCA(clade)
        output.append(luca.nameAbbreviation)
        
        collectiveMetabolismLength = len( luca.collectiveMetabolism().getECs() )
        
        #-     REPEAT for varying majority percentages
        for percentage in range(100, 0, -10):
            
            #-         calculate core metabolism
            ecNumbers = luca.substanceEcGraph( percentage ).getECs()
            
            #-         print number of ECs by majority percentage
            if percentage == 100:
                #-         IF majority percentage is 100%, i.e. consensus, also print the ECs
                output.append( str(percentage) + "% ECs: " + str(len( ecNumbers )) + " (" + ", ".join([str(x) for x in ecNumbers]) + ")"  + ' -> ' + Percent.getPercentStringShort(len( ecNumbers ), collectiveMetabolismLength, 0) + '% of collective')
            else:
                output.append( str(percentage) + "% ECs: " + str(len( ecNumbers )) + ' -> ' + Percent.getPercentStringShort(len( ecNumbers ), collectiveMetabolismLength, 0) + '% of collective')
            
        output.append("\n")
        
        luca.clade.group.freeHeap()
        luca = None
    
    for line in output:
        print(line)
