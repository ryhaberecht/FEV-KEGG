"""
Question
--------
Which EC numbers in Desulfobacterales, come from neofunctionalised enzymes?
How do they differ varying certain parameters of computation?
Which is the set of parameters most suited for investigating neofunctionalisation?

Method
------
- Get the clade
- Print all neofunctionalised EC numbers
- REPEAT for varying parameters of computation:
-     Different gene duplication models: SimpleGeneDuplication / SimpleGroupGeneDuplication
-     Different set of allowed gene duplicates: restricted to enzymes within the core metabolism / unrestricted
-     Different neofunctionalisation-majority values. This is the number of organisms (in percent of the group) required to have a neofunctionalisation involving an EC number fort this EC number to be reported as "neofunctionalised".

Result
------

::
    neofunctionalisation-majority: 70%    
    
    SimpleGeneDuplication + unrestricted: 29 (62%)
    1.17.1.9
    1.3.5.4
    1.4.3.16
    2.1.3.2
    2.1.3.3
    2.5.1.1
    2.5.1.10
    2.5.1.29
    2.5.1.90
    2.6.1.11
    2.6.1.17
    2.6.1.85
    4.1.3.27
    4.2.1.46
    5.1.3.2
    5.1.3.6
    5.4.2.10
    5.4.2.8
    5.4.3.8
    6.1.1.4
    6.1.1.5
    6.1.1.6
    6.1.1.9
    6.2.1.-
    6.2.1.1
    6.2.1.3
    6.3.1.-
    6.3.2.45
    6.3.2.8
    
    SimpleGeneDuplication + restriction: 7 (15%)
    2.1.3.2
    2.1.3.3
    4.2.1.46
    5.1.3.6
    6.1.1.4
    6.1.1.5
    6.1.1.9
    
    SimpleGeneDuplication + unrestricted: 56 (119%)
    1.1.1.205
    1.1.1.23
    1.1.1.308
    1.1.1.42
    1.1.1.85
    1.3.5.1
    1.3.5.4
    1.4.1.1
    1.4.3.16
    1.5.1.20
    1.6.1.2
    1.7.1.7
    2.1.1.10
    2.1.2.1
    2.1.2.10
    2.1.3.2
    2.1.3.3
    2.2.1.1
    2.2.1.7
    2.3.1.1
    2.5.1.1
    2.5.1.10
    2.5.1.29
    2.5.1.47
    2.5.1.54
    2.5.1.55
    2.5.1.90
    2.6.1.-
    2.6.1.83
    2.6.1.85
    2.6.1.9
    2.7.2.8
    3.5.2.2
    3.5.2.3
    4.1.1.3
    4.1.1.48
    4.1.1.81
    4.1.3.27
    4.2.1.12
    4.2.1.46
    4.2.1.9
    5.1.3.2
    5.1.3.6
    5.3.1.24
    5.4.2.10
    5.4.2.8
    5.4.99.18
    6.1.1.16
    6.1.1.17
    6.1.1.18
    6.1.1.4
    6.1.1.5
    6.1.1.6
    6.1.1.9
    6.3.4.13
    6.4.1.1

    
    SimpleGroupGeneDuplication + restriction: 29 (62%)
    1.1.1.42
    1.1.1.85
    2.1.3.2
    2.1.3.3
    2.3.1.1
    2.5.1.47
    2.5.1.54
    2.5.1.55
    2.6.1.1
    2.6.1.83
    2.6.1.9
    2.7.2.8
    2.7.9.1
    2.7.9.2
    4.1.1.3
    4.2.1.46
    5.1.3.2
    5.1.3.6
    6.1.1.12
    6.1.1.16
    6.1.1.17
    6.1.1.18
    6.1.1.4
    6.1.1.5
    6.1.1.6
    6.1.1.9
    6.2.1.1
    6.2.1.3
    6.4.1.1
    
    
    
    
    neofunctionalisation-majority: 20%
        
    SimpleGeneDuplication + restriction: 25 (53%)
    1.1.1.42
    1.1.1.85
    2.1.3.2
    2.1.3.3
    2.4.2.14
    2.5.1.47
    2.5.1.54
    2.5.1.55
    2.6.1.1
    2.6.1.16
    2.6.1.83
    2.6.1.9
    4.2.1.46
    5.1.3.2
    5.1.3.6
    6.1.1.12
    6.1.1.16
    6.1.1.4
    6.1.1.5
    6.1.1.6
    6.1.1.9
    6.2.1.1
    6.2.1.3
    6.3.2.8
    6.3.2.9

    
    SimpleGroupGeneDuplication + restriction: 34 (72%)
    1.1.1.42
    1.1.1.85
    2.1.3.2
    2.1.3.3
    2.3.1.1
    2.4.2.14
    2.5.1.47
    2.5.1.54
    2.5.1.55
    2.6.1.1
    2.6.1.16
    2.6.1.83
    2.6.1.9
    2.7.2.8
    2.7.9.1
    2.7.9.2
    4.1.1.3
    4.2.1.46
    5.1.3.2
    5.1.3.6
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
    6.3.2.8
    6.3.2.9
    6.4.1.1
    
    
    
    
    neofunctionalisation-majority: 0%
        
    SimpleGeneDuplication + restriction: 30 (64%)
    1.1.1.42
    1.1.1.85
    2.1.3.2
    2.1.3.3
    2.3.1.1
    2.4.2.14
    2.5.1.47
    2.5.1.54
    2.5.1.55
    2.6.1.1
    2.6.1.16
    2.6.1.83
    2.6.1.9
    2.7.2.8
    2.7.9.1
    2.7.9.2
    4.2.1.46
    5.1.3.2
    5.1.3.6
    6.1.1.12
    6.1.1.16
    6.1.1.22
    6.1.1.4
    6.1.1.5
    6.1.1.6
    6.1.1.9
    6.2.1.1
    6.2.1.3
    6.3.2.8
    6.3.2.9

Conclusion
----------
Restricting the set of allowed gene duplicates to enzymes within the core metabolism helps to keep secondary metabolism out, when this is wanted. As we usually handle core metabolisms, we certainly want to keep ties to non-core metabolism out of the results.

One problem could be too few cases of neofunctionalisation, which can be avoided using SimpleGeneDuplication with a low neofunctionalisation-majority value.

SimpleGroupGeneDuplication is far too slow for groups of significant size, even for the small group of Desulfobacterales (9 organisms), the download of all gene matching takes about one hour. With the increase of group size, the time to download grows almost exponentially.

In summary, we will only be using SimpleGeneDuplication + restriction for investigating neofunctionalisation.
"""

from FEV_KEGG.KEGG.File import cache
from FEV_KEGG.Evolution.Clade import Clade

@cache(folder_path='experiments/44', file_name='desulfobacterales_clade')
def getDesulfobacteralesClade():
    archaeaClade = Clade('Deltaproteobacteria/Desulfobacterales')
    # pre-fetch collective metabolism into memory
    archaeaClade.collectiveMetabolism(excludeMultifunctionalEnzymes=True)
    # pre-fetch collective enzyme metabolism into memory
    archaeaClade.collectiveMetabolismEnzymes(excludeMultifunctionalEnzymes=True)
    return archaeaClade

if __name__ == '__main__':

    #- Get the clade
    desulfoClade = getDesulfobacteralesClade()
    
    #- Print all neofunctionalised EC numbers
    neofunctionalisedMetabolismSet = desulfoClade.neofunctionalisedECs().getECs()
    
    sortedList = list(neofunctionalisedMetabolismSet)
    sortedList.sort()
    
    for ecNumber in sortedList:
        print( ecNumber )    