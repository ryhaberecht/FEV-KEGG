"""
Question
--------
Which "neofunctionalised" EC numbers exist in Deltaproteobacteria?
Does their existence increase redundancy? If yes, how much?

Method
------
- get Deltaproteobacteria
- get "neofunctionalised" ECs
- calculate redundancy
- REPEAT for each "neofunctionalised" EC
-     print number of ECs the "neofunctionalised" EC provides redundancy for (robustness and flexibility)
- print ECs which have "neofunctionalised" ECs contributing to their robustness

Result
------

::
    core metabolism majority: 80%
    neofunctionalisation majority: 0% (this means that gene duplication within a single organism is enough)
    
    Deltaproteobacteria:
    
    core metabolism ECs: 228
    
    "neofunctionalised" ECs: 36 (16%)
    1.1.1.42        (0 target-flexibility) (0 robustness)
    1.1.1.85        (0 target-flexibility) (0 robustness)
    2.1.3.2         (0 target-flexibility) (0 robustness)
    2.1.3.3         (1 target-flexibility) (0 robustness)
    2.2.1.1         (6 target-flexibility) (2 robustness)
    2.2.1.7         (7 target-flexibility) (2 robustness)
    2.4.2.14        (2 target-flexibility) (0 robustness)
    2.5.1.47        (0 target-flexibility) (0 robustness)
    2.6.1.1         (2 target-flexibility) (0 robustness)
    2.6.1.16        (2 target-flexibility) (0 robustness)
    2.6.1.62        (0 target-flexibility) (0 robustness)
    2.6.1.9         (0 target-flexibility) (0 robustness)
    4.1.3.-         (1 target-flexibility) (1 robustness)
    4.2.1.46        (0 target-flexibility) (0 robustness)
    5.1.1.1         (1 target-flexibility) (0 robustness)
    5.1.3.2         (2 target-flexibility) (0 robustness)
    5.1.3.6         (1 target-flexibility) (0 robustness)
    5.3.1.16        (0 target-flexibility) (0 robustness)
    5.4.3.8         (0 target-flexibility) (0 robustness)
    5.4.99.18       (1 target-flexibility) (0 robustness)
    6.1.1.12        (1 target-flexibility) (0 robustness)
    6.1.1.16        (0 target-flexibility) (0 robustness)
    6.1.1.17        (0 target-flexibility) (0 robustness)
    6.1.1.18        (1 target-flexibility) (0 robustness)
    6.1.1.22        (0 target-flexibility) (0 robustness)
    6.1.1.4         (0 target-flexibility) (0 robustness)
    6.1.1.5         (0 target-flexibility) (0 robustness)
    6.1.1.6         (0 target-flexibility) (0 robustness)
    6.1.1.9         (0 target-flexibility) (0 robustness)
    6.2.1.1         (3 target-flexibility) (0 robustness)
    6.2.1.3         (0 target-flexibility) (0 robustness)
    6.3.2.10        (0 target-flexibility) (0 robustness)
    6.3.2.13        (0 target-flexibility) (0 robustness)
    6.3.2.8         (0 target-flexibility) (0 robustness)
    6.3.2.9         (0 target-flexibility) (0 robustness)
    6.3.4.13        (0 target-flexibility) (0 robustness)
    
    Target-flexibility: 17.1%
    Target-flexibility contribution: 46.2%
    
    Robustness: 7.0%
    Robustness contribution: 31.2%
    
    Robust ECs contributed by "neofunctionalised" ECs:
    [1.2.1.12, 2.2.1.2, 2.4.2.-, 2.7.9.2, 5.1.3.1]
    

Conclusion
----------
31.2% of the core metabolism's robust EC numbers have at least one redundant path, making this EC number robust, which contains a "neofunctionalised" EC number.
7% of the core metabolism's EC numbers are robust.
This makes a total of 7% * 31.2% = 2.18% * 228 ECs in the core metabolism = 5 EC numbers in the core metabolism, which are robust with a "neofunctionalised" EC number on at least one of their alternative paths.
"""

from FEV_KEGG.KEGG.File import cache
from FEV_KEGG.Evolution.Clade import Clade
from FEV_KEGG.Statistics import Percent
from FEV_KEGG.Robustness.Topology.Redundancy import RedundancyType, Redundancy, RedundancyContribution

@cache(folder_path='experiments', file_name='deltaproteobacteria_clade')
def getCladeA():
    clade = Clade('Deltaproteobacteria')
    # pre-fetch collective metabolism into memory
    clade.collectiveMetabolism(excludeMultifunctionalEnzymes=True)
    # pre-fetch collective enzyme metabolism into memory
    clade.collectiveMetabolismEnzymes(excludeMultifunctionalEnzymes=True)
    return clade

if __name__ == '__main__':

    output = ['']

    #- get Deltaproteobacteria
    cladeA = getCladeA()
    majorityPercentageCoreMetabolism = 80
    majorityPercentageNeofunctionalisation = 0
    output.append( 'core metabolism majority: ' + str(majorityPercentageCoreMetabolism) + '%' )
    output.append( 'neofunctionalisation majority: ' + str(majorityPercentageNeofunctionalisation) + '% (this means that gene duplication within a single organism is enough)' )
    output.append('')
    output.append(', '.join(cladeA.ncbiNames) + ':')
    output.append('')
    cladeAEcGraph = cladeA.coreMetabolism(majorityPercentageCoreMetabolism)
    cladeAEcCount = len(cladeAEcGraph.getECs())
    output.append( 'core metabolism ECs: ' + str(cladeAEcCount) )
    output.append('')
    
    #- get "neofunctionalised" ECs
    cladeANeofunctionalisedMetabolismSet = cladeA.neofunctionalisedECs(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation).getECs()
    cladeASortedList = list(cladeANeofunctionalisedMetabolismSet)
    cladeASortedList.sort()
    
    #- calculate redundancy
    cladeARedundancy = Redundancy(cladeAEcGraph)
    cladeARedundancyContribution = RedundancyContribution(cladeARedundancy, cladeANeofunctionalisedMetabolismSet)
    cladeARobustECsForNeofunctionalisedEC = cladeARedundancyContribution.getContributedKeysForSpecial(RedundancyType.ROBUSTNESS)
    cladeAFlexibleECsForNeofunctionalisedEC = cladeARedundancyContribution.getContributedKeysForSpecial(RedundancyType.TARGET_FLEXIBILITY)
      
    #- REPEAT for each "neofunctionalised" EC
    output.append( '"neofunctionalised" ECs: ' + str(len(cladeASortedList)) + ' (' + str(Percent.getPercentStringShort(len(cladeASortedList), cladeAEcCount, 0)) + '%)' )
    for ecNumber in cladeASortedList:
        #-     print number of ECs the "neofunctionalised" EC provides redundancy for (robustness and flexibility)
        output.append( str(ecNumber)  + ' \t(' + str(len(cladeAFlexibleECsForNeofunctionalisedEC.get(ecNumber, []))) + ' target-flexibility) (' + str(len(cladeARobustECsForNeofunctionalisedEC.get(ecNumber, []))) + ' robustness)')
    
    output.append('')
    output.append( 'Target-flexibility: ' + Percent.floatToPercentString(cladeARedundancy.getRedundancyRatio(RedundancyType.TARGET_FLEXIBILITY) ) )
    output.append( 'Target-flexibility contribution: ' + Percent.floatToPercentString(cladeARedundancyContribution.getKeyContributionRatio(RedundancyType.TARGET_FLEXIBILITY) ) )
    output.append('')
    output.append( 'Robustness: ' + Percent.floatToPercentString(cladeARedundancy.getRedundancyRatio(RedundancyType.ROBUSTNESS) ) )
    output.append( 'Robustness contribution: ' + Percent.floatToPercentString(cladeARedundancyContribution.getKeyContributionRatio(RedundancyType.ROBUSTNESS) ) )
    output.append('')
    
    #- print ECs which have "neofunctionalised" ECs contributing to their robustness
    contributedECs = list(cladeARedundancyContribution.getContributingSpecialForKey(RedundancyType.ROBUSTNESS).keys())
    contributedECs.sort()
    output.append( 'Robust ECs contributed by "neofunctionalised" ECs:' )
    output.append( contributedECs ) 
    
    
    
    
    
    
    
    
    for line in output:
        print( line )
