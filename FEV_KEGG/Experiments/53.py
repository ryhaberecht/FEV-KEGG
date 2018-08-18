"""
Context
-------
In :mod:`52`, we did this before, but this time we include partial robustness.

Question
--------
Which neofunctionalised enzymes cause the core metabolism of Deltaproteobacteria to have increased redundancy? How much do they contribute?
Grouped by contributed functions, sorted lexicographically, annotated with links to KEGG, annotated with human-readable names (if possible), and exported into an HTML file.

Method
------
- get clade
- get core metabolism
- calculate "neofunctionalised" ECs
- calculate redundancy
- REPEAT for each "neofunctionalised" EC contributing to redundancy
-     report enzyme pairs of neofunctionalisations, which caused the EC to be considered "neofunctionalised", and are in return contributing to redundancy
-     print them into nice HTML

Result
------

::

    core metabolism majority: 80%
    neofunctionalisation majority: 0% (this means that gene duplication within a single organism is enough)
    
    Deltaproteobacteria:
    
    core metabolism ECs: 228
    
    "neofunctionalised" ECs: 36 (16%)
    
    Neofunctionalisations contributing to robustness: 61
    [see file Deltaproteobacteria_ROBUSTNESS_PARTIAL_Neofunctionalisations-For-Contributed-EC.html]
    
    "neofunctionalised" ECs: 36 (16%)
    
    Neofunctionalisations contributing to robustness: 91
    
    [see file Deltaproteobacteria_ROBUSTNESS_BOTH_Neofunctionalisations-For-Contributed-EC.html]
    

Conclusion
----------
In :mod:`52`, we saw 84 neofunctionalisations contributing to full robustness.
Additionally, 61 neofunctionalisations contribute to partial redundancy.

However, when adding these two sets of neofunctionalisations (ROBUSTNESS_BOTH), they result in only 91 different neofunctionalisations.
This means that some neofunctionalisations contribute to both, full and partial robustness.
This is to be expected, because EC numbers (of which a neofunctionalisation's function change has at least two) can, 
if they are in a central position, easily contribute to the redundancy of multiple other EC numbers.
Some of which may be fully redundant, some only partially redundant.

Note: the well-known neofunctionalisation 2.6.1.1 <-> 2.6.1.9 does not occur in :mod:`52`.
Possible causes were believed to be a high majority percentage, or that both of the two EC numbers can only contribute to partial redundancy.
In Deltaproteobacteria_ROBUSTNESS_PARTIAL_Neofunctionalisations-For-Contributed-EC.html we see seven neofunctionalisations involving 2.6.1.1 and 2.6.1.9, 
therefore, the cause of the non-appearance in :mod:`52` was indeed that they only contribute to partial redundancy.
This shows that it might be wise to choose partial types of redundancy, to receive more interesting results, especially when dealing with incomplete data such as KEGG's GENE and PATHWAY databases.
"""

from FEV_KEGG.KEGG.File import cache
from FEV_KEGG.Evolution.Clade import Clade
from FEV_KEGG.Statistics import Percent
from FEV_KEGG.Robustness.Topology.Redundancy import RedundancyType, Redundancy, RedundancyContribution
from FEV_KEGG import settings
from FEV_KEGG.Util.Util import dictToHtmlFile

@cache(folder_path='experiments', file_name='deltaproteobacteria_clade')
def getCladeA():
    clade = Clade('Deltaproteobacteria')
    # pre-fetch collective metabolism into memory
    clade.collectiveMetabolism(excludeMultifunctionalEnzymes=settings.defaultNoMultifunctional)
    # pre-fetch collective enzyme metabolism into memory
    clade.collectiveMetabolismEnzymes(excludeMultifunctionalEnzymes=settings.defaultNoMultifunctional)
    return clade

if __name__ == '__main__':

    output = ['']

    #- get clade
    cladeA = getCladeA()
    majorityPercentageCoreMetabolism = 80
    majorityPercentageNeofunctionalisation = 0
    
    redundancyTypes = [RedundancyType.ROBUSTNESS_PARTIAL, RedundancyType.ROBUSTNESS_BOTH]
    
    output.append( 'core metabolism majority: ' + str(majorityPercentageCoreMetabolism) + '%' )
    output.append( 'neofunctionalisation majority: ' + str(majorityPercentageNeofunctionalisation) + '% (this means that gene duplication within a single organism is enough)' )
    output.append('')
    output.append(', '.join(cladeA.ncbiNames) + ':')
    output.append('')
    
    #- get core metabolism
    cladeAEcGraph = cladeA.coreMetabolism()
    cladeAEcCount = len(cladeAEcGraph.getECs())
    output.append( 'core metabolism ECs: ' + str(cladeAEcCount) )
    output.append('')
    
    #- calculate "neofunctionalised" ECs
    cladeANeofunctionalisedMetabolismSet = cladeA.neofunctionalisedECs(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation).getECs()
    cladeANeofunctionalisationsForFunctionChange = cladeA.neofunctionalisationsForFunctionChange(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation)
    
    #- calculate redundancy
    cladeARedundancy = Redundancy(cladeAEcGraph)
    cladeARedundancyContribution = RedundancyContribution(cladeARedundancy, cladeANeofunctionalisedMetabolismSet)
    
    for redundancyType in redundancyTypes:
        
        cladeARobustnessContributedECsForContributingNeofunctionalisedEC = cladeARedundancyContribution.getContributedKeysForSpecial(redundancyType)
        cladeARobustnessContributingNeofunctionalisedECs = set(cladeARobustnessContributedECsForContributingNeofunctionalisedEC.keys())
        
        #- REPEAT for each function change consisting of "neofunctionalised" ECs, which also contribute to redundancy
        output.append( '"neofunctionalised" ECs: ' + str(len(cladeANeofunctionalisedMetabolismSet)) + ' (' + str(Percent.getPercentStringShort(len(cladeANeofunctionalisedMetabolismSet), cladeAEcCount, 0)) + '%)' )
        
        robustnessContributingNeofunctionalisations = dict()
        
        for functionChange, neofunctionalisations in cladeANeofunctionalisationsForFunctionChange.items():
            #-     report enzyme pairs of neofunctionalisations, which caused the EC to be considered "neofunctionalised", and are in return contributing to redundancy        
            
            if functionChange.ecA in cladeARobustnessContributingNeofunctionalisedECs or functionChange.ecB in cladeARobustnessContributingNeofunctionalisedECs: # function change contributes to robustness
                
                for neofunctionalisation in neofunctionalisations:
                    currentSetOfContributedECs = robustnessContributingNeofunctionalisations.get(neofunctionalisation, None)
                    
                    if currentSetOfContributedECs is None:
                        currentSetOfContributedECs = set()
                        robustnessContributingNeofunctionalisations[neofunctionalisation] = currentSetOfContributedECs
                    
                    for ec in functionChange.ecPair:
                        contributedECs = cladeARobustnessContributedECsForContributingNeofunctionalisedEC.get(ec, None)
                        if contributedECs is not None:
                            currentSetOfContributedECs.update(contributedECs)
        
        output.append('')
        output.append( 'Neofunctionalisations contributing to robustness: ' + str(len(robustnessContributingNeofunctionalisations)) )
        output.append('')
        output.append('')
        
        
        neofunctionalisationsForContributedEC = dict()
        for neofunctionalisation, contributedECs in robustnessContributingNeofunctionalisations.items():
            
            for contributedEC in contributedECs:
                
                currentSetOfNeofunctionalisations = neofunctionalisationsForContributedEC.get(contributedEC, None)
                
                if currentSetOfNeofunctionalisations is None:
                    currentSetOfNeofunctionalisations = set()
                    neofunctionalisationsForContributedEC[contributedEC] = currentSetOfNeofunctionalisations
                    
                currentSetOfNeofunctionalisations.add(neofunctionalisation)
        
        
        ecNumbers = set()
        for contributedEC in neofunctionalisationsForContributedEC.keys():
            ecNumbers.add( contributedEC )
        
        dictToHtmlFile(neofunctionalisationsForContributedEC, cladeA.ncbiNames[0] + '_' + redundancyType.name + '_Neofunctionalisations-For-Contributed-EC.html', byValueFirst=False, inCacheFolder=True, addEcDescriptions = ecNumbers)
  
    
    
    
    
    
    for line in output:
        print( line )
