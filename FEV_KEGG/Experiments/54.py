"""
Context
-------
Combining the list of all neofunctionalisations, sorted by their function change, from :mod:`51`, 
and the list of neofunctionalisations contributing to robustness, sorted by the EC number they provide robustness for, from :mod:`53`.

E-value also plays an important role in the number of neofunctionalisations reported, 
because the bigger the E-value is, the easier two enzymes match as gene-duplicates, the higher the probability of a discovered function change, together resulting in a neofunctionalisation.

Question
--------
Which neofunctionalisations exist in the core metabolism of Archaea, using a certain E-value?
Grouped by function change, sorted lexicographically, annotated with links to KEGG, annotated with human-readable names (if possible), and exported into an HTML file.

Which neofunctionalised enzymes cause the core metabolism of Archaea to have increased redundancy?
Grouped by contributed functions, sorted lexicographically, annotated with links to KEGG, annotated with human-readable names (if possible), and exported into an HTML file.

Method
------
- get clade
- get core metabolism
- calculate "neofunctionalised" ECs
- report neofunctionalisations
-     print them into nice HTML
- calculate redundancy
- REPEAT for each "neofunctionalised" EC contributing to redundancy
-     report enzyme pairs of neofunctionalisations, which caused the EC to be considered "neofunctionalised", and are in return contributing to redundancy
-         print them into nice HTML

Result
------

::

    core metabolism majority: 80%
    neofunctionalisation majority: 0% (this means that gene duplication within a single organism is enough)
    
    Archaea:
    
    core metabolism ECs: 114
    
    
    All neofunctionalisations: 1140
        [see Archaea_Neofunctionalisations-For-FunctionChange.html]
            
    "neofunctionalised" ECs: 16 (14%)
    
    Neofunctionalisations contributing to robustness: 93 (8%)
        [see Archaea_Neofunctionalisations-For-Contributed-EC.html]
    

Conclusion
----------
Archaea do not seem to have many neofunctionalisations contributing to robustness.
This might be limited by:
- the size of the core metabolism, i.e. the majority percentage.
- the E-value. The lower the E-value, the closer two enzymes' sequences have to match, for them to be recognised as gene-duplicated.
Using a higher E-value will greatly increase the number of neofunctionalisations reported. Of course, this will also increase the number of false-positives.
"""
from FEV_KEGG.KEGG.File import cache
from FEV_KEGG.Evolution.Clade import Clade
from FEV_KEGG.Statistics import Percent
from FEV_KEGG.Robustness.Topology.Redundancy import RedundancyType, Redundancy, RedundancyContribution
from FEV_KEGG import settings
from FEV_KEGG.Util.Util import dictToHtmlFile

@cache(folder_path='experiments', file_name='archaea_clade')
def getClade():
    clade = Clade('Archaea')
    # pre-fetch collective metabolism into memory
    clade.collectiveMetabolism(excludeMultifunctionalEnzymes=settings.defaultNoMultifunctional)
    # pre-fetch collective enzyme metabolism into memory
    clade.collectiveMetabolismEnzymes(excludeMultifunctionalEnzymes=settings.defaultNoMultifunctional)
    return clade

if __name__ == '__main__':

    output = ['']

    #- get clade
    clade = getClade()
    majorityPercentageCoreMetabolism = 80
    majorityPercentageNeofunctionalisation = 0
    eValue = 1e-15
    redundancyType = RedundancyType.ROBUSTNESS_BOTH
    
    output.append( 'core metabolism majority: ' + str(majorityPercentageCoreMetabolism) + '%' )
    output.append( 'neofunctionalisation majority: ' + str(majorityPercentageNeofunctionalisation) + '% (this means that gene duplication within a single organism is enough)' )
    output.append('')
    output.append(', '.join(clade.ncbiNames) + ':')
    output.append('')
    
    #- get core metabolism
    cladeEcGraph = clade.coreMetabolism(majorityPercentageCoreMetabolism)
    cladeEcCount = len(cladeEcGraph.getECs())
    output.append( 'core metabolism ECs: ' + str(cladeEcCount) )
    output.append('')
    
    #- report neofunctionalisations
    neofunctionalisationsForFunctionChange = clade.neofunctionalisationsForFunctionChange(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation, eValue=eValue)
    allNeofunctionalisations = set() # set of all neofunctionalisations, no matter which function change they belong to
    for valueSet in neofunctionalisationsForFunctionChange.values():
        allNeofunctionalisations.update( valueSet )

    output.append('')
    output.append( 'All neofunctionalisations: ' + str(len(allNeofunctionalisations)) )
    
    #-     print them into nice HTML
    ecNumbers = set()
    for functionChange in neofunctionalisationsForFunctionChange.keys():
        ecNumbers.update( functionChange.ecPair )
    dictToHtmlFile(neofunctionalisationsForFunctionChange, clade.ncbiNames[0] + '_Neofunctionalisations-For-FunctionChange.html', byValueFirst=False, inCacheFolder=True, addEcDescriptions=ecNumbers)
    output.append( '\t[see ' + clade.ncbiNames[0] + '_Neofunctionalisations-For-FunctionChange.html]' )
    output.append('')
    
    #- calculate "neofunctionalised" ECs
    cladeNeofunctionalisedMetabolismSet = clade.neofunctionalisedECs(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation, eValue=eValue).getECs()
    cladeNeofunctionalisationsForFunctionChange = clade.neofunctionalisationsForFunctionChange(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation, eValue=eValue)
    
    #- calculate redundancy
    cladeRedundancy = Redundancy(cladeEcGraph)
    cladeRedundancyContribution = RedundancyContribution(cladeRedundancy, cladeNeofunctionalisedMetabolismSet)
        
    cladeRobustnessContributedECsForContributingNeofunctionalisedEC = cladeRedundancyContribution.getContributedKeysForSpecial(redundancyType)
    cladeRobustnessContributingNeofunctionalisedECs = set(cladeRobustnessContributedECsForContributingNeofunctionalisedEC.keys())
    
    #- REPEAT for each function change consisting of "neofunctionalised" ECs, which also contribute to redundancy
    output.append( '"neofunctionalised" ECs: ' + str(len(cladeNeofunctionalisedMetabolismSet)) + ' (' + str(Percent.getPercentStringShort(len(cladeNeofunctionalisedMetabolismSet), cladeEcCount, 0)) + '%)' )
    
    robustnessContributingNeofunctionalisations = dict()
    
    for functionChange, neofunctionalisations in cladeNeofunctionalisationsForFunctionChange.items():
        #-     report enzyme pairs of neofunctionalisations, which caused the EC to be considered "neofunctionalised", and are in return contributing to redundancy        
        
        if functionChange.ecA in cladeRobustnessContributingNeofunctionalisedECs or functionChange.ecB in cladeRobustnessContributingNeofunctionalisedECs: # function change contributes to robustness
            
            for neofunctionalisation in neofunctionalisations:
                currentSetOfContributedECs = robustnessContributingNeofunctionalisations.get(neofunctionalisation, None)
                
                if currentSetOfContributedECs is None:
                    currentSetOfContributedECs = set()
                    robustnessContributingNeofunctionalisations[neofunctionalisation] = currentSetOfContributedECs
                
                for ec in functionChange.ecPair:
                    contributedECs = cladeRobustnessContributedECsForContributingNeofunctionalisedEC.get(ec, None)
                    if contributedECs is not None:
                        currentSetOfContributedECs.update(contributedECs)
    
    output.append('')
    output.append( 'Neofunctionalisations contributing to robustness: ' + str(len(robustnessContributingNeofunctionalisations)) + ' (' + str(Percent.getPercentStringShort(len(robustnessContributingNeofunctionalisations), len(allNeofunctionalisations), 0)) + '%)' )
    output.append( '\t[see ' + clade.ncbiNames[0] + '_Neofunctionalisations-For-Contributed-EC.html]' )
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
    
    #-         print them into nice HTML
    ecNumbers = set()
    for contributedEC in neofunctionalisationsForContributedEC.keys():
        ecNumbers.add( contributedEC )
    
    dictToHtmlFile(neofunctionalisationsForContributedEC, clade.ncbiNames[0] + '_Neofunctionalisations-For-Contributed-EC.html', byValueFirst=False, inCacheFolder=True, addEcDescriptions=ecNumbers)
    
    
    
    for line in output:
        print( line )
    