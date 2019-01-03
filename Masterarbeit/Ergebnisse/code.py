from FEV_KEGG.KEGG.File import cache
from FEV_KEGG.Evolution.Clade import Clade
from FEV_KEGG.Statistics import Percent
from FEV_KEGG.Robustness.Topology.Redundancy import RedundancyType, Redundancy, RedundancyContribution
from FEV_KEGG.Util.Util import dictToHtmlFile
from FEV_KEGG.Drawing import Export

@cache(folder_path='experiments', file_name='alphaproteobacteria_clade')
def getClade():
    clade = Clade('Alphaproteobacteria')
    # pre-fetch collective metabolism into memory
    clade.collectiveMetabolism(excludeMultifunctionalEnzymes=True)
    # pre-fetch collective enzyme metabolism into memory
    clade.collectiveMetabolismEnzymes(excludeMultifunctionalEnzymes=True)
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
    output.append(', '.join(clade.ncbiNames) + ': ' + str(clade.organismsCount) + ' organisms')
    output.append('')
    
    #- get core metabolism
    cladeEcGraph = clade.coreMetabolism(majorityPercentageCoreMetabolism)
    cladeEcCount = len(cladeEcGraph.getECs())
    output.append( 'core metabolism ECs: ' + str(cladeEcCount) )
    output.append('')

    #- calculate "neofunctionalised" ECs
    cladeNeofunctionalisedMetabolismSet = clade.neofunctionalisedECs(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation, eValue=eValue).getECs()
    cladeNeofunctionalisationsForFunctionChange = clade.neofunctionalisationsForFunctionChange(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation, eValue=eValue)
    
    #- calculate redundancy
    cladeRedundancy = Redundancy(cladeEcGraph)
    cladeRedundancyContribution = RedundancyContribution(cladeRedundancy, cladeNeofunctionalisedMetabolismSet)
        
    cladeRobustnessContributedECsForContributingNeofunctionalisedEC = cladeRedundancyContribution.getContributedKeysForSpecial(redundancyType)
    cladeRobustnessContributingNeofunctionalisedECs = set(cladeRobustnessContributedECsForContributingNeofunctionalisedEC.keys())
    
    #- export graph of core metabolism, colouring the edges of (contributing) neofunctionalised ECs
    edgesOfNeofunctionalisedECs = set()
    for key in cladeNeofunctionalisedMetabolismSet:
        for edge in cladeEcGraph.getEdgesFromKey(key):
            edgesOfNeofunctionalisedECs.add(edge)
    Export.addColourAttribute(cladeEcGraph, Export.Colour.BLUE, nodes = False, edges = edgesOfNeofunctionalisedECs)
    
    edgesOfContributingNeofunctionalisedECs = set()
    for key in cladeRobustnessContributingNeofunctionalisedECs:
        for edge in cladeEcGraph.getEdgesFromKey(key):
            edgesOfContributingNeofunctionalisedECs.add(edge)
    Export.addColourAttribute(cladeEcGraph, Export.Colour.GREEN, nodes = False, edges = edgesOfContributingNeofunctionalisedECs)
    
    Export.forCytoscape(cladeEcGraph, clade.ncbiNames[0], inCacheFolder = True, addDescriptions = True, totalNumberOfOrganisms = clade.organismsCount)
    
    #- REPEAT for each function change consisting of "neofunctionalised" ECs, which also contribute to redundancy
    output.append( '"neofunctionalised" ECs: ' + str(len(cladeNeofunctionalisedMetabolismSet)) + ' (' + str(Percent.getPercentStringShort(len(cladeNeofunctionalisedMetabolismSet), cladeEcCount, 0)) + '%)' )
    
    robustKeys = cladeRedundancy.getRedundantKeys(redundancyType)
    robustKeysCount = len(robustKeys)
    output.append( 'robust ECs: ' + str(robustKeysCount) + ' (' + str(Percent.getPercentStringShort(robustKeysCount, cladeEcCount, 0)) + '%)' )
    
    output.append( '    of which are robust due to "neofunctionalised" ECs: ' + str(Percent.getPercentStringShort(cladeRedundancyContribution.getKeyContributionRatio(redundancyType), 1, 0)) + '%' )
    
    
    
    #- report neofunctionalisations
    neofunctionalisationsForFunctionChange = clade.neofunctionalisationsForFunctionChange(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation, eValue=eValue)
    allNeofunctionalisations = set() # set of all neofunctionalisations, no matter which function change they belong to
    for valueSet in neofunctionalisationsForFunctionChange.values():
        allNeofunctionalisations.update( valueSet )

    output.append('')
    output.append('')
    output.append( 'All neofunctionalisations: ' + str(len(allNeofunctionalisations)) )
    
    #-     print them into nice HTML
    ecNumbers = set()
    for functionChange in neofunctionalisationsForFunctionChange.keys():
        ecNumbers.update( functionChange.ecPair )
    dictToHtmlFile(neofunctionalisationsForFunctionChange, clade.ncbiNames[0] + '_Neofunctionalisations-For-FunctionChange.html', byValueFirst=False, inCacheFolder=True, addEcDescriptions=ecNumbers)
    output.append( '\t[see ' + clade.ncbiNames[0] + '_Neofunctionalisations-For-FunctionChange.html]' )
    output.append('')
    
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
    