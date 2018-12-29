"""
Question
--------
How many EC numbers in Deltaproteobacteria are (partially) redundant?

Method
------
- get Deltaproteobacteria
- get core metabolism
- calculate redundancy
- print percentage of redundant ECs, for all types of redundancy

Result
------

::
    core metabolism majority: 80%

    Deltaproteobacteria:
    
    core metabolism ECs: 228
    
    
    Robustness fully: 7.0%
    Robustness partial: 4.8%
    Robustness both: 11.8%
    
    Flexibility fully: 17.1%
    Flexibility partial: 13.6%
    Flexibility both: 30.7%
    
    Target-flexibility fully: 21.9%
    Target-flexibility partial: 22.8%
    Target-flexibility both: 44.7%
    
    Source-flexibility fully: 40.8%
    Source-flexibility partial: 24.6%
    Source-flexibility both: 65.4%


Conclusion
----------
As expected, the occurence of redundancy increases with a more general definition.

However, it seems to be odd, at first, to have a much higher source-flexibility than target-flexibility.
This might be explained by the idea that, in general, metabolism has a main direction, because its purpose is to take one heap of molecules and turn it into a more useful heap of molecules.
This obviously only works if at least some enzymes only catalyse one direction of the reaction, the one towards the more useful molecule.
If we presume this is the case, there have to be substances which are mostly consumed, likely in the center of the metabolic network, and substances which are mostly produced, 
likely at the rim of the metabolic network. Substances mostly consumed are very likely to have redundant edges leading away from them, 
but any intermediate product might not be as likely to be the intermediate product of another path.
There is a similar case for end products: they are unlikely to be the end product of another path, but it might still be that their predecessor is the intermediate product of another path.
"""

from FEV_KEGG.KEGG.File import cache
from FEV_KEGG.Evolution.Clade import Clade
from FEV_KEGG.Statistics import Percent
from FEV_KEGG.Robustness.Topology.Redundancy import RedundancyType, Redundancy

@cache(folder_path='experiments', file_name='deltaproteobacteria_clade')
def getClade():
    clade = Clade('Deltaproteobacteria')
    # pre-fetch collective metabolism into memory
    clade.collectiveMetabolism(excludeMultifunctionalEnzymes=True)
    # pre-fetch collective enzyme metabolism into memory
    clade.collectiveMetabolismEnzymes(excludeMultifunctionalEnzymes=True)
    return clade

if __name__ == '__main__':

    output = ['']

    #- get Deltaproteobacteria
    clade = getClade()
    majorityPercentageCoreMetabolism = 80
    output.append( 'core metabolism majority: ' + str(majorityPercentageCoreMetabolism) + '%' )
    output.append('')
    output.append(', '.join(clade.ncbiNames) + ':')
    output.append('')
    
    #- get ECs
    cladeEcGraph = clade.coreMetabolism(majorityPercentageCoreMetabolism)
    cladeEcCount = len(cladeEcGraph.getECs())
    output.append( 'core metabolism ECs: ' + str(cladeEcCount) )
    output.append('')
    
    #- calculate redundancy
    cladeRedundancy = Redundancy(cladeEcGraph)
    
    #- print number of redundant ECs, including percentage of all ECs, for all types of redundancy
    output.append('')
    output.append( 'Robustness fully: ' + Percent.floatToPercentString(cladeRedundancy.getRedundancyRatio(RedundancyType.ROBUSTNESS) ) )
    output.append( 'Robustness partial: ' + Percent.floatToPercentString(cladeRedundancy.getRedundancyRatio(RedundancyType.ROBUSTNESS_PARTIAL) ) )
    output.append( 'Robustness both: ' + Percent.floatToPercentString(cladeRedundancy.getRedundancyRatio(RedundancyType.ROBUSTNESS_BOTH) ) )
    output.append('')
    output.append( 'Flexibility fully: ' + Percent.floatToPercentString(cladeRedundancy.getRedundancyRatio(RedundancyType.FLEXIBILITY) ) )
    output.append( 'Flexibility partial: ' + Percent.floatToPercentString(cladeRedundancy.getRedundancyRatio(RedundancyType.FLEXIBILITY_PARTIAL) ) )
    output.append( 'Flexibility both: ' + Percent.floatToPercentString(cladeRedundancy.getRedundancyRatio(RedundancyType.FLEXIBILITY_BOTH) ) )
    output.append('')
    output.append( 'Target-flexibility fully: ' + Percent.floatToPercentString(cladeRedundancy.getRedundancyRatio(RedundancyType.TARGET_FLEXIBILITY) ) )
    output.append( 'Target-flexibility partial: ' + Percent.floatToPercentString(cladeRedundancy.getRedundancyRatio(RedundancyType.TARGET_FLEXIBILITY_PARTIAL) ) )
    output.append( 'Target-flexibility both: ' + Percent.floatToPercentString(cladeRedundancy.getRedundancyRatio(RedundancyType.TARGET_FLEXIBILITY_BOTH) ) )
    output.append('')
    output.append( 'Source-flexibility fully: ' + Percent.floatToPercentString(cladeRedundancy.getRedundancyRatio(RedundancyType.SOURCE_FLEXIBILITY) ) )
    output.append( 'Source-flexibility partial: ' + Percent.floatToPercentString(cladeRedundancy.getRedundancyRatio(RedundancyType.SOURCE_FLEXIBILITY_PARTIAL) ) )
    output.append( 'Source-flexibility both: ' + Percent.floatToPercentString(cladeRedundancy.getRedundancyRatio(RedundancyType.SOURCE_FLEXIBILITY_BOTH) ) )
    output.append('')
    
    
    for line in output:
        print( line )
