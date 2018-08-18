"""
Context
-------
In :mod:`50`, we printed long lists of enzymes and EC numbers. This is hard to read and hard to use for further research, which is why we now export this information in HTML.

Question
--------
Which neofunctionalisations exist in Archaea?
Grouped by function change, sorted lexicographically, annotated with links to KEGG, annotated with human-readable names (if possible), and exported into an HTML file.

Method
------
- get clade
- get neofunctionalisations
- export to HTML, including links and further info

Result
------

::
    
    [see file Archaea.html]

Conclusion
----------
Much more useful than a plain text list of cryptic IDs.
"""

from FEV_KEGG.KEGG.File import cache
from FEV_KEGG.Evolution.Clade import Clade
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

    #- get clade
    clade = getClade()
    
    majorityPercentageCoreMetabolism = 80
    majorityPercentageNeofunctionalisation = 0
    
    #- get neofunctionalisations
    neofunctionalisationsForFunctionChange = clade.neofunctionalisationsForFunctionChange(majorityPercentageCoreMetabolism, majorityPercentageNeofunctionalisation)

    #- export to HTML, including links and further info
    dictToHtmlFile(neofunctionalisationsForFunctionChange, clade.ncbiNames[0], byValueFirst=False, inCacheFolder=True, addEcDescriptions = True)
    