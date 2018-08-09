"""
Question
--------
Which EC numbers of the core metabolism of Archaea come into the core at which majority percentage?

Method
------
- get Archaea clade
- calculate collective metabolism, containing the majority percentage thershold for each edge
- export graph annotating each edge with its majority threshold (can be interpreted as transparency)

Result
------

::
    See archaea.png


Conclusion
----------
The resulting graph, interpreting the majority threshold at which an edge appears as transparency, shows that minorities of Archaea usually add to the core metabolism on the outside.
This adds capabilities to metabolise under more conditions. Sometimes ths even connects components of the graph not connected in a majority of Archaea.
However, it also seems quite widespread to add redundancy to the central parts of the core, increasing robustness.
Although sometimes, the added edges in the central parts of the core might also be necessary to prevent breakage of the graph, because another edge which is part of a majority of Archaea might be missing in a certain minority, 
requiring replacement without actually providing redundancy. This can not be visualised in this kind of graph.
"""
from FEV_KEGG.Evolution.Clade import Clade
from FEV_KEGG.Drawing import Export
from FEV_KEGG.KEGG.File import cache

@cache(folder_path='experiments/43', file_name='archaea_clade')
def getArchaeaClade():
    archaeaClade = Clade('/Archaea')
    # pre-fetch collective metabolism into memory
    archaeaClade.collectiveMetabolism(excludeMultifunctionalEnzymes=True)
    return archaeaClade

if __name__ == '__main__':
    
    #- get Archaea clade
    archaeaClade = getArchaeaClade()
    
    #- calculate collective metabolism, containing the majority percentage thershold for each edge
    collectiveMetabolism = archaeaClade.collectiveMetabolism(excludeMultifunctionalEnzymes=True)
    
    #- export graph annotating each edge with its majority threshold (can be interpreted as transparancy)
    Export.forCytoscape(collectiveMetabolism, file = 'experiments/43/archaea', inCacheFolder=True, totalNumberOfOrganisms = archaeaClade.organismsCount)