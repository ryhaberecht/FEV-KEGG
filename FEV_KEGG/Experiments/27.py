"""
Context
-------
Some organisms belong to an 'unclassified' taxon. They are hard to compare against other taxa.

Question
--------
How does excluding 'unclassified' Archaea change the number of EC numbers new to Thaumarchaeota?
How does it change the number of neofunctionalised EC numbers?

Method
------
- get NCBI taxonomy tree
- get group of organisms 'Archaea/Thaumarchaeota'
- get supergroup of organisms 'Archaea'
- calculate new EC numbers occuring in group's core metabolism compared to supergroup's core metabolism
- calculate neofunctionalised EC numbers in group's core metabolism
- repeat all, but excluding taxa containing 'unclassified' in both group and supergroup
- repeat all, but excluding taxa containing 'Nitrososphaeria' in both group and supergroup


Result
------

::

    taxon exception: None
    new EC numbers: 110
    neofunctionalised EC numbers: 11
    
    taxon exception: unclassified
    new EC numbers: 172
    neofunctionalised EC numbers: 16
    
    taxon exception: Nitrososphaeria
    new EC numbers: 114
    neofunctionalised EC numbers: 9


Conclusion
----------
Excluding 'unclassfied' organisms has a significant impact on core metabolism. This might be due to
1) reduced number of organisms, leading to higher chances of a consensus in metabolism
2) the unclassified organisms here have an unusual metabolism

Excluding a sub-group of the same size (4 organisms), however, does not yield significant changes. This implies a strong preference for option 2) above.
"""
from FEV_KEGG.Evolution.Events import SimpleGeneDuplication, NeofunctionalisedECs, NeofunctionalisedEnzymes
from FEV_KEGG.Evolution.Taxonomy import NCBI
import FEV_KEGG.KEGG.Organism as Organism


if __name__ == '__main__':
    
    output = []
    
    #- get NCBI taxonomy tree
    taxonomy = NCBI.getTaxonomy()
    
    #- repeat all, but excluding taxa containing 'unclassified' in both group and supergroup
    exceptPathsList = [None, 'unclassified', 'Nitrososphaeria']
    
    for exceptPaths in exceptPathsList:
    
        output.append( '\ntaxon exception: ' + str(exceptPaths) )
        
        #- get group of organisms 'Archaea/Thaumarchaeota'
        group = Organism.Group( taxonomy.getOrganismAbbreviationsByPath('Archaea/Thaumarchaeota', exceptPaths = exceptPaths, oneOrganismPerSpecies=False) )
        
        #- get supergroup of organisms 'Archaea'
        supergroup = Organism.Group( taxonomy.getOrganismAbbreviationsByPath('Archaea', exceptPaths = exceptPaths, oneOrganismPerSpecies=False) )
        
        #- calculate new EC numbers occuring in group's core metabolism compared to supergroup's core metabolism
        newECs = group.consensusEcGraph(noMultifunctional = True).getECs().difference( supergroup.consensusEcGraph(noMultifunctional = True).getECs() )
        output.append( 'new EC numbers: ' + str( len( newECs ) ) )
        
        #- calculate neofunctionalised EC numbers in group's core metabolism
        descendantEnzymeGraph = group.collectiveEnzymeGraphByEcConsensus(True)
        descendantOrganismsCount = group.organismsCount
        output.append( 'neofunctionalised EC numbers: ' + str( len( NeofunctionalisedECs(NeofunctionalisedEnzymes(descendantEnzymeGraph.getEnzymes(), SimpleGeneDuplication)).getECs() ) ) )
        
    for line in output:
        print(line)