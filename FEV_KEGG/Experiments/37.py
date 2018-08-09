"""
Context
-------
The concept of Simple Gene Duplication can be specialised further, to only include possibly duplicated genes which actually have an orthologous gene in the specified supergroup of ancestors.
This is what Chevron Gene Duplication is for.
This experiment is based on the method from :mod:`27`. Keep in mind, though, that this experiment contains updated data from KEGG and, thus, yields different results.

Question
--------
How does changing the definition of gene duplication from Simple to Chevron change the size of the resulting set of neofunctionalised ECs and Enzymes?

Method
------
- get NCBI taxonomy tree
- get group of organisms 'Archaea/Thaumarchaeota'
- get supergroup of organisms 'Archaea'
- calculate new EC numbers occuring in group's core metabolism compared to supergroup's core metabolism
- calculate neofunctionalised EC numbers in group's core metabolism
-     1) for Simple Gene Duplication
-     2) for Chevron Gene Duplication


Result
------

::

    new EC numbers: 110
    
    Consensus neofunctionalisation
    Simple gene duplication neofunctionalised EC numbers: 11
    Chevron gene duplication neofunctionalised EC numbers: 11
    
    Majority neofunctionalisation
    Simple gene duplication neofunctionalised EC numbers: 18
    Chevron gene duplication neofunctionalised EC numbers: 18


Conclusion
----------
The Chevron Gene Duplication model yields the same number of neofunctionalised EC numbers as the Simple model.
Even when using the majority neofunctionalisation approach, the numbers rise, but stay equal.

This leads to the question of whether there actually are abundant cases of genes with a paralog, but without any ortholog. 
Considering the general abundance of HGT, there are bound to be such cases, but hardly will they stand out when researching a small subset of Archaea.
"""
from FEV_KEGG.Evolution.Events import SimpleGeneDuplication, ChevronGeneDuplication, NeofunctionalisedECs,\
    NeofunctionalisedEnzymes
from FEV_KEGG.Evolution.Taxonomy import NCBI
import FEV_KEGG.KEGG.Organism as Organism


if __name__ == '__main__':
    
    output = []
    
    #- get NCBI taxonomy tree
    taxonomy = NCBI.getTaxonomy()
    
    #- get group of organisms 'Archaea/Thaumarchaeota'
    group = Organism.Group( taxonomy.getOrganismAbbreviationsByPath('Archaea/Thaumarchaeota') )
    
    #- get supergroup of organisms 'Archaea'
    supergroup = Organism.Group( taxonomy.getOrganismAbbreviationsByPath('Archaea') )
    
    #- calculate new EC numbers occuring in group's core metabolism compared to supergroup's core metabolism
    newECs = group.consensusEcGraph(noMultifunctional = True).getECs().difference( supergroup.consensusEcGraph(noMultifunctional = True).getECs() )
    output.append( 'new EC numbers: ' + str( len( newECs ) ) )
    
    #- calculate neofunctionalised EC numbers in group's core metabolism
    descendantEnzymeGraph = group.collectiveEnzymeGraphByEcConsensus(noMultifunctional = True)
    output.append( '\nConsensus neofunctionalisation' )
    output.append( 'Simple gene duplication neofunctionalised EC numbers: ' + str( len( NeofunctionalisedECs(NeofunctionalisedEnzymes(descendantEnzymeGraph.getEnzymes(), SimpleGeneDuplication)).getECs() ) ) )
    output.append( 'Chevron gene duplication neofunctionalised EC numbers: ' + str( len( NeofunctionalisedECs(NeofunctionalisedEnzymes(descendantEnzymeGraph.getEnzymes(), ChevronGeneDuplication(group))).getECs() ) ) )
    descendantEnzymeGraph = group.collectiveEnzymeGraphByEcMajority(majorityPercentage = 90, noMultifunctional = True)
    output.append( '\nMajority neofunctionalisation' )
    output.append( 'Simple gene duplication neofunctionalised EC numbers: ' + str( len( NeofunctionalisedECs(NeofunctionalisedEnzymes(descendantEnzymeGraph.getEnzymes(), SimpleGeneDuplication)).getECs() ) ) )
    output.append( 'Chevron gene duplication neofunctionalised EC numbers: ' + str( len( NeofunctionalisedECs(NeofunctionalisedEnzymes(descendantEnzymeGraph.getEnzymes(), ChevronGeneDuplication(group))).getECs() ) ) )
    
    for line in output:
        print(line)