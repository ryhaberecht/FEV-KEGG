"""
Question
--------
How many core metabolism enzymes may have arisen due to gene duplication in Desulfobacterales, compared to Deltaproteobacteria?
How many of those performed neofunctionalisation?
How many new functions (EC numbers) arose in Desulfobacterales due to neofunctionalisation?

Method
------
- get NCBI taxonomy tree
- get group of organisms 'Proteobacteria/Deltaproteobacteria/Desulfobacterales'
- get supergroup of organisms 'Proteobacteria/Deltaproteobacteria'
- calculate number of enzymes in the group, including multifunctional enzymes
- calculate number of enzymes in the group, excluding multifunctional enzymes
- calculate number of enzymes in group's core metabolism
- calculate number of possible gene duplicates in the group
- calculate number of possible neofunctionalisations in the group
- calculate number of new EC numbers in the group
- calculate neofunctionalised EC numbers in the group

Result
------

::

    enzymes in group: 5245
    of which excluding multifunctional: 4661
    of which core metabolism: 2340
    of which gene duplicates: 770
    of which neofunctionalisations: 143
    EC numbers new in group: 166
    EC numbers due to neofunctionalisation: 23


Conclusion
----------
As we have seen before, the core metabolism of Desulfobacterales seems rather large. When applying a majority approach, instead of a consensus, this number is only going to rise.
Gene duplications, according to the simple model, are abundant. Many, but by far not all, gene duplications have also lead to new functions.  This can be well explained by patchwork evolution theory.
Only 9 new functions arose in Desulfobacterales core metabolism which can be explained by neofunctionalisation. Compared to 161 new functions which arose in total.
There may be several reasons for this:
1) gene duplication had been too obfuscated to be detected any longer, depends on E-value
2) horizontal gene transfer plays a dominating role in propagating new functions
3) majority core metabolism would yield much different results
""" 
from FEV_KEGG.Evolution import Events
from FEV_KEGG.Evolution.Taxonomy import NCBI
from FEV_KEGG.KEGG.Organism import Group


if __name__ == '__main__':
    
    output = []
    
    #- get NCBI taxonomy tree
    taxonomy = NCBI.getTaxonomy()
    
    #- get group of organisms 'Proteobacteria/Deltaproteobacteria/Desulfobacterales'
    group = Group( taxonomy.getOrganismAbbreviationsByPath('Bacteria/Proteobacteria/Deltaproteobacteria/Desulfobacterales') )
    
    #- get supergroup of organisms 'Proteobacteria/Deltaproteobacteria'
    supergroup = Group( taxonomy.getOrganismAbbreviationsByPath('Bacteria/Proteobacteria/Deltaproteobacteria') )
    
    #- calculate number of enzymes in the group, including multifunctional enzymes
    output.append( 'enzymes in group: ' + str( len( group.collectiveEnzymeGraph(noMultifunctional = False).getEnzymes() ) ) )
    
    #- calculate number of enzymes in the group, excluding multifunctional enzymes
    output.append( 'of which excluding multifunctional: ' + str( len( group.collectiveEnzymeGraph(noMultifunctional = True).getEnzymes() ) ) )
    
    #- calculate number of enzymes in group's core metabolism
    enzymeGraph = group.collectiveEnzymeGraphByEcConsensus(noMultifunctional = True)
    output.append( 'of which core metabolism: ' + str( len( enzymeGraph.getEnzymes() ) ) )
    
    #- calculate number of possible gene duplicates in the group
    geneDuplicationEnzymeGraph = Events.SimpleGeneDuplication.filterEnzymes(enzymeGraph)
    output.append( 'of which gene duplicates: ' + str( len( geneDuplicationEnzymeGraph.getEnzymes() ) ) )
    
    #- calculate number of possible neofunctionalisations in the group
    neofunctionalisationEnzymes = Events.NeofunctionalisedEnzymes(enzymeGraph.getEnzymes(), Events.SimpleGeneDuplication)
    output.append( 'of which neofunctionalisations: ' + str( len( neofunctionalisationEnzymes.getEnzymes() ) ) )
    
    #- calculate number of new EC numbers in the group
    newECs = group.consensusEcGraph(noMultifunctional = True).getECs().difference( supergroup.consensusEcGraph(noMultifunctional = True).getECs() )
    output.append('EC numbers new in group: ' + str( len( newECs ) ))
    
    #- calculate neofunctionalised EC numbers in the group    
    neoECs = Events.NeofunctionalisedECs(neofunctionalisationEnzymes).getECs()
    output.append( 'EC numbers due to neofunctionalisation: ' + str( len( neoECs ) ) )
    
    for line in output:
        print(line)
        