"""
A unit test for the creation of collective enzyme graphs by using an EC consensus.

Two methods of doing so are compared and must always return the same value.
At last, the identity of included enzymes is compared, which must result in *0*.

Warnings
--------
The first result is time-dependant! It will certainly change on future updates of KEGG. However, it should be somewhere around *87*.
"""

import unittest

from FEV_KEGG.KEGG import Organism
import FEV_KEGG


class Test(unittest.TestCase):


    def test_consensusEcGraph_difference(self):
        
        FEV_KEGG.startProcessPool()
        
        enterobacteriales_organisms_abbreviations = ['eco', 'ses', 'sfl', 'ent', 'esa', 'kpn', 'cko', 'ype', 'spe', 'buc']
        enterobacteriales_organisms = Organism.Group(organismAbbreviations = enterobacteriales_organisms_abbreviations)
        
        enterobacteriales_organisms_abbreviations = ['eco', 'ses', 'sfl', 'ent', 'esa', 'kpn', 'cko', 'ype', 'spe', 'buc']
        gammaproteobacteria_organisms_abbreviations = ['hin', 'mht', 'xcc', 'vch', 'pae', 'acb', 'son', 'pha', 'amc', 'lpn', 'ftu', 'aha']
        gammaproteobacteria_organisms_abbreviations.extend(enterobacteriales_organisms_abbreviations) # extend with the sub-set, because they are also part of the set
        gammaproteobacteria_organisms = Organism.Group(organismAbbreviations = gammaproteobacteria_organisms_abbreviations)
        
         
         
          
        enterobacteriales_EC_graph = enterobacteriales_organisms.consensusEcGraph(noMultifunctional = True)
        gammaproteobacteria_EC_graph = gammaproteobacteria_organisms.consensusEcGraph(noMultifunctional = True)
        enterobacteriales_EC_set = enterobacteriales_EC_graph.getECs()
        gammaproteobacteria_EC_set = gammaproteobacteria_EC_graph.getECs()
        only_enterobacteriales_EC_set = enterobacteriales_EC_set.difference(gammaproteobacteria_EC_set)
            
        output = []
        for ec in only_enterobacteriales_EC_set:
            output.append(ec.__str__())
        
        result = len(output)
        print(str(result) + ' results')
        self.assertEqual(result, 87)
     
         
         
         
        enterobacteriales_enzyme_graph = enterobacteriales_organisms.collectiveEnzymeGraphByEcConsensus(noMultifunctional = True)
        gammaproteobacteria_enzyme_graph = gammaproteobacteria_organisms.collectiveEnzymeGraphByEcConsensus(noMultifunctional = True)
        enterobacteriales_enzyme_graph.removeMultifunctionalEnzymes()
        gammaproteobacteria_enzyme_graph.removeMultifunctionalEnzymes()
        enterobacteriales_enzymes = enterobacteriales_enzyme_graph.getEnzymes()
        gammaproteobacteria_enzymes = gammaproteobacteria_enzyme_graph.getEnzymes()
        
        enterobacteriales_EC_set_2 = set()
        for enzyme in enterobacteriales_enzymes:
            ecNumbers = enzyme.ecNumbers
            enterobacteriales_EC_set_2.update(ecNumbers)
              
        gammaproteobacteria_EC_set_2 = set()
        for enzyme in gammaproteobacteria_enzymes:
            ecNumbers = enzyme.ecNumbers
            gammaproteobacteria_EC_set_2.update(ecNumbers)
        
        only_enterobacteriales_EC_set_2 = enterobacteriales_EC_set_2.difference(gammaproteobacteria_EC_set_2)
          
        output = []
        for ec in only_enterobacteriales_EC_set_2:
            output.append(ec.__str__())
            
        result2 = len(output)
        print(str(result2) + ' results')
        self.assertEqual(result2, result)
    
         
         
         
        output = []
        for ec in only_enterobacteriales_EC_set_2.symmetric_difference(only_enterobacteriales_EC_set):
            output.append(ec.__str__())
        result3 = len(output)
        print(str(result3) + ' results')
        self.assertEqual(result3, 0)
        
        for ecString in output:
            print( ecString )


if __name__ == "__main__":
    
    unittest.main()