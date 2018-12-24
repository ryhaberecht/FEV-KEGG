"""
Context
-------
Closer look at the possible reasons for occurence of EC numbers completely unknown to us in experiment :mod:`32`.
1. EC number is associated with the organism, but not listed in one of KEGG's hand-drawn pathways. For example 1.13.11.24 is associated with all our 15 organisms, but not present in any pathway.
2. As seen in experiment :mod:`30`, there may be EC numbers predicted by Oh et al. which are outdated.
3. Oh et al. used a compilation of several sources, some may have predicted EC numbers for B. subtilis which never made their way into KEGG at all, which is our only source.

Question
--------
Does the high number of EC numbers only in their set from experiment :mod:`32` result from outdated/faulty data?

Method
------
- get all EC numbers known to any organism in KEGG, using NUKA.
- take the EC numbers only in Oh's set at 1% majority (see :mod:`32`).
- keep only the ones occuring in NUKA.
- keep only the ones with a gene occuring in any of our 15 KEGG organisms found as "Bacillus subtilis"


Result
------

::

    ECs not in any pathway: 50
    Leaving ECs only in theirs: 58
    ---------------------------------
    ECs not in any of our organisms: 57
    Leaving ECs only in theirs: 1
    3.1.3.3


Conclusion
----------
1. EC number is associated with the organism, but not listed in one of KEGG's hand-drawn pathways. For example 1.13.11.24 is associated with all our 15 organisms, but not present in any pathway.

About half (50) of the EC numbers only in Oh's set are not listed in any of KEGG's pathways today. While B. subtilis may contain them, KEGG does not, and thus our model can not.


2. As seen in experiment :mod:`30`, there may be EC numbers predicted by Oh et al. which are outdated.
3. Oh et al. used a compilation of several sources, some may have predicted EC numbers for B. subtilis which never made their way into KEGG at all, which is our only source.

After removing all EC numbers unknown to our organisms in KEGG today, only one remains: 3.1.3.3. This seems to be another case of inconsistent data in KEGG, because 3.1.3.3 supposedly is part of pathway 00260 and 00680, 
and even has three associated genes in 'bsu', but it is nowhere to be found in the actual pathways bsu00260 or bsu00680.

Answering the question: yes, all EC numbers missing in our model can be explained by outdated or faulty data in either their set or in KEGG. Discriminating between these two possibilities is impossible for us.

While this experiment does not imply that our model is complete, it does imply that it is correct under the constraints imposed by the completeness and correctness of the KEGG database.
"""
from FEV_KEGG.Graph.Elements import EcNumber
from FEV_KEGG.Evolution.Taxonomy import NCBI
from FEV_KEGG.KEGG import Organism
from FEV_KEGG.KEGG.NUKA import NUKA
import FEV_KEGG.lib.Biopython.KEGG.REST

def isKnownToKEGG(ecNumberString, organisms):
    
    download = FEV_KEGG.lib.Biopython.KEGG.REST.kegg_get(ecNumberString).read()

    if download.find("PATHWAY     ") <= -1:
        return False
    
    foundAnOrganism = False
    for abbreviation in organisms:
        uppercase = abbreviation.upper() + ": "
        if download.find(uppercase) >= 0:
            foundAnOrganism = True
            break
    
    if foundAnOrganism == False:
        return False
    else:
        return True


if __name__ == '__main__':
    
    output = []
    
    #- get all EC numbers known to any organism in KEGG, using NUKA.
    nuka = NUKA()
    nukaECs = nuka.substanceEcGraph.getECs()
    
    #- take the EC numbers only in Oh's set at 1% majority (see :mod:`32`)
    theirECnumberStrings = ['5.3.1.23', '1.4.1.13', '6.3.5.5', '6.3.1.2', '2.3.1.1', '3.5.1.47', '1.4.1.1', '1.5.1.2', '2.1.3.3', '2.7.2.8', '2.6.1.11', '3.5.3.11', '1.2.1.38', '2.6.1.2', '4.1.1.19', '1.2.1.11', '6.3.5.4', '2.7.2.4', '2.6.1.1', '4.1.1.20', '5.1.1.7', '1.3.1.26', '4.2.1.52', '1.1.1.3', '2.7.1.39', '4.2.3.1', '1.2.1.41', '2.7.2.11', '4.1.3.18', '1.2.1.27', '1.4.1.2', '4.1.3.27', '2.4.2.18', '5.4.99.5', '4.2.3.5', '4.1.2.15', '4.2.1.10', '4.2.3.4', '4.1.1.48', '1.3.1.12', '4.2.1.51', '5.3.1.24', '2.5.1.19', '1.1.1.25', '2.7.1.71', '4.2.1.20', '2.6.1.5', '4.2.1.22', '4.1.99.1', '2.1.2.1', '2.7.1.25', '4.2.99.8', '4.4.1.1', '4.4.1.8', '1.1.1.95', '2.6.1.52', '2.3.1.30', '2.6.1.9', '2.6.1.58', '1.4.3.3', '2.5.1.6', '2.6.1.16', '2.1.1.13', '2.3.2.2', '1.5.3.1', '3.2.2.16', '4.1.1.50', '2.7.1.100', '2.5.1.16', '2.5.1.22', '2.6.1.19', '3.5.1.16', '4.2.99.8', '4.2.99.8', '5.1.1.1', '3.5.3.1', '4.3.2.1', '6.3.4.5', '2.4.2.17', '4.2.1.9', '4.2.1.9', '3.5.3.8', '3.5.1.31', '3.5.1.2', '5.1.1.3', '2.3.1.29', '1.1.1.31', '4.3.1.3', '1.1.1.23', '3.1.3.15', '1.5.1.12', '2.3.1.46', '4.2.1.19', '1.4.1.9', '2.6.1.42', '1.1.1.85', '4.2.1.33', '4.1.3.12', '3.5.2.7', '1.1.1.86', '1.1.1.86', '2.6.1.42', '1.4.1.9', '5.4.3.2', '4.1.1.18', '1.2.1.27', '1.2.1.27', '2.6.1.13', '2.3.1.35', '1.5.1.12', '3.5.4.19', '3.6.1.31', '5.3.1.16', '1.5.1.2', '3.1.3.3', '3.5.1.18', '4.3.1.18', '4.3.1.17', '4.2.99.9', '4.2.99.9', '4.2.99.9', '1.2.1.16', '1.1.1.103', '4.2.1.49', '3.5.1.5', '1.4.1.9', '2.6.1.42', '3.5.1.1', '4.3.1.1', '2.3.1.57', '2.3.1.57', '1.2.1.27', '2.1.1.37', '3.1.2.4', '3.1.3.3', '4.1.1.29', '4.2.99.9', '1.13.11.24', '1.1.1.2', '1.1.1.4', '4.1.1.5', '1.1.1.94', '2.7.1.12', '2.7.1.15', '1.1.1.14', '6.2.1.1', '6.4.1.1', '1.2.1.12', '1.1.1.44', '5.4.2.2', '5.1.3.1', '5.3.1.6', '1.1.99.5', '1.1.1.49', '2.7.2.1', '1.2.1.12', '4.2.1.11', '4.1.2.13', '1.1.1.17', '4.1.1.49', '1.1.1.40', '1.1.1.21', '3.2.1.108', '5.1.3.14', '4.2.1.42', '4.2.1.40', '4.1.1.2', '1.1.3.15', '5.4.2.10', '3.5.1.25', '1.1.1.267', '3.5.99.6', '3.2.1.55', '4.1.2.13', '4.2.1.7', '5.3.1.12', '1.1.1.57', '4.2.1.8', '1.1.1.58', '5.3.1.12', '2.7.1.2', '5.3.1.17', '1.1.1.127', '2.7.1.45', '4.1.2.14', '3.2.1.20', '3.2.1.48', '2.7.1.15', '1.1.1.21', '4.1.1.47', '1.1.1.60', '1.1.1.83', '4.1.1.73', '1.1.1.93', '1.1.1.93', '1.2.1.3', '5.1.3.2', '1.3.99.1', '6.2.1.5', '4.1.3.31', '4.1.3.30', '4.2.1.3', '4.1.3.7', '4.2.1.2', '1.1.1.42', '1.1.1.37', '3.2.1.21', '3.2.1.21', '3.2.1.26', '3.2.1.1', '3.2.1.86', '1.2.1.19', '3.6.1.13', '1.2.1.46', '1.2.1.3', '1.2.1.3', '1.2.1.3', '5.3.1.4', '3.2.1.21', '2.7.2.7', '4.2.1.41', '2.7.1.56', '3.1.3.11', '2.7.1.4', '1.1.1.118', '2.7.1.41', '5.3.1.9', '2.7.1.6', '3.2.1.22', '3.2.1.22', '3.2.1.22', '3.2.1.22', '2.7.7.10', '2.4.1.18', '1.2.1.21', '2.4.1.1', '2.4.1.21', '2.7.7.27', '2.7.1.30', '4.1.3.5', '3.2.1.26', '3.2.1.122', '2.3.1.79', '5.3.1.8', '4.2.3.3', '2.8.3.5', '2.3.1.19', '2.7.1.11', '3.2.1.86', '3.1.3.18', '5.4.2.6', '3.1.3.8', '5.4.2.8', '6.4.1.3', '2.7.9.2', '2.7.1.47', '2.7.1.16', '5.1.3.4', '2.7.1.5', '3.2.1.86', '2.2.1.2', '2.2.1.1', '2.2.1.1', '3.2.1.93', '5.1.3.14', '5.3.1.5', '5.3.1.5', '2.7.1.17', '2.3.1.8', '1.2.1.51', '2.7.1.11', '5.3.1.9', '5.4.2.1', '2.7.1.40', '5.3.1.1', '2.7.2.3', '5.3.1.14', '1.1.1.18', '4.1.3.19', '6.2.1.1', '5.3.1.3', '4.3.1.7', '3.5.1.49', '4.1.2.20', '2.7.1.31', '5.4.99.11', '1.2.1.22', '3.1.3.25', '2.4.1.8', '3.1.1.31', '4.1.2.19', '1.1.99.21', '4.1.2.40', '5.1.3.7', '3.5.3.19', '3.2.1.37', '2.7.8.13', '2.7.7.9', '2.7.7.39', '1.1.1.22', '2.7.7.23', '6.3.2.4', '3.5.1.28', '3.2.1.52', '6.3.2.13', '2.5.1.7', '6.3.2.9', '6.3.2.8', '1.1.1.158', '6.3.2.15', '2.3.1.157', '3.6.1.27', '1.5.1.5', '1.7.99.5', '2.7.7.1', '2.7.7.3', '2.7.1.24', '1.1.1.169', '2.5.1.30', '5.4.99.6', '2.1.2.11', '2.4.2.11', '3.5.1.19', '2.7.7.18', '4.1.3.36', '6.3.2.1', '2.7.1.34', '6.2.1.26', '2.7.7.58', '2.3.1.47', '1.1.1.193', '1.5.1.3', '6.3.2.12', '2.5.1.15', '2.5.1.9', '3.5.4.16', '3.5.4.9', '4.1.2.25', '2.6.1.62', '2.8.1.6', '6.2.1.14', '6.3.3.3', '3.5.4.26', '2.7.1.50', '2.7.6.3', '2.7.4.7', '2.5.1.9', '2.7.4.16', '2.5.1.3', '4.1.1.11', '1.4.3.16', '1.3.1.28', '1.5.1.3', '2.7.7.2', '3.5.1.10', '3.5.4.25', '3.3.2.1', '4.1.1.71', '2.7.1.33', '1.3.3.4', '4.2.3.12', '2.7.1.26', '2.4.2.19', '6.3.1.5', '3.5.99.2', '1.3.3.3', '4.99.1.1', '5.4.3.8', '4.3.1.8', '2.1.1.107', '4.2.1.24', '4.2.1.75', '4.1.1.37', '2.7.1.49', '3.1.3.5', '4.1.1.36', '6.3.2.5', '2.3.1.15', '3.1.3.27', '2.7.7.41', '4.1.1.65', '6.4.1.2', '1.2.1.25', '1.2.1.25', '1.2.1.25', '2.5.1.1', '5.3.3.2', '2.3.1.41', '2.3.1.41', '2.3.1.41', '2.3.1.41', '2.3.1.41', '2.3.1.41', '2.3.1.41', '2.3.1.41', '2.3.1.41', '2.3.1.41', '2.3.1.41', '2.3.1.41', '2.3.1.41', '2.3.1.41', '2.7.8.7', '6.2.1.3', '6.2.1.3', '6.2.1.3', '6.2.1.3', '6.2.1.3', '6.2.1.3', '6.2.1.3', '6.2.1.3', '6.2.1.3', '6.2.1.3', '6.2.1.3', '6.2.1.3', '6.2.1.3', '2.5.1.10', '2.7.8.8', '2.3.1.51', '2.3.1.16', '2.3.1.9', '2.3.1.16', '2.3.1.16', '2.3.1.16', '2.3.1.16', '2.3.1.16', '2.3.1.16', '1.3.99.2', '1.3.99.2', '1.3.99.3', '1.3.99.2', '1.3.99.3', '1.3.99.3', '1.3.99.3', '1.3.99.3', '1.3.99.3', '3.1.4.14', '2.7.1.107', '4.2.1.17', '4.2.1.17', '4.2.1.17', '4.2.1.17', '4.2.1.17', '4.2.1.17', '4.2.1.17', '4.2.1.17', '4.2.1.17', '6.2.1.3', '3.1.4.46', '3.1.4.46', '3.1.4.46', '1.1.1.35', '1.1.1.35', '1.1.1.35', '1.1.1.35', '1.1.1.35', '1.1.1.35', '1.1.1.35', '1.1.1.35', '1.1.1.35', '4.1.3.4', '2.5.1.31', '2.5.1.29', '2.5.1.33', '1.1.1.27', '1.6.5.3', '1.7.99.4', '3.6.3.14', '1.2.1.2', '1.5.1.30', '1.8.1.9', '1.9.3.1', '1.10.2.2', '1.5.1.29', '1.1.1.27', '3.1.3.5', '3.1.3.5', '3.1.3.5', '3.1.3.5', '3.1.3.5', '3.1.3.5', '3.1.3.5', '3.1.3.5', '3.1.3.5', '3.1.3.5', '3.1.4.16', '3.1.4.16', '3.1.4.16', '3.1.4.16', '3.1.3.6', '3.1.3.6', '3.1.3.6', '3.1.3.6', '2.7.1.23', '2.7.1.23', '2.7.1.23', '2.7.1.23', '2.7.1.23', '2.7.1.23', '3.2.2.9', '3.5.4.12', '3.5.4.12', '1.17.4.1', '1.17.4.1', '1.17.4.1', '1.17.4.1', '2.7.6.1', '3.5.4.3', '3.5.4.5', '2.4.2.2', '5.4.2.7', '5.4.2.7', '2.4.2.4', '2.4.2.2', '3.5.4.14', '2.7.1.74', '3.5.4.2', '2.7.1.20', '3.5.2.5', '2.7.4.14', '2.7.4.14', '2.7.4.11', '2.7.1.76', '2.7.1.113', '4.1.2.4', '2.7.4.9', '3.6.1.23', '2.7.4.8', '1.7.1.7', '2.7.6.5', '2.4.2.8', '2.7.4.6', '2.7.4.6', '2.7.4.6', '2.7.4.6', '2.7.4.6', '2.7.4.6', '2.7.4.6', '6.3.2.6', '2.4.2.1', '2.4.2.1', '2.4.2.1', '2.4.2.1', '2.4.2.1', '2.4.2.1', '2.4.2.1', '2.4.2.1', '2.4.2.2', '2.4.2.1', '2.7.1.21', '2.1.1.45', '1.7.3.3', '2.7.1.48', '2.7.1.48', '2.7.1.48', '1.1.1.204', '2.4.2.22', '6.3.4.13', '6.3.3.1', '6.3.5.3', '6.3.4.2', '2.7.4.3', '2.4.2.7', '4.3.2.2', '4.3.2.2', '6.3.4.4', '2.1.2.3', '2.1.3.2', '6.3.4.2', '1.3.3.1', '3.5.2.3', '2.1.2.2', '2.7.4.8', '2.4.2.14', '6.3.5.2', '3.5.4.10', '1.1.1.205', '2.7.4.6', '2.4.2.10', '2.4.2.9', '2.7.4.3', '3.5.3.4', '3.5.4.1', '2.7.4.6', '4.1.1.23', '2.7.4.4', '1.1.1.204', '1.2.1.8', '1.1.99.1', '3.5.2.6', '2.1.2.9', '1.11.1.9', '6.1.1.10', '3.5.1.11', '1.15.1.1', '6.1.1.17', '1.11.1.6', '3.6.1.15', '4.2.1.1', '2.7.7.4', '1.8.1.2', '3.1.3.1', '3.6.1.1', '3.1.3.1', '3.6.1.2']
    theirECnumbers = set()
    
    for string in theirECnumberStrings:
        theirECnumbers.add( EcNumber(string) )
    
    taxonomy = NCBI.getTaxonomy()
    organisms = taxonomy.getOrganismAbbreviationsByPath('Bacillus subtilis', oneOrganismPerSpecies=False)
    group = Organism.Group( organisms )
    ourECnumbers = group.majorityEcGraph(majorityPercentage = 1, noMultifunctional = False).getECs()
    ourECnumbers = EcNumber.removeWildcards(ourECnumbers)
    
    onlyInTheirs = theirECnumbers.difference( ourECnumbers )
    
    #- keep only the ones occuring in NUKA.
    existingECs = []
    for theirECnumber in theirECnumbers:
        if theirECnumber in nukaECs:
            existingECs.append(theirECnumber)
    
    outdatedEcNumbers = theirECnumbers.difference( existingECs )
            
    output.append( "ECs not in any pathway: " + str(len(outdatedEcNumbers)) )
    
    onlyInTheirs.difference_update( outdatedEcNumbers )
    
    output.append( "Leaving ECs only in theirs: " + str(len(onlyInTheirs)) )
    
    output.append( "---------------------------------" )
    
    #- keep only the ones with a gene occuring in any of our 15 KEGG organisms found as "Bacillus subtilis"
    outdatedEcNumbers = set()
    for ecNumber in onlyInTheirs:
        ecNumberString = ecNumber.__str__()
        if not isKnownToKEGG(ecNumberString, organisms):
            outdatedEcNumbers.add(ecNumber)
    
    
    output.append( "ECs not in any of our organisms: " + str(len(outdatedEcNumbers)) )
    
    onlyInTheirs.difference_update( outdatedEcNumbers )
    
    output.append( "Leaving ECs only in theirs: " + str(len(onlyInTheirs)) )
    
    for ecNumber in onlyInTheirs:
        output.append(ecNumber.__str__())
    
        
    for line in output:
        print(line)
