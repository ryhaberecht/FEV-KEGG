from builtins import str
import re
from typing import List, Iterable

from FEV_KEGG.KEGG.DataTypes import Gene
from FEV_KEGG.Util import Util
from FEV_KEGG import settings


class Element(object):
    def __init__(self, uniqueID: str):
        """
        Generic graph element with a `uniqueID`.
        
        Comparable (==, !=, <, >, <=, >=) and hashable by this unique ID. Converting to a string returns the `uniqueID`.
        
        Parameters
        ----------
        uniqueID : str
            String uniquely identifying this element among all other possible elements.
        
        Attributes
        ----------
        self.uniqueID : str
            Unique element ID.
        """
        self.uniqueID = uniqueID
    
    def getUrl(self):
        """
        Get the link to KEGG for this EC number.
        
        Returns
        -------
        str
            URL to KEGG.
        """
        return "http://kegg.jp/dbget-bin/www_bget?" + self.uniqueID
    
    def getRestUrl(self):
        """
        Get the link to KEGG's REST-API for this EC number.
        
        Essentially the same as :func:`getUrl`, but meant to be read by machines, therefore no eye-candy.
        
        Returns
        -------
        str
            URL to KEGG's REST-API
        """
        return "http://rest.kegg.jp/get/" + self.uniqueID
    
    def toHtml(self, short = False, noTd = False):
        """
        Get the Element's string representation surrounded by its URL as an HTML line.
        """
        if self.name is None or short is True:
            if noTd is True:
                return "<a target='_blank' href='" + self.getUrl() + "'>" + self.__str__() + "</a>"
            else:
                return "<td><a target='_blank' href='" + self.getUrl() + "'>" + self.__str__() + "</a></td><td></td>"
        else:
            if noTd is True:
                return "<a target='_blank' href='" + self.getUrl() + "'>" + self.__str__() + "</a>&nbsp;(" + self.name + ")"
            else:
                return "<td><a target='_blank' href='" + self.getUrl() + "'>" + self.__str__() + "</a></td><td>(" + self.name + ")</td>"
    
        
    def __str__(self):
        if settings.printElementUrl:
            return str(self.uniqueID) + ' (' + self.getUrl() + ')'
        else:
            return self.uniqueID
    
    def __repr__(self):
        return self.__str__()
        
    def __eq__(self, other):
        if isinstance(self, other.__class__):
            return self.uniqueID == other.uniqueID
        return False
        
    def __ne__(self, other):
        return not self == other
    
    def __hash__(self):
        return self.uniqueID.__hash__()
    
    def __lt__(self, other):
        return self.uniqueID < other.uniqueID
    
    def __gt__(self, other):
        return self.uniqueID > other.uniqueID
    
    def __le__(self, other):
        return self.uniqueID <= other.uniqueID
    
    def __ge__(self, other):
        return self.uniqueID >= other.uniqueID
    


class DrugIdError(Exception):
    """
    Raised if a :class:`SubstanceID` is created from a drug ID, because only compounds and glycans are useful in our model of metabolism.
    """
    pass

class SubstanceID(Element):
    REGEX_PATTERN = re.compile('^C|G[0-9]{5}$')
    
    def __init__(self, keggSubstanceID: 'C01102'):
        """
        Represents a substrate/product of metabolism by compound/glycan ID from KEGG, eg. 'C01102' or 'G00160'.
        
        Parameters
        ----------
        keggSubstanceID : str
            Unique ID of the compound or glycan.
        description : str, optional
            Descriptive chemical name of the compound/glycan.
        
        Attributes
        ----------
        self.keggCompoundID : str
            Unique compound/glycan ID.
        self.description : str
            Descriptive chemical name of the compound/glycan. May likely be *None*. Usually a list of synonymous names.
        self.name : str
            Short chemical name of the compound/glycan. May likely be *None*. Is the shortest name occuring in `description`.
        
        Raises
        ------
        DrugIdError
            Drug IDs, eg. D08603, raise a DrugIdError, because only compounds and glycans are useful in our model of metabolism. Use the synonymous Compound ID instead.
        
        Note
        ----
        This does not check if the compound/glycan actually exists in KEGG! You will find out eventually when trying to retrieve information about it.
        
        See Also
        --------
        FEV_KEGG.Graph.SubstanceGraphs.SubstanceGraph.addSubstanceDescriptions : The function to download and add `self.description`, and `self.name`.
        """
        if keggSubstanceID[0] == 'D':
            raise DrugIdError('Drug IDs are not accepted, as there are usually accompanied by a synonymous Compound ID.')
        
        if self.__class__.REGEX_PATTERN.match(keggSubstanceID) is None: # wrong format
            raise ValueError('Compound/Glycan ID not formatted correctly: ' + keggSubstanceID)
        
        Element.__init__(self, keggSubstanceID)
        self.keggCompoundID = self.uniqueID
        self.description = None
        self.name = None




class ReactionID(Element):
    
    def __init__(self, keggReactionID: 'R01899'):
        """
        Represents a reaction of metabolism by reaction ID from KEGG, eg. 'R01899'.
        
        Parameters
        ----------
        keggReactionID : str
            Unique ID of the reaction.
        
        Attributes
        ----------
        self.keggReactionID : str
            Unique reaction ID.
        
        Note
        ----
        This does not check if the reaction actually exists in KEGG! You will find out eventually when trying to retrieve information about it.
        """
        Element.__init__(self, keggReactionID)
        self.keggReactionID = self.uniqueID



class Enzyme(Element):

    def __init__(self, organismAbbreviation: 'eco', geneName: 'b0004', ecNumberStrings: List[str], name: 'thrC' = None, description: '(RefSeq) hydrogenase 4, subunit' = None):
        """
        Represents an enzyme of metabolism.
        
        It has exactly one GeneID, which is its unique identifier.
        
        Parameters
        ----------
        organismAbbreviation : str
            Abbreviation string of the organism this enzyme belongs to, as known to KEGG, e.g. 'eco'. Must obviously be unique and existant in KEGG.
        geneName : str
            Name of the gene which represents this enzyme, e.g. 'b0004'. Will be combined with `organismAbbreviation` to form the unique :class:`GeneID`. Thus, must be unique within the organism.
        ecNumberStrings : List[str]
            List of strings representing the EC numbers associated with this enzyme. Will be split and parsed into :class:`EcNumber` objects.
        name : str, optional
            Colloquial name of this enzyme, e.g. 'thrC'. This is not used for automatic identification, you may make it *None*.
        description : str, optional
            Full description of this enzyme from KEGG, e.g. '(RefSeq) hydrogenase 4, subunit'. This is not used for automatic identification, you may make it *None*.
        
        Attributes
        ----------
        self.organismAbbreviation : str
        self.geneName : str
        self.geneID : GeneID
        self.name : str
        self.ecNumbers : Set[EcNumber]
        self.description : str
            
        Raises
        ------
        ValueError
            If `organismAbbreviation` and `geneName` do not form a valid gene ID. Or if any of the EC numbers in `ecNumberStrings` is not a valid EC number.
        
        Note
        ----
        This does not check if the organism, gene ID, EC numbers, or anything else actually exist in KEGG! You will find out eventually when trying to retrieve information about them.
        """
        # build subclasses
        # GeneID
        geneID = GeneID(organismAbbreviation + ':' + geneName)
        # EcNumbers
        ecNumbers = set()
        for ecNumberString in ecNumberStrings:
            ecNumber = EcNumber(ecNumberString)
            ecNumbers.add(ecNumber)
        
        # determine unique ID
        Element.__init__(self, geneID.__str__())
        
        # save object attributes
        self.organismAbbreviation = organismAbbreviation
        self.geneID = geneID
        self.geneName = geneName
        if name is not None and name.__eq__(geneName):
            self.name = None
        else:
            self.name = name
        self.ecNumbers = ecNumbers
        # replace useless substrings
        if description is not None:
            description = description.replace('(RefSeq) ', '')
        self.description = description
    
    def getEcNumbersString(self):
        """
        EC numbers associated with this enzyme as a string.
        
        Returns
        -------
        str
            EC numbers associated with this enzyme in a string, eg. '1.2.3.4, 2.3.4.5'
        """
        strings = []
        for ecNumber in self.ecNumbers:
            strings.append(ecNumber.__str__())
            
        return ', '.join(strings)
    
    @classmethod
    def fromGene(cls, gene: Gene) -> 'Enzyme':
        """
        Creates an :class:`Enzyme` from a :class:`FEV_KEGG.KEGG.DataTypes.Gene`.
        
        Parameters
        ----------
        gene : Gene
            Gene object, retrieved and parsed from KEGG GENE at some point.
        
        Returns
        -------
        Enzyme
            An enzyme object.
        
        Raises
        ------
        ValueError
            If `organismAbbreviation` and `geneName` do not form a valid gene ID. Or if any of the EC numbers in `ecNumberStrings` is not a valid EC number.
        """
        return cls(organismAbbreviation = gene.organismAbbreviation, geneName = gene.number, ecNumberStrings = gene.ecNumbers, name = gene.symbol, description = gene.name)
    
    def __lt__(self, other):
        
        # sort by EC number first
        selfEcList = list(self.ecNumbers)
        otherEcList = list(other.ecNumbers)

        if selfEcList == otherEcList:
            # then by gene ID
            return self.uniqueID < other.uniqueID

        else:
            return selfEcList < otherEcList
    
    def __gt__(self, other):
        
        # sort by EC number first
        selfEcList = list(self.ecNumbers)
        otherEcList = list(other.ecNumbers)

        if selfEcList == otherEcList:
            # then by gene ID
            return self.uniqueID > other.uniqueID

        else:
            return selfEcList > otherEcList
    
    def __le__(self, other):
        
        # sort by EC number first
        selfEcList = list(self.ecNumbers)
        otherEcList = list(other.ecNumbers)

        if selfEcList == otherEcList:
            # then by gene ID
            return self.uniqueID <= other.uniqueID

        else:
            return selfEcList <= otherEcList
    
    def __ge__(self, other):
        
        # sort by EC number first
        selfEcList = list(self.ecNumbers)
        otherEcList = list(other.ecNumbers)

        if selfEcList == otherEcList:
            # then by gene ID
            return self.uniqueID >= other.uniqueID

        else:
            return selfEcList >= otherEcList
    

class EnzymeComplete(Enzyme):
    
    def __init__(self, gene: Gene):
        """
        Represents an enzyme of metabolism, saving the original underlying gene description `gene` for later manual use.
        
        The underlying gene description is usually not necessary, use the parent class to save memory space.
        
        Parameters
        ----------
        gene : Gene
            Gene object, retrieved and parsed from KEGG GENE at some point. Will be kept in memory in the *gene* attribute.
        
        Attributes
        ----------
        self.gene : :class:`FEV_KEGG.KEGG.DataTypes.Gene`
            Original underlying gene description.
        
        Raises
        ------
        ValueError
            See parent class.
        """
        super().__init__(gene.organismAbbreviation, gene.number, gene.symbol, gene.ecNumbers)
        self.gene = gene
    
    
class EcNumber(Element):
    WILDCARD = '-'
    REGEX_PATTERN = re.compile('^[1-7]\.(([1-9][0-9]{0,1})|\-)\.(((?<!\-\.)([1-9][0-9]{0,1}))|\-)\.(((?<!\-\.)([1-9][0-9]{0,2}))|\-)$')
    
    def __init__(self, ecNumberString: '4.2.3.1'):
        """
        Represents an enzyme of metabolism by EC number, e.g. '4.2.3.1'.
        
        Parameters
        ----------
        ecNumberString : str
            EC number represented as a string. Will be checked for correct formatting!
        
        
        Attributes
        ----------
        self.ecNumberString : str
            E.g. '4.2.3.-'.
        self.ecNumberLevels : List[str]
            E.g. ['4', '2', '3', '-'].
        self.ecNumberLevelsInteger : List[int]
            E.g. [4, 2, 3, -1]. A wildcard is translated to -1.
        self.description : str
            Descriptive name of the enzymes behind this EC number. May likely be *None*. Usually a list of synonymous names.
        self.name : str
            Short name of the enzymes behind this EC number. May likely be *None*. Is the shortest name occuring in `description`.
        self.reaction : str
            IUBMB string describing the reaction formula. May likely be *None*.
        
        Raises
        ------
        ValueError
            If EC number is not formatted correctly.
        
        See Also
        --------
        FEV_KEGG.Graph.SubstanceGraphs.SubstanceEcGraph.addEcDescriptions : The function to download and add `self.description`, `self.name`, and `self.reaction`.
        """
        if self.__class__.REGEX_PATTERN.match(ecNumberString) is None: # wrong format
            raise ValueError('EC number not formatted correctly: ' + ecNumberString)
        
        # determine unique ID
        Element.__init__(self, ecNumberString)
        
        # save object attributes
        self.ecNumberString = self.uniqueID
        self.ecNumberLevels = self.ecNumberString.split('.')
        self._ecNumberLevelsInteger = [-1 if level == EcNumber.WILDCARD else int(level) for level in self.ecNumberLevels]
        self.description = None
        self.name = None
        self.reaction = None
    
    @classmethod
    def fromArray(cls, ecNumberLevels: Iterable) -> 'EcNumber':
        """
        Creates EcNumber object from single EC number levels.
        
        Parameters
        ----------
        ecNumberLevels : Iterable
            Iterable of the EC number levels, can be int or str. For a wildcard, obviously only str is reasonable.
        
        Raises
        ------
        ValueError
            If the resulting EC number is not formatted correctly.
        """
        return cls('.'.join(ecNumberLevels))
    
    @property
    def ecNumberLevelsInteger(self) -> List[int]:
        if not hasattr(self, '_ecNumberLevelsInteger'):
            self._ecNumberLevelsInteger = [-1 if level == EcNumber.WILDCARD else int(level) for level in self.ecNumberLevels]
        return self._ecNumberLevelsInteger
    
    def contains(self, ecNumber: 'EcNumber') -> bool:
        """
        Check whether this EC number is a superset of `ecNumber`, made possibly by the wildcard.
        
        Parameters
        ----------
        ecNumber : EcNumber
            The EC number to compare against.
        
        Returns
        -------
        bool
            *True*, if the other EC number is part of the set of EC numbers defined by wildcard dashes in the levels of this EC number.
            For example 1.2.3.- contains 1.2.3.1 up to 1.2.3.999, but 1.2.3.4 can only contain itself.
        """
        selfLevels = self.ecNumberLevels
        otherLevels = ecNumber.ecNumberLevels
        
        for i in range(0, 4):
            selfNumber = selfLevels[i]
            otherNumber = otherLevels[i]
            if selfNumber != EcNumber.WILDCARD and selfNumber != otherNumber: # current level does not match AND is has no wildcard '-' in this EC number
                return False
        
        return True
    
    
    def matchingLevels(self, ecNumber: 'EcNumber', wildcardMatchesNumber = True) -> int:
        """
        Determines the number of levels which match between this EC number and `ecNumber`.
        
        This could act as a coarse distance measure for EC numbers.
        
        Parameters
        ----------
        ecNumber : EcNumber
            The EC number to compare against.
        wildcardMatchesNumber : bool, optional
            If *True*, a wildcard acts as a sure match: '1.-.-.-'.matchingLevels('1.2.3.4') = 4.
            If *False*, a wildcard only matches another wildcard.
        
        Returns
        -------
        int
            Number of consecutive levels that match, if any, starting with the first (leftmost).
            '1.2.3.4'.matchingLevels('1.2.6.7') = 2 because the first two levels match consecutively.
            '1.2.3.4'.matchingLevels('2.2.3.4') = 0 because the very first level does not match.
        """
        matchingLevels = 0
        
        selfLevels = self.ecNumberLevels
        otherLevels = ecNumber.ecNumberLevels
        
        for i in range(0, 4):
            selfNumber = selfLevels[i]
            otherNumber = otherLevels[i]
            
            if wildcardMatchesNumber == True:
                if selfNumber == EcNumber.WILDCARD or otherNumber == EcNumber.WILDCARD or selfNumber == otherNumber: # current level matches OR is a wildcard
                    matchingLevels += 1
                else:
                    return matchingLevels
            else:
                if selfNumber == otherNumber: # current level matches
                    matchingLevels += 1
                else:
                    return matchingLevels
        
        return matchingLevels
    
    def hasWildcard(self) -> bool:
        """
        Whether this EC number contains a wildcard.
        
        Returns
        -------
        bool
            *True* if this EC number contains a wildcard (-) at any level, otherwise, returns *False*.
        """
        for level in self.ecNumberLevels:
            if level == EcNumber.WILDCARD:
                return True
        return False
    
    @staticmethod
    def removeWildcards(ecNumbers: Iterable) -> Iterable:
        """
        Remove EC numbers containing wildcards from an Iterable.
        
        Parameters
        ----------
        ecNumbers : Iterable[EcNumber]
            The EcNumber objects to check for wildcards.
        
        Returns
        -------
        Iterable[EcNumber]
            A new Iterable of the same type, containing only EC numbers which do **not** have a wildcard (-) anywhere. This does not deduplicate EC numbers.
        """
        validECnumbers = []
        for ecNumber in ecNumbers:
            if not ecNumber.hasWildcard():
                validECnumbers.append(ecNumber)
        
        return ecNumbers.__class__(validECnumbers)
    
    @staticmethod
    def insertWildcards(ecNumbers: Iterable, keepLevels = 3, allowHigherWildcards = True, returnSet = True, deduplicateList = False) -> Iterable:
        """
        Turns EC numbers without wildcards into EC numbers with wildcards. 
        
        Returning them in a list preserves order.
        
        Parameters
        ----------
        ecNumbers : Iterable
            The EcNumber objects to abstract using wildcards.
        keepLevels : int, optional
            The first x levels of each EC number are kept intact. If `keepLevels` == 3, turns 1.2.3.4 into 1.2.3.-. Only 1, 2, 3, and 4 are allowed. EC numbers already containing wildcards are left unchanged.
        allowHigherWildcards : bool, optional
            If *False* and there is a wildcard in a level above 'keepLevels' (e.g. 3):, 1.2.3.4 -> 1.2.3.- and 2.3.4.- -> 2.3.4.-, but 3.4.-.- is removed completely.
        returnSet : bool, optional
            If *True*, returns results in a set. Takes precedence over 'deduplicateList', as sets automatically deduplicate.
        deduplicateList : bool, optional
            If *True*, result list is deduplicated before returning, preserving order.
        
        Returns
        -------
        Iterable
            Either a list or a set of abstracted EC numbers.
        
        Raises
        ------
        ValueError
            If `keepLevels` is not one of [1, 2, 3, 4].
        """
        if not keepLevels in [1, 2, 3, 4]:
            raise ValueError('Can not keep ' + str(keepLevels) + ' levels, there are only 1, 2, 3, or 4.')
        
        filtered = []
        
        for ecNumber in ecNumbers:
            
            levels = ecNumber.ecNumberLevels
            
            filteredLevels = []
            for i in range(0, keepLevels):
                
                level = levels[i]
                
                # check for higher wildcards
                if allowHigherWildcards is False and level == EcNumber.WILDCARD:
                    filteredLevels = None
                    break
                
                else:
                    filteredLevels.append(level)
                
            if filteredLevels is None: # higher wildcard found but disallowed
                continue
            
            else: # pad with wildcards
                for _ in range(4, keepLevels, -1):
                    filteredLevels.append(EcNumber.WILDCARD)
            
            filtered.append( EcNumber.fromArray(filteredLevels) )
        
        if returnSet is True:
            return set( filtered )
            
        if deduplicateList is True:
            filtered = Util.deduplicateList(filtered, preserveOrder = True)
        
        return filtered
    
    def addDescription(self):
        """
        Query KEGG and add further description to this EC number.
        
        Warnings
        --------
        Much slower than doing :func:`addEcDescriptions` for several EC numbers in bulk!
        """
        from FEV_KEGG.KEGG import Database
        
        ecNumberIdToEcEnzyme = Database.getEcEnzymeBulk([self])
        ecEnzyme = ecNumberIdToEcEnzyme.get(self.uniqueID)
        if ecEnzyme is not None:
            self.description = ecEnzyme.description
            self.name = ecEnzyme.name
            self.reaction = ecEnzyme.reaction
    
    @staticmethod
    def addEcDescriptions(ecNumbers: Iterable):
        """
        Query KEGG for further descriptions and add them to each EC number in `ecNumbers`.
        """
        from FEV_KEGG.KEGG import Database
        
        ecNumberIdToEcEnzyme = Database.getEcEnzymeBulk(ecNumbers)
        for ecNumber in ecNumbers:
            ecEnzyme = ecNumberIdToEcEnzyme.get(ecNumber.uniqueID)
            if ecEnzyme is not None:
                ecNumber.description = ecEnzyme.description
                ecNumber.name = ecEnzyme.name
                ecNumber.reaction = ecEnzyme.reaction
    
    def __lt__(self, other):
        return self.ecNumberLevelsInteger < other.ecNumberLevelsInteger
    
    def __gt__(self, other):
        return self.ecNumberLevelsInteger > other.ecNumberLevelsInteger
    
    def __le__(self, other):
        return self.ecNumberLevelsInteger <= other.ecNumberLevelsInteger
    
    def __ge__(self, other):
        return self.ecNumberLevelsInteger >= other.ecNumberLevelsInteger

    
    
class GeneID(Element):
    
    REGEX_PATTERN = re.compile('^[a-z]{3,4}:[a-zA-Z0-9_\-\.]+$')
    
    def __init__(self, geneIDString: 'eco:b0004'):
        """
        Represents am enzyme of metabolism by gene ID, e.g. 'eco:b0004'.
        
        Parameters
        ----------
        geneIDString : str
            Gene ID represented by a string, e.g. 'eco:b0004'. Will be checked for correct formatting!
        
        Attributes
        ----------
        self.geneIDString : str
        
        Raises
        ------
        ValueError
            If gene ID is not formatted correctly.
        """
        # check input
        if self.__class__.REGEX_PATTERN.match(geneIDString) is None: # wrong format
            raise ValueError('Gene ID not formatted correctly: ' + geneIDString)
        
        # determine unique ID
        Element.__init__(self, geneIDString)
        
        # save object attributes
        self.geneIDString = self.uniqueID
        
    @property
    def organismAbbreviation(self) -> str:
        """
        Returns
        -------
        str
            'eco' from 'eco:b0004'.
        """
        geneIDSplit = self.geneIDString.split(':')
        organismAbbreviation = geneIDSplit[0]
        return organismAbbreviation
    
    @property
    def geneName(self) -> str:
        """
        Returns
        -------
        str
            'b0004' from 'eco:b0004'.
        """
        geneIDSplit = self.geneIDString.split(':')
        geneName = geneIDSplit[1]
        return geneName
    
        
    
class KeggOrthologyID(Element):
    
    REGEX_PATTERN = re.compile('^K[0-9]{5}$')
    
    def __init__(self, keggOrthologyIDString: 'K01733'):
        """
        Represents an enzyme of metabolism by KEGG Orthology ID.
        
        Parameters
        ----------
        keggOrthologyIDString : str
            String representation of a KEGG Orthology ID. Will be checked for correct formatting!
        
        Attributes
        ----------
        self.keggOrthologyIDString : str
        
        Raises
        ------
        ValueError
            If KEGG Orthology ID is not formatted correctly.
        """
        # check input
        if self.__class__.REGEX_PATTERN.match(keggOrthologyIDString) is None: # wrong format
            raise ValueError('KEGG Orthology ID not formatted correctly: ' + keggOrthologyIDString)
        
        # determine unique ID
        Element.__init__(self, keggOrthologyIDString)
        
        # save object attributes
        self.keggOrthologyIDString = self.uniqueID
        