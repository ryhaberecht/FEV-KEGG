from _collections import defaultdict
import string

class Gene(object):
    
    _digit_keeper = defaultdict(type(None))
    _digit_keeper.update({ord(c): c for c in string.digits})
    
    def __init__(self, content):
        """
        Gene as defined by a gene description file in KEGG GENE database.
        
        All attributes might be *None*, depending on whether they actually occur in the gene description! Occurence varies between organisms and sources.
        
        Parameters
        ----------
        content : str
            A multi-line gene description.
        
        Attributes
        ----------
        self.number : str
            The name of the gene, e.g. 'Acav_0021'.
        self.name : str
            Colloquial name of the gene product, e.g. 'ThrC'.
        self.definition : str
            Long description of the gene, e.g. '(GenBank) Homoserine dehydrogenase'.
        
        self.isProtein : bool
        self.ecNumbers : List[str]
        self.isEnzyme : bool
            If `isProtein` and has EC numbers.
        
        self.keggOrthologyNames : List[str]
            Names for each KEGG Orthology ID associated with this gene, e.g. 'homoserine dehydrogenase'.
        self.keggOrthologyIDs : List[str]
            Each KEGG Orthology ID assocaited with this gene, e.g. 'K00003'.
        self.keggOrthologies : List[Tuple[str, str, List[str]]]
            List of each associated KEGG Orthology entry, represented as a tuple of (ID, name, EC numbers), e.g. ('K00003', 'homoserine dehydrogenase', ['1.1.1.3']).
        
        self.organismAbbreviation : str
        self.organismName : str
        self.organismTnumber : str
            KEGG Onthology ID of the organism this gene belongs to, e.g. 'T01445'.
        
        self.pathways : List[Tuple[str, str]]
            List of organism-specific pathways this gene occurs in, represented as a tuple of (ID, name), e.g. ('aaa00260', 'Glycine, serine and threonine metabolism').
        
        self.positionFrom : int
            Nucleotide sequence position this gene starts at.
        self.positionTo : int
            Nucleotide sequence positon this gene ends at.
        self.positionIsComplement : bool
            Whether the positions are on the complement strand.
        
        self.aaseqLength : int
        self.aaseq : str
        
        self.ntseqLength : int
        self.ntseq : str
        """
        
        if isinstance(content, str):
            linesList = content.splitlines()
        
        try:
            # default values
            self.number = None
            self.isProtein = False
            self.organismTnumber = None
            
            self.name = None
            
            self.definition = None
            
            self.ecNumbers = list()
            self.keggOrthologyNames = list()
            self.keggOrthologyIDs = list()
            self.keggOrthologies = list()
            self.isEnzyme = False
            
            self.organismAbbreviation = None
            self.organismName = None
            
            self.pathways = list()
            
            self.positionFrom = None
            self.positionTo = None
            self.positionIsComplement = False
            
            self.aaseqLength = None
            self.aaseq = None
            
            self.ntseqLength = None
            self.ntseq = None
            
            # parse file data
            currentSection = None
            currentContent = None
            
            for line in linesList:
                
                if len(line) == 0 or line[0] == ' ': # section content
                    
                    if currentSection is not None:
                        currentContent.append(line.lstrip())
                    
                else: # section beginning
                    
                    # process previous section
                    lastSection = currentSection
                    
                    if lastSection is not None:
                        if lastSection == 'ENTRY':
                            
                            nextWords = currentContent[0].split()
                            self.number = nextWords[0]
                            self.isProtein = nextWords[1] == 'CDS'
                            self.organismTnumber = nextWords[2]
                            
                        elif lastSection == 'NAME':
                            
                            self.name = currentContent[0]
                        
                        elif lastSection == 'DEFINITION':
                            
                            self.definition = currentContent[0]
                        
                        elif lastSection == 'ORTHOLOGY':
                            
                            for contentLine in currentContent:
                                
                                keggOrthologyID, rest = contentLine.split('  ') # two spaces
                                self.keggOrthologyIDs.append(keggOrthologyID)
                                
                                restSplit = rest.split(' [EC:')
                                if len(restSplit) > 1: # has EC number and long name
                                    
                                    self.isEnzyme = True
                                    ecNumbers = restSplit[1][:-1].split(' ')
                                    longName = restSplit[0]
                                    
                                    self.ecNumbers.extend(ecNumbers) # one space
                                    self.keggOrthologyNames.append(longName)
                                
                                elif len(restSplit) == 1: # has only long name
                                    
                                    ecNumbers = None
                                    longName = restSplit[0]
                                    
                                    self.keggOrthologyNames.append(longName)
                                    
                                else: # has nothing
                                    
                                    ecNumbers = None
                                    longName = None
                                
                                self.keggOrthologies.append( (keggOrthologyID, longName, ecNumbers) )
                        
                        elif lastSection == 'ORGANISM':
                            
                            split = currentContent[0].split('  ') # two spaces
                            self.organismAbbreviation = split[0]
                            self.organismName = split[1]
                        
                        elif lastSection == 'PATHWAY':
                            
                            for contentLine in currentContent:
                                
                                pathwayID, pathwayName = contentLine.split('  ') # two spaces
                                
                                self.pathways.append( (pathwayID, pathwayName) )
                        
                        elif lastSection == 'POSITION':
                            
                            split = currentContent[0].split(':')
                            
                            if len(split) > 1: # there was a colon
                                split = split[1]
                            else:
                                split = split[0]
                            
                            split = split.split('..')
                            if 'complement' in currentContent[0]:
                                self.positionIsComplement = True
                            
                            if len(split) == 0 or len(split) == 1 and (split[0] == 'X' or split[0] == 'Y' or split[0] == 'Unknown'):
                                self.positionFrom = None
                                self.positionTo = None
                            
                            elif len(split) == 1:
                                self.positionFrom = int( split[0].translate(self.__class__._digit_keeper) )    
                                self.positionTo = None
                            
                            else:
                                self.positionFrom = int( split[0].translate(self.__class__._digit_keeper) )
                                self.positionTo = int( split[1].translate(self.__class__._digit_keeper) )
                        
                        elif lastSection == 'AASEQ':
                            
                            self.aaseqLength = int(currentContent[0])
                            self.aaseq = ''.join(currentContent[1:])
                            
                        elif lastSection == 'NTSEQ':
                            
                            self.ntseqLength = int(currentContent[0])
                            self.ntseq = ''.join(currentContent[1:])
                            
                    
                    # begin reading next section
                    split = line.split(maxsplit = 1)
                    
                    firstWord = split[0]
                    if firstWord.startswith('///'):
                        break
                    
                    else:
                        if len(split) > 1:
                            restLine = split[1]
                            currentContent = [restLine]
                        else:
                            currentContent = []
                            
                        currentSection = firstWord
        
        except:
            print( "Error while parsing a gene description into a KEGG.DataTypes.Gene object:" )
            print( content )
            raise
        
    def getGeneID(self) -> 'GeneID':
        from FEV_KEGG.Graph.Elements import GeneID
        return GeneID(self.organismAbbreviation + ':' + self.number)



class Substance(object):
    
    def __init__(self, content):
        """
        A compound/glycan found in KEGG pathways.
        
        Depending on whether this is a compound or a glycan, there are more attributes than listed below.
        
        Attributes
        ----------
        self.uniqueID : str
            Unique string identifying a substance, e.g. 'C00084'.
        self.description : str
            Human-readable description of a substance, e.g. 'Acetaldehyde;Ethanal'.
        self.name : str
            Short human-readable description of a substance, e.g. 'Acetaldehyde'. The shortest of all words in `self.description`.
        """
        try:
            # determine whether this is a compound or glycan
            if not isinstance(content, str):
                raise ValueError('Substance content is not a string.')
            
            linesList = content.splitlines()
            
            if len(linesList) <= 1:
                raise ValueError('Substance content has only one line.')
            
            # default values
            self.description = ''
            self.shortestDescription = ''
            self.firstDescription = ''
            
            firstLine = linesList[0]
            if 'Glycan' in firstLine:
                self._parseGlycan(linesList)
            elif 'Compound' in firstLine:
                self._parseCompound(linesList)
            else:
                raise ValueError('Substance type unknown.')
            
            self.name = self.firstDescription
            
        except:
            print( "Error while parsing a substance description into a KEGG.DataTypes.Substance object:" )
            print( content )
            raise
    
    def _parseCompound(self, linesList):
        
        # parse file data
        currentSection = None
        currentContent = None
        
        for line in linesList:
                
            if len(line) == 0 or line[0] == ' ': # section content
                
                if currentSection is not None:
                    currentContent.append(line.lstrip())
                
            else: # section beginning
                
                # process previous section
                lastSection = currentSection
                
                if lastSection is not None:
                    if lastSection == 'ENTRY':
                        
                        nextWords = currentContent[0].split()
                        self.uniqueID = nextWords[0]
#                         self.entry = currentContent[0]
                        
                    elif lastSection == 'NAME':
                        
                        self.description = ' \n'.join(currentContent)
                        self.shortestDescription = min(currentContent, key=len).replace(';','')
                        self.firstDescription = currentContent[0].replace(';','')
                    
                    elif lastSection == 'FORMULA':
                        
                        self.formula = currentContent[0]
                    
                    elif lastSection == 'EXACT_MASS':
                        
                        self.exact_mass = currentContent[0]
                    
                    elif lastSection == 'MOL_WEIGHT':
                        
                        self.mol_weight = currentContent[0]
                    
                        
                
                # begin reading next section
                split = line.split(maxsplit = 1)
                
                firstWord = split[0]
                if firstWord.startswith('///'):
                    break
                
                else:
                    if len(split) > 1:
                        restLine = split[1]
                        currentContent = [restLine]
                    else:
                        currentContent = []
                        
                    currentSection = firstWord
    
    def _parseGlycan(self, linesList):
        
        # parse file data
        currentSection = None
        currentContent = None
        
        for line in linesList:
                
            if len(line) == 0 or line[0] == ' ': # section content
                
                if currentSection is not None:
                    currentContent.append(line.lstrip())
                
            else: # section beginning
                
                # process previous section
                lastSection = currentSection
                
                if lastSection is not None:
                    if lastSection == 'ENTRY':
                        
#                         self.entry = currentContent[0]
                        nextWords = currentContent[0].split()
                        self.uniqueID = nextWords[0]
                        
                    elif lastSection == 'COMPOSITION':
                        
                        self.description = ' \n'.join(currentContent)
                        self.shortestDescription = min(currentContent, key=len).replace(';','')
                        self.firstDescription = currentContent[0].replace(';','')
#                         self.composition = self.description
                    
                    elif lastSection == 'MASS':
                        
                        self.mass = currentContent[0]
                    
                        
                
                # begin reading next section
                split = line.split(maxsplit = 1)
                
                firstWord = split[0]
                if firstWord.startswith('///'):
                    break
                
                else:
                    if len(split) > 1:
                        restLine = split[1]
                        currentContent = [restLine]
                    else:
                        currentContent = []
                        
                    currentSection = firstWord


class EcEnzyme(object):
    
    def __init__(self, content):
        """
        An enzyme found in KEGG pathways, defined by its EC number.
        
        Attributes
        ----------
        self.uniqueID : str
            Unique string identifying the EC number, e.g. '4.1.2.48'.
        self.description : str
            Human-readable description of the EC number, e.g. 'low-specificity L-threonine aldolase;LtaE'.
        self.name : str
            Short human-readable name of the EC number, e.g. 'LtaE'. The shortest of all words in `self.description`.
        """
        try:
            # determine whether this is a compound or glycan
            if not isinstance(content, str):
                raise ValueError('Enzyme content is not a string.')
            
            linesList = content.splitlines()
            
            if len(linesList) <= 1:
                raise ValueError('Enzyme content has only one line.')
            
            # default values
            self.description = ''
            self.shortestDescription = ''
            self.firstDescription = ''
            self.reaction = ''
            
            # parse file data
            currentSection = None
            currentContent = None
            
            for line in linesList:
                    
                if len(line) == 0 or line[0] == ' ': # section content
                    
                    if currentSection is not None:
                        currentContent.append(line.lstrip())
                    
                else: # section beginning
                    
                    # process previous section
                    lastSection = currentSection
                    
                    if lastSection is not None:
                        if lastSection == 'ENTRY':
                            
    #                         self.entry = currentContent[0]
                            nextWords = currentContent[0].split()
                            self.uniqueID = nextWords[1]
                            
                        elif lastSection == 'NAME':
                            
                            self.description = ' \n'.join(currentContent)
                            self.shortestDescription = min(currentContent, key=len).replace(';','')
                            self.firstDescription = currentContent[0].replace(';','')
    #                         self.composition = self.description
                        
                        elif lastSection == 'REACTION':
                            
                            self.reaction = ' \n'.join(currentContent)
                        
                            
                    
                    # begin reading next section
                    split = line.split(maxsplit = 1)
                    
                    firstWord = split[0]
                    if firstWord.startswith('///'):
                        break
                    
                    else:
                        if len(split) > 1:
                            restLine = split[1]
                            currentContent = [restLine]
                        else:
                            currentContent = []
                            
                        currentSection = firstWord
                        
            self.name = self.firstDescription
            
        except:
            print( "Error while parsing an enzyme description into a KEGG.DataTypes.EcEnzyme object:" )
            print( content )
            raise
    