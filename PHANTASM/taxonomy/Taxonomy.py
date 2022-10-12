# Author: Joseph S. Wirth
# Last edit: September 29, 2022

from __future__ import annotations
import re, copy, random, textdistance, os
from Bio.Entrez import Parser
from PHANTASM.taxonomy.TaxRank import TaxRank
from PHANTASM. utilities import ncbiEfetchById, ncbiIdsFromSearchTerm, \
                      ncbiSummaryFromIdList, coerceEntrezToNormal, parseCsv


class Taxonomy:
    """ Taxonomy:
            This class is designed to store taxonomic information from NCBI and 
            LPSN in a convenient structure that can be queried both quickly and
            easily. Corresponding data from NCBI's assembly database is also 
            saved in the object. It is capable of automatically selecting refe-
            rence genomes to use for phylogenomic analyses.
    """
    ### constant used for comparing ranks to the taxonomy floor/ceiling
    SPECIES = TaxRank('species')
    DOMAIN  = TaxRank('domain')




    ### CONSTRUCTOR
    def __init__(self, taxid, ncbiName:str, rank, parentTax:Taxonomy=None, \
                 isExternal:bool=False, taxidSynonyms:dict=dict()) -> None:
        """ __init__:
                Accepts an NCBI taxonomy id as a string or an int, a string 
                containing the name from NCBI, a string or TaxRank object 
                indicating the taxonomic rank, the Taxonomy object of its par-
                ental taxonomic rank (default = None), and boolean indicating 
                if the object can be considered "external" to the focal clade
                (used for picking candidate outgroups) as inputs. Creates an 
                object with the specified information. Does not return.
        """
        # initialize member variables from inputs
        self.taxid:str = str(taxid)     # ensure taxid is a string
        self.sciName:str = str(ncbiName)  # protect against Bio.Entrez types
        self.ncbiName:str = str(ncbiName) # protect against Bio.Entrez types
        self.rank:TaxRank = rank
        self.descendantsD:dict = {}
        self.parent:Taxonomy = parentTax
        self.typeMaterial:Taxonomy = None
        self.isExternal:bool = isExternal
        self.ncbiNameIsCorrect:bool = True
        self.assemblyFtp:str = None
        self.assemblyAccn:str = None
        self.assemblyFromType:bool = False
        self.assemblyStrain:str = None
        self.refSeqCategory:str = None
        self.assCoverage:float = 0.0
        self.assIsComplete:bool = False
        self.allAssIds:list = []
        self.taxidSynonyms:dict = taxidSynonyms  # this dict is shared!!!
        
        # ensure self.rank is a TaxRank object
        if type(self.rank) is not TaxRank:
            self.rank = TaxRank(self.rank)




    ### OVERLOADED FUNCTIONS    
    def __repr__(self) -> str:
        """ __repr__:
                Accepts no inputs. Returns a string indicating the name, rank, 
                and NCBI taxonomy id of the object.
        """
        # construct and return a string containing the name, rank, and taxid
        out = self.sciName
        out += ' (rank: ' + str(self.rank)
        out += ', txid: ' + self.taxid + ')'
        return out

    
    def __str__(self) -> str:
        """ __str__:
                Accepts no inputs. Returns a string indicating the name, rank, 
                and NCBI taxonomy id of the object as well as the names, ranks,
                and NCBI taxonomy ids of the Taxonomy objects immediately below
                it (eg. the genera within a family).
        """
        # constants
        GAP_SIZE = 4
        TYPE_PREFIX = '* '

        # start by getting the base string from __repr__
        out = self.__repr__()

        if self.rank > Taxonomy.SPECIES:
            # add a new line if lower ranks are being included
            out += '\n'

            # add immediate descendant info if any are present
            if len(self.descendantsD) > 0:
                # report the taxonomic level below the current
                plural = self.rank.getRankBelow().plural()
                out += '  sub-taxa (' + plural + '):\n'

                # add a line for each descendant with its info
                for key in self.descendantsD.keys():
                    subEntry = self.descendantsD[key]

                    if self.typeMaterial is not None:
                        # annotate the type descendant
                        if subEntry == self.typeMaterial:
                            out += ' '*(GAP_SIZE - len(TYPE_PREFIX))
                            out += TYPE_PREFIX
                            out += subEntry.__repr__()
                        
                        # otherwise just print descendant
                        else:
                            out += ' '*GAP_SIZE
                            out += subEntry.__repr__()
                    
                    # otherwise just print descendant
                    else:
                        out += ' '*GAP_SIZE
                        out += subEntry.__repr__()

                    out += '\n'
                
                # remove the extra newline character from the end of the string 
                out = out[:-1]
            
            # if no descendant taxa present, then report this.
            else:
                out += '  Descendant taxa have not been populated.'
        
        return out
    

    def __len__(self) -> int:
        """ __len__:
                Accepts no inputs. Returns an integer indicating the total 
                number of taxa contained within the object. The minimum length
                of a Taxonomy object is '1'.
        """
        # constants
        SELF_LEN = 1

        # add 1 to the length of all the children and return
        return SELF_LEN + len(self.getAllDescendants())


    def __contains__(self, item) -> bool:
        """ __contains__:
                Accepts either a Taxonomy object, a taxid (int or str), or a
                scientific name (str) as input. Returns a boolean indicating if
                the input is present within the calling object.
        """
        # constants
        ERR_MSG_PREFIX = "'in <Taxonomy>' requires Taxonomy, taxid " + \
                         "(int or str), or sciName (str) as left operand, not "
        ERR_GREP_FIND = r"^<class '([^']+)'>$"
        ERR_GREP_REPL = r"\1"

        # if the item is a Taxonomy object
        if type(item) is Taxonomy:
            # attempt lookup by its taxid; True if found, False if not
            if self.getDescendantByTaxId(item.taxid):
                return self.getDescendantByTaxId(item.taxid) == item
            else:
                return False
        
        # if the item is an integer, then it must be a taxid
        elif type(item) is int:
            # attempt lookup by the taxid
            if self.getDescendantByTaxId(item):
                return True
            else:
                return False

        # if the item is a string
        elif type(item) is str:
            # first try to coerce it to an integer
            try:
                # if it can be coerced, then it must be a taxid; attempt lookup
                item = int(item)
                if self.getDescendantByTaxId(item):
                    return True
                else:
                    return False
            
            # if the item is not an integer, then it must be a sciName
            except ValueError:
                if self.getDescendantBySciName(item):
                    return True
                else:
                    return False
        
        else:
            error = re.sub(ERR_GREP_FIND, ERR_GREP_REPL, str(type(item)))
            raise ValueError(ERR_MSG_PREFIX + error)


    def __hash__(self) -> int:
        return int(self.taxid)


    def __eq__(self, other:Taxonomy) -> bool:
        """ __eq__:
                Directly compares two Taxonomy objects. Checks if all member 
                variables (except other Tax objects) match. Then recursively
                checks all of the children in each of the two objects. If 
                everything is identical, then returns True. Otherwise, returns
                False.

                Important note: inequivalent objects can still have equivalent
                subtaxa.
        """
        # if other is None, then False
        if other is None:
            return False
        
        # first check if the member variables are equal
        if not Taxonomy.__eqVars(self, other):
            return False
        
        # if at the species-level, then compare the two lists of strain names
        if self.rank == Taxonomy.SPECIES:
            if self.typeMaterial != other.typeMaterial:
                return False
        
        # otherwise if either's type material is None
        elif self.typeMaterial is None or other.typeMaterial is None:
            # if both aren't none, then false
            if self.typeMaterial is not None or other.typeMaterial is not None:
                return False
        
        # otherwise, if both types not none, then compare taxids (objects compared already)
        elif self.typeMaterial.taxid != other.typeMaterial.taxid: 
            return False
    
        # if either's parent is none
        if self.parent is None or other.parent is None:
            # if both aren't none, then false
            if self.parent is not None or other.parent is not None:
                return False

        # if parents not None, then compare taxids
        elif self.parent.taxid != other.parent.taxid:
            return False
        
        # check the children
        sChildren = self.getChildren()
        oChildren = other.getChildren()

        # if the number of children don't match, then false
        if len(sChildren) != len(oChildren):
            return False
    
        # recurse through the children
        for key in sChildren:
            # false if the key is not available
            if key not in oChildren.keys():
                return False

            # get the children           
            sChild = sChildren[key]
            oChild = oChildren[key]

            # return false if there is ever a mismatch
            if sChild != oChild:
                return False

        # if everything matched, then true
        return True


    def __eqVars(lhs:Taxonomy, rhs:Taxonomy) -> bool:
        """ eqVars:
                Accepts two Taxonomy objects as input. Compares all the member
                variables except those that are themselves Taxonomy objects. 
                Returns a boolean indicating whether or not the members contain
                equivalent values.
        """
        # constants
        DESC = 'descendantsD'
        PRNT = 'parent'
        TYPE = 'typeMaterial'
        TAX_OBJ_KEYS = [DESC, PRNT, TYPE, 'taxidSynonyms']
        RANK = 'rank'

        # get dicts of the member variables
        lhsVar = vars(lhs)
        rhsVar = vars(rhs)

        # go through each member variable (except those that are Taxonomy)
        for key in lhsVar:
            # don't compare dicts of taxonomy objects directly
            # doing so will exceed max recursive depth
            if key not in TAX_OBJ_KEYS:
                # return false if a mismatch in value
                if lhsVar[key] != rhsVar[key]:
                    return False

        # if all keys matched, then return true
        return True
  

    def __ne__(self, other:Taxonomy) -> bool:
        """ __ne__:
                Directly compares two Taxonomy objects. Uses __eq__ to compare 
                and returns the inverted the result.
        """
        return not self == other


    def __deepcopy__(self, memo:dict) -> Taxonomy:
        """ __deepcopy__:
                Accepts a dictionary as input. Returns a deepcopy of the calli-
                ng object. This function may seem redundant, but it allows dee-
                pcopying of collections of Taxonomies, such as lists, dicts, or
                sets. Although unused, the 'memo' input is required to for it 
                properly interface with the builtin function.
        """
        return self.copy()




    ### COPYING, SAVING, AND LOADING OBJECTS
    def copy(self) -> Taxonomy:
        """ copy:
                Accepts no inputs. Converts the calling object into the "memo"
                format and then rebuilds a new Taxonomy object from that memo.
                Returns a deep copy of the calling object.
        """
        # navigate to the root
        root = self.getRoot()

        # convert the root object to the memo format
        rootMemo = root.__makeMemo()

        # rebuild a full copy from the memo
        rootCopy = Taxonomy.__buildTaxonomyFromMemo(rootMemo)

        # navigate to the equivalent node of the calling object and return
        return rootCopy.getDescendantByTaxId(self.taxid, resolveSynonym=True)


    def save(self, filename:str) -> None:
        """ save:
                Accepts a filename as input. Writes the full taxonomy object to
                the specified file. Does not return.
        """
        # navigate to the root
        root = self.getRoot()

        # convert the root object to the memo format
        rootMemo = root.__makeMemo()

        # write the object to file
        Taxonomy.__writeMemoToFile(rootMemo, filename)


    def load(filename:str) -> Taxonomy:
        """ load:
                Accepts a filename as input. Calls helper functions to make the
                Taxonomy object from the saved file. Returns the object at its
                root.
        """
        # read the file in as the memo format
        memo = Taxonomy.__makeMemoFromFile(filename)

        # build a Taxonomy object from the memo and return
        return Taxonomy.__buildTaxonomyFromMemo(memo)


    def _loadFromTable(table:list) -> Taxonomy:
        """ loadFromTable:
                Accepts a list in the form of a parsed Taxonomy save file as 
                input. Constructs a Taxonomy object from the table. Returns the
                Taxonomy object.
        """
        # convert the table into a memo
        memo = Taxonomy.__makeMemoFromTable(table)

        # build a Taxonomy object from the memo and return
        return Taxonomy.__buildTaxonomyFromMemo(memo)
        

    def __makeMemo(self) -> dict:
        """ makeMemo:
                Accepts a dictionary as input. Calls a helper function to con-
                vert each descendant and the calling object to the "memo" for-
                mat and add it to the memo. Returns the memo dictionary.
        """
        # initialize the memo
        memo = {}

        # get all the descendant objects
        allDesc = self.getAllDescendants(set)

        # add each descendant to the memo
        desc:Taxonomy
        for desc in allDesc:
            desc.__addTaxonToMemo(memo)

        # add the calling object to the memo and return
        self.__addTaxonToMemo(memo)
        return memo


    def __addTaxonToMemo(self, memo:dict) -> None:
        """ addTaxonToMemo:
                Accepts a dictionary as input. Converts the calling object into
                the "memo" format. Makes a new Taxonomy object that is a deep-
                copy of the calling object with one exception: the member vari-
                ables that are themselves Taxonomy objects retain their default
                initialized values. In order to save the structural data, a 3-
                item tuple is created consisting of the copy, the taxid of its
                parent, and a boolean indicating whether or not the calling ob-
                ject is the type material of that parent. The tuple is then ad-
                ded to the dictionary with its taxid as the key. Does not retu-
                rn.
        """
        # constants
        DESC = 'descendantsD'
        PRNT = 'parent'
        TYPE = 'typeMaterial'
        TAX_OBJ_KEYS = [DESC, PRNT, TYPE]
        RANK = 'rank'

        # if not in the memo then make a new thing and add it
        if self.taxid not in memo.keys():
            # make the copy (call str to ensure no tether to original)
            dup = Taxonomy(str(self.taxid), str(self.ncbiName), str(self.rank))

            # populate the copy
            dupVars = vars(dup)
            selfVars = vars(self)
            for key in dupVars:
                # deepcopy works for non-Taxonomy member variables
                if key not in TAX_OBJ_KEYS:
                    dupVars[key] = copy.deepcopy(selfVars[key])
                
            # for species, typeMaterial is a list of str; can safely deepcopy
            if dupVars[RANK] == Taxonomy.SPECIES:
                dupVars[TYPE] = copy.deepcopy(selfVars[TYPE])

            # determine the parent's taxid
            if self.parent is not None:
                parTaxId = self.parent.taxid

                # determine if the calling object is type material
                if self.parent.typeMaterial is not None:
                    isTypeMat = self.taxid == self.parent.typeMaterial.taxid
                else:
                    isTypeMat = False
            
            # set the values for the root (no parent)
            else:
                parTaxId = None
                isTypeMat = False

            # add the data to the memo
            memo[dup.taxid] = (dup, parTaxId, isTypeMat)


    def __buildTaxonomyFromMemo(memo:dict) -> Taxonomy:
        """ buildTaxonomyFromMemo:
                Accepts a dictionary in the "memo" format of a Taxonomy object
                as input. Builds a Taxonomy object from the objects contained
                within the memo. Returns the root of the newly built Taxonomy
                object.
        """
        # constants
        TAX_OBJ   = 0
        PARENT_ID = 1
        TYPE_MAT  = 2

        # for each item in the memo,
        for key in memo:
            # retrieve the object and its structural data
            taxO:Taxonomy = memo[key][TAX_OBJ]
            parTxid:str = memo[key][PARENT_ID]
            isTypeMat:bool = memo[key][TYPE_MAT]

            # retrieve its parent (if possible)
            if parTxid is not None:
                parent:Taxonomy = memo[parTxid][TAX_OBJ]

                # add the object to its parental object
                parent._importDirectDescendant(taxO)

                # if the object is type material, then update its parent
                if isTypeMat:
                    parent.typeMaterial = taxO
        
        # return the root of the newly constructed object
        root = taxO.getRoot()
        return root


    def __writeMemoToFile(memo:dict, filename:str) -> None:
        """ writeMemoToFile:
                Accepts a filename as input. Converts a Taxonomy object in memo
                format to a table and writes the table to the specified file.
                Does not return.
        """
        # constants
        DESC = 'descendantsD'
        PARENT = 'parent'
        TYPE_MAT = 'typeMaterial'
        ASS_IDS  = 'allAssIds'
        TXID_SYN = 'taxidSynonyms'
        VARS_TO_SKIP = [DESC, PARENT, TYPE_MAT, ASS_IDS]

        RANK = 'rank'
        IS_TYPE = 'isTypeMat'
        PAR_TXID = 'parentTaxid'

        VAR_TO_INDEX = {'taxid':             0,
                        'sciName':           1,
                        'ncbiName':          2,
                        'rank':              3,
                        'isExternal':        4,
                        'ncbiNameIsCorrect': 5,
                        'assemblyFtp':       6,
                        'assemblyAccn':      7,
                        'assemblyFromType':  8,
                        'assemblyStrain':    9,
                        'refSeqCategory':    10,
                        'assCoverage':       11,
                        'assIsComplete':     12,
                        'allAssIds':         13,
                        'typeMaterial':      14,
                        'taxidSynonyms':     15,
                        'parentTaxid':       16,
                        'isTypeMat':         17}
        
        EMPTY_ROW = [None] * len(VAR_TO_INDEX)

        LIST_SEP_CHAR = "|~|"
        LEN_LIST_SEP = len(LIST_SEP_CHAR)
        TAB_CHAR = "\t"
        LEN_TAB = len(TAB_CHAR)
        NEW_LINE   = "\n"

        # for each item in the memo, create a 2d array
        ### each item in the memo becomes a row (list)
        ### each row has exactly len(VAR_TO_INDEX) columns
        ### each column is a string
        outStr = ''
        for mKey in memo.keys():
            # parse the tuple
            taxon, parTxid, isTypeMat = memo[mKey]

            # for each member variable (except the variables to skip)
            row = EMPTY_ROW
            mVars = vars(taxon)
            for var in mVars:
                if var not in VARS_TO_SKIP:
                    # save the variable as a string
                    row[VAR_TO_INDEX[var]] = str(mVars[var])

            # populate the parentTaxid and isTypeMat fields
            row[VAR_TO_INDEX[PAR_TXID]] = str(parTxid)
            row[VAR_TO_INDEX[IS_TYPE]]  = str(isTypeMat)

            # handle type material and allAssIds for species
            typeIdx  = VAR_TO_INDEX[TYPE_MAT]
            assIdIdx = VAR_TO_INDEX[ASS_IDS]
            if mVars[RANK] == Taxonomy.SPECIES:
                # convert the type strain list to a string
                if mVars[TYPE_MAT] is not None and len(mVars[TYPE_MAT]) > 0:
                    typeStrainStr = ''
                    for strain in mVars[TYPE_MAT]:
                        typeStrainStr += strain + LIST_SEP_CHAR
                    typeStrainStr = typeStrainStr[:-LEN_LIST_SEP]
                else:
                    typeStrainStr = str(None)

                # convert the assembly id list to a string
                assIdStr = ''
                for assId in mVars[ASS_IDS]:
                    assIdStr += str(assId) + LIST_SEP_CHAR
                assIdStr = assIdStr[:-LEN_LIST_SEP]

                # save the strings
                row[typeIdx]  = typeStrainStr
                row[assIdIdx] = assIdStr
            
            else:
                row[typeIdx]  = None
                row[assIdIdx] = ''
            
            # convert the row to a string
            rowStr = ""
            for i in range(len(row)):
                rowStr += str(row[i]) + TAB_CHAR
            rowStr = rowStr[:-LEN_TAB]

            # add the row to the outstring
            outStr += rowStr + NEW_LINE
        
        # write the string to a file
        handle = open(filename, "w")
        handle.write(outStr)
        handle.close()


    def __makeMemoFromFile(filename:str) -> dict:
        """ makeMemoFromFile:
                Accepts a filename as input. Reads the file into memory as a 
                table with all fields as strings. Converts the table to a dict
                in the "memo" format. Returns the dictionary.
        """
        TAB_CHAR = "\t"

        # read the file into memory
        table = parseCsv(filename, TAB_CHAR)
        return Taxonomy.__makeMemoFromTable(table)


    def __makeMemoFromTable(table:list) -> dict:
        """ makeMemoFromTable:
                Accepts a list in the form of a parsed Taxonomy save file as 
                input. Decodes the strings at each column into their appropria-
                tely typed values. Uses those values to construct a dictionary
                in the "memo" format. Returns the dictionary.
        """
        # constants
        TXID = 'taxid'
        RANK = 'rank'
        NCBI = 'ncbiName'
        IS_TYPE = 'isTypeMat'
        TYPE_MAT = 'typeMaterial'
        PARENT_TXID = 'parentTaxid'
        VARS_TO_SKIP = [IS_TYPE, PARENT_TXID, RANK, NCBI, TXID]

        TAX_SYN = 'taxidSynonyms'
        ASS_COV = 'assCoverage'
        ASS_IDS = 'allAssIds'
        SCI_NAME = 'sciName'

        GREP_F = r'^(.+) \(invalid name\)$'
        GREP_R = r'"\1" (invalid name)'

        INDEX_TO_VAR = {0:  'taxid',
                        1:  'sciName',
                        2:  'ncbiName',
                        3:  'rank',
                        4:  'isExternal',
                        5:  'ncbiNameIsCorrect',
                        6:  'assemblyFtp',
                        7:  'assemblyAccn',
                        8:  'assemblyFromType',
                        9:  'assemblyStrain',
                        10: 'refSeqCategory',
                        11: 'assCoverage',
                        12: 'assIsComplete',
                        13: 'allAssIds',
                        14: 'typeMaterial',
                        15: 'taxidSynonyms',
                        16: 'parentTaxid',
                        17: 'isTypeMat'}

        SEP_CHAR_1 = "|~|"

        # initialize the memo
        memo = {}
        # for each row
        for row in table:
            # parse the row into a dictionary object keyed by the var name
            rowDict = dict()
            for idx in range(len(row)):
                # decode any lists
                if SEP_CHAR_1 in row[idx]:
                    row[idx] = row[idx].split(SEP_CHAR_1)
                
                # decode none strings
                elif row[idx] == 'None':
                    row[idx] = None
                
                # decode true strings
                elif row[idx] == 'True':
                    row[idx] = True
                
                # decode false strings
                elif row[idx] == 'False':
                    row[idx] = False
                
                # save all values to the dictionary
                rowDict[INDEX_TO_VAR[idx]] = row[idx]
            
            # correctly format sci names containing '(invalid name)'
            rowDict[SCI_NAME] = re.sub(GREP_F, GREP_R, rowDict[SCI_NAME])

            # handle empty or length==1 typeMaterial for species
            if rowDict[RANK] == 'species':
                if rowDict[TYPE_MAT] is None:
                    rowDict[TYPE_MAT] = []
                if type(rowDict[TYPE_MAT]) is str:
                    rowDict[TYPE_MAT] = [rowDict[TYPE_MAT]]

            # make ass coverage a float
            rowDict[ASS_COV] = Taxonomy.__coverageToFloat(rowDict[ASS_COV])

            # process the 'taxidSynonyms' dictionary string
            taxidSynRaw = rowDict[TAX_SYN]   # obtain the string
            taxidSynRaw = re.sub("'", '', taxidSynRaw)   # remove single-quotes
            taxidSynRaw = taxidSynRaw[1:-1]   # remove '{' and '}'

            # decode the string into the taxidSynonyms dictionary
            taxidSynonymsD = dict()
            if taxidSynRaw != '':   # already done if the dictionary was empty
                # split the string into key:val pairs
                taxidSynRaw = taxidSynRaw.split(", ")

                # split the pairs into key/val and add to a dictionary
                for pair in taxidSynRaw:
                    key, val = pair.split(": ")
                    taxidSynonymsD[key] = val
            
            # update the value in rowDict with the decoded 'taxidSynonyms' dictionary
            rowDict[TAX_SYN] = taxidSynonymsD

            # make sure allAssIds is a list and not a string
            if rowDict[ASS_IDS] == '':
                rowDict[ASS_IDS] = []
            elif type(rowDict[ASS_IDS]) is str:
                rowDict[ASS_IDS] = [rowDict[ASS_IDS]]

            # make a new Taxonomy object
            taxO = Taxonomy(rowDict[TXID], rowDict[NCBI], rowDict[RANK])

            # transfer the data stored in rowDict to the new object
            mVars = vars(taxO)
            for var in rowDict:
                if var not in VARS_TO_SKIP:
                    mVars[var] = rowDict[var]
                
            # make the "memo" format tuple and save it
            memo[rowDict[TXID]] = (taxO, rowDict[PARENT_TXID], rowDict[IS_TYPE])
        
        return memo




    ### TESTING OBJECT CONSISTENCY
    def _isConsistent(self) -> bool:
        """ isConsistent:
                Accepts no inputs. Navigates to the root of the calling object
                and then uses a recursive helper function to check for consist-
                ency between the parents and their children. Returns a boolean.
        """
        # navigate to the root
        root = self.getRoot()

        # check the root for consistency
        return root.__isNodeConsistent()


    def __isNodeConsistent(self) -> bool:
        """ isNodeConsistent:
                Accepts no inputs. Calls recursive helper functions to evaluate
                the consistency of the calling object. Returns a boolean indic-
                ating whether or not the calling object is consistent.
        """
        # check that each child lists self as parent
        if not self.__parentChecker():
            print("parent checker failed")
            return False
        
        # check that the parent lists self as child
        if not self.__childChecker():
            print("child checker failed")
            return False

        # check that the ranks are valid
        if not self.__rankChecker():
            print("rank checker failed")
            return False
        
        # check for any taxids that are shared
        if not self.__taxidChecker(set()):
            print("taxid checker failed")
            return False
        
        # check for any inconsistencies with type material
        if not self.__typeChecker():
            print("type checker failed")
            return False
        
        # check for any scientific names that are shared
        if not self.__nameChecker():
            print("name checker failed")
            return False
        
        return True

    
    def __parentChecker(self) -> bool:
        """ parentChildChecker:
                Accepts no inputs. For each of the calling object's children, 
                checks that they list the calling object as their parent.
                Returns a boolean.
        """
        # a species has no descendants, so we are done
        if self.rank == Taxonomy.SPECIES:
            return True
        
        # for each child
        children = self.getChildren(set)
        child:Taxonomy
        for child in children:
            # check that self is the child's parent
            if self != child.parent:
                return False
            
            # recurse to check the children's children
            if not child.__parentChecker():
                # False if a discrepancy occurs
                return False
        
        # if all the children were correct, then true
        return True
    
    
    def __childChecker(self) -> bool:
        """ childChecker:
                Accepts no inputs. For each of the calling objects children, 
                checks if the calling object's parent lists the calling object
                as its child. Returns a boolean.
        """
        # get all the children
        allChildren = self.getAllDescendants()

        # for each child
        for key in allChildren.keys():
            child:Taxonomy = allChildren[key]

            # check that the child's parent lists it as a child
            parent = child.parent
            if child not in parent.getChildren(set):
                return False
        
        # if all children passed, then true
        return True
        

    def __rankChecker(self) -> bool:
        """ rankChecker:
                Accepts no inputs. Recursively evaluates if the taxonomic rank
                of the calling object is exactly one rank below its parent and
                exactly one rank above its children. Returns a boolean.
        """
        # get the current rank and those that flank it
        curRank = self.rank
        parRank = curRank.getRankAbove()
        kidRank = curRank.getRankBelow()

        # make sure that the parent is the expected rank
        if self.parent is not None:
            if self.parent.rank != parRank:
                return False
        
        # make sure the children are the expected rank
        child:Taxonomy
        for child in self.getChildren(set):
            if child.rank != kidRank:
                return False
            
            # recurse on the child objects
            child.__rankChecker()
        
        # if parents and children have appropriate ranks then true
        return True


    def __taxidChecker(self, uids:set=set()) -> bool:
        """ taxidChecker:
                Accepts a set of unique taxids as input. Recursively searches
                for any instances where there are shared taxids. Returns false
                if any taxids are shared, otherwise returns true.
        """
        # add the calling object to the set
        if self.taxid not in uids:
            uids.add(str(self.taxid))
        
        # if duplicate keys exist then false
        else:
            return False
        
        # recurse on the children
        child:Taxonomy
        for child in self.getChildren(set):
            # return false if a non-unique id is found
            if not child.__taxidChecker(uids):
                print(repr(child))
                return False
        
        # if all the taxids are unique, then true
        return True


    def __typeChecker(self) -> bool:
        """ typeChecker:
                Accepts no inputs. Recursively searches the calling object for
                any instances where the type material object is not also a dir-
                ect descendant of the calling object. If such an instance is 
                found, returns false. Otherwise, returns true.
        """
        # for ranks higher than species
        if self.rank > Taxonomy.SPECIES:
            # if the type material field is populated
            if self.typeMaterial is not None:
                # then the type should be found in a set of its children
                if self.typeMaterial not in self.getChildren(set):
                    return False
                # and the type should list it as the type's parent
                if self != self.typeMaterial.parent:
                    return False
            
            # recurse on the children
            child:Taxonomy
            for child in self.getChildren(set):
                if not child.__typeChecker():
                    return False
        
        # true only if all children and self are also true
        return True
    

    def __nameChecker(self) -> bool:
        """ nameChecker:
                Accepts no inputs. Recursively searches the calling object for
                any objects with shared scientific names. If one or more shared
                names are found, returns false. Otherwise, returns true.
        """
        sharedNames = self.__findSubTaxaWithSharedSciNames()
        
        if len(sharedNames) > 0:
            return False
        
        else: return True




    ### BUILDING AND MODIFYING OBJECTS
    def _initializeDescendants(self, ncbiTaxaParsed:dict=None, \
                              maxDepth='species') -> None:
        """ initializeDescendants:
                Accepts a dictionary containing the parsed NCBI taxonomy data 
                and a TaxRank indicating the maximum recursive depth as input. 
                Recursively adds Taxonomy objects that are descendants of the 
                current object until it reaches the maximum depth. Does not re-
                turn.
        """
        # constants used for querying NCBI
        TAXONOMY_DB = 'taxonomy'
        BAD_CHAR_GREP = r'[\[\]\(\)]'  # brackets or parentheses
        SPECIES_GREP = r'^\S+ \S+$'    # 2 words; sep=' '
        HIGHER_TAX_GREP = r'^\S+$'     # 1 word; no whitespace characters

        # make sure maxDepth is a TaxRank object
        if type(maxDepth) is not TaxRank:
            maxDepth = TaxRank(maxDepth)

        # only proceed if the calling object is at higher rank than maxDepth
        if self.rank > maxDepth:
            # only query NCBI if a dictionary has not been passed in
            if ncbiTaxaParsed is None:
                # attain a list of all taxids from NCBI that are contained
                # within the taxonomic rank of the calling object 
                allTxids = self.__getChildTaxIdsFromNcbi(maxDepth)
                
                # download the summaries from NCBI (much faster than efetch!!)
                ncbiSums = ncbiSummaryFromIdList(allTxids, TAXONOMY_DB)
                
                # go through the summaries and save txid for well-named taxa
                goodTxids = list()
                for sumry in ncbiSums:
                    # extract data of interest from the summary
                    txid = sumry['Id']
                    name = sumry['ScientificName']
                    rank = TaxRank(sumry['Rank'])

                    # save entries ranks higher than species
                    if rank > Taxonomy.SPECIES:
                        # higher ranks should be exactly one word
                        if re.match(HIGHER_TAX_GREP, name):
                            goodTxids.append(txid)
                    
                    # save species-level entries
                    else:
                        # species should have binomial names
                        # speices should not contain the string 'sp.'
                        if re.match(SPECIES_GREP, name) and 'sp.' not in name:
                            goodTxids.append(txid)

                # if good taxids were found
                if len(goodTxids) > 0:
                    # get the detailed records for the good taxids
                    # (this step is very slow!)
                    ncbiTaxa = ncbiEfetchById(goodTxids, TAXONOMY_DB)
                
                # otherwise, return as there is nothing left to-do
                else:
                    return
                
                # initialize a dictionary to hold the parsed data
                ncbiTaxaParsed = dict()
                
                # parse the important information from the NCBI taxa
                taxon:dict
                for taxon in ncbiTaxa:
                    # extract data of interest from the record
                    txid = taxon['TaxId']
                    name = taxon['ScientificName']
                    rank = TaxRank(taxon['Rank'])

                    # extract the full lineage and reverse it (child -> parent)
                    lineage:Parser.ListElement = taxon['LineageEx']
                    lineage.reverse()

                    # get the taxid one valid rank above the current taxon
                    parentId = None
                    ancestor:dict
                    for ancestor in lineage:
                        # extract the ancestor's rank
                        ancRank = ancestor['Rank']

                        # if it is a valid rank and exactly one level higher
                        if TaxRank.isValidRank(ancRank):
                            if TaxRank(ancRank) - rank == 1:
                                # stop looping once the parent id is found
                                parentId = ancestor['TaxId']
                                break
                    
                    # make sure parentId gets assigned if nothing was found
                    if parentId is None:
                        parentId = taxon['ParentTaxId']

                    # remove bad characters from the names
                    name = re.sub(BAD_CHAR_GREP, '', name)

                    # create the entry for the taxon
                    entry = dict()
                    entry['name'] = name
                    entry['txid'] = txid
                    entry['rank'] = rank

                    # add the entry to the dictionary
                    # key the entry by its parental taxid
                    if parentId in ncbiTaxaParsed.keys():
                        ncbiTaxaParsed[parentId].append(entry)
                    else:
                        ncbiTaxaParsed[parentId] = [entry]
            
            # populate descendants recursively
            try:
                ncbiDesc = ncbiTaxaParsed[self.taxid]
                skip = False
            
            # problems in NCBI's database will cause some entries to be invalid
            except KeyError:
                skip = True
                self.parent.__removeDirectDescendant(self.taxid)

            if not skip:
                for descD in ncbiDesc:
                    # get data from the dictionary
                    txid = descD['txid']
                    name = descD['name']
                    rank = descD['rank']

                    # populate the direct descendants
                    if rank == self.rank.getRankBelow():
                        self.__addDescendant(txid, name, rank)

                        # for each subrank, populate with recursion
                        if rank > Taxonomy.SPECIES:
                            tmp = self.getDescendantByTaxId(txid)
                            tmp._initializeDescendants(ncbiTaxaParsed, maxDepth)
                            

    def __getChildTaxIdsFromNcbi(self, maxDepth='species') -> list:
        """ getChildTaxIdsFromNcbi:
                Accepts no inputs. Queries NCBI taxonomy for the taxonomy ids 
                of the rank immediately below (and contained by) the current 
                rank. Returns a list of NCBI taxonomy ids that match the query.
        """
        # constants used for generating the NCBI search term
        FIELD1_PREFIX = 'txid'
        FIELD1_SUFFIX = '[Organism:exp]'
        FIELD2_SUFFIX = '[Rank]'
        AND = ' AND '
        OR = ' OR '
        DATABASE = 'taxonomy'

        # ensure that maxDepth is a TaxRank object
        if type(maxDepth) is not TaxRank:
            maxDepth = TaxRank(maxDepth)

        # make the search string to find descendant taxa
        search = FIELD1_PREFIX + self.taxid + FIELD1_SUFFIX
        search += AND + '('

        # obtain the taxonomic rank directly below the calling object
        curRank = copy.deepcopy(self.rank)
        curRank.decrement()
        while(curRank >= maxDepth):
            if curRank != Taxonomy.DOMAIN:
                search += str(curRank) + FIELD2_SUFFIX
            else:
                search += "superkingdom" + FIELD2_SUFFIX

            if curRank == maxDepth:
                search += ')'
                break
            else:
                search += OR
                curRank.decrement()
        
        # return the search result
        return ncbiIdsFromSearchTerm(search, DATABASE)


    def __addDescendant(self, childTxid:str, childName:str, childRank:TaxRank) \
                        -> None:
        """ addDescendant:
                Accepts a taxid, name, and taxonomic rank for the descendant 
                taxon, as well as the taxonomic rank indicating the maximum 
                recursive depth as inputs. Constructs a new Taxonomy object 
                for the descendant with the current Taxonomy object set as its
                parent. Adds this object to the dictionary of descendant taxa 
                for the current Taxonomy object. Does not return.
        """
        # construct a new Taxonomy object and set 'self' as the parentTax
        # make sure the external status matches and possess shared dicts
        child = Taxonomy(childTxid, childName, childRank, parentTax=self, \
                         isExternal=self.isExternal)

        # add the child to the descendants dict with the child taxid as the key
        self.descendantsD[child.taxid] = child


    def __removeDirectDescendant(self, childTxid:str) -> None:
        """ removeDirectDescendant:
                Accepts an NCBI taxonomy id as input. Removes the corresponding
                descendant from the calling object. Can only remove descendant 
                taxa that are exactly one taxonomic rank below the calling 
                object. Does not return.
        """
        # constants
        ERR_MSG = " (taxid) is not a direct descendant of the calling object."

        # make sure the input is a string
        childTxid = str(childTxid)
        
        # make sure that the input is an NCBI taxonomy id and is a descendant
        if not childTxid in self.getChildren():
            raise ValueError(childTxid + ERR_MSG)
        
        # delete the matching object
        del self.descendantsD[childTxid]


    def _importDirectDescendant(self, other:Taxonomy) -> None:
        """ importDirectDescendant:
                Accepts a Taxonomy object as input. Nests the input Taxonomy 
                object within the current Taxonomy object. Modifies both the 
                input and calling objects. Does not return.
        """
        #constants
        ERR_MSG_1 = "Input is not an instance of type 'Taxonomy'."
        ERR_MSG_2 = "Input is not exactly one taxonomic rank below the current object!"

        # make sure the input is a Taxonomy object
        if type(other) is not Taxonomy:
            raise TypeError(ERR_MSG_1)

        # make sure that the input is exactly one rank below current
        if self.rank.getRankBelow() != other.rank:
            raise AttributeError(ERR_MSG_2)
    
        # save the original parent
        oldParent = other.parent
        
        # add self as the parent of other
        other.parent = self

        # add other as a descendant of self
        self.descendantsD[other.taxid] = other

        # remove other from its old parent (if necessary)
        if oldParent is not None and oldParent != self:
            oldParent.__removeDirectDescendant(other.taxid)

        # if there is typeMaterial present, then further investigation required
        if self.typeMaterial is not None:
            # if the names match ...
            if other.sciName == self.typeMaterial.sciName:
                # if the taxids match, then update the field
                if other.taxid == self.typeMaterial.taxid:
                    self.typeMaterial = other

                # if the taxids don't match
                else:
                    # if current type's ncbi name is wrong but other's is right
                    if not self.typeMaterial.ncbiNameIsCorrect and \
                                                       other.ncbiNameIsCorrect:
                        # then update the field
                        self.typeMaterial = other
        
        # resolve any duplicate names that have been introduced
        self.__resolveDuplicateNames()

    def _importExistingSubTax(self, other:Taxonomy, lpsnD:dict) -> None:
        """ importExistingSubTax:
                Accepts a Taxonomy object and the LPSN dictionary as inputs.
                Finds a suitable node to nest the input to the calling object
                and does so, modifying the Taxonomy objects in the process. Do-
                es not return.
        """
        # constants
        TEMP_TAXID = '0'
        ERR_MSG_1 = "The calling object's rank does not exceed the input's rank."
        ERR_MSG_2 = "failed to find a suitable parent"

        # other must be at a rank below self
        if self.rank <= other.rank:
            raise AttributeError(ERR_MSG_1)

        # get the root of tax1
        root1 = self.getRoot()

        # follow the other object's parent until it is within root1
        root2 = other.getRoot()
        merged = False
        while not root1.getDescendantByTaxId(root2.taxid):
            # raise an error if we're at the rank ceiling
            if root2.rank == Taxonomy.DOMAIN:
                raise RuntimeError(ERR_MSG_2)
            
            # find a new parent
            root2 = root2._findNewParent(lpsnD)

            # handle the introduction of artificial ids
            if root2.taxid == TEMP_TAXID:
                # first check if name is already present in root1
                if root2.sciName in root1:
                    newId = root1.getDescendantBySciName(root2.sciName).taxid

                else:
                    # get the next unused id from root1
                    newId = root1.__getUnusedArtificialTaxId()

                    # make sure it hasn't already been used in root2
                    while newId in root2:
                        newId = str(int(newId) - 1)

                # replace the temporary id with the next unused artificial one
                root2.taxid = newId

            # a merge is required if the ranks are equal but taxids differ
            if root1.rank == root2.rank and root1.taxid != root2.taxid:
                root1 = Taxonomy._mergeMultipleRoots([root1, root2], lpsnD)
                merged = True
                break

        # if not already merged after looping ...
        if not merged:
            # ... find a node in the calling object to be the new parent
            newPar1 = root1.getDescendantByTaxId(root2.taxid)
            
            # ... find an ancestor that can be a child of the new parent ...
            newDesc2 = other.getAncestorAtRank(newPar1.rank.getRankBelow())

            # ... newDesc2 is now a direct descendant of newPar1; import it
            newPar1._importDirectDescendant(newDesc2)


    def _mergeMultipleRoots(allRoots:list, lpsnD:dict) -> Taxonomy:
        """ mergeMultipleRoots:
                Accepts a list of Taxonomy objects at the root of their object
                and the LPSN dictionary as input. Identifies a taxonomic mrca
                ancestral to all the objects in the input list. Removes all
                taxa higher than the input objects' rank that are not ancestral
                to the input objects. Nests each object into the newly identif-
                ied taxonomic mrca. Returns the new Taxonomy object.
        """
        # constants
        TEMP_TAXID = '0'
        ERR_1 = "All Taxonomy objects in the list must be at their root"
        ERR_2 = "Cannot merge objects with shared taxids. Aborting."
        ERR_3 = "Cannot merge different domains (eg. Archaea and Bacteria)."
    
        # sort inputs so the highest ranks come first (important for looping)
        taxaL = sorted(allRoots, key=lambda x: x.rank, reverse=True)

        # loop through the list, comparing a pair of objects each iteration
        for idx in range(1,len(taxaL)):
            # get the previous index
            prv = idx - 1

            # extract the objects to compare
            tax1:Taxonomy = taxaL[prv]
            tax2:Taxonomy = taxaL[idx]

            # make sure they are at the root
            if not tax1.isRoot() or not tax2.isRoot():
                raise AttributeError(ERR_1)

            # if the taxids are identical, a safe merge cannot be performed
            if tax1.taxid == tax2.taxid:
                raise AttributeError(ERR_2)

            # make sure the objects start at the same rank
            # only works if the highest ranks are encountered first (sorted)
            while tax1.rank > tax2.rank:
                # get the new parent
                tax2 = tax2._findNewParent(lpsnD)

                # handle the introduction of artificial ids
                if tax2.taxid == TEMP_TAXID:
                    tax2.taxid = tax2.__getUnusedArtificialTaxId()

                # update the list so that they contain only roots 
                taxaL[idx] = tax2

        # get a set of all the unique root taxids
        uniqRoots = set([x.getRoot().taxid for x in taxaL])

        # if there is only one unique taxid
        if len(uniqRoots) == 1:
            # then the mrca has already been found
            found = True
        
        # otherwise, mrca has not yet been found
        else:
            found = False

        # keep track of the current rank
        curRank = copy.deepcopy(tax1.rank)

        # if the current rank is domain or mrca is already found
        # then root needs to be defined.
        if curRank == Taxonomy.DOMAIN or found:
            # get the highest ranking object from the list (first item!)
            highestRankingObject:Taxonomy = taxaL[0]

            # get the root of the highest ranking object
            root = highestRankingObject.getRoot()

            # indicate if the elements in the list are also the roots
            # this is important for how _importExistingSubTax functions
            taxaAreMrca = True
        
        # if not a domain or not already found, then the taxa will not be roots
        else:
            taxaAreMrca = False

        # initialize variables for looping
        mrca = None
        taxO:Taxonomy

        # continue to deepen the root of each object until an mrca is found
        while not found and curRank != Taxonomy.DOMAIN:
            # for each taxonomy object in the list
            for taxO in taxaL:
                # ensure we are still at the root of the object
                root = taxO.getRoot()

                # deepen one rank at a time so all roots remain equally deep
                if root.rank == curRank:
                    # get the new parent
                    root = root._findNewParent(lpsnD)

                    # handle the introduction of artificial ids
                    if root.taxid == TEMP_TAXID:
                        root.taxid = root.__getUnusedArtificialTaxId()
                
                # make sure that all roots are at the same rank
                # it is possible that the 'break' will cause some to be lower
                # this else/while block allows unprocessed taxa to catch up
                else:
                    while root.rank <= curRank:
                        root = root._findNewParent(lpsnD)

                        if root.taxid == TEMP_TAXID:
                            root.taxid = root.__getUnusedArtificialTaxId()

                # for the first object, save the root as the mrca
                if mrca is None:
                    mrca = root.taxid
                
                # if a mismatch is found, then restart the for-loop
                elif mrca != root.taxid:
                    mrca = None
                    found = False
                    break
                
                # if a match is found, keep checking the other objects
                else: found = True
            
            # increment to the next rank before looping again
            curRank.increment()

        # at this point, all roots should be congruent with one another.
        # we just need to nest the children of each root into one of the roots
        # root was defined in the loop above (or in the conditional before it)
        taxO:Taxonomy
        for taxO in taxaL:

            # handle any collisions from artificial taxids
            # get a list of all artificial ids in the root object
            artificialIds = root.__getAllArtificialTaxIds()

            if len(artificialIds) > 0:
                # save the current minimum value for artificial ids in root
                curMin = min(map(int, artificialIds))

            # replace any artificial ids in taxO with unused ones in root
            for txid in artificialIds:
                # attempt to extract a descendant from an in-use art. id
                desc = taxO.getDescendantByTaxId(txid)

                # if a descendant was found, then its taxid needs to change
                if desc:
                    # move curMin to the next available taxid
                    curMin -= 1

                    # update the descendant's taxid with the safe value
                    desc.taxid = str(curMin)

                    # add the descendant to its parent's dict under new id
                    desc.parent.descendantsD[str(curMin)] = desc

                    # remove the old id from its parent's dict
                    del desc.parent.descendantsD[txid]

            # handle when the items in the list are also roots
            if taxaAreMrca:
                # two different taxa cannot be merged
                if root.sciName != taxO.sciName:
                    raise AttributeError(ERR_3)
                
                # taxO is a domain, get the children of taxO
                children = taxO.getChildren(list)

                # nest each child beneath the root (also a domain)
                child:Taxonomy
                for child in children:
                    root._importExistingSubTax(child, lpsnD)

            # if the items in the list are not roots, then import them
            else:
                # nest taxO beneath the root
                root._importExistingSubTax(taxO, lpsnD)

        return root


    def __externalHelper(self, isExternal:bool) -> None:
        """ externalHelper:
                Accepts a boolean indicating the new value for 'isExternal'.
                Modifies the calling object's member variable, and recursively
                modifies all of the calling object's children the same way.
                Does not return.
        """
        # set the ingroup status for the calling object
        self.isExternal = isExternal

        # recurse through each rank below the calling object
        if self.rank > Taxonomy.SPECIES:
            descL = self.getChildren(list)

            desc:Taxonomy
            for desc in descL:
                desc.__externalHelper(isExternal)
    

    def __setAsExternal(self) -> None:
        """ setAsExternal:
                Accepts no inputs. Modifies the calling object's 'isExternal' 
                member variable to be 'True' and recursively modifies all of
                the calling object's children the same way. Does not return.
        """
        self.__externalHelper(isExternal=True)


    def __setAsInternal(self, modifyChildren:bool=False) -> None:
        """ setAsInternal:
                Accepts a boolean (optional) indicating whether or not to 
                recursively modify all of the calling object's children as 
                input. Modifies the calling object's 'isExternal' member 
                variable to be 'False' and, if requested, recursively modifies 
                all of the calling object's children the same way. Does not 
                return.
        """
        # constants
        NEW_VALUE = False

        # recursively modify children only if requested
        if modifyChildren:
            self.__externalHelper(isExternal=NEW_VALUE)
        
        # otherwise, change only the calling object (default)
        else:
            self.isExternal = NEW_VALUE




    ### QUERYING OBJECTS
    def getDescendantByTaxId(self, taxid, resolveSynonym:bool=False) \
                             -> Taxonomy:
        """ getDescendantByTaxId:
                Accepts an NCBI taxonomy id as input. Recursively searches the
                current Taxonomy object for a Taxonomy object whose taxonomy id
                matches the input. Returns the matching Taxonomy object if one 
                is found. Otherwise, returns False.
        """
        # convert taxid to a string (enables int as input)
        taxid = str(taxid)

        # if the taxid is a synonym, then replace it with the preferred taxid
        if resolveSynonym and taxid in self.taxidSynonyms.keys():
            taxid = self.taxidSynonyms[taxid]
        
        # check self for the taxid
        if self.taxid == taxid:
            return self

        # get all the children
        allDesc = self.getAllDescendants()

        # return the matching child, or false if none found
        if taxid in allDesc.keys():
            return allDesc[taxid]
        
        return False


    def getDescendantBySciName(self, sciName:str) -> Taxonomy:
        """ getDescendantBySciName:
                Accepts a scientific name as input. Recursively searches the
                current Taxonomy object for a Taxonomy object whose scientific
                name matches the input. Returns the matching Taxonomy object if
                one is found. Otherwise, returns False.
        """
        # check self for the sci name
        if self.sciName == sciName:
            return self

        # get all the descendants of the calling object
        allDesc = self.getAllDescendants()

        # return the descendant that has the requested name
        for key in allDesc.keys():
            desc:Taxonomy = allDesc[key]
            if desc.sciName == sciName:
                return desc
        
        # return false if no descendant was found
        return False


    def getDescendantByNcbiName(self, ncbiName:str) -> Taxonomy:
        """ getDescendantByNcbiName:
                Accepts a scientific name as input. Recursively searches the
                current Taxonomy object for a Taxonomy object whose NCBI name
                matches the input. Returns the matching Taxonomy object if one
                is found. Otherwise, returns False.
        """
        # check self for the ncbi name
        if self.ncbiName == ncbiName:
            return self

        # get all the descendants of the calling object
        allDesc = self.getAllDescendants()

        # return the descendant that has the requested name
        for key in allDesc.keys():
            child = allDesc[key]
            if child.ncbiName == ncbiName:
                return child
        
        # return false if no descendant was found
        return False


    def getLineage(self) -> list:
        """ getLineage:
                Accepts no inputs. Recursively retrieves and returns the taxo-
                nomic lineage for the current object and its parents.
        """
        # base case: no parent
        if self.parent is None:
            return []
        
        # recursive case: self.parent exists
        else:
            return self.parent.getLineage() + [self.parent]
    

    def printLineage(self) -> None:
        """ printLineage:
                Accepts no inputs. Constructs and prints the lineage of the 
                current Taxonomy object. Does not return.
        """
        # constants
        GAP = ' '*2

        # obtain the lineage
        lineage = self.getLineage()

        # add the calling object to the lineage
        lineage.append(self)

        # print initial statement
        print("Lineage for " + self.sciName + ":")

        # initialize a counter to dynamically increase the indentation
        count = 1

        # print an entry for each taxon within the lineage
        for taxon in lineage:
            # extract info from taxon
            name = taxon.sciName
            rank = str(taxon.rank)
            txid = taxon.taxid

            # determine the indentation
            indent = GAP*count

            # construct and print the string
            print(indent + name + ' (' + rank + ', ' + txid + ')')

            # increment the count
            count += 1


    def getAncestorAtRank(self, rank) -> Taxonomy:
        """ getAncestorAtRank:
                Accepts a taxonomic rank (string or TaxRank) as input. Returns 
                the Taxonomy object of the calling object's ancestor at the sp-
                ecified rank.
        """
        # make sure rank is a TaxRank object
        if type(rank) is not TaxRank:
            rank = TaxRank(rank)

        # raise an error if the rank is below self
        if rank < self.rank:
            raise AttributeError("Provided rank is lower than the calling object's rank.")
        
        # base case: rank matches self
        if rank == self.rank:
            return self
        
        # raise an error if the rank exceeds the root's rank
        if self.parent is None:
            raise AttributeError("Provided rank exceeds the root object's rank.")
        
        # recursive case: rank exceeds self
        else:
            return self.parent.getAncestorAtRank(rank)
    

    def getRoot(self) -> Taxonomy:
        """ getRoot:
                Accepts no inputs. Recursively follows self.parent until the 
                calling object's parent is None, then returns that object.
        """
        # base case: calling object is the root
        if self.parent is None:
            return self
        
        # recursive case: calling object has a parent
        else:
            return self.parent.getRoot()


    def numChildren(self) -> int:
        """ numChildren:
                Accepts no inputs. Returns the number of immediate taxonomic 
                descendants (eg. the number of genera in a family).
        """
        return len(self.descendantsD)


    def getChildren(self, retType:type=dict):
        """ getChildren:
                Accepts either dict (type) or list (type) as inputs (others may
                work but have not been tested; default = dict). Returns an 
                object of the input type containing all of the descendant taxa
                as Taxonomy objects.
        """
        if retType == dict:
            return self.descendantsD
        else:
            return retType(self.descendantsD.values())
    

    def descendantsKeyedBySciName(self) -> dict:
        """ descendantsKeyedByName:
                Accepts no inputs. Returns a dictionary keyed by the sciName of
                the descendants of the calling object with the taxids as the 
                values.
        """
        # get a set of the child objects
        childrenS = self.getChildren(set)

        # make and return a dict keyed by child names
        namesD = {}
        for child in childrenS:
            namesD[child.sciName] = child.taxid
        return namesD
    

    def getDescendantsAtRank(self, rank, rettype:type=dict):
        """ getDescendantsAtRank:
                Accepts a taxonomic rank (string or TaxRank) and a return type
                as input. Obtains a collection of all descedants of the calling
                object at the input rank. Converts the collection to the input
                return type and returns the collection.
        """
        # ensure rank is a TaxRank object
        if type(rank) is not TaxRank:
            rank = TaxRank(rank)
        
        # raise an error if rank is not less than current
        if rank >= self.rank:
            raise AttributeError(str(rank) + " is not lower than the current object (" + str(self.rank) + ')')
        
        # get all the descendants
        allDesc = self.getAllDescendants()

        # make a new dictionary to hold the descendants at the desired rank
        descAtRank = dict()
        for key in allDesc:
            desc = allDesc[key]
            if desc.rank == rank:
                descAtRank[key] = desc
        
        if rettype is dict:
            return descAtRank
        else:
            return rettype(descAtRank.values())


    def numDescendantsAtRank(self, rank) -> int:
        """ numDescendantsAtRank:
                Accepts a taxonomic rank (string or TaxRank) as input. Returns 
                the total number of descendant taxa of the calling object that
                are at the specified rank.
        """
        return len(self.getDescendantsAtRank(rank))


    def getAllDescendants(self, rettype:type=dict):
        """ getAllDescendants:
                Accepts a type as input. Recursively generates a dictionary
                that contains all of the taxonomic objects that are descendants
                of the calling object. Converts the dictionary to the desired
                return type and returns that object.
                
                Rettype can be dict (default), set, or list.
        """
        # initialize output
        allDescendantsD = dict()
        
        # get the direct descendants
        children = self.getChildren()

        # add children to the descendants list
        allDescendantsD.update(children)

        # recurse on each child
        for key in children.keys():
            child:Taxonomy = children[key]
            allDescendantsD.update(child.getAllDescendants())

        return Taxonomy.__applyReturnType(allDescendantsD, rettype)
    

    def __applyReturnType(collection:dict, rettype:type):
        """ applyReturnType:
                Accepts a dictionary and a return type as inputs. Converts the
                dictionary to the specified return type and returns the new co-
                llection.
        """
        # if the return type is dict, then nothing to do.
        if rettype == dict:
            return collection
        
        # otherwise, convert the dict to the desired return type
        else:
            return rettype(collection.values())


    def getSiblings(self, rettype:type=list):
        """ getSiblings:
                Accepts a return type as input. Returns a collection of Taxono-
                my objects that are siblings of the calling object.
        """
        # if there is no parent, then the object has no siblings
        if self.parent is not None:
            # add all the parent's children (except self) to the collection
            allChildren = self.parent.getChildren(set)

            # if self is the only child, then no siblings
            if len(allChildren) > 1:
                # make a dictionary of the siblings keyed by taxid
                siblings = dict()
                child:Taxonomy
                for child in allChildren:
                    # make sure not to include the calling object as a sibling
                    if child != self:
                        siblings[child.taxid] = child
        
            # the calling object has no siblings if it's parent has ≤1 children
            else:
                siblings = dict()

        # the calling object has no siblings if it does not have a parent
        else:
            siblings = dict()

        # convert the collection to the desired type and return
        return Taxonomy.__applyReturnType(siblings, rettype)


    def numSiblings(self) -> int:
        """ numSiblings:
                Accepts no inputs. Returns the number of siblings of the calli-
                ng object.
        """
        return len(self.getSiblings(set))


    def getExternalChildren(self, rettype:type=dict):
        """ getExternalChildren:
                Accepts a return type as an input. Constructs a collection of
                the children of the calling object which are already marked as
                external. Converts the collection to the input return type and
                returns the collection.
        """
        # get all the descendants as a list
        allKids = self.getChildren()

        # populate a dictionary with only the external kids
        extKids = dict()
        for key in allKids:
            kid = allKids[key]
            if kid.isExternal:
                extKids[key] = kid

        # convert the dictionary to the desired type and return
        if rettype is dict:
            return extKids
        else:
            return rettype(extKids.values())


    def getExternalDescendantsAtRank(self, rank, rettype:type=list) -> list:
        """ getExternalDescendantsAtRank:
                Accepts a TaxRank object and a return type as input. Constructs
                a collection of descendants that at the provided rank that are 
                already marked as external. Returns the collection as the requ-
                ested return type.
        """
        # make sure the input is a TaxRank
        if type(rank) is not TaxRank:
            rank = TaxRank(rank)

        # get all descendants at the specified rank
        allKidsD = self.getDescendantsAtRank(rank)

        # populate and return a list with the internal kids
        intKidsD = dict()
        for taxid in allKidsD.keys():
            kid:Taxonomy = allKidsD[taxid]
            if not kid.isExternal:
                intKidsD[taxid] = kid

        return Taxonomy.__applyReturnType(intKidsD, rettype)


    def getInternalChildren(self, rettype:type=dict):
        """ getInternalChildren:
                Accepts a return type as an input. Constructs a collection of
                the children of the calling object which are already marked as
                internal. Converts the collection to the input return type and
                returns the collection.
        """
        # get all the descendants as a list
        allKids = self.getChildren()

        # populate a new dict with only the internal kids
        intKids = dict()
        for key in allKids:
            kid = allKids[key]
            if not kid.isExternal:
                intKids[key] = kid

        # convert collection to the desired return type
        if rettype is dict:
            return intKids
        else:
            return rettype(intKids.values())


    def getInternalDescendantsAtRank(self, rank, rettype:type=list):
        """ getInternalDescendantsAtRank:
                Accepts a TaxRank object (or equivalent string) and a return
                type as input. Constructs a collection of descendants that at
                the provided rank that are already marked as internal. Returns
                the collection as the requested return type.
        """
        # make sure the input is a TaxRank
        if type(rank) is not TaxRank:
            rank = TaxRank(rank)

        # get all descendants at the specified rank
        allKidsD = self.getDescendantsAtRank(rank)

        # populate and return a list with the internal kids
        intKidsD = dict()
        for taxid in allKidsD.keys():
            kid:Taxonomy = allKidsD[taxid]
            if not kid.isExternal:
                intKidsD[taxid] = kid

        return Taxonomy.__applyReturnType(intKidsD, rettype)


    def numSavedAssemblies(self) -> int:
        """ numSavedAssemblies:
                Accepts no inputs. Returns the number of species with an assem-
                bly from NCBI.
        """
        # if the calling object is a species, return 1 or 0
        if self.rank == Taxonomy.SPECIES:
            if self.assemblyFtp is not None:
                return 1
            else:
                return 0
                
        # otherwise get all the species
        kids = self.getDescendantsAtRank('species')

        # count and return the number of species with assemblies
        outSum:int = 0
        for key in kids.keys():
            kid:Taxonomy = kids[key]

            if kid.assemblyFtp is not None:
                outSum += 1
    
        return outSum


    def getTypeAtRank(self, rank) -> Taxonomy:
        """ getTypeAtRank:
                Accepts a taxonomic rank (or coercible str) as input. Recur-
                sively follows self.typeMaterial until the self.rank is equal
                to the input, then returns that object.
        """
        # make sure input is a rank (or can be coerced to one)
        if type(rank) is not TaxRank:
            rank = TaxRank(rank)
        
        # make sure rank is not greater than self
        if rank > self.rank:
            raise ValueError("Rank exceeds that of the calling object.")

        # follow the type material until the desired rank is reached
        if self.rank > rank:
            if self.typeMaterial is not None:
                return self.typeMaterial.getTypeAtRank(rank)
            else:
                return None
        else:
            return self


    def getAllSpeciesWithAssemblies(self, rettype:type=list):
        """ getAllSpeciesWithAssemblies:
                Accepts a return type as input. Creates a collection of Taxono-
                my objects that are descendants of the calling object, and pos-
                sess an assembly. Returns the collection as an object of the
                requested type.
        """
        # get a list of all descendant species
        descSpecies:list = self.getDescendantsAtRank(Taxonomy.SPECIES, list)

        # get a dictionary for only those species that possess assemblies
        outD = Taxonomy.__filterListForSpeciesWithAssemblies(descSpecies)

        # convert to the requested type before returning
        return Taxonomy.__applyReturnType(outD, rettype)
    

    def __filterListForSpeciesWithAssemblies(speciesL:list) -> dict:
        """ filterListForSpeciesWithAssemblies:
                Accepts a list of Taxonomy objects as input. Makes a dictionary
                keyed by taxids whose values are the corresponding Taxonomy ob-
                jects. Only saves those taxa that possess an assembly. Returns
                the dictionary.
        """
        # for each species, add those with an assembly to the output
        outD = dict()
        spe:Taxonomy
        for spe in speciesL:
            if spe.assemblyFtp is not None:
                outD[spe.taxid] = spe
        
        return outD


    def getAllIngroupCandidateSpecies(self, rettype:type=list):
        """ getIngroupCandidates:
                Accepts a return type as input. Creates a collection of Taxono-
                my objects that are descendants of the calling object, marked
                as internal, and possess an assembly. Returns the collection as
                an object of the requested type.
        """
        # get all internal species
        intSpecies:list = self.getInternalDescendantsAtRank(Taxonomy.SPECIES)
        
        # get a dictionary for only those species that possess assemblies
        outD = Taxonomy.__filterListForSpeciesWithAssemblies(intSpecies)

        # convert to the requested type before returning
        return Taxonomy.__applyReturnType(outD, rettype)
    

    def isRoot(self) -> bool:
        """ isRoot:
                Accepts no inputs. Returns a boolean indicating if the calling
                object is also the root of the taxonomy.
        """
        return self.rank == self.getRoot().rank

    
    def getTaxaWithIncorrectNcbiNames(self, rettype:type=dict):
        """ getTaxaWithIncorrectNcbiNames:
                Accepts a return type as input. Creates a collection of Taxono-
                my objects whose NCBI names are flagged as being incorrect. Re-
                turns the collection as an object of the requested type.
        """
        # initialize output
        allConflictingTaxaD = dict()

        # navigate to the root
        taxO = self.getRoot()

        # check self
        if not taxO.ncbiNameIsCorrect:
            allConflictingTaxaD.update({taxO.taxid: taxO})
        
        # check all descendants of the root
        desc:Taxonomy
        for desc in taxO.getAllDescendants(list):
            if not desc.ncbiNameIsCorrect:
                allConflictingTaxaD.update({desc.taxid: desc})
        
        # apply the return type and return
        return Taxonomy.__applyReturnType(allConflictingTaxaD, rettype)




    ### REMODELING OBJECTS ACCORDING TO LPSN DATA
    def _addLpsnData(self, lpsnD:dict, removeEmptyTaxa:bool=True) -> None:
        """ addLpsnData:
                Accepts the LPSN dictionary and a boolean indicating whether or
                not to remove taxa higher than species with no descendants as
                inputs. Calls helper functions to restructure the calling obje-
                ct so that it is congruent with the LPSN taxonomy. Does not re-
                turn.
        """
        # update the names before doing anything else
        self.__renameSynonyms(lpsnD)

        # resolve any duplicate names introduced by the renaming
        self.__resolveDuplicateNames()

        # get the original number of siblings of the calling object
        originalNumSiblings = self.numSiblings()

        # update the taxonomic structure to match LPSN
        self.__updateStructureByLpsn(lpsnD, removeEmptyTaxa)


        # if new siblings were introduced, they need to be processed
        if self.numSiblings() > originalNumSiblings:
            self.__processSiblings(lpsnD)

        # finally, import the type material from LPSN
        self.__importTypeMaterial(lpsnD)


    def __renameSynonyms(self, lpsnD:dict) -> None:
        """ renameSynonyms:
                Accepts an LPSN dictionary as input. Uses the dictionary to re-
                cursively replace the synonymous names with the preferred ones.
                Does not return.
        """
        # constants
        SUFFIX = '_syn'
        INVALID = '" (invalid name)'

        # there is no lpsn dict for bacteria or archaea (no type; no parent)
        if self.rank < Taxonomy.DOMAIN:
            # determine which lpsnD dict to use
            rnkKey = str(self.rank)
            synKey = rnkKey + SUFFIX

            # rename self if necessary
            if self.sciName in lpsnD[synKey].keys():
                self.sciName = lpsnD[synKey][self.sciName]
                self.ncbiNameIsCorrect = False
            
            # mark taxa missing from LPSN as invalid
            if self.sciName not in lpsnD[rnkKey].keys():
                if INVALID not in self.sciName:
                    self.sciName = '"' + self.sciName + INVALID
            
        # for ranks higher than species
        if self.rank > Taxonomy.SPECIES:
            # recurse on the children
            child:Taxonomy
            for child in self.getChildren(list):
                child.__renameSynonyms(lpsnD)


    def __resolveDuplicateNames(self) -> None:
        """ resolveDuplicateNames:
                Accepts no inputs. Iterates through the taxa with shared names.
                For each, identifies the preferred taxon and migrates all of 
                the synonyms' children underneath it. Deletes the synonym Taxo-
                nomy objects once their children have been depopulated. Does 
                not return.
        """
        # get the duplicate named taxa
        problemTax = self.__findSubTaxaWithSharedSciNames()

        # sort the keys of problemTax from lowest taxonomic rank to highest
        sortMethod = lambda x: problemTax[x][0].rank
        sortedKeys = sorted(problemTax.keys(), key=sortMethod)

        # for each scientific name needing to be resolved
        for dupName in sortedKeys:
            # extract the list of Taxonomy objects that are colliding
            probList:list = problemTax[dupName]

            # track the ncbi names while looping (prevents looping again later)
            ncbiNamesL:list = list()

            # find the preferred node, if possible
            probTax:Taxonomy
            preferredNode = None
            for probTax in probList:
                # if the ncbiName is correct, then we found the preferred node
                if probTax.ncbiNameIsCorrect:
                    preferredNode = probTax
                    break
            
                # save the ncbiName in case preferred node cannot be found
                ncbiNamesL.append(probTax.ncbiName)

            # if none of the names were preferred, then find one another way
            if preferredNode is None:
                preferredNode = self.__findPreferredNode(probList, ncbiNamesL \
                                                                     , dupName)
                
            # remove the preferred node from the list of problem taxa
            if preferredNode in probList:
                probList.remove(preferredNode)
            
            # for each of the remaining problem taxa
            synonomousTax:Taxonomy
            for synonomousTax in probList:
                # get the parent
                parent = synonomousTax.parent
                # synonomousTax = parent.getDescendantByTaxId(synonomousTax.taxid)

                # get a list of the kids that may become orphaned
                kidsToMove = synonomousTax.getChildren(list)

                # move the children from the synonym to the preferred parent
                kid:Taxonomy
                for kid in kidsToMove:
                    # a collision occurs if:
                    #    a preferred node has a descendant with a sciName that
                    #    matches the sciName of the current kid.
                    # this is rare, but occurs in tests 12, 17, and 18
                    collision = preferredNode.getDescendantBySciName(kid.sciName)

                    # collision will be false if the look-up failed
                    # if a collision occured, then a resolution is required
                    if collision:
                        preferredNode.__resolveDuplicateNames()

                    # only import the kid if there won't be a collision
                    else:
                        preferredNode._importDirectDescendant(kid)
 
                # add the bad taxid to the synonym dictionary
                preferredNode.taxidSynonyms[synonomousTax.taxid] = \
                                                            preferredNode.taxid

                # finally, remove the bad taxon from the object
                parent.__removeDirectDescendant(synonomousTax.taxid)


    def __findSubTaxaWithSharedSciNames(self) -> dict:
        """ findSubTaxaWithSharedSciNames:
                Accepts no inputs. Uses a helper function to identify any Taxo-
                nomy objects that share the same scientific name. Returns a di-
                ctionary of the identified taxa.
        """
        # initialize a dictionary to store the duplicate taxa
        dupTaxa = dict()
        uniqueNames = self.__getUniqueNamesDict(dict())

        # for each unique name
        for key in uniqueNames:
            # get a list of the taxa
            taxaList = uniqueNames[key]
            
            # save the lists where duplicate taxa exists
            if len(taxaList) > 1:
                dupTaxa[key] = taxaList

    
        return dupTaxa


    def __getUniqueNamesDict(self, uniqueNames:dict=dict()) -> dict:
        """ getUniqueNamesDic:
                Accepts a dictionary as input. Recursively populates the input
                dictionary with a list of all the scientific names within the
                Taxonomy object that are unique. Returns a dictionary keyed by
                the unique scientific names whose values are a list of Taxonomy
                objects with that share that name.
        """
        # save the name if it hasn't been seen before
        if self.sciName not in uniqueNames.keys():
            uniqueNames[self.sciName] = [self]
        
        # otherwise, if the name has been seen before
        else:
            # initialize a boolean indicating whether to include self
            include = True

            # determine if the calling object should be added to the list
            taxon:Taxonomy
            for taxon in uniqueNames[self.sciName]:
                # if the taxid matches, then this name isn't being shared
                if self.taxid == taxon.taxid:
                    include = False
                    break

                # if the ranks don't match, then this name isn't being shared
                # (shared names across ranks happens between classes and phyla)
                if self.rank != taxon.rank:
                    include = False
                    break
            
            # if the object should be included, then add it to the list
            if include:
                uniqueNames[self.sciName].append(self)

        # recurse on each child then return
        child:Taxonomy
        for child in self.getChildren(set):
            child.__getUniqueNamesDict(uniqueNames)
        
        return uniqueNames


    def __findPreferredNode(self, collidingTaxa:list, ncbiNames:list, \
                                                      sciName:str) -> Taxonomy:
        """ findPreferredNode:
                Accepts a list of Taxonomy objects at the same rank whose scie-
                ntific names are colliding, a list of NCBI names whose indices
                correspond to those in the first list, and the scientific name
                responsible for the collision as inputs. Attempts to find which
                taxid NCBI prefers for the scientific name. If this fails, then
                selects the NCBI name that is more similar to the scientific
                name as determined via the Ratcliff/Obershelp method. Returns
                the newly found preferred node.
        """
        # constants
        AND = ' AND '
        RANK_SUFFIX = '[Rank]'

        # determine the rank for the collision; will be used in the ncbi search
        firstTaxon:Taxonomy = collidingTaxa[0]
        rankStr = str(firstTaxon.rank)

        # search ncbi for the scientific name at the desired rank
        searchStr = sciName + AND + rankStr + RANK_SUFFIX
        goodTaxid:list = ncbiIdsFromSearchTerm(searchStr, 'taxonomy')

        # initialize preferred Node
        preferredNode = None

        # ideally, only one taxid was found. Use it to find the preferred node
        if len(goodTaxid) == 1:
            goodTaxid = goodTaxid.pop()
            preferredNode = self.getDescendantByTaxId(goodTaxid)

        # if more than one taxid was found, then raise an error (untested)
        elif len(goodTaxid) > 1:
            raise RuntimeError('untested scenario')

        # If preferred node still not found, then ncbiNamesL indices will match
        # those of probList. Pick the most similar string (arbitrary).
        else:
            # initialize variables for tracking the highest scoring ncbi name
            maxScore = -float('Inf')
            maxIndex = 0

            # go through each of the ncbiNames
            for idx in range(len(ncbiNames)):
                # compare the ncbiName to the sciName
                curName = ncbiNames[idx]
                score = textdistance.ratcliff_obershelp(curName, sciName)

                # update variables if a better score is observed
                if score > maxScore:
                    maxScore = score
                    maxIndex = idx
            
            # set the preferred node to the highest scoring ncbi name
            preferredNode:Taxonomy = collidingTaxa[maxIndex]

        # raise an exception if this function failed to find a preferred node
        if preferredNode is None or preferredNode is False:
            raise RuntimeError('Failed to find a preferred node.')

        return preferredNode


    def __processSiblings(self, lpsnD:dict) -> None:
        """ processSiblings:
                Accepts the LPSN dictionary as input. Ensures that each sibling
                of the calling object is nested beneath the preferred parent
                taxon. Modifies the calling object. Does not return.
        """
        # constants
        TEMP_TAXID = '0'
        INVALID_STR = '" (invalid name)'

        # get a list of all of the siblings of the calling object
        allSiblings = self.getSiblings(list)

        # get the LPSN dictionary for the sibling's rank
        rankD:dict = lpsnD[str(self.rank)]

        # for each sibling
        sibling:Taxonomy
        for sibling in allSiblings:
            # get the current parent
            parent = sibling.parent

            # if possible, retrieve the preferred parent name
            if sibling.sciName in rankD.keys():
                preferredParentName = rankD[sibling.sciName]['parent']
                
                # compare to the current parent name to the preferred name
                # find the appropriate parent if a mismatch occurs
                if parent.sciName != preferredParentName:
                    # find a new parent for the sibling
                    siblingParent = sibling._findNewParent(lpsnD)

                    # handle the introduction of any artificial taxids
                    if siblingParent.taxid == TEMP_TAXID:
                        siblingParent.taxid = \
                                     siblingParent.__getUnusedArtificialTaxId()

                    # merge the sibling with the calling object
                    # this modifies both parent and siblingParent
                    Taxonomy._mergeMultipleRoots([parent,siblingParent], lpsnD)
            
            # if the sibling's parent's name is absent from LPSN ...
            # ... keep the same parent, because a new one cannot be located
            else:
                # ... set 'ncbiNameIsCorrect' to false
                sibling.ncbiNameIsCorrect = False

                # ... mark the name as invalid (unless it is already marked)
                if INVALID_STR not in sibling.sciName:
                    sibling.sciName = '"' + sibling.sciName + INVALID_STR


    def __updateStructureByLpsn(self, lpsnD:dict, removeEmpties:bool) -> None:
        """ updateStructureByLpsn:
                Accepts the LPSN dictionary and a boolean indicating whether or
                not to remove taxa higher than species with no descendants as
                inputs. Starts at the species level and works up. Updates the
                taxonomic position of each descendant at the current rank, and
                then does the same for the children at the next highest rank.
                Continues until all descendant taxa have been appropiately res-
                tructured. Does not return.
        """
        # start at the species level
        curRank = TaxRank('species')

        # work from species up
        while curRank < self.rank:
            # get a list of all the children at the current rank
            children = self.getDescendantsAtRank(curRank, list)

            # for each child
            child:Taxonomy
            for child in children:
                # update the tax structure to conform to LPSN
                child.__updateStructureHelper(lpsnD, removeEmpties)

            # increment the rank
            curRank.increment()
  
    
    def __updateStructureHelper(self, lpsnD:dict, removeEmpties:bool) -> None:
        """ updateStructureHelper:
                Accepts the LPSN dictionary and a boolean indicating whether or
                not to remove taxa higher than species with no descendants as
                inputs. Checks if the taxonomic structure of the calling object
                conforms with the LPSN taxonomy and modifies it if necessary.
                Depending on the conditions, the reconciliation may vary. May
                call a recursive helper function  to assist with the restructu-
                ring. Does not return.
        """
        # constants
        PARENT  = 'parent'
        TEMP_TAXID = '0'
        SEP     = "\n\t"
        ERR_MSG = "attempting to merge inequivalent objects" + \
                  " that share keys! (DANGEROUS COLLISION!)" + SEP
        
        # get lpsnD key
        lpsnKey = str(self.rank)

        # if the object has no descendants (and is not a species) ...
        if self.rank > Taxonomy.SPECIES and len(self) == 1:
            # ... remove the calling object from its parent (if requested)
            if removeEmpties:
                self.parent.__removeDirectDescendant(self.taxid)

        # otherwise, retrieve the parent name from lpsn
        elif self.sciName in lpsnD[lpsnKey].keys():
            # get the current parent and its name
            curParent = self.parent
            curParentName = curParent.sciName

            # get the lpsn parent name
            lpsnParentName = lpsnD[lpsnKey][self.sciName][PARENT]

            # if the parent name is incorrect
            if lpsnParentName != curParentName:
                # retrieve the Taxonomy object for the new parent
                newParent:Taxonomy = self._findNewParent(lpsnD, lpsnParentName)

                # protect from failure due to NoneType parent
                if newParent is not None:
                    # findNewParent can make nodes w/multiple taxonomic ranks
                    # this ensures that we reference the immediate parent
                    newParent = newParent.getDescendantBySciName(lpsnParentName)

                    # if getDescendantBySciName fails, newParent will be False
                    if newParent:
                        # if an artificial taxid has been created, then replace it
                        if newParent.taxid == TEMP_TAXID:
                            newParent.taxid = self.__getUnusedArtificialTaxId()
                        
                        # nest self under the new parent
                        newParent._importDirectDescendant(self)

                        # get the roots of the two objects
                        newRoot = newParent.getRoot()
                        curRoot = curParent.getRoot()

                        # extract booleans used in comparisons
                        sameTxid = curRoot.taxid == newRoot.taxid
                        oneChild = len(newRoot.getChildren(list)) == 1
                        inequivalentObjects = curRoot != newRoot
                        sameRank = curRoot.rank == newRoot.rank
                        newRootBelowCurRoot = newRoot.rank < curRoot.rank

                        # newRoot might be functionally identical to curRoot
                        if sameTxid and oneChild:
                            # get the correct node to import
                            newParent = self.getAncestorAtRank(curRoot.rank - 1)

                            # import the new parent beneath the current root
                            curRoot._importDirectDescendant(newParent)

                        # if the roots are siblings
                        elif inequivalentObjects and sameRank:
                            # abort if inequivalent objects share a name
                            if curRoot.sciName == newRoot.sciName:
                                # bad LPSN data can produce this error!!
                                raise RuntimeError(ERR_MSG + str(curRoot) + \
                                                    SEP + str(newRoot))

                            # find the parent of the curRoot
                            ### findNewParent will nest curRoot automatically
                            ### findNewParent will handle lpsn reconciliation
                            grandParent = curRoot._findNewParent(lpsnD)

                            # handle the introduction of aritificial taxids
                            if grandParent.taxid == TEMP_TAXID:
                                grandParent.taxid = grandParent.__getUnusedArtificialTaxId()

                            # nest the newRoot inside the grandParent
                            grandParent._importDirectDescendant(newRoot)
                
                        # otherwise, if the new root is not currently nested, but it could be
                        elif inequivalentObjects and newRootBelowCurRoot and not curRoot.getDescendantByTaxId(newRoot.taxid):
                            futureParent = curRoot.getDescendantBySciName(newRoot.sciName)
                            # find a new parent until it can be nested within the current root
                            while not futureParent and newRoot.rank < curRoot.rank:
                                newRoot = newRoot._findNewParent(lpsnD)

                                # handle the introduction of aritificial taxids
                                if newRoot.taxid == TEMP_TAXID:
                                    newRoot.taxid = newRoot.__getUnusedArtificialTaxId()

                                futureParent = curRoot.getDescendantBySciName(newRoot.sciName)
                            
                            # if a suitable parent was found, then migrate the object
                            if futureParent:
                                curParent = curRoot.getDescendantBySciName(newRoot.sciName)
                                curParent._importDirectDescendant(newParent.getAncestorAtRank(curParent.rank.getRankBelow()))
                            
                            # if a suitable parent was not found, then attempt to merge roots
                            else:
                                curRoot = Taxonomy._mergeMultipleRoots([curRoot, newRoot], lpsnD)


    def _findNewParent(self, lpsnD:dict, newParName:str=None) -> Taxonomy:
        """ findNewParent:
                Accepts an LPSN dictionary and a string for the new parent's 
                name as inputs. Uses recursion to find the best parent to use 
                for the calling object. Attempts the follow this fundamental
                strategy:
                    * If a suitable parent already exists somewhere below or 
                      above the calling object, then use that node as the 
                      parent and nest the calling object underneath.
                    
                    * Otherwise, make a suitable parent object and nest the
                      calling object underneath.

                    * Recursion may occur if multigenenerational parent objects
                      are needed to map the LPSN taxonomy onto the existing 
                      NCBI taxonomy.
                
                After creating a suitable parent that can be successfully inte-
                grated into the calling object's structure, returns a handle to
                the newly created parent Taxonomy object.
        """
        # constants
        PARENT = 'parent'
        SUPERKINGDOM = 'superkingdom'
        DATABASE = 'taxonomy'
        AND = ' AND '
        QUERY_SUFFIX = "[rank]"
        TEMP_TAXID = '0'
        PLACEHOLDER_STR = "; not assigned to "

        # get the root
        root = self.getRoot()

        # get the parent rank
        parentRank = self.rank.getRankAbove()

        # determine the new parent name (if not provided)
        if newParName is None:
            # attempt to look up lpsn parent name
            lpsnKey = str(self.rank)
            if self.sciName in lpsnD[lpsnKey]:
                newParName = lpsnD[lpsnKey][self.sciName][PARENT]
        
        # attempt to use an existing node
        newParent = root.getDescendantBySciName(newParName)

        # make sure the new parent is at the correct rank
        # (eg. 'Actinobacteria' is both a class and a phylum)
        if newParent:
            if newParent.rank != parentRank:
                newParent = False

        # if existing parent still not found
        if not newParent:
            # make sure not to look up LPSN placeholder names
            if PLACEHOLDER_STR not in newParName:
                # rename 'domain' to facilitate communication with ncbi
                if parentRank == Taxonomy.DOMAIN:
                    rankStr = SUPERKINGDOM
                else:
                    rankStr = str(parentRank)

                # attempt to look up parent on ncbi taxonomy
                query = newParName + AND + rankStr + QUERY_SUFFIX
                
                # get a list of taxids, or make an empty one if the query fails
                try:
                    taxid = ncbiIdsFromSearchTerm(query, DATABASE)
                except:
                    taxid = []
            else:
                taxid = []

            # if a taxid was found
            if len(taxid) == 1:
                # make a new node
                taxid = taxid.pop()
                newParent = Taxonomy(taxid, newParName, parentRank)
                
                # if the new parent is a descendant (indirect ok) of the root
                if newParent.rank < root.rank:
                    # use recursion to find the correct parent
                    newParent = newParent._findNewParent(lpsnD)

                    # handle the introduction of artificial taxids
                    if newParent.taxid == TEMP_TAXID:
                        newParent.taxid = \
                                         newParent.__getUnusedArtificialTaxId()
                
                # if the new parent is ancestral to the current root
                elif newParent.rank > root.rank:
                    # nest the root under the new parent
                    newParent._importDirectDescendant(root)
                
                # if the same rank, but different names, then migrate self
                elif newParent.sciName != root.sciName:
                    newParent._importDirectDescendant(self)
            
            # it's possible the name wasn't found because ncbi uses a synonym
            # it is not safe to assume that the parent is also a synonym.
            ## eg. a species was moved from a valid genus to a new genus.
            # because of this, an arificial NCBI taxid must be created
            ## findNewParent will always use 0 and expect the calling funxn
            ## to handle the assignment of an appropriate taxid
            ## all "real" NCBI taxids are integers >= 1

            # force a new parent into existence (REQUIRES AN ARTIFICIAL TAXID!)
            else:
                newParent = Taxonomy(TEMP_TAXID, newParName, \
                                     self.rank.getRankAbove())
                newParent._importDirectDescendant(self)
        
        # return the newParent
        return newParent


    def __importTypeMaterial(self, lpsnD:dict) -> None:
        """ importTypeMaterial:
                Accepts an LPSN dictionary as input. Migrates to the root and
                then calls a recursive helper function to import the type mat-
                erial data into the calling object. Does not return.
        """
        # pass the root to the helper
        root = self.getRoot()

        # the helper will recurse through the entire data structure
        root.__importTypeHelper(lpsnD)


    def __importTypeHelper(self, lpsnD:dict) -> None:
        """ importTypeHelper:
                Accepts an LPSN dictionary as input. Recursively populates the
                typeMaterial member variable for self and all its descendants.
                Does not return.
        """
        # constants
        TYPE = 'type'
        NCBI_DB = 'taxonomy'

        # no type material (or lpsnD dict) for domains
        if self.rank == Taxonomy.DOMAIN:
            # still need to recurse on their children
            for child in self.getChildren(set):
                child.__importTypeHelper(lpsnD)
            return
        
        # get the rank dictionary
        rankD:dict = lpsnD[str(self.rank)]
        
        # get the entry from the type field (if possible)
        if self.sciName in rankD.keys():
            lpsnType = rankD[self.sciName][TYPE]
        
        # if not possible, then set it to None
        else:
            lpsnType = None

        # save the list of strains if the calling object is a species
        if self.rank == Taxonomy.SPECIES:
            # empty list instead of None for species
            if lpsnType is None:
                lpsnType = []
            
            # update the member variable
            self.typeMaterial = lpsnType
        
        # for ranks between domain and species
        else:
            # search for the type material in an existing child first
            child:Taxonomy
            for child in self.getChildren(list):
                if child.sciName == lpsnType:
                    self.typeMaterial = child
                
                # recurse on all the children
                child.__importTypeHelper(lpsnD)

            # if the type material field is still empty
            if self.typeMaterial is None and lpsnType is not None:
                # attempt to make a node for the type material
                newTxid:list = ncbiIdsFromSearchTerm(lpsnType, NCBI_DB)

                # only do something if exactly one id was retrieved
                if len(newTxid) == 1:
                    newTxid = newTxid.pop()

                    # make sure no collisions will occur
                    if self.getRoot().getDescendantByTaxId(newTxid, resolveSynonym=True):
                        return

                    newRank = self.rank.getRankBelow()
                    newChild = Taxonomy(newTxid, lpsnType, newRank)

                    # save the new node as a child and update type material
                    self._importDirectDescendant(newChild)                      
                    self.typeMaterial = newChild

                    # recurse on the new child
                    newChild.__importTypeHelper(lpsnD)




    ### UPDATING OBJECTS ACCORDING TO BLAST TABLES
    def _updateInternalExternalStatus(self, blastnHits:set, lpsnD:dict) -> None:
        """ updateInternalExternalStatus:
                Accepts a set of taxids (blastn hits) and an LPSN dictionary as
                inputs. Determines the deepest node (furthest from the root)
                that encompasses all of the input taxids. Modifies that object
                so that it and all its descendants are marked as internal. Does
                not return.
        """
        # map the blastn hits to species in the calling object
        blastnSpecies = self.__mapBlastHitsToSpecies(blastnHits, lpsnD)

        # determine the deepest shared parent of the blastn species
        taxonomicMrca = self.__findTaxonomicMrca(blastnSpecies)

        # the highest ranked inteneral ancestor is at least the taxonomicMrca
        deepestInternalAncestor = taxonomicMrca

        # lineage is listed from highest tax to lowest
        # need to reverse it to go through each parent in order
        lineage = taxonomicMrca.getLineage()
        lineage.reverse()
       
        ancestor:Taxonomy
        for ancestor in lineage:
            # get a list of the ancestor's children and loop through them
            allChildren:list = ancestor.getChildren(list)
            child:Taxonomy
            for child in allChildren:
                # if any children are external, then the ancestor is external
                if child.isExternal:
                    isInternal = False
                    break
                
                # if all children are internal, then the ancestor is internal
                else:
                    isInternal = True
            
            # if the ancestor is internal, then it is the deepest internal anc.
            if isInternal: deepestInternalAncestor = ancestor
            
            # otherwise, the previous value is the deepest internal anc.
            else: break

        # set the entire tree as external
        self.getRoot().__setAsExternal()

        # change only the deepest internal ancestor
        deepestInternalAncestor.__setAsInternal(modifyChildren=True)


    def __mapBlastHitsToSpecies(self, blastnHits:set, lpsnD:dict) -> set:
        """ mapBlastHitsToSpecies:
                Accepts a set of taxids (blastn hits) and an LPSN dictionary as
                inputs. Finds the Taxonomy objects that directly correspond to
                the input ids. Returns these objects as a set.
        """
        # constants
        TAX_DB = 'taxonomy'
        TAXID_KEY = 'Id'
        GENUS_KEY = 'Genus'
        SPECIES_KEY = 'Species'
        SPECIES_SYN = 'species_syn'
        INVALID_STR = "(invalid name)"
        GREP_FIND_1 = r'^(\S+) (\S+) .+$'
        GREP_REPL_1 = r'\1 \2'
        GREP_FIND_2 = r' .+$'
        GREP_REPL_2 = r''
        TEMP_TAXID = '0'

        # get the root of the object
        root = self.getRoot()

        # retrieve summaries from NCBI
        try:
            summaries = ncbiSummaryFromIdList(blastnHits, TAX_DB)
        except:
            efetch = ncbiEfetchById(blastnHits, TAX_DB)

            validIds = set()
            for record in efetch:
                validIds.add(record['TaxId'])
            
            summaries = ncbiSummaryFromIdList(validIds, TAX_DB)


        # find species in the object that are present in the blastn table
        speciesPresentInBlastn = set()
        for i in range(len(summaries)):
            # for each summary retrieve the name and id
            summary = summaries[i]
            txid = summary[TAXID_KEY]
            name = summary[GENUS_KEY] + ' ' + summary[SPECIES_KEY]

            # make sure all species names are only two words long
            name = re.sub(GREP_FIND_1, GREP_REPL_1, name)

            # make sure the name is not a synonym!
            if name in lpsnD[SPECIES_SYN].keys():
                name = lpsnD[SPECIES_SYN][name]
            
            # determine what the invalid string would look like
            invalidName = '"' + name + '" ' + INVALID_STR

            # attempt to retrieve a matching object
            foundByTxid = root.getDescendantByTaxId(txid)
            foundByName = root.getDescendantBySciName(name)
            foundByInvalid = root.getDescendantBySciName(invalidName)

            # if an object is found by taxid, then use it
            if foundByTxid:
                newSpecies = foundByTxid
            
            # otherwise, if an object was found by name, then use that one
            # this is primarily used to protect against subranks (eg. 'strain')
            elif foundByName:
                newSpecies = foundByName
            
            elif foundByInvalid:
                newSpecies = foundByInvalid
            
            # otherwise
            else:
                # create a new object
                newSpecies = Taxonomy(txid, name, Taxonomy.SPECIES)

                # update the species with the lpsn data
                newSpecies._addLpsnData(lpsnD)

                # find the parental genus name ...
                ### ... find the genus in the NCBI name if the name is invalid
                if INVALID_STR in newSpecies.sciName:
                    genusName = re.sub(GREP_FIND_2, GREP_REPL_2, \
                                                        newSpecies.ncbiName)
                
                ### ... find the genus in scientific name if the name is valid
                else:
                    genusName = re.sub(GREP_FIND_2, GREP_REPL_2, \
                                                        newSpecies.sciName)
                
                # look for an existing object with the genus name
                parentalGenus = root.getDescendantBySciName(genusName)

                # if a parentalGenus was found ...
                if parentalGenus:
                    # ... then nest the new species under the its parent
                    parentalGenus._importDirectDescendant(newSpecies)


                # if a genus was not found, then parentalGenus will be False
                else:
                    # need to find a new parent for the species
                    try:
                        parentalGenus = newSpecies._findNewParent(lpsnD)
                    except:
                        parentalGenus = False
                        newSpecies = None

                    if parentalGenus:
                        # using a species without a parent means an artificial
                        # taxid may have been created. '0' is reserved, so we
                        # need replace it with an unused artificial taxid
                        if parentalGenus.taxid == TEMP_TAXID:
                            parentalGenus.taxid = \
                                              self.__getUnusedArtificialTaxId()

                        # it is also possible that the taxid is already in use
                        # if that is the case, then we need an artifical taxid
                        ### get all taxids currently in use
                        allTxids = self.getRoot().getAllDescendants(dict)
                        allTxids = set(allTxids.keys())
                        allTxids.add(self.getRoot().taxid)

                        ### if the id is already in use, get an artificial id
                        if parentalGenus.taxid in allTxids:
                            parentalGenus.taxid = \
                                              self.__getUnusedArtificialTaxId()

                        # need to find a safe node to connect this lineage
                        mergeRequired = False
                        newRoot:Taxonomy = parentalGenus.getRoot()
                        while not newRoot.sciName in root:
                            # keep track of the previous deepest ancestor
                            prevRoot = newRoot
                            newRoot  = newRoot._findNewParent(lpsnD)

                            # handle the introduction of artificial taxids
                            if newRoot.taxid == TEMP_TAXID:
                                # get a the next available id from calling obj
                                newId = self.__getUnusedArtificialTaxId()

                                # make sure it is not already in use in newRoot
                                while newId in newRoot.getRoot():
                                    newId = str(int(newId) - 1)

                                newRoot.taxid = newId
                            
                            # a merge is required if the ranks ever match
                            if newRoot.rank == root.rank:
                                mergeRequired = True
                                break
                        
                        # merge the two objects if necessary
                        if mergeRequired:
                            root = Taxonomy._mergeMultipleRoots( \
                                                        [root, newRoot], lpsnD)
                        
                        # otherwise, nest prevRoot beneath root
                        else:
                            root._importExistingSubTax(prevRoot, lpsnD)
            
            # add newSpecies to the growing set (no need to save 'None')
            if newSpecies is not None:
                speciesPresentInBlastn.add(newSpecies)
            
        return speciesPresentInBlastn
    

    def __getUnusedArtificialTaxId(self) -> str:
        """ getUnusedArtificialTaxId:
                Accepts no input. Looks for the presence of any artificial (ne-
                gative integers) taxids in the calling object. Returns the next
                artificial taxid (artificial ids start at -1; decrement by 1)
                that is not currently in use by the calling object.
        """
        # get the set of all ids that are currently in use
        allTxids = set(self.getRoot().getAllDescendants(dict).keys())
        allTxids.add(self.getRoot().taxid)

        # find the minimum value in use
        minTaxId = float('Inf')
        for txid in allTxids:
            if int(txid) < minTaxId:
                minTaxId = int(txid)
    
        # if there are no artificial ids in use, then -1 is unused
        if minTaxId > 0:
            return str(-1)

        # otherwise, return the next lowest integer
        return str(minTaxId - 1)


    def __getAllArtificialTaxIds(self) -> list:
        """ getAllArtificialTaxIds:
                Accepts no inputs. Constructs a list of all taxids in the call-
                ing object that are artificial (<1). Returns the list.
        """
        # get the set of all ids that are currently in use
        allTxids = set(self.getRoot().getAllDescendants(dict).keys())
        allTxids.add(self.getRoot().taxid)

        # go through each id and save those that are artificial (<1)
        artificialTxids = list()
        for txid in allTxids:
            if int(txid) < 1:
                artificialTxids.append(txid)
        
        return artificialTxids


    def __findTaxonomicMrca(self, species:set) -> Taxonomy:
        """ findTaxonomicMrca:
                Accepts a set of Taxonomy objects at the species-level. Finds
                the deepest (lowest rank; furthest from the root) object that
                contains all of the input objects as descendants. Returns the
                identified Taxonomy object.
        """
        # constants
        ERR_1 = "Cannot find taxonomic MRCA. '"
        ERR_2 = "' not within the root of the calling object"

        # get the root
        root = self.getRoot()

        # for each species
        curParentInitialized = False
        spe:Taxonomy
        for spe in species:
            # quit if a species is not witin the calling object's root
            if not root.getDescendantByTaxId(spe.taxid):
                raise BaseException(ERR_1 + str(spe.taxid) + ERR_2)
            
            # if this is the first species, then initialize the parent
            if not curParentInitialized:
                curParent = spe
                curParentInitialized = True
            
            # for every species, update the earliest shared parent
            else:
                # advance the parent until the species is a descendant
                while not curParent.getDescendantByTaxId(spe.taxid):
                    curParent = curParent.parent
        
        return curParent




    ### IMPORTING NCBI ASSEMBLY DATA
    def _importAssemblies(self) -> None:
        """ importAssemblies:
                Accepts no inputs. Retrieves a list of all assemblies for the 
                species contained within the calling object. For each assembly,
                adds the information to the correct species (linked by taxonomy
                id). Does not return.
        """
        # if no species available, then nothing to do
        if self.numDescendantsAtRank(Taxonomy.SPECIES) < 1:
            return

        # obtain a list of assembly ids for the species within the object
        allAssIds  = self.__getAssemblyIds()

        # nothing to do if no assemblies were found
        if allAssIds == []:
            return
        
        # for each assembly ...
        allAssSums = Taxonomy.__getAssemblySummary(allAssIds)
        for assSum in allAssSums:
            # look up the species
            species = self.getDescendantByTaxId(assSum['Taxid'])

            # if the species was not found, then try a different field
            if not species:
                # this handles atypical taxonomy structures in NCBI
                # (sometimes strains are the entries within a species)
                species = self.getDescendantByTaxId(assSum['SpeciesTaxid'])

            # only proceed if a species in the object matches the assSum
            if species:
                strain, strainIsType = species.__getStrainDataFromSummary(assSum)

                if strain is not None:                    
                    # only save if the strain is type material
                    if strainIsType:
                        species._updateAssemblyInfo(assSum)

    
    def __getAssemblyIds(self) -> list:
        """ getAssemblyIds:
                Retrieves a list of NCBI assembly ids for genomes that meet the 
                following criteria:
                    * Organism is within the taxonomic scope of the provided id
                    * It is the latest refSeq version
                    * The assembly was generated from type material
                Returns a list of NCBI assembly ids matching those criteria.
        """
        # constants
        RANK = 'species'
        DATABASE = 'assembly'
        
        # get a list of the taxids for all of the species
        if self.rank > Taxonomy.SPECIES:
            allSpecies = self.getDescendantsAtRank(RANK)
            speciesTxids = list(allSpecies.keys())
        else:
            speciesTxids = [self.taxid]

        # handle an empty id list
        if speciesTxids == []:
            return []

        # construct string to search the assembly database
        searchStr = Taxonomy.__makeAssemblySearchString(speciesTxids)
        
        # search NCBI for valid assemblies ids and return
        return ncbiIdsFromSearchTerm(searchStr, DATABASE)
    
    
    def __makeAssemblySearchString(txids:list) -> str:
        """ makeAssemblySearchString:
                Accepts a list of NCBI taxonomy ids as input. Generates and 
                returns the search string required for identifying the NCBI 
                assemblies that meet the following criteria:
                    * Organism is within the taxonomic scope of the provided id
                    * It is the latest refSeq version
                    * It has annotations
                    * It has full genome representation
                    * It was generated from type material
        """
        # constants
        FIELD_1_PREFIX = 'txid'
        FIELD_1_SUFFIX = '[Organism]'
        FIELD_2 = 'latest_refseq[filter]'
        FIELD_3 = '"has annotation"[Properties]'
        FIELD_4 = '"full genome representation"[Properties]'
        FIELD_5 = '("assembly from type material"[FromType] OR ' \
                    '"assembly from synonym type material"[FromType] OR ' \
                    '"assembly from pathotype material"[FromType] OR ' \
                    '"assembly designated as neotype"[FromType] OR ' \
                    '"assembly designated as reftype"[FromType])'
        OR = ' OR '
        AND = ' AND '
        SEP_LEN  = len(OR)
        
        # convert the list of ids to an 'OR-separated' search term
        searchStr = '('
        for txid in txids:
            searchStr += FIELD_1_PREFIX + str(txid) + FIELD_1_SUFFIX
            searchStr += OR
        searchStr = searchStr[:-SEP_LEN]
        searchStr += ')'

        # add other qualifiers to limit the search results
        searchStr += AND
        searchStr += FIELD_2
        searchStr += AND
        searchStr += FIELD_3
        searchStr += AND
        searchStr += FIELD_4
        searchStr += AND
        searchStr += FIELD_5

        # return
        return searchStr
            

    def __getAssemblySummary(assId) -> list:
        """ getAssemblySummary:
                Accepts an assembly id (str/int) OR a list of assembly ids as
                input. Queries NCBI Assembly via esummary. Returns the parsed
                NCBI result.
        """
        # constants
        DATABASE = 'assembly'
        RES_F1   = 'DocumentSummarySet'
        RES_F2   = 'DocumentSummary'

        # query NCBI
        result = ncbiSummaryFromIdList(assId, DATABASE)

        # parse the result
        return result[RES_F1][RES_F2]
    

    def __getStrainDataFromSummary(self, assSum:Parser.DictionaryElement) \
                                                                      -> tuple:
        """ getStrainDataFromSummary:
                Accepts an assembly summary result (Entrez.esummary) as input.
                Extracts the type strain information for a given assembly summ-
                ary. Returns a tuple composed of the strain name (str) and a
                boolean indicating whether or not the strain is type material.

                This function should only ever be called on a species object.
        """
        # constants
        ERR_MSG = 'calling object is not a species'
        SEP_CHAR = "; "
        STRN_F1 = 'Biosource'
        STRN_F2_A = 'InfraspeciesList'
        STRN_F2_B = 'Isolate'
        STRN_F3 = 'Sub_value'
        TAX_F1 = 'OtherNames'
        TAX_F2 = 'Synonym'
        TAX_F3 = 'ScientificName'

        # ensure that this function has been called on a species object
        if self.rank != Taxonomy.SPECIES:
            raise Exception(ERR_MSG)

        # extract the strain name from the assembly
        allStrainsL = list()
        for strains in assSum[STRN_F1][STRN_F2_A]:
            # make sure that strains is a list of strains (split on SEP_CHAR)
            strains = strains[STRN_F3].split(SEP_CHAR)

            # add the strains to the growing list
            allStrainsL.extend(strains)

        # grab strain from alternate field (split on SEP_CHAR)
        allStrainsL.extend(assSum[STRN_F1][STRN_F2_B].split(SEP_CHAR))
        
        # if strain data still not found
        if allStrainsL == []:
            # query NCBI taxonomy for the full entry
            taxRecord = ncbiEfetchById(assSum['Taxid'], 'taxonomy')

            # only examine one record
            taxRecord = taxRecord[0]
            
            # parse strain name from the record
            try:
                # if 'OtherNames' field exists, parse it
                strains = taxRecord[TAX_F1][TAX_F2]
            except KeyError:
                # otherwise, try the 'ScientificName' field
                strains = list(taxRecord[TAX_F3])

            for strain in strains:
                # remove the species name from the front of the str
                strain = re.sub(self.sciName + ' ', '', strain)
                allStrainsL.append(strain)

        # if still no strains identified, then do not proceed further
        if allStrainsL == []:
            return None, False

        # if there is type material, check if the strains match
        elif self.typeMaterial is not None:
            for strain in allStrainsL:
                # done if a type strain was identified
                if strain in self.typeMaterial:
                    return str(strain), True
        
        # if we reach this point, then return the first strain and False
        return allStrainsL[0], False


    def _updateAssemblyInfo(self, assSummary:Parser.DictionaryElement) -> None:
        """ updateAssemblyInfo:
                Accepts an assembly summary result (Entrez.esummary) as input.
                Updates the assembly information of the calling object with the
                values found in the summary. Does not return.
        """
        # constants
        STATUS  = 'AssemblyStatus'
        COMPLETE = 'Complete Genome'
        REF_SEQ ='RefSeq_category'
        REF_GEN = 'representative genome'
        COVERAGE = 'Coverage'

        # extract some helpful information from the summary
        isComp = assSummary[STATUS]  == COMPLETE
        isRef  = assSummary[REF_SEQ] == REF_GEN
        cover  = Taxonomy.__coverageToFloat(assSummary[COVERAGE])
        strain, strainIsType = self.__getStrainDataFromSummary(assSummary)
        
        # save the assembly info if nothing to overwrite
        if len(self.allAssIds) == 0:
            self.__populateAssemblyVars(assSummary)
        
        # save if a type assembly can replace a nontype assembly
        elif not self.assemblyFromType and strainIsType:
            self.__populateAssemblyVars(assSummary)
        
        # save if a complete genome can replace a draft
        elif not self.assIsComplete and isComp:
            self.__populateAssemblyVars(assSummary)

        # save if a reference genome can replace a non-ref
        elif self.refSeqCategory != REF_GEN and isRef:
            self.__populateAssemblyVars(assSummary)
        
        # save if the coverage is better and its completeness is better (True > False == True)
        elif self.assCoverage < cover and self.assIsComplete <= isComp:
            self.__populateAssemblyVars(assSummary)


    def __populateAssemblyVars(self, assSum:Parser.DictionaryElement) -> None:
        """ populateAssemblyVars:
                Accepts an assembly summary result (Entrez.esummary) as input.
                Uses it to populate the assembly-related member variables of 
                the calling object. Does not return.
        """
        # constants
        STATUS   = 'AssemblyStatus'
        COMPLETE = 'Complete Genome'
        REF_SEQ  = 'RefSeq_category'
        COVERAGE = 'Coverage'

        # get the assembly database's uid
        assId = copy.deepcopy(assSum.attributes['uid'])
        
        # remove all of Bio.Entrez's custom classes
        # this absolutely necessary to maintain copy functionality
        assSum = coerceEntrezToNormal(assSum)

        # get the ftp path
        ftpPath = Taxonomy.__getFtpPathFromAssSummary(assSum)

        # get the accession number
        accn = Taxonomy.__getAccessionFromFtpPath(ftpPath)
    
        # get the coverage as a float
        coverage = Taxonomy.__coverageToFloat(assSum[COVERAGE])

        # get the strain information from the summary
        strain, strainIsType = self.__getStrainDataFromSummary(assSum)

        # populate member variables with assembly data
        self.assemblyFtp = ftpPath
        self.assemblyAccn = accn
        self.assemblyStrain = str(strain)
        self.assemblyFromType = strainIsType
        self.refSeqCategory = assSum[REF_SEQ]
        self.allAssIds.append(assId)
        self.assCoverage = coverage

        if assSum[STATUS] == COMPLETE:
            self.assIsComplete = True


    def __getFtpPathFromAssSummary(assSummary:Parser.DictionaryElement) -> str:
        """ getFtpPathFromAssSummary:
                Accepts an assembly summary result (Entrez.esummary) as input.
                Returns the FTP path for downloading the gbff.gz file from NCBI
                (string).
        """
        # constants
        SLASH = '/'
        FTP_REF_SEQ = 'FtpPath_RefSeq'
        FTP_GENBANK = 'FtpPath_GenBank'
        FTP_SUFFIX = '_genomic.gbff.gz'
        GREP_FIND = r'^.+/([^/]+)/$'
        GREP_REPLACE = r'\1'

        # extract the ftp path from the entry (folder)
        ftpPath = assSummary[FTP_REF_SEQ]

        # if refseq doesn't have it, then get it from genbank
        if ftpPath == "":
            ftpPath = assSummary[FTP_GENBANK]

        # make sure it ends the correct character
        if ftpPath[-1] != SLASH:
            ftpPath += SLASH
        
        # add the file name of the gbff.gz to the path and return
        ftpPath += re.sub(GREP_FIND, GREP_REPLACE, ftpPath)
        ftpPath += FTP_SUFFIX

        return ftpPath


    def __getAccessionFromFtpPath(ftpPath:str) -> str:
        """ getAccessionFromFtpPath:
                Accepts a string indicating the ftp path for the assembly as
                input. Returns a string corresponding to the accession number.
        """
        # expected pattern for the accession number's format
        EXPECT_GREP = r'^GC[AF]_\d{9}$'

        # get the file name of the assembly
        accession = os.path.basename(ftpPath)

        # keep removing "extensions"
        while "." in accession:
            accession = os.path.splitext(accession)[0]
        
        # make sure that the accession has the expected format
        if not re.match(EXPECT_GREP, accession):
            raise RuntimeError('incorrect accession number')

        return accession


    def __coverageToFloat(coverage:str) -> float:
        """ coverageToFloat:
                Accepts a string (assembly coverage) as input. Converts the st-
                ring to a float or raises an error. Returns the float.
        """
        # constants
        ERR_MSG = "assembly coverage cannot be coerced to a numeric"

        # try to convert the string directly
        try:
            coverage = float(coverage)

        # if that's not possible ...
        except:
            # if the coverage is absent (empty string)
            if coverage == "":
                # then set it to the default value (0.0)
                coverage = 0.0

            # if first character is '>', then coerce remaining chars and add 1
            elif coverage[0] == ">":
                coverage = float(coverage[1:]) + 1
            
            # if first character is '<', then coerce remaining chars and sub 1
            elif coverage[0] == "<":
                coverage = float(coverage[1:]) - 1
            
            # otherwise, something unexpected is happening...raise an error
            else:
                raise ValueError(ERR_MSG)
        
        return coverage




    ### PICKING GENOMES FOR PHYLOGENOMIC ANALYSES
    def _pickOutgroup(self, lpsnD:dict) -> tuple:
        """ pickOutgroup:
                Accepts the LPSN dictionary as input. Finds a species that pos-
                sesses an assembly and is external to the ingroup candidates.
                Retrieves the up-to-date handle for the calling object in case
                the calling object was a synononym that became resolved while
                picking the outgroup. Returns a tuple of two Taxonomy objects:
                the first is the up-to-date handle of the calling object and
                the second is the outgroup itself.
        """
        # constants
        BACT_NAME = "Bacteria"
        ARCH_NAME = "Archaea"

        # attempt to retrieve the outgroup
        try:
            outgroup = self.__pickOutgroupWrapper(lpsnD)

        # if this fails, then we must pick from another domain
        except:
            # navigate to the root of the object
            root = self.getRoot()

            # pick a bacterial outgroup if self is archaeal
            if root.sciName == ARCH_NAME:
                outgroup = Taxonomy.__outgroupFromDiffDomain(BACT_NAME, lpsnD)
            
            # pick an archaeal outgroup if self is bacterial
            elif root.sciName == BACT_NAME:
                outgroup = Taxonomy.__outgroupFromDiffDomain(ARCH_NAME, lpsnD)
            
            # raise an error if the root isn't at a domain (should not happen)
            else:
                raise RuntimeError("Could not find an outgroup")

        # get the root
        root = self.getRoot()

        # protect against a discarded taxO (can happen to synonyms)
        if self.parent is not None:
            parent = root.getDescendantByTaxId(self.parent.taxid)
            parent._importDirectDescendant(self)
        
        taxO = root.getDescendantByTaxId(self.taxid, resolveSynonym=True)

        return taxO, outgroup


    def __pickOutgroupWrapper(self, lpsnD:dict) -> Taxonomy:
        """ pickOutgroupWrapper:
                Accepts the LPSN dictionary as input. Finds an external species
                that possesses an assembly and returns it. Relies on a helper
                function to find the correct outgroup. If an outgroup cannot be
                identified, then it will deepen the root until a sibling is av-
                ailable. It will then attempt to find an outgroup within the
                sibling(s). Returns a Taxonomy object that can be used as an
                outgroup.
        """
        # constants
        TEMP_TAXID = '0'
        OR = ' OR '
        
        # get the root
        root = self.getRoot()

        # attempt to pick an outgroup from an existing species
        outgroup = root.__pickOutgroupRecursiveHelper()

        # if the attempt succeeded, then stop and return
        if outgroup is not None:
            return outgroup
        
        # initialize the collection of siblings
        allSiblings = set()

        # if the root is a domain, then it is not possible to go deeper
        # must assume that unrepresented phyla will be suitable for an outgroup
        if root.rank == Taxonomy.DOMAIN:
            # save handles to the current children of the root
            children = root.getChildren(list)

            # remove the children from the root object
            root.descendantsD = dict()

            # repopulate root's descendants with all children taxa
            root._initializeDescendants(maxDepth=root.rank.getRankBelow())
            root._addLpsnData(lpsnD, removeEmptyTaxa=False)

            # mark the root, and all its descendants, as external
            root.__setAsExternal()

            # get a set of the siblings
            allSiblings:set = root.getChildren(set)

            # get a set of the taxids in allSiblings
            allSiblingTaxids = {sib.taxid for sib in allSiblings}

            # get a set of the taxids in children
            childrenTaxids = {kid.taxid for kid in children}

            # only keep the taxids not present in children
            allSiblingTaxids.difference_update(childrenTaxids)

            # revise allSiblings to reflect those in allSiblingTaxids
            allSiblings = {sib for sib in allSiblings if sib.taxid in allSiblingTaxids}

            # nest each of the original children back into the root object
            for child in children:
                root._importDirectDescendant(child)
    
        # keep looping until an outgroup is found. exit via return or exception
        while True:
            # continue to the deepen the root until it has sibling(s)
            while allSiblings == set():
                # save a handle to the previous root
                prevRoot = root

                # deepen the root and add its immediate descendants
                root = root._findNewParent(lpsnD)

                # handle artificial taxids
                if root.taxid == TEMP_TAXID:
                    # initialize a boolean indicating if a synonym was found
                    foundSynonym = False

                    # get the synonym dictionary
                    synToPrefD:dict = lpsnD[str(root.rank) + '_syn']

                    # make the query for searching NCBI taxonomy
                    query = str(root.rank) + "[Rank] AND ("
                    for syn in synToPrefD.keys():
                        # get the preferred name for the synonym
                        prefName = synToPrefD[syn]

                        # if the preferred name matches the object ...
                        if prefName == root.sciName:                    
                            # ... add the synonym to the search string
                            query += syn + OR

                            # ... and indicate that a synonym was found
                            foundSynonym = True
                    
                    # if a synonym was found
                    if foundSynonym:
                        # remove trailing OR and close the parentheses
                        query = query[:-len(OR)] + ")"
                        
                        # attempt to retrieve possible alternative taxids
                        synonymTaxIds = ncbiIdsFromSearchTerm(query, 'taxonomy')

                        # use one of the possible taxids
                        root.taxid = str(synonymTaxIds.pop())

                        # if there are additional taxids ...
                        if len(synonymTaxIds) > 0:
                            # ... then add each one to the taxidSynonyms dictionary
                            for txid in synonymTaxIds:
                                root.taxidSynonyms[str(txid)] = root.taxid

                        # look up the ncbi name for the new taxid
                        summary = ncbiSummaryFromIdList([root.taxid], 'taxonomy')
                        root.ncbiName = str(summary[0]['ScientificName'])

                        # set ncbiNameIsCorrect to False
                        root.ncbiNameIsCorrect = False
                    
                    # if a synonym was not found
                    else:
                        # use an artificial taxid
                        root.taxid = root.__getUnusedArtificialTaxId()

                        # loop again
                        continue

                # get all of the direct descendants of the new root
                root._initializeDescendants(maxDepth=root.rank.getRankBelow())

                # make sure the imported descendants are congruent with LPSN
                root._addLpsnData(lpsnD, removeEmptyTaxa=False)

                # the previous steps displaced prevRoot from root's descendants
                # reinsert prevRoot into the root and get its siblings
                root._importDirectDescendant(prevRoot)

                # make sure duplicate siblings were not introduced by NCBI
                # example from test case 15 (GCF_900109565.1):
                #   Sphingomonadales (taxid: 204457) -> Flavobacteriales
                #   findNewParent retrieved: Flavobacteriales (taxid: 200644)
                root.__resolveDuplicateNames()

                # get a set of the new siblings of the previous root
                allSiblings = prevRoot.getSiblings(set)
    
                # addLpsnData could merge sibling with prevRoot
                # let that happen before exiting from while-loop
                if len(allSiblings) > 0:
                    # initialize a collection of faux siblings
                    fauxSiblings = set()

                    # for each of the new siblings
                    sibling:Taxonomy
                    for sibling in allSiblings:
                        # if the sibling's sciName matches the prevRoot
                        if sibling.sciName == prevRoot.sciName:
                            # this is a faux sibling, add it to the set
                            fauxSiblings.add(sibling)
                        
                        # if the sibling doesn't have a type species
                        if sibling.numDescendantsAtRank('species') < 1:
                            # try to populate the sibling first
                            sibling._initializeDescendants(maxDepth='species')
                            sibling._addLpsnData(lpsnD, removeEmptyTaxa=False)

                            # check again
                            if sibling.numDescendantsAtRank('species') < 1:
                                # this is a faux sibling, add it to the set
                                fauxSiblings.add(sibling)
                        
                        # otherwise, mark the sibling as external
                        else:
                            sibling.__setAsExternal()
                    
                    # remove any faux siblings from the set of siblings
                    allSiblings.difference_update(fauxSiblings)

                    # addLpsnData can produce duplicate sibling names
                    # make sure these names are resolved before proceeding
                    root.__resolveDuplicateNames()

            # attempt to use the type species of the sibling as an outgroup
            # type material was populated to the species-level by addLpsnData
            sibling:Taxonomy
            for sibling in allSiblings:
                sibling._importAssemblies()

                # attempt to find an outgroup; stop looping upon success
                outgroup = root.__pickOutgroupRecursiveHelper()
                if outgroup is not None:
                    return outgroup
        
            # if none of the existing type species were suitable as outgroups,
            # then try any of the species of an external type genus
            for sibling in allSiblings:
                typeGen = sibling.getTypeAtRank('genus')
                if typeGen is not None:
                    typeGen._initializeDescendants()
                    typeGen._addLpsnData(lpsnD)
                    typeGen.__setAsExternal()
                    typeGen._importAssemblies()

                    # attempt to find an outgroup; stop looping upon success
                    outgroup = typeGen.__pickOutgroupRecursiveHelper()
                    if outgroup is not None:
                        return outgroup
            
            # if the type genus could not be used, 
            # then try any species of the type family
            for sibling in allSiblings:
                typeFam = sibling.getTypeAtRank('family')
                if typeFam is not None:
                    typeFam._initializeDescendants()
                    typeFam._addLpsnData(lpsnD)
                    typeFam.__setAsExternal()
                    typeFam._importAssemblies()

                    # attempt to find an outgroup; stop looping upon success
                    outgroup = typeFam.__pickOutgroupRecursiveHelper()
                    if outgroup is not None:
                        return outgroup
            
            # if the type family could not be used,
            # then try any species descendant to the sibling (last resort)
            for sibling in allSiblings:
                sibling._initializeDescendants()
                sibling._addLpsnData(lpsnD)
                sibling.__setAsExternal()
                sibling._importAssemblies()

                # attempt to find an outgroup; stop looping upon success
                outgroup = sibling.__pickOutgroupRecursiveHelper()
                if outgroup is not None:
                    return outgroup


            
            # if we're here, the siblings weren't good enough. loop
            allSiblings = set()

            # if we're here and the root is a domain, then the function failed
            if root.rank == Taxonomy.DOMAIN:
                raise RuntimeError("Could not find an outgroup")


    def __pickOutgroupRecursiveHelper(self) -> Taxonomy:
        """ pickRecursiveOutgroupHelper:
                Accepts no inputs. Requires that the calling object has been 
                successfully reconciled with LPSN, and that NCBI assembly info
                has been mapped onto the species. Gets a list of children which
                are marked as external. Goes through the list, starting with
                the type material, and recurses on the first child object that
                contains one or more assembly. The base case occurs when the 
                calling object is a species. Returns a Taxonomy object that can
                be used as an outgroup (with respect to those objects which are
                not "external").
        """
        # get a list of external children
        extKids:list = self.getExternalChildren(list)
        
        # make sure the list begins with typeMaterial, if possible
        if self.typeMaterial is not None:
            if self.typeMaterial in extKids:
                idx = extKids.index(self.typeMaterial)
                typ = extKids.pop(idx)
                extKids.insert(0, typ)
            
        # for each external child (starting with type if possible)
        child:Taxonomy
        for child in extKids:
            # if there are assemblies
            if child.numSavedAssemblies() > 0:
                # if the calling object is a species, then this is the outgroup
                if child.rank == Taxonomy.SPECIES:
                    return child
                
                # otherwise, recurse on the children, 1-by-1
                outgroup = child.__pickOutgroupRecursiveHelper()

                # return if an outgroup is identified, otherwise keep looping
                if outgroup is not None:
                    return outgroup
                
        return None


    def __outgroupFromDiffDomain(name:str, lpsnD:dict) -> Taxonomy:
        """ outgroupFromDiffDomain:
                Accepts one of two domain names ('Bacteria' or 'Archaea') and 
                the lpsnD dictionary as inputs. Randomly selects an order from
                the requested domain that contains at least one type species
                with an assembly available. Returns one of the species from th-
                at order as a Taxonomy object.
        """
        # constants
        SEARCH_SUFFIX = " AND superkingdom[RANK]"
        DATABASE = 'taxonomy'
        ORDER = 'order'

        # determine the taxid ny searching NCBI
        result = ncbiIdsFromSearchTerm(name + SEARCH_SUFFIX, DATABASE)
        taxid  = result.pop()
        
        # initialize the descendants to the order
        taxO = Taxonomy(taxid, name, Taxonomy.DOMAIN)
        taxO._initializeDescendants(maxDepth=ORDER)

        # get a list of all the orders
        allOrders = taxO.getDescendantsAtRank(ORDER, list)

        # get a shuffled list of indices
        indices = list(range(len(allOrders)))
        random.shuffle(indices)

        # initialize a boolean to handle failed attempts
        found = False

        # go through the list of shuffled indices
        for idx in indices:
            # extract the order from the list and make it a root
            order:Taxonomy = allOrders[idx]
            order.parent = None

            # import the lpsn data and the assembly data
            order._addLpsnData(lpsnD, removeEmptyTaxa=False)
            order._importAssemblies()

            # if there is an assembly available, then we're done searching
            if len(order.getAllSpeciesWithAssemblies()) > 0:
                found = True
                break
        
        # if an assembly was not found, then raise an error
        if not found:
            raise RuntimeError("Could not find an outgroup")
        
        # set the order as external
        order.__setAsExternal()

        # pick an outgroup
        # All species are external and at least one assembly is present
        # because of this, we can just call the recursive helper
        outgroup = order.__pickOutgroupRecursiveHelper()

        return outgroup


    def _pickIngroup(self, maxNumSeqs:int) -> list:
        """ pickIngroup:
                Accepts an integer indicating the maximum number of ingroup sp-
                ecies to pick as input. Returns a list containing the ingroup
                species.
        """
        # constants
        MIN_NUM_SEQS = 2

        # make sure the max number is greater than the minimum
        if maxNumSeqs < MIN_NUM_SEQS:
            raise ValueError('max sequences must be greater than ' + \
                                                             str(MIN_NUM_SEQS))

        # evenly distribute the sequences across the "internal" Taxonomy
        ingroupRaw:list
        candidates:list
        ingroupRaw, candidates = self.__evenlyDistributeIngroup(maxNumSeqs)

        # unpack the raw ingroup result into lists of taxa that have contribut-
        # ed all possible members
        ingroupSpecies:list
        ingroupHigherTaxa:list
        ingroupSpecies, ingroupHigherTaxa = \
                                        Taxonomy.__unpackIngroupRaw(ingroupRaw)

        # redistribute the surplus to the remaining candidates
        Taxonomy.__redistributeSurplusToCandidates(ingroupSpecies, \
                                                 ingroupHigherTaxa, candidates)

        # check that all ingroup species are unique
        if len(ingroupSpecies) != len(set(ingroupSpecies)):
            raise RuntimeError('ingroup contains redundant species!')

        # get a set of all the genera represented by the ingroup species
        ingroupGenera = set()
        spe:Taxonomy
        for spe in ingroupSpecies:
            ingroupGenera.add(spe.parent)

        # find all internal genera with assemblies
        allCandidateGenera:list = self.getInternalDescendantsAtRank('genus')
        canGenIndices = list(range(len(allCandidateGenera)))
        canGenIndices.reverse()  # tail-to-head allows on-the-fly popping
        for genIdx in canGenIndices:
            genTaxon:Taxonomy = allCandidateGenera[genIdx]
            if genTaxon.numSavedAssemblies() < 1:
                allCandidateGenera.pop(genIdx)
        
        # find the candidate genera that are not represented by the ingroup
        missingGenera:set = set(allCandidateGenera)
        missingGenera.difference_update(ingroupGenera)

        # if there are missing genera AND some genera have multiple species
        if len(missingGenera) > 0 and len(ingroupGenera) < len(ingroupSpecies):
            # redistribute some sequences to some unrepresented genera
            ingroupSpecies = Taxonomy.__improveIngroupGenusDiversity(\
                                   candidates, ingroupSpecies, ingroupGenera, \
                                              missingGenera, ingroupHigherTaxa)

        return ingroupSpecies


    def __evenlyDistributeIngroup(self, maxNumSeqs:int) -> tuple:
        """ evenlyDistributeIngroup:
                Accepts an integer indicating the maximum number of ingroup
                species to pick as input. Iteratively distributes the the sequ-
                ences across all possible ingroup candidates. The iterations
                continue until the ingroup being created stops changing. Retur-
                ns a tuple consisting of two lists of dictionaries: one for the
                ingroup and one for the candidates that could not be evenly
                distributed across without supervision. These lists will enable
                downstream functions to keep track of both the parental taxa 
                and the total number of species that they should contribute to
                the ingroup.
        """
        # constants
        TAX = 'tax'
        SURPLUS = 'surplus'
        DONE = 'done'

        # evenly distribute the sequences across the taxonomy
        ingroup:list
        candidates:list
        ingroup, candidates = self.__getIngroupDistribution(maxNumSeqs)

        # make sure we will try to redistribute at least once
        prevIngroup = float('Inf')

        # loop until the number of ingroup assemblies stops changing
        # (unpacking the ingroup reveils the number of assemblies)
        while prevIngroup != ingroup:
            # however, if ever there are no candidates, we're done
            if candidates == []:
                break

            # update the previous ingroup
            prevIngroup = copy.deepcopy(ingroup)
            
            # transfer all of the surplus to a dictionary of the parents
            parentD = dict()
            Taxonomy.__transferSurplusToParentD(candidates, parentD)
            Taxonomy.__transferSurplusToParentD(ingroup, parentD)
        
            # make a copy we can modify without affecting the original
            taxCopy = self.copy()

            # for each item in the ingroup, remove it from the copy
            # this will allow the distribution functions to completely ignore
            # the taxa that have already been added to the ingroup
            for ingMemb in ingroup:
                # retrieve the taxon in the copy that matches the ingroup taxon
                match = taxCopy.getDescendantByTaxId(ingMemb[TAX].taxid)

                # find the parent and use it to delete the ingroup taxon
                parent = match.parent
                parent.__removeDirectDescendant(match.taxid)

                # if we just removed a taxon present in parentD,
                # then we need to move that surplus elsewhere
                if match.taxid in parentD.keys():
                    # update an existing parent if possible
                    if parent.taxid in parentD.keys():
                        parentD[parent.taxid] += parentD[match.taxid]
                    
                    # otherwise, create a new parent
                    else:
                        parentD[parent.taxid] = parentD[match.taxid]
                    
                    # either way, delete the match from parentD
                    del parentD[match.taxid]

            # sort the parent keys to control the order of iteration
            # we want lower taxa to be processed before their ancestors
            rootCopy = taxCopy.getRoot()
            sortByTaxRank = lambda x: rootCopy.getDescendantByTaxId(x).rank
            parentIdsSorted = sorted(list(parentD.keys()), key=sortByTaxRank)

            # for each parent, redistribute the surplus
            newIngroup = list()
            newCandidates = list()
            for parentId in parentIdsSorted:
                # find the equivalent parent in the copied object
                parentCopy = rootCopy.getDescendantByTaxId(parentId)

                # redistribute its surplus
                result = parentCopy.__getIngroupDistribution(parentD[parentId])

                # remove new ingroup members to protect against duplicates
                newEntry:dict
                for newEntry in result[0]:
                    newTax:Taxonomy = newEntry[TAX]
                    newParent = newTax.parent
                    newParent.__removeDirectDescendant(newTax.taxid)

                # split the resulting tuple and add it to the growing lists
                newIngroup += result[0]
                newCandidates += result[1]

            # update the ingroup with the new results
            for idx in range(len(newIngroup)):
                # the new list contains copies of the taxa; get the originals
                realTax = self.getDescendantByTaxId(newIngroup[idx][TAX].taxid)

                # build the new ingroup entry
                entry = realTax.__makeIngroupCandidateEntry(
                                   newIngroup[idx][SURPLUS], newIngroup[idx][DONE])

                # get the indices backwards so that list modification can occur
                indices = list(range(len(ingroup)))
                indices.reverse()

                # remove any of the new taxon's children from the ingroup
                for idx in indices:
                    # descendants would be redundant, so we need to remove them
                    if realTax.getDescendantByTaxId(ingroup[idx][TAX].taxid):
                        ingroup.pop(idx)

                # add the entry to the ingroup list
                ingroup.append(entry)
        
            # completely replace the candidate list with the new result
            oldCandidates = copy.deepcopy(candidates)
            candidates = list()
            for idx in range(len(newCandidates)):
                # the new list contains copies of the taxa; get the originals
                realTax = self.getDescendantByTaxId(
                                    newCandidates[idx][TAX].taxid)

                # build the new candidate entry
                entry = realTax.__makeIngroupCandidateEntry(
                            newCandidates[idx][SURPLUS], newCandidates[idx][DONE])

                # add it to the list
                candidates.append(entry)
            
            # make sure any old candidates with surplus weren't lost
            for oldCand in oldCandidates:
                if oldCand[SURPLUS] > 0:
                    candidates.append(oldCand)

        # return once the ingroup has stopped changing
        return ingroup, candidates


    def __getIngroupDistribution(self, maxNumSeqs:int) -> tuple:
        """ getIngroupDistribution:
                Accepts an integer indicating the maximum number of ingroup
                species to pick as input. Calls a recursive helper function to
                distribute the input number across the calling object's ingroup
                candidate species. Returns a tuple consisting of two lists of
                dictionaries: one for the ingroup and one for the candidates
                that could not be evenly distributed across without any superv-
                ision.
        """
        # constants
        DONE = 'done'

        # get the full list
        combinedList = self.__getIngroupDistributionHelper(maxNumSeqs)

        # use reversed indices so popping can happen on-the-fly
        indices = list(range(len(combinedList)))
        indices.reverse()

        # split the full list into the separate ingroup and candidate lists
        ingroup = list()
        candidates = list()
        for i in indices:
            # if marked as 'done', then add it to the ingroup
            if combinedList[i][DONE]:
                ingroup.append(combinedList.pop(i))
            
            # otherwise, it is unfinished so add it to the candidates
            else: candidates.append(combinedList.pop(i))
        
        return ingroup, candidates


    def __getIngroupDistributionHelper(self, maxNumSeqs:int) -> list:
        """ getIngroupDistributionHelper:
                Accepts an integer indicating the maximum number of ingroup
                species to pick as input. Recursively distributes the value un-
                til the calling object cannot evenly distribute the input to
                its children, or until all of the surplus has been distributed
                to all of the calling object's species. Returns a list of dict-
                ionaries indicating how the ingroup should be selected from the
                Taxonomy object.
        """
        # constants
        SURPLUS = 'surplus'
        
        # intiialize values for looping
        intChildren = self.getInternalChildren(set)
        intChild:Taxonomy
        ingroupCandidates = set()
        numIngroupAssemblies = 0

        # if there are no internal children, 
        # then keep trying lower ranks until internal descendants are found
        # start exactly two ranks below current
        curRank:TaxRank = self.rank - 2
        while intChildren == set() and curRank is not None:
            intChildren:set = self.getInternalDescendantsAtRank(curRank, set)

            # exit the loop if curRank is at the rank floor
            if curRank == Taxonomy.SPECIES:
                break
                
            # otherwise, move on to the next rank
            else:
                curRank.decrement()

        # for each internal child of the calling object
        for intChild in intChildren:
            # if the there are assemblies, then the child is a candidate
            if intChild.numSavedAssemblies() > 0:
                ingroupCandidates.add(intChild)
                numIngroupAssemblies += intChild.numSavedAssemblies()
        
        # if there are less than (or equal to) the number requested, we're done
        if numIngroupAssemblies <= maxNumSeqs:
            # remember, the expected output is list
            return [self.__makeIngroupCandidateEntry(
                                    (maxNumSeqs - numIngroupAssemblies), True)]
        
        # get the share of sequences that each candidate is allotted
        share = maxNumSeqs // len(ingroupCandidates)
        leftover = maxNumSeqs % len(ingroupCandidates)

        # recurse on the ingroup candidates with their share before returning
        if share > 0:
            outL = list()
            igc:Taxonomy
            for igc in ingroupCandidates:
                outL+= igc.__getIngroupDistributionHelper(share)
            
            # we need to save the surplus leftover from the quotient division
            # arbitrarily add it to the last candidate in the list
            outL[-1][SURPLUS] += leftover

            return outL
        
        # otherwise there is not enough for all children; stop and handle later
        # remember, the expected output is list
        else: return [self.__makeIngroupCandidateEntry(maxNumSeqs, False)]


    def __makeIngroupCandidateEntry(self:Taxonomy, surplus:int, 
                                                    done:bool) -> dict:
        """ makeIngroupCandidateEntry:
                Accepts an integer indicating the maximum number of ingroup
                species to pick and a boolean indicating whether nor not distr-
                ibution to the candidate is done as inputs. Creates and returns
                the entry as a dictionary.
        """
        # constants (keys)
        TAX = 'tax'
        SURPLUS = 'surplus'
        DONE = 'done'

        # create and return the entry
        entry = dict()
        entry[TAX] = self
        entry[SURPLUS] = surplus
        entry[DONE] = done

        return entry


    def __unpackIngroupRaw(ingroupRaw:list) -> tuple:
        """ unpackIngroupRaw:
                Accepts a list of dictionaries for those candidates taht have
                been selected for the ingroup. Returns two lists of Taxonomy 
                objects: one is a list of all the species that have currently
                been selected; the other is a list of all the higher taxa that
                have had all of their ingroup species selected.
        """
        # constants
        TAX = 'tax'

        # initialize variables to hold unpacked data
        ingroupSpecies = list()
        ingroupHigherTaxa = list()

        # for each item in the list
        for member in ingroupRaw:
            # extract the Taxonomy object
            taxO:Taxonomy = member[TAX]

            # update the list of ingroup species
            ingroupSpecies += taxO.getAllIngroupCandidateSpecies(list)

            # update the list of ingroup higher taxa
            if taxO.rank > Taxonomy.SPECIES:
                ingroupHigherTaxa.append(taxO)
        
        return ingroupSpecies, ingroupHigherTaxa


    def __transferSurplusToParentD(entryL:list, parentD:dict) -> None:
        """ transferSurplusToParentD:
                Accepts a list of dictionaries (candidate entries) and a dicti-
                onary as inputs. Transfers the surplus stored in the entries to
                the corresponding parent in the parent dictionary. Modifies the
                inputs but does not return.
        """
        # constants
        TAX = 'tax'
        SURPLUS = 'surplus'

        # for each entry
        for idx in range(len(entryL)):
            # parse data into more descriptive variable names
            entry = entryL[idx]
            taxO:Taxonomy = entry[TAX]
            surplus:int   = entry[SURPLUS]

            # if at the root, then no parent to transfer to
            if not taxO.isRoot():
                # get the parent's taxid
                parentId = taxO.parent.taxid

                # if the parent doesn't exist yet
                if parentId not in parentD.keys() and surplus > 0:
                    # make a new entry with the surplus as the value
                    parentD[parentId] = surplus
                
                # otherwise update the surplus
                elif surplus > 0:
                    parentD[parentId] += surplus
                
                # reset entry's surplus (use index to ensure modification)
                entryL[idx][SURPLUS] = 0


    def __redistributeSurplusToCandidates(ingroupSpecies:list, \
                              ingroupHigherTaxa:list, candidates:list) -> None:
        """ redistributeSurplusToCandidates:
                Accepts three lists of inputs. Two are lists of Taxonomy objec-
                ts, and the third is a list of candidate entries. Iteratively 
                distributes the surplus to the candidates. The iteration conti-
                nues until all of the surplus in the candidates has been fully
                distributed. Modifies the inputs but does not return.
        """
        # constants
        TAX = 'tax'
        SURPLUS = 'surplus'
        MAX_CHANCES_BEFORE_AVOIDING = 24

        # calculate total surplus contained in the candidates list
        sumCand = 0
        for can in candidates:
            sumCand += can[SURPLUS]
        
        # calculate the maximum number of sequences based on the input
        maxNumSeqs = sumCand + len(ingroupSpecies)
        
        # tracks the number of times a taxid has been seen
        # avoids the taxon if it exceeds the limit
        timesTried = dict()
        avoid = set()

        # while candidates still remain, continue picking species
        while candidates != []:
            # try to pick the type species's assembly if possible
            Taxonomy.__pickTypeAssembly(ingroupSpecies, ingroupHigherTaxa, 
                                        candidates)

            # if that depleted the candidates, then we're done
            if candidates == []: break

            # track the taxids of the ingroup for faster Taxonomy look-ups
            ### track the species
            ingroupTaxIds = set()
            ingSpe:Taxonomy
            for ingSpe in ingroupSpecies:
                ingroupTaxIds.add(ingSpe.taxid)

            ### track the higher taxa as well
            ingTax:Taxonomy
            for ingTax in ingroupHigherTaxa:
                ingroupTaxIds.add(ingTax.taxid)

            # prepare to replace existing candidates with their children
            newCandidates = list()

            # track the taxids of the new candidates for index look-ups
            newCandidateTaxIds = list()

            # begin looping through the candidates (tail to head allows pop)
            candidateIndices = list(range(len(candidates)))
            candidateIndices.reverse()
            for canIdx in candidateIndices:
                # extract the data into smaller variable names
                candTaxon:Taxonomy = candidates[canIdx][TAX]
                candSurplus:int    = candidates[canIdx][SURPLUS]

                # update the 'timesTried' dictionary
                if candTaxon.taxid in timesTried.keys():
                    timesTried[candTaxon.taxid] += 1

                    # if it appears that we are infinitely looping
                    if timesTried[candTaxon.taxid] > \
                                                   MAX_CHANCES_BEFORE_AVOIDING:

                        # we can't avoid the root! we need a new solution
                        if candTaxon.isRoot():
                            # get a list of all the remaining candidates
                            allCandSpe:set = candTaxon.getAllIngroupCandidateSpecies(set)
                            remCandSpe = list(allCandSpe.difference(set(ingroupSpecies)))

                            # pick candidates until the surplus is consumed
                            for numIter in range(candSurplus):
                                randIdx = random.choice(range(len(remCandSpe)))
                                ingroupSpecies.append(remCandSpe.pop(randIdx))
                            
                            # remove the taxon from the candidate list and move on
                            candidates.pop(canIdx)
                            continue
                        
                        # otherwise, make sure that the taxon is avoided in the future
                        else:
                            avoid.add(candTaxon.taxid)
                
                # the first time we see an id, it hasn't been tried yet.
                else: timesTried[candTaxon.taxid] = 0

                # if the taxon is avoided ...
                if candTaxon.taxid in avoid:
                    # pass the surplus to an existing parental 'new candidate'
                    if candTaxon.parent.taxid in newCandidateTaxIds:
                        newCanIdx = newCandidateTaxIds.index(
                                                        candTaxon.parent.taxid)
                        newCandidates[newCanIdx][SURPLUS] += candSurplus
            
                    # transfer surplus to a new parental 'new candidate'
                    else:
                        entry = candTaxon.parent.__makeIngroupCandidateEntry(
                                                            candSurplus, False)
                        newCandidates.append(entry)
                        newCandidateTaxIds.append(candTaxon.parent.taxid)
                
                # if the candidate is a species
                elif candTaxon.rank == Taxonomy.SPECIES:
                    # if the candidate isn't already in the ingroup
                    if candTaxon.taxid not in ingroupTaxIds:
                        # move it to the ingroup
                        ingroupSpecies.append(candTaxon)
                        ingroupTaxIds.add(candTaxon.taxid)
                    
                    # otherwise, if its parent is a new candidate
                    elif candTaxon.parent.taxid in newCandidateTaxIds:
                        # transfer the surplus to the existing new candidate
                        newCanIdx = newCandidateTaxIds.index(
                                                        candTaxon.parent.taxid)
                        newCandidates[newCanIdx][SURPLUS] += candSurplus
                    
                    # otherwise, make the parent a new candidate
                    else:
                        entry = candTaxon.parent.__makeIngroupCandidateEntry(
                                                            candSurplus, False)
                        newCandidates.append(entry)
                        newCandidateTaxIds.append(candTaxon.parent.taxid)

                # if the candidate is not a species
                else:
                    # get all the children that would make suitable candidates
                    allCandChildren:list = candTaxon.getInternalChildren(list)

                    # curate the list so only suitable candidates are present
                    # for each of the children (tail to head allows pop)
                    childIndices = list(range(len(allCandChildren)))
                    childIndices.reverse()
                    for childIdx in childIndices:
                        # get the child object
                        child:Taxonomy = allCandChildren[childIdx]
        
                        # candidates must have at least 1 assembly saved
                        if child.numSavedAssemblies() < 1:
                            allCandChildren.pop(childIdx)

                        # candidates must not be current ingroup members
                        elif child.taxid in ingroupTaxIds:
                            allCandChildren.pop(childIdx)

                        # candidates must not be currently avoided
                        elif child.taxid in avoid:
                            allCandChildren.pop(childIdx)
                        
                        # for candidates at higher ranks ...
                        elif child.rank > Taxonomy.SPECIES:
                            # get a set of species not yet in the ingroup
                            potentialIngroupSpecies:set = \
                                                child.getAllIngroupCandidateSpecies(set)
                            potentialIngroupSpecies.difference_update(
                                                                ingroupSpecies)

                            # candidates must contain at least 1 new species
                            if len(potentialIngroupSpecies) < 1:
                                allCandChildren.pop(childIdx)
                                ingroupHigherTaxa.append(child)
                                ingroupTaxIds.add(child.taxid)

                    # if there are no suitable candidates
                    if allCandChildren == []:
                        # then return the surplus to its parent
                        if candTaxon.parent is not None:
                            entry = \
                                candTaxon.parent.__makeIngroupCandidateEntry(
                                                            candSurplus, False)
                            newCandidates.append(entry)
                            newCandidateTaxIds.append(candTaxon.parent.taxid)
                        
                        # if there is no parent, then pick random internal 
                        # assemblies not already in the ingroup
                        else:
                            # find all the species still available to select
                            allCandidateSpecies:set = \
                                            candTaxon.getAllIngroupCandidateSpecies(set)
                            allCandidateSpecies.difference_update(
                                                                ingroupSpecies)

                            # pick 'surplus' species at random
                            selectedCandidateSpecies = random.choices(
                                      list(allCandidateSpecies), k=candSurplus)

                            # make new entries and add them to the list
                            spe:Taxonomy
                            for spe in selectedCandidateSpecies:
                                entry = spe.__makeIngroupCandidateEntry(1, 
                                                                        False)
                                newCandidates.append(entry)
                                newCandidateTaxIds.append(spe.taxid)
                
                    # if there are suitable candidates
                    else:
                        # randomly select 'surplus' children
                        newCandTaxa:list = random.choices(allCandChildren, 
                                                                k=candSurplus)

                        # build the new candidate entry and add it to the list
                        newCandTaxon:Taxonomy
                        for newCandTaxon in newCandTaxa:
                            entry = newCandTaxon.__makeIngroupCandidateEntry(1, 
                                                                        False)
                            newCandidates.append(entry)

                # remove the candidate from the candidates list
                candidates.pop(canIdx)
            
            # add the new candidates to the candidates list
            candidates += newCandidates
                
        # check that the maxNumSeqs has not been modified by this function
        if len(ingroupSpecies) != maxNumSeqs:
            raise RuntimeError("'ingroupSpecies' is not the expected length")


    def __pickTypeAssembly(ingroupSpecies:list, ingroupHigherTaxa:list, \
                                                      candidates:list) -> None:
        """ pickTypeAssembly:
                Accepts three lists of inputs. Two are lists of Taxonomy objec-
                ts, and the third is a list of candidate entries. Retrieves the
                type species for the Taxonomy objects in the candidates list.
                Adds it to the ingroup if it is not there already. Modifies the
                inputs but does not return.
        """
        # constants
        TAX = 'tax'
        SURPLUS = 'surplus'

        # track the ingroup species taxids for faster 'in' Taxonomy look-ups
        ingroupTaxIds = set()
        ingSpe:Taxonomy
        for ingSpe in ingroupSpecies:
            ingroupTaxIds.add(ingSpe.taxid)

        # prepare to loop through the list from the tail to head (allows pop)
        canIndices = list(range(len(candidates)))
        canIndices.reverse()

        # for each candidate
        for canIdx in canIndices:
            # get the type species of the candidate
            canTax:Taxonomy = candidates[canIdx][TAX]
            typSpe = canTax.getTypeAtRank(Taxonomy.SPECIES)

            if typSpe is not None:
                # if the type species has an assembly and isn't in the ingroup
                if typSpe.taxid not in ingroupTaxIds \
                                           and typSpe.numSavedAssemblies() > 0:
                    # ... add it to the ingroup and reduce surplus by one
                    ingroupSpecies.append(typSpe)
                    ingroupTaxIds.add(typSpe.taxid)
                    candidates[canIdx][SURPLUS] -= 1
                
            # if the surplus of a candidate is 0, then remove it from the list
            if candidates[canIdx][SURPLUS] == 0:
                candidates.pop(canIdx)
                if canTax.rank > Taxonomy.SPECIES:
                    ingroupHigherTaxa.append(canTax)


    def __improveIngroupGenusDiversity(candidates:list, ingroupSpecies:list, \
                ingroupGenera:set, missingGenera:set, ingroupHigherTaxa:list) \
                                                                       -> list:
        """ improveIngroupGenusDiversity:
                Accepts a list of candidate entries (dictionaries), two lists
                of Taxonomy objects, a set of the genera represented by the in-
                group species, and a set of the genera not represented by the
                ingroup species as inputs. Attempts to increase the number of
                unique, ingroup genera represented by the ingroup species. Mod-
                ifies the inputs. However, one list (ingroupSpecies) is reassi-
                gned within the function, so this list is returned.
        """
        # constants
        TAX = 'tax'
        SURPLUS = 'surplus'

        # get a list of the higher taxa's taxids (for faster TaxO lookups)
        ingroupHigherTaxIds = list()
        ingTax:Taxonomy
        for ingTax in ingroupHigherTaxa:
            ingroupHigherTaxIds.append(ingTax.taxid)

        # initialize variables for looping
        ingGen:Taxonomy
        newSurplus = 0
        multiSpeGenera = dict()

        # for each ingroup genus
        for ingGen in ingroupGenera:
            # get a list of all its children currently in the ingroup
            children:set = ingGen.getChildren(set)
            children.intersection_update(ingroupSpecies)

            # if this genus harbors multiple ingroup species, then redistribute
            if len(children) > 1:
                # save the genus and the remaining number of children
                multiSpeGenera[ingGen.taxid] = len(children) - 1

                # remove it from ingroupHigherTaxa (if present)
                if ingGen.taxid in ingroupHigherTaxIds:
                    idx1 = ingroupHigherTaxa.index(ingGen)
                    ingroupHigherTaxa.pop(idx1)
                    idx2 = ingroupHigherTaxIds.index(ingGen.taxid)
                    ingroupHigherTaxIds.pop(idx2)

                # create a new entry for the genus and add it to the candidates
                entry = ingGen.__makeIngroupCandidateEntry(1, False)
                candidates.append(entry)

                # update the accumulating surplus
                newSurplus += len(children) - 1  # we added one species already

                # remove all of the children from the ingroup
                temp = set(ingroupSpecies)
                temp.difference_update(children)
                ingroupSpecies = list(temp)

        # if there is more surplus than missing genera
        if newSurplus > len(missingGenera):
            # add one species for each missing genus
            allNewEntries = dict()
            misGen:Taxonomy
            for misGen in missingGenera:
                entry = misGen.__makeIngroupCandidateEntry(1, False)
                allNewEntries[misGen.taxid] = entry
                newSurplus -= 1
            
            # keep adding back to the multi-species genera until full
            while sum(multiSpeGenera.values()) > 0 and newSurplus > 0:
                for mulSpeGenTaxId in multiSpeGenera.keys():
                    # if there are still more species left for this genus
                    if multiSpeGenera[mulSpeGenTaxId] >= 1:
                        # find the candidate index
                        canIdx = None
                        for idx in range(len(candidates)):
                            if candidates[idx][TAX].taxid == mulSpeGenTaxId:
                                canIdx = idx
                                break
                        
                        # transfer one species back to the genus
                        candidates[canIdx][SURPLUS] += 1
                        multiSpeGenera[mulSpeGenTaxId] -= 1
                        newSurplus -=1
                    
                    # once full, stop looping
                    if sum(multiSpeGenera.values()) == 0:
                        break
                    
                    # once newSurplus is drained, stop looping
                    if newSurplus == 0:
                        break
            
            # add the new entries to the candidates
            for key in allNewEntries:
                candidates.append(allNewEntries[key])

        # if there are more missing genera than the surplus allows for
        else:
            # then randomly select some missing genera as candidates
            randomGenera = random.choices(list(missingGenera), k=newSurplus)
            ranGen:Taxonomy
            for ranGen in randomGenera:
                # one species per genus
                entry = ranGen.__makeIngroupCandidateEntry(1, False)
                candidates.append(entry)

        # we just made a genus-deep candidate list; redistribute them!
        Taxonomy.__redistributeSurplusToCandidates(ingroupSpecies,
                                                ingroupHigherTaxa, candidates)

        # ingroupSpecies was reassigned (near lines 543-546). As a result,
        # ingroupSpecies in the stack frame above has not been changed.
        # Instead, the ingroupSpecies of this stack frame has simply shifted 
        # its 'pointer' to a new list. We will have to return ingroupSpecies to
        # compensate. The other inputs were only modified, not reassigned, so 
        # they do not need to be returned.
        return ingroupSpecies




    ### MISCELLANEOUS METHODS
    def makeNewickStr(self, maxDepth:TaxRank=SPECIES) -> str:
        """ newickStr:
                Accepts a TaxRank (or equivalent str) as optional input. Calls
                a recursive helper function on the root of the object to gener-
                ate a textual representation of the full Taxonomy in Newick fo-
                rmat. Returns the string.
        """
        # constants
        SEP_CHAR = ','
        END_CHAR = ';'
        GREP_FIND = "\(\)"
        GREP_REPL = ''

        # convert max depth to a rank
        if type(maxDepth) is str:
            maxDepth = TaxRank(maxDepth)

        # start at root
        root = self.getRoot()

        # call recursive helper function on the root of the object
        nwkStr = root.__taxonomyToNewickHelper(maxDepth)

        # replace the trailing comma with a semi-colon
        if nwkStr[-1] == SEP_CHAR:
            nwkStr = nwkStr[:-1]
            nwkStr += END_CHAR

        # remove empty nodes and return
        return re.sub(GREP_FIND, GREP_REPL, nwkStr)


    def __taxonomyToNewickHelper(self, maxDepth:TaxRank) -> str:
        """ taxonomyToNewickHelper:
                Accepts a TaxRank object as input. Recursively constructs a te-
                xtual representation of the calling object in Newick format as
                a string. Returns the string.
        """
        # constants
        OPEN = '('
        CLOSE = ')'
        COMMA = ','
        QUOTE = '"'
        INVALID_STR = ' (invalid name)'
        GREP_FIND = r'\"'
        GREP_REPL = r''

        nwkStr = ''
        # if taxO is greater than maxDepth, then it is an internal node
        if self.rank > maxDepth:
            # open the node
            nwkStr += OPEN

            # recurse on all the children
            child:Taxonomy
            for child in self.getChildren(set):
                nwkStr += child.__taxonomyToNewickHelper(maxDepth)
            
            # remove the trailing comma
            if nwkStr[-1] == COMMA:
                nwkStr = nwkStr[:-1]
            
            # close the node
            nwkStr += CLOSE

        # get the scientific name
        name = self.sciName

        # if it is an invalid name, then remove quotation marks
        if INVALID_STR in name:
            name = re.sub(GREP_FIND, GREP_REPL, name)
        
        # add quotation marks around the full name followed by a comma
        name = QUOTE + name + QUOTE + COMMA

        # append and return
        nwkStr += name

        return nwkStr


