# Author: Joseph S. Wirth

import re
from PHANTASM.utilities import parseCsv


def makeLpsnD(genSpeCsv:str, genFamCsv:str, \
               phyClsOrdFamCsv:str, phyClsOrdFamSynonymsCsv:str) -> dict:
    """ importLPSN:
            Accepts the paths to four CSVs as input: the gss file from LPSN, a
            custom built genus-family csv, a custom built phylum-class-order-
            family csv, and a custom built synonyms table for the phyla, class-
            es, orders, and families. Imports the tables as parsed lists, and
            then relies on helper functions to convert them into their specific
            dictionaries. Adds these to an umbrella dictionary and returns it.
            
            For each rank between species and phylum, two dictionaries are 
            created:
                * str(rank)             (eg. 'species')
                * str(rank) + '_syn''   (eg. 'species_syn')
            
            The '[rank]_syn' dictionaries are keyed by LPSN's synonymous names
            with the preferred name as the values. The '[rank]' dictionaries are
            keyed by the preferred name and their values are dictionaries that
            are keyed by the following strings (extra fields allowed):
                * 'type'
                * 'status'
                * 'parent'
            Some examples follow below:
            
            Synonym dictionaries for all ranks follow this model behavior
            >>> lpsnD['genus_syn']['Silicibacter']
            >>> 'Ruegeria'

            Rank dictionaries for higher ranks have their typeMaterial as str
            >>> lpsnD['genus']['Ruegeria']
            >>> {'status': 'correct name',
                 'type': 'Ruegeria atlantica',
                 'parent': 'Rhodobacteraceae'}
            
            The species rank dictionary has its typeMaterial field as a list of strains
            >>> lpsnD['species']['Ruegeria pomeroyi']
            >>> {'type': ['ATCC 700808', 'DSM 15171', 'DSS-3', 'ATCC700808', 'DSM15171'],
                 'status': 'correct name',
                 'parent': 'Ruegeria'}
    """
    # read csv into parsed lists of rows
    genSpeData = parseCsv(genSpeCsv)
    genFamData = parseCsv(genFamCsv)
    phyClsOrdFamData = parseCsv(phyClsOrdFamCsv)
    synonymsData = parseCsv(phyClsOrdFamSynonymsCsv)

    # convert parsed data into individual dictionaries
    genusD, genSynD, speciesD, speSynD = __importGenusSpeciesData(genSpeData,\
                                                                  genFamData)
    phylumD, classD, orderD, familyD = __importHigherRankData(phyClsOrdFamData)
    phySynD, clsSynD, ordSynD, famSynD = __importSynonymsData(synonymsData)

    # replace any type genera in familyD that are synonyms
    __replaceSynonymTypeGenusWithPreferredName(familyD, genSynD)
    __replaceSynonymParentGenusWithPreferredName(speciesD, genSynD)

    # remove any synonyms whose preferred name cannot be successfully looked up
    __dropBadLpsnEntriesFromSynonymsDict(phySynD, phylumD)
    __dropBadLpsnEntriesFromSynonymsDict(clsSynD, classD)
    __dropBadLpsnEntriesFromSynonymsDict(ordSynD, orderD)
    __dropBadLpsnEntriesFromSynonymsDict(famSynD, familyD)
    __dropBadLpsnEntriesFromSynonymsDict(genSynD, genusD)
    __dropBadLpsnEntriesFromSynonymsDict(speSynD, speciesD)

    # remove any taxa from higher ranks (genus and above) without children
    __removeTaxaWithoutChildren(genusD, genSynD, speciesD)
    __removeTaxaWithoutChildren(familyD, famSynD, genusD)
    __removeTaxaWithoutChildren(orderD, ordSynD, familyD)
    __removeTaxaWithoutChildren(classD, clsSynD, orderD)
    __removeTaxaWithoutChildren(phylumD, phySynD, classD)

    # taxa higher than family may list a genus as the type; this needs to be fixed
    __fixHigherRankTypeMaterial(genusD, genSynD, familyD, famSynD, orderD, ordSynD, classD, clsSynD, phylumD)

    # make final dictionary
    lpsnD = {}
    lpsnD['phylum'] = phylumD
    lpsnD['class'] = classD
    lpsnD['order'] = orderD
    lpsnD['family'] = familyD
    lpsnD['genus'] = genusD
    lpsnD['species'] = speciesD
    lpsnD['phylum_syn'] = phySynD
    lpsnD['class_syn'] = clsSynD
    lpsnD['order_syn'] = ordSynD
    lpsnD['family_syn'] = famSynD
    lpsnD['genus_syn'] = genSynD
    lpsnD['species_syn'] = speSynD

    return lpsnD


def __importGenusSpeciesData(genSpeParsed:list, genFamParsed:list, \
                            headers:bool=True) -> tuple:
    """ importGenusSpeciesData:
            Accepts the parsed content of two specific csv files and a boolean
            indicating whether or not the file contained headers. Assumes the 
            first csv file possesses the structure indicated by the indices 
            (int) in the constants list. The genFamParsed file is used to link
            the genera to their parental families. Returns the following dict-
            ionaries in their final forms:
                * genus
                * genus_syn
                * species
                * species_syn
    """
    # constants for genSpeParsed
    GENUS = 0
    SPECIES = 1
    SUBSPECIES = 2
    REFERENCE = 3
    STATUS = 4
    AUTHORS = 5
    URL = 6
    RISK_GRP = 7
    NOM_TYPE = 8
    REC_NO = 9
    REC_LNK = 10
    BAD_NOM_TYPE = 'PENDING'
    SPLT_CHR = '; '
    CORRECT = 'correct name'
    SYNONYM = ['synonym', 'misspelling']
    STAT_IDX = -1

    # check for headers
    if headers:
        startRow = 1
    else:
        startRow = 0
    
    # create dictionaries for easy data access (first row == headers)
    allFullRecords = {}
    generaD = {}
    genSynD = {}
    speciesD = {}
    speSynD = {}

    for row in genSpeParsed[startRow:]:
        # make a dictionary containing all of the field info
        fullRec = {}
        fullRec['genus'] = row[GENUS]
        fullRec['species'] = row[SPECIES]
        fullRec['subspecies'] = row[SUBSPECIES]
        fullRec['reference'] = row[REFERENCE]
        fullRec['status'] = row[STATUS].split(SPLT_CHR)
        fullRec['authors'] = row[AUTHORS]
        fullRec['url'] = row[URL]
        fullRec['risk_grp'] = row[RISK_GRP]
        fullRec['nom_type'] = row[NOM_TYPE]
        fullRec['rec_no'] = row[REC_NO]
        fullRec['rec_lnk'] = row[REC_LNK]

        # add the entry to the dictionary
        allFullRecords[row[REC_NO]] = fullRec

        stat = fullRec['status'][STAT_IDX]

        # populate a separate dictionary to hold the genus info
        if row[SPECIES] == "":
            if row[NOM_TYPE] != BAD_NOM_TYPE:
                # make the entry
                genusEntry = {}
                key = row[GENUS]

                # save synonyms in a separate dictionary
                if stat in SYNONYM:
                    genSynD[key] = row[REC_LNK]

                # save the correct names only
                elif stat == CORRECT:
                    genusEntry['lpsn_no'] = row[REC_NO]
                    genusEntry['nom_type'] = row[NOM_TYPE]
                    genusEntry['status'] = stat
                    generaD[key] = genusEntry

        # populate a separate dictionary to hold the species info
        elif row[SUBSPECIES] == "":
            speciesEntry = {}

            # make the key the binomial name
            key = row[GENUS] + ' ' + row[SPECIES]

            # default override to off
            addAnyway = False

            # save synonyms separately
            if stat in SYNONYM:
                # if LPSN has not defined a synonym's referrant
                if row[REC_LNK] == '':
                    # then add the synonym to the rank dictionary instead
                    addAnyway = True
                
                # else if LPSN has defined the referrant
                else:
                    # then save the record link as the value
                    speSynD[key] = row[REC_LNK]

            # save the correct names only (or if the override is active)
            elif stat == CORRECT or addAnyway:
                speciesEntry['lpsn_no'] = row[REC_NO]
                speciesEntry['type'] = row[NOM_TYPE].split(SPLT_CHR)
                speciesEntry['status'] = stat

                speciesD[key] = speciesEntry

                # reset the override boolean
                addAnyway = False
    
    # add type species names to the genus dictionary
    __addTypeToGenera(generaD, allFullRecords, speciesD)

    # add the parents to the genus and species dictionaries
    __addParentsToGenera(generaD, genFamParsed, headers=headers)
    __addParentsToSpecies(speciesD)

    # fix type strain information
    __fixTypeStrains(speciesD)

    # import species referenced but missing --- THIS MAY BECOME OBSOLETE
    __importMissingSpecies(speciesD, speSynD, genFamParsed, headers=headers)

    # replace rec_no in synonym dictionaries with the full name
    __addPreferredNamesToSynonyms(genSynD, allFullRecords)
    __addPreferredNamesToSynonyms(speSynD, allFullRecords)

    # remove any residual homonyms in the synonym table
    # (eg. Pararhodobacter marinus)
    __removeHomonymsFromSynonymsDict(genSynD, generaD)
    __removeHomonymsFromSynonymsDict(speSynD, speciesD)

    # replace synonyms in the type material fields with their preferred names
    __replaceTypeSynonyms(speSynD, generaD)

    # remove extra fields that are no longer needed
    __removeExtraDataFromGenusSpecies(generaD, speciesD)

    #return allFullRecords
    return generaD, genSynD, speciesD, speSynD


def __addTypeToGenera(genRecs:dict, allRecs:dict, speRecs:dict) -> None:
    """ addTypeToGenera:
            Accepts the genus dictionary and a dictionary containing all the
            records for all of the genera and species as inputs. Identifies the
            type species for each genus and adds its name to the 'type' field 
            of the genus dictionary. Does not return.
    """
    # for each genus
    for genus in genRecs.keys():
        # get the record num of the nomenclatural type
        typeRecNo = genRecs[genus]['nom_type']

        # extract the record for the type species
        try:
            typeRec = allRecs[typeRecNo]
        
        # if that isn't possible, then just pick the first species available
        except:
            for species in speRecs.keys():
                if genus in species:
                    typeRec = allRecs[speRecs[species]['lpsn_no']]
                    break

        # use the record to make the binomial name and add it to the genus dict
        typeName = typeRec['genus'] + ' ' + typeRec['species']
        genRecs[genus]['type'] = typeName


def __addParentsToGenera(genusD:dict, genFamParsed:list, headers:bool=True) \
                        -> None:
    """ addParentsToGenera:
            Accepts the genus dictionary, the parsed genus-family csv, and a
            boolean indicating if the genus-family csv contains headers. Uses
            the genus-family table to determine the family for each genus. Adds
            the family name to the genus dictionary. Does not return.
    """
    # constants
    GENUS_IDX    = 0
    TYPE_SPE_IDX = 1
    PARENT_IDX   = 2
    LPSN_NUM_IDX = 3
    NOM_TYPE_IDX = 4
    STATUS_IDX   = 5
    SYNONYM = ['synonym', 'misspelling']
    ERROR_MSG = "No parent designated for the genus "

    # check for headers
    if headers:
        startRow = 1
    else:
        startRow = 0

    # for each row in genFamParsed
    for row in genFamParsed[startRow:]:
        # get the genus and its parental family from the table
        genus   = row[GENUS_IDX]
        family  = row[PARENT_IDX]
        typeSpe = row[TYPE_SPE_IDX]
        status  = row[STATUS_IDX]

        # look up (or add) the corresponding genus in the dictionary
        if genus in genusD.keys():
            # add the family information
            genusD[genus]['parent'] = family

        elif status not in SYNONYM:
            # make a new entry for the missing genus
            entry = {}
            entry['type'] = typeSpe
            entry['parent'] = family
            entry['status'] = status
            genusD[genus] = entry
    
    # make sure that all genera have a parent!
    missingParents = set()
    for genus in genusD.keys():
        if 'parent' not in genusD[genus].keys():
            missingParents.add(genus)
    
    # alert user if their data sets are not fully linked.
    if len(missingParents) > 0:
        parStr = ''
        for missing in missingParents:
            parStr += "'" + missing + "', "
        parStr = parStr[:-2]

        raise BaseException(ERROR_MSG + parStr)


def __importMissingSpecies(speciesD:dict, speSynD:dict, genFamParsed:list, \
                                                    headers:bool=True) -> None:
    """ importMissingSpecies:
            Accepts the species dictionary, the species synonym dictionary, the
            parsed genus-family csv (list), and a boolean as inputs. Searches
            the list for any species that are not currently present in either
            of the dictionaries. Adds any missing species to the species dicti-
            onary. Modifies the dictionary but does not return.
    """
    # constants
    GENUS_IDX    = 0
    TYPE_SPE_IDX = 1
    PARENT_IDX   = 2
    LPSN_NUM_IDX = 3
    NOM_TYPE_IDX = 4
    STATUS_IDX   = 5
    SYNONYM = ['synonym', 'misspelling']
    NEW_SPE_STATUS = "status unknown; placeholder taxon"
    NEW_SPE_TYPEMAT = []

    # check for headers
    if headers:
        startRow = 1
    else:
        startRow = 0

    # for each row in genFamParsed
    for row in genFamParsed[startRow:]:
        # get the type species and its parent
        typeSpe = row[TYPE_SPE_IDX]
        genus   = row[GENUS_IDX]

        # make a new entry if it is not found in the species dict
        if typeSpe not in speciesD.keys() and typeSpe not in speSynD.keys() \
                          and typeSpe != "" and row[STATUS_IDX] not in SYNONYM:
            # assign values for parent, type, and status
            entry = dict()
            entry['parent'] = genus
            entry['type'] = NEW_SPE_TYPEMAT
            entry['status'] = NEW_SPE_STATUS

            speciesD[typeSpe] = entry


def __addParentsToSpecies(speciesD:dict) -> None:
    """ addParentsToSpecies:
            Accepts the species dictionary as input. For each species, extracts
            the genus name from the binomial name and saves it in the 'parent'
            field. Does not return.
    """
    # constants
    GREP_FIND = r' \S+$'
    GREP_REPL = r''

    # for each species
    for species in speciesD.keys():
        # extract the genus name
        genus = re.sub(GREP_FIND, GREP_REPL, species)

        # save the genus name as parent
        speciesD[species]['parent'] = genus


def __fixTypeStrains(speciesD:dict) -> None:
    """ fixTypeStrains:
            Accepts the species dictionary as input. Removes 'strain ' from the
            front of any strain name. Adds alternate spellings of applicable
            strains such that spaces have been ommitted (improves accuracy when
            mapping NCBI assembly data later). Does not return.
    """
    # constants
    TYPE = 'type'
    GREP_FIND = r'^strain '
    GREP_REPL = r''
    SPACE = ' '

    # for each species
    for spe in speciesD.keys():
        # initialize alt-spelling strain list
        altSpellings = list()

        # for each type strain
        for i in range(len(speciesD[spe][TYPE])):
            strain = speciesD[spe][TYPE][i]

            # remove 'strain ' from the beginning of any strain name
            strain = re.sub(GREP_FIND, GREP_REPL, strain)
            speciesD[spe][TYPE][i] = strain

            # make duplicate entry where spaces are ommitted as alt. spellings
            if SPACE in strain:
                altSpell = re.sub(SPACE, GREP_REPL, strain)
                altSpellings.append(altSpell)
        
        # add the alternate spellings to the list of strains
        speciesD[spe][TYPE] += altSpellings


def __addPreferredNamesToSynonyms(synD:dict, fullRecs:dict) -> None:
    """ addPreferredNamesToSynonyms:
            Accepts a synonym dictionary and all of the genus-species records 
            as input. Assumes the value of the synonym dictionary is an LPSN
            record number and uses it to look up the preferred name. Replaces
            the value with the preferred name. Does not return.
    """
    # initialize the set of keys that need to be removed from the dictionary
    badKeys = set()

    # for each synonym
    for key in synD.keys():
        # get the record for the preferred name
        prefRec = fullRecs[synD[key]]  # synD's values are lpsn_record_links

        # make the binomial name
        genName = prefRec['genus']
        speName = prefRec['species']
    
        # if the species name is empty, then the synonym is a genus
        if speName == '':
            # replace the lpsn_record_link with the preferred name
            if key != genName:
                synD[key] = genName
            
            # mark synonyms for removal if their replacement has the same name
            else:
                badKeys.add(key)
        
        # otherwise, the synonym is a species name
        else:
            # replace the lpsn_record_link with the preferred name
            sciName = genName + " " + speName
            if key != sciName:
                synD[key] = sciName
            
            # mark synonyms for removal if their replacement has the same name
            else:
                badKeys.add(key)
    
    # remove any unneeded keys from the synonym dictionary
    for key in badKeys:
        synD.pop(key)


def __removeHomonymsFromSynonymsDict(synD:dict, rankD:dict) -> None:
    """ removeHomonymsFromSynonymsDict:
            Accepts a synonym dictionary and a rank dictionary as input. Finds
            and removes any entries in the synonym dictionary whose key is also
            a key in the rank dictionary. Does not return.
    """
    # for each name in the rank dictionary
    for sciName in rankD.keys():
        # remove any homonyms in the synonym dictionary
        if sciName in synD.keys():
            synD.pop(sciName)


def __replaceSynonymTypeGenusWithPreferredName(familyD:dict, genSynD:dict) \
                                                                       -> None:
    """ replaceSynonymTypeGenusWithPreferredName:
            Accepts the family dictionary and the genus synonym dictionary as
            inputs. Searches the family dictionary for any instances where the
            type genus is a synonym and replaces the synonym with the preferred
            name. Modifies the family dictionary but does not return.
    """
    # for each family in the dictionary
    for fam in familyD.keys():
        # get the type genus
        typeGen = familyD[fam]['type']

        # if the type genus is a synonym, then replace with the preferred name
        if typeGen in genSynD.keys():
            familyD[fam]['type'] = genSynD[typeGen]


def __replaceSynonymParentGenusWithPreferredName(speciesD:dict, genSynD:dict) \
                                                                       -> None:
    """ replaceSynonymParentGenusWithPreferredName:
            Accepts the species dictionary and the genus synonym dictionary as
            inputs. Searches the species dictionary for any instances where the
            parental genus is a synonym and replaces the synonym with the pref-
            erred name. Modifies the species dictionary but does not return.
    """
    # for each species in the dictionary
    for spe in speciesD.keys():
        # get the parental genus
        parentGen = speciesD[spe]['parent']

        # if the parental genus is a synonym, then replace with preferred name
        if parentGen in genSynD.keys():
            speciesD[spe]['parent'] = genSynD[parentGen]


def __dropBadLpsnEntriesFromSynonymsDict(synD:dict, rankD:dict) -> None:
    """ dropBadLpsnEntriesFromSynonymsDict:
            Accepts a synonym dictionary and a rank dictionary as input. Finds
            and drops any synonyms whose preferred name is missing in the rank
            dictionary. Does not return.
    """
    # initialize a set of keys to drop
    keysToDrop = set()

    # stage synonyms for removal if their preferred name isn't in the rank dict
    for synKey in synD.keys():
        if synD[synKey] not in rankD.keys():
            keysToDrop.add(synKey)
    
    # drop the marked keys
    for key in keysToDrop:
        synD.pop(key)


def __replaceTypeSynonyms(speSynD:dict, genusD:dict) -> None:
    """ replaceTypeSynonyms:
            Accepts two dictionaries as inputs. Uses the species synonym dicti-
            onary to check the type material field of the genus dictionary for
            the presence of synonyms. Uses the species synonym dictionary to
            replace any synonyms in the genus dictionary's type material field.
            Does not return; modifies the genus dictionary.
    """
    # for each scientific genus name
    for sciName in genusD.keys():
        # get the type material
        typeMat = genusD[sciName]['type']

        # if the type is a synonym, then replace it
        if typeMat in speSynD.keys():
            genusD[sciName]['type'] = speSynD[typeMat]


def __removeExtraDataFromGenusSpecies(genusD:dict, speciesD:dict) -> None:
    """ removeExtraDataFromGenusSpecies:
            Accepts two dictionaries as inputs. Removes all extra fields from
            the genus/species dictionaries keeping only those keys that are
            shared between all of the dictionaries. Does not return.
    """
    # for each genus
    for genName in genusD.keys():
        __removeExtraDataHelper(genName, genusD)
        
    
    for speName in speciesD.keys():
        __removeExtraDataHelper(speName, speciesD)
        

def __removeExtraDataHelper(name:str, rankD:dict) -> None:
    """ removeExtraDataHelper:
            Accepts an scientific name and its rank dictionary as input. Looks
            for and removes extra keys from the entry in their rank dictionary.
            Does not return.
    """
    # constant
    KEYS_TO_KEEP = {'type', 'status', 'parent'}

    # obtain a list of keys that should be removed
    keysToRemove = set()
    for key in rankD[name].keys():
        if key not in KEYS_TO_KEEP:
            keysToRemove.add(key)
    
    # remove the marked keys
    for badKey in keysToRemove:
        rankD[name].pop(badKey)


def __importHigherRankData(phyClsOrdFamData:list, headers:bool=True) -> tuple:
    """ importHigherRankData:
            Accepts a parsed phylum-class-order-family csv file and a boolean
            indicating whether or not the file contains headers as input. Assu-
            mes the csv file uses the structure as indicated by the indices in
            the constants list. Converts the parsed csv into rank dictionaries.
            Returns the following rank dictionaries as tuple (listed in order)
            in their final forms:
                * phylum
                * class
                * order
                * family
    """
    # constants
    NAME = 0
    TYPE_NAME = 1
    PARENT = 2
    STATUS = 3
    RANK = 4

    # check headers
    if headers:
        startRow = 1
    else:
        startRow = 0

    # use parsed data to populate dictionaries
    phylumD = {}
    classD = {}
    orderD = {}
    familyD = {}
    for row in phyClsOrdFamData[startRow:]:
        # parse fields into convenient variable names
        key = row[NAME]
        typeName = row[TYPE_NAME]
        parent = row[PARENT]
        status = row[STATUS]
        rank = row[RANK]

        # convert missing type material to None (happens for placeholders)
        if typeName == '':
            typeName = None

        # make the entry
        entry = {}
        entry['type'] = typeName
        entry['parent'] = parent
        entry['status'] = status

        # populate phylum dictionary
        if rank == 'phylum':
            phylumD[key] = entry

        # populate class dictionary
        elif rank == 'class':
            classD[key] = entry

        # populate  order dictionary
        elif rank == 'order':
            orderD[key] = entry
        
        # populate family dictionary
        elif rank == 'family':
            familyD[key] = entry
        
    return phylumD, classD, orderD, familyD


def __importSynonymsData(synonymsParsed:list, headers:bool=True) -> tuple:
    """ importSynonymsData:
            Accepts a parsed phylum-class-order-family synonyms csv file and a
            boolean indicating whether or not the file contains headers as inp-
            uts. Assumes the csv file has the structure indicated by the values 
            in the constants list. Converts the parsed file into the synonyms
            dictionaries. Returns the following dictionaries in their final fo-
            rms as a tuple (listed in order):
                * phylum_syn
                * class_syn
                * order_syn
                * family_syn
    """
    # constants
    SYNONYM = 0
    PREFERRED = 1
    RANK = 2

    # check for headers
    if headers:
        startRow = 1
    else:
        startRow = 0

    # use parsed data to populate dictionaries
    phySynD = {}
    clsSynD = {}
    ordSynD = {}
    famSynD = {}
    for row in synonymsParsed[startRow:]:
        key = row[SYNONYM]
        val = row[PREFERRED]
        rank = row[RANK]

        # populate phylum_syn dict
        if rank == 'phylum':
            phySynD[key] = val

        # populate class_syn dict
        if rank == 'class':
            clsSynD[key] = val

        # populate order_syn dict
        if rank == 'order':
            ordSynD[key] = val
        
        # populate family_syn dict
        if rank == 'family':
            famSynD[key] = val
    
    return phySynD, clsSynD, ordSynD, famSynD


def __removeTaxaWithoutChildren(parentRankD:dict, parentSynD:dict, \
                                childRankD:dict) -> None:
    """ removeTaxaWithoutChildren:
            Accepts a three dictionaries as input: a rank dictionary, the syno-
            nym dictionary for that rank, and the dictionary for the rank below
            as inputs. Looks for and removes for taxa that are not listed as a
            parent for one or more taxa in the rank below. Does not return.
    """
    # get a set of all the unique parents that are referenced by the children
    parentNames = set()
    for child in childRankD:
        parentNames.add(childRankD[child]['parent'])
    
    # for each parent, mark those without children for removal
    barrenParents = set()
    for parent in parentRankD:
        if parent not in parentNames:
            barrenParents.add(parent)

    # for each barrent parent, remove them it the rank dictionary
    for barPar in barrenParents:
        # remove from the rank dictionary
        parentRankD.pop(barPar)

    # mark synonyms for removal if they reference a barren parent
    badSynonyms = set()
    for syn in parentSynD:
        if parentSynD[syn] in barrenParents:
            badSynonyms.add(syn)

    # for each bad synonym, remove it from the synonym dictionary
    for badSyn in badSynonyms:
        parentSynD.pop(badSyn)


def __fixHigherRankTypeMaterial(genusD:dict, genSynD:dict, familyD:dict, famSynD:dict, orderD:dict, ordSynD:dict, classD:dict, clsSynD:dict, phylumD:dict) -> None:
    # start by replacing any synonyms in the type material
    __replaceSynonymsInTypeMaterial(familyD, 'family', genSynD, famSynD, ordSynD, clsSynD)
    __replaceSynonymsInTypeMaterial(orderD,  'order',  genSynD, famSynD, ordSynD, clsSynD)
    __replaceSynonymsInTypeMaterial(classD,  'class',  genSynD, famSynD, ordSynD, clsSynD)
    __replaceSynonymsInTypeMaterial(phylumD, 'phylum', genSynD, famSynD, ordSynD, clsSynD)

    # check orders first
    for ord in orderD.keys():
        # if the order's type is a genus
        oType = orderD[ord]['type']
        if oType in genusD.keys():
            # then replace it with that genus's parental family
            parentalFam = genusD[oType]['parent']
            orderD[ord]['type'] = parentalFam
    
    # check classes next
    for cls in classD.keys():
        # if the class's type is a genus
        cType = classD[cls]['type']
        if cType in genusD.keys():
            # then replace it with that genus's parental order
            parentalFam = genusD[cType]['parent']
            parentalOrd = familyD[parentalFam]['parent']
            classD[cls]['type'] = parentalOrd
    
    # check phyla last
    for phy in phylumD.keys():
        # if the phylum's type is a genus
        pType = phylumD[phy]['type']
        if pType in genusD.keys():
            # then replace it with that genus's parental class
            parentalFam = genusD[pType]['parent']
            parentalOrd = familyD[parentalFam]['parent']
            parentalCls = orderD[parentalOrd]['parent']
            phylumD[phy]['type'] = parentalCls
        
        # if the phylum's type is an order
        if pType in orderD.keys():
            parentalCls = orderD[pType]['parent']
            phylumD[phy]['type'] = parentalCls


def __replaceSynonymsInTypeMaterial(rankD:dict, rank:str, genSynD:dict, famSynD:dict, ordSynD:dict, clsSynD:dict):
    RANK_TO_NUM = {"family": 0,
                   "order":  1,
                   "class":  2,
                   "phylum": 3}

    for k in rankD.keys():
        typeName = rankD[k]['type']
        
        # check for genus synonyms
        if typeName in genSynD.keys():
            rankD[k]['type'] = genSynD[typeName]
        
        # check for family synonyms
        elif typeName in famSynD.keys() and RANK_TO_NUM[rank] > RANK_TO_NUM['family']:
            rankD[k]['type'] = famSynD[typeName]

        # check for order synonyms
        elif typeName in ordSynD.keys() and RANK_TO_NUM[rank] > RANK_TO_NUM['order']:
            rankD[k]['type'] = ordSynD[typeName]
        
        # check for class synonyms
        elif typeName in clsSynD.keys() and RANK_TO_NUM[rank] > RANK_TO_NUM['class']:
            rankD[k]['type'] = clsSynD[typeName]