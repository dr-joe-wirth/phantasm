import re
from PHANTASM.utilities import ncbiIdsFromSearchTerm, isFloat, parseCsv


def getTaxIdsFromRnaBlast(blastResultFile:str) -> list:
    """ getTaxIdsFromRnaBlast:
            Accepts the filename of a blastn file as input. Reads the file and
            extracts the NCBI taxonomy IDs as a list. Returns the list.
    """
    # constants
    ERR_MSG_1 = "The file "
    ERR_MSG_2 = " does not contain any valid blastn hits."

    # get subjects from file
    subjects = __readRnaBlastFile(blastResultFile)

    # extract taxonomy ids from hits
    taxids = __getUniqueTaxIdsFromSubjects(subjects)

    # handle incorrectly formatted subject entries
    __fixTaxIds(taxids)

    if len(taxids) == 0:
        raise Exception(ERR_MSG_1 + blastResultFile + ERR_MSG_2)
    
    return taxids


def __readRnaBlastFile(blastFile:str) -> list:
    """ readRnaBlastFile:
            Accepts the filename of a blastn file as input. Returns a list con-
            taining the strings found in the stitle column.
    """
    # constants
    PERC_ID_CUTOFF = 94
    PID_IDX = 2
    STITLE_IDX = 1

    # read the file as a list of lists
    blastParsed = parseCsv(blastFile, '\t')
    
    # extract the hits from the blast table
    subjects = set()
    for row in blastParsed:
        if isFloat(row[PID_IDX]):
            # ensure that the hit exceeds the percent id cut-off:
            pid = float(row[PID_IDX])
            if pid >= PERC_ID_CUTOFF:
                subjects.add(row[STITLE_IDX])
    
    return subjects


def __getUniqueTaxIdsFromSubjects(subjects:list) -> set:
    """ getUniqueTaxIdsFromSubjects:
            Accepts a list of subject names (stitle) as input. Extracts the NC-
            BI taxonomic IDs from the strings. Returns a set containing the un-
            ique IDs present in the list.
    """
    # initialize the list to hold unique taxids
    uids = set()

    # for each subject in the set
    for subj in subjects:
        # get the NCBI taxonomy id and add it to the set
        taxid = __getTaxIdFromSubject(subj)
        uids.add(taxid)
    
    return uids


def __getTaxIdFromSubject(subject:str) -> str:
    """ getTaxIdFromSubject:
            Accepts a subject name (stitle) string as input. Uses a regular 
            expression to extract the NCBI taxonomy ID from the string. Returns
            the extracted ID.
    """
    # constants
    RE_FIND = r'^.+taxon:(\d+)$'
    RE_REPLACE = r'\1'

    return re.sub(RE_FIND, RE_REPLACE, subject)


def __fixTaxIds(taxids:set) -> None:
    """ fixTaxIds:
            Accepts a set of NCBI Taxonomy ids as input. Looks for any ids that
            are not correctly formatted. Attempts to coerce their format and
            adds them to the set. Removes any ids that cannot be properly form-
            atted. Does not return.
    """
    # constants
    GREP_F = r'^(\S+ [^\|]+)\|.+$'
    GREP_R = r'\1'
    DATABASE = 'taxonomy'
    SEP = ' OR '
    SEP_LEN = len(SEP)
    ERR = 'The following search to NCBI failed:\n  database:    ' + DATABASE \
                                                          + '\n  search term: '

    # make a search term to find missing taxids and remove bad names from set
    searchTerm = ''
    badTaxIds = set()
    txid:str
    for txid in taxids:
        # if the entry is not numeric, then it likely has the species name
        if not txid.isnumeric():
            # save the taxid for removal (cannot change set size during iter)
            badTaxIds.add(txid)

            # attempt to get the species name from the string
            speName = re.sub(GREP_F, GREP_R, txid)

            # add the name to the search term
            searchTerm += speName + SEP
    
    # remove the trailing SEP from the search term
    searchTerm = searchTerm[:-SEP_LEN]
    
    # only do something if there is a problem to fix
    if len(badTaxIds) > 0:
        # remove the bad ids from the current set
        for bad in badTaxIds:
            taxids.remove(bad)
        # look up new taxids via NCBI
        try:
            newTaxIds = ncbiIdsFromSearchTerm(searchTerm, DATABASE)
        except:
            raise RuntimeError(ERR + searchTerm)
        
        # for each new taxid retrieved, add it to the set
        for newId in newTaxIds:
            taxids.add(newId)


