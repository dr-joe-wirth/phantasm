# Author: Joseph S. Wirth
# Last edit: September 27, 2022

import csv, ftplib, glob, gzip, math, os, re, shutil
from Bio import Entrez
from Bio.Entrez import Parser
from PHANTASM.taxonomy.TaxRank import TaxRank
from PHANTASM.Parameter import Parameters


def getParamO_1(email:str) -> Parameters:
    """ getParamO_1:
            Helper function to facilitate the creation of the first Parameters
            object. Accepts a valid email address as input. Returns the Parame-
            ters object.
    """
    # import relevant data from param.py and specify the anlysis directory
    from param import BLASTPLUS_DIR, MUSCLE_EXE, FASTTREE_EXE, IQTREE_EXE, \
                                     NUM_PROCESSORS, MAX_LEAVES, NUM_BOOTSTRAPS
    ANALYSIS_DIR = './initialAnalysis'

    # build the parameter object
    parameterO = Parameters(email,
                            ANALYSIS_DIR,
                            BLASTPLUS_DIR,
                            MUSCLE_EXE,
                            FASTTREE_EXE,
                            IQTREE_EXE,
                            NUM_PROCESSORS,
                            MAX_LEAVES,
                            NUM_BOOTSTRAPS)

    return parameterO


def getParamO_2(email:str) -> Parameters:
    """ getParamO_2:
            Helper function to facilitate the creation of the second Parameters
            object. Accepts a valid email address as input. Returns the Parame-
            ters object.
    """
    # import relevant data from param.py and specify the anlysis directory
    from param import BLASTPLUS_DIR, MUSCLE_EXE, FASTTREE_EXE, IQTREE_EXE, \
                                     NUM_PROCESSORS, MAX_LEAVES, NUM_BOOTSTRAPS
    ANALYSIS_DIR = './finalAnalysis'

    # build the parameter object
    parameterO = Parameters(email,
                            ANALYSIS_DIR,
                            BLASTPLUS_DIR,
                            MUSCLE_EXE,
                            FASTTREE_EXE,
                            IQTREE_EXE,
                            NUM_PROCESSORS,
                            MAX_LEAVES,
                            NUM_BOOTSTRAPS)
    
    return parameterO


def getTaxidsFromFile(taxidsFN:str) -> list[str]:
    """ getTaxidsFromFile:
            Accepts the filename of the taxids file. Expects the file to cont-
            ain one taxid per line. Returns a list of the taxids present.
    """
    # constants
    ERR_MSG_1 = "'" + taxidsFN + "' is empty."
    ERR_MSG_2 = "'" + taxidsFN + "' is improperly formatted."
    EOL = " ... "
    PRNT_1 = "Extracting taxids from the specified file" + EOL
    DONE = "Done.\n"
    
    # print status
    print(PRNT_1, end="", flush=True)

    # open the file
    fh = open(taxidsFN, "r")

    # go through each row
    outL = list()
    for row in fh:
        # drop the new line character if present
        if row[-1] == "\n":
            row = row[:-1]
        
        # append the extracted taxid to the list (no empty strings allowed)
        if len(row) > 1:
            outL.append(row)
    
    fh.close()

    # check that the file was not empty
    if len(outL) == 0:
        raise BaseException(ERR_MSG_1)
    
    for txid in outL:
        # all ids should be coerceable to int
        try:
            int(txid)
        except:
            raise BaseException(ERR_MSG_2)

    # print status
    print(DONE)

    return outL


def validEmailAddress(email:str) -> bool:
    """ validEmailAddress:
            Accepts an email address (str) as input. Ensures that it matches an
            expected pattern as defined by two general regular expressions. Re-
            turns a boolean indicating if the email is valid.
    """
    GREP_FIND_1 = r'\s'
    GREP_FIND_2 = r'^[^@]+@[^@]+\.[^@]+$'
    
    # if whitespace present, then invalid
    if re.search(GREP_FIND_1, email):
        return False
    
    # one '@' followed by at least one '.'
    elif re.search(GREP_FIND_2, email):
        return True
    
    # if here, then false
    return False


def checkEntrezEmail(email:str) -> None:
    """ checkEntrezEmail:
            Accepts an email address as input. Checks if Entrez.email has been
            set. If not, then it attempts to set Entrez.email to the provided
            email address. If not possible, then it raises an exception. Does 
            not return.
    """
    # make sure Entrez.email can be set
    if Entrez.email is None:
        if email is not None:
            Entrez.email = email
        else:
            raise BaseException("Entrez.email has not been set.")


def ncbiEfetchById(id, database:str, mode:str='text', rettype:str='xml', \
                                             retmode:str=None, email:str=None):
    """ ncbiEfetchById:
            Accepts an NCBI id (or a list of ids) and an NCBI database as 
            input. Uses NCBI's Efetch tool to retrieve full records from NCBI
            for the provided id(s). Returns the retrieved records.
    """
    # make sure Entrez.email has been set
    checkEntrezEmail(email)

    # retrieve the result from NCBI and return it
    handle = Entrez.efetch(db=database, id=id, mode=mode, rettype=rettype, retmode=retmode)
    result = Entrez.read(handle)
    handle.close()

    return result


def ncbiIdsFromSearchTerm(searchStr:str, database:str, retmax:int=0, \
                                                               email:str=None):
    """ ncbiIdsFromSearchTerm:
            Accepts a search term and an NCBI database as input. Uses NCBI's 
            Esearch tool to retrieve a list of ids that were found. By default,
            retrieves all possible records, but modifying the 'retmax' para-
            meter to a value greater than zero will limit the number of NCBI 
            ids that will be retrieved.
    """
    # constants
    BIG_NUMBER = 10000

    # make sure Entrez.email has been set
    checkEntrezEmail(email)

    # if retmax is less than 1, then retrieve all results
    if retmax < 1:
        retmax = BIG_NUMBER
        retrieveAll = True
    else:
        retrieveAll = False

    # retrieve result from NCBI
    handle = Entrez.esearch(db=database, term=searchStr, retmax=retmax)
    result = Entrez.read(handle)
    handle.close()

    # grab all results if requested (default)
    if retrieveAll:
        # determine the number of records found
        numFound = int(result['Count'])

        # repeat search if the amount retrieved is less than the amount found
        if numFound > retmax:
            retmax = numFound

            # retrieve result from NCBI
            handle = Entrez.esearch(db=database, term=searchStr, retmax=retmax)
            result = Entrez.read(handle)
            handle.close()
    
    # return the list of ids that were found
    return result['IdList']

    
def ncbiSummaryFromIdList(idList:list, database:str, email:str=None):
    """ ncbiSummaryFromIdList:
            Accepts a list of NCBI ids and an NCBI database as input. Uses 
            NCBI's Esummary tool to retrieve the summaries of the records for 
            the provided ids. Returns the search result.
    """
    # constants
    SEP_CHAR  = ', '
    SEP_LEN   = len(SEP_CHAR)
    RET_MAX   = 10000    # this is NCBI's default; cannot be exceeded
    RESULT_K1 = 'DocumentSummarySet'
    RESULT_K2 = 'DocumentSummary'
    VALIDATE  = False   # Bio.Entrez may throw a tantrum if modified

    # check that Entrez.email has been set
    checkEntrezEmail(email)

    # if not a list, assume it is an int or a string
    if type(idList) is str:
        try:
            int(idList)
        except ValueError:
            raise ValueError("invalid input: expected a list of ids or a single id")

    if type(idList) in [int, str]:
        queryL = [str(idList)]

    # otherwise assume it is a list (should work with stupid Entrez formats, too)
    else:
        # convert the list into a comma-separated string
        queryL = list()
        idStr  = str()
        count  = 0
        for id in idList:
            # make sure the ids are valid
            try:
                int(id)
            except ValueError:
                raise ValueError("invalid id present in the input list: " + str(id))

            # separate each id with a comma
            idStr += str(id) + SEP_CHAR
            count += 1

            # make sure each string contains no more than max allowed ids
            if count % RET_MAX == 0:
                # remove the trailing SEP_CHAR from the string
                idStr = idStr[:-SEP_LEN]

                # add the query to the list and begin making the next query
                queryL.append(idStr)
                idStr = str()
            
            # once all ids have been added to a query string
            elif count == len(idList):
                # remove the trailing SEP_CHAR from the string
                idStr = idStr[:-SEP_LEN]
                # add the query to the list
                queryL.append(idStr)

    # search NCBI with the first query
    handle = Entrez.esummary(id=queryL[0], db=database)
    result = Entrez.read(handle, validate=VALIDATE)
    handle.close()

    # if there are multiple queries, then keep querying NCBI until done
    for i in range(1,len(queryL)):
        query = queryL[i]

        handle = Entrez.esummary(id=query, db=database)
        nextResult = Entrez.read(handle, validate=VALIDATE)
        handle.close()

        # something I don't wuote understand is happening here.
        # sometimes, Entrez returns a list, other times its a dict
        # (technically, they are Entrez's stupid format, and not list or dict)
        # the problem may be database-specific:
        ### I think taxonomy returns list and assembly returns dict
        ### it's not worth resolving since the try-except works
        try:
            # append the results if they are returned as lists
            result.extend(nextResult)
        except:
            # append the results if they are returned as dicts
            result[RESULT_K1][RESULT_K2] += nextResult[RESULT_K1][RESULT_K2]

    return result


def ncbiELinkFromIdList(idL:list, dbFrom:str, dbTo:str) -> list:
    """ ncbiELinkFromIdList:
            Accepts a list of uids, a string indicating which database the uids
            correspond to, and a string indicating which database to link to as
            inputs. Retrieves and returns a list of the link records for the q-
            uery.
    """
    # max num ids per elink request is roughly 200 
    # (https://www.ncbi.nlm.nih.gov/books/NBK25499/)
    MAX_NUM_IDS = 200
    
    # calculate the number of times we will need to loop to not exceed the max 
    numIters = len(idL) / MAX_NUM_IDS
    numIters = math.ceil(numIters)

    # get a list of nuccore ids for the protein ids that are represented
    result = list()
    for iter in range(numIters):
        # determine which portion of the list to slice
        start = iter * MAX_NUM_IDS
        end = start + MAX_NUM_IDS

        # get a handle to the elink result
        handle = Entrez.elink(id=idL[start:end], dbfrom=dbFrom, db=dbTo)

        # first time through, save the result
        result += Entrez.read(handle)
        
        # close the handle
        handle.close()

    return result


def extractIdsFromELink(ncbiELink:list) -> dict:
    """ extractIdsFromELink:
            Accepts a list of ncbiELink results as input. Extracts the uids for
            both the original database and the linked database. Returns a dict-
            ionary where the linked database uids are keys and the correspondi-
            ng original database uids are the values.
    """
    # go through each link record
    idsD = dict()
    for link in ncbiELink:
        # if one or more link is present
        if link['LinkSetDb'] != []:
            # get the the id for the original database
            idFrom = link['IdList'][0]

            # get the id for the linked database
            idTo = link['LinkSetDb'][0]['Link'][0]['Id']

            # make a dictionary linking the dbTo to dbFrom
            idsD[idTo] = idFrom
    
    return idsD


def getParentalTaxIds(taxids, targetRank) -> list:
    """ getParentalTaxId:
            Accepts a taxonomy id (or a list of ids) and a target taxonomic
            rank as input. Returns a list of dictionaries with the keys 'name'
            and 'txid'.
    """
    # constants
    DATABASE = "taxonomy"
    
    # ensure the target rank is valid (let TaxRank constructor check)
    if type(targetRank) is str:
        targetRank = TaxRank(targetRank)
    
    # get the string of the input rank
    targetRank = str(targetRank)
    
    # import data from the taxonomy database
    taxRecords = ncbiEfetchById(taxids, DATABASE)
    
    # for each record retrieved
    parentTaxid = []
    for tax in taxRecords:
        # get the lineage
        allLineages = tax['LineageEx']

        # go through the lineage until target rank is found
        for lineage in allLineages:
            if lineage['Rank'] == targetRank:
                newEntry = {'txid': lineage['TaxId'], 
                            'name': lineage['ScientificName']}
                if newEntry not in parentTaxid:
                    parentTaxid.append(newEntry)
    
    return parentTaxid


def coerceEntrezToNormal(entrezElement):
    """ coerceEntrezToNormal:
            Accepts an object whose type is a sub-type of the Bio.Entrez class
            as input. Creates an equivalent object with the corresponding "nor-
            mal" type. Returns the "normal" object.
    """
    # determine the type
    curType = type(entrezElement)
    newType = curType

    # convert to string
    if curType is Parser.StringElement:
        newType = str
    
    # convert to list
    elif curType is Parser.ListElement:
        newType = list
    
    # convert to dictionary
    elif curType is Parser.DictionaryElement:
        newType = dict
    
    # if checking a list or dict, then recurse through the items within
    if newType is list:
        for i in range(len(entrezElement)):
            entrezElement[i] = coerceEntrezToNormal(entrezElement[i])
    elif newType is dict:
        for key in entrezElement:
            entrezElement[key] = coerceEntrezToNormal(entrezElement[key])

    # return the modified object
    return newType(entrezElement)


def parseCsv(csvFN:str, delim:str=',') -> list:
    """ parseCsv:
            Accepts a csv filename and delimiter (default = ,) as input. Opens
            the file, converts each row to a list of columns, and then adds 
            that list to a list of rows (list of lists). Returns the parsed da-
            ta.
    """
    handle = open(csvFN, 'r')
    parsed = list(csv.reader(handle, delimiter=delim))
    handle.close()

    return parsed


def isFloat(string:str) -> bool:
    """ isFloat:
            Accepts a string as input. Returns a boolean indicating whether or 
            not the string can be safely converted to a float.
    """
    try:
        float(string)
        return True
    except ValueError:
        return False


def changeFileExtension(inFile:str, newExtension:str) -> str:
    """ changeFileExtension:
            Accepts two strings: a filename and the new extension to add to it.
            Returns the modified filename with the new extension added.
    """
    # remove the original file extension
    out = removeFileExtension(inFile)

    # make sure newExtension will follow a period
    if newExtension[0] != '.' and out[-1] != '.':
        out += '.'

    # add the extension and return
    out += newExtension
    return out


def removeFileExtension(inFile:str) -> str:
    # remove the original file extension
    return re.sub(r'\.[^\.]+$', '', inFile)


def extractFilenameFromPath(filePath:str) -> str:
    """ extractFilenameFromPath:
            Accepts a string containing the path to a file. Returns only the 
            filename of the file (removes the path to the file).
    """
    # remove all characters except those after the last '/' and return
    return re.sub(r'^.+/([^/]+)$', r'\1', filePath)


def downloadFileFromFTP(ftpHost:str, ftpFilePath:str, filename:str) -> None:
    """ downloadFileFromFTP:
            Accepts three strings: the ftp server, the full path to the file to
            be downloaded, and the name of the file to write locally. Downloads
            the file and writes it. Does not return.
    """
    # initialize the ftp object, then connect and login to the server
    ftp = ftplib.FTP()
    ftp.connect(ftpHost)
    ftp.login()

    # construct the ftp retrieve command
    retrCmd = "RETR " + ftpFilePath
    
    # download the file and write it locally
    with open(filename, 'wb') as fileHandle:
        ftp.retrbinary(retrCmd, fileHandle.write)


def extractContentsFromGZ(filename:str) -> str:
    """ extractContentsFromGZ:
            Accepts a string containing the path to a gzip file. Reads the 
            contents as text and returns the string.
    """
    # Open the gzip file, extract the contents, close the file, and return
    fileHandle  = gzip.open(filename, 'rt')
    fileContent = fileHandle.read()
    fileHandle.close()
    return fileContent


def decompressGZ(gzFileName:str) -> None:
    """ decompressGZ:
            Accepts a compressed gbff file as input. Extracts the gbff file and
            deletes the compressed file. Does not return.
    """
    # get gbff file name
    gbffFileName = os.path.splitext(gzFileName)[0]

    # read gz file
    gzFH = gzip.open(gzFileName, "rt") # read text
    
    # create empty gbff file
    gbffFH = open(gbffFileName, 'w')

    # copy contents of decompressed file
    shutil.copyfileobj(gzFH, gbffFH)

    # close files
    gzFH.close()
    gbffFH.close()

    # remove the gzip file
    os.remove(gzFileName)


def cleanup(paramO:Parameters) -> None:
    """ cleanup:
            Accepts a Parameters object as input. Removes unnecessary
    """
    # constants
    FILE_PATTERN = "*.dtd"

    # get the working directory
    workdir = os.path.dirname(paramO.workdir)

    # get a list of files to remove
    filesToRemove = glob.glob(os.path.join(workdir, FILE_PATTERN))

    # remove each file in the list
    for file in filesToRemove:
        os.remove(file)

