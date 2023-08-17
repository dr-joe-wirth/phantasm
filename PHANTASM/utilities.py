# Author: Joseph S. Wirth

import csv, ftplib, getopt, glob, gzip, logging, math, os, re, shutil, string, sys, time
from http.client import HTTPResponse
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.Entrez import Parser
from PHANTASM.taxonomy.TaxRank import TaxRank
from PHANTASM.Parameter import Parameters


def parseArgs() -> tuple[list, list, Parameters]:
    """ parseArgs:
            Accepts no inputs. Parses arguments from the command line and
            handles invalid syntax/inputs. Returns a three-ple containing:
                1. a list of input genomes
                2. a list of locus tags
                3. a Parameters object
    """
    # import location of dependencies and other useful information
    from phantasm import JOB_1, JOB_2, JOB_3, JOB_4, JOB_5
    from param import BLASTPLUS_DIR, MUSCLE_EXE, FASTTREE_EXE, IQTREE_EXE
    
    # command line flags
    SHORT_OPTS = "i:e:t:m:O:N:L:B:FS:E:h"
    LONG_OPTS = ["input=",
                 "email=",
                 "locus_tags=",
                 "map_file=",
                 "out=",
                 "num_threads=",
                 'max_leaves=',
                 'bootstrap=',
                 'fewer_coregenes',
                 'skip_16S=',
                 'exclude_taxids=',
                 'help']
    INPUT_FLAGS = ("-i","--input")
    EMAIL_FLAGS = ("-e","--email")
    LOCUS_FLAGS = ("-t", "--locus_tags")
    MAP_FLAGS = ("-m", "--map_file")
    OUT_DIR_FLAGS = ("-O", "--out")
    CORES_FLAGS = ("-N", "--num_threads")
    LEAF_FLAGS = ("-L", "--max_leaves")
    BOOTS_FLAGS = ("-B", "--bootstrap")
    REDUCE_FLAGS = ("-F", "--fewer_coregenes")
    SKIP_16S_FLAGS = ("-S", "--skip_16S")
    EXCLUDE_FLAGS = ("-E", "--exclude_taxa")
    
    # default values for some parameters
    DEFAULT_DIR_1 = os.path.join(os.getcwd(), "initialAnalysis")
    DEFAULT_DIR_2 = os.path.join(os.getcwd(), "finalAnalysis")
    DEFAULT_THREADS = 1
    DEFAULT_BOOTS = 0
    DEFAULT_LEAVES = 50
    DEFAULT_REDUCE = False
    
    # misc
    MIN_NUM_GENOMES = 4
    TAG_SEP = ","
    
    # messages
    INVALID_MSG = "ignoring invalid option: "
    UNUSED_MSG = "ignoring unused option: "
    ERR_MSG_1 = "invalid (or missing) email address"
    ERR_MSG_2 = "locus tag(s) (-t or --locus_tag) required for "
    ERR_MSG_3 = "the number of genes is not a multiple of the number of input genomes."
    ERR_MSG_4 = "the specified genome directory is not a directory."
    ERR_MSG_5 = "cannot analyze fewer than " + str(MIN_NUM_GENOMES) + " genomes"
    ERR_MSG_6 = "output file exists"
    ERR_MSG_7 = "could not find specified directory: "
    ERR_MSG_8 = "input is not a directory"
    
    # initialize the logger
    logger = logging.getLogger(__name__ + "." + parseArgs.__name__)
    
    # extract the job name
    job = sys.argv[1]
    
    # parse the additional arguments (start 2 to skip the job name)
    opts,args = getopt.getopt(sys.argv[2:], SHORT_OPTS, LONG_OPTS)
    
    # set default values for a few variables
    threads = DEFAULT_THREADS
    boots = DEFAULT_BOOTS
    leaves = DEFAULT_LEAVES
    reduce = DEFAULT_REDUCE
    
    # set initial empty values for a few other variables
    genomesL = []
    genomeDir = None
    email = ""
    tagsL = []
    mapFN = ""
    outDir = ""
    taxidsFN = None
    excludeFN = None
    
    # extract the arguments from the command
    for opt,arg in opts:
        # input genomes
        if opt in INPUT_FLAGS:
            if os.path.isdir(arg):
                genomeDir = arg
                genomesL = glob.glob(os.path.join(genomeDir, "*"))
            else:
                genomesL = [arg]
        
        # email address
        elif opt in EMAIL_FLAGS:
            email = arg
        
        # locus tags
        elif opt in LOCUS_FLAGS:
            tagsL = arg.split(TAG_SEP)
        
        # human map file
        elif opt in MAP_FLAGS:
            mapFN = os.path.abspath(arg)
        
        # out directory
        elif opt in OUT_DIR_FLAGS:
            outDir = os.path.abspath(arg)
        
        # processors
        elif opt in CORES_FLAGS:
            threads = int(arg)
        
        # max number of leaves
        elif opt in LEAF_FLAGS:
            leaves = int(arg)
        
        # num bootstraps
        elif opt in BOOTS_FLAGS:
            boots = int(arg)
        
        # reduce core genes
        elif opt in REDUCE_FLAGS:
            reduce = True
        
        # skip 16S
        elif opt in SKIP_16S_FLAGS:
            taxidsFN = os.path.abspath(arg)
        
        # excluded taxids
        elif opt in EXCLUDE_FLAGS:
            excludeFN = os.path.abspath(arg)
        
        # ignore any other flags
        else:
            print(INVALID_MSG + opt)
    
    # process getPhyloMarker
    if job == JOB_1:
        # email address is required
        if not validEmailAddress(email):
            logger.critical(ERR_MSG_1)
            raise ValueError(ERR_MSG_1)
        
        # set the output directory
        outDir = DEFAULT_DIR_1
        
        # report any unused arguments
        for opt,arg in opts:
            if opt in (OUT_DIR_FLAGS + LOCUS_FLAGS + MAP_FLAGS + BOOTS_FLAGS):
                print(UNUSED_MSG + opt)

    # process refinePhylogeny
    elif job in JOB_2:
        # email address is required
        if not validEmailAddress(email):
            logger.critical(ERR_MSG_1)
            raise ValueError(ERR_MSG_1)
        
        # set the output directory
        outDir = DEFAULT_DIR_2
        
        # ensure that locus tags have been provided
        if tagsL == []:
            logger.critical(ERR_MSG_2)
            raise RuntimeError(ERR_MSG_2 + job)
        
        # report any unused arguments
        for opt,arg in opts:
            if opt in (OUT_DIR_FLAGS + MAP_FLAGS + SKIP_16S_FLAGS):
                print(UNUSED_MSG  + opt)
    
    # process knownPhyloMarker
    elif job in JOB_3:
        # email address is required
        if not validEmailAddress(email):
            logger.critical(ERR_MSG_1)
            raise ValueError(ERR_MSG_1)
        
        # set the output directory if one was not specified
        if outDir == "":
            outDir = DEFAULT_DIR_2
        
        # ensure that locus tags have been provided
        if tagsL == []:
            logger.critical(ERR_MSG_2)
            raise RuntimeError(ERR_MSG_2 + job)
        
        # raise error if num locus tags is not a multiple of num genomes
        if len(tagsL) % len(genomesL) != 0:
            logger.critical(ERR_MSG_3)
            raise ValueError(ERR_MSG_3)
        
        # report any unused arguments
        for opt,arg in opts:
            if opt in (MAP_FLAGS + SKIP_16S_FLAGS):
                print(UNUSED_MSG + opt)
    
    # process analyzeGenomes
    elif job == JOB_4:
        # email address is required
        if not validEmailAddress(email):
            logger.critical(ERR_MSG_1)
            raise ValueError(ERR_MSG_1)
        
        # set the output directory if one was not specified
        if outDir == "":
            outDir = DEFAULT_DIR_2
        
        # make sure genomeDir is a directory
        if not os.path.isdir(genomeDir):
            logger.critical(ERR_MSG_4)
            raise NotADirectoryError(ERR_MSG_4)
        
        # make sure that enough genomes are present
        if len(genomesL) < MIN_NUM_GENOMES:
            logger.critical(ERR_MSG_4)
            raise ValueError(ERR_MSG_5)
        
        # report any unused arguments
        for opt,arg in opts:
            if opt in (LOCUS_FLAGS + LEAF_FLAGS + SKIP_16S_FLAGS + EXCLUDE_FLAGS):
                print(UNUSED_MSG + opt)

    # process rankPhyloMarkers
    elif job == JOB_5:
        # if set, the output file will be set to `outDir`; move to `outFN`
        if outDir != "":
            outFN = os.path.abspath(outDir)
            
            # make sure the file doesn't exist
            if os.path.exists(outFN):
                logger.critical(ERR_MSG_6)
                raise FileExistsError(ERR_MSG_6)
            
        # set outFN to False to ensure that the file replacement below works
        else:
            outFN = False
        
        # the working directory will be set to `genomeDir`; move to `outDir`
        outDir = genomeDir
        
        # make sure output directory exists 
        if not os.path.exists(outDir):
            logger.critical(ERR_MSG_7 + outDir)
            raise FileNotFoundError(ERR_MSG_7 + outDir)
        
        # make sure output directory is a directory
        if not os.path.isdir(outDir):
            logger.critical(ERR_MSG_8 + outDir)
            raise NotADirectoryError(ERR_MSG_8 + outDir)
        
        # report any unused arguments
        for opt,arg in opts:
            if opt in (EMAIL_FLAGS + LOCUS_FLAGS + MAP_FLAGS + LEAF_FLAGS + BOOTS_FLAGS + REDUCE_FLAGS + EXCLUDE_FLAGS + SKIP_16S_FLAGS):
                print(UNUSED_MSG + opt)
    
    # create a Parameters object using the extracted values:
    paramO = Parameters(email,
                        outDir,
                        BLASTPLUS_DIR,
                        MUSCLE_EXE,
                        FASTTREE_EXE,
                        IQTREE_EXE,
                        threads,
                        leaves,
                        boots,
                        reduce)
    
    # add the skip16S taxids file if requested
    if job == JOB_1 and taxidsFN is not None:
        if not os.path.isfile(taxidsFN):
            raise FileNotFoundError(taxidsFN)
        paramO.taxidsFN = taxidsFN

    # add the excluded taxids file if requested
    if job in (JOB_1, JOB_2, JOB_3) and excludeFN is not None:
        if not os.path.isfile(excludeFN):
            raise FileNotFoundError(excludeFN)
        paramO.excludedTaxidsFN = excludeFN
    
    # replace the genomes directory and map file if analyzing genomes
    if job == JOB_4:
        paramO.genbankFilePath = os.path.join(genomeDir, "*")
        paramO.fileNameMapFN = os.path.abspath(mapFN)
        
        # ensure that the provided map file is valid
        checkForValidHumanMapFile(paramO)
    
    # replace the output file in the paramO object
    elif job == JOB_5 and outFN:
        paramO.phyloMarkersFN = outFN
    
    # make sure that the specified executables are accessible
    checkForValidExecutables(paramO)

    return genomesL, tagsL, paramO


def getHelpMessage(task:str) -> str:
    """ getHelpMessage
            Accepts a task (job name) as an input. Returns the appropriate help
            message for the given task.
    """
    # constants
    from phantasm import PHANTASM_PY, JOB_0A, JOB_0B, JOB_1, JOB_2, JOB_3, JOB_4, JOB_5
    GAP = 4*" "

    # messages
    HELP_0 = "Usage: " + PHANTASM_PY + " TASK [OPTIONS]\n\n" + \
             "Available tasks:\n" + \
             GAP + f"{JOB_1:20}" + "identify phylogenetic markers (option 1, step 1)\n" + \
             GAP + f"{JOB_2:20}" + "refine phylogeny and perform phylogenomic analyses (option 1, step 2)\n" + \
             GAP + f"{JOB_3:20}" + "idenfify genomes from known phylogenetic marker and perform phylogenomic analyses (option 2)\n" + \
             GAP + f"{JOB_4:20}" + "perform phylogenomic analyses on a user-specified set of genomes (option 3)\n" + \
             GAP + f"{JOB_5:20}" + "identify core genes and rank phylogenetic markers\n" + \
             GAP + f"{JOB_0A:20}" + "print this message\n" + \
             GAP + f"{JOB_0B:20}" + "print the version\n\n" + \
             "Run '" + PHANTASM_PY + " TASK --help' for more information on a task.\n\n"
    HELP_1 = "Usage: " + PHANTASM_PY + " " + JOB_1 + " [-ieNLF]\n\n" + \
             "Required arguments:\n" + \
             GAP + "-i, --input <file>           gbff file or a directory containing gbff files\n" + \
             GAP + "-e, --email <email>          email address\n\n" + \
             "Optional arguments:\n" + \
             GAP + "-N, --num_threads <int>      number of processors to use                    [default: 1]\n" + \
             GAP + "-L, --max_leaves <int>       maximum number of leaves in the species tree   [default: 50]\n" + \
             GAP + "-F, --fewer_coregenes        limit the core genes to those with ≤5% gaps    [default: no limiting]\n" + \
             GAP + "-S, --ski6_16S <file>        specify related taxids instead of relying on 16S homology\n" + \
             GAP + "-E, --exclude_taxa <file>    specify taxids to exclude\n" + \
             GAP + "-h, --help                   print this message\n\n" + \
             "Results:\n" + \
             GAP + "'initialAnalysis/putativePhylogeneticMarkers.txt'\n\n"
    HELP_2 = "Usage: " + PHANTASM_PY + " " + JOB_2 + " [-ietNLBF]\n\n" + \
             "Required arguments:\n" + \
             GAP + "-i, --input <file>           gbff file or a directory containing gbff files\n" + \
             GAP + "-e, --email <email>          email address\n" + \
             GAP + "-t, --locus_tags <str>       a comma-separated list of locus tags to use a phylogenetic markers\n\n" + \
             "Optional arguments:\n" + \
             GAP + "-N, --num_threads <int>      number of processors to use                    [default: 1]\n" + \
             GAP + "-L, --max_leaves <int>       maximum number of leaves in the species tree   [default: 50]\n" + \
             GAP + "-B, --bootstrap <int>        number of bootstraps to perform                [default: no bootstrapping]\n" + \
             GAP + "-F, --fewer_coregenes        limit the core genes to those with ≤5% gaps    [default: no limiting]\n" + \
             GAP + "-E, --exclude_taxa <file>    specify taxids to exclude\n" + \
             GAP + "-h, --help                   print this message\n\n" + \
             "Results:\n" + \
             GAP + "'finalAnalysis/aai_matrix.txt'\n" + \
             GAP + "'finalAnalysis/aai_heatmap.pdf'\n" + \
             GAP + "'finalAnalysis/ani_matrix.txt'\n" + \
             GAP + "'finalAnalysis/ani_heatmap.pdf'\n" + \
             GAP + "'finalAnalysis/coreGenesSummary.txt\n" + \
             GAP + "'finalAnalysis/speciesTree.nwk'\n" + \
             GAP + "'finalAnalysis/speciesTree_outgroupPruned.nwk'\n\n"
    HELP_3 = "Usage: " + PHANTASM_PY + " " + JOB_3 + " [-ietONLBF]\n\n" + \
             "Required arguments:\n" + \
             GAP + "-i, --input <file>           gbff file or a directory containing gbff files\n" + \
             GAP + "-e, --email <email>          email address\n" + \
             GAP + "-t, --locus_tags <str>       a comma-separated list of locus tags to use a phylogenetic markers\n\n" + \
             "Optional arguments:\n" + \
             GAP + "-O, --out <dir>              output directory                               [default: '" + os.path.join(".", "finalAnalysis") + "']\n" + \
             GAP + "-N, --num_threads <int>      number of processors to use                    [default: 1]\n" + \
             GAP + "-L, --max_leaves <int>       maximum number of leaves in the species tree   [default: 50]\n" + \
             GAP + "-B, --bootstrap <int>        number of bootstraps to perform                [default: no bootstrapping]\n" + \
             GAP + "-F, --fewer_coregenes        limit the core genes to those with ≤5% gaps    [default: no limiting]\n" + \
             GAP + "-E, --exclude_taxa <file>    specify taxids to exclude\n" + \
             GAP + "-h, --help                   print this message\n\n" + \
             "Results (directory may vary if '-O' or '--out' used):\n" + \
             GAP + "'finalAnalysis/aai_matrix.txt'\n" + \
             GAP + "'finalAnalysis/aai_heatmap.pdf'\n" + \
             GAP + "'finalAnalysis/ani_matrix.txt'\n" + \
             GAP + "'finalAnalysis/ani_heatmap.pdf'\n" + \
             GAP + "'finalAnalysis/coreGenesSummary.txt\n" + \
             GAP + "'finalAnalysis/speciesTree.nwk'\n" + \
             GAP + "'finalAnalysis/speciesTree_outgroupPruned.nwk'\n\n"
    HELP_4 = "Usage: " + PHANTASM_PY + " " + JOB_4 + " [-iemONLBF]\n\n" + \
             "Required arguments:\n" + \
             GAP + "-i, --input <file>         gbff file or a directory containing gbff files\n" + \
             GAP + "-e, --email <email>        email address\n" + \
             GAP + "-m, --map_file <file>      a file with two tab-separated columns (no headers): filename, taxon name\n\n" + \
             "Optional arguments:\n" + \
             GAP + "-O, --out <dir>            output directory                               [default: '" + os.path.join(".", "finalAnalysis") + "']\n" + \
             GAP + "-N, --num_threads <int>    number of processors to use                    [default: 1]\n" + \
             GAP + "-B, --bootstrap <int>      number of bootstraps to perform                [default: no bootstrapping]\n" + \
             GAP + "-F, --fewer_coregenes      limit the core genes to those with ≤5% gaps    [default: no limiting]\n" + \
             GAP + "-h, --help                 print this message\n\n" + \
             "Results (directory may vary if '-O' or '--out' used):\n" + \
             GAP + "'finalAnalysis/aai_matrix.txt'\n" + \
             GAP + "'finalAnalysis/aai_heatmap.pdf'\n" + \
             GAP + "'finalAnalysis/ani_matrix.txt'\n" + \
             GAP + "'finalAnalysis/ani_heatmap.pdf'\n" + \
             GAP + "'finalAnalysis/coreGenesSummary.txt\n" + \
             GAP + "'finalAnalysis/speciesTree.nwk'\n" + \
             GAP + "'finalAnalysis/speciesTree_outgroupPruned.nwk'\n\n"
    HELP_5 = "Usage: " + PHANTASM_PY + " " + JOB_5 + " [-iON]\n\n" + \
             "Required arguments:\n" + \
             GAP + "-i, --input <dir>          directory containing phantasm results\n\n" + \
             "Optional arguments:\n" + \
             GAP + "-O, --out <file>           output file                    [default: '" + os.path.join("<dir>", "putativePhylogeneticMarkers.txt") + "']\n" + \
             GAP + "-N, --num_threads <int>    number of processors to use    [default: 1]\n" + \
             GAP + "-h, --help                 print this message\n\n" + \
             "Results:\n" + \
             GAP + "'<dir>/putativePhylogeneticMarkers.txt' (or at the specified file location)\n\n"
    ADDITIONAL_HELP = "see https://github.com/dr-joe-wirth/phantasm for additional documentation.\n"
    
    # return the appropriate help message
    if task == JOB_0A:
        return HELP_0 + ADDITIONAL_HELP
    elif task == JOB_1:
        return HELP_1 + ADDITIONAL_HELP
    elif task == JOB_2:
        return HELP_2 + ADDITIONAL_HELP
    elif task == JOB_3:
        return HELP_3 + ADDITIONAL_HELP
    elif task == JOB_4:
        return HELP_4 + ADDITIONAL_HELP
    elif task == JOB_5:
        return HELP_5 + ADDITIONAL_HELP


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
    
    # initialize logger
    logger = logging.getLogger(__name__ + "." + getTaxidsFromFile.__name__)
    
    # print status
    print(PRNT_1, end="", flush=True)
    logger.info(PRNT_1)

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
        logger.error(ERR_MSG_1)
        raise BaseException(ERR_MSG_1)
    
    for txid in outL:
        # all ids should be coerceable to int
        try:
            int(txid)
        except:
            logger.error(ERR_MSG_2)
            raise BaseException(ERR_MSG_2)

    # print status
    print(DONE)
    logger.info(DONE)

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
    # constant
    ERR_MSG = "Entrez.email has not been set."
    
    # initialize logger
    logger = logging.getLogger(__name__ + "." + checkEntrezEmail.__name__)
    
    # make sure Entrez.email can be set
    if Entrez.email is None:
        if email is not None:
            Entrez.email = email
        else:
            logger.error(ERR_MSG)
            raise BaseException(ERR_MSG)


def loadHumanMap(humanMapFN:str) -> dict:
    """ loadHumanMap:
            Accepts a string indicating the filename of the human map file as
            input. Constructs a dictionary keyed by gbff filenames with the co-
            rresponding human names as the values. Returns the dictionary.
    """
    # constants
    ERR_MSG = "human map file is improperly formatted"
    FILE_NAME_IDX  = 0
    HUMAN_NAME_IDX = 1
    
    # initialize logger
    logger = logging.getLogger(__name__ + "." + loadHumanMap.__name__)

    # read the file into memory
    parsed = parseCsv(humanMapFN, '\t')

    # remove any trailing empty lines from the file
    while parsed[-1] == []:
        parsed = parsed[:-1]

    # create the dict
    humanMapD = dict()

    # go through each row in the file and add the data to the dictionary
    for row in parsed:
        # invalid format if the row does not have exactly two columns
        if len(row) != 2:
            logger.error(ERR_MSG)
            raise ValueError(ERR_MSG)

        # keys are filenames; values are human names
        humanMapD[row[FILE_NAME_IDX]] = row[HUMAN_NAME_IDX]
    
    return humanMapD


def ncbiEfetchById(id, database:str, mode:str='text', rettype:str='xml', \
                       retmode:str=None, email:str=None) -> Parser.ListElement:
    """ ncbiEfetchById:
            Accepts an NCBI id (or a list of ids) and an NCBI database as 
            input. Uses NCBI's Efetch tool to retrieve full records from NCBI
            for the provided id(s). Returns the retrieved records.
    """
    # make sure Entrez.email has been set
    checkEntrezEmail(email)

    result = __efetch(database, id, mode, rettype, retmode)

    return result


def __efetch(db:str, id:str, mode:str, rettype:str, retmode:str) -> \
                                                            Parser.ListElement:
    """ efetch:
            Accepts a database, a uid, a query mode, a return type, and return
            mode as inputs. Obtains efetch records from NCBI. Handles errors
            associated with failed requests to NCBI. Returns the records as a
            "list" (Entrez format).
    """
    # constants
    MAX_ATTEMPTS = 5
    PAUSE = 0.5
    ERR_MSG = "failed to retrieve records from NCBI after " + str(MAX_ATTEMPTS) + " tries"

    logger = logging.getLogger(__name__ + "." + __efetch.__name__)

    # keep trying until it works or fails too many times
    result = None
    numAttempts = 0
    while result is None and numAttempts < MAX_ATTEMPTS:
        try:
            # retrieve the result from NCBI
            handle = Entrez.efetch(db=db, id=id, mode=mode, rettype=rettype, retmode=retmode)
            result = Entrez.read(handle)
            handle.close()
        except:
            # reset if this fails; pause so NCBI is not overwhelmed
            result = None
            numAttempts += 1
            time.sleep(PAUSE)
    
    # raise an error if the result was not retrieved
    if result is None:
        logger.critical("database: " + db)
        logger.critical("query:    " + id)
        logger.critical(ERR_MSG)
        raise RuntimeError(ERR_MSG)
    
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

    result = __esearch(database, searchStr, retmax)

    # grab all results if requested (default)
    if retrieveAll:
        # determine the number of records found
        numFound = int(result['Count'])

        # repeat search if the amount retrieved is less than the amount found
        if numFound > retmax:
            retmax = numFound

            result = __esearch(database, searchStr, retmax)
    
    # return the list of ids that were found
    return result['IdList']


def __esearch(db:str, term:str, retmax:int) -> Parser.DictionaryElement:
    """ esearch:
            Accepts a database, a search term, and a max number of records to
            return as inputs. Obtains esearch records from NCBI. Handles errors
            associated with failed requests to NCBI. Returns the records as a
            "dict" (Entrez format) containing the search result.
    """
    # constants
    MAX_ATTEMPTS = 5
    PAUSE = 0.5
    ERR_MSG = "failed to connect to NCBI after " + str(MAX_ATTEMPTS) + " tries"

    logger = logging.getLogger(__name__ + "." + __esearch.__name__)

    # keep trying until it works or fails too many times
    result = None
    numAttempts = 0
    while result is None and numAttempts < MAX_ATTEMPTS:
        # retrieve the result from NCBI
        try:
            handle = Entrez.esearch(db=db, term=term, retmax=retmax)
            result = Entrez.read(handle)
            handle.close()
        
        # reset if it fails; pause to protect NCBI
        except:
            result = None
            numAttempts += 1
            time.sleep(PAUSE)
    
    # raise an error if a result was not retrieved
    if result is None:
        logger.critical(ERR_MSG)
        logger.critical("database: " + db)
        logger.critical("query:    " + term)
        raise RuntimeError(ERR_MSG)

    return result


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
    ERR_MSG_1 = "invalid input: expected a list of ids or a single id"
    ERR_MSG_2 = "invalid id present in the input list: "
    
    logger = logging.getLogger(__name__ + "." + ncbiSummaryFromIdList.__name__)

    # check that Entrez.email has been set
    checkEntrezEmail(email)

    # if not a list, assume it is an int or a string
    if type(idList) is str:
        try:
            int(idList)
        except ValueError:
            logger.error(ERR_MSG_1)
            raise ValueError(ERR_MSG_1)

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
                logger.error(ERR_MSG_2 + str(id))
                raise ValueError(ERR_MSG_2 + str(id))

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
    result = __esummary(queryL[0], database)

    # if there are multiple queries, then keep querying NCBI until done
    for i in range(1,len(queryL)):
        query = queryL[i]

        nextResult = __esummary(query, database)

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
            result[RESULT_K1][RESULT_K2].append(nextResult[RESULT_K1][RESULT_K2])

    return result


def __esummary(query:str, db:str):
    """ esummary:
            Accepts a query string and a database as inputs. Obtains summary
            records from NCBI. Handles errors associated with failed requests
            to NCBI. Returns the records as either a Parser.ListElement or,
            depending on the database, a Parser.DictionaryElement.
    """
    # constants
    VALIDATE  = False   # Bio.Entrez may throw a tantrum if modified
    MAX_ATTEMPTS = 5
    PAUSE = 0.5
    ERR_MSG = "failed to retrieve the result from NCBI after " + str(MAX_ATTEMPTS) + " tries"

    logger = logging.getLogger(__name__ + "." + __esummary.__name__)

    # keep trying until it works or fails too many times
    result = None
    numAttempts = 0
    while result is None and numAttempts < MAX_ATTEMPTS:
        # retrieve the result from NCBI
        try:
            handle = Entrez.esummary(id=query, db=db)
            result = Entrez.read(handle, validate=VALIDATE)
            handle.close()

        # reset on a failure; pause to protect NCBI
        except:
            result = None
            numAttempts += 1
            time.sleep(PAUSE)
    
    # raise an error if a result was not obtained
    if result is None:
        logger.error(ERR_MSG)
        logger.error("database: " + db)
        logger.error("query:    " + query)
        raise RuntimeError(ERR_MSG)

    return result


def ncbiELinkFromIdList(idL:list, dbFrom:str, dbTo:str, keepTrying=True) -> list:
    """ ncbiELinkFromIdList:
            Accepts a list of uids, a string indicating which database the uids
            correspond to, a string indicating which database to link to, and a
            boolean indicating if it should repeat a query after a failed query
            to NCBI as inputs. Retrieves and returns a list of the link records
            for the query.
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

        # save the result
        result.extend(__elink(handle, keepTrying))
        
        # close the handle
        handle.close()

    return result


def __elink(handle:HTTPResponse, keepTrying:bool) -> Parser.ListElement:
    """ elink:
            Accepts a handle to an HTTP request (ncbi query) and a boolean ind-
            icating if it should repeat a query after a failed queryto NCBI as
            inputs. Obtains elink records from NCBI. Handles errors associated
            with failed requests to NCBI, but does not log the error message.
            Returns the records as a "list" (Entrez format).
    """
    # constants
    MAX_ATTEMPTS = 5
    PAUSE = 0.5
    ERR_MSG = "failed to retreive the result from NCBI after " + str(MAX_ATTEMPTS) + " tries"

    # keep trying until it works
    result = None
    numAttempts = 0
    while result is None and numAttempts < MAX_ATTEMPTS:
        try:
            result = Entrez.read(handle)

        # reset on a failure
        except:
            result = None
            numAttempts += 1
            
            # if we are to retry query, then pause to protect ncbi
            if keepTrying:
                time.sleep(PAUSE)
            
            # stop looping if we are not retrying the query
            else: break
    
    # raise an error if a result was not obtained
    if result is None:
        raise RuntimeError(ERR_MSG)

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
    DOMAIN = 'domain'
    DOMAIN_NCBI = 'superkingdom'
    LINEAGE = 'LineageEx'
    RANK = 'Rank'
    TAXID = 'TaxId'
    SCI_NAME = 'ScientificName'
    
    # ensure the target rank is valid (let TaxRank constructor check validity)
    if type(targetRank) is str:
        targetRank = TaxRank(targetRank)
    
    # get the string of the input rank
    targetRank = str(targetRank)

    # make sure to use 'superkingdom' if requesting a domain from NCBI
    if targetRank == DOMAIN:
        targetRank = DOMAIN_NCBI
    
    # import data from the taxonomy database
    taxRecords = ncbiEfetchById(taxids, DATABASE)
    
    # for each record retrieved
    parentTaxid = []
    for tax in taxRecords:
        # get the lineage
        allLineages = tax[LINEAGE]

        # go through the lineage until target rank is found
        for lineage in allLineages:
            if lineage[RANK] == targetRank:
                newEntry = {'txid': lineage[TAXID], 
                            'name': lineage[SCI_NAME]}
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
    with open(csvFN, 'r') as fh:
        parsed = list(csv.reader(fh, delimiter=delim))

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
    out = os.path.splitext(inFile)[0]

    # make sure newExtension will follow a period
    if newExtension[0] != '.' and out[-1] != '.':
        out += '.'

    # add the extension and return
    out += newExtension
    return out


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


def checkForValidInputGenomes(gbffL:list) -> None:
    """ checkForValidInputGenomes:
            Accepts a list of filenames (str) as input. Raises an error if any
            of the expected features of a suitable genbank file are missing.
            Does not return.
    """
    # constants
    FORMAT = 'genbank'
    CDS = "CDS"
    LOCUS_TAG = 'locus_tag'
    ERR_MSG_1 = "no input genome files were found"
    ERR_MSG_2 = "' is not in genbank format"
    ERR_MSG_3 = "' is missing 'locus_tag' for one or more of its CDS features"
    MSG_START = "'"
    
    # initialize logger
    logger = logging.getLogger(__name__ + "." + checkForValidInputGenomes.__name__)

    # invalid if an empty list is provided
    if len(gbffL) == 0:
        logger.error(ERR_MSG_1)
        raise FileNotFoundError(ERR_MSG_1)
    
    # for each file in the list
    for gbFN in gbffL:
        # open it as a genbank file
        parsed = SeqIO.parse(gbFN, FORMAT)

        # if there are no records, then the file is invalid
        if len(list(parsed)) == 0:
            logger.error(MSG_START + gbFN + ERR_MSG_2)
            raise ValueError(MSG_START + gbFN + ERR_MSG_2)

        # for each record in the file
        rec:SeqRecord
        for rec in parsed:
            feature:SeqFeature
            for feature in rec.features:
                # if the feature is a CDS
                if feature.type == CDS:
                    # try to find a locus tag
                    try:
                        feature.qualifiers[LOCUS_TAG]
                    
                    # raise an error if there are no locus tags
                    except:
                        logger.error(MSG_START + gbFN + ERR_MSG_3)
                        raise ValueError(MSG_START + gbFN + ERR_MSG_3)


def checkForValidExecutables(paramO:Parameters) -> None:
    """ checkForValidExecutables:
            Accepts a Parameters object as input. Raises an error if any of the
            necessary executable files are missing or not executable. Does not
            return.
    """
    # constants
    ERR_MSG_1 = "invalid blast+ directory specified: '"
    ERR_MSG_2 = "MUSCLE is not executable: '"
    ERR_MSG_3 = "FastTree is not executable: '"
    ERR_MSG_4 = "IQTree is not executable: '"
    MSG_END = "'"
    
    # initialize logger
    logger = logging.getLogger(__name__ + "." + checkForValidExecutables.__name__)

    # get a list of the blast+ executables
    blastExeL = glob.glob(os.path.join(paramO.blastExecutDirPath, "*"))

    # check each blast+ executable
    for exe in blastExeL:
        if not os.access(exe, os.X_OK):
            logger.error(ERR_MSG_1 + paramO.blastExecutDirPath + MSG_END)
            raise ValueError(ERR_MSG_1 + paramO.blastExecutDirPath + MSG_END)

    # check muscle
    if not os.access(paramO.musclePath, os.X_OK):
        logger.error(ERR_MSG_2 + paramO.musclePath + MSG_END)
        raise ValueError(ERR_MSG_2 + paramO.musclePath + MSG_END)

    # check fasttree
    if not os.access(paramO.fastTreePath, os.X_OK):
        logger.error(ERR_MSG_2 + paramO.musclePath + MSG_END)
        raise ValueError(ERR_MSG_3 + paramO.fastTreePath + MSG_END)

    # check iqtree
    if not os.access(paramO.iqTreePath, os.X_OK):
        logger.error(ERR_MSG_2 + paramO.musclePath + MSG_END)
        raise ValueError(ERR_MSG_4 + paramO.iqTreePath + MSG_END)


def checkForValidHumanMapFile(paramO:Parameters) -> None:
    """ checkForValidHumanMapFile:
            Accepts a Parameters object as input. Raises an error if the human
            map file is inproperly formatted. Does not return.
    """
    # constants
    ERR_MSG_1 = "filenames in the map file are not unique"
    ERR_MSG_2 = "human names in the map file are not unique"
    ERR_MSG_3 = "input genbank filenames do not match those in the human map file"
    ERR_MSG_4 = "human names may only contain alphanumeric and the following characters _-|."
    ALLOWED_CHARS = string.ascii_letters + string.digits + "_-|."

    # initialize logger
    logger = logging.getLogger(__name__ + "." + checkForValidHumanMapFile.__name__)

    # extract the relevant data from paramO
    gbffL = glob.glob(paramO.genbankFilePath)
    mapFN = paramO.fileNameMapFN

    # load the human map file
    mapD = loadHumanMap(mapFN)

    # make sure that the filenames in the human map file are unique
    if len(mapD.keys()) != len(set(mapD.keys())):
        logger.error(ERR_MSG_1)
        raise ValueError(ERR_MSG_1)
    
    # make sure that the human names in the human map file are unique
    if len(mapD.values()) != len(set(mapD.values())):
        logger.error(ERR_MSG_2)
        raise ValueError(ERR_MSG_2)

    # get a set of the gbff basenames that were found
    foundGbk = set([os.path.basename(fn) for fn in gbffL])

    # get a set of genome basenames specified in the human map file
    specifiedGbk = set(mapD.keys())

    # make sure the specified files match those that were found
    if not foundGbk == specifiedGbk:
        logger.error(ERR_MSG_3)
        raise ValueError(ERR_MSG_3)

    # make sure there are no illegal characters in the human names
    for humanName in mapD.values():
        for char in humanName:
            if char not in ALLOWED_CHARS:
                logger.error(ERR_MSG_4)
                raise ValueError(ERR_MSG_4)


def getLpsnAge() -> str:
    """ getLpsnAge:
            Accepts no input. Returns a string (yyyy-mm-dd) indicating the date
            that the LPSN data was last updated.
    """
    # constant
    DATE_FN = "date.txt"

    # get the path to the date file
    from param import PHANTASM_DIR
    dateFN = os.path.join(PHANTASM_DIR, 'PHANTASM', 'lpsn_data', DATE_FN)
    
    # get the date from the file
    fh = open(dateFN, 'r')
    date = fh.read()
    fh.close()

    return date
    

def redactEmailAddress() -> str:
    """ redactEmailAddress:
            Accepts no inputs. Redacts the email address from the input command
            to ensure that the user's email remains anonymous in the resulting
            log files. Returns the command as a string with the email address
            redacted.
    """
    # constants
    SEP = r'|~|'
    GREP_FIND = r'([^\@]+)\@(.+)(\..+)$'
    GREP_REPL = r'\1' + SEP + r'\2' + SEP + r'\3'
    X = 'x'
    
    # extract the email address
    for idx in range(len(sys.argv)):
        if sys.argv[idx] in ("-e", "--email"):
            idx += 1
            email = sys.argv[idx]
            break
    
    # extract the address, domain, and extension from the email address
    address,domain,ext = re.sub(GREP_FIND, GREP_REPL, email).split(SEP)
    
    # replace all characters in the address (except the first) with `x`
    address = address[0] + (len(address)-1)*X
    
    # replace all characters in domain (except the first) with `x`
    domain = domain[0] + (len(domain)-1)*X
    
    # construct the obscured email address
    obscured =  address + '@' + domain + ext
    
    # convert the list of arguments to a string command
    cmdL = sys.argv[:idx] + [obscured] + sys.argv[idx+1:]
    
    return " ".join(cmdL)
    

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

