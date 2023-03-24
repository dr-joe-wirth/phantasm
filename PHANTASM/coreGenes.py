# Author: Joseph S. Wirth
# Last edit: September 28, 2022

from __future__ import annotations
import csv, glob, logging, os, re, scipy.stats, string, subprocess, sys
from PHANTASM.utilities import parseCsv, loadHumanMap
from PHANTASM.Parameter import Parameters
from Bio import Phylo, SeqIO
from Bio.Phylo import Newick
from Bio.SeqRecord import SeqRecord
from PHANTASM.downloadGbff import downloadGbffsForRootTaxonomy, _makeHumanMapString, _makeTaxonName
from PHANTASM.taxonomy.Taxonomy import Taxonomy
from param import XENOGI_DIR
sys.path.insert(0,os.path.join(sys.path[0], XENOGI_DIR))
import xenoGI.xenoGI, xenoGI.scores, xenoGI.Tree, xenoGI.genomes, xenoGI.trees

# global constant for separating the data of multiple input genomes
SEP_CHAR = ","

def xenogiInterfacer_1(taxO:Taxonomy, allQryGbksL:list, paramsO:Parameters) \
                                                                   -> Taxonomy:
    """ xenogiInterfacer_1:
            Accepts a Taxonomy object, a list of paths (str) to the query genb-
            ank files, and a Parameters object as inputs. Downloads the ingroup
            and outgroup sequences and puts them into a file structure that is
            expected by xenoGI. Creates the human map file required by xenoGI.
            Returns the outgroup species as a Taxonomy object.
    """
    # extract relevant data from paramD
    gbffFN = paramsO.genbankFilePath
    humanMapFN = paramsO.fileNameMapFN
    taxObjFilePath = paramsO.taxonomyObjectFilePath
    maxNumSeqs = paramsO.maxNumTreeLeaves - len(allQryGbksL) - 1 # outgroup
    
    # determine the directory for the gbff file
    gbffDir = os.path.dirname(gbffFN)

    # identify and download gbff files from NCBI
    outgroup:Taxonomy = downloadGbffsForRootTaxonomy(taxO, maxNumSeqs, paramsO)

    # determine paths to taxonomy files and the extensions
    oldTaxFN = glob.glob(taxObjFilePath).pop()
    dir = os.path.dirname(taxObjFilePath)
    ext = os.path.splitext(taxObjFilePath)[1]

    # save and replace existing taxonomy file
    taxO = taxO.getRoot()
    newTaxFN = os.path.join(dir, taxO.sciName + ext)
    os.remove(oldTaxFN)
    taxO.save(newTaxFN)

    # open the human map file
    filehandle = open(humanMapFN, "a")

    # for each query genome in the list
    for queryGbff in allQryGbksL:
        # get the file name
        basename = os.path.basename(queryGbff)

        # make a symlink to the user's input file
        oldFN = os.path.abspath(queryGbff)
        newFN = os.path.join(gbffDir, basename)
        os.symlink(oldFN, newFN)

        # get the human name for the query file name
        humanName = _humanNameFromQueryGenbankFN(queryGbff)

        # make the human map string and append it to the human map file
        humanMapStr = _makeHumanMapString(humanName, basename)
        filehandle.write(humanMapStr)
    
    # close the human map file
    filehandle.close()

    return outgroup


def _humanNameFromQueryGenbankFN(queryGbff:str) -> str:
    """ humanNameeFromQueryGenbankFN:
            Accepts a filename (str) of a query genome as input. Returns the
            corresponding human name as a string for use in the human map file.
    """
    # constants
    ALLOWED_CHARS = string.ascii_letters + string.digits + "_"
    
    # get the filename without the extension or path
    noext = os.path.splitext(os.path.basename(queryGbff))[0]

    # initialize the human name string
    humanName = ""

    # for each character in the filename
    for char in noext:
        # keep allowed characters
        if char in ALLOWED_CHARS:
            humanName += char
        
        # replace illegal characters with underscores
        else:
            humanName += "_"

    return humanName


def parseGenbank(paramD:dict) -> None:
    """ parseGenbank:
            Accepts the parameter dictionary as input. Calls the 'parseGenbank'
            command of xenoGI. Does not return.
    """
    # constants
    PRINT = 'Parsing genbank files ... '
    DONE = 'Done.'

    # intialize logger
    logger = logging.getLogger(__name__ + "." + parseGenbank.__name__)

    # parse the genbank files
    print(PRINT, end='', flush=True)
    logger.info(PRINT)
    xenoGI.xenoGI.parseGenbankWrapper(paramD)
    print(DONE)
    logger.info(DONE)


def allVsAllBlast(paramD:dict) -> None:
    """ allVsAllBlast:
            Accepts the parameter dictionary as input. Calls the 'runBlast'
            command of xenoGI. Does not return.
    """
    # constants
    PRINT = 'Running all pairwise blastp comparisons ... '
    DONE = 'Done.'
    
    # initialize logger
    logger = logging.getLogger(__name__ + "." + allVsAllBlast.__name__)

    # run all pairwise blast comparisons
    print(PRINT, end='', flush=True)
    logger.info(PRINT)
    xenoGI.xenoGI.runBlastWrapper(paramD)
    print(DONE)
    logger.info(DONE)


def copyExistingBlastFiles(oldParamO:Parameters, newParamO:Parameters) -> None:
    """ copyExistingBlastFiles:
            Accepts two Parameters objects as inputs. Identifies any existing
            blast tables that could be used in the new folder. Creates modifi-
            ed files that are equivalent to the original but have replaced the
            old gene number with the new gene number. Does not return.
    """
    # constants
    PRINT = "Copying existing blastp comparisons ... "
    DONE = 'Done.'
    
    logger = logging.getLogger(__name__ + "." + copyExistingBlastFiles.__name__)

    # print job-start statement
    print(PRINT, end='', flush=True)
    logger.info(PRINT)

    # extract the new blast directory from newParamD
    newBlastDir = os.path.splitext(newParamO.blastFilePath)[0][:-1]
    blastJoinStr = newParamO.blastFileJoinStr

    # create the new blast directory (if it doesn't already exist)
    if glob.glob(newBlastDir) == []:
        os.mkdir(newBlastDir)

    # determine which blastp tables need to be made
    filesToCopyD = __getBlastFilesToCopy(oldParamO, newParamO)

    # get a dictionary to link locus tags to their new names
    locusTagToNewNameD = __getLocusTagToNameD(newParamO)

    # for each file to copy
    for oldFilename in filesToCopyD:
        # get the new
        newFilename = filesToCopyD[oldFilename]

        # modify the query and subject names to match the new gene numbers
        # save the modified string to the new filename
        __updateQuerySubjectNames(oldFilename, newFilename, blastJoinStr, \
                                                            locusTagToNewNameD)
    
    print(DONE)
    logger.info(DONE)


def __getBlastFilesToCopy(oldParamO:Parameters, newParamO:Parameters) -> dict:
    """ getBlastFilesToCopy:
            Accepts the old and new Parameters objects as inputs. Determines
            which blastp comparisons have already been performed. Creates a di-
            ctionary where the key is the path to the existing blastp files and
            the value is the path where the existing file should be copied. Re-
            turns the dictionary.
    """
    # constant
    DELIM = "\t"

    # extract relevant data from oldParamD
    oldBlastPath = oldParamO.blastFilePath # path/*

    # extract relevant data from newParamD
    blastJoin = newParamO.blastFileJoinStr
    wgsMapFN = newParamO.fileNameMapFN
    blastSplt = os.path.splitext(newParamO.blastFilePath)

    # get the blast file extension and the new blast dicrectory
    blastExt =  blastSplt[1]
    newBlastDir = blastSplt[0][:-1]

    # parse human map file into a list of lists
    parsedWgsMap = parseCsv(wgsMapFN, DELIM)

    # create a set of all the filenames needed
    allFilesNeeded = set()
    for row1 in parsedWgsMap:
        for row2 in parsedWgsMap:
            blastFN = row1[1] + blastJoin + row2[1] + blastExt
            allFilesNeeded.add(blastFN)
    
    # get the full path to all existing blastp tables
    oldBlastPaths = glob.glob(oldBlastPath)

    # make a dictionary to enable look-up old filepath from basename
    allOldFilesD = dict()
    for filename in oldBlastPaths:
        allOldFilesD[os.path.basename(filename)] = filename

    # determine which needed files already exist and need to be copied
    allOldFiles = set(allOldFilesD.keys())
    filesToCopy = allFilesNeeded.intersection(allOldFiles)

    # create a dict where key is an old filepath and val is the new destination
    filesToCopyD = dict()
    for filename in filesToCopy:
        oldFilePath = allOldFilesD[filename]
        newFilePath = os.path.join(newBlastDir, filename)
        filesToCopyD[oldFilePath] = newFilePath

    return filesToCopyD


def __getLocusTagToNameD(paramO:Parameters) -> dict:
    """ getLocusTagToNameD:
            Accecpts a Parameters object as input. Uses genesO.geneInfoD to cr-
            eate a new dictionary keyed by locus tag with the gene blast string
            as the value (eg. "999_Methanosarcina_acetivorans-MA_RS05270"). Re-
            turns the newly created dictionary.
    """
    # constants
    LOCUS_TAG_IDX  = 2
    BLAST_NAME_IDX = 0
    GREP_FIND_PREFIX = r'\d+_(.+)-'
    GREP_REPL = r'\1'

    # get the genesO object for the new folder and initialize its geneInfoD
    paramD = paramO.toDict()
    newGenesO = xenoGI.xenoGI.loadGenomeRelatedData(paramD)[1] # genesO object
    newGenesO.initializeGeneInfoD(paramO.geneInfoFN)

    # create a new dictionary: key = human name + locus_tag, val = blast name
    locusTagToNameD = dict()
    for geneNum in newGenesO.geneInfoD.keys():
        # extract relevant data from genesInfoD
        locusTag  = newGenesO.geneInfoD[geneNum][LOCUS_TAG_IDX]
        blastName = newGenesO.geneInfoD[geneNum][BLAST_NAME_IDX]
        humanName = re.sub(GREP_FIND_PREFIX+locusTag, GREP_REPL, blastName)

        # store the data in the dictionary
        locusTagToNameD[(humanName, locusTag)] = blastName    
    
    return locusTagToNameD


def __updateQuerySubjectNames(existingFile:str, newFile:str, blastJoinStr:str,\
                                              locusTagToNewNameD:dict) -> None:
    """ updateQuerySubjectNames:
            Accepts a path to an existing blast file, the path to save the mod-
            ified file, the blast join string, and a dictionary as inputs. Rea-
            ds the blast file and changes the names of the query and subject so
            that their gene numbers are compatible with the current folder. Do-
            es not return.
    """
    # constants
    DELIM = "\t"
    GREP_FIND_1 = r'^\d+_[^\|]+\|[^\|]+\|[^-]+-(\S+)$'
    GREP_FIND_2 = r"^\d+_[^-]+-(\S+)$"
    GREP_REPL = r"\1"
    QRY_COL = 0
    SBJ_COL = 1
    
    # open the old file for reading
    oldHandle = open(existingFile, "r", newline='')
    oldReader = csv.reader(oldHandle, delimiter=DELIM)
    
    # open the new file for writing
    newHandle = open(newFile, "w", newline='')
    newWriter = csv.writer(newHandle, delimiter=DELIM)
    
    # get the human names from the filename
    comparison = os.path.splitext(os.path.basename(newFile))[0]
    qryName,sbjName = comparison.split(blastJoinStr)

    # for each row in the blast table ...
    for row in oldReader:
        # ... extract the locus tag for the subject and the query
        locusTag_qry = re.sub(GREP_FIND_1, GREP_REPL, row[QRY_COL])
        locusTag_sbj = re.sub(GREP_FIND_1, GREP_REPL, row[SBJ_COL])

        # if grep didn't work for query (eg user_input), then try GREP_FIND_2
        if locusTag_qry == row[QRY_COL]:
            locusTag_qry = re.sub(GREP_FIND_2, GREP_REPL, row[QRY_COL])
        
        # if grep didn't work for subject (eg user_input), then try GREP_FIND_2
        if locusTag_sbj == row[SBJ_COL]:
            locusTag_sbj = re.sub(GREP_FIND_2, GREP_REPL, row[SBJ_COL])

        # ... replace the query and subject names with the new names
        row[QRY_COL] = locusTagToNewNameD[(qryName,locusTag_qry)]
        row[SBJ_COL] = locusTagToNewNameD[(sbjName,locusTag_sbj)]

        # ... write the row to the new file
        newWriter.writerow(row)
    
    # close the files
    oldHandle.close()
    newHandle.close()


def calculateCoreGenes(paramO:Parameters) -> None:
    """ calculateCoreGenes:
            Accepts a Parameters object as input. Calculates the hardcore genes
            via the 'createAabrhL' function from xenoGI. Does not return.
    """
    # constants
    PRINT_1 = 'Calculating core genes ... '
    DONE = 'Done.'
    
    logger = logging.getLogger(__name__ + "." + calculateCoreGenes.__name__)

    # parse parameters into shorter variable names
    strainInfoFN = paramO.strainInfoFN
    blastFileJoinStr = paramO.blastFileJoinStr
    blastDir,blastExt = paramO.blastFilePath.split("*")
    evalueThresh = paramO.evalueThresh
    alignCoverThresh = paramO.alignCoverThresh
    percIdentThresh = paramO.percIdentThresh
    aabrhFN = paramO.aabrhFN

    # read the strain info into a list
    strainNamesL = xenoGI.xenoGI.readStrainInfoFN(strainInfoFN)

    # determine the all-against-all best reciprocal hits (hardcore genes)
    print(PRINT_1, end='', flush=True)
    logger.info(PRINT_1)
    xenoGI.scores.createAabrhL(strainNamesL,
                               blastFileJoinStr,
                               blastDir,
                               blastExt,
                               evalueThresh,
                               alignCoverThresh,
                               percIdentThresh,
                               aabrhFN)
    print(DONE)
    logger.info(DONE)


def makeSpeciesTree(allQryGbksL:list, paramO:Parameters, outgroup:Taxonomy, \
                                      program:str, makeGeneTrees:bool) -> None:
    """ makeSpeciesTree:
            Accepts a list of the query genbank files, a Parameters object,
            an outgroup (Taxonomy), and a string indicating which program to
            use as inputs. Uses xenoGI to make gene trees, and then builds ano-
            ther tree based on the concatenated alignment of the hardcore genes
            and roots the tree on the specified outgroup. Does not return.
    """
    # constants
    PRINT_1A = 'Aligning core genes and making gene trees ... '
    PRINT_1B = 'Aligning core genes ... '
    PRINT_2A = 'Using FastTree to make the species tree from a ' + \
                                'concatenated alignment of the core genes ... '
    PRINT_2B = 'Using IQTree to make the species tree with bootstrap ' + \
                 'supports from a concatenated alignment of the core genes ...'
    DONE = 'Done.'
    FAST_TREE = 'fasttree'
    IQTREE = 'iqtree'
    ERR_MSG = 'Invalid program specified.'
    
    logger = logging.getLogger(__name__ + "." + makeSpeciesTree.__name__)

    # align hardcore genes and build gene trees with fasttree
    if makeGeneTrees:
        print(PRINT_1A, end='', flush=True)
        logger.info(PRINT_1A)
        __makeGeneTreesWrapper(paramO)
        print(DONE)
        logger.info(DONE)
    
    else:
        print(PRINT_1B, end='', flush=True)
        logger.info(PRINT_1B)
        __alignGenesWrapper(paramO)
        print(DONE)
        logger.info(DONE)

    # parse data from paramO
    speTreWorkDir = paramO.makeSpeciesTreeWorkingDir
    catAlnFN = paramO.concatenatedAlignmentFN
    keyFN = paramO.famToGeneKeyFN
    mapFN = paramO.fileNameMapFN
    speTreeFN = paramO.speciesTreeFN

    # load the human map file
    humanMapD = loadHumanMap(mapFN)

    # make a list of the human names for the query genomes
    queryHumanNamesL = list()
    for qryGbk in allQryGbksL:
        basename = os.path.basename(qryGbk)
        humanName = humanMapD[basename]
        queryHumanNamesL.append(humanName)

    # concatenate the alignments
    __concatenateAlignments(queryHumanNamesL,
                            speTreWorkDir,
                            catAlnFN,
                            keyFN,
                            mapFN)

    # print the summary
    __printSummary(paramO)

    # save more detailed core gene summary file
    __saveCoreGenesDetails(paramO)

    # outgroup will not be the root if phantasm identified reference genomes
    if not outgroup.isRoot():
        # determine the outgroup's name
        outgroupTaxonName = _makeTaxonName(outgroup)
    
    # the outgroup will be the root if user is specifying their own genomes
    else:
        outgroupTaxonName = outgroup.ncbiName

    # use fast tree if requested
    if program == FAST_TREE:
        print(PRINT_2A, end='', flush=True)
        logger.info(PRINT_2A)
        __runFastTree(paramO)
    
    # use iqtree if requested
    elif program == IQTREE:
        logger.info(PRINT_2B)
        print(PRINT_2B, end='', flush=True)
        __runIqTree(outgroupTaxonName, paramO)
    
    # raise an error if an unexpected program was requested
    else:
        logger.error(ERR_MSG)
        raise ValueError(ERR_MSG)

    # root the tree on the specified outgroup
    __rootTree(speTreeFN, [outgroupTaxonName])

    print(DONE)
    logger.info(DONE)


def __makeGeneTreesWrapper(paramO:Parameters) -> None:
    """ makeGeneTreeWrapper:
            Accepts a Parameters object as input. Uses xenoGI to make a tree
            for each core gene identified. Does not return.
    """
    # extract necessary data from paramO
    aabrhFN = paramO.aabrhFN
    workDir = paramO.makeSpeciesTreeWorkingDir
    gtFileStem = paramO.aabrhHardCoreGeneTreeFileStem

    # make parameter dictionary for xenoGI
    paramD = paramO.toDict()

    # load genesO and aarbHardCoreL using xenoGI
    genesO = xenoGI.xenoGI.loadGenomeRelatedData(paramD)[1]
    aabrhHardCoreL = xenoGI.scores.loadOrthos(aabrhFN)

    ## set up
    # create work dir if it doesn't already exist
    if glob.glob(workDir)==[]:
        os.mkdir(workDir)

    # get the pattern for identifying gene trees
    allGtFilePath = os.path.join(workDir,gtFileStem+'*.tre')
   
    # delete any pre-existing hard core gene trees
    for fn in glob.glob(allGtFilePath):
        os.remove(fn)
    
    ## make gene tree for each aabrh hard Core set
    # add numbering to list
    newAabrhHardCoreL = []
    for orthoNum,orthoT in enumerate(aabrhHardCoreL):
        newAabrhHardCoreL.append((orthoNum,orthoT))
    
    # make the gene trees
    xenoGI.trees.makeGeneTreesFastTree(paramD,
                                       True,
                                       genesO,
                                       workDir,
                                       gtFileStem,
                                       newAabrhHardCoreL)


def __alignGenesWrapper(paramO:Parameters) -> None:
    """ alignGenesWrapper:
            Accepts a Parameters object as input. Aligns core genes in parallel
            using MUSCLE via xenoGI. Does not return.
    """
    # extract necessary data from paramO
    aabrhFN = paramO.aabrhFN
    workDir = paramO.makeSpeciesTreeWorkingDir
    gtFileStem = paramO.aabrhHardCoreGeneTreeFileStem
    musclePath = paramO.musclePath
    numProcesses = paramO.numProcesses

    # make parameter dictionary for xenoGI
    paramD = paramO.toDict()

    # get the alignStem
    alignStem = "align"+"-"+gtFileStem

    # load genesO and aarbHardCoreL using xenoGI
    genesO = xenoGI.xenoGI.loadGenomeRelatedData(paramD)[1]
    aabrhHardCoreL = xenoGI.scores.loadOrthos(aabrhFN)

    ## set up
    # create work dir if it doesn't already exist
    if glob.glob(workDir)==[]:
        os.mkdir(workDir)

    newAabrhHardCoreL = []
    for orthoNum,orthoT in enumerate(aabrhHardCoreL):
        newAabrhHardCoreL.append((orthoNum,orthoT))

    xenoGI.trees.createAlignments(paramD,
                                  alignStem,
                                  True,
                                  genesO,
                                  workDir,
                                  musclePath,
                                  numProcesses,
                                  newAabrhHardCoreL)


def __concatenateAlignments(qryHumanNamesL:list, speciesTreeWorkDir:str, \
                                alnOutFN:str, keyFN:str, wgsMapFN:str) -> None:
    """ concatenateAlignments:
            Accepts a list of human names for the query genomes, a string indi-
            cating the make species tree working directory, a string indicating
            the output filename, a string indicating the famToGeneKey filename,
            and a string indicating the wgsHumanMap filename as inputs. Concat-
            enates the alignments of all the hardcore genes and then saves the
            concatenated alignment in fasta format. Does not return.
    """
    # constants
    from param import REDUCE_NUM_CORE_GENES
    FILE_NAME_PATTERN = 'align*afa'
    FORMAT = 'fasta'
    GREP_FIND_1 = r'^.+align-aabrhHardCoreFam-(\d+)\.afa$'
    GREP_FIND_2 = r'^\S+ (\d+)$'
    GREP_REPL = r'\1'
    DELIM = "\t"
    EOL = "\n"
    ERR_MSG = "untested condition; please report this at https://github.com/dr-joe-wirth/phantasm"

    logger = logging.getLogger(__name__ + "." + __concatenateAlignments.__name__)

    # get the files sring
    fileString = os.path.join(speciesTreeWorkDir, FILE_NAME_PATTERN)

    # get a list of all the alignment files
    allAfa = glob.glob(fileString)

    if REDUCE_NUM_CORE_GENES:
        __reduceCoreGenes(allAfa)

    # concatenate sequences in dictionary format
    seqD = dict()
    famKeyD = dict()
    for filename in allAfa:
        # get the aabrhHardCoreFam number from the filename
        famNum = str(int(re.sub(GREP_FIND_1, GREP_REPL, filename)))

        # parse the sequence file
        parsed = SeqIO.parse(filename, FORMAT)

        # for each entry in the file
        record:SeqRecord
        for record in parsed:
            # make a new entry for each species that is absent
            if record.id not in seqD.keys():
                seqD[record.id] = record.seq
            
            # concatenate sequences for species that are present
            else:
                seqD[record.id] += record.seq
            
            # save the gene number for the user's input genomes
            if record.id in qryHumanNamesL:
                # get the gene number and save it in the number file
                genNum = re.sub(GREP_FIND_2, GREP_REPL, record.description)

                # make a new dictionary for each famNum that is absent
                if famNum not in famKeyD.keys():
                    famKeyD[famNum] = {record.id: genNum}

                # add subsequent entries to nested dictionary
                else:
                    famKeyD[famNum].update({record.id: genNum})

    # make sure keyFN is empty
    keyFH = open(keyFN, 'w')
    keyFH.close()

    # open keyFN to write famKeyD data to file
    keyFH = open(keyFN, 'a')

    # for each family number
    for famNum in famKeyD.keys():
        # initialize the line with the family number
        lineStr = famNum + DELIM

        # loop through the query names so that ordering is consistent
        for queryName in qryHumanNamesL:
            # add the associated gene numbers separated by the SEP character
            lineStr += famKeyD[famNum][queryName] + SEP_CHAR
        
        # remove the trailing SEP character and add EOL to the line
        lineStr = lineStr[:-1] + EOL

        # write the line to the file
        keyFH.write(lineStr)
    
    # close the file
    keyFH.close()

    # load the wgsHumanMap file into memory
    mapD = loadHumanMap(wgsMapFN)

    # convert the dictionary to a list of SeqRecord objects
    allConcatenatedRecords = list()
    for key in seqD.keys():
        rec = SeqRecord(seqD[key])

        # get a list of all the taxa containing the key
        taxaL = [x for x in mapD.values() if key in x]
        
        # if only one taxon matched, then use it
        if len(taxaL) == 1:
            rec.id = taxaL.pop()
        
        # if multiple taxa matched, then find the proper taxon
        else:
            # go through each taxon and see if the key matches
            found = False
            for taxonName in taxaL:
                if key == taxonName:
                    # use the name that matched and stop looping
                    rec.id = taxonName
                    found = True
                    break
            
            # if a key couldn't be found, then raise an error
            if not found:
                logger.error(ERR_MSG)
                raise RuntimeError(ERR_MSG)
            
        # description field should be empty
        rec.description = ""
        allConcatenatedRecords.append(rec)
    
    # write the file in fasta format
    SeqIO.write(allConcatenatedRecords, alnOutFN, FORMAT)


def __reduceCoreGenes(allAlignmentsL:list) -> int:
    """ reduceCoreGenes:
            Accepts a list of alignment filenames in fasta format as input.
            Removes any alignment files from the list that contain greater than
            5% gaps for one or more taxa. Does not return; only modifies the
            input list.
    """
    # get a list of indices in reverse order for on-the-fly popping
    indices = list(range(len(allAlignmentsL)))
    indices.reverse()
    
    # for each index in the list
    for idx in indices:
        # get the associated alignment
        alnFN = allAlignmentsL[idx]

        # check if it should be kept; remove it if necessary
        if not __lessThanFivePercentGaps(alnFN):
            allAlignmentsL.pop(idx)


def __lessThanFivePercentGaps(alignmentFN:str) -> bool:
    """ lessThanFivePercentGaps:
            Accepts the filename of an alignment file in fasta format as input.
            Returns True if the all taxa in the alignment have less than or eq-
            ual to 5% gaps relative to the alignment length. Otherwise returns
            False.
    """
    # constants
    FORMAT = 'fasta'
    INITIAL_SEQ_LEN = -1
    GAP_CHAR = "-"
    MAX_PERC = 0.05
    ERR_MSG = "alignment length inconsistent for "
    
    logger = logging.getLogger(__name__ + "." + __lessThanFivePercentGaps.__name__)

    # parse the file
    parsed = SeqIO.parse(alignmentFN, FORMAT)

    # initialize the sequence length to an impossible value
    seqLen = INITIAL_SEQ_LEN

    # for each sequence record
    rec:SeqRecord
    for rec in parsed:
        # get the new sequence length
        newSeqLen = len(rec.seq)
        # evaluate if the sequence length is changing
        notFirstRecord = seqLen != INITIAL_SEQ_LEN
        seqLenChanged = newSeqLen != seqLen
    
        # raise an error if seq length changed and it's not the first record
        if notFirstRecord and seqLenChanged:
            logger.error(ERR_MSG + alignmentFN)
            raise RuntimeError(ERR_MSG + alignmentFN)
        
        seqLen = newSeqLen

        # count the number of gaps in each character of the sequence
        numGaps = 0
        for char in rec.seq:
            if char == GAP_CHAR:
                numGaps += 1

        # if a single seq has > than max allowed gaps, then do not keep it
        if numGaps / seqLen > MAX_PERC:
            return False
    
    # only keep an alignment if all sequences do not exceed max allowed gaps
    return True


def __printSummary(paramO:Parameters) -> None:
    """ printSummary:
            Accepts a Parameters object as input. Prints the number of core ge-
            nes used to generate the alignment and and the sequence length of
            each concatenated alignment. Does not return.
    """
    # constants
    GAP = " " * 4
    PRINT_1 = GAP + "Number of core genes used to construct the tree: "
    PRINT_2 = GAP + "Sequence length for each concatenated alignment: "
    PRINT_3A = GAP*2 + "See '"
    PRINT_3B = "' for more details."

    # extract data from the alignment
    numCoreGenes, lenAlignment = __getSummary(paramO)

    # print the data
    print(PRINT_1 + str(numCoreGenes))
    print(PRINT_2 + str(lenAlignment))

    # determine where the full results will be saved
    workdir = os.path.basename(paramO.workdir)
    deetsFN = os.path.basename(paramO.coreGenesSummaryFN)
    deetsFN = os.path.join(workdir, deetsFN)

    print(PRINT_3A + deetsFN + PRINT_3B)


def __getSummary(paramO:Parameters) -> tuple[int,int]:
    """ getSummary:
            Accepts a Parameters object as input. Uses existing files to deter-
            mine the number of core genes used to construct the concatenated a-
            lignment and the sequence length of each individual alignment. Ret-
            urns the values as a tuple.
    """
    # get the filenames for the famToGeneKey and the concatated alignment files
    famKeyFN = paramO.famToGeneKeyFN
    concatAlnFN = paramO.concatenatedAlignmentFN
    
    # get the number of core genes (nrow of coreGenesFN)
    parsed = parseCsv(famKeyFN, delim="\t")
    numCoreGenes = len(parsed)

    # parse the concatenated alignment
    parsed = SeqIO.parse(concatAlnFN, 'fasta')

    # begin iterating through the records
    record:SeqRecord
    for record in parsed:
        # get the length of the first aligned sequence in the file
        lenAlignment = len(record.seq)

        # all other sequences will be the same length, so no need to check them
        break
    
    return numCoreGenes, lenAlignment


def __saveCoreGenesDetails(paramO:Parameters) -> None:
    """ saveCoreGenesDetails:
            Accepts a Parameters object as input. Extracts relevant data about
            the core genes and saves them to the file specified in the Parame-
            ters object. Does not return.
    """
    # constants
    DELIM = "\t"
    EOL = "\n"
    FORMAT = "fasta"
    ALIGN_PRE = "align-aabrhHardCoreFam-"
    ALIGN_EXT = ".afa"
    HEADER_STR = "alignment_file" + DELIM + "alignment_len" + DELIM + \
                 "locus_tag" + DELIM + "protein_len" + DELIM + "gene_name" + \
                 DELIM + "annotation" + DELIM + "gene_number" + EOL
    TAIL_1 = "Number of core genes used to construct the tree: "
    TAIL_2 = "Sequence length for each concatenated alignment: "

    # extract relevant data from paramO
    workdir = paramO.workdir
    famToGeneKeyFN = paramO.famToGeneKeyFN
    geneInfoFN = paramO.geneInfoFN
    alnDir = paramO.makeSpeciesTreeWorkingDir
    summaryFN = paramO.coreGenesSummaryFN

    # the famToGene dict has famNum key with a list of geneNum as value
    #### Most importantly, it only contains entries for the genes to construct
    #### the concatenated alignment. This allows it function properly even if
    #### REDUCE_NUM_CORE_GENES is True
    famToGeneD = __loadFamToGeneKey(famToGeneKeyFN)

    # load genesO object
    genesO = xenoGI.genomes.genes(geneInfoFN)
    genesO.initializeGeneInfoD(geneInfoFN)

    # open the output file
    summaryFH = open(summaryFN, 'w')

    # write the header string
    summaryFH.write(HEADER_STR)

    # for each family number
    for famNum in famToGeneD.keys():
        # get the alignment filename
        # the number needs to be six digits long; add zeros to make it so
        alnNum = "0" * (6 - len(famNum)) + famNum
        basename = ALIGN_PRE + alnNum + ALIGN_EXT

        # determine the absolute path to the alignment file
        alignFN = os.path.join(alnDir, basename)

        # determine the abbreviated alignment file name
        work = os.path.basename(workdir)
        aln = os.path.basename(alnDir)
        alignFileStr = os.path.join(work, aln, basename)

        # get the length of the alignment for this gene family
        parsed = SeqIO.parse(alignFN, FORMAT)
        record:SeqRecord
        for record in parsed:
            # the length is the same for all records in the file
            alignmentLen = len(record.seq)
            break
    
        # start making the string for the current line
        lineStr = alignFileStr + DELIM + str(alignmentLen) + DELIM

        # get the data for the gene numbers and add them to the line string
        geneNumL = famToGeneD[famNum]
        lineStr += __geneNumToSummaryStr(geneNumL, genesO)

        # write the line to the summary file
        summaryFH.write(lineStr)
    
    # get the basic details about the concatenated alignment
    numCoreGenes, lenConcatAln = __getSummary(paramO)
    
    # construct the final lines of the file
    # matches the string produced by __printSummary
    tail = EOL + EOL + TAIL_1 + str(numCoreGenes) + EOL + \
                                               TAIL_2 + str(lenConcatAln) + EOL

    # write the tail string then close the file
    summaryFH.write(tail)
    summaryFH.close()


def __runFastTree(paramO:Parameters) -> None:
    """ runFastTree:
            Accepts a Parameters object as input. Calls FastTree on the concat-
            enated alignment. Does not return.
    """
    # constants
    ENV_VAR = "OMP_NUM_THREADS"

    # extract the necessary data from paramO
    fastTree = paramO.fastTreePath
    speTreeFN = paramO.speciesTreeFN
    catAlnFN = paramO.concatenatedAlignmentFN
    numThreads = paramO.numProcesses

    # set the number of threads for fasttree
    os.putenv(ENV_VAR, str(numThreads))

    # run fasttree on the concatenated alignment
    cmd = [fastTree, "-quiet", "-out", speTreeFN, catAlnFN]
    subprocess.run(cmd)


def __runIqTree(outgroupName:str, paramO:Parameters) -> None:
    """ runIqTree:
            Accepts an outgroup name (str) and a Parameters object as inputs.
            Calls IQTree on the concatenated alignment with the specified numb-
            er of bootstrap supports. Does not return.
    """
    # constants
    IQTREE_OUT_EXT = ".treefile"

    # extract necessary data from paramO
    alnFN = paramO.concatenatedAlignmentFN
    numBootstraps = str(paramO.numBootstraps)
    iqtree = paramO.iqTreePath
    numThreads = str(paramO.numProcesses)
    humanMapFN = paramO.fileNameMapFN
    speTreeFN = paramO.speciesTreeFN

    # reformat the outgroup name
    outgroupName = __formatNamesForIqTree(outgroupName)

    # build iqtree command
    cmd = [iqtree, '-quiet', \
           '-b', numBootstraps, \
           '-nt', numThreads, \
           '-o', outgroupName, \
           '-s', alnFN]
    
    # run iqtree
    subprocess.run(cmd)
    
    # determine the output filename
    iqtreeOutFN = alnFN + IQTREE_OUT_EXT
    
    # load the tree file as a string
    iqtreeFH = open(iqtreeOutFN, 'r')
    treeStr = iqtreeFH.read()
    iqtreeFH.close()

    # load the human map file
    humanMapD = loadHumanMap(humanMapFN)

    # use the human map file to replace the names in the tree
    for origName in humanMapD.values():
        treeName = __formatNamesForIqTree(origName)
        treeStr = re.sub(treeName, origName, treeStr)
    
    # write the tree string to the species tree file
    treeFH = open(speTreeFN, 'w')
    treeFH.write(treeStr)
    treeFH.close()


def __rootTree(treeFN:str, outGroupTaxaL:list) -> None:
    """ rootTree:
            Accepts a tree file name and a list of outgroup taxa as inputs. Lo-
            ads the tree and roots it on the specified taxa. Overwites the the
            original file with the rooted tree. Does not return.
    """
    # constants
    FORMAT = 'newick'
    OUTGROUP_IDX = -1
    INGROUP_IDX = 1

    # read tree
    tree:Newick.Tree = Phylo.read(treeFN, FORMAT)

    # root tree
    tree.root_with_outgroup(outGroupTaxaL)

    # get a list of clades
    cladesL = list(tree.find_clades())

    # get the outgroup and ingroup clades
    outgroup:Newick.Clade = cladesL[OUTGROUP_IDX]
    ingroup:Newick.Clade = cladesL[INGROUP_IDX]

    # get the total distance between the two clades
    totalDist = outgroup.branch_length + ingroup.branch_length
    
    # evenly distribute the distance between the two nodes
    outgroup.branch_length = totalDist / 2
    ingroup.branch_length  = totalDist / 2

    # write tree
    Phylo.write(tree, treeFN, FORMAT)


def __formatNamesForIqTree(taxonName:str) -> str:
    """ formatNamesForIqtree:
            Accepts a taxon name (str) as input. Replaces all disallowed chara-
            cters with underscores in order to determine the IQTree will use
            for the specified taxon name. Returns the modified name.
    """
    # constants
    ALLOWED_CHARS = string.ascii_letters + string.digits + "_-"

    # initialize the new name
    newName = ""

    # for each character in the taxon name
    for char in taxonName:
        # keep the character if it is allowed
        if char in ALLOWED_CHARS:
            newName += char
        
        # otherwise replace the character with an underscore
        else:
            newName += '_'
    
    return newName


def rankPhylogeneticMarkers(paramO:Parameters) -> None:
    """ rankPhylogeneticMarkers:
            Accepts a Parameters object as input. Calculates the cophenetic co-
            rrelation coefficient for all hardcore genes and writes the results
            to the file specified in the Parameters object. Does not return.
    """
    # constants
    PRINT_1 = 'Ranking phylogenetic markers ... '
    DONE = 'Done.'
    GENE_NUM_IDX = 0
    COPH_COR_IDX = 1 
    
    DELIM = "\t"
    EOL = "\n"
    HEADER_STRING = "cophenetic_corr_coef" + DELIM + "locus_tag" + DELIM + \
                    "protein_len" + DELIM + "gene_name" + DELIM + \
                    "annotation" + DELIM + "gene_number" + EOL

    logger = logging.getLogger(__name__ + "." + rankPhylogeneticMarkers.__name__)

    # print status
    print(PRINT_1, end='', flush=True)
    logger.info(PRINT_1)
    
    # extract necessary data from paramO
    geneInfoFN = paramO.geneInfoFN
    phyloMarkersFN = paramO.phyloMarkersFN

    # load genesO
    genesO = xenoGI.genomes.genes(geneInfoFN)
    genesO.initializeGeneInfoD(geneInfoFN)

    # obtain the cophenetic correlation coefficients for each core gene
    cophCorrCoefs = __calculateCopheneticCorrelations(paramO)

    # make sure the file is already empty
    filehandle = open(phyloMarkersFN, "w")
    filehandle.close()

    # save the results to the specified file
    filehandle = open(phyloMarkersFN, 'a')
    filehandle.write(HEADER_STRING)
    entry:tuple
    for entry in cophCorrCoefs:
        # parse the tuple into gene number and cophenetic correlation coef.
        geneNumStr:str = entry[GENE_NUM_IDX]
        cophCor = str(entry[COPH_COR_IDX])

        lineStr = __geneNumToSummaryStr(geneNumStr.split(SEP_CHAR), genesO)
        lineStr = cophCor + DELIM + lineStr
        filehandle.write(lineStr)
    
    # close the file
    filehandle.close()
    print(DONE)
    logger.info(DONE)


def __calculateCopheneticCorrelations(paramO:Parameters) -> list:
    """ calculateCopheneticCorrelations:
            Accepts a Parameters object as input. Calculates the cophenetic co-
            rrelation coefficient for each hardcore gene and constructs a list
            of tuples where the first item is the gene number and the next item
            is the cophenetic correlation coefficient. Returns the list with
            the items sorted from highest correlation to lowest.
    """
    # constants
    TREE_EXT = ".tre"

    # extract values from paramO
    geneTreeFileStem = paramO.aabrhHardCoreGeneTreeFileStem
    speciesTreeWorkDir = paramO.makeSpeciesTreeWorkingDir
    speciesTreeFN = paramO.speciesTreeFN
    famToGeneKeyFN = paramO.famToGeneKeyFN

    # get the distance matrix of the species tree
    speDistMat = __distanceMatrixFromNewickFile(speciesTreeFN)

    # load the fam to gene number key as a dictionary
    famToGeneKey = __loadFamToGeneKey(famToGeneKeyFN)

    # for each gene tree in the famToGeneKey
    result = dict()
    for famNum in famToGeneKey.keys():
        # the gene tree number will be exactly six char long; fill out with 0
        geneTreeNum = '0' * (6 - len(famNum)) + famNum

        # add the filestem and its path to the gene number to get the filename
        geneTreeFN = geneTreeFileStem + geneTreeNum + TREE_EXT
        geneTreeFN = os.path.join(speciesTreeWorkDir, geneTreeFN)

        # get the distance matrix of the gene tree
        genDistMat = __distanceMatrixFromNewickFile(geneTreeFN)

        # get the list of gene numbers for this gene family
        geneNumL = famToGeneKey[famNum]

        # convert the list to a string of separated numbers
        geneNumStr = ""
        for geneNum in geneNumL:
            geneNumStr += geneNum + SEP_CHAR
        
        # truncate trailing sep character
        geneNumStr = geneNumStr[:-len(SEP_CHAR)]

        # calculate the cophenetic correlation coefficient and save the result
        result[geneNumStr] = __cophenetic(speDistMat, genDistMat)
    
    # sort the output from highest value to lowest value
    result = sorted(result.items(), key=lambda x: x[1], reverse=True)

    return result


def __geneNumToSummaryStr(geneNumL:list, genesO:xenoGI.genomes.genes) -> str:
    """ geneNumToSummaryStr:
            Accepts a list of gene numbers (int or str) for a single core gene
            and a genes object as inputs. Extracts data from genesO for each
            gene number, and constructs a string with all the information that
            can be written to a file. Returns the newly constructed string.
    """
    # constants
    GENE_NAME_IDX = 1
    LOCUS_TAG_IDX = 2
    ANNOTATXN_IDX = 4
    CDS_START_IDX = 6
    CDS_END_IDX = 7
    ANN_SEP_CHAR = " | "
    PRO_SEP_CHAR = "|"
    DELIM = "\t"
    EOL = "\n"

    # initialize empty strings for each field
    locusTag = ""
    geneName = ""
    annotation = ""
    proteinLen = ""
    geneNumStr = ""

    # for each gene number 
    for geneNum in geneNumL:
        # add the gene number to the gene num string
        geneNumStr += str(geneNum) + SEP_CHAR

        # get a tuple of all the gene information
        geneInfo = genesO.numToGeneInfo(int(geneNum))

        # parse the desired strings from the geneInfo dictionary
        # append the data onto the growing strings
        # multiple genomes will have multiple values at each position
        locusTag += geneInfo[LOCUS_TAG_IDX] + SEP_CHAR
        geneName += geneInfo[GENE_NAME_IDX] + SEP_CHAR
        annotation += geneInfo[ANNOTATXN_IDX] + ANN_SEP_CHAR

        # calculate the length of the protein
        cdsStart = int(geneInfo[CDS_START_IDX])
        cdsEnd = int(geneInfo[CDS_END_IDX])
        proteinLen += str(abs(cdsEnd - cdsStart) // 3) + PRO_SEP_CHAR

    # remove the trailing sep character(s)
    locusTag = locusTag[:-len(SEP_CHAR)]
    geneName = geneName[:-len(SEP_CHAR)]
    annotation = annotation[:-len(ANN_SEP_CHAR)]
    proteinLen = proteinLen[:-len(PRO_SEP_CHAR)]
    geneNumStr = geneNumStr[:-len(SEP_CHAR)]

    # construct the string and return it
    return locusTag + DELIM + proteinLen + DELIM + geneName + DELIM + \
                                          annotation + DELIM + geneNumStr + EOL


def __cophenetic(speciesDistMat:dict, geneDistMat:dict) -> float:
    """ cophenetic:
            Accepts two distance matrices in dictionary format as inputs whose
            keys are identical. The keys are expected to be pairs of tips as a
            tuple and the values are expected to be the distance between the
            two tips in the tuple. Calculates and returns the cophenetic corre-
            lation coefficient based on the input matrices.
    """
    # constant
    ERR_MSG = "The distance matrices must use identical names."
    
    logger = logging.getLogger(__name__ + "." + __cophenetic.__name__)
    
    # make sure the tip names match
    if speciesDistMat.keys() != geneDistMat.keys():
        logger.error(ERR_MSG)
        raise ValueError(ERR_MSG)

    # flatten the matrices into 1d arrays that are linked by their indices
    speDistArray = list()
    genDistArray = list()
    for key in speciesDistMat.keys():
        speDistArray.append(speciesDistMat[key])
        genDistArray.append(geneDistMat[key])
    
    # calculate the pearson correlation coeffient and its p-value
    r, p = scipy.stats.pearsonr(speDistArray, genDistArray)

    # return the cophenetic correlation coefficient
    return r


def __distanceMatrixFromNewickFile(treeFN:str) -> dict:
    """ distanceMatrixFromNewickFile:
            Accepts a string indicating the filename to a tree file as input.
            Calculates a distance matrix for all pairwise combinations of the
            tips as dictionary where the keys are tuples containing a pair of
            tips and the values are the cophenetic distances between the pair
            of tips. Returns the dictionary.
    """
    # import the file as a Utree object
    tree = xenoGI.Tree.Utree()
    tree.fromNewickFile(treeFN)

    # make and return the distance matrix in dictionary format
    return tree.makeDistanceMatrix()


def __loadFamToGeneKey(famToGeneKeyFN:str) -> dict:
    """ loadFamToGeneKey:
            Accepts a string indicating a filename as input. Expects the input
            file to be tab-delimited and contain exactly two columns where the
            first column is the aabrhHardCore family number and the second col-
            umn is the gene number in the user's input. Loads the file into a 
            dictionary where the family number is the key and the values are
            the gene numbers. Returns the dictionary.
    """
    # constants
    DELIM = "\t"
    FAM_NUM_IDX = 0
    GENE_NUM_IDX = 1

    # load the file as a list of lists
    raw = parseCsv(famToGeneKeyFN, delim=DELIM)

    # for each item, make a dictionary keyed by fam number
    outD = dict()
    for line in raw:
        outD[line[FAM_NUM_IDX]] = line[GENE_NUM_IDX].split(SEP_CHAR)
    
    return outD

