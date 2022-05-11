# Author: Joseph S. Wirth

from __future__ import annotations
import csv, glob, os, re, scipy.stats, string, subprocess, sys
from PHANTASM.utilities import parseCsv
from PHANTASM.Parameter import Parameters
from Bio import SeqIO
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
    PRINT = 'Parsing genbank files ... '
    DONE = 'Done.'

    # parse the genbank files
    print(PRINT, end='', flush=True)
    xenoGI.xenoGI.parseGenbankWrapper(paramD)
    print(DONE)


def allVsAllBlast(paramD:dict) -> None:
    """ allVsAllBlast:
            Accepts the parameter dictionary as input. Calls the 'runBlast'
            command of xenoGI. Does not return.
    """
    # constants
    PRINT = 'Running all pairwise blastp comparisons ... '
    DONE = 'Done.'

    # run all pairwise blast comparisons
    print(PRINT, end='', flush=True)
    xenoGI.xenoGI.runBlastWrapper(paramD)
    print(DONE)


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

    # print job-start statement
    print(PRINT, end='', flush=True)

    # extract the new blast directory from newParamD
    newBlastDir = os.path.splitext(newParamO.blastFilePath)[0][:-1]

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
        __updateQuerySubjectNames(oldFilename, newFilename, locusTagToNewNameD)
    
    print(DONE)


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

    # get the genesO object for the new folder and initialize its geneInfoD
    paramD = paramO.toDict()
    newGenesO = xenoGI.xenoGI.loadGenomeRelatedData(paramD)[1] # genesO object
    newGenesO.initializeGeneInfoD(paramO.geneInfoFN)

    # create a new dictionary: key = locus_tag, val = blast name
    locusTagToNameD = dict()
    for geneNum in newGenesO.geneInfoD.keys():
        locusTag  = newGenesO.geneInfoD[geneNum][LOCUS_TAG_IDX]
        blastName = newGenesO.geneInfoD[geneNum][BLAST_NAME_IDX]
        locusTagToNameD[locusTag] = blastName
    
    return locusTagToNameD


def __updateQuerySubjectNames(existingFile:str, newFile:str, \
                                              locusTagToNewNameD:dict) -> None:
    """ updateQuerySubjectNames:
            Accepts a path to an existing blast file, the path to save the mod-
            ified file, and a dictionary as inputs. Reads the blast file and
            changes the names of the query and subject so that their gene numb-
            ers are compatible with the current folder. Does not return.
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
        row[QRY_COL] = locusTagToNewNameD[locusTag_qry]
        row[SBJ_COL] = locusTagToNewNameD[locusTag_sbj]

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
    xenoGI.scores.createAabrhL(strainNamesL,
                               blastFileJoinStr,
                               blastDir,
                               blastExt,
                               evalueThresh,
                               alignCoverThresh,
                               percIdentThresh,
                               aabrhFN)
    print(DONE)


def makeSpeciesTree(allQryGbksL:list, paramO:Parameters, outgroup:Taxonomy) \
                                                                       -> None:
    """ makeSpeciesTree:
            Accepts a list of the query genbank files, a Parameters object and
            an outgroup (Taxonomy) as inputs. Uses xenoGI to make a species tr-
            ee, and then builds another tree based on the concatenated alignme-
            nts of the hardcore genes and roots the tree on the specified outg-
            roup. Does not return.
    """
    # constants
    PRINT_1 = 'Aligning core genes and making gene trees ... '
    PRINT_2 = 'Making the species tree from a concatenated alignment ' + \
                                                       'of the core genes ... '
    DONE = 'Done.'

    # align hardcore genes and build gene trees with fasttree
    print(PRINT_1, end='', flush=True)
    __makeGeneTreesWrapper(paramO)
    print(DONE)
    
    # parse data from paramO
    speTreWorkDir = paramO.makeSpeciesTreeWorkingDir
    catAlnFN = paramO.concatenatedAlignmentFN
    keyFN = paramO.famToGeneKeyFN
    wgsFN = paramO.fileNameMapFN
    fastTree = paramO.fastTreePath
    speTreeFN = paramO.speciesTreeFN

    # make a list of the human names for the query genomes
    queryHumanNamesL = list()
    for qryGbk in allQryGbksL:
        basename = os.path.basename(qryGbk)
        humanName = os.path.splitext(basename)[0]
        queryHumanNamesL.append(humanName)

    # concatenate the alignments
    __concatenateAlignments(queryHumanNamesL,
                            speTreWorkDir,
                            catAlnFN,
                            keyFN,
                            wgsFN)

    # run fasttree on the concatenated alignment
    print(PRINT_2, end='', flush=True)
    cmd = [fastTree, "-quiet", "-out", speTreeFN, catAlnFN]
    subprocess.run(cmd)

    # replace the invalid string for Taxonomy objects with that of trees
    outgroupTaxonName = _makeTaxonName(outgroup)

    # root the tree on the specified outgroup
    __rootTree(speTreeFN, [outgroupTaxonName])

    print(DONE)

    # print summary for the concatenated alignment
    __printSummary(paramO)


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
    xenoGI.trees.makeGeneTrees(paramD,
                               True,
                               genesO,
                               workDir,
                               gtFileStem,
                               newAabrhHardCoreL)


def __concatenateAlignments(qryHumanNamesL:list, speciesTreeWorkDir:str, \
                                 alnOutFN:str, keyFN:str,wgsMapFN:str) -> None:
    """ concatenateAlignments:
            Accepts a list of human names for the query genomes, a string indi-
            cating the make species tree working directory, a string indicating
            the output filename, a string indicating the famToGeneKey filename,
            and a string indicating the wgsHumanMap filename as inputs. Concat-
            enates the alignments of all the hardcore genes and then saves the
            concatenated alignment in fasta format. Does not return.
    """
    # constants
    FILE_NAME_PATTERN = 'align*afa'
    FORMAT = 'fasta'
    GREP_FIND_1 = r'^.+align(\d+)\.afa$'
    GREP_FIND_2 = r'^\S+ (\d+)$'
    GREP_REPL = r'\1'
    DELIM = "\t"
    EOL = "\n"

    # get the files sring
    fileString = os.path.join(speciesTreeWorkDir, FILE_NAME_PATTERN)

    # get a list of all the alignment files
    allAfa = glob.glob(fileString)

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
    wgsMapL = parseCsv(wgsMapFN, delim=DELIM)

    # convert the dictionary to a list of SeqRecord objects
    allConcatenatedRecords = list()
    for key in seqD.keys():
        rec = SeqRecord(seqD[key])

        # find the name that contains the key
        for row in wgsMapL:
            # taxon name is the second column
            taxonName = row[1]

            # if the key is within the taxon name
            if key in taxonName:
                # then use this as the id
                rec.id = taxonName
                break
            
        # description field should be empty
        rec.description = ""
        allConcatenatedRecords.append(rec)
    
    # write the file in fasta format
    SeqIO.write(allConcatenatedRecords, alnOutFN, FORMAT)


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

    # extract data on the alignment
    numCoreGenes, lenAlignment = __getSummary(paramO)

    # print the data
    print(PRINT_1 + str(numCoreGenes))
    print(PRINT_2 + str(lenAlignment))


def __getSummary(paramO:Parameters) -> tuple:
    """ getSummary:
            Accepts a Parameters object as input. Uses existing files to deter-
            mine the number of core genes used to construct the concatenated a-
            lignment and the sequence length of each individual alignment. Ret-
            urns the values as a tuple.
    """
    # get the filename for core genes file
    coreGenesFN = paramO.aabrhFN

    # get the filename for the concatenated alignment
    concatAlnFN = paramO.concatenatedAlignmentFN
    
    # get the number of core genes (nrow of coreGenesFN)
    parsed = parseCsv(coreGenesFN, delim="\t")
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


def __rootTree(treeFN, outGroupTaxaL) -> None:
    """ rootTree:
            Accepts a tree file name and a list of outgroup taxa as inputs. Lo-
            ads the tree and roots it on the specified taxa. Overwites the the
            original file with the rooted tree. Does not return.
    """
    # constants
    ROOT_NODE = 's0'

    # load the tree into memory as an Rtree object
    speciesRtreeO = xenoGI.Tree.Rtree()
    speciesRtreeO.fromNewickFileLoadSpeciesTree(treeFN, outGroupTaxaL, \
                                                             includeBrLen=True)

    # initialize the total distance and a list of keys that need to be revised
    totalRootDist = 0
    keysToChange = list()

    # for each branch length (keys are pairs of node names as a tuple) ...
    for key in speciesRtreeO.branchLenD.keys():
        # ... if the branch is directly connected to the root ...
        if ROOT_NODE in key:
            # ... update the distance and add the key to the list
            totalRootDist += speciesRtreeO.branchLenD[key]
            keysToChange.append(key)
    
    # determine what the new "to-root" distance should be
    newRootDist = totalRootDist / len(keysToChange)

    # for each key that needs to change, update with the "to-root" distance
    for key in keysToChange:
        speciesRtreeO.branchLenD[key] = newRootDist
    
    # convert the object to a newick string
    treeStr = speciesRtreeO.toNewickStr(includeBrLength=True)

    # ensure the newick string ends in a semicolon (important for heatmap!)
    if treeStr[-1] != ";":
        treeStr += ";"

    # write the rooted tree to file
    handle = open(treeFN, "w")
    handle.write(treeStr)
    handle.close()


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
    GENE_NAME_IDX = 1
    LOCUS_TAG_IDX = 2
    ANNOTATXN_IDX = 4
    CDS_START_IDX = 6
    CDS_END_IDX = 7
    ANN_SEP_CHAR = " | "
    PRO_SEP_CHAR = "|"
    DELIM = "\t"
    EOL = "\n"
    HEADER_STRING = "cophenetic_corr_coef" + DELIM + "protein_len" + DELIM + \
             "gene_num" + DELIM + "locus_tag" + DELIM + "gene_name" + DELIM + \
                                                             "annotation" + EOL

    # print status
    print(PRINT_1, end='', flush=True)

    # extract necessary data from paramO
    geneInfoFN = paramO.geneInfoFN
    phyloMarkersFN = paramO.phyloMarkersFN

    # load the gene info file as a dictionary
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

        # initialize empty strings for each field
        locusTag = ""
        geneName = ""
        annotation = ""
        proteinLen = ""

        # for each gene number 
        for geneNum in geneNumStr.split(SEP_CHAR):
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

        # construct the string and write it to file
        lineStr = cophCor + DELIM + proteinLen + DELIM + geneNumStr + DELIM + \
                         locusTag + DELIM + geneName + DELIM + annotation + EOL
        filehandle.write(lineStr)
    
    # close the file
    filehandle.close()
    print(DONE)


def __calculateCopheneticCorrelations(paramO:Parameters) -> list:
    """ calculateCopheneticCorrelations:
            Accepts a Parameters object as input. Calculates the cophenetic co-
            rrelation coefficient for each hardcore gene and constructs a list
            of tuples where the first item is the gene number and the next item
            is the cophenetic correlation coefficient. Returns the list with
            the items sorted from highest correlation to lowest.
    """
    # constants
    TREE_EXT = "*.tre"
    GREP_REPL = r"\1"
    GREP_FIND_PREFIX = r"^.+"
    GREP_FIND_SUFFIX = r"0{0,}(\d+)\.tre"

    # extract values from paramD
    geneTreeFileStem = paramO.aabrhHardCoreGeneTreeFileStem
    speciesTreeWorkDir = paramO.makeSpeciesTreeWorkingDir
    speciesTreeFN = paramO.speciesTreeFN
    famToGeneKeyFN = paramO.famToGeneKeyFN
    
    # construct the grep find pattern
    grepFind = GREP_FIND_PREFIX + geneTreeFileStem + GREP_FIND_SUFFIX

    # construct the pattern that matches all tree files
    filePattern = os.path.join(speciesTreeWorkDir, TREE_EXT)

    # get the distance matrix of the species tree
    speDistMat = __distanceMatrixFromNewickFile(speciesTreeFN)

    # load the fam to gene number key as a dictionary
    famToGeneKey = __loadFamToGeneKey(famToGeneKeyFN)

    # get a list of gene tree filenames
    allGeneTreeFiles = glob.glob(filePattern)

    # for each gene tree
    result = dict()
    for geneTreeFN in allGeneTreeFiles:
        # get the distance matrix of the gene tree
        genDistMat = __distanceMatrixFromNewickFile(geneTreeFN)

        # extract the fam number from the file name
        famNum = re.sub(grepFind, GREP_REPL, geneTreeFN)

        # replace the number with the gene number
        geneNum = famToGeneKey[famNum]

        # calculate the cophenetic correlation coefficient and save the result
        result[geneNum] = __cophenetic(speDistMat, genDistMat)
    
    # sort the output from highest value to lowest value
    result = sorted(result.items(), key=lambda x: x[1], reverse=True)

    return result


def __cophenetic(speciesDistMat:dict, geneDistMat:dict) -> float:
    """ cophenetic:
            Accepts two distance matrices in dictionary format as inputs whose
            keys are identical. The keys are expected to be pairs of tips as a
            tuple and the values are expected to be the distance between the
            two tips in the tuple. Calculates and returns the cophenetic corre-
            lation coefficient based on the input matrices.
    """
    # make sure the tip names match
    if speciesDistMat.keys() != geneDistMat.keys():
        raise ValueError("The distance matrices must use identical names.")

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
        outD[line[FAM_NUM_IDX]] = line[GENE_NUM_IDX]
    
    return outD

