# Author: Joseph S. Wirth

from __future__ import annotations
import sys, os, re, glob, scipy.stats, csv, subprocess
from PHANTASM.utilities import parseCsv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from PHANTASM.downloadGbff import downloadGbffsForRootTaxonomy, _makeHumanMapString, _makeTaxonName
from PHANTASM.taxonomy.Taxonomy import Taxonomy
from param import XENOGI_DIR
sys.path.insert(0,os.path.join(sys.path[0], XENOGI_DIR))
import xenoGI.xenoGI, xenoGI.scores, xenoGI.Tree, xenoGI.genomes, xenoGI.trees


def xenogiInterfacer_1(taxO:Taxonomy, queryGbff:str, paramD:dict) -> Taxonomy:
    """ xenogiInterfacer_1:
            Accepts a Taxonomy object, a string indicating the path to the qu-
            ery genbank, and the parameter dictionary as inputs. Downloads the
            ingroup and outgroup sequences and puts them into a file structure
            that is expected by xenoGI. Creates the human map file required by 
            xenoGI. Returns the outgroup species as a Taxonomy object
    """
    # constants
    USER_INPUT = "user_input"

    # extract relevant data from paramD
    gbffFN = paramD['genbankFilePath']
    humanMapFN = paramD['fileNameMapFN']
    taxObjFilePath = paramD['taxonomyObjectFilePath']
    
    # determine the directory for the gbff file
    gbffDir = os.path.dirname(gbffFN)

    # identify and download gbff files from NCBI
    outgroup:Taxonomy = downloadGbffsForRootTaxonomy(taxO, paramD)

    # determine paths to taxonomy files and the extensions
    oldTaxFN = glob.glob(taxObjFilePath).pop()
    dir = os.path.dirname(taxObjFilePath)
    ext = os.path.splitext(taxObjFilePath)[1]

    # save and replace existing taxonomy file
    taxO = taxO.getRoot()
    newTaxFN = os.path.join(dir, taxO.sciName + ext)
    os.remove(oldTaxFN)
    taxO.save(newTaxFN)

    # make a symlink to the user's input file
    oldFN = os.path.abspath(queryGbff)
    newFN = os.path.join(gbffDir, os.path.basename(queryGbff))
    os.symlink(oldFN, newFN)

    # add the query to the human map file
    humanMapStr = _makeHumanMapString(USER_INPUT, os.path.basename(queryGbff))
    filehandle = open(humanMapFN, "a")
    filehandle.write(humanMapStr)
    filehandle.close()

    return outgroup


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


def copyExistingBlastFiles(oldParamD:dict, newParamD:dict) -> None:
    """ copyExistingBlastFiles:
            Accepts two parameter dictionaries as inputs. Identifies any exist-
            ing blast tables that could be used in the new folder. Creates mod-
            ified files that are equivalent to the original but have replaced
            the old gene number with the new gene number. Does not return. 
    """
    # constants
    PRINT = "Copying existing blastp comparisons ... "
    DONE = 'Done.'

    # print job-start statement
    print(PRINT, end='', flush=True)

    # extract the new blast directory from newParamD
    newBlastDir = os.path.splitext(newParamD['blastFilePath'])[0][:-1]

    # create the new blast directory (if it doesn't already exist)
    if glob.glob(newBlastDir) == []:
        os.mkdir(newBlastDir)

    # determine which blastp tables need to be made
    filesToCopyD = __getBlastFilesToCopy(oldParamD, newParamD)

    # get a dictionary to link locus tags to their new names
    locusTagToNewNameD = __getLocusTagToNameD(newParamD)

    # for each file to copy
    for oldFilename in filesToCopyD:
        # get the new
        newFilename = filesToCopyD[oldFilename]

        # modify the query and subject names to match the new gene numbers
        # save the modified string to the new filename
        __updateQuerySubjectNames(oldFilename, newFilename, locusTagToNewNameD)
    
    print(DONE)


def __getBlastFilesToCopy(oldParamD, newParamD) -> dict:
    """ getBlastFilesToCopy:
            Accepts the old and new parameter dictionarys as inputs. Determines
            which blastp comparisons have already been performed. Creates a di-
            ctionary where the key is the path to the existing blastp files and
            the value is the path where the existing file should be copied. Re-
            turns the dictionary.
    """
    # constant
    DELIM = "\t"

    # extract relevant data from oldParamD
    oldBlastPath = oldParamD['blastFilePath'] # path/*

    # extract relevant data from newParamD
    blastJoin = newParamD['blastFileJoinStr']
    wgsMapFN = newParamD['fileNameMapFN']
    blastSplt = os.path.splitext(newParamD['blastFilePath'])

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


def __getLocusTagToNameD(paramD:dict) -> dict:
    """ getLocusTagToNameD:
            Accecpts a parameter dictionary as input. Uses genesO.geneInfoD to
            create a new dictionary keyed by locus tag with the gene blast str-
            ing as the value (eg. "999_Methanosarcina_acetivorans-MA_RS05270").
            Returns the newly created dictionary.
    """
    # constants
    LOCUS_TAG_IDX  = 2
    BLAST_NAME_IDX = 0

    # get the genesO object for the new folder and initialize its geneInfoD
    newGenesO = xenoGI.xenoGI.loadGenomeRelatedData(paramD)[1] # genesO object 
    newGenesO.initializeGeneInfoD(paramD['geneInfoFN'])

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
    GREP_FIND = r"^\d+_[^-]+-(\S+)$"
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
        locusTag_qry = re.sub(GREP_FIND, GREP_REPL, row[QRY_COL])
        locusTag_sbj = re.sub(GREP_FIND, GREP_REPL, row[SBJ_COL])

        # ... replace the query and subject names with the new names
        row[QRY_COL] = locusTagToNewNameD[locusTag_qry]
        row[SBJ_COL] = locusTagToNewNameD[locusTag_sbj]

        # ... write the row to the new file
        newWriter.writerow(row)
    
    # close the files
    oldHandle.close()
    newHandle.close()


def calculateCoreGenes(paramD) -> None:
    """ calculateCoreGenes:
            Accepts the parameter dictionary as input. Calculates the hardcore
            genes via the 'createAabrhL' function from xenoGI. Does not return.
    """
    # parse parameters into shorter variable names
    strainInfoFN = paramD['strainInfoFN']
    blastFileJoinStr = paramD['blastFileJoinStr']
    blastDir,blastExt = paramD['blastFilePath'].split("*")
    evalueThresh = paramD['evalueThresh']
    alignCoverThresh = paramD['alignCoverThresh']
    percIdentThresh = paramD['percIdentThresh']
    aabrhFN = paramD['aabrhFN']

    # read the strain info into a list
    strainNamesL = xenoGI.xenoGI.readStrainInfoFN(strainInfoFN)

    # determine the all-against-all best reciprocal hits (hardcore genes)
    print('Calculating core genes ... ', end='', flush=True)
    xenoGI.scores.createAabrhL(strainNamesL,
                               blastFileJoinStr,
                               blastDir,
                               blastExt,
                               evalueThresh,
                               alignCoverThresh,
                               percIdentThresh,
                               aabrhFN)
    print('Done.')


def makeSpeciesTree(paramD:dict, outgroup:Taxonomy) -> None:
    """ makeSpeciesTree:
            Accepts the parameters dictionary and an outgroup (Taxonomy) as in-
            puts. Uses xenoGI to make a species tree, and then builds another
            tree based on the concatenated alignments of the hardcore genes and
            roots the tree on the specified outgroup. Does not return.
    """
    # constants
    PRINT_1 = 'Aligning core genes and making gene trees ... '
    PRINT_2 = 'Making the species tree from a concatenated alignment ' + \
                                                       'of the core genes ... '
    DONE = 'Done.'

    # align hardcore genes and build gene trees with fasttree
    print(PRINT_1, end='', flush=True)
    __makeGeneTreesWrapper(paramD)
    print(DONE)
    
    # parse data from paramD
    speTreWorkDir = paramD['makeSpeciesTreeWorkingDir']
    catAlnFN = paramD['concatenatedAlignmentFN']
    keyFN = paramD['famToGeneKeyFN']
    wgsFN = paramD['fileNameMapFN']
    fastTree = paramD['fastTreePath']
    speTreeFN = paramD['speciesTreeFN']

    # concatenate the alignments
    __concatenateAlignments(speTreWorkDir, catAlnFN, keyFN, wgsFN)

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
    __printSummary(paramD)


def __makeGeneTreesWrapper(paramD) -> None:
    """ makeGeneTreeWrapper:
            Accepts the parameter dictionary as input. Uses xenoGI to make a
            tree for each core gene identified. Does not return.
    """
    # load genesO and aarbHardCoreL using xenoGI
    strainNamesT,genesO,geneOrderD = xenoGI.xenoGI.loadGenomeRelatedData(paramD)
    aabrhHardCoreL = xenoGI.scores.loadOrthos(paramD['aabrhFN'])

    ## set up
    # create work dir if it doesn't already exist
    workDir = paramD['makeSpeciesTreeWorkingDir']
    if glob.glob(workDir)==[]:
        os.mkdir(workDir)

    # get the pattern for identifying gene trees
    gtFileStem = paramD['aabrhHardCoreGeneTreeFileStem']
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
    xenoGI.trees.makeGeneTrees(paramD,True,genesO,workDir,gtFileStem,newAabrhHardCoreL)


def __concatenateAlignments(speciesTreeWorkDir:str, alnOutFN:str, keyFN:str, wgsMapFN:str) -> None:
    """ concatenateAlignments:
            Accepts a string indicating the make species tree working director-
            y, a string indicating the output filename, a string indicating the
            famToGeneKey filename, and a string indicating the wgsHumanMap fil-
            ename as inputs. Concatenates the alignments of all the hardcore
            genes and then saves the concatenated alignment in fasta format.
            Does not return.
    """
    # constants
    FILE_NAME_PATTERN = 'align*afa'
    FORMAT = 'fasta'
    GREP_FIND_1 = r'^.+align(\d+)\.afa$'
    GREP_FIND_2 = r'^\S+ (\d+)$'
    GREP_REPL = r'\1'
    USER_INPUT = "user_input"

    # get the files sring
    fileString = os.path.join(speciesTreeWorkDir, FILE_NAME_PATTERN)

    # get a list of all the alignment files
    allAfa = glob.glob(fileString)

    # make sure the fam number to gene number file is empty
    filehandle = open(keyFN, "w")
    filehandle.close()

    # open the file for appending data during the loop
    filehandle = open(keyFN, "a")

    # concatenate sequences in dictionary format
    seqD = dict()
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
            
            # save the gene number for the user's input
            if record.id == USER_INPUT:
                # get the gene number and save it in the number file
                genNum = re.sub(GREP_FIND_2, GREP_REPL, record.description)

                # make the tab delimited string and append it to the file
                filehandle.write(famNum + "\t" + genNum + "\n")

    # close the key file
    filehandle.close()

    # load the wgsHumanMap file into memory
    wgsMapL = parseCsv(wgsMapFN, delim="\t")

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


def __printSummary(paramD:dict) -> None:
    """ printSummary:
            Accepts the parameter dictionary as input. Prints the number of co-
            re genes used to generate the alignment and and the sequence length
            of each concatenated alignment. Does not return.
    """
    # constants
    GAP = " " * 4
    PRINT_1 = GAP + "Number of core genes used to construct the tree: "
    PRINT_2 = GAP + "Sequence length for each concatenated alignment: "

    # extract data on the alignment
    numCoreGenes, lenAlignment = __getSummary(paramD)

    # print the data
    print(PRINT_1 + str(numCoreGenes))
    print(PRINT_2 + str(lenAlignment))


def __getSummary(paramD:dict) -> tuple:
    """ getSummary:
            Accepts the parameter dictionary as input. Uses existing files to
            determine the number of core genes used to construct the concatena-
            ted alignment and the sequence length of each individual alignment.
            Returns the values as a tuple.
    """
    # get the filename for core genes file
    coreGenesFN = paramD["aabrhFN"]

    # te the filename for the concatenated alignment
    concatAlnFN = paramD["concatenatedAlignmentFN"]
    
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


def rankPhylogeneticMarkers(paramD:dict) -> None:
    """ rankPhylogeneticMarkers:
            Accepts the parameter dictionary as input. Calculates the cophenet-
            ic correlation coefficient for all hardcore genes and writes the
            results to the file specified in the parameter dictionary. Does not
            return.
    """
    # constants
    GENE_NUM_IDX = 0
    COPH_COR_IDX = 1 
    GENE_NAME_IDX = 1
    LOCUS_TAG_IDX = 2
    ANNOTATXN_IDX = 4
    DELIM = "\t"
    EOL = "\n"
    HEADER_STRING = "cophenetic_corr_coef" + DELIM + "gene_num" + DELIM + \
                    "locus_tag" + DELIM + "gene_name" + DELIM + "annotation" + EOL

    print('Ranking phylogenetic markers ... ', end='', flush=True)
    # load the gene info file as a dictionary
    geneInfoFN = paramD['geneInfoFN']
    genesO = xenoGI.genomes.genes(geneInfoFN)
    genesO.initializeGeneInfoD(geneInfoFN)

    # obtain the cophenetic correlation coefficients for each core gene
    cophCorrCoefs = __calculateCopheneticCorrelations(paramD)

    # make sure the file is already empty
    phyloMarkersFN = paramD['phyloMarkersFN']
    filehandle = open(phyloMarkersFN, "w")
    filehandle.close()

    # save the results to the specified file
    filehandle = open(phyloMarkersFN, 'a')
    filehandle.write(HEADER_STRING)
    entry:tuple
    for entry in cophCorrCoefs:
        # parse the tuple into gene number and cophenetic correlation coef.
        geneNum = entry[GENE_NUM_IDX]
        cophCor = str(entry[COPH_COR_IDX])

        # get a tuple of all the gene information
        geneInfo = genesO.numToGeneInfo(int(geneNum))

        # parse the desired strings from the geneInfo dictionary
        locusTag = geneInfo[LOCUS_TAG_IDX]
        geneName = geneInfo[GENE_NAME_IDX]
        annotation = geneInfo[ANNOTATXN_IDX]

        # construct the string and write to file
        lineStr = cophCor + DELIM + geneNum + DELIM + locusTag + DELIM + \
                                    geneName + DELIM + annotation + "\n"
        filehandle.write(lineStr)
    
    # close the file
    filehandle.close()
    print('Done.')


def __calculateCopheneticCorrelations(paramD:dict) -> list:
    """ calculateCopheneticCorrelations:
            Accepts the parameter dictionary as input. Calculates the cophenet-
            ic correlation coefficients for each hardcore gene and constructs a
            list of tuples where the first item is the gene number and the next
            item is the cophenetic correlation coefficient. Returns the list 
            with the items sorted from highest correlation to lowest.
    """
    # constants
    GREP_REPL = r"\1"
    GREP_FIND_PREFIX = r"^.+"
    GREP_FIND_SUFFIX = r"0{0,}(\d+)\.tre"

    # extract values from paramD
    geneTreeFileStem = paramD['aabrhHardCoreGeneTreeFileStem']
    speciesTreeWorkDir = paramD['makeSpeciesTreeWorkingDir']
    speciesTreeFN = paramD['speciesTreeFN']
    famToGeneKeyFN = paramD['famToGeneKeyFN']
    

    grepFind = GREP_FIND_PREFIX + geneTreeFileStem + GREP_FIND_SUFFIX

    filePattern = os.path.join(speciesTreeWorkDir, "*.tre")

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


