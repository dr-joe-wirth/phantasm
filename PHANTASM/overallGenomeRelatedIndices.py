import glob, math, os, sys
import rpy2.robjects as robjects
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from param import XENOGI_DIR, PHANTASM_DIR
sys.path.insert(0,os.path.join(sys.path[0],XENOGI_DIR))
import xenoGI.xenoGI, xenoGI.analysis
from PHANTASM.rRNA.runRnaBlast import __makeblastdbFromFna, __makeOutfmtString
from PHANTASM.utilities import parseCsv


###############################################################################
def overallGenomeRelatedIndices(paramD:dict) -> None:
    """ overallGenomeRelatedIndices:
            Accepts the parameter dictionary as input. Calculates the average
            amino-acid identity (AAI) and the average nucleotide identity (ANI)
            for the data-set specified in the parameter dictionary. Writes the
            results to file. Does not return.
    """
    # constants
    INDENT  = ' '*4
    PRINT_1 = 'Calculating average amino acid identity (AAI) ... '
    SAVE_PRINT = INDENT + 'Saving result to '
    PRINT_2 = 'Calculating average nucleotide identity (ANI) ... '
    DONE = 'Done.'

    # extract relevant data from paramD
    aaiFN = paramD['aaiFN']
    aniFN = paramD['aniFN']
    aaiHeatmapFN = paramD['aaiHeatmapFN']
    treeFN = paramD['speciesTreeFN']

    # calculate AAI
    print(PRINT_1, end='', flush=True)
    aaiD = _calculateAAI(paramD)
    print(DONE)

    # save the result to file
    print(SAVE_PRINT + aaiFN)
    __saveOgriMatrix(aaiD, aaiFN)


def __saveOgriMatrix(ogriD:dict, filename:str) -> None:
    """ saveOgriMatrix:
            Accepts a matrix of OGRI values in dictionary format and a string 
            indicating the where to save the matrix as inputs. Saves the matrix
            to the file as a tab-delimited matrix with both column names and
            row names. Does not return.
    """
    # constants
    DELIM_CHAR = '\t'
    EOL_CHAR = '\n'
    NUM_DECIMAL_PLACES = 3

    # get a set of all the unique species names
    speciesNamesS = set()
    for key in ogriD.keys():
        # stop looping once all the species have been added to the set
        if len(speciesNamesS) == int(math.sqrt(len(ogriD.keys()))):
            break

        # otherwise, keep adding species names to the set
        else:
            speciesNamesS.add(key[0])
            speciesNamesS.add(key[1])

    # impose an order by converting the set to a list
    speciesNamesL = list(speciesNamesS)

    # add the header to the out string
    outStr = ''
    for spe1 in speciesNamesL:
        # start with a tab, then each species name separated by a tab
        outStr += DELIM_CHAR + spe1
    
    # add a new line at the end of the headers (column names)
    outStr += EOL_CHAR
    
    # now start adding rows to the string
    for spe1 in speciesNamesL:
        # add the species name at the beginning of a row
        outStr += spe1

        # include all pairwise combinations (spe2 matches current column name)
        for spe2 in speciesNamesL:
            # round to only 3 decimal points
            ogriVal = round(ogriD[(spe1, spe2)], NUM_DECIMAL_PLACES)

            # add the value to the string
            outStr += DELIM_CHAR + str(ogriVal)
        
        # add a new line to the end of each row
        outStr += EOL_CHAR

    # write the string to specified file
    filehandle = open(filename, 'w')
    filehandle.write(outStr)
    filehandle.close()
###############################################################################


###############################################################################
def _calculateAAI(paramD:dict) -> dict:
    """ calculateAAI:
            Accepts the parameter dictionary as input. Uses xenoGI's amino acid
            identity calculator. Returns a dictionary whose keys are tuples co-
            ntaining pairs of species names and whose values are the AAI for 
            the pair.
    """
    # extract relevant data from paramD
    humanMapFN:str = paramD['fileNameMapFN']
    blastJoinStr:str = paramD['blastFileJoinStr']
    evalThresh:str = paramD['evalueThresh']
    blastFilePath:str = paramD['blastFilePath']

    # get a list of strain names from the human map file
    humanMapD  = __loadHumanMap(humanMapFN)
    taxaNamesL = list(humanMapD.values())
    
    # get the blast directory and the blast extension from the path
    blastDir, blastExt = blastFilePath.split("*")

    # calculate the AAI
    return xenoGI.analysis.aminoAcidIdentity(taxaNamesL,
                                             blastJoinStr,
                                             blastDir,
                                             blastExt,
                                             evalThresh)


def makeAaiHeatmap(paramD:dict, outgroup:str) -> None:
    """ makeAaiHeatmap:
            Accepts the parameter dictionary and an outgroup as inputs. Calls a
            custom R-script to generate a heatmap of the AAI values. Does not
            return.
    """
    # get the filename of the R script and source it
    rScript = os.path.join(PHANTASM_DIR, "PHANTASM", "aaiHeatmap.R")
    robjects.r['source'](rScript)

    # extract data from paramD
    treeFN:str = paramD['speciesTreeFN']
    aaiFN:str = paramD['aaiFN']
    pdfOutFN:str = paramD['aaiHeatmapFN']

    # create the python function that can call the R function
    heatmapRunner = robjects.globalenv['heatmapRunner']

    # make the outgroup an R vector
    outgroupVec = robjects.StrVector([outgroup])

    # call the function to generate the heatmap
    heatmapRunner(treeFN, aaiFN, outgroupVec, pdfOutFN)

###############################################################################


###############################################################################
def _calculateANI(paramD:dict) -> dict:
    """ calculateANI:
            Accepts the parameter dictionary as input. Calculates the average
            nucleotide identity for the dataset specified by the input. Returns
            a dictionary whose keys are pairs of species names and whose values
            are the corresponding ANI value.
    """
    # extract relevant data from paramD
    blastJoinStr = paramD['blastFileJoinStr']
    aniBlastPath = paramD['aniBlastFilePath']

    # get all the blast files
    blastFiles = glob.glob(aniBlastPath)

    # if the blastn files don't exist, then fix that
    if blastFiles == []:
        __aniBlastnPrep(paramD)
        __aniBlastnWrapper(paramD)

    # for each file, begin calculating the ANI
    aniD = dict()
    for blastFN in blastFiles:
        # extract the query and subject names from the blast filename
        basename:str = os.path.basename(blastFN)
        basename = os.path.splitext(basename)[0]
        query, subject = basename.split(blastJoinStr)

        # get the ANI value for the pair and save it in the dictionary
        aniD[(query, subject)] = __calculateAniFromBlastnFile(blastFN)

    return aniD


def __calculateAniFromBlastnFile(blastFN:str) -> float:
    """ calculateAniFromBlastnFile:
            Accepts a string indicating the path to a blastn table with the fo-
            rmat indicated by the constants below as input. Calculates the ave-
            rage nucleotide identity (ANI) for that blast file. Returns the ca-
            lculated value as a float.
    """
    # constants
    DELIM_CHAR = '\t'
    QSEQID = 0
    SSEQID = 1
    PIDENT = 2
    LENGTH = 3
    QLEN   = 4
    EVALUE = 5
    BITSCO = 6

    # parse the table
    parsedBlastnTable = parseCsv(blastFN, DELIM_CHAR)

    # initialize variables for looping
    row:list
    sumPercId = 0
    sumLength = 0

    # for each row in the blast table
    for row in parsedBlastnTable:
        # extract specific values
        alignLen = int(row[LENGTH])
        percId   = float(row[PIDENT])
        queryLen = int(row[QLEN])

        # recalculate the overall sequence identity
        overallSeqId = percId * alignLen / queryLen

        # determine the percentage of the alignable region
        percAligned = alignLen / queryLen * 100

        # only proceed >30% overall identity and >70% alignable region
        if overallSeqId >= 30 and percAligned >= 70:
            # weight the percent id by the alignment length
            sumPercId += percId * alignLen
            sumLength += alignLen
    
    # return the mean average nucleotide identity
    return sumPercId / sumLength


def __aniBlastnWrapper(paramD:dict) -> None:
    """ aniBlastnWrapper:
            Accepts the parameter dictionary as input. Uses the dictionary to
            determine which blast comparisons are required and then uses blastn
            to make those comparison files. Does not return.
    """
    # constants
    FRAG_END = '_fragmented.fna'

    # extract relevant data from paramD
    blastJoinStr = paramD['blastFileJoinStr']
    aniBlastPath = paramD['aniBlastFilePath']
    blastExeDir  = paramD['blastExecutDirPath']
    workdir = paramD['aniWorkDir']

    # get the file extension for the blastn results
    fileExt = os.path.splitext(aniBlastPath)[1]

    # get the path to the blastn executable
    blastnExe = os.path.join(blastExeDir, 'blastn')

    # make the blast directory
    blastdir = os.path.dirname(aniBlastPath)
    os.mkdir(blastdir)

    # get query fna files
    queryFnaFiles = glob.glob(os.path.join(workdir, '*' + FRAG_END))

    # run all pairwise blast comparisons
    for queryFN in queryFnaFiles:
        for q2 in queryFnaFiles:
            # extract the name of the blastdb for the comparison
            subjectDb = q2[:-len(FRAG_END)]

            # get the human map names for both the query and subject
            qname = os.path.basename(queryFN)[:-len(FRAG_END)]
            sname = os.path.basename(subjectDb)

            # construct the filename for the blastp result
            outFN = qname + blastJoinStr + sname + fileExt
            outFN = os.path.join(blastdir, outFN)

            # run blastn
            __runAniBlastn(queryFN, subjectDb, outFN, blastnExe)


def __runAniBlastn(queryFN:str, subjectDb:str, outFN:str, blastnExe:str) \
                                                                       -> None:
    """ runAniBlastn:
            Accepts four strings as inputs: the filename for a query fasta file
            (fragmented), the filename for a blastdb, a filename to save the 
            blastn results, and the full path the blastn executatble. Calls the
            executable to and saves the result to the specified file. Does not
            return.
    """
    # constants
    HEADERS = ['qseqid', 'sseqid', 'pident', 
               'length', 'qlen', 'evalue', 'bitscore']
    OUTFMT = '6'

    # make the outfmt string
    outfmtStr = __makeOutfmtString(OUTFMT, HEADERS)

    # construct the blastn command
    blastnCmd = NcbiblastnCommandline(cmd=blastnExe,
                                      query=queryFN,
                                      db=subjectDb,
                                      out=outFN,
                                      outfmt=outfmtStr)

    # execute the blastn command
    blastnCmd()


def __aniBlastnPrep(paramD:dict) -> None:
    """ aniBlastnPrep:
            Accepts the parameter dictionary as input. Makes nucleotide fasta
            files (.fna) and blastn databases required to run the ANI blast co-
            mparisons. Does not return.
    """
    # constants
    FNA_EXT = '.fna'
    FRAG_STR = '_fragmented'

    # extract relevant data from paramD
    gbkFilePath = paramD['genbankFilePath']
    aniWorkDir  = paramD['aniWorkDir']
    humanMapFN  = paramD['fileNameMapFN']
    blastExePath = paramD['blastExecutDirPath']

    # get genbank files
    genbankFiles = glob.glob(gbkFilePath)

    # make ANI working directory
    os.mkdir(aniWorkDir)

    # load the human map file
    humanMapD = __loadHumanMap(humanMapFN)

    # for each genbank file
    for gbkFN in genbankFiles:
        # remove the path information from the filename
        base = os.path.basename(gbkFN)
        
        # use the filename to look up the human name
        humanName = humanMapD[base]

        # add the working directory to the human name
        humanName = os.path.join(aniWorkDir, humanName)

        # make the filenames for the full and fragmented fasta files
        fullFastaFN = humanName + FNA_EXT
        fragFastaFN = humanName + FRAG_STR + FNA_EXT

        # make the regular fasta file (used to make blast db)
        __genbankToFna(gbkFN, fullFastaFN)

        # make the fragmented fasta (used as blastn queries)
        __makeFragmentedFnaForAni(gbkFN, fragFastaFN)

        # make a blastdb for the regular fasta
        __makeblastdbFromFna(fullFastaFN, humanName, blastExePath)


def __loadHumanMap(humanMapFN:str) -> dict:
    """ loadHumanMap:
            Accepts a string indicating the filename of the human map file as
            input. Constructs a dictionary keyed by gbff filenames with the co-
            rresponding human names as the values. Returns the dictionary.
    """
    # constants
    FILE_NAME_IDX  = 0
    HUMAN_NAME_IDX = 1

    # read the file into memory
    parsed = parseCsv(humanMapFN, '\t')

    # create the dict
    humanMapD = dict()
    for row in parsed:
        # keys are filenames; values are human names
        humanMapD[row[FILE_NAME_IDX]] = row[HUMAN_NAME_IDX]
    
    return humanMapD


def __genbankToFna(genbankFN:str, fastaFN:str) -> None:
    """ genbankToFna:
            Accepts two strings as input: one indicating the path to a genbank
            file and another indicating where to write the nucleotide fasta fi-
            le (.fna). Converts the genbank into fasta format. Does not return.
    """
    # constants
    IN_FORMAT = 'genbank'
    OUT_FORMAT = 'fasta'

    # parse the genbank file
    parsed:SeqIO.InsdcIO.GenBankIterator = SeqIO.parse(genbankFN, IN_FORMAT)

    # get all of the records as a list
    allRecords = list(parsed.records)

    # write the records to specified file
    SeqIO.write(allRecords, fastaFN, OUT_FORMAT)
    

def __makeFragmentedFnaForAni(genbankFN:str, fastaFN:str) -> None:
    """ makeFragmentedFnaForAni:
            Accepts a string indicating the path to a genbank file and a string
            indicating where to save the resulting fasta file as inputs. Fragm-
            ents the genome into 1020bp pieces, and writes a fasta containing
            these fragments to the specified file. Does not return.
    """
    # constants
    MAX_FRAGMENT_SIZE = 1020  # bp
    IN_FORMAT = 'genbank'
    OUT_FORMAT = 'fasta'
    SEP_1 = '|'
    SEP_2 = '..'

    # import the genbank into memory
    parsed:SeqIO.InsdcIO.GenBankIterator = SeqIO.parse(genbankFN, IN_FORMAT)

    # make the fragmented fasta file
    fragRecords = list()
    record:SeqIO.SeqRecord
    for record in parsed.records:
        # determine the number of iterations for the record
        numIter = math.ceil(len(record.seq) / MAX_FRAGMENT_SIZE)

        for iter in range(numIter):
            # determine the sequence boundaries based on the iteration
            start = iter * MAX_FRAGMENT_SIZE
            end = start + MAX_FRAGMENT_SIZE
            
            # slice the sequence to reflect the boundaries
            fragSeq = record.seq[start:end]

            # the starting base pair is always 1 larger than 'start'
            start += 1

            # then end base pair matches except at the end of a sequence
            if len(fragSeq) < MAX_FRAGMENT_SIZE:
                end = len(record.seq)
        
            # build the new record for the fragment
            fragRec = SeqIO.SeqRecord(fragSeq)
            fragRec.id = record.id + SEP_1 + str(start) + SEP_2 + str(end)
            fragRec.description = ''

            # add the new record to the list
            fragRecords.append(fragRec)

    # save the fragmented fasta to the specified file
    SeqIO.write(fragRecords, fastaFN, OUT_FORMAT)
###############################################################################

