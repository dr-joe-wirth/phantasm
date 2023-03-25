# Author: Joseph S. Wirth

import glob, logging, math, os, sys, subprocess
import rpy2.robjects as robjects
from Bio import SeqIO
from PHANTASM.Parameter import Parameters
from PHANTASM.utilities import parseCsv, loadHumanMap
from PHANTASM.taxonomy.Taxonomy import Taxonomy
from PHANTASM.downloadGbff import _makeTaxonName
from param import XENOGI_DIR, PHANTASM_DIR
sys.path.insert(0,os.path.join(sys.path[0],XENOGI_DIR))
import xenoGI.analysis


###############################################################################
def overallGenomeRelatedIndices(paramO:Parameters, outgroup:Taxonomy) -> None:
    """ overallGenomeRelatedIndices:
            Accepts a Parameters object and an outgroup Taxonomy as inputs.
            Calculates the average amino acid identity (AAI) and the average
            nucleotide identity (ANI) for the data-set specified in the Parame-
            ters object and omits the outgroup from these calculations. Writes
            the results to file. Does not return.
    """
    # constants
    GAP  = ' '*4
    PRINT_1 = 'Calculating average amino acid identity (AAI) ... '
    SAVE_PRINT = GAP + 'Saving result to '
    PRINT_2 = 'Calculating average nucleotide identity (ANI) ... '
    DONE = 'Done.'

    logger = logging.getLogger(__name__ + "." + overallGenomeRelatedIndices.__name__)

    # extract relevant data from paramO
    aaiFN = paramO.aaiFN
    aniFN = paramO.aniFN

    # calculate AAI
    print(PRINT_1, end='', flush=True)
    logger.info(PRINT_1)
    aaiD = _calculateAAI(paramO, outgroup)
    print(DONE)
    logger.info(DONE)

    # save the result to file
    print(SAVE_PRINT + aaiFN)
    __saveOgriMatrix(aaiD, aaiFN)

    # calculate ANI
    print(PRINT_2, end='', flush=True)
    logger.info(PRINT_2)
    aniD = _calculateANI(paramO, outgroup)
    print(DONE)
    logger.info(DONE + "\n")

    # save the result to file
    print(SAVE_PRINT + aniFN)
    __saveOgriMatrix(aniD, aniFN)


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
def _calculateAAI(paramO:Parameters, outgroup:Taxonomy) -> dict:
    """ calculateAAI:
            Accepts a Parameters object and an outgroup Taxonomy as inputs.
            Calculates the average amino acid identity (AAI) for each pair of
            strains in the analysis. Constructs and returns a dictionary whose
            keys are strain pairs (tuple) and whose values are the respective
            AAI values.
    """
    # extract relevant data from paramO
    wgsMapFN = paramO.fileNameMapFN

    # the outgroup will not be the root if phantasm selected reference genomes
    if not outgroup.isRoot():
        # get the outgroup name
        outgroupName = _makeTaxonName(outgroup)
    
    # the outgroup will be the root if the user specified reference genomes
    else:
        outgroupName = outgroup.ncbiName

    # get a list of all the strains
    allStrainsL = list(loadHumanMap(wgsMapFN).values())

    # remove the outgroup from the list
    allStrainsL.remove(outgroupName)

    # convert the list to a tuple
    allStrainsT = tuple(allStrainsL)

    # calculate and return the aai
    return xenoGI.analysis.calculateAAI(paramO.toDict(), allStrainsT)


def makeAaiHeatmap(paramO:Parameters, outgroup:Taxonomy) -> None:
    """ makeAaiHeatmap:
            Accepts a Parameters object and an outgroup (Taxonomy) as inputs.
            Creates a heatmap of average amino acid identity. Does not return.
    """
    # extract data from paramD
    treeFN:str = paramO.speciesTreeFN
    aaiFN:str = paramO.aaiFN
    pdfOutFN:str = paramO.aaiHeatmapFN

    __makeHeatmap(outgroup, treeFN, aaiFN, pdfOutFN)

###############################################################################


###############################################################################
def _calculateANI(paramO:Parameters, outgroup:Taxonomy) -> dict:
    """ calculateANI:
            Accepts a Parameters object as input. Uses pyani's average nucleot-
            ide calculator. Returns a dictionary whose keys are tuples contain-
            ing pairs of species names and whose values are the ANI for the sp-
            ecified pair.
    """
    # the outgroup will not be the root if phantasm selected reference genomes
    if not outgroup.isRoot():
        # determine the outgroup name
        outgroupName = _makeTaxonName(outgroup)
    
    # the outgroup will be the root if the user specified reference genomes
    else:
        outgroupName = outgroup.ncbiName

    # prepare all files needed to pyani to run
    fnaDir = __aniPrep(paramO, outgroupName)

    # calculate and return ani
    return __aniRunner(fnaDir, paramO)
    

def __aniPrep(paramO:Parameters, outgroupName:str) -> str:
    """ aniPrep:
            Accepts a Parameters object as input. Creates the ANI working dire-
            ctory and makes a nucleotide fasta for each whole genome sequence.
            Returns a string indicating the directory containing the newly cre-
            ated fasta files.
    """
    # constants
    FNA_EXT = ".fna"

    # extract necessary data from paramO
    gbFilePath = paramO.genbankFilePath
    aniDir = paramO.aniWorkDir
    wgsMapFN = paramO.fileNameMapFN

    # make the ani directory if it does not exist
    if not os.path.exists(aniDir):
        os.mkdir(aniDir)
    
    # make the fasta directory if it does not exist
    fnaDir = os.path.join(aniDir, "fasta")
    if not os.path.exists(fnaDir):
        os.mkdir(fnaDir)

    # get a list of all gbff files
    gbFilesL = glob.glob(gbFilePath)

    # look up the outgroup filename
    humanMapD = loadHumanMap(wgsMapFN)
    for filename in humanMapD.keys():
        # assign the variable and stop looping once found
        if humanMapD[filename] == outgroupName:
            outgroupFN = filename
            break

    # go through the list of gbff files and convert to fna
    for gbFN in gbFilesL:
        # only proceed if the file is not the outgroup
        if outgroupFN not in gbFN:
            # drop file path and extension
            basename = os.path.basename(gbFN)
            basename = os.path.splitext(basename)[0]

            # make fasta file name
            fnaFN = os.path.join(fnaDir, basename + FNA_EXT)

            # create the fasta file
            __genbankToFna(gbFN, fnaFN)
    
    return fnaDir


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


def __aniRunner(fnaDir:str, paramO:Parameters) -> dict:
    """ aniRunner:
            Accepts a string indicating the directory containing the nucleotide
            fasta files and a Parameters object as inputs. Uses pyani to calcu-
            late the average nucleotide identities for all genomes in the prov-
            ided directory. Returns a dictionary whose keys are pairs of speci-
            es and whose values are the ANI scores for that pair.
    """
    # constants
    PYANI_O = "pyani"
    OUTFILE = "ANIm_percentage_identity.tab"

    # extract necessary data from paramO
    aniDir = paramO.aniWorkDir
    numThreads = paramO.numProcesses
    humanMapFN = paramO.fileNameMapFN

    # determine the path for pyani's output
    outDir = os.path.join(aniDir, PYANI_O)

    # build the ANI command
    cmd = ['average_nucleotide_identity.py',
            '-i', fnaDir,
            '-o', outDir,
            '--workers', str(numThreads)]
    
    # run the command and suppress the standard error
    subprocess.run(cmd, stderr=subprocess.DEVNULL)

    # convert the file to a dictionary with pairs of human names as keys
    pyaniOutFN = os.path.join(outDir, OUTFILE)
    return __createAniD(humanMapFN, pyaniOutFN)


def __createAniD(humanMapFN:str, aniFN:str) -> dict:
    """ createAniD:
            Accepts a string indicating the path to the human map file and a
            string indicating the path to the output of pyani's analysis. Crea-
            tes and returns a dictionary whose keys are pairs of species and
            whose values are the average nucleotide identities for the pair of
            species.
    """
    # load the human map as a dictionary
    humanMapD:dict = loadHumanMap(humanMapFN)

    # keys contain file extensions; remove them from the keys
    tempD = dict()
    for oldKey in humanMapD.keys():
        newKey = os.path.splitext(oldKey)[0]
        tempD[newKey] = humanMapD[oldKey]

    humanMapD = tempD

    # parse the ANI results into a list of lists
    aniL = parseCsv(aniFN, delim='\t')

    # the first row contains the headers; all sequence names
    header = aniL[0]

    # initialize the output
    aniD = dict()

    # for each row in the file (excluding the header)
    for row in aniL[1:]:
        # for each column (first column is the sequence name)
        for idx in range(1,len(row)):
            # determine the two sequences that are being compared
            seqA = row[0]
            seqB = header[idx]

            # use humanMapD to get the human names for each sequence
            nameA = humanMapD[seqA]
            nameB = humanMapD[seqB]
            
            # save the result as a float percentage in the dictionary
            aniD[(nameA, nameB)] = float(row[idx]) * 100
    
    return aniD


def makeAniHeatmap(paramO:Parameters, outgroup:Taxonomy) -> None:
    """ makeAaiHeatmap:
            Accepts a Parameters object and an outgroup (Taxonomy) as inputs.
            Creates a heatmap of average nucleotide identity. Does not return.
    """
    # extract data from paramD
    treeFN:str = paramO.speciesTreeFN
    aniFN:str = paramO.aniFN
    pdfOutFN:str = paramO.aniHeatmapFN

    __makeHeatmap(outgroup, treeFN, aniFN, pdfOutFN)

###############################################################################
def __makeHeatmap(outgroup:Taxonomy, treeFN:str, matFN:str, pdfFN:str) -> None:
    """ makeHeatmap:
            Accepts an outgroup (Taxonomy) and three strings indicating the pa-
            ths to the species tree, a similariy matrix (ie. ANI/AAI), and the
            output pdf file as inputs. Calls a custom R-script to generate a
            heatmap of the values in the provided matrix. Does not return.
    """
    # get the filename of the R script and source it
    rScript = os.path.join(PHANTASM_DIR, "PHANTASM", "heatmap.R")
    robjects.r['source'](rScript)

    # create the python function that can call the R function
    heatmapRunner = robjects.globalenv['heatmapRunner']

    # the outgroup will not be the root if phantasm selected reference genomes
    if not outgroup.isRoot():
        # get the outgroup name
        outgroupTaxonName = _makeTaxonName(outgroup)
    
    # the outgroup will be the root if the user specified reference genomes
    else:
        outgroupTaxonName = outgroup.ncbiName

    # make the outgroup an R vector
    outgroupVec = robjects.StrVector([outgroupTaxonName])

    # call the function to generate the heatmap
    heatmapRunner(treeFN, matFN, outgroupVec, pdfFN)

