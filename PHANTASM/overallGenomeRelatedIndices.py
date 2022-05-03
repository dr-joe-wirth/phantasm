# Author: Joseph S. Wirth

import glob, math, os, sys, subprocess
import rpy2.robjects as robjects
from Bio import SeqIO
from param import XENOGI_DIR, PHANTASM_DIR
sys.path.insert(0,os.path.join(sys.path[0],XENOGI_DIR))
import xenoGI.xenoGI, xenoGI.analysis
from PHANTASM.utilities import parseCsv
from PHANTASM.taxonomy.Taxonomy import Taxonomy
from PHANTASM.downloadGbff import _makeTaxonName


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

    # calculate AAI
    print(PRINT_1, end='', flush=True)
    aaiD = _calculateAAI(paramD)
    print(DONE)

    # save the result to file
    print(SAVE_PRINT + aaiFN)
    __saveOgriMatrix(aaiD, aaiFN)

    # calculate ANI
    print(PRINT_2, end='', flush=True)
    aniD = _calculateANI(paramD)
    print(DONE)

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


def makeAaiHeatmap(paramD:dict, outgroup:Taxonomy) -> None:
    """ makeAaiHeatmap:
            Accepts the parameter dictionary and an outgroup as inputs. Calls a
            custom R-script to generate a heatmap of the AAI values. Does not
            return.
    """
    # extract data from paramD
    treeFN:str = paramD['speciesTreeFN']
    aaiFN:str = paramD['aaiFN']
    pdfOutFN:str = paramD['aaiHeatmapFN']

    __makeHeatmap(outgroup, treeFN, aaiFN, pdfOutFN)

###############################################################################


###############################################################################

def _calculateANI(paramD:dict) -> dict:
    """ calculateANI:
            Accepts the parameter dictionary as input. Uses pyani's average nu-
            cleotide calculator. Returns a dictionary whose keys are tuples co-
            ntaining pairs of species names and whose values are the ANI for 
            the pair.
    """
    # prepare all files needed to pyani to run
    fnaDir = __aniPrep(paramD)

    # calculate and return ani
    return __aniRunner(fnaDir, paramD)
    

def __aniPrep(paramD:dict) -> str:
    """ aniPrep:
            Accepts the parameter dictionary as input. Creates the ANI working
            directory and makes a nucleotide fasta for each whole genome seque-
            nce. Returns a string indicating the directory containing the newly
            created fasta files.
    """
    # constants
    FNA_EXT = ".fna"

    # extract necessary data from paramD
    gbFilePath = paramD['genbankFilePath']
    aniDir = paramD['aniWorkDir']

    # make the ani directory if it does not exist
    if not os.path.exists(aniDir):
        os.mkdir(aniDir)
    
    # make the fasta directory if it does not exist
    fnaDir = os.path.join(aniDir, "fasta")
    if not os.path.exists(fnaDir):
        os.mkdir(fnaDir)

    # get a list of all gbff files
    gbFilesL = glob.glob(gbFilePath)

    # go through the list of gbff files and convert to fna
    for gbFN in gbFilesL:
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


def __aniRunner(fnaDir:str, paramD:dict) -> dict:
    """ aniRunner:
            Accepts a string indicating the directory containing the nucleotide
            fasta files and the parameter dictionary as inputs. Uses pyani to
            calculate the average nucleotide identities for all genomes in the
            provided directory. Returns a dictionary whose keys are pairs of
            species and whose values are the ANI scores for that pair.
    """
    # constants
    PYANI_O = "pyani"
    OUTFILE = "ANIb_percentage_identity.tab"

    # extract necessary data from paramD
    aniDir = paramD['aniWorkDir']
    numThreads = paramD['numProcesses']
    humanMapFN = paramD['fileNameMapFN']

    # determine the path for pyani's output
    outDir = os.path.join(aniDir, PYANI_O)

    # call ANI function via bash
    cmd = ['average_nucleotide_identity.py',
            '-m', 'ANIb',
            '-i', fnaDir,
            '-o', outDir,
            '--workers', str(numThreads)]
    subprocess.run(cmd)

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
    humanMapD:dict = __loadHumanMap(humanMapFN)

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
            
            # save the result as a float in the dictionary
            aniD[(nameA, nameB)] = float(row[idx])
    
    return aniD


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


def makeAniHeatmap(paramD:dict, outgroup:Taxonomy) -> None:
    """ makeAaiHeatmap:
            Accepts the parameter dictionary and an outgroup as inputs. Calls a
            custom R-script to generate a heatmap of the ANI values. Does not
            return.
    """
    # extract data from paramD
    treeFN:str = paramD['speciesTreeFN']
    aniFN:str = paramD['aniFN']
    pdfOutFN:str = paramD['aniHeatmapFN']

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

    # get the taxon name for the outgroup
    outgroupTaxonName = _makeTaxonName(outgroup)

    # make the outgroup an R vector
    outgroupVec = robjects.StrVector([outgroupTaxonName])

    # call the function to generate the heatmap
    heatmapRunner(treeFN, matFN, outgroupVec, pdfFN)

