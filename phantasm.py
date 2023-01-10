# Author: Joseph S. Wirth

from PHANTASM.utilities import validEmailAddress, getParamO_1, getParamO_2, getParamO_3, checkForValidInputGenomes, checkForValidExecutables, checkForValidHumanMapFile, getLpsnAge
from PHANTASM.main import getPhyloMarker, refinePhylogeny, knownPhyloMarker, analyzeSpecifiedGenomes
from PHANTASM.findMissingNeighbors import _locusTagToGeneNum
from param import PHANTASM_DIR
import glob, os, sys

# constants
PHANTASM_PY = "python3 " + os.path.join(PHANTASM_DIR, "phantasm.py")
VERSION = "v1.0.0"
JOB_0A = 'help'
JOB_0B = "-h"
JOB_0C = ["-v", "version"]
JOB_1 = 'getPhyloMarker'
JOB_2 = 'refinePhylogeny'
JOB_3 = 'knownPhyloMarker'
JOB_4 = 'analyzeGenomes'
SEP = ","
GAP = " "*4
LOCUS_TAG = "--locus_tag"
GENE_NUM  = "--gene_num"

# error messages
ERR_MSG_1 = "Incorrect syntax used.\nType '" + PHANTASM_PY + " help' for more information."
ERR_MSG_2 = "Invalid email address"
ERR_MSG_3 = "The number of genes is not a multiple of the number of input genomes."
ERR_MSG_4 = "Invalid flag: "
ERR_MSG_5 = "One or more specified genes could not be extracted."
ERR_MSG_6 = "The specified genome directory is not a directory."
ERR_MSG_7 = "Less than 2 files were found in the specified genome directory."
ERR_MSG_8 = "The specified output directory already exists."
ERR_MSG_10A = "Invalid task: "
ERR_MSG_10B = "\n\ntype '" + PHANTASM_PY + " help' for information."

# reference message
REF_MSG = "\nIf you use this software in your research, please cite our paper:\n" + \
          GAP + "Automating microbial taxonomy workflows with PHANTASM: PHylogenomic\n" + \
          GAP + "ANalyses for the TAxonomy and Systematics of Microbes\n" + \
          GAP*2 + "Joseph S. Wirth and Eliot C. Bush, 2023\n" + \
          GAP*2 + "https://www.biorxiv.org/content/10.1101/2022.10.18.512716v1\n"

# version message
VERSION_MSG = "PHANTASM " + VERSION + "\n"

# no arguments message
NO_ARGS_MSG = "\ntype '" + PHANTASM_PY + " -h' for help\n"

# LPSN age message
AGE_MSG = "PHANTASM relies on data manually acquired from the LPSN.\n" + \
          GAP + "These data were last retrieved on "

# detailed help message (for JOB_0A)
DETAILED_HELP_MSG = "\nPHANTASM: PHylogenomic ANalyses for the TAxonomy and Systematics of Microbes\n\n" + \
                    "In order to run PHANTASM without a known phylogenetic marker, follow these steps:\n" + \
                    GAP + "1. cd into the desired working directory\n" + \
                    GAP + "2. identify a suitable phylogenetic marker:\n\n" + \
                    GAP*2 + "'" + PHANTASM_PY + " " + JOB_1 + " <gbff file> <email address>'\n\n" + \
                    GAP + "3. examine the file 'initialAnalysis/putativePhylogeneticMarkers.txt'\n" + \
                    GAP + "4. refine the phylogeny and perform phylogenomic analyses:\n\n" + \
                    GAP*2 + "'" + PHANTASM_PY + " " + JOB_2 + " --locus_tag <locus_tag> <gbff file> <email address>'\n" + \
                    GAP*3 + "OR\n" + \
                    GAP*2 + "'" + PHANTASM_PY + " " + JOB_2 + " --gene_num <gene_number> <gbff file> <email address>'\n\n" + \
                    GAP + "5. results can be found in the following files:\n\n" + \
                    GAP*2 + "'finalAnalysis/aai_matrix.txt'\n" + \
                    GAP*2 + "'finalAnalysis/aai_heatmap.pdf'\n" + \
                    GAP*2 + "'finalAnalysis/ani_matrix.txt'\n" + \
                    GAP*2 + "'finalAnalysis/ani_heatmap.pdf'\n" + \
                    GAP*2 + "'finalAnalysis/coreGenesSummary.txt\n" + \
                    GAP*2 + "'finalAnalysis/speciesTree.nwk'\n" + \
                    GAP*2 + "'finalAnalysis/speciesTree_outgroupPruned.nwk'\n\n" + \
                    "In order to run PHANTASM with a known phylogenetic marker, follow these steps:\n" + \
                    GAP + "1. cd into the desired working directory\n" + \
                    GAP + "2. call PHANTASM on the known marker(s) and input genome(s):\n\n" + \
                    GAP*2 + "'" + PHANTASM_PY + " " + JOB_3 + " <locus tag> <gbff file> <email address>'\n\n" + \
                    GAP + "3. results can be found in the following files:\n\n" + \
                    GAP*2 + "'finalAnalysis/aai_matrix.txt'\n" + \
                    GAP*2 + "'finalAnalysis/aai_heatmap.pdf'\n" + \
                    GAP*2 + "'finalAnalysis/ani_matrix.txt'\n" + \
                    GAP*2 + "'finalAnalysis/ani_heatmap.pdf'\n" + \
                    GAP*2 + "'finalAnalysis/coreGenesSummary.txt\n" + \
                    GAP*2 + "'finalAnalysis/speciesTree.nwk'\n" + \
                    GAP*2 + "'finalAnalysis/speciesTree_outgroupPruned.nwk'\n\n" + \
                    "In order to run PHANTASM to analyze user-specified genomes, follow these steps:\n" + \
                    GAP + "1. cd into the desired working directory\n" + \
                    GAP + "2. make a directory containing ONLY the genomes you wish to analyze\n" + \
                    GAP + "3. make a 'human map' file with the outgroup as the last entry in file (see readme for more details).\n" + \
                    GAP + "4. call PHANTASM on the specified genomes:\n\n" + \
                    GAP*2 + "'" + PHANTASM_PY + " " + JOB_4 + " <gbff directory> <human map file> <output directory> <email address>'\n\n" + \
                    GAP + "5. results can be found in the following files:\n\n" + \
                    GAP*2 + "'<output dir>/aai_matrix.txt'\n" + \
                    GAP*2 + "'<output dir>/aai_heatmap.pdf'\n" + \
                    GAP*2 + "'<output dir>/ani_matrix.txt'\n" + \
                    GAP*2 + "'<output dir>/ani_heatmap.pdf'\n" + \
                    GAP*2 + "'<output dir>/coreGenesSummary.txt\n" + \
                    GAP*2 + "'<output dir>/speciesTree.nwk'\n" + \
                    GAP*2 + "'<output dir>/speciesTree_outgroupPruned.nwk'\n\n" + \
                    "Optional Features\n" + \
                    GAP + "multiple input genomes:\n" + \
                    GAP*2 + "replace '<gbff file>' with '<gbff directory>'\n" + \
                    GAP*2 + "the specified directory must contain all of the input genome files\n\n" + \
                    GAP + "multiple phylogenetic markers:\n" + \
                    GAP*2 + "list each desired marker separated by a comma (no spaces)\n" + \
                    GAP*2 + "the number of specified markers must be a multiple of the number input genomes\n\n" + \
                    GAP + "excluding specific taxa from the final analysis:\n" + \
                    GAP*2 + "create a file in the working directory named 'excludedTaxids.txt'\n" + \
                    GAP*2 + "the file should contain exactly one NCBI Taxonomy ID to be excluded per line\n" + \
                    GAP*2 + "WARNING: this feature is experimental and should be used with caution.\n\n" + \
                    GAP + "bypassing the requirement for 16S rRNA gene sequences in the '" + JOB_1 + "' step:\n" + \
                    GAP*2 + "create a file in the working directory named 'taxids.txt'\n" + \
                    GAP*2 + "the file should contain exactly one NCBI Taxonomy ID per line\n" + \
                    GAP*3 + " * the taxids should be somewhat related to the input genome(s)\n" + \
                    GAP*3 + " * the taxids should not exceed the taxonomic rank of Family\n"

# short help message (for JOB_0B)
SHORT_HELP_MSG = "Getting detailed help\n" + \
                  GAP + PHANTASM_PY + " help\n\n\n" + \
                  "Getting this message\n" + \
                  GAP + PHANTASM_PY + " -h\n\n\n" + \
                  "Getting the version\n" + \
                  GAP + PHANTASM_PY + " -v\n" + \
                  GAP*2 + "OR\n" +\
                  GAP + PHANTASM_PY + " version\n\n\n" + \
                  "Option 1: unknown reference genomes and unknown phylogenetic markers\n" + \
                  GAP + "Step 1: rank phylogenetic markers\n" + \
                  GAP*2 + PHANTASM_PY + " " + JOB_1 + " <input genome(s)> <email>\n\n" + \
                  GAP + "Step 2: refine phylogeny using a phylogenetic marker\n" + \
                  GAP*2 + PHANTASM_PY + " " + JOB_2 + " --locus_tag <locus tag(s)> <input genome(s)> <email>\n" + \
                  GAP*3 + "OR\n" + \
                  GAP*2 + PHANTASM_PY + " " + JOB_2 + " --gene_num <gene number(s)> <input genome(s)> <email>\n\n\n" + \
                  "Option 2: unkonwn reference genomes and known phylogenetic marker\n" + \
                  GAP + PHANTASM_PY + " " + JOB_3 + " <locus tag(s)> <input genome(s)> <email>\n\n\n" + \
                  "Option 3: known reference genomes\n" + \
                  GAP + PHANTASM_PY + " " + JOB_4 + " <genome directory> <human map file> <output directory> <email>\n"


# begin main function
if __name__ == "__main__":
    # print the reference message first
    print(REF_MSG)

    # print help message if no arguments are provided
    if len(sys.argv) == 1:
        print(VERSION_MSG + NO_ARGS_MSG)
    
    else:
        # extract the job name
        job = sys.argv[1]

        # print the detailed help message if requested
        if job == JOB_0A:
            print(VERSION_MSG)
            print(DETAILED_HELP_MSG)

        # print the short help message if requested
        elif job == JOB_0B:
            print(VERSION_MSG)
            print(SHORT_HELP_MSG)
        
        elif job in JOB_0C:
            print(VERSION_MSG)

        # if JOB_1 specified
        elif job == JOB_1:
            # check that the expected number of arguments has been provided
            if len(sys.argv) != 4:
                raise SyntaxError(ERR_MSG_1)
            
            # extract the data from the arguments
            inputFile = sys.argv[2]
            email = sys.argv[3]

            # ensure that a valid email address is being used
            if not validEmailAddress(email):
                raise ValueError(ERR_MSG_2)

            # if a directory was provided, then get the list of files within it
            if os.path.isdir(inputFile):
                gbffL = glob.glob(os.path.join(inputFile, "*"))
            
            # if a single file was provided, then convert it to a one-item list
            else:
                gbffL = [inputFile]
            
            # construct the Parameters object
            paramO = getParamO_1(email)

            # check inputs
            checkForValidInputGenomes(gbffL)
            checkForValidExecutables(paramO)

            # print the LPSN age data
            age = getLpsnAge()
            print(AGE_MSG + age + "\n\n")

            # execute job 1
            getPhyloMarker(gbffL, paramO)
        
        # if JOB_2 requested
        elif job == JOB_2:
            # check that the expected number of arguments has been provided
            if len(sys.argv) != 6:
                raise SyntaxError(ERR_MSG_1)
            
            # extract the data from the arguments
            flag = sys.argv[2]
            gene = sys.argv[3]
            inputFile = sys.argv[4]
            email = sys.argv[5]

            # ensure that a valid email address is being used
            if not validEmailAddress(email):
                raise ValueError(ERR_MSG_2)
            
            # if a directory was provided
            if os.path.isdir(inputFile):
                # get the list of files within it
                gbffL = glob.glob(os.path.join(inputFile, "*"))
            
            # if a single file was provided, then convert it to a one-item list
            else:
                gbffL = [inputFile]

            # split the gene argument into a list of genes
            genesL = gene.split(SEP)

            # raise error if num genes is not a multiple of num genomes
            if len(genesL) % len(gbffL) != 0:
                raise ValueError(ERR_MSG_3)

            # construct the Parameters objects
            paramO_1 = getParamO_1(email)
            paramO_2 = getParamO_2(email)

            # check inputs
            checkForValidInputGenomes(gbffL)
            checkForValidExecutables(paramO_1)
            checkForValidExecutables(paramO_2)

            # initialize geneNumsL
            geneNumsL = list()

            # if a locus tag was provided, convert it to an integer
            if flag == LOCUS_TAG:
                for gene in genesL:
                    geneNum = _locusTagToGeneNum(gene, paramO_1.geneInfoFN)
                    geneNumsL.append(geneNum)
                    
            # otherwise, a gene number should have been specified
            elif flag == GENE_NUM:
                for gene in genesL:
                    geneNum = int(gene)
                    geneNumsL.append(geneNum)
            
            # if an invalid flag was specified, then raise an error
            else:
                raise ValueError(ERR_MSG_4 + flag)
            
            # make sure the number of gene numbers matches the number of genes
            if len(genesL) != len(geneNumsL):
                raise ValueError(ERR_MSG_5)

            # print the LPSN age data
            age = getLpsnAge()
            print(AGE_MSG + age + "\n\n")

            # execute JOB_2
            refinePhylogeny(geneNumsL, gbffL, paramO_1, paramO_2)
        
        # if JOB_3 requested
        elif job == JOB_3:
            # check that the expected number of arguments has been provided
            if len(sys.argv) != 5:
                raise SyntaxError(ERR_MSG_1)
            
            # extract the data from the arguments
            locusTag = sys.argv[2]
            inputFile = sys.argv[3]
            email = sys.argv[4]

            # ensure that a valid email address is being used
            if not validEmailAddress(email):
                raise ValueError(ERR_MSG_2)
            
            # if a directory was provided
            if os.path.isdir(inputFile):
                # get the list of files within it
                gbffL = glob.glob(os.path.join(inputFile, "*"))
            
            # if a single file was provided, then convert it to a one-item list
            else:
                gbffL = [inputFile]
            
            # split the locusTag argument into a list of locus tags
            locusTagsL = locusTag.split(SEP)

            # raise error if num locus tags is not a multiple of num genomes
            if len(locusTagsL) % len(gbffL) != 0:
                raise ValueError(ERR_MSG_3)
            
            # create the Parameters object
            paramO = getParamO_2(email)

            # check inputs
            checkForValidInputGenomes(gbffL)
            checkForValidExecutables(paramO)

            # print the LPSN age data
            age = getLpsnAge()
            print(AGE_MSG + age + "\n\n")
            
            # execute JOB_3
            knownPhyloMarker(gbffL, locusTagsL, paramO)

        # if JOB_4 requested
        elif job == JOB_4:
            # check that the expected number of arguments has been provided
            if len(sys.argv) != 6:
                raise SyntaxError(ERR_MSG_1)
            
            # parse command line arguments
            genomeDir = os.path.abspath(sys.argv[2])
            humanMapFN = os.path.abspath(sys.argv[3])
            outdir = sys.argv[4]
            email = sys.argv[5]

            # ensure that a valid email address is being used
            if not validEmailAddress(email):
                raise ValueError(ERR_MSG_2)

            # make sure genomeDir is a directory
            if not os.path.isdir(genomeDir):
                raise ValueError(ERR_MSG_6)

            # make a Parameters object
            paramO = getParamO_3(email, outdir)

            # update the parameters object with the genomes and map locations
            paramO.genbankFilePath = os.path.join(genomeDir, "*")
            paramO.fileNameMapFN = humanMapFN

            # get a list of all the genbank files
            gbffL = glob.glob(paramO.genbankFilePath)

            # make sure there are at least two files in the list
            if len(gbffL) < 2:
                raise ValueError(ERR_MSG_7)

            # check for invalid inputs
            checkForValidInputGenomes(gbffL)
            checkForValidExecutables(paramO)
            checkForValidHumanMapFile(paramO)

            # make the output directory; error if it already exists
            if not os.path.exists(paramO.workdir):
                os.mkdir(paramO.workdir)
        
            else:
                raise FileExistsError(ERR_MSG_8)

            # analyze the specified genomes
            analyzeSpecifiedGenomes(gbffL, paramO)

        # raise an error if an invalid job was specified
        else:
            raise ValueError(ERR_MSG_10A + job + ERR_MSG_10B)

