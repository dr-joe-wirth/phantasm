# Author: Joseph S. Wirth


from PHANTASM.main import getPhyloMarker, refinePhylogeny, knownPhyloMarker, analyzeSpecifiedGenomes
from PHANTASM.utilities import parseArgs, getLpsnAge, getCmdWithRedactedEmail
from PHANTASM.findMissingNeighbors import _locusTagToGeneNum
from PHANTASM.Parameter import Parameters
from param import PHANTASM_DIR
import copy, glob, logging, os, sys

# constants
PHANTASM_PY = "python3 " + os.path.join(PHANTASM_DIR, "phantasm.py")
VERSION = "v1.0.4"
JOB_0A = 'help'
JOB_0B = "-h"
JOB_0C = ("-v", "version")
JOB_1 = 'getPhyloMarker'
JOB_2 = 'refinePhylogeny'
JOB_3 = 'knownPhyloMarker'
JOB_4 = 'analyzeGenomes'
SEP = ","
GAP = " "*4

# error messages
ERR_MSG_1 = "One or more specified genes could not be extracted."
ERR_MSG_2 = "The specified output directory already exists."
ERR_MSG_3A = "Invalid task: "
ERR_MSG_3B = "\n\ntype '" + PHANTASM_PY + " help' for information."

# reference message
REF_MSG = "\nIf you use this software in your research, please cite our paper:\n" + \
          GAP + "Automating microbial taxonomy workflows with PHANTASM: PHylogenomic\n" + \
          GAP + "ANalyses for the TAxonomy and Systematics of Microbes\n" + \
          GAP*2 + "Joseph S. Wirth and Eliot C. Bush, 2023\n" + \
          GAP*2 + "https://doi.org/10.1093/nar/gkad196\n"

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
            # parse the parameters
            gbffL, locusTagsL, paramO = parseArgs()

            # print the LPSN age data
            age = getLpsnAge()
            print(AGE_MSG + age + "\n\n")

            # initialize logger
            logging.basicConfig(filename=paramO.logFN, level=logging.INFO)
            logger = logging.getLogger(__name__)
            
            # save details about the run
            logger.info(getCmdWithRedactedEmail())
            logger.info(VERSION)
            logger.info('num cpus:        ' + str(paramO.numProcesses))
            logger.info('max leaves:      ' + str(paramO.maxNumTreeLeaves))
            logger.info('reduce num core: ' + str(paramO.reduceNumCoreGenes))
            logger.info('bootstrap tree:  ' + str(paramO.numBootstraps > 0))
            logger.info('num bootstraps:  ' + str(paramO.numBootstraps) + "\n")

            # execute job 1
            logger.info("start " + JOB_1 + "\n")
            getPhyloMarker(gbffL, paramO)
            logger.info("end " + JOB_1 + "\n")            
        
        # if JOB_2 requested
        elif job == JOB_2:
            # parse the parameters
            gbffL, locusTagsL, paramO_2 = parseArgs()
            
            # create paramO_1 from paramO_2
            tmp = paramO_2.toDict()
            tmp['workdir'] = os.path.join(os.path.dirname(paramO_2.workdir), "initialAnalysis")
            paramO_1 = Parameters.fromDict(tmp)
            
            # initialize geneNumsL
            geneNumsL = list()
            
            # convert the locus tags to gene numbers
            for tag in locusTagsL:
                geneNum = _locusTagToGeneNum(tag, paramO_1.geneInfoFN)
                geneNumsL.append(geneNum)
            
            # make sure the number of gene numbers matches the number of genes
            if len(locusTagsL) != len(geneNumsL):
                raise ValueError(ERR_MSG_1)

            # print the LPSN age data
            age = getLpsnAge()
            print(AGE_MSG + age + "\n\n")
            
            # initialize logger
            logging.basicConfig(filename=paramO_1.logFN, level=logging.INFO)
            logger = logging.getLogger(__name__)
            
            # save details about the run
            logger.info(getCmdWithRedactedEmail())
            logger.info(VERSION)
            logger.info('num cpus:        ' + str(paramO_1.numProcesses))
            logger.info('max leaves:      ' + str(paramO_1.maxNumTreeLeaves))
            logger.info('reduce num core: ' + str(paramO_1.reduceNumCoreGenes))
            logger.info('bootstrap tree:  ' + str(paramO_1.numBootstraps > 0))
            logger.info('num bootstraps:  ' + str(paramO_1.numBootstraps) + "\n")

            # execute JOB_2
            logger.info("start " + JOB_2 + "\n")
            refinePhylogeny(geneNumsL, gbffL, paramO_1, paramO_2)
            logger.info('end ' + JOB_2 + "\n")
        
        # if JOB_3 requested
        elif job == JOB_3:
            # parse the parameters
            gbffL, locusTagsL, paramO = parseArgs()
            
            # print the LPSN age data
            age = getLpsnAge()
            print(AGE_MSG + age + "\n\n")
            
            # initialize logger
            logging.basicConfig(filename=paramO.logFN, level=logging.INFO)
            logger = logging.getLogger(__name__)
            
            # save details about the run
            logger.info(getCmdWithRedactedEmail())
            logger.info(VERSION)
            logger.info('num cpus:        ' + str(paramO.numProcesses))
            logger.info('max leaves:      ' + str(paramO.maxNumTreeLeaves))
            logger.info('reduce num core: ' + str(paramO.reduceNumCoreGenes))
            logger.info('bootstrap tree:  ' + str(paramO.numBootstraps > 0))
            logger.info('num bootstraps:  ' + str(paramO.numBootstraps) + "\n")
            
            # execute JOB_3
            logger.info("start " + JOB_3 + "\n")
            knownPhyloMarker(gbffL, locusTagsL, paramO)
            logger.info('end ' + JOB_3 + "\n")

        # if JOB_4 requested
        elif job == JOB_4:
            # parse the parameters
            gbffL, locusTagsL, paramO = parseArgs()

            # make the output directory; error if it already exists
            if not os.path.exists(paramO.workdir):
                os.mkdir(paramO.workdir)
        
            else:
                raise FileExistsError(ERR_MSG_2)

            # initialize logger
            logging.basicConfig(filename=paramO.logFN, level=logging.INFO)
            logger = logging.getLogger(__name__)
            
            # save details about the run
            logger.info(getCmdWithRedactedEmail())
            logger.info(VERSION)
            logger.info('num cpus:        ' + str(paramO.numProcesses))
            logger.info('max leaves:      ' + str(paramO.maxNumTreeLeaves))
            logger.info('reduce num core: ' + str(paramO.reduceNumCoreGenes))
            logger.info('bootstrap tree:  ' + str(paramO.numBootstraps > 0))
            logger.info('num bootstraps:  ' + str(paramO.numBootstraps) + "\n")
            
            # analyze the specified genomes
            logger.info("start " + JOB_4 + "\n")
            analyzeSpecifiedGenomes(gbffL, paramO)
            logger.info("end " + JOB_4 + "\n")

        # raise an error if an invalid job was specified
        else:
            raise ValueError(ERR_MSG_3A + job + ERR_MSG_3B)

