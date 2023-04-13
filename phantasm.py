# Author: Joseph S. Wirth


from PHANTASM.main import getPhyloMarker, refinePhylogeny, knownPhyloMarker, analyzeSpecifiedGenomes
from PHANTASM.utilities import parseArgs, getLpsnAge, redactEmailAddress
from PHANTASM.findMissingNeighbors import _locusTagToGeneNum
from PHANTASM.Parameter import Parameters
from param import PHANTASM_DIR
import copy, glob, logging, os, sys

# constants
PHANTASM_PY = "python3 " + os.path.join(PHANTASM_DIR, "phantasm.py")
VERSION = "v1.0.4"
JOB_0A = '--help'
JOB_0B = "-h"
JOB_0C = ("-v", "--version")
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
DETAILED_HELP_MSG = "# Getting detailed help (this message):\n\n" + \
                    GAP + PHANTASM_PY + " " + JOB_0A + "\n\n\n" + \
                    "# Getting a short help message:\n\n" + \
                    GAP + PHANTASM_PY + " " + JOB_0B + "\n\n\n" + \
                    "# Getting the version:\n\n" + \
                    GAP + PHANTASM_PY + " " + JOB_0C[0] + "\n" + \
                    GAP*2 + " OR \n" + \
                    GAP + PHANTASM_PY + " " + JOB_0C[1] + "\n\n\n" + \
                    "# Identifying phylogenetic markers (option 1, step 1)\n\n" + \
                    GAP + PHANTASM_PY + " " + JOB_1 + " [-ieNLF]\n\n" + \
                    GAP + "required arguments:\n" + \
                    GAP*2 + "-i, --input <file>         gbff file or a directory containing gbff files\n" + \
                    GAP*2 + "-e, --email <email>        email address\n\n" + \
                    GAP + "optional arguments:\n" + \
                    GAP*2 + "-N, --num_threads <int>    number of processors to use                    [default: 1]\n" + \
                    GAP*2 + "-L, --max_leaves <int>     maximum number of leaves in the species tree   [default: 50]\n" + \
                    GAP*2 + "-F, --fewer_coregenes      limit the core genes to those with ≤5% gaps    [default: no limiting]\n\n" + \
                    GAP + "results:\n" + \
                    GAP*2 + "'initialAnalysis/putativePhylogeneticMarkers.txt'\n\n\n" + \
                    "# Refining phylogeny and performing phylogenomic analyses (option 1, step 2)\n\n" + \
                    GAP + PHANTASM_PY + " " + JOB_2 + " [-ietNLBF]\n\n" + \
                    GAP + "required arguments:\n" + \
                    GAP*2 + "-i, --input <file>         gbff file or a directory containing gbff files\n" + \
                    GAP*2 + "-e, --email <email>        email address\n" + \
                    GAP*2 + "-t, --locus_tags <str>     a comma-separated list of locus tags to use a phylogenetic markers\n\n" + \
                    GAP + "optional arguments:\n" + \
                    GAP*2 + "-N, --num_threads <int>    number of processors to use                    [default: 1]\n" + \
                    GAP*2 + "-L, --max_leaves <int>     maximum number of leaves in the species tree   [default: 50]\n" + \
                    GAP*2 + "-B, --bootstrap <int>      number of bootstraps to perform                [default: no bootstrapping]\n" + \
                    GAP*2 + "-F, --fewer_coregenes      limit the core genes to those with ≤5% gaps    [default: no limiting]\n\n" + \
                    GAP + "results:\n" + \
                    GAP*2 + "'finalAnalysis/aai_matrix.txt'\n" + \
                    GAP*2 + "'finalAnalysis/aai_heatmap.pdf'\n" + \
                    GAP*2 + "'finalAnalysis/ani_matrix.txt'\n" + \
                    GAP*2 + "'finalAnalysis/ani_heatmap.pdf'\n" + \
                    GAP*2 + "'finalAnalysis/coreGenesSummary.txt\n" + \
                    GAP*2 + "'finalAnalysis/speciesTree.nwk'\n" + \
                    GAP*2 + "'finalAnalysis/speciesTree_outgroupPruned.nwk'\n\n\n" + \
                    "# Identifying genomes using a known phylogenetic marker and performing phylogenomic analyses (option 2)\n\n" + \
                    GAP + PHANTASM_PY + " " + JOB_3 + " [-ietONLBF]\n\n" + \
                    GAP + "required arguments:\n" + \
                    GAP*2 + "-i, --input <file>         gbff file or a directory containing gbff files\n" + \
                    GAP*2 + "-e, --email <email>        email address\n" + \
                    GAP*2 + "-t, --locus_tags <str>     a comma-separated list of locus tags to use a phylogenetic markers\n\n" + \
                    GAP + "optional arguments:\n" + \
                    GAP*2 + "-O, --out <dir>            output directory                               [default: '" + os.path.join(".", "finalAnalysis") + "']\n" + \
                    GAP*2 + "-N, --num_threads <int>    number of processors to use                    [default: 1]\n" + \
                    GAP*2 + "-L, --max_leaves <int>     maximum number of leaves in the species tree   [default: 50]\n" + \
                    GAP*2 + "-B, --bootstrap <int>      number of bootstraps to perform                [default: no bootstrapping]\n" + \
                    GAP*2 + "-F, --fewer_coregenes      limit the core genes to those with ≤5% gaps    [default: no limiting]\n\n" + \
                    GAP + "results (directory may vary if '-O' or '--out' used):\n" + \
                    GAP*2 + "'finalAnalysis/aai_matrix.txt'\n" + \
                    GAP*2 + "'finalAnalysis/aai_heatmap.pdf'\n" + \
                    GAP*2 + "'finalAnalysis/ani_matrix.txt'\n" + \
                    GAP*2 + "'finalAnalysis/ani_heatmap.pdf'\n" + \
                    GAP*2 + "'finalAnalysis/coreGenesSummary.txt\n" + \
                    GAP*2 + "'finalAnalysis/speciesTree.nwk'\n" + \
                    GAP*2 + "'finalAnalysis/speciesTree_outgroupPruned.nwk'\n\n\n" + \
                    "# Performing phylogenomic analyses on a user-specified set of genomes (option 3)\n\n" + \
                    GAP + PHANTASM_PY + " " + JOB_4 + " [-iemONLBF]\n\n" + \
                    GAP + "required arguments:\n" + \
                    GAP*2 + "-i, --input <file>         gbff file or a directory containing gbff files\n" + \
                    GAP*2 + "-e, --email <email>        email address\n" + \
                    GAP*2 + "-m, --map_file <file>      a file with two tab-separated columns (no headers): filename, taxon name\n\n" + \
                    GAP + "optional arguments:\n" + \
                    GAP*2 + "-O, --out <dir>            output directory                               [default: '" + os.path.join(".", "finalAnalysis") + "']\n" + \
                    GAP*2 + "-N, --num_threads <int>    number of processors to use                    [default: 1]\n" + \
                    GAP*2 + "-L, --max_leaves <int>     maximum number of leaves in the species tree   [default: 50]\n" + \
                    GAP*2 + "-B, --bootstrap <int>      number of bootstraps to perform                [default: no bootstrapping]\n" + \
                    GAP*2 + "-F, --fewer_coregenes      limit the core genes to those with ≤5% gaps    [default: no limiting]\n\n" + \
                    GAP + "results (directory may vary if '-O' or '--out' used):\n" + \
                    GAP*2 + "'finalAnalysis/aai_matrix.txt'\n" + \
                    GAP*2 + "'finalAnalysis/aai_heatmap.pdf'\n" + \
                    GAP*2 + "'finalAnalysis/ani_matrix.txt'\n" + \
                    GAP*2 + "'finalAnalysis/ani_heatmap.pdf'\n" + \
                    GAP*2 + "'finalAnalysis/coreGenesSummary.txt\n" + \
                    GAP*2 + "'finalAnalysis/speciesTree.nwk'\n" + \
                    GAP*2 + "'finalAnalysis/speciesTree_outgroupPruned.nwk'\n\n\n" + \
                    "# Optional Features\n" + \
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
                  GAP + PHANTASM_PY + " " + JOB_0A + "\n\n\n" + \
                  "Getting this message\n" + \
                  GAP + PHANTASM_PY + " " + JOB_0B + "\n\n\n" + \
                  "Getting the version\n" + \
                  GAP + PHANTASM_PY + " -v\n\n\n" + \
                  "Option 1: unknown reference genomes and unknown phylogenetic markers\n" + \
                  GAP + "Step 1: rank phylogenetic markers\n" + \
                  GAP*2 + PHANTASM_PY + " " + JOB_1 + " -i <input genome(s)> -e <email>\n\n" + \
                  GAP + "Step 2: refine phylogeny using a phylogenetic marker\n" + \
                  GAP*2 + PHANTASM_PY + " " + JOB_2 + " -t <locus tag(s)> -i <input genome(s)> -e <email>\n\n\n" + \
                  "Option 2: unkonwn reference genomes and known phylogenetic marker\n" + \
                  GAP + PHANTASM_PY + " " + JOB_3 + " -t <locus tag(s)> -i <input genome(s)> -e <email>\n\n\n" + \
                  "Option 3: known reference genomes\n" + \
                  GAP + PHANTASM_PY + " " + JOB_4 + " -i <genome directory> -m <map file> -o <output directory> -e <email>\n"


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
            logger.info(redactEmailAddress())
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
            tmp = paramO_2.toDict().copy()
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
            logger.info(redactEmailAddress())
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
            logger.info(redactEmailAddress())
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
            logger.info(redactEmailAddress())
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

