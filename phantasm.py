# Author: Joseph S. Wirth


from PHANTASM.main import getPhyloMarker, refinePhylogeny, knownPhyloMarker, analyzeSpecifiedGenomes, rankPhyloMarkers
from PHANTASM.utilities import parseArgs, getLpsnAge, redactEmailAddress, getHelpMessage
from PHANTASM.findMissingNeighbors import _locusTagToGeneNum
from PHANTASM.Parameter import Parameters
from param import PHANTASM_DIR
import logging, os, sys

# constants
PHANTASM_PY = "python3 " + os.path.join(PHANTASM_DIR, "phantasm.py")
VERSION = "v1.1.0"
JOB_0A = "help"
JOB_0B = "version"
JOB_1 = 'getPhyloMarker'
JOB_2 = 'refinePhylogeny'
JOB_3 = 'knownPhyloMarker'
JOB_4 = 'analyzeGenomes'
JOB_5 = 'rankPhyloMarkers'
ALL_JOBS = (JOB_0A, JOB_0B, JOB_1, JOB_2, JOB_3, JOB_4, JOB_5)
SEP = ","
GAP = " "*4

# reference message
REF_MSG = "\nIf you use this software in your research, please cite our paper:\n" + \
          GAP + "Automating microbial taxonomy workflows with PHANTASM: PHylogenomic\n" + \
          GAP + "ANalyses for the TAxonomy and Systematics of Microbes\n" + \
          GAP*2 + "Joseph S. Wirth and Eliot C. Bush, 2023\n" + \
          GAP*2 + "https://doi.org/10.1093/nar/gkad196\n"

# version message
VERSION_MSG = "PHANTASM " + VERSION + "\n"

# no arguments message
NO_ARGS_MSG = "\ntype '" + PHANTASM_PY + " help' for help\n"

# LPSN age message
AGE_MSG = "PHANTASM relies on data manually acquired from the LPSN.\n" + \
          GAP + "These data were last retrieved on "

# error messages
ERR_MSG_1A = "Invalid task: "
ERR_MSG_1B = "\n" + NO_ARGS_MSG
ERR_MSG_2 = "One or more specified genes could not be extracted."
ERR_MSG_3 = "The specified output directory already exists."

# begin main function
if __name__ == "__main__":
    # print the reference message first
    print(REF_MSG)

    # print a message if no arguments are provided
    if len(sys.argv) == 1:
        print(VERSION_MSG + NO_ARGS_MSG)
    
    else:
        # extract the job name
        job = sys.argv[1]

        # make sure a valid job was specified
        if job not in ALL_JOBS:
            raise ValueError(ERR_MSG_1A + job + ERR_MSG_1B)

        # print the help message for the `help` job
        elif job == JOB_0A:
            print(getHelpMessage(job))

        # print the version if requested
        elif job == JOB_0B:
            print(VERSION_MSG)

        # print the help message for the specified job
        elif "-h" in sys.argv or "--help" in sys.argv:
            print(getHelpMessage(job))

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
            logger.info('lpsn data:       ' + age)
            logger.info('num cpus:        ' + str(paramO.numProcesses))
            logger.info('max leaves:      ' + str(paramO.maxNumTreeLeaves))
            logger.info('reduce num core: ' + str(paramO.reduceNumCoreGenes))

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
                raise ValueError(ERR_MSG_2)

            # print the LPSN age data
            age = getLpsnAge()
            print(AGE_MSG + age + "\n\n")
            
            # initialize logger
            logging.basicConfig(filename=paramO_1.logFN, level=logging.INFO)
            logger = logging.getLogger(__name__)
            
            # save details about the run
            logger.info(redactEmailAddress())
            logger.info(VERSION)
            logger.info('lpsn data:       ' + age)
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
            logger.info('lpsn data:       ' + age)
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
                raise FileExistsError(ERR_MSG_3)

            # initialize logger
            logging.basicConfig(filename=paramO.logFN, level=logging.INFO)
            logger = logging.getLogger(__name__)
            
            # save details about the run
            logger.info(redactEmailAddress())
            logger.info(VERSION)
            logger.info('num cpus:        ' + str(paramO.numProcesses))
            logger.info('reduce num core: ' + str(paramO.reduceNumCoreGenes))
            logger.info('bootstrap tree:  ' + str(paramO.numBootstraps > 0))
            logger.info('num bootstraps:  ' + str(paramO.numBootstraps) + "\n")
            
            # analyze the specified genomes
            logger.info("start " + JOB_4 + "\n")
            analyzeSpecifiedGenomes(gbffL, paramO)
            logger.info("end " + JOB_4 + "\n")

        elif job == JOB_5:
            gbffL, locusTagsL, paramO = parseArgs()
            
            # initialize logger
            logging.basicConfig(filename=paramO.logFN, level=logging.INFO)
            logger = logging.getLogger(__name__)
            
            # save details about the run
            logger.info(VERSION)
            logger.info('num cpus: ' + str(paramO.numProcesses))
            
            # rank the phylogenetic markers
            logger.info("start " + JOB_5 + "\n")
            rankPhyloMarkers(paramO)
            logger.info("end " + JOB_5 + "\n")
