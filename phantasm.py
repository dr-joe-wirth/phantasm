# Author: Joseph S. Wirth

from PHANTASM.main import getPhyloMarker, refinePhylogeny
from PHANTASM.findMissingNeighbors import _locusTagToGeneNum
from PHANTASM.utilities import validEmailAddress, getParamO_1, getParamO_2
import glob, os, sys

# constants
GAP = " "*4
ERR_MSG_1 = "Incorrect syntax used.\nType 'python3 <path>/phantasm.py help' for more information."
ERR_MSG_2 = "Invalid email address"
ERR_MSG_3A = "Invalid task: "
ERR_MSG_3B = "\ntype '<path>/phantasm.py help' for information."
HELP_MSG = "\nPHANTASM: PHylogenomic ANalyses for the TAxonomy and Systematics of Microbes\n\n" + \
           "In order to run PHANTASM, first follow these steps:\n" + \
           GAP + "1. cd into the desired working directory\n" + \
           GAP + "2. identify a suitable phylogenetic marker:\n" + \
           GAP*2 + "'python3 <path>/phantasm.py getPhyloMarker <gbff file> <email address>'\n\n" + \
           GAP + "3. examine the file 'initialAnalysis/putativePhylogeneticMarkers.txt'\n" + \
           GAP + "4. refine the phylogeny and perform phylogenomic analyses:\n" + \
           GAP*2 + "'python3 <path>/phantasm.py refinePhylogeny --locus_tag <locus_tag> <gbff file> <email address>'\n" + \
           GAP*3 + "OR\n" + \
           GAP*2 + "'python3 <path>/phantasm.py refinePhylogeny --gene_num <gene_number> <gbff file> <email address>'\n" + \
           GAP + "5. results can be found in the following files:\n" + \
           GAP*2 + "'finalAnalysis/aai_matrix.txt'\n" + \
           GAP*2 + "'finalAnalysis/aai_heatmap.pdf'\n" + \
           GAP*2 + "'finalAnalysis/speciesTree.nwk'\n" + \
           GAP*2 + "'finalAnalysis/speciesTree_outgroupPruned.nwk'\n"
JOB_0 = 'help'
JOB_1 = 'getPhyloMarker'
JOB_2 = 'refinePhylogeny'

# begin main function
if __name__ == "__main__":
    # print help message if no arguments are provided
    if len(sys.argv) == 1:
        print(HELP_MSG)
    
    else:
        # extract the job name
        job = sys.argv[1]

        # print help message if requested
        if job == JOB_0:
            print(HELP_MSG)

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
            
            # construct the Parameters object and execute JOB_1
            paramO = getParamO_1(email)
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

            # construct the Parameters objects
            paramO_1 = getParamO_1(email)
            paramO_2 = getParamO_2(email)

            # if a locus tag was provided, convert it to an integer
            if flag == "--locus_tag":
                geneNum = _locusTagToGeneNum(gene, paramO_1.geneInfoFN)
            
            # otherwise, a gene number should have been specified
            elif flag == "--gene_num":
                geneNum = int(gene)
            
            else:
                raise ValueError("Invalid flag: " + flag)
            
            # if a directory was provided, then get the list of files within it
            if os.path.isdir(inputFile):
                gbffL = glob.glob(os.path.join(inputFile, "*"))
            
            # if a single file was provided, then convert it to a one-item list
            else:
                gbffL = [inputFile]

            refinePhylogeny(geneNum, gbffL, paramO_1, paramO_2)
        
        # raise an error if an invalid job was specified
        else:
            raise ValueError(ERR_MSG_3A + job + ERR_MSG_3B)

