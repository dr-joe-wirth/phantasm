# Author: Joseph S. Wirth

from PHANTASM.main import getPhyloMarker, refinePhylogeny, knownPhyloMarker
from PHANTASM.findMissingNeighbors import _locusTagToGeneNum
from PHANTASM.utilities import validEmailAddress, getParamO_1, getParamO_2
import glob, os, sys

# constants
JOB_0 = 'help'
JOB_1 = 'getPhyloMarker'
JOB_2 = 'refinePhylogeny'
JOB_3 = 'knownPhyloMarker'
SEP = ","
GAP = " "*4
LOCUS_TAG = "--locus_tag"
GENE_NUM  = "--gene_num"

# error messages
ERR_MSG_1 = "Incorrect syntax used.\nType 'python3 <path>/phantasm.py help' for more information."
ERR_MSG_2 = "Invalid email address"
ERR_MSG_3 = "The number of genes is not a multiple of the number of input genomes."
ERR_MSG_4 = "Invalid flag: "
ERR_MSG_5A = "Invalid task: "
ERR_MSG_5B = "\ntype '<path>/phantasm.py help' for information."

# help message
HELP_MSG = "\nPHANTASM: PHylogenomic ANalyses for the TAxonomy and Systematics of Microbes\n\n" + \
           "In order to run PHANTASM without a known phylogenetic marker, follow these steps:\n" + \
           GAP + "1. cd into the desired working directory\n" + \
           GAP + "2. identify a suitable phylogenetic marker:\n\n" + \
           GAP*2 + "'python3 <path>/phantasm.py getPhyloMarker <gbff file> <email address>'\n\n" + \
           GAP + "3. examine the file 'initialAnalysis/putativePhylogeneticMarkers.txt'\n" + \
           GAP + "4. refine the phylogeny and perform phylogenomic analyses:\n\n" + \
           GAP*2 + "'python3 <path>/phantasm.py refinePhylogeny --locus_tag <locus_tag> <gbff file> <email address>'\n" + \
           GAP*3 + "OR\n" + \
           GAP*2 + "'python3 <path>/phantasm.py refinePhylogeny --gene_num <gene_number> <gbff file> <email address>'\n\n" + \
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
           GAP*2 + "'python3 <path>/phantasm.py knownPhyloMarker <locus tag> <gbff file> <email address>'\n\n" + \
           GAP + "3. results can be found in the following files:\n\n" + \
           GAP*2 + "'finalAnalysis/aai_matrix.txt'\n" + \
           GAP*2 + "'finalAnalysis/aai_heatmap.pdf'\n" + \
           GAP*2 + "'finalAnalysis/ani_matrix.txt'\n" + \
           GAP*2 + "'finalAnalysis/ani_heatmap.pdf'\n" + \
           GAP*2 + "'finalAnalysis/coreGenesSummary.txt\n" + \
           GAP*2 + "'finalAnalysis/speciesTree.nwk'\n" + \
           GAP*2 + "'finalAnalysis/speciesTree_outgroupPruned.nwk'\n\n" + \
           "Optional Features\n" + \
           GAP + "multiple input genomes:\n" + \
           GAP*2 + "replace '<gbff file>' with '<gbff directory>'\n" + \
           GAP*2 + "the specified directory must contain all of the input genome files\n\n" + \
           GAP + "multiple phylogenetic markers:\n" + \
           GAP*2 + "list each desired marker separated by a comma (no spaces)\n" + \
           GAP*2 + "the number of specified markers must be a multiple of the number input genomes\n\n" + \
           GAP + "excluding specific taxa from the final analysis:\n" + \
           GAP*2 + "create a file in the working directory named 'excludedTaxids.txt'\n" + \
           GAP*2 + "the file should contain exactly one NCBI Taxonomy taxid to be excluded per line\n" + \
           GAP*2 + "WARNING: this feature is experimental and should be used with caution.\n"

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
            
            # execute JOB_3
            knownPhyloMarker(gbffL, locusTagsL, paramO)

        # raise an error if an invalid job was specified
        else:
            raise ValueError(ERR_MSG_5A + job + ERR_MSG_5B)

