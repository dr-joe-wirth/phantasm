from PHANTASM.Parameter import Parameters
from PHANTASM.main import runPt1, runPt2
from PHANTASM.findMissingNeighbors import _locusTagToGeneNum
from PHANTASM.utilities import validEmailAddress, getParamD_1, getParamD_2
import sys


GAP = " "*4
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

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print(HELP_MSG)
    
    else:
        job = sys.argv[1]

        if job == JOB_0:
            print(HELP_MSG)

        elif job == JOB_1:
            if len(sys.argv) != 4:
                raise SyntaxError("Incorrect syntax used.\nType 'python3 <path>/phantasm.py help' for more information.")
            gbffFN = sys.argv[2]
            email = sys.argv[3]

            if not validEmailAddress(email):
                raise ValueError("Invalid email address")

            paramD = getParamD_1(email)
            runPt1(gbffFN, paramD)
        
        elif job == JOB_2:
            if len(sys.argv) != 6:
                raise SyntaxError("Incorrect syntax used.\nType 'python3 <path>/phantasm.py help' for more information.")
            flag = sys.argv[2]
            gene = sys.argv[3]
            gbffFN = sys.argv[4]
            email = sys.argv[5]

            if not validEmailAddress(email):
                raise ValueError("Invalid email address")

            paramD_1 = getParamD_1(email)
            paramD_2 = getParamD_2(email)

            if flag == "--locus_tag":
                geneNum = _locusTagToGeneNum(gene, paramD_1['geneInfoFN'])
            
            elif flag == "--gene_num":
                geneNum = gene
            
            else:
                raise ValueError("Invalid flag: " + flag)
            
            runPt2(geneNum, gbffFN, paramD_1, paramD_2)
        
        else:
            raise ValueError("Invalid task: " + job + "\ntype '<path>/phantasm.py help' for information.")

