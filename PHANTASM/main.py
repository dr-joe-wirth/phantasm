# Author: Joseph S. Wirth

from PHANTASM.rRNA.runRnaBlast import rnaBlastRunner
from PHANTASM.rRNA.processRnaBlast import getTaxIdsFromRnaBlast
from PHANTASM.taxonomy.taxonomyConstruction import Taxonomy, constructTaxonomy, __getLpsnData
from PHANTASM.coreGenes import rankPhylogeneticMarkers, xenogiInterfacer_1, \
                        parseGenbank, allVsAllBlast, copyExistingBlastFiles, \
                        calculateCoreGenes, makeSpeciesTree
from PHANTASM.findMissingNeighbors import phyloMarkerBlastRunner, xenogiInterfacer_2
from PHANTASM.overallGenomeRelatedIndices import overallGenomeRelatedIndices, makeAaiHeatmap
from Bio import Entrez


def runPt1(queryGbff:str, paramD:dict) -> None:
    """ run the entire first portion from start to finish
    """
    outgroup = taxonomyWrapper(queryGbff, paramD)
    coreGenesWrapper_1(paramD)
    findPhylogeneticMarkersWrapper(outgroup, paramD)


def runPt2(geneNum:int, queryGbff:str, paramD_1:dict, paramD_2:dict) -> None:
    """ run the entire second portion from start to finish
    """
    outgroup = findMissingRelativesWrapper(geneNum, queryGbff, paramD_1, paramD_2)
    coreGenesWrapper_2(paramD_1, paramD_2)
    finalAnalysesWrapper(outgroup, paramD_2)


def taxonomyWrapper(queryGenbank:str, paramD_1:dict) -> Taxonomy:
    """ creates a Taxonomy object, downloads gbffs, and makes the human map.
        returns the outgroup species as a Taxonomy object.
    """
    # set the entrez email address
    Entrez.email = paramD_1['email']

    # get 16S rRNA sequences, create 16S db, and run BLASTn
    blastResultsFile = rnaBlastRunner(queryGenbank, paramD_1['workdir'], paramD_1['blastExecutDirPath'])

    # get the taxids from blastn table
    taxids = getTaxIdsFromRnaBlast(blastResultsFile)

    # construct a taxonomy object for the taxids
    taxO = constructTaxonomy(taxids, saveTax=True, dir=paramD_1['workdir'])

    # make/download all files required for the first pass of xenoGI
    outgroup = xenogiInterfacer_1(taxO, queryGenbank, paramD_1)

    return outgroup


def coreGenesWrapper_1(paramD_1:dict) -> None:
    """ runs the first set of all-vs-all blastp comparisons and calculates the
        core genes
    """
    parseGenbank(paramD_1)
    allVsAllBlast(paramD_1)
    calculateCoreGenes(paramD_1)


def findPhylogeneticMarkersWrapper(outgroup:Taxonomy, paramD_1:dict) -> None:
    """ makes the first species tree and ranks the core genes as phylogenetic
        markers.
    """
    makeSpeciesTree(paramD_1, outgroup)
    rankPhylogeneticMarkers(paramD_1)


def findMissingRelativesWrapper(geneNum:int, queryGenbank:str, paramD_1:dict,\
                                                    paramD_2:dict) -> Taxonomy:
    """ uses blastp to find missing relatives, downloads a new set of genomes,
        and makes the new human map file. returns the outgroup species as a str
        (scientific name).
    """
    lpsnD = __getLpsnData()
    phyloMarkerBlastRunner(geneNum, paramD_1)
    outgroup = xenogiInterfacer_2(queryGenbank, paramD_1, paramD_2, lpsnD)

    return outgroup


def coreGenesWrapper_2(paramD_1:dict, paramD_2:dict) -> None:
    """ runs the first set of all-vs-all blastp comparisons and calculates the
        core genes
    """
    parseGenbank(paramD_2)
    copyExistingBlastFiles(paramD_1, paramD_2)
    allVsAllBlast(paramD_2)
    calculateCoreGenes(paramD_2)


def finalAnalysesWrapper(outgroup:Taxonomy, paramD_2:dict) -> None:
    """ makes the second species tree and calculates OGRIs.
    """
    makeSpeciesTree(paramD_2, outgroup)
    overallGenomeRelatedIndices(paramD_2)
    makeAaiHeatmap(paramD_2, outgroup)



