
from PHANTASM.rRNA.runRnaBlast import rnaBlastRunner
from PHANTASM.rRNA.processRnaBlast import getTaxIdsFromRnaBlast
from PHANTASM.taxonomy.taxonomyConstruction import constructTaxonomy, __getLpsnData
from PHANTASM.coreGenes import rankPhylogeneticMarkers, xenogiInterfacer_1, \
                        parseGenbank, allVsAllBlast, copyExistingBlastFiles, \
                        calculateCoreGenes, makeSpeciesTree
from PHANTASM.findMissingNeighbors import phyloMarkerBlastRunner, xenogiInterfacer_2
from PHANTASM.overallGenomeRelatedIndices import overallGenomeRelatedIndices, makeAaiHeatmap
from Bio import Entrez


def runPt1(queryGbff:str, paramD:dict) -> None:
    """ run the entire first portion from start to finish
    """
    outgroupSciName = taxonomyWrapper(queryGbff, paramD)
    coreGenesWrapper_1(paramD)
    findPhylogeneticMarkersWrapper(outgroupSciName, paramD)


def runPt2(locusTag:str, queryGbff:str, paramD_1:dict, paramD_2:dict) -> None:
    """ run the entire second portion from start to finish
    """
    outgroupSciName = findMissingRelativesWrapper(locusTag, queryGbff, paramD_1, paramD_2)
    coreGenesWrapper_2(paramD_1, paramD_2)
    finalAnalysesWrapper(outgroupSciName)


def taxonomyWrapper(queryGenbank:str, paramD_1:dict) -> str:
    """ creates a Taxonomy object, downloads gbffs, and makes the human map.
        returns the outgroup species as a string (scientific name).
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
    outgroupSciName = xenogiInterfacer_1(taxO, queryGenbank, paramD_1)

    return outgroupSciName


def coreGenesWrapper_1(paramD_1:dict) -> None:
    """ runs the first set of all-vs-all blastp comparisons and calculates the
        core genes
    """
    parseGenbank(paramD_1)
    allVsAllBlast(paramD_1)
    calculateCoreGenes(paramD_1)


def findPhylogeneticMarkersWrapper(outgroupSciName:str, paramD_1:dict) -> None:
    """ makes the first species tree and ranks the core genes as phylogenetic
        markers.
    """
    makeSpeciesTree(paramD_1, outgroupSciName)
    rankPhylogeneticMarkers(paramD_1)


def findMissingRelativesWrapper(locusTag:str, queryGenbank:str, paramD_1:dict,\
                                                         paramD_2:dict) -> str:
    """ uses blastp to find missing relatives, downloads a new set of genomes,
        and makes the new human map file. returns the outgroup species as a str
        (scientific name).
    """
    lpsnD = __getLpsnData()
    phyloMarkerBlastRunner(locusTag, queryGenbank, paramD_1)
    outgroupSciName = xenogiInterfacer_2(queryGenbank, paramD_1, paramD_2, lpsnD)

    return outgroupSciName


def coreGenesWrapper_2(paramD_1:dict, paramD_2:dict) -> None:
    """ runs the first set of all-vs-all blastp comparisons and calculates the
        core genes
    """
    parseGenbank(paramD_2)
    copyExistingBlastFiles(paramD_1, paramD_2)
    allVsAllBlast(paramD_2)
    calculateCoreGenes(paramD_2)


def finalAnalysesWrapper(outgroupSciName:str, paramD_2:dict) -> None:
    """ makes the second species tree and calculates OGRIs.
    """
    makeSpeciesTree(paramD_2, outgroupSciName)
    overallGenomeRelatedIndices(paramD_2)
    makeAaiHeatmap(paramD_2, outgroupSciName)



