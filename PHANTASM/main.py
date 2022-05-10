# Author: Joseph S. Wirth

from PHANTASM.rRNA.runRnaBlast import rnaBlastRunner
from PHANTASM.rRNA.processRnaBlast import getTaxIdsFromRnaBlast
from PHANTASM.Parameter import Parameters
from PHANTASM.taxonomy.taxonomyConstruction import Taxonomy, constructTaxonomy, __getLpsnData
from PHANTASM.coreGenes import rankPhylogeneticMarkers, xenogiInterfacer_1, \
                        parseGenbank, allVsAllBlast, copyExistingBlastFiles, \
                        calculateCoreGenes, makeSpeciesTree
from PHANTASM.findMissingNeighbors import phyloMarkerBlastRunner, xenogiInterfacer_2
from PHANTASM.overallGenomeRelatedIndices import overallGenomeRelatedIndices, makeAaiHeatmap, makeAniHeatmap
from Bio import Entrez


def getPhyloMarker(allQueryGenbanksL:list, paramO:Parameters) -> None:
    """ run the entire first portion from start to finish
    """

    outgroup = taxonomyWrapper(allQueryGenbanksL, paramO)
    coreGenesWrapper_1(paramO)
    findPhylogeneticMarkersWrapper(allQueryGenbanksL, outgroup, paramO)


def refinePhylogeny(geneNum:int, allQueryGenbanksL:list, paramO_1:Parameters, \
                                                  paramO_2:Parameters) -> None:
    """ run the entire second portion from start to finish
    """
    outgroup = findMissingRelativesWrapper(geneNum, allQueryGenbanksL, \
                                                            paramO_1, paramO_2)
    coreGenesWrapper_2(paramO_1, paramO_2)
    finalAnalysesWrapper(outgroup, paramO_2)


def taxonomyWrapper(allQueryGenbanksL:list, paramO_1:Parameters) -> Taxonomy:
    """ creates a Taxonomy object, downloads gbffs, and makes the human map.
        returns the outgroup species as a Taxonomy object.
    """
    # set the entrez email address
    Entrez.email = paramO_1.email

    # get 16S rRNA sequences, create 16S db, and run BLASTn
    blastResultsFile = rnaBlastRunner(allQueryGenbanksL, paramO_1.workdir, \
                                                   paramO_1.blastExecutDirPath)

    # get the taxids from blastn table
    taxids = getTaxIdsFromRnaBlast(blastResultsFile)

    # construct a taxonomy object for the taxids
    taxO = constructTaxonomy(taxids, saveTax=True, dir=paramO_1.workdir)

    # make/download all files required for the first pass of xenoGI
    outgroup = xenogiInterfacer_1(taxO, allQueryGenbanksL, paramO_1)

    return outgroup


def coreGenesWrapper_1(paramO_1:Parameters) -> None:
    """ runs the first set of all-vs-all blastp comparisons and calculates the
        core genes
    """
    parseGenbank(paramO_1.toDict())
    allVsAllBlast(paramO_1.toDict())
    calculateCoreGenes(paramO_1)


def findPhylogeneticMarkersWrapper(allQryGbksL:list, outgroup:Taxonomy, \
                                                  paramO_1:Parameters) -> None:
    """ makes the first species tree and ranks the core genes as phylogenetic
        markers.
    """
    makeSpeciesTree(allQryGbksL, paramO_1, outgroup)
    rankPhylogeneticMarkers(paramO_1)


def findMissingRelativesWrapper(geneNumsL:list, allQueryGenbanksL:list, \
                         paramO_1:Parameters, paramO_2:Parameters) -> Taxonomy:
    """ uses blastp to find missing relatives, downloads a new set of genomes,
        and makes the new human map file. returns the outgroup species as a str
        (scientific name).
    """
    lpsnD = __getLpsnData()
    phyloMarkerBlastRunner(geneNumsL, paramO_1)
    outgroup = xenogiInterfacer_2(allQueryGenbanksL, paramO_1, paramO_2, lpsnD)

    return outgroup


def coreGenesWrapper_2(paramO_1:Parameters, paramO_2:Parameters) -> None:
    """ runs the first set of all-vs-all blastp comparisons and calculates the
        core genes
    """
    parseGenbank(paramO_2.toDict())
    copyExistingBlastFiles(paramO_1, paramO_2)
    allVsAllBlast(paramO_2.toDict())
    calculateCoreGenes(paramO_2)


def finalAnalysesWrapper(outgroup:Taxonomy, paramO_2:Parameters) -> None:
    """ makes the second species tree and calculates OGRIs.
    """
    makeSpeciesTree(paramO_2, outgroup)
    overallGenomeRelatedIndices(paramO_2)
    makeAaiHeatmap(paramO_2, outgroup)
    makeAniHeatmap(paramO_2, outgroup)



