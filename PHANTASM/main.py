# Author: Joseph S. Wirth
# Last edit: September 27, 2022

import os
from Bio import Entrez
from PHANTASM.Parameter import Parameters
from PHANTASM.rRNA.runRnaBlast import rnaBlastRunner
from PHANTASM.rRNA.processRnaBlast import getTaxIdsFromRnaBlast
from PHANTASM.utilities import cleanup, getTaxidsFromFile, parseCsv
from PHANTASM.taxonomy.taxonomyConstruction import Taxonomy, constructTaxonomy, _getLpsnData
from PHANTASM.findMissingNeighbors import phyloMarkerBlastRunner, xenogiInterfacer_2, xenogiInterfacer_3
from PHANTASM.overallGenomeRelatedIndices import overallGenomeRelatedIndices, makeAaiHeatmap, makeAniHeatmap
from PHANTASM.coreGenes import rankPhylogeneticMarkers, xenogiInterfacer_1, parseGenbank, allVsAllBlast, copyExistingBlastFiles, calculateCoreGenes, makeSpeciesTree


def getPhyloMarker(allQueryGenbanksL:list, paramO_1:Parameters) -> None:
    """ run the getPhyloMarker process
    """
    # skip the 16s sequences if taxids file is present
    if os.path.exists(paramO_1.taxidsFN):
        outgroup = skip16sWrapper(allQueryGenbanksL, paramO_1)

    # use 16s if taxids file is absent
    else:
        outgroup = taxonomyWrapper(allQueryGenbanksL, paramO_1)

    # rank the phylogenetic markers
    coreGenesWrapper_1(paramO_1)
    findPhylogeneticMarkersWrapper(allQueryGenbanksL, outgroup, paramO_1)
    cleanup(paramO_1)


def refinePhylogeny(geneNumsL:list, allQueryGenbanksL:list, \
                             paramO_1:Parameters, paramO_2:Parameters) -> None:
    """ run the refinePhylogeny job (only used after getPhyloMarker)
    """
    outgroup = findMissingRelativesWrapper(geneNumsL, allQueryGenbanksL, \
                                                            paramO_1, paramO_2)
    coreGenesWrapper_2(paramO_1, paramO_2)
    finalAnalysesWrapper(allQueryGenbanksL, outgroup, paramO_2)
    cleanup(paramO_2)


def knownPhyloMarker(allQueryGenbanksL:list, locusTagsL:list, \
                                                    paramO:Parameters) -> None:
    """ run phantasm with known markers in a single step
    """
    lpsnD = _getLpsnData()
    outgroup = xenogiInterfacer_3(allQueryGenbanksL, locusTagsL, paramO, lpsnD)
    coreGenesWrapper_1(paramO)
    finalAnalysesWrapper(allQueryGenbanksL, outgroup, paramO)
    cleanup(paramO)


def analyzeSpecifiedGenomes(allGbksL:list, paramO:Parameters) -> None:
    """ use phantasm to analyze a set of user-specified genomes
    """
    coreGenesWrapper_1(paramO)
    outgroup = makeOutgroupObject(paramO)
    finalAnalysesWrapper(allGbksL, outgroup, paramO)
    cleanup(paramO)
###############################################################################


###############################################################################
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


def skip16sWrapper(allQueryGenbanksL:list, paramO_1:Parameters) -> Taxonomy:
    """ creates a Taxonomy object, downloads gbffs, and makes the human map.
        returns the outgroup species as a Taxonomy object.
    """
    # set the entrez email address
    Entrez.email = paramO_1.email

    # create the working directory
    os.mkdir(paramO_1.workdir)

    # extract the taxids from the file
    taxids = getTaxidsFromFile(paramO_1.taxidsFN)

    # constrcut a taxonomy object for the taxids
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
    makeSpeciesTree(allQryGbksL, paramO_1, outgroup, 'fasttree')
    rankPhylogeneticMarkers(paramO_1)


def findMissingRelativesWrapper(geneNumsL:list, allQueryGenbanksL:list, \
                         paramO_1:Parameters, paramO_2:Parameters) -> Taxonomy:
    """ uses blastp to find missing relatives, downloads a new set of genomes,
        and makes the new human map file. returns the outgroup species as a str
        (scientific name).
    """
    lpsnD = _getLpsnData()
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


def makeOutgroupObject(paramO:Parameters) -> Taxonomy:
    """ makeOutgroupObject:
            Accepts a Parameters object as input. Constructs a Taxonomy object
            for the last entry in the human map file. The object uses the rese-
            rved artificial taxid. This function should only be used when users
            specify their own reference genomes instead of allowing PHANTASM to
            pick for them. Returns the newly constructed Taxonomy object.
    """
    # constants
    ARTIFICIAL_TAXID = "0"
    DELIM = "\t"
    
    # read the human map file
    humanMapL = parseCsv(paramO.fileNameMapFN, delim=DELIM)

    # truncate empty lines from the human map file
    while len(humanMapL[-1]) == 0:
        humanMapL.pop()
    
    # the outgroup will always be the last entry in the human map file
    outgroupName = humanMapL[-1][1]

    # create and return a Taxonomy object for the outgroup
    return Taxonomy(ARTIFICIAL_TAXID, outgroupName, 'species')


def finalAnalysesWrapper(allQryGbksL:list, outgroup:Taxonomy, \
                                                  paramO_2:Parameters) -> None:
    """ makes the second species tree and calculates OGRIs.
    """
    from param import BOOTSTRAP_FINAL_TREE

    if BOOTSTRAP_FINAL_TREE:
        makeSpeciesTree(allQryGbksL, paramO_2, outgroup, 'iqtree')
    else:
        makeSpeciesTree(allQryGbksL, paramO_2, outgroup, 'fasttree')
    overallGenomeRelatedIndices(paramO_2, outgroup)
    makeAaiHeatmap(paramO_2, outgroup)
    makeAniHeatmap(paramO_2, outgroup)
