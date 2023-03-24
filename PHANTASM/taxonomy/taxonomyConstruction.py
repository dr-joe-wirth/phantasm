# Author: Joseph S. Wirth

from PHANTASM.taxonomy.Taxonomy import *
from PHANTASM.taxonomy.processLpsn import makeLpsnD
from PHANTASM.utilities import getParentalTaxIds
import logging


def constructTaxonomy(taxids:list, saveTax:bool=False, dir:str='./') -> Taxonomy:
    """ constructTaxonomy:
            Accepts a list of NCBI taxonomy ids, a boolean indicating whether
            or not to save the resulting Taxonomy object to a file, and the di-
            rectory where the file should be saved as inputs. Identifies the
            order(s) represented by the input ids, and makes a Taxonomy object
            for each order. Populates the order(s) to the species-level and re-
            conciles their taxonomy with LPSN data. If multiple orders were id-
            entified, then finds their taxonomic MRCA and nests the orders wit-
            hin it. Uses the taxids to find and label the taxa that could be
            phylogenomic outgroups (external) or ingroups (internal). Imports
            data from NCBI assembly. Returns the newly built Taxonomy object.
    """
    # constants
    GAP = " " * 4
    PRNT_1 = "Constructing taxonomy for the query genome(s) ... "
    PRNT_2 = GAP + "Reconciling taxonomy with LPSN taxonomy ... "
    PRNT_3 = GAP + "Merging multiple orders into a single taxonomy ... "
    PRNT_4 = GAP + "Mapping blast hits onto the taxonomy ... "
    PRNT_5 = GAP + "Importing data from NCBI Assembly ... "
    DONE   = "Done."
    ERR_MSG_1 = "Taxonomic order could not be identified from provided taxids."
    ERR_MSG_2 = "Failed to create a Taxonomy object. Aborting"
    ERR_MSG_3 = "Resulting taxonomy is not consistent. Aborting."
    
    # initialize logger
    logger = logging.getLogger(__name__ + ".constructTaxonomy")

    # retrieve a list of Taxonomy objects at the order-level 
    # that are represented by the input taxids
    print(PRNT_1)
    logger.info(PRNT_1)
    allOrdersL = __initializeTaxonomy(taxids)

    # raise an error if no orders were found
    if len(allOrdersL) == 0:
        logger.error(ERR_MSG_1)
        raise RuntimeError(ERR_MSG_1)

    print(PRNT_2, end='', flush=True)
    logger.info(PRNT_2)
    lpsnD = _getLpsnData()

    # initialize a dictionary to resolve colliding root taxids
    rootTaxD = dict()

    # for each order ...
    for idx in range(len(allOrdersL)):
        # reconcile the object with LPSN data
        #### note: this step can produce roots that will collide
        order:Taxonomy = allOrdersL[idx]
        order._addLpsnData(lpsnD)

        # create a new entry for the taxid if it is not in the dictionary
        if order.getRoot().taxid not in rootTaxD.keys():
            rootTaxD[order.getRoot().taxid] = []

        # store the current index under the taxid 
        rootTaxD[order.getRoot().taxid].append(idx)

        # make sure all objects are at their roots
        allOrdersL[idx] = order.getRoot()

    # initialize a list of indices that will need to be removed from allOrdersL
    indicesToRemoveL = list()
    
    # check the taxids in in the dictionary
    for taxid in rootTaxD.keys():
        # get the last index seen for this taxid
        idx = rootTaxD[taxid].pop()

        # use the index to get the first seen taxonomy object
        taxO:Taxonomy = allOrdersL[idx]

        # resolve any roots whose taxids are colliding
        # this loop will only be entered if collisions will occur
        for idx in rootTaxD[taxid]:
            # get the offending  taxonomy
            nextO:Taxonomy = allOrdersL[idx]

            # resolve the collisions
            __resolveCollidingRoots(taxO, nextO)

            # mark the index as flagged for removal
            indicesToRemoveL.append(idx)

    # sort the indices from high-to-low for on-the-fly popping
    indicesToRemoveL = sorted(indicesToRemoveL, reverse=True)

    # remove the specificied objects from allOrdersL
    for idx in indicesToRemoveL:
        allOrdersL.pop(idx) 

    print(DONE)
    logger.info(GAP + DONE)
    
    # if only one order was made, then use it
    if len(allOrdersL) == 1:
        taxonomy = allOrdersL.pop()
    
    # if multiple orders were made, then merge them into a single object
    elif len(allOrdersL) > 1:
        print(PRNT_3, end='', flush=True)
        logger.info(PRNT_3)
        taxonomy = Taxonomy._mergeMultipleRoots(allOrdersL, lpsnD)
        print(DONE)
        logger.info(GAP + DONE)
    
    # otherwise, no orders were made; raise an error
    else:
        logger.error(ERR_MSG_2)
        raise BaseException(ERR_MSG_2)
    
    # use the taxids to update which taxa are considered interal/external
    print(PRNT_4, end="", flush=True)
    logger.info(PRNT_4)
    taxonomy._updateInternalExternalStatus(taxids, lpsnD)
    print(DONE)
    logger.info(GAP + DONE)

    # import data from NCBI assembly for the object
    print(PRNT_5, end='', flush=True)
    logger.info(PRNT_5)
    taxonomy._importAssemblies()
    print(DONE)
    logger.info(GAP + DONE)

    # make sure the resulting object is consistent
    if not taxonomy._isConsistent():
        logger.error(ERR_MSG_3)
        raise RuntimeError(ERR_MSG_3)

    # write object to file (if requested)
    if saveTax:
        __saveTaxonomy(taxonomy, dir)
    
    print(DONE + "\n")
    logger.info(DONE)

    return taxonomy


def __initializeTaxonomy(taxids:list) -> list:
    """ initializeTaxonomy:
            Accepts a list of NCBI taxonomy ids as input. Finds the taxonomic 
            order(s) represented by the input ids. For each order found, const-
            ructs a Taxonomy object at the order level and fully populates it to
            the species-level. Adds each order to a list
    """
    # constants
    GAP = " " * 4
    PRNT_1 = GAP + "Identifying order(s) represented in the blast results ... "
    PRNT_2 = GAP + "Importing data from NCBI Taxonomy for "
    PRNT_3 = " orders:"
    DONE = 'Done.'
    ORDER  = 'order'
    DOMAIN = 'domain'
    ERR_MSG = "Blast results contain taxids from both Bacteria and Archaea" + \
              ". Please select a different phylogenetic marker(s)."

    # initialize a logger
    logger = logging.getLogger(__name__ + ".__initializeTaxonomy")

    # make sure that all orders are in the same domain
    domInfo = getParentalTaxIds(taxids, DOMAIN)
    if len(domInfo) > 1:
        logger.error(ERR_MSG)
        raise RuntimeError(ERR_MSG)

    # find the order(s) represented by the blastn
    print(PRNT_1, end='', flush=True)
    logger.info(PRNT_1)
    ordInfo = getParentalTaxIds(taxids, ORDER)
    print(DONE)
    logger.info(GAP + DONE)

    # print the status
    status = PRNT_2 + str(len(ordInfo)) + PRNT_3
    logger.info(status)
    print(status)

    # for each order found ...
    allOrders = list()
    for order in ordInfo:
        # make a taxonomy object and add it to the list of orders
        taxO = __getOrder(order)
        allOrders.append(taxO)
    
    return allOrders


def __getOrder(ordInfo:dict) -> Taxonomy:
    """ getOrder:
            Accepts a single dictionary with the keys 'name' and 'txid' for a
            taxonomic order as input. Makes a new Taxonomy object for the order
            and populates it to the species-level. Returns the Taxonomy
    """
    # constants
    GAP = " " * 8
    PRNT_2 = " (taxid: "
    PRNT_3 = ") ... "
    DONE = "Done."
    RANK = 'order'
    
    logger = logging.getLogger(__name__ + ".__getOrder")
    
    # print the status
    status = GAP + ordInfo['name'] + PRNT_2 + ordInfo['txid'] + PRNT_3
    logger.info(status)
    print(status, end="", flush=True)

    # create a taxonomy object for the order
    order =  Taxonomy(ordInfo['txid'], ordInfo['name'], RANK, isExternal=False)
    order._initializeDescendants()
    print(DONE)
    logger.info(GAP + DONE)
    
    return order


def _getLpsnData() -> dict:
    """ getLpsnData:
            Accepts no inputs. Constructs and returns the LPSN data dictionary.

            This needs to be retired once an API is available.
    """
    # constants
    from param import CSV_1, CSV_2, CSV_3, CSV_4

    # import LPSN data
    return makeLpsnD(CSV_1, CSV_2, CSV_3, CSV_4)


def __resolveCollidingRoots(taxO:Taxonomy, otherO:Taxonomy) -> None:
    """ resolveCollidingRoots:
            Accepts two taxonomy objects whose taxids are equal. Resolves any
            collisions by coercing the second object's data into the first obj-
            ect. Accomplishes this by recursing along shared children of the 
            two objects. Modifies both inputs. Does not return.
    """
    # constants
    ERR_MSG_1 = "objects do not share taxids"
    ERR_MSG_2 = "untested condition: artificial taxids are colliding"
    ERR_MSG_3 = "untested condition: taxid present, but name doesn't exist"
    ERR_MSG_4 = "untested condition: synonym names don't match"

    logger = logging.getLogger(__name__ + ".__resolveCollidingRoots")

    # ensure that the taxids of the two objects are identical
    if taxO.taxid != otherO.taxid:
        logger.error(ERR_MSG_1)
        raise RuntimeError(ERR_MSG_1)

    # only do work if objects are inequivalent
    if taxO != otherO:
        # for each child of the second object
        otherChild:Taxonomy
        for otherChild in otherO.getChildren(list):
            # if the taxid is not already present in taxO
            if otherChild.taxid not in taxO:
                # if the scientific name is not already present in taxO
                if otherChild.sciName not in taxO:
                    # simply import the descendant; no collisions will occur
                    taxO._importDirectDescendant(otherChild)
                
                # if the scientific name is already present in taxO
                else:
                    # get the synonymous child from taxO
                    taxChild = taxO.getDescendantBySciName(otherChild.sciName)

                    # artificial ids colliding here is untested. raise an error
                    if int(taxChild.taxid) < 0 or int(otherChild.taxid) < 0:
                        logger.error(ERR_MSG_2)
                        raise RuntimeError(ERR_MSG_2)

                    # import the descendant, then resolve the duplicate names
                    taxO._importDirectDescendant(otherChild)

            # taxid present with the name absent is an untested condition
            elif otherChild.sciName not in taxO:
                logger.error(ERR_MSG_3)
                raise RuntimeError(ERR_MSG_3)

            # if taxid is present in taxO 
            else:
                # get the synonymous child from taxO
                taxChild = taxO.getDescendantByTaxId(otherChild.taxid)

                # make sure the names match; error if not
                if taxChild.sciName != otherChild.sciName:
                    logger.error(ERR_MSG_4)
                    raise RuntimeError(ERR_MSG_4)
                
                # recurse on the two synonymous children
                __resolveCollidingRoots(taxChild, otherChild)


def __saveTaxonomy(taxO:Taxonomy, dir:str) -> None:
    """ saveTaxonomy:
            Accepts a Taxonomy object and a string indicating the desired save
            directory as inputs. Saves the object to a file with the scientific
            name of its root as the file name. Does not return.
    """
    # constants
    EXTENSION = ".tax"

    # make sure directory ends in a slash
    if dir[-1] != '/':
        dir += '/'

    # determine the file name
    filename = dir + taxO.getRoot().sciName + EXTENSION

    # save the file
    taxO.save(filename)


