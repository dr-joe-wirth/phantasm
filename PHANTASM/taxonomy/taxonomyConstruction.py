# Author: Joseph S. Wirth

from PHANTASM.taxonomy.Taxonomy import *
from PHANTASM.taxonomy.processLpsn import makeLpsnD
from PHANTASM.utilities import getParentalTaxIds


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
    PRNT_1 = "Constructing taxonomy for the query genome ... "
    PRNT_2 = GAP + "Reconciling taxonomy with LPSN taxonomy ... "
    PRNT_3 = GAP + "Merging multiple orders into a single taxonomy ... "
    PRNT_4 = GAP + "Mapping blast hits onto the taxonomy ... "
    PRNT_5 = GAP + "Importing data from NCBI Assembly ... "
    DONE   = "Done."

    # retrieve a list of Taxonomy objects at the order-level 
    # that are represented by the input taxids
    print(PRNT_1)
    allOrders = __initializeTaxonomy(taxids)

    # for each order ...
    print(PRNT_2, end='', flush=True)
    lpsnD = __getLpsnData()
    for idx in range(len(allOrders)):
        # ... reconcile the object with LPSN data
        order:Taxonomy = allOrders[idx]
        order._addLpsnData(lpsnD)

        # ... make sure the object is at its root
        allOrders[idx] = order.getRoot()
    print(DONE)
    
    # if only one order was made, then use it
    if len(allOrders) == 1:
        taxonomy = allOrders.pop()
    
    # if multiple orders were made, then merge them into a single object
    elif len(allOrders) > 1:
        print(PRNT_3, end='', flush=True)
        taxonomy = Taxonomy._mergeMultipleRoots(allOrders, lpsnD)
        print(DONE)
    
    # otherwise, no orders were made; raise an error
    else: raise BaseException("Failed to create a Taxonomy object; aborting")
    
    # use the taxids to update which taxa are considered interal/external
    print(PRNT_4, end="", flush=True)
    taxonomy._updateInternalExternalStatus(taxids, lpsnD)
    print(DONE)

    # import data from NCBI assembly for the object
    print(PRNT_5, end='', flush=True)
    taxonomy._importAssemblies()
    print(DONE)

    # make sure the resulting object is consistent
    if not taxonomy._isConsistent():
        raise RuntimeError("Resulting taxonomy is not consistent. Aborting.")

    # write object to file (if requested)
    if saveTax:
        __saveTaxonomy(taxonomy, dir)
    
    print(DONE + "\n")

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

    # find the order(s) represented by the blastn
    print(PRNT_1, end='', flush=True)
    ordInfo = getParentalTaxIds(taxids, ORDER)
    print(DONE)

    # print the status
    status = PRNT_2 + str(len(ordInfo)) + PRNT_3
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
    
    # print the status
    status = GAP + ordInfo['name'] + PRNT_2 + ordInfo['txid'] + PRNT_3
    print(status, end="", flush=True)

    # create a taxonomy object for the order
    order =  Taxonomy(ordInfo['txid'], ordInfo['name'], RANK, isExternal=False)
    order._initializeDescendants()
    print(DONE)
    
    return order


def __getLpsnData() -> dict:
    """ __getLpsnData:
            Accepts no inputs. Constructs and returns the LPSN data dictionary.

            This needs to be retired once an API is available.
    """
    # constants
    from param import CSV_1, CSV_2, CSV_3, CSV_4

    # import LPSN data
    return makeLpsnD(CSV_1, CSV_2, CSV_3, CSV_4)


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


