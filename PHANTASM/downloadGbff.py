# Author: Joseph S. Wirth

import re, os
from PHANTASM.utilities import downloadFileFromFTP, removeFileExtension, decompressGZ
from PHANTASM.taxonomy.Taxonomy import Taxonomy
from PHANTASM.taxonomy.taxonomyConstruction import __getLpsnData


def downloadGbffsForRootTaxonomy(taxO:Taxonomy, paramD:dict) -> Taxonomy:
    """ downloadGbffsForRootTaxonomy:
            Accepts a Taxonomy object and the parameters dictionary as inputs.
            Uses the Taxonomy object to determine which assemblies to download
            and then downloads those assemblies. Returns the outgroup as a Tax-
            onomy object.
    """
    # constants
    PRINT_1 = 'Identifying a suitable set of whole genome sequences ... '
    PRINT_2 = 'Downloading genbank files from NCBI ... '
    DONE = 'Done.'

    maxNumSeqs = paramD['maxNumTreeLeaves'] - 2  # user_input + outgroup = 2
    humanMapFN = paramD['fileNameMapFN']
    gbffFN = paramD['genbankFilePath']

    gbffDir = os.path.dirname(gbffFN)

    # get lpsnD
    lpsnD = __getLpsnData()

    # get a list of species for phylogenetic anlayses
    print(PRINT_1, end='', flush=True)
    speciesList = __selectSpeciesFromTaxonomyObject(taxO, maxNumSeqs, lpsnD)
    print(DONE)

    # download the genomes
    print(PRINT_2, end='', flush=True)
    __downloadGbffFromSpeciesList(speciesList, humanMapFN, gbffDir)
    print(DONE)

    # outgroup is always the last element of speciesList; return it
    return speciesList[-1]


def __selectSpeciesFromTaxonomyObject(taxO:Taxonomy, maxNumSeqs:int, \
                                                           lpsnD:dict) -> list:
    """ selectSpeciesFromTaxonomyObject:
            Accepts a Taxonomy object, an integer indicating the max number of
            ingroup sequences to download, and the LPSN dictionary as inputs.
            Uses methods from the Taxonomy class to create a list of species
            Taxonomy objects that possess assemblies available for downloading.
            Returns the list.
    """
    # constants
    ERR_PREFIX = "Could not pick an ingroup with exactly "
    ERR_SUFFIX = " sequences. Please try a larger number."
    
    # get the outgroup first (taxO may be modified in the process)
    outgroup:Taxonomy
    taxO, outgroup = taxO._pickOutgroup(lpsnD)

    # get an ingroup based on the object's composition (result can be random)
    speciesL = taxO._pickIngroup(maxNumSeqs)

    if speciesL == []:
        raise RuntimeError(ERR_PREFIX + str(maxNumSeqs) + ERR_SUFFIX)

    # append the outgroup to the end of the list
    speciesL.append(outgroup)

    return speciesL


def __downloadGbffFromSpeciesList(allTaxa:list, humanMapFN:str, gbffDir:str) \
                                                                       -> None:
    """ downloadGbffFromSpeciesList:
            Accepts a list of species Taxonomy objects, a string indicating the
            filename of the human map, and a string indicating the directory
            where assemblies should be downloaded as inputs. Downloads the seq-
            uences for each of the species and populates the human map file wi-
            th the necessary data. Does not return.
    """
    # delete any existing humanMap file
    if os.path.exists(humanMapFN):
        os.remove(humanMapFN)

    # download the ingroup sequences and make human map file
    filehandle = open(humanMapFN, "a")
    taxon:Taxonomy
    for taxon in allTaxa:
        # download gbff and get the filename
        gbffFN = __downloadGbff(taxon.assemblyFtp, gbffDir)

        # create the human map string
        humanMapStr = __makeHumanMapString(taxon, gbffFN)
        
        # add the string to the file
        filehandle.write(humanMapStr)
    
    # close the file
    filehandle.close()


def __downloadGbff(ftpPath:str, gbffDir:str) -> str:
    """ downloadGbff:
            Accepts two strings as inputs: one indicating the ftp path of the
            assembly and one indicating the directory where the downloaded file
            should be saved. Returns the filename of the newly downloaded file
            as a string.
    """
    # constants
    NCBI_FTP  = "ftp.ncbi.nlm.nih.gov"
    LEAD_STR  = "ftp://" + NCBI_FTP
    
    # make the folder if it doesn't yet exist
    if not os.path.exists(gbffDir):
        os.mkdir(gbffDir)

    # remove the leading string from the ftp path
    ftpPath = ftpPath[len(LEAD_STR):]

    # extract the filename from the ftp path
    gzFileName = os.path.basename(ftpPath)

    # add the folder info to the path
    gzFileName = os.path.join(gbffDir, gzFileName)

    # download the file
    downloadFileFromFTP(NCBI_FTP, ftpPath, gzFileName)

    # get the filename for the unpacked gbff
    gbffFileName = removeFileExtension(gzFileName)

    # unpack the gz
    decompressGZ(gzFileName)

    # return the gbff filename (without directory information)
    return os.path.basename(gbffFileName)


def __makeHumanMapString(speciesO:Taxonomy, filename:str) -> str:
    """ makeHumanMapString:
            Accepts a Taxonomy object (or a string) and a string indicating the
            filename of the genome for that object (or string) as inputs. Cons-
            tructs and returns the string for the human map file.
    """
    # if the speciesO is a string, then use that string as the species name
    if type(speciesO) is str:
        speciesName = speciesO

    # otherwise if the species name is invalid, then reflect this in the name
    elif ' (invalid name)' in speciesO.sciName:
        speciesName = speciesO.ncbiName + '__invalid'
    
    # otherwise, use the scientific name as the species name
    else:
        speciesName = speciesO.sciName

    # replace any spaces in the species name with underscores
    speciesName = re.sub(' ', '_', speciesName)

    # combine inputs to make map string and return it
    return filename + '\t' + speciesName + '\n'


