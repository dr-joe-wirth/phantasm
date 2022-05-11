# Author: Joseph S. Wirth

import sys, os, re, glob, copy
from PHANTASM.Parameter import Parameters
from PHANTASM.taxonomy.Taxonomy import Taxonomy
from PHANTASM.rRNA.runRnaBlast import __makeOutfmtString
from PHANTASM.utilities import parseCsv, ncbiIdsFromSearchTerm, ncbiSummaryFromIdList, ncbiELinkFromIdList, extractIdsFromELink, removeFileExtension
from PHANTASM.downloadGbff import __downloadGbffFromSpeciesList, _makeHumanMapString
from PHANTASM.coreGenes import _humanNameFromQueryGenbankFN
from Bio import Entrez, SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaIterator
from Bio.Blast.Applications import NcbiblastpCommandline
from param import XENOGI_DIR
sys.path.insert(0,os.path.join(sys.path[0],XENOGI_DIR))
import xenoGI.genomes, Bio.Entrez.Parser



###############################################################################
def phyloMarkerBlastRunner(geneNumsL:list, paramO:Parameters) -> None:
    """ phyloMarkerBlastRunner:
            Accepts a list of xenoGI gene numbers and a Parameters object as
            inputs. Creates an amino acid fasta file (faa) for the CDS in the
            query genomes corresponding to the provided gene number. Uses the
            fasta file to run a blastp against NCBI's nr database. If that bla-
            stp was not successful (often due to exceeding CPU limit), then it
            tries to run a blastp against NCBI's refseq_protein database. Does
            not return.
    """
    # constants
    PRINT_1 = 'Using the phylogenetic marker(s) to search for closely-related genomes ... '
    DONE = 'Done.'
    TEMP_FAA_FN = "temp.faa"
    TEMP_BLAST_FN = 'temp.blastp'
    FORMAT = 'fasta'
    DB_1 = "nr"
    DB_2 = "refseq_protein"

    # print status
    print(PRINT_1, end='', flush=True)
    
    # extract relevant parameters
    faaFN = paramO.phyloMarkerFaaFN
    blastFN = paramO.blastpResultFN
    blastExecutDirPath = paramO.blastExecutDirPath
    
    # get a SeqRecord object for the protein
    allSeqRecs = list()
    for geneNum in geneNumsL:
        seqRec = __getSeqRecordFromGeneNum(geneNum, paramO)
        allSeqRecs.append(seqRec)

    # save the file as a fasta
    SeqIO.write(allSeqRecs, faaFN, FORMAT)

    # determine the temporary filenames
    tempFaaFN = os.path.join(os.path.dirname(faaFN), TEMP_FAA_FN)
    tempBlastFN = os.path.join(os.path.dirname(blastFN), TEMP_BLAST_FN)

    # make sure the blast file is empty
    blastFH = open(blastFN, 'w')
    blastFH.close()

    

    # for each seq record
    #### doing it this way should ease the strain on remote blast computation
    #### this way is more likely to successfully get a result back
    for seqRec in allSeqRecs:
        # write the record to file
        SeqIO.write(seqRec, tempFaaFN, FORMAT)

        # run the blastp with the new fasta against nr


        print('against nr')


        __blastFaaAgainstDb(tempFaaFN, tempBlastFN, blastExecutDirPath, DB_1)

        # check that the blast was successful
        if os.path.getsize(tempBlastFN) == 0:


            print('against refseq_protein')


            # if not, then run a blast against refseq_protein
            __blastFaaAgainstDb(tempFaaFN, tempBlastFN, blastExecutDirPath, DB_2)
        
        # check that the blast was successful
        if os.path.getsize(tempBlastFN) == 0:
            # if not, then raise an error
            raise RuntimeError("blastp failed.")
        
        # open the temp file
        tempBlastFH = open(tempBlastFN, 'r')

        # move the data from the temp file into the blast file
        # open the blast file to begin saving results
        blastFH = open(blastFN, 'a')
        for line in tempBlastFH:
            blastFH.write(line)
        
        # close the blast file
        blastFH.close()
        
        # close the temp blast file
        tempBlastFH.close()

        # remove the temp files
        os.remove(tempFaaFN)
        os.remove(tempBlastFN)

        print('finished first blastp')

    print(DONE)


def __getSeqRecordFromGeneNum(geneNum:int, paramO:Parameters) -> SeqRecord:
    """ getSeqRecordFromGeneNum:
            Accepts an integer indicating a xenoGI gene number and a Parameters
            object as inputs. Retrieves and returns a SeqRecord object for the
            desired gene.
    """
    # constants
    GREP_FIND_1 = r'^\d+_([^\|]+\|[^\|]+\|[^-]+)-\S+$'
    GREP_FIND_2 = r"^\d+_([^-]+)-\S+$"
    GREP_REPL = r"\1"
    GENE_NAME_IDX = 0
    FORMAT = "fasta"

    # get the necessary data from paramO
    geneInfoFN = paramO.geneInfoFN
    fastaFilePath = paramO.fastaFilePath

    # load the geneInfo data for the given gene number
    genesO = xenoGI.genomes.genes(geneInfoFN)
    genesO.initializeGeneInfoD(geneInfoFN)
    geneInfo = genesO.numToGeneInfo(int(geneNum))

    # get the name for the provided gene
    geneName = geneInfo[GENE_NAME_IDX]

    # extract the species name from the gene name
    speciesName = re.sub(GREP_FIND_1, GREP_REPL, geneName)

    # if the grep didn't work
    if geneName == speciesName:
        speciesName = re.sub(GREP_FIND_2, GREP_REPL, geneName)
    
    # use the species name to determine which file to load
    allFastasL = glob.glob(fastaFilePath)

    # navigate to the correct fasta file
    for fastaFN in allFastasL:
        if speciesName in fastaFN:
            break
    
    # read the fasta file
    parsed:FastaIterator = SeqIO.parse(fastaFN, FORMAT)

    # identify and return the desired record
    record:SeqRecord
    for record in parsed.records:
        if record.name == geneName:
            return record
    
    # if record was not found, then raise an exception
    raise Exception("invalid gene number provided")


def _locusTagToGeneNum(locusTag:str, geneInfoFN:str) -> int:
    """ geneNumToLocusTag:
            Accepts an string indicating a locus tag and a string indicating
            the gene info filename as inputs. Retrieves and returns a xenoGI
            gene number for the provided locus tag as an integer.
    """
    # constants
    LOCUS_TAG_IDX = 2

    # load the genesO object
    genesO = xenoGI.genomes.genes(geneInfoFN)
    genesO.initializeGeneInfoD(geneInfoFN)

    # get the gene number for the provided gene and return it
    for geneNum in genesO.geneInfoD.keys():
        if genesO.geneInfoD[geneNum][LOCUS_TAG_IDX] == locusTag:
            return geneNum
    
    # if a number was not found then raise an exception
    raise Exception("Invalid locus tag provided")


def __blastFaaAgainstDb(faaFN:str, outFN:str, blastExecutDirPath:str, \
                                                         database:str) -> None:
    """ blastFaaAgainstDb:
            Accepts a string indicating the amino acid fasta file to be used as
            the query, a string indicating where the save the blastp results,
            a string indicating the folder containing the blast+ executables,
            and a string indicating the database to query against. Runs blastp
            against the provided database remotely. Does not return.
    """
    # constants
    REMOTE = True
    HEADERS = ['sseqid', 'stitle', 'staxid', 'sacc', 'pident', 'length',
               'mismatch', 'gaps', 'gapopen', 'evalue', 'bitscore']
    
    # make the outfmt string
    outfmt = __makeOutfmtString('6', HEADERS)

    # construct the blastp command
    blastpExe = os.path.join(blastExecutDirPath, 'blastp')
    blastpCmd = NcbiblastpCommandline(cmd=blastpExe, query=faaFN, db=database,
                                       remote=REMOTE, outfmt=outfmt, out=outFN)
    
    # execute the blastp command
    blastpCmd()
###############################################################################


###############################################################################
def getRelatives(oldParamO:Parameters, newParamO:Parameters, lpsnD:dict, \
                                                       maxNumSeqs:int) -> list:
    """ getRelatives:
            Accepts the old and new Parameters objects and the LPSN dictionary
            as inputs. Uses the blastp results for the phylogenetic marker to
            determine which species should be used for downstream phylogenomic
            analyses. Returns a list of the selected species.
    """

    # extract relevant data from parameter dictionaries
    oldTaxFN = glob.glob(oldParamO.taxonomyObjectFilePath).pop()
    blastFN = oldParamO.blastpResultFN
    newTaxDir = os.path.dirname(newParamO.taxonomyObjectFilePath)
    newTaxExt = os.path.splitext(newParamO.taxonomyObjectFilePath)[1]

    # set entrez email
    Entrez.email = oldParamO.email

    # load the taxonomy object from file
    taxO = Taxonomy.load(oldTaxFN)

    # extract the assembly and blastp data from the blast file
    assemblyData:tuple = __findMissingAssemblies(blastFN, taxO, lpsnD)

    # complete assembly and blastp data keyed by assembly id
    assembliesD:dict = assemblyData[0]

    # key = taxid of objects in taxO w/o assemblies; val = list of assembly ids
    existNoAssD:dict = assemblyData[1]

    # key = valid species name absent from taxO; val = list of assembly ids
    missingD:dict = assemblyData[2]
    
    # add any missing taxa to the existing Taxonomy object
    # once added, make an entry in existNoAssD
    for speciesName in missingD.keys():
        # get the taxonomy id for the missing species
        taxid = assembliesD[missingD[speciesName][0]]['taxid']

        # create a new object for the missing species (internal by default)
        newTaxon = Taxonomy(taxid, speciesName, 'species')

        # try to add the newTaxon to the existing object
        try:
            taxO._importExistingSubTax(newTaxon, lpsnD)
        
        # if import fails, then skip the newTaxon
        # should really only happen if the newTaxon is in a different domain
        except: continue

        # make sure taxO is still referencing the root!
        taxO = taxO.getRoot()

        # make a corresponding entry in existNoAssD
        existNoAssD[taxid] = missingD[speciesName]

    # update assembly info for species w/o assemblies but present in the blastp
    for taxId in existNoAssD.keys():
        # get a list of all the assembly ids associated with that taxon
        assmIdsL:list = existNoAssD[taxId]

        # get a handle to the species of interest
        speciesO = taxO.getRoot().getDescendantByTaxId(taxId)

        # add the assembly to taxO; automatically saves the best assembly
        for assmId in assmIdsL:
            # update object with the assembly
            assembly:dict = assembliesD[assmId]
            speciesO._updateAssemblyInfo(assembly['summary'])

    # save a taxonomy file with the new taxonomy
    taxO = taxO.getRoot()
    newTaxFN = os.path.join(newTaxDir, taxO.sciName + newTaxExt)
    taxO.save(newTaxFN)

    # get the unique genera with their min and max bitscores
    uniqueGeneraD = __getUniqueGeneraFromAssembliesD(assembliesD, lpsnD)

    # sort the unique genera from highest max bitscore to lowest max bitscore
    sortMethod = lambda genName: uniqueGeneraD[genName]['max']
    generaSortedL = sorted(uniqueGeneraD.keys(), key=sortMethod, reverse=True)

    # get the minimum score for the genus with the highest maximum score
    closestGenusMinScore = uniqueGeneraD[generaSortedL[0]]['min']
    
    # track the index while iterating through genera from best hit to worst
    relatives = list()
    for idx in range(len(generaSortedL)):
        # use the index to extract the genus name from the list
        genusName = generaSortedL[idx]

        # if the max score for the genus >= the lowest score for the best genus
        if uniqueGeneraD[genusName]['max'] >= closestGenusMinScore:
            # get a handle to that genus
            genus = taxO.getDescendantBySciName(genusName)

            # add all of the species with assemblies to the relatives list
            relatives.extend(genus.getAllSpeciesWithAssemblies(list))
        
        # otherwise, the list is sorted; it is safe to stop looping
        else: break

    # it is possible that too many relatives were selected in the previous step
    # if this is the case, then we need to select a subset of them
    if len(relatives) > maxNumSeqs:
        # copy relatives and then empty the list
        speciesToPickL = copy.deepcopy(relatives)
        relatives = list()

        # reverse the indices for on-the-fly popping
        speIndicesL = list(range(len(speciesToPickL)))
        speIndicesL.reverse()

        # go through the indices (in reverse order)
        for speIdx in speIndicesL:
            # extract the current relative
            relative:Taxonomy = speciesToPickL[speIdx]

            # if the current relative is type material ...
            if relative.parent.typeMaterial == relative:
                # ... then pop it from the old and add to relatives
                relatives.append(speciesToPickL.pop(speIdx))
        
        # next, begin adding all species from the related genera
        # speciesToPickL is already in order from closest genus to furthest
        # speciesToPickL will not have any overlap in relatives
        while len(relatives) < maxNumSeqs and len(speciesToPickL) > 0:
            # this will functionally loop through the remaining species
            relatives.append(speciesToPickL.pop(0))
        
    # the current index can be used to remove the already processed genera
    generaSortedL = generaSortedL[idx:]
    
    # for the remaining genera in the list, add one species from the genus
    # loop until the max num seqs is reached or the genus list is depleted
    while len(relatives) < maxNumSeqs and generaSortedL != []:
        # get the current genus name and a handle to the corresponding object
        curGenusName = generaSortedL.pop(0)
        curGenus = taxO.getDescendantBySciName(curGenusName)
        
        # only proceed if a genus was found (curGenus=False on a failed search)
        if curGenus:
            # get a list of all possible candidate species
            allCandidates = curGenus.getAllIngroupCandidateSpecies()

            # sort the candidate species by assembly coverage
            sortMethod = lambda speciesO: speciesO.assCoverage
            allCandidates = sorted(allCandidates, key=sortMethod, reverse=True)

            # if the type material is a candidate, then this is the selection
            if curGenus.typeMaterial in allCandidates:
                relatives.append(curGenus.typeMaterial)
            
            # if type unavailable, select candidates by the following criteria:
                # bin1:  complete and representative
                # bin2:  complete and not representative
                # bin3:  incomplete and representative
                # bin4:  incomplete and not representative
            else:
                # for each candidate species add to the correct bin
                cand:Taxonomy
                bin1 = list()
                bin2 = list()
                bin3 = list()
                bin4 = list()
                for cand in allCandidates:
                    # complete ...
                    if cand.assIsComplete:
                        # ... and representative
                        if cand.refSeqCategory == 'representative genome':
                            bin1.append(cand)
                        # ... and not representative
                        else:
                            bin2.append(cand)
                    # incomplete ...
                    else:
                        # ... and representative
                        if cand.refSeqCategory == 'representative genome':
                            bin3.append(cand)
                        # ... and not representative
                        else:
                            bin4.append(cand)
                
                # genomes are sorted by assembly coverage in each bin
                # prioritize bins first, then prioritize assembly coverage
                # append only one species per genus
                if bin1 != []:
                    relatives.append(bin1[0])
                elif bin2 != []:
                    relatives.append(bin2[0])
                elif bin3 != []:
                    relatives.append(bin3[0])
                else:
                    relatives.append(bin4[0])
    
    return relatives


def __findMissingAssemblies(blastFN:str, taxO:Taxonomy, lpsnD:dict) -> tuple:
    """ findMissingAssemblies:
            Accepts a string indicating a blastp result file in outfmt 6 with
            custom headers (see parseBlastpFile for more info), a Taxonomy obj-
            ect, and the LPSN dictionary as inputs. Creates and returns three 
            dictionaries: the first is the dictionary returned by linkAssembli-
            esWithBlastpResults. The second is a dictionary keyed by the taxids
            of objects already present in the Taxonomy object that currently
            lack assemblies. The third is a dictionary keyed by the names of 
            validly published species present in the blast result but absent
            from the Taxonomy object. Both the second and third dictionaries
            have the corresponding assembly uid as their values.
    """
    # get a dictionary with both the assembly and blastp data
    assembliesD = __linkAssembliesWithBlastpResults(blastFN)

    # populate additional diciontaries with hits that need to be resolved
    existNoAssD = dict()
    missingD = dict()
    for assId in assembliesD.keys():
        # extract values from the assembly dictionary
        taxid = assembliesD[assId]['taxid']
        name = assembliesD[assId]['name']
    
        # make sure name is not a synonym
        sciName = __renameSpeciesByLpsn(name, lpsnD)
        
        # attempt to lookup an existing object by taxid
        desc = taxO.getDescendantByTaxId(taxid)

        # if lookup by taxid failed, then try lookup by sciName
        if not desc:
            desc = taxO.getDescendantBySciName(sciName)
        
        # if lookup by taxid and sciName failed, then try lookup by ncbiName
        if not desc:
            desc = taxO.getDescendantByNcbiName(name)
        
        # if an existing descendant was found 
        if desc:
            # make sure the descendant is marked as internal
            desc.isExternal = False

            # if the existing descendant lacks an assembly
            if desc.assemblyFtp is None:
                # make a list of assembly ids if the taxid isn't a key
                if desc.taxid not in existNoAssD.keys():
                    existNoAssD[desc.taxid] = [assId]
                
                # update the list of assembly ids if the taxid is already a key
                else:
                    existNoAssD[desc.taxid].append(assId)
        
        # if no existing descendant was found in taxO, but it is a valid name
        elif sciName is not None:
            # make a list of assembly ids if the sciName isn't a key
            if sciName not in missingD.keys():
                missingD[sciName] = [assId]
            
            # update the list of assembly ids if the sciName is already a key
            else:
                missingD[sciName].append(assId)
    
    return assembliesD, existNoAssD, missingD


def __getUniqueGeneraFromAssembliesD(assembliesD:dict, lpsnD:dict) -> dict:
    """ getUniqueGeneraFromAssembliesD:
            Accepts the assembly and LPSN dictionaries as inputs. Extracts the
            unique genera for all the blastp hits to species known by the LPSN
            as well as each genus's min and max bitscores. Creates a dictionary
            of dictionaries. The nested dictionary is keyed by 'min' and 'max'
            and contains the corresponding bitscores. The parental dictionary 
            is keyed by genus name with the corresponding min/max dictionary as
            its value. Returns the newly created dictionary.
    """
    # constants
    GREP_FIND = r'^(\S+) .+$'
    GREP_REPL = r'\1'
    MIN = 'min'
    MAX = 'max'

    # for each assembly id
    outD = dict()
    for assmId in assembliesD.keys():
        # get the entry from the dictionary
        entry:dict = assembliesD[assmId]

        # extract the name and the blast bitscore from the entry
        fullName:str  = entry['name']
        score = float(entry['blast']['bitscore'])

        # rename the hit with its LPSN name
        fullName = __renameSpeciesByLpsn(fullName, lpsnD)

        # fullName will be None if the input was absent from LPSN
        if fullName is not None:
            # extract the genus name from the binomial name
            genusName = re.sub(GREP_FIND, GREP_REPL, fullName)

            # make a new entry for a genus that has not yet been encountered
            if genusName not in outD.keys():
                # save the current score as both the min and the max values
                entry = {MIN: score,
                         MAX: score}
                outD[genusName] = entry
            
            # update an entry if the genus name has already been encountered
            else:
                # update the min score if it is greater than the current score
                if outD[genusName][MIN] > score:
                    outD[genusName][MIN] = score
                
                # update the max score if it is less than the current score
                if outD[genusName][MAX] < score:
                    outD[genusName][MAX] = score

    return outD


def __linkAssembliesWithBlastpResults(blastFN:str) -> dict:
    """ linkAssembliesWithBlastpResults:
            Accepts a string indicating a blastp result file in outfmt 6 with
            custom headers (see parseBlastpFile for more info) as input. Finds
            and links the assembly information with the blast result via a dic-
            tionary. Returns the dictionary.
    """
    # constants
    PROT_DB = 'protein'
    NUCL_DB = 'nuccore'
    ASSM_DB = 'assembly'

    # parse the blastp file into a dictionary of hits and search string
    parsedBlast:tuple = __parseBlastpFile(blastFN)
    hitsD:dict = parsedBlast[0]
    searchStr:str = parsedBlast[1]

    # use the search string to get the protein ids
    protIds = ncbiIdsFromSearchTerm(searchStr, PROT_DB)

    # get a summary of the protein ids (used to link assemblies later)
    protSummaries = ncbiSummaryFromIdList(protIds, PROT_DB)

    # make a dictionary linking the protein ids to the blast hits
    prot2hitsD = dict()
    for protSum in protSummaries:
        # extract the import values from the summary
        accession = protSum['Caption']
        uid = protSum['Id']

        # handle the rare instance where 'accession' is not in 'hitsD'
        if accession not in hitsD.keys():
            continue

        # update hitD so it includes the protein id
        hitsD[accession]['protein id'] = uid

        # add the hit to the dictionary
        prot2hitsD[uid] = hitsD[accession]

    # use the protein ids to link to nuccore
    try:
        prot2nuclLink = ncbiELinkFromIdList(protIds, PROT_DB, NUCL_DB)
    except:
        prot2nuclLink = list()
        for idx in range(len(protIds)):
            try:
                proId = protIds[idx]
                prot2nuclLink += ncbiELinkFromIdList([proId], PROT_DB, NUCL_DB)
            except:
                continue
    
    # extract the nuccore ids from the link
    nucl2protD = extractIdsFromELink(prot2nuclLink)
    nuclIds = list(nucl2protD.keys())

    # use the nuccore ids to get link to assembly
    nucl2assmLink = ncbiELinkFromIdList(nuclIds, NUCL_DB, ASSM_DB)

    # get a dictionary keyed by assembly ids whos values are nuccore ids
    assm2nuclD = extractIdsFromELink(nucl2assmLink)

    # refine the assemblies to those with annotations and full genome rep
    __refineAssemblySelection(assm2nuclD)

    # get the assembly summaries for the ids
    assmSummaries = ncbiSummaryFromIdList(list(assm2nuclD.keys()), 'assembly')
    assmSummaries = assmSummaries['DocumentSummarySet']['DocumentSummary']

    # construct the final dictionary
    assembliesD = dict()
    assmSum:Bio.Entrez.Parser.DictionaryElement
    for assmSum in assmSummaries:
        # get the assembly id
        assmId = assmSum.attributes['uid']

        # use the assembly id to find protein id via the nuccore id
        protId = nucl2protD[assm2nuclD[assmId]]

        # use the protein id to look up the blast hit
        blast = prot2hitsD[protId]

        # construct and save the entry to the assemblies dictionary
        assembliesD[assmId] = {'name': assmSum['SpeciesName'],
                               'taxid': assmSum['Taxid'],
                               'blast': blast,
                               'summary': assmSum}
    
    return assembliesD


def __renameSpeciesByLpsn(ncbiName:str, lpsnD:dict) -> str:
    """ renameSpeciesByLpsn:
            Accepts a species name from NCBI as a string and the LPSN dictiona-
            ry as inputs. Renames the string to coincide with the name saved in
            the LPSN. Returns the name. If the provided name is absent from the
            LPSN, then returns None.
    """
    # extract the species and species synonym dictionaries from lpsnD
    speciesSynD:dict = lpsnD['species_syn']
    speciesD:dict = lpsnD['species']

    # if the ncbiName is a synonym, then rename it
    if ncbiName in speciesSynD.keys():
        sciName = speciesSynD[ncbiName]

    # if the ncbiName is the correct name, then no modifications necessary
    elif ncbiName in speciesD.keys():
        sciName = ncbiName
    
    # if the ncbiName is absent from the LPSN, then renaming is not possible
    else:
        sciName = None
    
    return sciName


def __parseBlastpFile(blastFN:str) -> tuple:
    """ parseBlastpFile:
            Accepts a string indicating a blastp result file in outfmt 6 with
            headers as specified below. Parses the file into a dictionary and
            creates a string that can be used to search for the corresponding
            protein ids. Returns a tuple of the dictionary and the search str-
            ing.
    """
    # expected blastp headers (and their order)
    SSEQID = 0
    STITLE = 1
    STAXID = 2
    SACC = 3
    PIDENT = 4
    LENGTH = 5
    MISMATCH = 6
    GAPS = 7
    GAPOPEN = 8
    EVALUE = 9
    BITSCORE = 10

    # constants
    BLAST_DELILM = '\t'
    SEARCH_SUFFIX = '[ACCN]'
    SEARCH_SEP = ' OR '

    # read the file into memory as a list of rows
    parsed = parseCsv(blastFN, BLAST_DELILM)

    # make the search string and parse the blast hits into a dictionary
    hitsD = dict()
    searchStr = str()
    for row in parsed:
        # extract the protein accession number
        protAccn = row[SACC]

        # save the blast hit in dictionary format
        hitsD[protAccn] = {'sseqid': row[SSEQID],
                           'stitle': row[STITLE],
                           'staxid': row[STAXID],
                           'sacc': protAccn,
                           'pident': row[PIDENT],
                           'length': row[LENGTH],
                           'mismatch': row[MISMATCH],
                           'gaps': row[GAPS],
                           'gapopen': row[GAPOPEN],
                           'evalue': row[EVALUE],
                           'bitscore': row[BITSCORE]}

        # construct the search string for getting the uids for the protein db
        searchStr += protAccn + SEARCH_SUFFIX + SEARCH_SEP
    
    # remove the trailing sep character(s) from the search string
    searchStr = searchStr[:-len(SEARCH_SEP)]

    # return the blast hit dictionary and the search string
    return hitsD, searchStr


def __refineAssemblySelection(assmD:dict) -> None:
    """ refineAssemblySelection:
            Accepts the assembly dictionary as input. Removes any assemblies
            from the dictionary that without annotations or full genome repres-
            entation. Does not return.
    """
    # constants
    FILTER_1 = '"has annotation"[Properties]'
    FILTER_2 = '"full genome representation"[Properties]'
    SUFFIX = '[UID]'
    OR = ' OR '
    AND = ' AND '
    DATABASE = 'assembly'

    # add the filters to the search string
    searchStr = FILTER_1 + AND + FILTER_2

    # add each assembly id to the search string
    searchStr += AND + ' ('
    for assmId in assmD.keys():
        searchStr += assmId + SUFFIX + OR
    
    # remove the trailing OR character(s) and close the parentheses
    searchStr = searchStr[:-len(OR)]
    searchStr += ')'

    # get a list of assembly ids that match the search string
    refinedIds = ncbiIdsFromSearchTerm(searchStr, DATABASE)

    # make a set of all current assembly ids
    allAssmIds = set(list(assmD.keys()))

    # if the set of refined assembly ids differs from the current set of ids
    if set(refinedIds) != allAssmIds:
        # get a set of the current ids that are not also refined ids
        badIds = allAssmIds.difference(set(refinedIds))

        # refine the dictionary by removing the "unrefined" ids from it
        for id in badIds:
            assmD.pop(id)
###############################################################################


###############################################################################
def xenogiInterfacer_2(allQryGbksL:list, oldParamO:Parameters, newParamO:Parameters, \
                                                       lpsnD:dict) -> Taxonomy:
    """ xenogiInterfacer_2:
            Accepts a list of strings indicating the filenames of the query ge-
            nbanks, a Parameters object for the first analysis, the Parameters
            object for the final analysis, and the LPSN dictionary as inputs.
            Uses the phylogenetic marker's blastp results to determine which
            genomes to download. Downloads (or makes symlinks for) these genom-
            es and makes the human map file in the process. Returns the outgro-
            up as a Taxonomy object.
    """
    # constants
    PRINT_1 = 'Finding missing relatives and updating the Taxonomy object ... '
    PRINT_2 = 'Downloading genbank files from NCBI ... '
    DONE = 'Done.'

    # extract relevant data from the Parameters objects
    maxNumSeqs = oldParamO.maxNumTreeLeaves - len(allQryGbksL) - 1 # 1 outgroup
    taxObjectFilePath = newParamO.taxonomyObjectFilePath
    oldGenbankFilePath = oldParamO.genbankFilePath
    newGenbankFilePath = newParamO.genbankFilePath
    newWorkDir = newParamO.workdir
    newHumanMapFN = newParamO.fileNameMapFN

    # make the new working directory
    if not os.path.exists(newWorkDir):
        os.mkdir(newWorkDir)

    # get the new relatives based on phylogenetic marker blastp
    print(PRINT_1, end='', flush=True)
    speciesL = getRelatives(oldParamO, newParamO, lpsnD, maxNumSeqs)
    print(DONE)

    # load the taxonomy object created by getRelatives
    taxFN = glob.glob(taxObjectFilePath).pop()
    taxO = Taxonomy.load(taxFN)

    # add the outgroup to the list
    taxO:Taxonomy
    outgroup:Taxonomy
    taxO, outgroup = taxO._pickOutgroup(lpsnD)
    speciesL.append(outgroup)

    # save the modified taxonomy object and remove previously loaded version
    dir = os.path.dirname(taxObjectFilePath)
    ext = os.path.splitext(taxObjectFilePath)[1]
    os.remove(taxFN)
    taxO.save(os.path.join(dir, taxO.getRoot().sciName + ext))

    # initialize variables to populate while looping
    fileToSpeD = dict()
    spe:Taxonomy

    # go through the list back-to-front to enable on-the-fly popping
    idxL = list(range(len(speciesL)))
    idxL.reverse()
    for idx in idxL:
        # get the species object from the list
        spe = speciesL.pop(idx)

        # extract the filename from the ftp path
        file = os.path.basename(spe.assemblyFtp)
        file = removeFileExtension(file)

        # add the entry to the dictionary
        fileToSpeD[file] = spe
    
    # get a set of the genomes that have already been downloaded
    filesDownloaded = glob.glob(oldGenbankFilePath)
    for idx in range(len(filesDownloaded)):
        filesDownloaded[idx] = os.path.basename(filesDownloaded[idx])
    filesDownloaded = set(filesDownloaded)

    # determine which files still need to be downloaded and which exist
    filesRequired = set(fileToSpeD.keys())
    filesToDownload = filesRequired.difference(filesDownloaded)
    filesToCopy = filesRequired.intersection(filesDownloaded)

    # make the new genbank directory
    newGbkDir = os.path.dirname(newGenbankFilePath)
    os.mkdir(newGbkDir)

    # make sure the genbank directories are absolute paths
    newGbkDir = os.path.abspath(newGbkDir)
    oldGbkDir = os.path.dirname(oldGenbankFilePath)
    oldGbkDir = os.path.abspath(oldGbkDir)

    # initialize the human map string
    #### cannot open file and append
    #### __downloadGbffFromSpeciesList will delete the existing human map file
    humanMapStr = ""

    # make symlinks for files and simultaneously make the human map string
    for filename in filesToCopy:
        # extract the Taxonomy object from the dictionary
        taxon:Taxonomy = fileToSpeD[filename]

        # add the human map string data to the growing string
        humanMapStr += _makeHumanMapString(taxon, filename)

        # make a symlink for the existing file
        oldFN = os.path.join(oldGbkDir, filename)
        newFN = os.path.join(newGbkDir, filename)
        os.symlink(oldFN, newFN)
    
    # for each of the user's query genomes
    for queryGbff in allQryGbksL:
        # get the file name
        basename = os.path.basename(queryGbff)

        # make a symlink to the user's input file
        oldFN = os.path.abspath(queryGbff)
        newFN = os.path.join(newGbkDir, basename)
        os.symlink(oldFN, newFN)

        # get the human name from the query file name
        humanName = _humanNameFromQueryGenbankFN(queryGbff)

        # make the human map string and append it to the human map file
        humanMapStr += _makeHumanMapString(humanName, basename)
    
    # make a list of species objects that still need to be downloaded
    speciesL = list()
    for file in filesToDownload:
        speciesL.append(fileToSpeD[file])
    
    # download the species in the list
    print(PRINT_2, end='', flush=True)
    __downloadGbffFromSpeciesList(speciesL, newHumanMapFN, newGbkDir)
    print(DONE)

    # append the human map string to the file
    filehandle = open(newHumanMapFN, 'a')
    filehandle.write(humanMapStr)
    filehandle.close()

    return outgroup
