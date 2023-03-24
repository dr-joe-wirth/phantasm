# Author: Joseph S. Wirth

from PHANTASM.utilities import downloadFileFromFTP, extractContentsFromGZ, changeFileExtension
from pathlib import Path
from io import TextIOWrapper
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.SeqIO.InsdcIO import GenBankIterator
import logging, os


def rnaBlastRunner(allQryGbksL:list, workingDir:str, blastExecutDirPath:str) \
                                                                        -> str:
    """ rnaBlastRunner:
            Accepts a list of all the query genbank files, the path to a worki-
            ng directory, and the directory containing the blast+ executables
            as inputs. Downloads and constructs a blast database from NCBI for
            the 16S rRNA gene sequences of Bacteria and Archaea. Uses blastn to
            compare the query 16S rRNA gene sequences against those in the dat-
            abase. Returns the filename containing the blastn results (str).
    """
    # constants
    RNA_FOLDER = "16S/"
    BLAST_FOLDER = "blastdb/"
    
    # make folders for analyses
    rnaDir = os.path.join(workingDir, RNA_FOLDER)
    blastDir = os.path.join(rnaDir, BLAST_FOLDER)
    Path(rnaDir).mkdir(parents=True, exist_ok=True)
    Path(blastDir).mkdir(parents=True, exist_ok=True)

    # Construct query fasta and blast database
    queryFna, blastdb, blastFileOut = __blastPrep(allQryGbksL, rnaDir, blastDir, blastExecutDirPath)

    # Run the blast
    __runRnaBlast(queryFna, blastdb, blastFileOut, blastExecutDirPath)

    # return the filename for accessing the results
    return blastFileOut


def __blastPrep(allQryGbksL:list, rnaDir:str, blastDir:str, \
                                              blastExecutDirPath:str) -> tuple:
    """ blastPrep:
            Accepts a list of query genbank files, and three directories as in-
            puts. Extracts all sequences from the genbanks and writes them to a
            file in fasta format. Downloads the curated 16S rRNA gene sequences
            from NCBI for Bacteria and Archaea and constructs a blast database
            made for those sequences. Returns a three-ple composed of the file-
            name of the newly created 16S rRNA fasta (query), the filename of
            the blast database, and the filename where the blast results will
            be written to.
    """
    # constants
    QUERY_FNA_FILENAME = "query.16SrRNA.fna"
    SUBJECT_GBFF = "prok.16SrRNA.gbff"
    RESULT_EXTENSION = '.blastn'
    EOL = " ... "
    PRNT_1 = "Constructing the the query 16S rRNA fasta file" + EOL
    PRNT_2 = "Downloading subject sequences from NCBI" + EOL
    PRNT_3 = "Constructing the subject 16S rRNA fasta file" + EOL
    PRNT_4 = "Making blast database from '"
    DONE = "Done.\n"
    
    # create logger
    logger = logging.getLogger(__name__ + "." + __blastPrep.__name__)
    
    # add path information to queryFnaFile
    queryFnaFile = os.path.join(rnaDir, QUERY_FNA_FILENAME)

    # add path information to SUBJECT_GBFF
    subjectGbff = os.path.join(rnaDir, SUBJECT_GBFF)

    # make query fasta from the input gbffs
    print(PRNT_1, end="", flush=True)
    logger.info(PRNT_1)
    
    # open the file as a filehandle
    queryFnaFH = open(queryFnaFile, "w")

    # get 16s rRNA sequences from each query gbff and save in a single fasta
    for queryGenbank in allQryGbksL:
        __rnaFastaFromGbk(queryGenbank, queryFnaFH)
    
    # close the fasta file
    queryFnaFH.close()
    print(DONE)
    logger.info(DONE)

    # download the sequence files from NCBI and save the gbff
    print(PRNT_2)
    logger.info(PRNT_2)
    __downloadRnaSequences(subjectGbff, rnaDir)
    print(DONE)
    logger.info(DONE)

    # construct an fna with all of NCBI's curated 16S rRNA gene sequences
    print(PRNT_3, end='', flush=True)
    logger.info(PRNT_3)
    subjectFna = __makeFastaForBlastDB(subjectGbff)
    print(DONE)
    logger.info(DONE)

    # get the blastn database filename
    dbFilename = os.path.splitext(subjectFna)[0]
    dbFilename = os.path.basename(dbFilename)
    dbFilename = os.path.join(blastDir, dbFilename)

    # create the blastn database from the subject fna
    print(PRNT_4 + subjectFna + EOL, end="", flush=True)
    logger.info(PRNT_4 + subjectFna)
    __makeblastdbFromFna(subjectFna, dbFilename, blastExecutDirPath)
    print(DONE)
    logger.info(DONE)

    # generate the filename for saving the blast result
    blastResultFilename = changeFileExtension(queryFnaFile, RESULT_EXTENSION)

    # return
    return queryFnaFile, dbFilename, blastResultFilename


def __rnaFastaFromGbk(filename:str, outFH:TextIOWrapper) -> None:
    """ rnaFastaFromGbk:
            Accepts the filename of a query genbank and a filhandle for saving
            the fasta file. Parses the genbank file and retrieves the 16S rRNA 
            gene sequences using the annotations. Writes the 16S rRNA gene seq-
            uences in fasta format. Does not return.
    """
    # constants
    ERR_MSG = "Could not extract 16S rRNA gene sequences from the " + \
                                                       "provided genbank file."
    
    logger = logging.getLogger(__name__ + "." + __rnaFastaFromGbk.__name__)
    
    # get a list of 16S rRNA SeqRecord objects
    recs = __rnaSeqGrabber(filename)

    if len(recs) == 0:
        logger.error(ERR_MSG)
        raise BaseException(ERR_MSG)

    # write the records to file
    SeqIO.write(recs, outFH, 'fasta')


def __rnaSeqGrabber(filename: str) -> list:
    """ rnaSeqGrabber:
            Accepts the filename of a whole genome sequence in genbank file 
            format as input. Parses the records and features from the file and 
            returns only those that contain either the string "16S" or the 
            string "16s".
    """
    # file format needs to be genbank
    FILE_FORMAT = "genbank"

    # read genome
    genome_parsed = SeqIO.parse(filename, FILE_FORMAT)

    # initalize variables
    rnaRecords = list()

    # for each record (loop allows multi-genbanks)
    record:SeqRecord
    for record in genome_parsed:
        # find and store rRNAs that also have '16S' in their product field
        feature:SeqFeature
        for feature in record.features:
            # make sure is16S is false each time through the loop
            is16S = False
            
            if feature.type == "rRNA":
                # check the product list for a string containing '16S'
                products = feature.qualifiers["product"]

                # for each string in the product field, check if 16S is present
                for prod in products:
                    if "16S" in prod or '16s' in prod:
                        is16S = True
                        break
                
                # save the 16S feature and its record as a tuple in the list
                if is16S:
                    # make a new SeqRecord for the 16S rRNA
                    rnaRec = SeqRecord(feature.extract(record.seq))
                    rnaRec.id = feature.qualifiers['locus_tag'][0]
                    rnaRec.description = feature.qualifiers['product'][0]

                    # add the record to the list
                    rnaRecords.append(rnaRec)
    
    # return the list of tuples
    return rnaRecords


def __downloadRnaSequences(outFile:str, outDir:str) -> None:
    """ downloadRnaSequences:
            Accepts a file name and the directory that contains it as input.
            Downloads the 16S rRNA gene sequences from NCBI's curated 'targeted
            loci' project. Combines the two sequence files into a single file 
            and writes this file. Does not return.
    """
    # constants
    NCBI_FTP = "ftp.ncbi.nlm.nih.gov"
    BACTERIA_FTP_PATH = "/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.gbff.gz"
    ARCHAEA_FTP_PATH = "/refseq/TargetedLoci/Archaea/archaea.16SrRNA.gbff.gz"
    GAP = " " * 4
    SUFFIX = "'."
    PRNT_1 = GAP + "Saving bacterial 16S rRNA sequences to '"
    PRNT_2 = GAP + "Saving archaeal  16S rRNA sequences to '"
    PRNT_3 = GAP + "Writing subject gbff to '"
    
    logger = logging.getLogger(__name__ + "." + __downloadRnaSequences.__name__)

    # get the filenames so that they can be written with their original names
    bacteriaFile = os.path.basename(BACTERIA_FTP_PATH)
    archaeaFile  = os.path.basename(ARCHAEA_FTP_PATH)

    # add the path information to the filenames
    bacteriaFile = os.path.join(outDir, bacteriaFile)
    archaeaFile  = os.path.join(outDir, archaeaFile)

    # download the individual files
    print(PRNT_1 + bacteriaFile + SUFFIX)
    logger.info(PRNT_1 + bacteriaFile + SUFFIX)
    downloadFileFromFTP(NCBI_FTP, BACTERIA_FTP_PATH, bacteriaFile)
    print(PRNT_2 + archaeaFile + SUFFIX)
    logger.info(PRNT_2 + archaeaFile + SUFFIX)
    downloadFileFromFTP(NCBI_FTP,  ARCHAEA_FTP_PATH, archaeaFile)

    # extract the contents of the files
    bacteriaFileContents = extractContentsFromGZ(bacteriaFile)
    archaeaFileContents  = extractContentsFromGZ(archaeaFile)

    # combine the two files
    print(PRNT_3 + outFile + SUFFIX)
    logger.info(PRNT_3 + outFile + SUFFIX)
    allContents = bacteriaFileContents + archaeaFileContents

    # write the file
    fH = open(outFile, 'w')
    fH.write(allContents)
    fH.close()


def __makeFastaForBlastDB(genbankFile:str) -> str:
    """ makeFastaForBlastDB
            Takes a genbank file containing 16S rRNA gene sequences of both 
            bacteria and archaea. Converts it to amino acid fasta format. Retu-
            rns the name of the newly written fasta file.
    """
    # constants
    OUT_FILE_EXTENSION = "fna"

    # construct output file name
    outFile = changeFileExtension(genbankFile, OUT_FILE_EXTENSION)

    # make the fasta file
    __makeFastaForBlastDBHelper(genbankFile, outFile)

    return outFile


def __makeFastaForBlastDBHelper(genbankFN:str, fastaFN:str) -> None:
    """ makeFastaHelper:
            Accepts a string indicating a genbank filename and a string indica-
            ting where to save the resulting fasta file as inputs. Converts the
            genbank file into a fasta of 16S rRNA gene sequences. Does not ret-
            urn.
    """
    # constants
    SEP = "|"

    # get an iterator for each genbank record
    parsed:GenBankIterator = SeqIO.parse(genbankFN, 'genbank')

    # initialize a list of records to write to file
    recordsL = list()

    # for each record in the genbank file
    rec:SeqRecord
    feature:SeqFeature
    for rec in parsed.records:

        # get the scientific name and its annotation from the record
        sciName = rec.annotations.get('organism')
        accession = rec.id

        # for each feature
        for feature in rec.features:
            # only the 'source' types will contain the fields of interest
            if feature.type == 'source':
                # get the taxid and the strain information
                taxNum = feature.qualifiers.get("db_xref")[0]
                strain = feature.qualifiers.get("strain")

                # format strain into a str
                if strain is not None:
                    strain = strain[0]
                else:
                    strain = ''
                
                # once found, stop iterating through any remaining features
                break
        
        # build the new record and add it to the list
        idStr = accession + SEP + sciName + SEP + strain + SEP + taxNum
        newRec = SeqRecord(rec.seq, id=idStr, description='')
        recordsL.append(newRec)
    
    # write the records to file
    SeqIO.write(recordsL, fastaFN, 'fasta')




def __makeblastdbFromFna(fnaFN:str, dbFN:str, blastExecutDirPath:str) -> None:
    """ makeblastdbFromFna:
            Accepts a filename for a fasta containing the query nucleotide seq-
            uences, a filename indicating where to write a blast database, and
            the directory containing the blast+ executables as inputs. Constru-
            cts a database from the fasta and saves it at the specified locati-
            on. Does not return.
    """
    # constants
    EXECUTABLE = 'makeblastdb'
    NT = "nucl"

    # construct the 'makeblastdb' command
    executable = os.path.join(blastExecutDirPath, EXECUTABLE)
    makeblastdb = NcbimakeblastdbCommandline(cmd=executable, dbtype=NT, \
                                                    input_file=fnaFN, out=dbFN)

    # execute the command
    makeblastdb()


def __runRnaBlast(queryFna:str, blastdb:str, outFile:str, \
                                               blastExecutDirPath:str) -> None:
    """ runBlast:
            Accepts a nucleotide fasta file, a nucleotide blastdb file, the lo-
            cation to save the blastn results, and the directory containing the
            blast+ executables as inputs. Blasts the query fasta against the
            database, and saves the results at the specified location. Does not
            return.
    """
    # constants
    PRNT_1 = "Using blastn to search for related 16S rRNA gene sequences ... "
    DONE = "Done.\n"
    HEADERS = ['qseqid', 'stitle', 'pident', 
               'length', 'mismatch', 'gaps', 
               'gapopen', 'evalue', 'bitscore']
    OUTFMT = '6'        # unsafe to adjust value
    EXECUTABLE = 'blastn'
    QCOV_HSP_PERC = '70'
    EVALUE_CUTOFF = '1e-5'

    # make logger
    logger = logging.getLogger(__name__ + "." + __runRnaBlast.__name__)

    # make the outfmt string
    outfmtStr = __makeOutfmtString(OUTFMT, HEADERS)

    # construct the 'blastn' command
    executable = os.path.join(blastExecutDirPath, EXECUTABLE)
    blastn = NcbiblastnCommandline(cmd = executable,
                                   db = blastdb,
                                   query = queryFna,
                                   qcov_hsp_perc = QCOV_HSP_PERC,
                                   evalue = EVALUE_CUTOFF,
                                   outfmt = outfmtStr,
                                   out = outFile)

    # execute the command with the shell (runs blastn)
    print(PRNT_1, end='', flush=True)
    logger.info(PRNT_1)
    blastn()
    print(DONE)
    logger.info(DONE)


def __makeOutfmtString(outfmt:str, headers:list) -> str:
    """ makeOutfmtString:
            Accepts an outfmt number (as a string) and a list of headers to be
            included in the blast output. Returns a string that can be used wi-
            th the '-outfmt' flag of blastn.
    """
    if len(headers) == 0:
        # if headers is empty, then just return outfmt
        outStr = outfmt
    else:
        # open quote and add the outfmt number
        outStr = '"'
        outStr += outfmt

        # add each header separated by a space
        for head in headers:
            outStr += " " + head
        
        # close quote
        outStr += '"'
    
    # return
    return outStr


def __makeHeaderString(headers:list) -> str:
    """ makeHeaderString:
            Accepts a list of headers and converts it into a tab-delimited 
            string. Returns the string.
    """
    # initialize output string
    outStr = ''

    # add each header to the string separated by a tab then return
    for head in headers:
        outStr += head
        outStr += "\t"
    return outStr


