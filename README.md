# PHANTASM: PHylogenomic ANalyses for the TAxonomy and Systematics of Microbes

## A docker image for this software can be found at https://hub.docker.com/r/jwirth/phantasm

## If you use this software, please email me at jwirth \<at> g \<dot> hmc \<dot> edu

## Running the software
### To run the software, perform the following steps:
  * Make a working directory containing either a single input genome in genbank file format or a directory of multiple input genomes in genbank file format
  * Modify ```<path to phantasm>/param.py``` appropriately
  * From within the working directory, call the following command on the commandline:

        $ python3 <path to phantasm>/phantasm.py help

## Identifying a suitable phylogenetic marker for your input genome
### Run the following command:

        $ python3 <path to phantasm>/phantasm.py getPhyloMarker <gbff filename> <email address>
 
  * After it finishes running, open ```initialAnalysis/putativePhylogeneticMarkers.txt``` and determine which gene you wish to use as the phylogenetic marker of the group. Make a note of the gene number or the locus tag as this will be used as input in the next step.
  * If using PHANTASM on multiple input genomes, replace ```<gbff filename>``` with ```<gbff directory>```

    * example for a single genome:

          $ cd ~/workdir
          $ ls
          R_pomeroyi.gbff
          $ python3 ~/phantasm/phantasm.py getPhyloMarker R_pomeroyi.gbff email@address.org      
          $ less ~/workdir/initialAnalysis/putativePhylogeneticMarkers.txt
    
    * example for multiple genomes:

          $ cd ~/workdir
          $ ls
          input_genomes
          $ ls input_genomes
          R_atlantica.gbff        R_pomeroyi.gbff
          $ python3 ~/phantasm/phantasm.py getPhyloMarker input_genomes email@address.org
          $ less ~/workdir/initialAnalysis/putativePhylogeneticMarkers.txt

## Refining the phylogeny and performing phylogenomic analyses
### Run one of the following commands:

        $ python3 <path to phantasm>/phantasm.py refinePhylogeny --locus_tag <locus tag> <gbff file> <email address>

  OR

        $ python3 <path to phantasm>/phantasm.py refinePhylogeny --gene_num <gene number> <gbff file> <email address>
  
  * The results can be found in the folder ```finalAnalysis```.
  
  * example using one genome and one phylogenetic marker:

        $ python3 ~/phantasm/phantasm.py refinePhylogeny --locus_tag SPO_RS17765 R_pomeroyi.gbff email@address.org
  
  * example using one genome and multiple (two) phylogenetic markers:

        $ python3 ~/phantasm/phantasm.py refinePhylogeny --locus_tag SPO3507,SPO0155 R_pomeroyi.gbff email@address.org

  * example using multiple genomes and one (per genome) phylogenetic marker:

        $ python3 ~/phantasm/phantasm.py refinePhylogeny --locus_tag RUA4292_04770,SPO3507 input_genomes email@address.org
  
  * example using multiple genomes and multiple (two per genome) phylogenetic markers:

        $ python3 ~/phantasm/phantasm.py refinePhylogeny --locus_tag RUA4292_04770,SPO3507,RUA4292_00884,SPO0155 input_genomes email@address.org

### Optional: Excluding taxa from the refined phylogeny
It is possible to exclude specific taxa from the refined phylogeny. To do so, create a file named ```excludedTaxids.txt``` containing exactly one NCBI Taxonomy id per line. The excluded ids can represent any taxonomic rank (eg. species, genus, family, etc.). This file should be in the same directory in which PHANTASM is called (```~/workdir``` in the examples above). This feature is largely experimental and should be used with caution.

## Description of the workflow
### Part 1: Identifying a suitable phylogenetic marker for your input genome
  * Extracts the 16S rRNA gene sequences from the input genbank.
  * Uses the 16S rRNA gene sequences and BLASTn to determine which taxonomic orders are related to the input.
  * Retreives data from NCBI's Taxonomy database to create a single taxonomy encompassing all the identified orders.
  * Reconciles the NCBI-based taxonomy with LPSN's taxonomic structure.
  * Imports data from NCBI's Assembly database for the species represented in the taxonomy.
  * Uses the BLASTn results and the taxonomic structure to determine which taxa could serve as the ingroup/outgroup.
  * Downloads assemblies for the ingroup and the outgroup the number to download is specified in ```<path to phatasm>/param.py``` in the ```MAX_LEAVES``` variable.
  * Calculates core genes for the set of gbff files.
  * Creates a species tree based on a concatenated alignment of the core genes.
  * Creates an individual gene tree for each core gene.
  * Uses the species tree and gene trees to calculate a score of how well correlated each individual core gene is to the species tree's topology. Saves this information in the file ```initialAnalysis/putativePhylogeneticMarkers.txt```.

### Part 2: Refining the phylogeny and performing phylogenomic analyses
  * Using the phylogenetic marker(s) specified by the user, finds the most closely-related taxa to the input genome.
  * Uses this information to select genomes from:
      * many species from the closest genus/genera
      * several species from more distant genera
  * Downloads the assemblies for the new ingroup and the outgroup
  * Calculates the core genes for the refined set of genomes.
  * Creates a species tree based on a concatenated alignment of the core genes.
  * Calculates the average amino acid identity (AAI) for the species in the tree (excluding the outgroup).
  * Calculates the average nucleotide identity (ANI) for the species in the tree (excluding the outgroup).

## Dependencies
### External software
  * FastTree (http://www.microbesonline.org/fasttree/)
  * MUSCLE (https://www.drive5.com/muscle/)
  * blast+ (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
  * IQTree (http://www.iqtree.org/)
  * mummer (https://github.com/mummer4/mummer)

### Python modules (python3.7 or above)
  * xenoGI (https://github.com/ecbush/xenoGI)
  * scipy
  * Bio
  * textdistance
  * numpy
  * parasail
  * rpy2
  * pyani

### R packages (R version 4.1.1 or above)
  * ape
  * gplots
  * dendextend
  * DECIPHER
