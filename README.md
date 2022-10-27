
# PHANTASM: PHylogenomic ANalyses for the TAxonomy and Systematics of Microbes

## A docker image for this software can be found at https://hub.docker.com/r/jwirth/phantasm

## This software is still under active development.

## If you use our software, please cite our paper
**Automating microbial taxonomy workflows with PHANTASM: PHylogenomic ANalyses for the TAxonomy and Systematics of Microbes**

Joseph S. Wirth & Eliot C. Bush, 2022

This manuscript is currently in revision. A preprint can be found on [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.10.18.512716v1).

##### Note: PHANTASM requires an email address for communicating with NCBI as described in python's [Bio.Entrez](https://biopython.org/docs/latest/api/Bio.Entrez.html) package. The email address is not stored or used for any other purposes.

## Installing PHANTASM
### Using the Docker image
The easiest way to get started is is to use the [Docker image](https://hub.docker.com/r/jwirth/phantasm). There are detailed instructions for using phantasm as a Docker container on [dockerhub](https://hub.docker.com/r/jwirth/phantasm). The Docker image contains all of the dependencies pre-installed and only requires that Docker Desktop is installed (or `docker` if using a linux server).


### Native installation
PHANTASM requires `python3` (v. 3.9+). This can be installed via the command line with the following command:

    sudo apt-get -y install python3.9 python3-pip python3-setuptools python3-dev
          
`mummer` is needed for calculating ANI using [pyani](https://github.com/widdowquinn/pyani). It can be installed with the following command:

    sudo apt-get -y install mummer

`curl` is necessary for certain R libraries to function properly. It can be installed with the following command:

    sudo apt-get -y install curl libcurl4-openssl-dev

`blast+` is required for calculating core genes and average amino acid identity. It can be installed with the following command:

    sudo apt-get -y install ncbi-blast+

`R` needs to be installed for plotting heatmaps of the overall genome-related indices. It can be installed with the following command:

    sudo apt-get -y install r-base

The following R libraries need to be installed: `ape`, `gplots`, `dendextend`, and `BiocManager`. It also needs the `DECIPHER` bioconductor package. These can be installed using the following commands. _**Warning**: These packages can take several minutes to install._

    R -e "install.packages('ape',dependencies=TRUE, repos='http://cran.rstudio.com/')"
    R -e "install.packages('gplots',dependencies=TRUE, repos='http://cran.rstudio.com/')"
    R -e "install.packages('dendextend',dependencies=TRUE, repos='http://cran.rstudio.com/')"
    R -e "install.packages('BiocManager',dependencies=TRUE, repos='http://cran.rstudio.com/')"
    R -e "BiocManager::install('DECIPHER')"

The following python packages need to be installed: `numpy`, `scipy`, `Bio`, `textdistance`, `parasail`, `rpy2`, `pyani`. They can be installed with `pip3` using the following commands:

    pip3 install scipy
    pip3 install Bio
    pip3 install textdistance
    pip3 install numpy
    pip3 install parasail
    pip3 install rpy2
    pip3 install pyani

PHANTASM requires [FastTree](https://www.microbesonline.org/fasttree/) (or FastTreeMP), [muscle v5.1](https://github.com/rcedgar/muscle/releases/), and [iqtree](https://github.com/Cibiv/IQ-TREE/releases) to be installed. It also needs the source code for [xenoGI](https://github.com/ecbush/xenoGI/releases) v3.1.0.

Finally, PHANTASM and its dependency, [xenoGI](https://github.com/ecbush/xenoGI), need to be cloned onto your machine.

In order for PHANTASM to function properly, you will need to modify the file `param.py` found in the `phantasm` directory. Specifically, you need to modify the following fields as shown in the table below:

    | variable name |                    variable value                     |
    |---------------|-------------------------------------------------------|
    | BLASTPLUS_DIR | absolute path to the blast+ executable directory (bin)|
    | MUSCLE_EXE    | absolute path the the muscle executable file          |
    | FASTTREE_EXE  | absolute path to the fasttree executable file         |
    | IQTREE_EXE    | absolute path to the IQTree executable file           |
    | PHANTASM_DIR  | absolute path to the phantasm directory               |
    | XENOGI_DIR    | absolute path to the xenoGI directory                 |



## Running the software
### To get set up, perform the following steps:
  * Make a working directory containing either a single input genome in genbank file format or a directory of multiple input genomes in genbank file format
  * From within the working directory, call the following command on the commandline (modify ```<path to phantasm>/param.py``` appropriately):

        $ python3 <path to phantasm>/phantasm.py help

## When a phylogenetic marker and reference genomes are unknown
If you have an annotated genome sequence but you do not know a suitable phylogenetic marker and you do not know which reference genomes would be appropriate to use, then PHANTASM can be run in two steps:

1. It uses 16S rRNA gene sequences to estimate a taxonomic placement, then prioritizes taxonomic breadth while in order to find genes that coevolve with the input genome(s) and its relatives.
2. The user selects a phylogenetic marker and uses it to prioritize the taxonomic placement of the input genome(s). It strives to select suitable reference genomes and then uses them to perform phylogenomic analyses.

### Optional: Bypassing the requirement for annotated 16S rRNA gene sequences in your input genome(s)
By default, PHANTASM will attempt to extract the 16S rRNA gene sequence(s) from your input genome(s) and blastn against the bacterial and archaeal 16S sequences from  [NCBI's Targetted Loci database](https://www.ncbi.nlm.nih.gov/refseq/targetedloci/). The blast results are then used to approximate the taxonomic placement of the input genome(s).

Annotated 16S rRNA gene sequences are not always available. It is possible for users to bypass this requirement by providing a file with NCBI Taxonomy IDs for the taxa that the input genomes may be related to. The specified taxa can be either any NCBI Taxonomy classification as long as it is of rank Family or below (eg. Genus and Species are acceptable, but Order and Class are not).

In order for PHANTASM to find these Taxonomy IDs, you must provide them in a file named  `taxids.txt`  in the directory where you are calling PHANTASM (ie. `./`). The file must contain exactly one NCBI Taxonomy ID per line. The taxids provided by the user cannot exceed the rank Family. This methodology allows the user to specify a collection of NCBI Taxonomy IDs that are somewhat related to the input genome(s). The user is free to obtain these IDs in whatever method they best see fit.

Keep in mind that annotated 16S rRNA gene sequences are  **not**  required if a phylogenetic marker and/or suitable reference genomes are already known (see other sections below).

### To identify putative phylogenetic markers, run the following command:

    $ python3 <path to phantasm>/phantasm.py getPhyloMarker <gbff filename> <email address>
 
  * After it finishes running, open ```initialAnalysis/putativePhylogeneticMarkers.txt``` and determine which gene(s) you wish to use as the phylogenetic marker(s) of the group. Make a note of the gene number or the locus tag as this will be used as input in the next step.
  * If using PHANTASM on multiple input genomes, replace ```<gbff filename>``` with ```<gbff directory>```

    * example for a single genome:

          $cd ~/workdir
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

### To refining the phylogeny and perform phylogenomic analyses, run one of the following commands:

    $ python3 <path to phantasm>/phantasm.py refinePhylogeny --locus_tag <locus tag> <gbff file> <email address>

  OR

    $ python3 <path to phantasm>/phantasm.py refinePhylogeny --gene_num <gene number> <gbff file> <email address>
  
  * The results can be found in the folder `finalAnalysis`.
  
  * example using one genome and one phylogenetic marker:

        $ python3 ~/phantasm/phantasm.py refinePhylogeny --locus_tag SPO_RS17765 R_pomeroyi.gbff email@address.org
  
  * example using one genome and multiple (two) phylogenetic markers:

        $ python3 ~/phantasm/phantasm.py refinePhylogeny --locus_tag SPO3507,SPO0155 R_pomeroyi.gbff email@address.org

  * example using multiple genomes and one (per genome) phylogenetic marker:

        $ python3 ~/phantasm/phantasm.py refinePhylogeny --locus_tag RUA4292_04770,SPO3507 input_genomes email@address.org
  
  * example using multiple genomes and multiple (two per genome) phylogenetic markers:

        $ python3 ~/phantasm/phantasm.py refinePhylogeny --locus_tag RUA4292_04770,SPO3507,RUA4292_00884,SPO0155 input_genomes email@address.org

## When phylogenetic marker(s) are known but suitable reference genomes are not
If you have an annotated genome sequence and you know suitable phylogenetic markers, but you do not know which reference genomes would be appropriate to use, then PHANTASM can be ran in a single step

    $ python3 <path to phantasm>/phantasm.py knownPhyloMarker <locus tag> <email address>
The locus tag(s) provided must be present within the input genome(s).
  * The results can be found in the folder ```finalAnalysis```.
  * If using PHANTASM on multiple input genomes, replace ```<gbff filename>``` with ```<gbff directory>```
  * If using PHANTASM with multiple phylogenetic markers:
      * List each locus tag separated by a comma (no spaces allowed)
      * Include the same number of locus tags for each input genome

  * example using one genome and one phylogenetic marker:

        $ python3 ~/phantasm/phantasm.py knownPhyloMarker SPO_RS17765 R_pomeroyi.gbff email@address.org

  * example using one genome and multiple (two) phylogenetic markers:

        $ python3 ~/phantasm/phantasm.py knownPhyloMarker SPO3507,SPO0155 R_pomeroyi.gbff email@address.org

  * example using multiple genomes and one (per genome) phylogenetic marker:

        $ python3 ~/phantasm/phantasm.py knownPhyloMarker RUA4292_04770,SPO3507 input_genomes email@address.org
  
  * example using multiple genomes and multiple (two per genome) phylogenetic markers:

        $ python3 ~/phantasm/phantasm.py knownPhyloMarker RUA4292_04770,SPO3507,RUA4292_00884,SPO0155 input_genomes email@address.org

## When suitable reference genomes are already known
If you already know which genomes you want to analyze, then you can tell PHANTASM to skip the reference genome selection step and only do the phylogenomic analyses. First, you must create a "human map file". This file tells PHANTASM what the human-readable names are (used in the trees and the alignments) for each genome file. It should be formatted with exactly two columns separated by a single tab per line. The first column is the filename of the genbank file (not including the path), and the second column is the human readable name. For example:

    GCF_009496005.1_ASM949600v1_genomic.gbff	Tritonibacter_litoralis__invalid|SM1979|GCF_009496005
    GCF_007923355.1_ASM792335v1_genomic.gbff	Phaeobacter_marinintestinus__invalid|UB-M7|GCF_007923355
    GCF_900106805.1_IMG-taxon_2693429872_annotated_assembly_genomic.gbff	Ruegeria_halocynthiae|DSM_27839_T|GCF_900106805
    GCF_000511355.1_ASM51135v1_genomic.gbff	Leisingera_methylohalidivorans|DSM_14336__MB2|GCF_000511355
    GCF_900172225.1_Tr.litoreus_CECT7639_Spades_Prokka_genomic.gbff	Ruegeria_litorea|CECT_7639_T|GCF_900172225
    GCF_900102795.1_IMG-taxon_2622736501_annotated_assembly_genomic.gbff	Epibacterium_ulvae|U95_T|GCF_900102795
    GCF_001458295.1_Ruegeria_sp._CECT5091_Spades_Prokka_genomic.gbff	Ruegeria_denitrificans|CECT_5091_T|GCF_001458295
    GCF_001507545.1_ASM150754v1_genomic.gbff	Ruegeria_profundi|ZGT108_T|GCF_001507545
    GCF_009617595.1_ASM961759v1_genomic.gbff	Tritonibacter_aquimaris__invalid|SM1969|GCF_009617595
    my_assembly.gbk	my_input_genome

Once you have made this file, then you can call the following command to run PHANTASM:


    $ python3 <path>/phantasm.py analyzeGenomes <input genomes dir> <human map file> <output directory> <email address>
##### Note: the specified output directory must not already exist.

For example:

    $ ls
    input_genomes humanMap.txt
    $ python3 ~/phantasm/phantasm.py analyzeGenomes ./input_genomes ./humanMap.txt ./results email@address.org

* The results can be found in the folder you specified (`./results` in the example above).

## Optional: Excluding taxa from phylogenomic analyses
It is possible to exclude specific taxa from the refined phylogeny. To do so, create a file named ```excludedTaxids.txt``` containing exactly one NCBI Taxonomy id per line. The excluded ids can represent any taxonomic rank (eg. species, genus, family, etc.). This file should be in the same directory in which PHANTASM is called (```~/workdir``` in the examples above). This feature is largely experimental and should be used with caution.



## Detailed descriptions of the workflows
### Option 1a: Identifying a suitable phylogenetic marker for your input genome
  * Extracts the 16S rRNA gene sequences from the input genbank.
  * Uses the 16S rRNA gene sequences and BLASTn to search [NCBI's Targetted Loci database](https://www.ncbi.nlm.nih.gov/refseq/targetedloci/) in order to determine which taxonomic orders are related to the input genome(s).
  * Retreives data from NCBI's Taxonomy database to create a single taxonomy encompassing all the identified orders.
  * Reconciles the NCBI-based taxonomy with LPSN's taxonomic structure.
  * Imports data from NCBI's Assembly database for the species represented in the taxonomy.
  * Uses the BLASTn results and the taxonomic structure to determine which taxa could serve as the ingroup/outgroup.
  * Downloads assemblies for the ingroup and the outgroup the number to download is specified in ```<path to phatasm>/param.py``` in the ```MAX_LEAVES``` variable.
  * Calculates core genes for the set of gbff files.
  * Creates a species tree based on a concatenated alignment of the core genes.
  * Creates an individual gene tree for each core gene.
  * Uses the species tree and gene trees to calculate a score of how well correlated each individual core gene is to the species tree's topology. Saves this information in the file ```initialAnalysis/putativePhylogeneticMarkers.txt```.

### Option 1b: Refining the phylogeny and performing phylogenomic analyses
  * Using the phylogenetic marker(s) specified by the user, finds the most closely-related taxa to the input genome.
  * Uses this information to select genomes from:
      * many species from the closest genus/genera
      * several species from more distant genera
  * Downloads the assemblies for the new ingroup and the outgroup
  * Calculates the core genes for the refined set of genomes.
  * Creates a species tree based on a concatenated alignment of the core genes.
  * Calculates the average amino acid identity (AAI) for the species in the tree (excluding the outgroup).
  * Calculates the average nucleotide identity (ANI) for the species in the tree (excluding the outgroup).

### Option 2: Using a known phylogenetic marker(s)
  * Using the phylogenetic marker(s) specified by the user, finds the most closely-related taxa to the input genome via BLASTp against either NCBI's nr or refseq_protein databases.
  * Uses the BLASTp results to determine which taxonomic orders are related to the input.
  * Retreives data from NCBI's Taxonomy database to create a single taxonomy encompassing all the identified orders.
  * Reconciles the NCBI-based taxonomy with LPSN's taxonomic structure.
  * Imports data from NCBI's Assembly database for the species represented in the taxonomy.
  * Uses the BLASTp results and the taxonomic structure to determine which taxa could serve as the ingroup/outgroup.
  * Downloads assemblies for the ingroup and the outgroup the number to download is specified in `<path to phantasm>/param.py` in the ```MAX_LEAVES``` variable.
  * Calculates core genes for the set of gbff files.
  * Creates a species tree based on a concatenated alignment of the core genes.
  * Calculates the average amino acid identity (AAI) for the species in the tree (excluding the outgroup).
  * Calculates the average nucleotide identity (ANI) for the species in the tree (excluding the outgroup).

### Option 3: Analyzing user-specified genomes
  * Calculates core genes for the set of gbff files.
  * Creates a species tree based on a concatenated alignment of the core genes.
  * Calculates the average amino acid identity (AAI) for the species in the tree (excluding the outgroup).
  * Calculates the average nucleotide identity (ANI) for the species in the tree (excluding the outgroup).


## Dependencies
### External software
  * [FastTree](http://www.microbesonline.org/fasttree/)
  * [MUSCLE](https://www.drive5.com/muscle/)
  * [blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
  * [IQTree](http://www.iqtree.org/)
  * [mummer](https://github.com/mummer4/mummer)

### Python modules (python3.9 or above)
  * [xenoGI v3.1.0](https://github.com/ecbush/xenoGI/releases)
  * scipy
  * Bio
  * textdistance
  * numpy
  * parasail
  * rpy2
  * [pyani](https://github.com/widdowquinn/pyani)

### R packages (R version 4.1.1 or above)
  * ape
  * gplots
  * dendextend
  * DECIPHER
