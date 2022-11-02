# PHANTASM: PHylogenomic ANalyses for the TAxonomy and Systematics of Microbes
## Developed by Joseph S. Wirth

## Source code is available at https://github.com/dr-joe-wirth/phantasm

## This software is still under active development.

## If you use our software, please cite our paper
**Automating microbial taxonomy workflows with PHANTASM: PHylogenomic ANalyses for the TAxonomy and Systematics of Microbes**

Joseph S. Wirth & Eliot C. Bush, 2022

This manuscript is currently in revision. A preprint can be found on [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.10.18.512716v1).
‎
##### Note: PHANTASM requires an email address for communicating with NCBI as described in python's [Bio.Entrez](https://biopython.org/docs/latest/api/Bio.Entrez.html) package. The email address is not stored or used for any other purposes.
##### Note: PHANTASM does not always run properly with the current version of `podman` (3.4.7). It is unclear what is causing this issue.
‎
 
 
# Table of Contents
1. [Mounting the Docker image](#1-mounting-the-docker-image)

    1.1. Running Docker as root

    1.2 Preparing Docker settings

    1.3 Mounting the image as a container

2. [Getting started with PHANTASM](#2-getting-started-with-phantasm)

    2.1. Getting help

    2.2. Modifying PHANTASM's settings (optional)

    2.3. Excluding specific taxa from phylogenomic analyses (optional)

3. [Running PHANTASM](#3-running-phantasm)

    3.1. Option 1: unknown reference genomes and unknown phylogenetic marker(s)

    3.2. Option 2: unknown reference genomes and known phylogenetic marker(s)

    3.3. Option 3: known reference genomes

4. [Analyzing the results](#4-analyzing-the-results)

5. [Detailed descriptions of the workflows](#5-detailed-descriptions-of-the-workflows)

    5.1. Option 1: Identifying suitable phylogenetic markers for your input genome

    5.2. Option 1: Refining the phylogeny and performing phylogenomic analyses

    5.3. Option 2: Using a known phylogenetic marker(s)

    5.4. Option 3: Analyzing a set of user-specified genomes

6. [Common error messages](#6-common-error-messages)

    6.1. Problems with 16S rRNA gene sequence annotation (or lack thereof) in your input genome(s)

    6.2. Failed to connect to NCBI database(s)

# 1. Mounting the Docker image
## 1.1. Running Docker as `root`
Docker must be run as `root`. This should not be a problem if [Docker Desktop](https://www.docker.com/products/docker-desktop/) is installed on your computer. However, running with the command line may require the use of `sudo` (`su` to switch to `root`). To circumvent this issue, add yourself to the group `docker` as described [here](https://docs.docker.com/engine/install/linux-postinstall/):

    $ sudo groupadd docker
    $ sudo usermod -aG docker $USER

Performing these steps will allow Docker to run as `root` without requiring `sudo` or `su`.

## 1.2. Preparing Docker settings
Be sure to allocate at least 8gb of memory to docker.
  * If you cannot allocate 8gb, then you will need to reduce the number of taxa (see section 2.2)
  * If you see the following error, then you have not allocated enough memory:

        Calculating core genes ... Killed

Be sure to allocate more processors to docker if you wish to take advantage of parallel processing
  * This is optional. PHANTASM can run successfully with only 1 CPU, but will run significantly faster with more.
  * If you wish to use parallel processing, you will also need to modify the `param.py` file as described in section 2.2.

## 1.3. Mounting the image as a container:
Pull the docker image from docker hub

    $ docker pull jwirth/phantasm:latest

Make a working directory containing an input gbff file or a directory containing multiple gbff files (note: all gbff files MUST be annotated!). For example:

    $ mkdir ~/myworkdir
    $ cp my_genome.gbff ~/myworkdir
    $ ls ~/myworkdir
    my_genome.gbff

Mount the image as a container (named `myContainer` below) and designate your working directory as a volume within it (named `/mydata` below). You must provide the __ABSOLUTE PATH__ to your working directory.

    $ docker run -it --name myContainer -v /Users/<USERNAME>/myworkdir:/mydata jwirth/phantasm

The following names __cannot__ be used for the volume (`/mydata` in the example above):
  * `/bin`
  * `/boot`
  * `/dev`
  * `/etc`
  * `/exec`
  * `/home`
  * `/lib`
  * `/lib64`
  * `/media`
  * `/mnt`
  * `/opt`
  * `/phantasm`
  * `/proc`
  * `/root`
  * `/run`
  * `/sbin`
  * `/srv`
  * `/sys`
  * `/tmp`
  * `/usr`
  * `/var`
  * `/xenoGI-3.1.0`

If you have successfully mounted the image, then you should see something like this:

    root@fe46a3c61f0b:/#

You should also be able to see your genome file of interest within the container. For example:

    root@fe46a3c61f0b:/# ls /mydata
    my_genome.gbff

# 2. Getting started with PHANTASM
## 2.1. Getting help:
In order to get a help message, use the following command:

    root@fe46a3c61f0b:/# phantasm help

## 2.2. Modifying PHANTASM's settings (optional)
Use the `nano` command to open the `param.py` file

    root@fe46a3c61f0b:/# nano /phantasm/param.py

Change the values for `NUM_PROCESSORS`, `MAX_LEAVES`, `REDUCE_NUM_CORE_GENES`, `BOOTSTRAP_FINAL_TREE`, and/or `NUM_BOOTSTRAPS` in the file. For your convenience, the relevant contents of the file are shown below:

    # specify the number of processors to use
    NUM_PROCESSORS:int = 1

    # specify the maximum number of taxa in a given analysis
    # do not set this value below 10
    MAX_LEAVES:int = 50

    # specify if the number of core genes used to calculate the final tree should be reduced
    #### `True` indicates yes
    #### `False` indicates no
    REDUCE_NUM_CORE_GENES:bool = False

    # specify if the final tree should have bootstrap supports
    #### WARNING: Bootstrapping trees will significantly increase run times (several hours -> several days)
    #### `True` indicates yes
    #### `False` indicates no
    BOOTSTRAP_FINAL_TREE:bool = False

    # specify the number of bootstrap supports for the final species tree
    #### this is only relevant if BOOTSTRAP_FINAL_TREE is True
    NUM_BOOTSTRAPS:int = 100
    . . .

To close the file, press `CTRL + X`. You will be asked if you wish to save the file; press `Y`. Press `ENTER` to confirm that the saved changes will overwrite the existing `param.py` file.

### Important caveats for modifying these values:
  * Although you can set `NUM_PROCESSORS` to any positive integer value, PHANTASM will ultimately be limited by the number of processors allocated to Docker.
  * `MAX_LEAVES` can be any positive integer, but very small values (<20) are not recommended. Keep in mind that run-times and RAM scale _exponentially_ with the number of leaves.
  * If the runtime of your phylogenetic tree is a concern, then setting `REDUCE_NUM_CORE_GENES` may be useful. It omits any core genes where one or more taxa has >5% gaps in its alignment.
  * Changing `BOOTSTRAP_FINAL_TREE` to `True` will likely increase the run time from a few hours to several days. This is due to the long run times of IQTree as compared to FastTree.
  * `NUM_BOOTSTRAPS` is only relevant if `BOOTSTRAP_FINAL_TREE` is set to `True`.

## 2.3. Excluding specific taxa from phylogenomic analyses (optional)
It is possible to exclude specific taxa from the refined phylogeny. To do so, create a file named `excludedTaxids.txt` containing exactly one NCBI Taxonomy id per line in the directory where you call PHANTASM. The excluded ids can represent any taxonomic rank (eg. species, genus, family, etc.). This option is only relevant is using option 1 or option 2 (see below). This file should be in the same directory in which PHANTASM is called (ie. `/mydata` in the example shown in section 2.1). This feature is largely experimental and should be used with caution.

# 3. Running PHANTASM
There are three ways of running PHANTASM. As show in the table below, the best option to use largely depends on the answer to two questions: 1) do you know a suitable set of genomes? and 2) do you know a suitable phylogenetic marker(s) for your genome(s) of interest?

     _____________________________________________________________________
    | reference genomes? | phylogenetic marker(s)? | best PHANTASM option |
    |--------------------|-------------------------|----------------------|
    |         no         |           no            |       option 1       |
    |         no         |           yes           |       option 2       |
    |         yes        |        irrelevant       |       option 3       |
    |____________________|_________________________|______________________|


## 3.1. Option 1: unknown reference genomes and unknown phylogenetic marker(s)
If you have an annotated genome sequence but you do not know a suitable phylogenetic marker and you do not know which reference genomes would be appropriate to use, then PHANTASM can be run in two steps:

1. PHANTASM uses 16S rRNA gene sequences to estimate a taxonomic placement, then prioritizes taxonomic breadth in order to find genes that coevolve with the input genome(s) and its relatives.
2. The user Looks at the results of this analysis and selects a suitable phylogenetic marker(s).
3. Using the selected marker(s), PHANTASM automatically selects a set of reference genomes and uses them to perform phylogenomic analyses.

### Optional: Bypassing the requirement for annotated 16S rRNA gene sequences in your input genome(s)

By default, PHANTASM will attempt to extract the 16S rRNA gene sequence(s) from your input genome(s) by relying on its accurate annotation in your input genome. It searches this gene with blastn against the bacterial and archaeal 16S sequences from [NCBI's Targetted Loci database](https://www.ncbi.nlm.nih.gov/refseq/targetedloci/). The blast results are then used to approximate the taxonomic placement of the input genome(s).

Annotated 16S rRNA gene sequences are not always available. It is possible for users to bypass this requirement by providing a file with NCBI Taxonomy IDs for the taxa that the input genomes may be related to. The specified taxa can be either any NCBI Taxonomy classification as long as it is of rank Family or below (eg. Genus and Species are acceptable, but Order and Class are not).

In order to bypass the 16S rRNA gene sequence requirement, you must provide these NCBI Taxonomy IDs a file named `taxids.txt` in the directory where you are calling PHANTASM (eg. `/mydata` in the examples below). The file must contain exactly one NCBI Taxonomy ID per line. The taxids provided by the user cannot exceed the rank Family. This methodology allows you to specify a collection of NCBI Taxonomy IDs that are somewhat related to the input genome(s). You are free to obtain these IDs in whatever method you best see fit. Ultimately, PHANTASM will use a suitable phylogenetic marker to refine the set of reference genomes.

Keep in mind that annotated 16S rRNA gene sequences are **not** required if a phylogenetic marker and/or suitable reference genomes are already known (see sections 3.2 and 3.3 for more information).

### Step 1: Identifying putative phylogenetic markers. . .
Make sure you are currently within the directory containing the gbff file (called `mydata` in this example)

    root@fe46a3c61f0b:/# cd /mydata

You are now ready to identify phylogenetic markers. Be sure to replace `my_genome.gbff` and `email@address.org` with appropriate values. If using multiple genomes, replace `my_genome.gbff` with the name of the directory containing the genome files (this directory must be nested within the working directory).

    root@fe46a3c61f0b:/mydata# phantasm getPhyloMarker my_genome.gbff email@address.org

If using PHANTASM on multiple input genomes, replace `<gbff filename>` with `<gbff directory>`. For example:

    root@fe46a3c61f0b:/mydata# ls input_genomes
    R_atlantica.gbff R_pomeroyi.gbff
    root@fe46a3c61f0b:/mydata# phantasm getPhyloMarker input_genomes email@address.org

##### Note: this step will take some time, especially when performing blastp.
‎
### Step 2: Selecting phylogenetic marker(s)
 After PHANTASM finishes running, examine the file `initialAnalysis/putativePhylogeneticMarkers.txt` and determine which gene(s) you wish to use as the phylogenetic marker(s).  In this example, it will also exist locally at `~/myworkdir/initialAnalysis/putativePhylogeneticMarkers.txt`. Core genes are ranked by their cophenetic correlation coefficients ("r<sub>ccc</sub>") and are listed from highest r<sub>ccc</sub> to lowest r<sub>ccc</sub>. Record either the `locus_tag` or the `gene_num` for the gene you wish to use as a phylogenetic marker. __As a general rule, an r<sub>ccc</sub> value of at least 0.9 is recommended for any gene used as a phylogenetic marker.__

Here is an example of what the first few lines of `putativePhylogeneticMarkers.txt` looks like:

    cophenetic_corr_coef    protein_len     gene_num        locus_tag       gene_name       annotation
    0.9891495196945915      1188    43248   SPO_RS03315     dnaE    dnaE - DNA polymerase III subunit alpha
    0.987983296592394       1379    46112   SPO_RS17770     rpoB    rpoB - DNA-directed RNA polymerase subunit beta
    0.9856682366935172      732     45731   SPO_RS15845             primosomal protein N'
    0.9854566492450753      1151    44689   SPO_RS10520     mfd     mfd - transcription-repair coupling factor
    0.9816288307386607      972     44830   SPO_RS11250     uvrA    uvrA - excinuclease ABC subunit UvrA
    0.9791067341269494      776     45773   SPO_RS16055     clpA    clpA - ATP-dependent Clp protease ATP-binding subunit ClpA
    0.9784173031511715      913     44666   SPO_RS10405     gyrA    gyrA - DNA gyrase subunit A
    . . .
‎

### Step 3: Refining the initial phylogeny and perform phylogenomic analyses. . .
First, make sure you are currently within the directory containing the gbff file (called `mydata` in this example)

    root@fe46a3c61f0b:/# cd /mydata

Next, run one of the following two commands depending on what information you recorded.
  * If you are using a locus tag, then run the following command. Be sure to replace `<locus tag>`, `my_genome.gbff`, and `email@address.org` with appropriate values.

        root@fe46a3c61f0b:/mydata# phantasm refinePhylogeny --locus_tag <locus tag> my_genome.gbff email@address.org

  * If you are using a gene number, then run the following command. Be sure to replace `<integer>`, `my_genome.gbff`, and `email@address.org` with appropriate values:

        root@fe46a3c61f0b:/mydata# phantasm refinePhylogeny --gene_num <integer> my_genome.gbff email@address.org

PHANTASM can handle multiple phylogenetic markers and/or multiple input genomes. **The examples below use locus tags, but PHANTASM works analogously when specifying gene numbers.**

  * example using one genome and one phylogenetic marker:

        root@fe46a3c61f0b:/mydata# phantasm refinePhylogeny --locus_tag SPO_RS17765 R_pomeroyi.gbff email@address.org

  * example using one genome and multiple (two) phylogenetic markers:

        root@fe46a3c61f0b:/mydata# python3 ~/phantasm/phantasm.py refinePhylogeny --locus_tag SPO3507,SPO0155 R_pomeroyi.gbff email@address.org

  * example using multiple genomes and one (per genome) phylogenetic marker:

        root@fe46a3c61f0b:/mydata# python3 ~/phantasm/phantasm.py refinePhylogeny --locus_tag RUA4292_04770,SPO3507 input_genomes email@address.org

  * example using multiple genomes and multiple (two per genome) phylogenetic markers:

        root@fe46a3c61f0b:/mydata# python3 ~/phantasm/phantasm.py refinePhylogeny --locus_tag RUA4292_04770,SPO3507,RUA4292_00884,SPO0155 input_genomes email@address.org

##### Notes:
  * this step will take some time as there are many slow steps, in particular:
    * blastp
    * calculating AAI
    * calculating ANI
    * making species trees (especially if bootstraps are involved)
  * If using multiple genomes:
    * replace `my_genome.gbff` with a directory containing the genome files
    * list at least one gene per genome as a comma-separated list (**no spaces allowed**)
  * If using multiple phylogenetic markers
    * list each gene as a comma-separated list (**no spaces allowed**)
    * can be used in conjunction with multiple genomes, but the you must specify a gene for each genome


## 3.2. Option 2: unknown reference genomes and known phylogenetic marker(s)
If you have an annotated genome sequence and you know suitable phylogenetic markers, but you do not know which reference genomes would be appropriate to use, then you can run PHANTASM in a single step.

First, make sure you are currently within the directory containing the gbff file (called `mydata` in this example)

    root@fe46a3c61f0b:/# cd /mydata

Next, run the following command (be sure to replace `<locus_tag>`, `my_genome.gbff`, and `email@address.org` with appropriate values):

    root@fe46a3c61f0b:/mydata# phantasm knownPhyloMarker <locus_tag> my_genome.gbff email@address.org

The locus tag(s) provided must be present within the input genome(s) and indicate the phylogenetic marker(s) that you wish to use. PHANTASM can handle multiple phylogenetic markers and/or multiple input genomes. If using multiple markers, they should be a comma-separated list **without spaces**. If using multiple input genomes, the same number of phylogenetic markers must be specified for each input genome.

### Examples:
  * one genome and one phylogenetic marker:

        root@fe46a3c61f0b:/mydata# phantasm knownPhyloMarker SPO_RS17765 my_genome.gbff email@address.org

  * one genome and multiple (two) phylogenetic markers:

        root@fe46a3c61f0b:/mydata# phantasm knownPhyloMarker SPO3507,SPO0155 my_genome.gbff email@address.org

  * multiple genomes and one (per genome) phylogenetic marker:

        root@fe46a3c61f0b:/mydata# phantasm knownPhyloMarker RUA4292_04770,SPO3507 input_genomes email@address.org
  
  * multiple genomes and multiple (two per genome) phylogenetic markers:

        root@fe46a3c61f0b:/mydata# phantasm knownPhyloMarker RUA4292_04770,SPO3507,RUA4292_00884,SPO0155 input_genomes email@address.org

### Notes:
  * this step will take some time as there are many slow steps, in particular:
    * blastp
    * calculating AAI
    * calculating ANI
    * making species trees (especially if bootstraps are involved)
  * If using multiple genomes:
    * replace `my_genome.gbff` with a directory containing the genome files
    * list at least one gene per genome as a comma-separated list (**no spaces allowed**)
  * If using multiple phylogenetic markers
    * list each gene as a comma-separated list
    * can be used in conjunction with multiple genomes, but the you must specify the same number of genes for each genome

## 3.3. Option 3: known reference genomes
It is also possible to run PHANTASM on a set of user-specified genomes. This can be accomplished with a single command.

First, make sure you are currently within the directory containing the gbff file (called `mydata` in this example)

    root@fe46a3c61f0b:/# cd /mydata

Next, you must create a "human map file". This file tells PHANTASM what the human-readable names are (used in the trees and the alignments) for each genome file. It should be formatted with exactly two columns separated by a single tab. The first column is the filename of the genbank file (not including the path), and the second column is the human readable name. **Important note: _the desired outgroup must be listed last in the human map file_.**  For example:

    GCF_009496005.1_ASM949600v1_genomic.gbff	Tritonibacter_litoralis__invalid|SM1979|GCF_009496005
    GCF_007923355.1_ASM792335v1_genomic.gbff	Phaeobacter_marinintestinus__invalid|UB-M7|GCF_007923355
    GCF_900106805.1_IMG-taxon_2693429872_annotated_assembly_genomic.gbff	Ruegeria_halocynthiae|DSM_27839_T|GCF_900106805
    GCF_000511355.1_ASM51135v1_genomic.gbff	Leisingera_methylohalidivorans|DSM_14336__MB2|GCF_000511355
    GCF_900172225.1_Tr.litoreus_CECT7639_Spades_Prokka_genomic.gbff	Ruegeria_litorea|CECT_7639_T|GCF_900172225
    GCF_900102795.1_IMG-taxon_2622736501_annotated_assembly_genomic.gbff	Epibacterium_ulvae|U95_T|GCF_900102795
    GCF_001458295.1_Ruegeria_sp._CECT5091_Spades_Prokka_genomic.gbff	Ruegeria_denitrificans|CECT_5091_T|GCF_001458295
    GCF_001507545.1_ASM150754v1_genomic.gbff	Ruegeria_profundi|ZGT108_T|GCF_001507545
    GCF_009617595.1_ASM961759v1_genomic.gbff	Tritonibacter_aquimaris__invalid|SM1969|GCF_009617595
    my_assembly.gbk	Ruegeria_sp
    outgroup_genome_seq.gbff	outgroup_name

Once this file has been created, run the following command to analyze the provided genomes:

    root@fe46a3c61f0b:/mydata# phantasm analyzeGenomes <input genomes directory> <human map file> <output directory> <email address>
##### Note: the specified output directory must not already exist.

# 4. Analyzing the results
At this point, all results should have been generated and it is safe to close and remove the docker container (named `myContainer` in this example:

          root@fe46a3c61f0b:/# exit
          $ docker rm myContainer

All data can be found within the volume that you mounted into the container (`~/myworkdir` in this example). Unless the command `analyzeGenomes` was used, the results will be found in the folder `finalAnalysis`. If `analyzeGenomes` was used, then the results will be found in the specified output directory.
  * The following files are likely to be the most useful for phylogenomic and taxonomic analyses:
    * `speciesTree.nwk`: this is the species tree, it is rooted on the outgroup, and the outgroup is present in the tree
    * `speciesTree_outgroupPruned.nwk`: this is the same species tree as above, but the outgroup has been pruned to allow better resolution of the relevant phylogenomic relationships
    * `aai_matrix.txt`: this is a raw text file containing the average amino acid identities for all of the taxa in `speciesTree_outgroupPruned.nwk`
    * `aai_heatmap.pdf`: this is a visualization of the data in `aai_matrix.txt` and the taxa are ordered to match the order of the taxa in `speciesTree_outgroupPruned.nwk`. Due to a bug in the R package `gplots`, the tree cannot be plotted alongside the heatmap at this time.
    * `ani_matrix.txt`: this is a raw text file containing the average nucleotide identities for all of the taxa in `speciesTree_outgroupPruned.nwk`
    * `ani_heatmap.pdf`: this is a visualization of the data in `ani_matrix.txt` and the taxa are ordered to match the order of the taxa in `speciesTree_outgroupPruned.nwk`. Due to a bug in the R package `gplots`, the tree cannot be plotted alongside the heatmap at this time.
    * `coreGenesSummary.txt`: this is a tab-delimited file containing detailed information on the core genes used to construct the species tree. Locus tags, gene numbers, gene names, and annotations are only specified for the input genomes. The indicated alignment file can be used to determine these data for other reference genomes used in the analysis.
    * `wgsHumanMap.txt`: this file lists all of the accession numbers for the genomes used in the analysis as well as their human readable names found in the alignments, trees, and heatmaps. The equivalent of this file was provided by the user if `analyzeGenomes` was used.

# 5. Detailed descriptions of the workflows

## 5.1. Option 1: Identifying suitable phylogenetic markers for your input genome

* Extracts the 16S rRNA gene sequences from the input genbank.

* Uses the 16S rRNA gene sequences and BLASTn to search [NCBI's Targetted Loci database](https://www.ncbi.nlm.nih.gov/refseq/targetedloci/) in order to determine which taxonomic orders are related to the input genome(s).

* Retreives data from NCBI's Taxonomy database to create a single taxonomy encompassing all the identified orders.

* Reconciles the NCBI-based taxonomy with LPSN's taxonomic structure.

* Imports data from NCBI's Assembly database for the species represented in the taxonomy.

* Uses the BLASTn results and the taxonomic structure to determine which taxa could serve as the ingroup/outgroup.

* Downloads assemblies for the ingroup and the outgroup the number to download is specified in `<path to phatasm>/param.py` in the `MAX_LEAVES` variable.

* Calculates core genes for the set of gbff files.

* Creates a species tree based on a concatenated alignment of the core genes.

* Creates an individual gene tree for each core gene.

* Uses the species tree and gene trees to calculate a score of how well correlated each individual core gene is to the species tree's topology. Saves this information in the file `initialAnalysis/putativePhylogeneticMarkers.txt`.

## 5.2. Option 1: Refining the phylogeny and performing phylogenomic analyses
  * Using the phylogenetic marker(s) specified by the user, finds the most closely-related taxa to the input genome.
  * Uses this information to select genomes from:
  * many species from the closest genus/genera
  * several species from more distant genera
  * Downloads the assemblies for the new ingroup and the outgroup
  * Calculates the core genes for the refined set of genomes.
  * Creates a species tree based on a concatenated alignment of the core genes.
  * Calculates the average amino acid identity (AAI) for the species in the tree (excluding the outgroup).
  * Calculates the average nucleotide identity (ANI) for the species in the tree (excluding the outgroup).

## 5.3. Option 2: Using a known phylogenetic marker(s)
  * Using the phylogenetic marker(s) specified by the user, finds the most closely-related taxa to the input genome via BLASTp against either NCBI's nr or refseq_protein databases.
  * Uses the BLASTp results to determine which taxonomic orders are related to the input.
  * Retrieves data from NCBI's Taxonomy database to create a single taxonomy encompassing all the identified orders.
  * Reconciles the NCBI-based taxonomy with LPSN's taxonomic structure.
  * Imports data from NCBI's Assembly database for the species represented in the taxonomy.
  * Uses the BLASTp results and the taxonomic structure to determine which taxa could serve as the ingroup/outgroup.
  * Downloads assemblies for the ingroup and the outgroup the number to download is specified in `<path to phantasm>/param.py` in the `MAX_LEAVES` variable.
  * Calculates core genes for the set of gbff files.
  * Creates a species tree based on a concatenated alignment of the core genes.
  * Calculates the average amino acid identity (AAI) for the species in the tree (excluding the outgroup).
  * Calculates the average nucleotide identity (ANI) for the species in the tree (excluding the outgroup).

## 5.4. Option 3: Analyzing a set of user-specified genomes
  * Calculates core genes for the set of gbff files.
  * Creates a species tree based on a concatenated alignment of the core genes.
  * Calculates the average amino acid identity (AAI) for the species in the tree (excluding the outgroup).
  * Calculates the average nucleotide identity (ANI) for the species in the tree (excluding the outgroup).

# 6. Common error messages
## 6.1. Problems with 16S rRNA gene sequence annotation (or lack thereof) in your input genome(s)
### Possible error messages
    BaseException: Could not extract 16S rRNA gene sequences from the provided genbank file.

    RuntimeError: The file '/mydata/initialAnalysis/16S/query.16SrRNA.blastn' does not contain any valid blastn hits

### Possible solutions
__Option 1.__ Create a file named `taxids.txt` as described in the optional section above.

__Option 2.__ Use a known phylogenetic marker as described above.

## 6.2. Failed to connect to NCBI database(s)
### Possible error messages
    urllib.error.HTTPError: HTTP Error 400: Bad Request

### Possible solutions
This error message only occurs when PHANTASM fails to connect to one or more NCBI databases. Currently, the only option is to rerun PHANTASM from the beginning