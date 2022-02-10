# prokTaxBuilder

## Running the software
### To run the software, perform the following steps:
  * Make a working directory containing an input gbff
  * Modify ```<path to prokTaxBuilder>/param.py``` appropriately
  * From within the working directory, call the following command on the commandline:

        $ python3.7 -i <path to prokTaxBuilder>/main.py

  * Within the interactive python session, run the following command:

        >>> runPt1("<genbank filename>.gbff")

  * After it finishes running, open ```initialAnalysis/putativePhylogeneticMarkers.txt``` and determine which CDS you wish to use as the phylogenetic marker of the group. Make a note of the gene number as this will be used as input in the next step.
  * Within the interactive python session, run the following command:

        >>> runPt2("<genbank filename>.gbff", <gene number:int>)
  
  * The results can be found in the folder ```finalAnalysis```.

## Description of the workflow
### Part 1: Determining a rough taxonomic estimate and performing initial phylogenomic analyses
  * Extract 16S rRNA gene sequences from the input genbank.
  * Use the 16S rRNA gene sequences and BLASTn to determine which taxonomic orders are related to the input.
  * Retreive data from NCBI's Taxonomy database to create a single taxonomy encompassing all the identified orders.
  * Reconcile the NCBI-based taxonomy with LPSN's taxonomic structure.
  * Import data from NCBI's Assembly database for the species represented in the taxonomy.
  * Use the BLASTn results to determine which taxa could serve as the ingroup/outgroup.
  * Download assemblies, the number to download is specified in ```<path to prokTaxBuilder>/param.py```.
  * Calculate core genes for the set of gbff files.
  * Create a species tree based on a concatenated alignment of the core genes and individual gene trees for each of the core genes.
  * Use the species tree and gene trees to calculate a score of how well correlated each individual core gene is to the species tree topology. Save this information in the file ```initialAnalysis/putativePhylogeneticMarkers.txt```.

### Part 2: Refining genome selection and performing final phylogenomic analyses
  * Using the gene specified by the user, find the most closely-related taxa to the input genome.
  * Use this knowledge to select genomes from:
      * many species from the closest genus/genera
      * several species from more distant genera
  * Download the assemblies.
  * Calculate the core genes for the refined set of genomes.
  * Create a species tree based on a concatenated alignment of the core genes.
  * Calculate the average amino acid identity (AAI) for the species in the tree.
  * Calculate the average nucleotide identity (ANI) for the species in the tree.
