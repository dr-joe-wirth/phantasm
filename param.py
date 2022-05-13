# specify the number of processors to use
NUM_PROCESSORS:int = 1

# specify the maximum number of taxa in a given analysis
MAX_LEAVES:int = 50

# specify the number of bootstrap supports for the final species tree
NUM_BOOTSTRAPS:int = 100

# specify the path to the folder containing the blast+ executables
BLASTPLUS_DIR:str = '/usr/local/ncbi/blast/bin'

# specify the path to the MUSCLE executable
MUSCLE_EXE:str = '/usr/local/bin/muscle'

# specify the path to the FastTree executable (FastTreeMP is acceptable)
FASTTREE_EXE:str = '/usr/local/bin/FastTree'

# specify the path to the IQTree exectuable
IQTREE_EXE:str = '/usr/local/bin/iqtree'

# specify the path to the phantasm directory
PHANTASM_DIR:str = '/phantasm'

# specify the path to the xenoGI directory
XENOGI_DIR:str = '/xenoGI'

# specify the paths to each of the required lpsn csv files
# MODIFYING THIS MAY CAUSE PHANTASM TO FAIL!
import os
CSV_1:str = os.path.join(PHANTASM_DIR, 'PHANTASM', 'lpsn_data', 'lpsn_gss_2022-02-01.csv')
CSV_2:str = os.path.join(PHANTASM_DIR, 'PHANTASM', 'lpsn_data', 'genusParents.csv')
CSV_3:str = os.path.join(PHANTASM_DIR, 'PHANTASM', 'lpsn_data', 'validPhylumClassOrderFamily.csv')
CSV_4:str = os.path.join(PHANTASM_DIR, 'PHANTASM', 'lpsn_data', 'synonymsPhylumClassOrderFamily.csv')
