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
#### WARNING: Bootstrapping trees will significantly increase run times
#### `True` indicates yes
#### `False` indicates no
BOOTSTRAP_FINAL_TREE:bool = False

# specify the number of bootstrap supports for the final species tree
#### this is only relevant if BOOTSTRAP_FINAL_TREE is True
NUM_BOOTSTRAPS:int = 100

# modifying these variables may cause phantasm to fail
BLASTPLUS_DIR:str = '/usr/bin'
MUSCLE_EXE:str = '/exec/muscle'
FASTTREE_EXE:str = '/exec/FastTreeMP'
IQTREE_EXE:str = '/exec/iqtree/bin/iqtree'
PHANTASM_DIR:str = '/phantasm'
XENOGI_DIR:str = '/xenoGI-3.1.1'

# modifying these variables may cause phantasm to fail
import os
CSV_1:str = os.path.join(PHANTASM_DIR, 'PHANTASM', 'lpsn_data', 'lpsn_gss_2022-02-01.csv')
CSV_2:str = os.path.join(PHANTASM_DIR, 'PHANTASM', 'lpsn_data', 'genusParents.csv')
CSV_3:str = os.path.join(PHANTASM_DIR, 'PHANTASM', 'lpsn_data', 'validPhylumClassOrderFamily.csv')
CSV_4:str = os.path.join(PHANTASM_DIR, 'PHANTASM', 'lpsn_data', 'synonymsPhylumClassOrderFamily.csv')
