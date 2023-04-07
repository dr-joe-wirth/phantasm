# modifying these variables may cause phantasm to fail
BLASTPLUS_DIR:str = '/usr/bin'
MUSCLE_EXE:str = '/exec/muscle'
FASTTREE_EXE:str = '/exec/FastTreeMP'
IQTREE_EXE:str = '/exec/iqtree/bin/iqtree'
PHANTASM_DIR:str = '/phantasm'
XENOGI_DIR:str = '/xenoGI-3.1.1'

# modifying these variables may cause phantasm to fail
import os
CSV_1:str = os.path.join(PHANTASM_DIR, 'PHANTASM', 'lpsn_data', 'speciesGenus.csv')
CSV_2:str = os.path.join(PHANTASM_DIR, 'PHANTASM', 'lpsn_data', 'genusParents.csv')
CSV_3:str = os.path.join(PHANTASM_DIR, 'PHANTASM', 'lpsn_data', 'validPhylumClassOrderFamily.csv')
CSV_4:str = os.path.join(PHANTASM_DIR, 'PHANTASM', 'lpsn_data', 'synonymsPhylumClassOrderFamily.csv')
