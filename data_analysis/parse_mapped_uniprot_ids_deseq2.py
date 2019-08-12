import glob
import os

import pandas as pd

dirpath = '../data/deseq2'

for dr, subdr, files in os.walk(dirpath):
    
    map_files = [x for x in files if x.endswith('_map.txt')]
    if len(map_files) == 1:
        print(dr, map_files)
        fpath = os.path.join(dr, map_files[0])
        
        try:
            df = pd.read_csv(fpath, sep='\t')
        except pd.errors.EmptyDataError:
            continue
        organism='Homo sapiens (Human)'
        status='reviewed'
        uniprot_genes = set(df[(df['Organism']==organism) & (df['Status']==status)].iloc[:, 1].values)
        
        with open(fpath.replace('_map.txt', '_uniprot_ids.txt'), 'w') as o:
            o.write('\n'.join(uniprot_genes))