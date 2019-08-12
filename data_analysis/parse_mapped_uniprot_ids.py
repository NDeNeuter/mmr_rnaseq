import glob

import pandas as pd

def parse_directory(dirpath):
    
    dct = {'name': [], 'mapped_count': [], 'unmapped_count': []}
    
    for fpath in glob.glob(dirpath+'/*_list.txt'):
        print(fpath)
        
        map_path = fpath.replace('/degs/', '/uniprot_ids/').replace('_list.txt', '_map.txt')
        
        mapped_genes_path = fpath.replace('_list.txt', '_mapped_genes.txt')
        unmapped_genes_path = fpath.replace('/degs/', '/non-coding_genes/').replace('_list.txt', '_unmapped_genes.txt')
        protein_path = map_path.replace('_map.txt', '_protein_ids.txt')
        
        mapped_genes, unmapped_genes = check_mapped(fpath, map_path)
        protein_ids = parse_mapping_file(map_path)
        
        with open(mapped_genes_path, 'w') as o:
            o.write('\n'.join(mapped_genes))

        with open(unmapped_genes_path, 'w') as o:
            o.write('\n'.join(unmapped_genes))
            
        with open(protein_path, 'w') as o:
            o.write('\n'.join(protein_ids))
        
        dct['name'].append(fpath.split('/')[-1])
        dct['mapped_count'].append(len(mapped_genes))
        dct['unmapped_count'].append(len(unmapped_genes))
        
        df = pd.DataFrame(dct)
        df.to_csv('../data/go_enrichment/non-coding_genes/summary.tsv', index=False, sep='\t')
        
        
def parse_mapping_file(fpath, organism='Homo sapiens (Human)', status='reviewed'):
    
    try:
        df = pd.read_csv(fpath, sep='\t')
    except pd.errors.EmptyDataError:
        return []
    filtered_df = df[(df['Organism']==organism) & (df['Status']==status)]
    
    return filtered_df['Entry'].values

    
def check_mapped(inputpath, mapped_path, organism='Homo sapiens (Human)', status='reviewed'):
    
    with open(inputpath, 'r') as f:
        input_genes = set(f.read().split('\n'))
    
    try:
        df = pd.read_csv(mapped_path, sep='\t')
    except pd.errors.EmptyDataError:
        return [], []
    
    uniprot_genes = set(df[(df['Organism']==organism) & (df['Status']==status)].iloc[:, 0].values)
    
    mapped_genes = {x for x in input_genes.intersection(uniprot_genes) if x != ''}
    unmapped_genes = {x for x in input_genes.difference(uniprot_genes) if x != ''} # input_genes - uniprot_genes 
    
    return mapped_genes, unmapped_genes

        
if __name__ == '__main__':
    
    parse_directory('../data/go_enrichment/degs')
    
    