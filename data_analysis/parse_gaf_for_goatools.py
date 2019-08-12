import pandas as pd

gafpath = '../original_data/goa_human.gaf'

df = pd.read_csv(gafpath, sep='\t', skiprows=12, header=None)
df = df[[1, 4]]
df = df.groupby(1)[4].apply(lambda x: ';'.join(x))
df.to_csv('../data/go_enrichment/parsed_goa_human.gaf', sep='\t')