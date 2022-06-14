import pandas as pd

count_tables = []
for sample,fname in zip(snakemake.params.samples, snakemake.input.tsv):
    count_tables.append(pd.read_table(fname, dtype=dict(geno='str'))
                        .set_index('geno')
                        .rename(columns=lambda x: sample))
    
(pd.concat(count_tables, axis=1)
 .fillna(0)
 .astype(int)
 .sort_index()
 .to_csv(snakemake.output.tsv, sep='\t')
)
