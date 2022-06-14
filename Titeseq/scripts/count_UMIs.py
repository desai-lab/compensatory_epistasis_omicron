import pandas as pd
import numpy as np

reads_df = pd.read_csv(snakemake.input.tsv, sep='\t', dtype=dict(geno='str'))

unique_reads = (
    reads_df
    .drop_duplicates(subset=['UMI_1','UMI_2'], keep='first')
    .drop(columns=['UMI_1','UMI_2'])
)

counts_per_geno = (
    unique_reads
    .assign(count=1)
    .groupby('geno')
    .sum()
    .sort_values(by='count', ascending=False)
)

# Write to file
counts_per_geno.to_csv(snakemake.output.tsv, sep='\t')

# Write stats
(pd.DataFrame({
    'sample' : [snakemake.wildcards.sample],
    'num_reads' : [len(reads_df)],
    'num_unique_UMIs' : [len(unique_reads)],
    'num_reads_dropped' : [len(reads_df) - len(unique_reads)],
    'frac_reads_dropped' : [(len(reads_df) - len(unique_reads)) / len(reads_df)] if len(reads_df) != 0 else np.nan,
    'num_genos_seen' : [len(counts_per_geno)]
})
 .set_index('sample')
 .to_csv(snakemake.output.stats, sep='\t')
)
