import pandas as pd

(pd.concat([
    pd.read_table(statsfile, index_col=[0]) for statsfile in snakemake.input.stats])
 .to_csv(snakemake.output.stats, sep='\t')
)