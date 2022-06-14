import regex
import pandas as pd
import numpy as np

n_subs = snakemake.config['num_substitutions_tolerated_in_inline_idx']
col_idx_regex = regex.compile('(%s){s<=%d}' % (snakemake.params.col_inline_idx, n_subs))
row_idx_regex = regex.compile('(%s){s<=%d}' % (snakemake.params.row_inline_idx, n_subs))


nb_reads = 0
reads_thrown_out = 0
reads_with_no_index_err = 0
reads_with_corrected_index_err = 0

with open(snakemake.input.tsv) as f, open(snakemake.output.tsv, 'w') as fw:
    f.readline() # remove header
    fw.write("\t".join(["UMI_1", "UMI_2", "read_1", "read_2"]) + "\n")
    for line in f:
        col_idx, row_idx, UMI_1, UMI_2, read_1, read_2 = line.strip().split("\t")
        nb_reads += 1
        # Match the row and column inline indices to the fuzzy regex
        # for the correct row/column index for this sample
        match_1 = col_idx_regex.match(col_idx)
        match_2 = row_idx_regex.match(row_idx)

        if match_1 is None or match_2 is None:
            reads_thrown_out += 1
            continue

        if match_1.fuzzy_counts[0] == 0 and match_2.fuzzy_counts[0] == 0:
            reads_with_no_index_err += 1
        else:
            reads_with_corrected_index_err += 1

        fw.write("\t".join([UMI_1, UMI_2, read_1, read_2]) + "\n")

# Write stats
(pd.DataFrame({
    'sample' : [snakemake.wildcards.sample],
    'reads_with_corrected_index_err' : [reads_with_corrected_index_err],
    'reads_with_no_index_err' : [reads_with_no_index_err],
    'reads_thrown_out' : [reads_thrown_out],
    'fraction_of_index_accepted' : [(reads_with_corrected_index_err +
                                     reads_with_no_index_err) / nb_reads
                                    ] if nb_reads != 0 else np.nan,
})
 .set_index('sample')
 .to_csv(snakemake.output.stats, sep='\t')
)
