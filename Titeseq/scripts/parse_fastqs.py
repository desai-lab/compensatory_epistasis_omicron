import pandas as pd
import numpy as np
import gzip
from Bio import SeqIO
import regex

read_1_regex = regex.compile(open(snakemake.input.read_1_regex).read())
read_2_regex = regex.compile(open(snakemake.input.read_2_regex).read())


total_reads_considered = 0
reads_matching_pattern = 0

with open(snakemake.output.tsv, 'w') as fw:
    fw.write("\t".join(
                     ["col_idx", "row_idx", "UMI_1", "UMI_2",
                      "read_1", "read_2"]) + "\n")
    # Go through both input fastq files simultaneously line by line
    with gzip.open(snakemake.input.read_1_fastq, 'rt') as read_1_handle:
        with gzip.open(snakemake.input.read_2_fastq, 'rt') as read_2_handle:
            for (read_1, read_2) in zip(
                    SeqIO.parse(read_1_handle, 'fastq'),
                    SeqIO.parse(read_2_handle, 'fastq')):
                if snakemake.config['debug_mode'] and total_reads_considered >= snakemake.config['max_reads_per_file_if_debugging']:
                    break
                total_reads_considered += 1

                match_1 = read_1_regex.match(str(read_1.seq))
                match_2 = read_2_regex.match(str(read_2.seq))

                if match_1 is None or match_2 is None:
                    continue
                reads_matching_pattern += 1

                UMI_1 = match_1.group('UMI')
                inline_idx_1 = match_1.group('index')
                read_1_seq = match_1.group('read_seq')

                UMI_2 = match_2.group('UMI')
                inline_idx_2 = match_2.group('index')
                read_2_seq = match_2.group('read_seq')

                # write to file
                fw.write("\t".join(
                    [inline_idx_1, inline_idx_2, UMI_1,
                     UMI_2, read_1_seq, read_2_seq]
                ) + "\n")


# Write stats
try:
    frac_accepted = reads_matching_pattern / total_reads_considered
except ZeroDivisionError:
    frac_accepted = np.nan
(pd.DataFrame({
    'sample' : [snakemake.wildcards.sample],
    'fastq_file' : [snakemake.input.read_1_fastq],
    'total_reads_considered': [total_reads_considered],
    'reads_matching_pattern': [reads_matching_pattern], 
    'fraction_of_fastq_accepted' : [frac_accepted]
})
 .set_index('sample')
 .to_csv(snakemake.output.stats, sep='\t')
)
