import regex
import pandas as pd
from Bio.Seq import Seq
import numpy as np
import gzip

geno_for_regex = regex.compile(open(snakemake.input.geno_for_regex).read())
geno_rev_regex = regex.compile(open(snakemake.input.geno_rev_regex).read())


mutations_df = pd.read_table(snakemake.input.mutation_list)
codon_0_list = mutations_df['codon_0']
codon_1_list = mutations_df['codon_1']

def extract_mutations_from_match_pair(for_match, rev_match):
    for_mutations = extract_mutations(for_match)
    rev_mutations = extract_mutations(rev_match)
    # Take the reverse complement of the mutations on reverse read
    # and append to the forward read mutations
    return for_mutations + [rc(codon) for codon in reversed(rev_mutations)]

def extract_mutations(match):
    return [codon for (groupname,codon) in match.groupdict().items()
            if 'mut_' in groupname]

class CodonDoesntMatchError(ValueError):
    pass

def rc(x):
    return str(Seq(x).reverse_complement())

def mutations_to_binary_geno(muts):
    s = ''
    for mut,codon_0,codon_1 in zip(muts,codon_0_list,codon_1_list):
        if mut == codon_0:
            s += '0'
        elif mut == codon_1:
            s += '1'
        else:
            raise CodonDoesntMatchError
    return s

unmatched_reads = 0
reads_with_err_in_mutation = 0
reads_with_no_err = 0
reads_with_indels = 0
reads_with_one_sub = 0
reads_with_two_or_more_subs = 0
total_accepted = 0
nb_reads = 0

with \
     open(snakemake.input.tsv, 'r') as f, \
     open(snakemake.output.tsv, 'w') as fw, \
     gzip.open(snakemake.output.bad_reads, 'wt') as fbad:

    f.readline() # remove header line
    fw.write("\t".join(['UMI_1', 'UMI_2', 'geno']) + "\n")
    fbad.write("\t".join(['read_1', 'read_2']) + "\n")
    for line in f:
        UMI_1, UMI_2, read_1, read_2 = line.strip().split("\t")

        match_1 = geno_for_regex.match(read_1)
        match_2 = geno_rev_regex.match(read_2)
        nb_reads += 1

        if match_1 is None or match_2 is None:
            unmatched_reads += 1
            fbad.write(f"{read_1}\t{read_2}\n")
            continue

        mutation_codon_list = extract_mutations_from_match_pair(match_1,
                                                                match_2)

        try:
            geno = mutations_to_binary_geno(mutation_codon_list)
        except CodonDoesntMatchError:
            reads_with_err_in_mutation += 1
            # TODO we can fix the error if it is closer to one codon than the other
            continue

        sub1, in1, del1 = match_1.fuzzy_counts
        sub2, in2, del2 = match_2.fuzzy_counts
        num_indels = in1 + in2 + del1 + del2
        num_subs = sub1 + sub2

        if num_indels != 0:
            reads_with_indels += 1
        elif num_subs == 1:
            reads_with_one_sub += 1
        elif num_subs >= 2:
            reads_with_two_or_more_subs += 1
        else:
            reads_with_no_err += 1

        total_accepted += 1
        fw.write("\t".join((UMI_1, UMI_2, geno)) + "\n")



# Write stats
(pd.DataFrame({
    'sample' : [snakemake.wildcards.sample],
    'unmatched_reads' : [unmatched_reads],
    'reads_with_err_in_mutation' : [reads_with_err_in_mutation],
    'reads_with_indels' : [reads_with_indels],
    'reads_with_one_sub' : [reads_with_one_sub],
    'reads_with_two_or_more_sub' : [reads_with_two_or_more_subs],
    'reads_with_no_err' : [reads_with_no_err],
    'total_accepted' : [total_accepted],
    'fraction_of_geno_accepted' : [(total_accepted) / nb_reads] if nb_reads != 0 else np.nan,
})
 .set_index('sample')
 .to_csv(snakemake.output.stats, sep='\t')
)
