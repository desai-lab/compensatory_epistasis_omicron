from Bio.Seq import Seq
import pandas as pd
import re
import numpy as np


def rc(x):
    return str(Seq(x).reverse_complement())


def translate(x):
    return str(Seq(x).translate())


seq_0 = snakemake.config['ancestral_seq'].strip().upper()
seq_1 = snakemake.config['descendant_seq'].strip().upper()
starting_resnum = snakemake.config['sequence_begins_with_residue']

for seq in [seq_0, seq_1]:
    assert len(seq_0) % 3 == 0, 'Ancestral seq length should be multiple of 3'
    assert len(seq_1) % 3 == 0, 'Descendant seq length should be multiple of 3'
    assert len(seq_0) == len(
        seq_1), 'Ancestral and Descendant seqs should be the same length'

mutation_list_file = snakemake.output.mutation_list

synonymously_coded_mutations = snakemake.config['synonymously_coded_mutations']


codons_0 = re.findall('...', seq_0)
codons_1 = re.findall('...', seq_1)

mutations = []
num_mutations = 0
for i, (codon_0, codon_1) in enumerate(zip(codons_0, codons_1)):
    if codon_0 == codon_1:
        continue

    aa_0 = translate(codon_0)
    aa_1 = translate(codon_1)
    if aa_0 != aa_1:
        resnum = starting_resnum + i
        mut_name = f'{aa_0}{resnum}{aa_1}'
        num_mutations += 1
        mutation = (mut_name, i*3, codon_0, codon_1, aa_0, aa_1, resnum)
        mutations.append(mutation)


mutations_df = (pd.DataFrame(mutations,
                             index=range(1, len(mutations) + 1),
                             columns=['mut_name', 'nt_pos', 'codon_0', 'codon_1', 'aa_0', 'aa_1', 'aa_pos'])
                .set_index('mut_name')
                )


# Check that the mutation names provided are consistent with
# the ones present between the sequences
assert set([a.strip() for a in mutations_df.index]) == set([a.strip() for a in snakemake.config['mutation_names']]
                                                           ), 'Mutations detected in sequences are not consistent with those in config file.\nMutations detected:' + '\n'.join(mutations_df.index) + '\nConfig file mutations:' + '\n'.join(snakemake.config['mutation_names'])
assert num_mutations == int(
    snakemake.config['num_mutations']), f'Number of mutations in detected sequences ({num_mutations}) is not consistent with config `num_mutations`'

mutations_df['synonymous'] = False

for syn_mut in synonymously_coded_mutations:
    mut_name = syn_mut['mutation_name']
    assert mut_name in mutations_df.index, "Synonymous mutation {mut_name} is not in list of mutations"

    nt_pos = 3 * (syn_mut['synonomous_position'] - starting_resnum)

    mutations_df.loc[mut_name, 'nt_pos'] = nt_pos
    mutations_df.loc[mut_name, 'codon_0'] = syn_mut['codon_0']
    mutations_df.loc[mut_name, 'codon_1'] = syn_mut['codon_1']
    mutations_df.loc[mut_name, 'aa_0'] = translate(syn_mut['codon_0'])
    mutations_df.loc[mut_name, 'aa_1'] = translate(syn_mut['codon_1'])
    mutations_df.loc[mut_name, 'aa_pos'] = syn_mut['synonomous_position']
    mutations_df.loc[mut_name, 'synonymous'] = True

mutations_df.to_csv(mutation_list_file, sep='\t')
