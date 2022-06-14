import math
from Bio.Seq import Seq
import pandas as pd
import regex
import numpy as np


def rc(x):
    return str(Seq(x).reverse_complement())


def translate(x):
    return str(Seq(x).translate())


seq_0 = snakemake.config['ancestral_seq'].strip().upper()
seq_1 = snakemake.config['descendant_seq'].strip().upper()
starting_resnum = snakemake.config['sequence_begins_with_residue']
num_mutations = snakemake.config['num_mutations']

num_mutations_on_forward_read = snakemake.config['num_mutations_on_forward_read']
num_mutations_on_reverse_read = snakemake.config['num_mutations_on_reverse_read']
for_primer = snakemake.config['foward_read_annealing_seq']
rev_primer = snakemake.config['reverse_read_annealing_seq']
forward_read_starting_offset = snakemake.config['forward_read_starting_offset']
reverse_read_starting_offset = snakemake.config['reverse_read_starting_offset']

UMI_length = snakemake.config['umi_length']
col_inline_idx = snakemake.config['col_inline_indices']
row_inline_idx = snakemake.config['row_inline_indices']

num_idx_subs = snakemake.config['num_substitutions_tolerated_in_inline_idx']
num_substitutions_tolerated_per_ten_bp = snakemake.config['num_substitutions_tolerated_per_ten_bp']
are_indels_tolerated_in_primers = snakemake.config['are_indels_tolerated_in_primers']

mutations_df = pd.read_table(snakemake.input.mutation_list, index_col=[0])


def make_or_regex(seqs, name, dist=0):
    '''Returns a regex string that is an OR of all the seqs in seqs.
    The capture is given a name `name`.
    Substitutions only are accepted within edit distance `dist`.
    '''
    ret = f"(?P<{name}>{'|'.join(seqs)})"
    if dist != 0:
        ret += '{s<=%d}' % dist
    return ret


UMI_regex = "(?P<UMI>[ACGT]{%d})" % UMI_length

cols_regex = make_or_regex(col_inline_idx, name='index', dist=num_idx_subs)
rows_regex = make_or_regex(row_inline_idx, name='index', dist=num_idx_subs)

read_seq_regex = "(?P<read_seq>[ACGTN]+)"

read_1_regex = UMI_regex + cols_regex + read_seq_regex
read_2_regex = UMI_regex + rows_regex + read_seq_regex

with open(snakemake.output.read_1_regex, 'w') as f:
    f.write(read_1_regex)

with open(snakemake.output.read_2_regex, 'w') as f:
    f.write(read_2_regex)

# we will trim or elongate the `amplicon` string to match what is amplified based on the for/rev primers and their provided offsets
amplicon = seq_0

if forward_read_starting_offset < 0:
    # amplicon begins before the provided gene sequence; we should extend the amplicon to the left
    # aaaaaa
    #     gggggggggg
    anneal_len = len(for_primer) - np.abs(forward_read_starting_offset)
    assert anneal_len > 0, "Forward primer is too short to anneal with gene sequence at provided offset"
    assert for_primer[-anneal_len:] == seq_0[:anneal_len], "Forward primer does not anneal provided gene sequence at provided offset"

    amplicon = for_primer[:-anneal_len] + amplicon
else:
    # amplicon begins somewhere in the provided gene sequence (I haven't tested this code)
    # we should trim amplicon from the left
    anneal_range = slice(forward_read_starting_offset,
                         forward_read_starting_offset + len(for_primer))
    assert seq_0[anneal_range] == for_primer, "Forward primer does not anneal provided gene sequence at provided offset"

    amplicon = amplicon[forward_read_starting_offset:]


if reverse_read_starting_offset < 0:
    # amplicon ends after the provided gene sequence; we should extend the amplicon to the right
    # gggggggggg
    #        aaaaa

    # The omicron project doesn't do this, so I haven't implemented it
    raise NotImplemented
else:
    # amplicon ends somewhere in the provided gene sequence; we should trim from the right
    anneal_range = slice(-reverse_read_starting_offset -
                         len(rev_primer), -reverse_read_starting_offset)
    assert seq_0[anneal_range] == rc(rev_primer), "Reverse primer does not anneal provided gene sequence at provided offset, reverse primer: " + \
        snakemake.config['reverse_read_annealing_seq'] + ", reverse primer (rc): " + rc(
            rev_primer) + ", gene sequence: " + seq_0[anneal_range]

    amplicon = amplicon[:-reverse_read_starting_offset]


mut_pos_in_amplicon = mutations_df['nt_pos'] - forward_read_starting_offset

constant_regions = []

for i in range(num_mutations + 1):
    if i == 0:
        start = 0
    else:
        start = mut_pos_in_amplicon[i - 1] + 3

    if i == num_mutations:
        end = len(amplicon)
    else:
        end = mut_pos_in_amplicon[i]

    c = amplicon[start:end]
    print(c)
    constant_regions.append(c)

constant_regions_rc = list(map(rc, reversed(constant_regions)))

for_read_constant_regions = constant_regions[:num_mutations_on_forward_read]
rev_read_constant_regions = constant_regions_rc[:num_mutations_on_reverse_read]

assert num_mutations_on_forward_read + \
    num_mutations_on_reverse_read == num_mutations

print(for_read_constant_regions)
print(rev_read_constant_regions)


def make_fuzzy_regex(search, num_err, name='', bestmatch=True, allow_indels=False):
    if search == '':
        return ''
    r = ''
    if bestmatch and num_err != 0:
        r += '(?e)'
    r += '('
    if name != '':
        r += f'?P<{name}>'
    r += search
    r += ')'
    if num_err != 0:
        r += '{%c<=%d}' % ('e' if allow_indels else 's', num_err)
    return r


def assemble_regex(const_regions):
    regex_parts = []
    for i, const_region in enumerate(const_regions):
        num_mut_tolerated = math.ceil(
            (len(const_region) / 10) * num_substitutions_tolerated_per_ten_bp)
        is_indel_tolerated = (i == 0) and are_indels_tolerated_in_primers

        regex_parts.append(make_fuzzy_regex(const_region,
                                            num_err=num_mut_tolerated,
                                            allow_indels=is_indel_tolerated,
                                            name=f'const_{i}'))

        regex_parts.append(make_fuzzy_regex('...',
                                            num_err=0,
                                            name=f'mut_{i}'))

    print('\n'.join(regex_parts))
    return ''.join(regex_parts)


geno_for_regex = assemble_regex(for_read_constant_regions)

geno_rev_regex = assemble_regex(rev_read_constant_regions)

with open(snakemake.output.geno_for_regex, 'w') as f:
    f.write(geno_for_regex)

with open(snakemake.output.geno_rev_regex, 'w') as f:
    f.write(geno_rev_regex)
