# importing things
import pandas as pd
import numpy as np
import scipy.optimize
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300

conf = snakemake.config

replicate_names = snakemake.params.replicate_names
reps = {}
for r, tsv in zip(replicate_names, snakemake.input.tsv):
    reps[r] = pd.read_csv(tsv, sep="\t", dtype={'geno': str})\
        .set_index("geno")\
        .loc[:,['Mean fluorescence expression', 'Std fluorescence expression']]\
        .add_suffix('_' + r)

# merging the replicates
df = pd.concat([reps[r] for r in replicate_names], axis=1)

# normalize expression
for rep in replicate_names:
    df[f"expression_norm_{rep}"] = df[f"Mean fluorescence expression_{rep}"]\
        / np.mean(df[f"Mean fluorescence expression_{rep}"])

# compute mean expression
df["expression_norm"] = df.apply(lambda r:
                                 np.nanmean([r[f"expression_norm_{rep}"] for rep in replicate_names
                                             if rep in replicate_names]),
                                 axis=1)

# save tsv
df.to_csv(snakemake.output.tsv, sep="\t")

# plot the comparison between replicate for expression
try:
    rep1, rep2 = [(a, b)
                  for a in replicate_names for b in replicate_names if a != b][0]
except IndexError:
    # if there's only one replicate, still plot the "jointplot",
    # but they won't be very informative
    rep1, rep2 = replicate_names[0], replicate_names[0]

if conf['do_plots']:
    sns.jointplot(x=f"Mean fluorescence expression_{rep1}",
                  y=f"Mean fluorescence expression_{rep2}",
                  data=df, kind="hex", bins='log')
    plt.title(f"Expression (mean bin) of {rep1} vs {rep2}")
else:
    plt.subplots()
plt.savefig(snakemake.output.plot_replicate_comparison_expr)
