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
        .add_suffix('_' + r)
    # remove the bad Kds
    reps[r].loc[(reps[r]['sigma_' + r] > conf["sigma_max"]) |
                (reps[r]['r2_' + r] < conf["r2_min"]), "log10Kd_" + r] = np.nan


# merging the replicates
df = pd.concat([reps[r] for r in replicate_names], axis=1)

# compute mean & ste
df["all_Kds"] = df.apply(lambda r:
                         [r[f"log10Kd_{rep}"] for rep in replicate_names
                          if not np.isnan(r[f"log10Kd_{rep}"])],
                         axis=1)
df["log10Kd"] = df.all_Kds.apply(np.mean)  # mean Kd
# unknown error if only one value
df.loc[df.all_Kds.apply(len) == 1, "err_log10Kd"] = np.nan
df.loc[df.all_Kds.apply(len) > 1, "err_log10Kd"] = df.all_Kds.apply(
    np.std) / np.sqrt(df.all_Kds.apply(len) - 1)

df["mean_sigma"] = df.apply(lambda r:
                            np.nanmean([r[f"sigma_{rep}"]
                                        for rep in replicate_names]),
                            axis=1)

# pin boundary Kds
df["log10Kd_pinned"] = df.log10Kd.apply(lambda x: max(
    min(x, conf['max_log10Kd']), conf['min_log10Kd']))

# drop useless columns
df = df.drop(["all_Kds"], axis=1)

# save tsv
df.to_csv(snakemake.output.tsv, sep="\t")

# plot the comparison between replicate for log10Kd
try:
    rep1, rep2 = [(a, b)
                  for a in replicate_names for b in replicate_names if a != b][0]
except IndexError:
    # if there's only one replicate, still plot the "jointplot",
    # but they won't be very informative
    rep1, rep2 = replicate_names[0], replicate_names[0]

if conf['do_plots']:
    sns.jointplot(x=f"log10Kd_{rep1}", y=f"log10Kd_{rep2}",
                  data=df, kind="hex", bins='log')
    plt.title(f"$\log_{{10}}(Kd)$ of {rep1} vs {rep2}")
else:
    plt.subplots()
plt.savefig(snakemake.output.plot_replicate_comparison)

# plot the comparison between Kd error inter and intra replicates

if conf['do_plots'] and rep1 != rep2:
    sns.jointplot(x=f"mean_sigma", y=f"err_log10Kd",
                  data=df, kind="hex", bins='log').set_axis_labels(
                      'Mean error on Kd for one replicate', 'Error between replicates')
    plt.title(f"Comparison between intra and inter-replicates errors")
else:
    plt.subplots()
plt.savefig(snakemake.output.plot_comparison_errors)
