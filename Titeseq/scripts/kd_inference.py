# importing things
import pandas as pd
import numpy as np
import scipy.optimize
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
mpl.rcParams['figure.dpi'] = 300

conf = snakemake.config


# reading in sample info
sample_info = snakemake.params.sample_info

# finding construct and replicate
construct = sample_info.construct.unique()[0]
replicate = sample_info.replicate.unique()[0]

# reading in counts
df = pd.read_csv(snakemake.input.count_table, sep="\t")
df["geno"] = df.geno.apply(lambda x: f"{int(x):0{conf['num_mutations']}d}")
# reading in fluorescence data
fluor = pd.read_csv(snakemake.input.fluorescence, sep="\t")
fluor = fluor[fluor.replicate == replicate]

def mean_expression(sample_info, fluor, df):
    sample_info = sample_info[sample_info.concentration == 'F'].copy()
    fluor = fluor[fluor.concentration == 'F'].copy()
    nb_bins = sample_info.bin.nunique()
    nb_genos = len(df)
    probas = np.zeros((nb_bins, nb_genos))
    counts = np.zeros((nb_bins, nb_genos))
    cells = np.zeros(nb_bins)
    meanfluor, stdfluor = np.zeros((2, nb_bins))

    for bb, gate in enumerate(range(1, nb_bins+1)):
        counts[bb, :] = df[f"{construct}_{replicate}_F_{gate}"]
        cells[bb] = sample_info[sample_info.bin == gate]["cell count"].iloc[0]
        meanfluor[bb] = fluor[fluor.bin == gate]["mean_log10_FITCs"].iloc[0]
        stdfluor[bb] = fluor[fluor.bin == gate]["std_log10_FITCs"].iloc[0]
    probas = counts / (counts.sum(axis=1)[:, None]) * cells[:, None]
    probas = probas / probas.sum(axis=0)[None, :]
    mean_log10_fluor = (probas * meanfluor[:, None]).sum(axis=0)
    std_log10_fluor = np.sqrt((stdfluor[:, None]**2 * probas**2
                               + meanfluor[:, None]**2 * probas**2 / counts).sum(axis=0))
    return mean_log10_fluor, std_log10_fluor

# plotting with error
# define function


def extractKd(concentrations, bins_mean, bins_std):
    """
        Arguments: concentrations and exponential mean bin numbers (bin number: 1, 2, 3, 4)
        Return: -log10(Kd), and the r^2
    """
    popt, pcov = scipy.optimize.curve_fit(sigmoid, concentrations,
                                          bins_mean,
                                          p0=[(-9), 10**(4), 10**(2)],
                                          sigma=bins_std, absolute_sigma=True,
                                          bounds=[(conf['bounds_log10Kd_min'],
                                                   conf['bounds_A_min'],
                                                   conf['bounds_B_min']),
                                                  (conf['bounds_log10Kd_max'],
                                                   conf['bounds_A_max'],
                                                   conf['bounds_B_max'])],
                                          maxfev=400000)
    return(-1*popt[0], popt[1], popt[2], 1 - np.sum((sigmoid(concentrations, *popt) - bins_mean)**2)/np.sum((bins_mean - bins_mean.mean())**2), np.sqrt(np.diag(pcov))[0])


def sigmoid(c, Kd, A, B):
    return np.log10(A * (10**c/((10**c)+(10**Kd))) + B)


def compute_Kds(sample_info, fluor, df):
    sample_info = sample_info[~sample_info.concentration_float.isna()].copy()
    nb_bins = sample_info.bin.nunique()
    nb_concs = sample_info.concentration.nunique()
    concentrations = sample_info.concentration_float.unique()
    nb_genos = len(df)
    probas = np.zeros((nb_bins, nb_genos, nb_concs))
    counts = np.zeros((nb_bins, nb_genos, nb_concs))
    cells = np.zeros((nb_bins, nb_concs))
    meanfluor, stdfluor = np.zeros((2, nb_bins, nb_concs))
    for bb, gate in enumerate(range(1, nb_bins+1)):
        for cc, conc in enumerate(sample_info.concentration.unique()):
            counts[bb, :, cc] = df[f"{construct}_{replicate}_{conc}_{gate}"]
            cells[bb, cc] = sample_info[(sample_info.concentration == conc)
                                        & (sample_info.bin == gate)]["cell count"].iloc[0]
            meanfluor[bb, cc] = fluor[(fluor.concentration == conc)
                                      & (fluor.bin == gate)]["mean_log10_PEs"].iloc[0]
            stdfluor[bb, cc] = fluor[(fluor.concentration == conc)
                                     & (fluor.bin == gate)]["std_log10_PEs"].iloc[0]

    probas = counts / (counts.sum(axis=1)[:, None, :]) * cells[:, None, :]
    probas = probas / probas.sum(axis=0)[None, :, :]
    mean_log10_fluor = (probas * meanfluor[:, None, :]).sum(axis=0)
    std_log10_fluor = np.sqrt((stdfluor[:, None, :]**2 * probas**2
                               + meanfluor[:, None, :]**2 * probas**2 / (1e-22 + counts)).sum(axis=0))

    # replace the "0" concentration by an arbitrary large value (here 20) and invert the values
    concentrations = -concentrations
    concentrations[concentrations ==
                   0] = concentrations[concentrations == 0] - 20
    # fit of Kd
    Kds, A, B, err, cov = np.zeros((5, nb_genos))
    for s in range(nb_genos):
        notnanindex = [ii for ii in range(nb_concs)
                       if not np.isnan(mean_log10_fluor[s, ii] + std_log10_fluor[s, ii])]
        if len(notnanindex) < 4:
            Kds[s], A[s], B[s], err[s], cov[s] = [np.nan]*5
        # not enough reads
        if np.sum(counts.sum(axis=0) > snakemake.config['min_number_counts']) < 4:
            Kds[s], A[s], B[s], err[s], cov[s] = [np.nan]*5
        else:
            Kds[s], A[s], B[s], err[s], cov[s] = extractKd(concentrations[notnanindex],
                                                           mean_log10_fluor[s,
                                                                            notnanindex],
                                                           std_log10_fluor[s, notnanindex])

    return Kds, A, B, err, cov, mean_log10_fluor, std_log10_fluor, concentrations


df_subset = df.sample(n=8)
xs = np.linspace(-14, -6, 100)
Kds, As, Bs, errs, covs, mlog10, slog10, concentrations = compute_Kds(
    sample_info, fluor, df_subset)

xlim = (-max(sample_info.concentration_float)-0.5,
        -min(sample_info[sample_info.concentration_float != 0].concentration_float)+0.5)
fig, ax = plt.subplots(4, 2, figsize=(10, 10), sharex=True)
for ii, ax in enumerate(ax.flatten()):
    ax.errorbar(x=concentrations,
                y=mlog10[ii, :], yerr=slog10[ii, :], label=df_subset.geno.iloc[ii])
    ax.plot(xs, sigmoid(xs, -Kds[ii], As[ii], Bs[ii]))
    ax.set_xlim(xlim)
    ax.set_xlabel("$\log_{10}(\mathrm{Concentration})$")
    ax.set_ylabel("Est. Mean Fluorescence")
plt.savefig(snakemake.output.plot_test_curve)


# Compute Kds for the full dataset
Kds, As, Bs, errs, covs, mean_log10_PE, std_log10_PE, concs = compute_Kds(
    sample_info, fluor, df)
df["log10Kd"] = Kds
df["A"] = As
df["B"] = Bs
df["r2"] = errs
df["sigma"] = covs

# also save the mean fluorescence values
for cc in range(0, mean_log10_PE.shape[1]):
    df[f"mean_log10PE{cc}"] = mean_log10_PE[:, cc-1]
    df[f"std_log10PE{cc}"] = std_log10_PE[:, cc-1]

if 'F' in sample_info.concentration.unique():
    m10s, s10s = mean_expression(sample_info, fluor, df)
    df["Mean fluorescence expression"] = m10s
    df["Std fluorescence expression"] = s10s

# save the Kds values
df[["geno", "log10Kd", "A", "B", "r2", "sigma"] +
   [f"mean_log10PE{cc}" for cc in range(0, mean_log10_PE.shape[1])] +
   [f"std_log10PE{cc}" for cc in range(0, mean_log10_PE.shape[1])] +
   (["Mean fluorescence expression", "Std fluorescence expression"]
    if ('F' in sample_info.concentration.unique()) else [])].to_csv(
    snakemake.output.tsv,
    sep="\t", index=False)


with pd.option_context('mode.use_inf_as_na', True):
    fig, ax = plt.subplots()
    if conf['do_plots']:
        # Kd distribution
        sns.histplot(x="log10Kd", data=df.dropna(subset=["log10Kd"]),
                     ax=ax, label="All intermediates")
        # No mutations
        if df.geno.apply(lambda x: '1' not in x).sum() > 0:
            ax.axvline(x=df[df.geno.apply(lambda x: '1' not in x)
                            ].log10Kd.iloc[0], label="Original", c="g")
        # All mutations
        if df.geno.apply(lambda x: '0' not in x).sum() > 0:
            ax.axvline(x=df[df.geno.apply(lambda x: '0' not in x)
                            ].log10Kd.iloc[0], label="Variant", c="r")
        ax.legend()
    plt.savefig(snakemake.output.plot_Kd_distribution)

    # Mean expression distribution
    fig, ax = plt.subplots()
    if "Mean fluorescence expression" in df and conf['do_plots']:
        sns.histplot(x="Mean fluorescence expression",
                     data=df.dropna(subset=["Mean fluorescence expression"]),
                     label="All intermediates", ax=ax)
        if df.geno.apply(lambda x: '1' not in x).sum() > 0:
            ax.axvline(x=df[df.geno.apply(lambda x: '1' not in x)]["Mean fluorescence expression"].iloc[0],
                       label="Original", c="g")
        if df.geno.apply(lambda x: '0' not in x).sum() > 0:
            ax.axvline(x=df[df.geno.apply(lambda x: '0' not in x)]["Mean fluorescence expression"].iloc[0],
                       label="Variant", c="r")
        ax.legend()
    plt.savefig(snakemake.output.plot_fluo_distribution)

    # Correlation btw expression and Kd
    if "Mean fluorescence expression" in df and conf['do_plots']:
        sns.jointplot(x="Mean fluorescence expression",
                      y="log10Kd",
                      data=df.dropna(
                          subset=["Mean fluorescence expression", "log10Kd"]),
                      kind="hex", bins='log')
    else:
         plt.subplots()
    plt.savefig(snakemake.output.plot_corr_fluo_Kd)

    # Correlation btw Kd error and Kd
    if conf['do_plots']:
        df["$\log_{10}(\mathrm{error})$"] = np.log(1 - df["r2"])
        sns.jointplot(x="$\log_{10}(\mathrm{error})$",
                      y="log10Kd",
                      data=df.dropna(
                          subset=["$\log_{10}(\mathrm{error})$", "log10Kd"]),
                      kind="hex", bins='log')
    else:
        plt.subplots()
    plt.savefig(snakemake.output.plot_corr_Kd_error)
