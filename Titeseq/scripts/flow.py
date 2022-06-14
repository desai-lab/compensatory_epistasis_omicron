#!/usr/bin/env python3

import pandas as pd
import numpy as np

sample_table = snakemake.params.sample_table.copy()
sample_table = sample_table[sample_table.concentration != "unsorted"]


def list_FITCs(csv):
    df = pd.read_csv(csv)
    try:
        return df["Comp-FITC-A"].to_list()
    except KeyError:
        return df["FITC-A"].to_list()


def list_PEs(csv):
    df = pd.read_csv(csv)
    try:
        return df["Comp-PE-A"].to_list()
    except KeyError:
        return df["PE-A"].to_list()


# load the data
sample_table["FITCs"] = sample_table.flow_file.apply(list_FITCs)
sample_table["PEs"] = sample_table.flow_file.apply(list_PEs)
sample_table["min_FITCs"] = sample_table.FITCs.apply(np.min)
sample_table["min_PEs"] = sample_table.PEs.apply(np.min)

# compute the shift needed for each experiment (same everywhere)
sample_table["experiment"] = (sample_table["construct"]
                              + sample_table["replicate"])
minFITCs = sample_table.groupby("experiment").min_FITCs.min()
minPEs = sample_table.groupby("experiment").min_PEs.min()

sample_table["shift_FITCs"] = sample_table.experiment.map(minFITCs)
sample_table["shift_PEs"] = sample_table.experiment.map(minPEs)
sample_table["log10FITCs"] = sample_table.apply(lambda r:
                                                [np.log10(
                                                    a-r["shift_FITCs"]+1
                                                ) for a in r["FITCs"]],
                                                axis=1)
sample_table["log10PEs"] = sample_table.apply(lambda r:
                                              [np.log10(
                                                  a-r["shift_PEs"]+1)
                                               for a in r["PEs"]],
                                              axis=1)

# compute the mean and standard deviation
sample_table["mean_log10_PEs"] = sample_table.log10PEs.apply(np.mean)
sample_table["std_log10_PEs"] = sample_table.log10PEs.apply(np.std)
sample_table["mean_log10_FITCs"] = sample_table.log10FITCs.apply(np.mean)
sample_table["std_log10_FITCs"] = sample_table.log10FITCs.apply(np.std)

# save the information
sample_table[["sample_id", "construct", "replicate", "concentration",
              "bin", "mean_log10_PEs", "std_log10_PEs",
              "mean_log10_FITCs", "std_log10_FITCs"]].to_csv(snakemake.output.tsv, sep="\t")
