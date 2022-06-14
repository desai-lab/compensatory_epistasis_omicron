# Running titeseq

The code in here fits a dissociation constant for each of the variants using the sequencing and flow cytometry data.

## About the pipeline

This analysis is done with the `snakemake` pipeline described in the [Snakefile](Snakefile). To run this, first make sure you have installed [conda](https://docs.conda.io/en/latest/), and build the environment with

```
conda create -n omicron
```

Then activate conda with

```
source activate omicron
```

And load the necessary packages with
```
conda install -c conda-forge mamba
mamba env update --file environment.yaml
```

You can run the analysis with

```
snakemake -j 1
```
Where you can replace `1` by the number of cores you want to use.

## How to run the analysis

### Configure the pipeline

1. Go over the settings in the files [metadata/config.yaml](metadata/config.yaml) about how to parse and analyze the data.
2. Enter information about each sample into [metadata/sample_info.tsv](metadata/sample_info.tsv).
3. Enter information about each Tite-seq run into [metadata/replicate_info.tsv](metadata/replicate_info.tsv).

The file [metadata/README.md](metadata/README.md) contains more information about each of the fields.

### Prepare the data input

1. Enter the cell counts into [metadata/sample_info.tsv](metadata/sample_info.tsv).
2. Transfer the NGS read data to the correct input directory. 
    1. Create a parent directory for all the NGS data `sequencing_data_dir`, specified at the top of [metadata/config.yaml](metadata/config.yaml). Place each separate NGS run in a separate subdirectory. The name of the subdirectory for each Titeseq run is specified in [metadata/replicate_info.tsv](metadata/replicate_info.tsv) under the `sequencing_data_subdir` column.
    2. Make sure the files are named `{sample_name}_S{sample_id}_L00?_R?_001.fastq.gz`, where `sample_name` and `sample_id` are the columns in [metadata/sample_info.tsv](metadata/sample_info.tsv). These should be identical to the columns in the NGS sample sheet.
3. Transfer the fcs files (flow cytometry information) to the correct input directory. 
    1. Use Flow-Jo to read the fcs files. There's already a compensation applied. Select the "Sorted" tubes and export them to csv format (Export or Concatenate → CSV - Scale values, All compensated parameters). (You can also rerun the compensation with flow-Jo but that shouldn't be necessary.)S
    2. Create a parent directory for all the cytometry data `cytometry_data_dir`, specified at the top of [metadata/config]. Place each Titeseq run in a separate subdirectory. The name of the subidrectory for each run is specified in the `cytometry_data_subdir` column of [metadata/replicate_info.tsv](metadata/replicate_info.tsv).
    3. Filenames should be in the format
        1. `export_Sorted_{construct}{replicate}_{concentration}_PE{bin}` for samples testing for binding and 
        2. `export_Sorted_{construct}{replicate}_{concentration}_mycFITC_FITC{bin}` for the expression samples.


### Run snakemake

Snakemake will automatically perform the following steps:

1. Parse the fastq files into a table of the number of reads per genotype, for each the samples.
    1. Parse fastq files and separate each read into fields `col_idx, row_idx, UMI_1, UMI_2, read_1, read_2`
    2. Discard reads with incorrect inline indices, leaving fields `UMI_1, UMI_2, read_1, read_2`
    3. Parse read sequences into binary string genotype, leaving fields `UMI_1, UMI_2, geno`
    4. Count the number of unique UMIs appearing for each genotype, leaving fields `geno, count`
    5. Collect all samples into a single count table, with a column for each sample and a row for each genotype
2. Fit the Kd for each of the variants
    1. Read the csv fluorescence files, correct for negative values and extract the mean and standard variance of the log-fluorescence.
    2. Calculate the mean bin from the cell counts, reads, and flourescence values
    3. Fit the binding curve to the mean bin.

The relevant output files will be the following:

* The count table is written to [results_{experiment_name}/count_table.tsv](results/count_table.tsv).
* The inferred Kds are written to `results_{experiment_name}/Kds/Kds_{construct}{replicate}.tsv'`
* The fluorescence information is written to [results_{experiment_name}/fluorescence.tsv](results/fluorescence.tsv).
* Statistics about sequencing errors, number of duplicate UMI's, etc. are written to the directory [results/stats/](results/stats)
