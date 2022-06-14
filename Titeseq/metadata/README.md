# Settings and sequences

There are three metadata files:

## `config.yaml`

These are the various settings for snakemake. Detailed information is included as comments in the file. These include

* the sequences, primers, number of mutations, error tolerance, etc.
* Locations of various files:
    * a directory for all the NGS sequencing data
        * Separate sequencing runs can be placed under different subdirectories
    * a directory for all the cytometry data
        * Separate Titeseq runs can be placed under different subdirectories
    * locations for the other two metadata files, `sample_info.tsv` and `replicate_info.tsv`

## `sample_info.tsv`

- `concentration`: Concentration of the solute, string used to label the tubes
- `concentration_float`: Measured concentration of the solute
- `bin`: bin name used to label the tubes

TODO finish this

## `replicate_info.tsv`

TODO finish this
