# sRNA-seq-analysis

Not untill I have to analyze sRNA-seq data for a trillion times that I finally started to learn **Snakemake** language (with the help of [**Gemini**](https://gemini.google.com/app) of course)... better late than never :)

This is a simply workflow for the **clean data** from Novogene as an input. Given that Novogene gives out clean fasta files has no quality, the first thing that needs to be done is to fix the fasta files t, then use **bowtie** to align to the genome (attached), then extract the length and first nucleotide of each read to visualize.

## Install [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html)

Snakemake can be installed with all goodies needed to run in any environment and for creating interactive reports via

```{bash}
conda create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
snakemake --help
snakemake --version
```

## Run Snakefile

**Snakefile** should NOT be changed, you can change the folder path (it should be the full path) in config.yaml. Three fasta files are shown in test/clean_dir

```{bash}
# Get today's date in YYYYMMDD format
TODAY=$(date +%Y%m%d)

# Run Snakemake, specifying the config file and redirecting all output to a log file
snakemake \
-j 12 \
--configfile config.yaml > "${TODAY}_snakemake.log" 2>&1
```

The output should be three files: CEA17_aligned_merged.fa, CEA17_aligned_merged_1st_nt.txt, CEA17_aligned_merged_length.txt

## Visualization via R
```{bash}
# change the "aligned_dir" in *visualization.R* if necessary
# make sure you have the required R packages
Rscript visualization.R
```
Output should look like this:

![get_1st_nt](../test/output/1st_nt.pdf)
![length_distribution](../test/ouptut/length_distribution.pdf)
