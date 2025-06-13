# sRNA-seq-analysis

not untill I have to analyze sRNA-seq data for a trillion times that I finally started to learn **Snakemake** language (with the help of [**Gemini**](https://gemini.google.com/app) of course)... better late than never :)


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
