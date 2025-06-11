# Project: sRNA-seq Alignment and Analysis
# Author: 潘晓清 Xiaoqing Pan
# Date: 2025-06-11
# Description: This Snakemake workflow handles cleaning FASTA reads (from Novogene),
#              aligning them to the pabA reference, merging aligned reads
#              from replicates, and extracting read length and first nucleotide information.
# Version: 8.10.8

from os.path import join, expanduser
import glob
from datetime import datetime
import yaml # Added this import for YAML parsing

# Load configuration from config.yaml manually
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

# Define input directory from config
FASTA_DIR = config["FASTA_DIR"]

# New output directory for pabA aligned files from config
PABA_ALIGNED_DIR = config["PABA_ALIGNED_DIR"] 

# Define the pabA index prefix from config
INDEX_PABA = config["INDEX_PABA"] 

# Get list of samples from the clean FASTA files that actually exist
SAMPLES, = glob_wildcards(join(FASTA_DIR, '{sample}_1_clean.fa'))

# Define the rule all to specify all final outputs
rule all:
    input:
        expand(join(PABA_ALIGNED_DIR, "{sample}_{rep}_aligned.fa"),
               sample=SAMPLES,
               rep=[1, 2, 3]),
        expand(join(PABA_ALIGNED_DIR, "{sample}_aligned_merged.fa"),
               sample=SAMPLES),
        expand(join(PABA_ALIGNED_DIR, "{sample}_aligned_merged_length.txt"), # New target: read lengths
               sample=SAMPLES),
        expand(join(PABA_ALIGNED_DIR, "{sample}_aligned_merged_1st_nt.txt"), # New target: first nucleotide
               sample=SAMPLES)

# Rule to fix FASTA headers (remains unchanged)
rule fix_fa:
    input:
        join(FASTA_DIR, "{sample}_{rep}_clean.fa")
    output:
        join(FASTA_DIR, "{sample}_{rep}_clean_fixed.fa")
    shell:
        """
        echo "Processing {input}"
        echo "Input file exists: $(test -f {input} && echo 'yes' || echo 'no')"
        echo "Input is symlink: $(test -L {input} && echo 'yes' || echo 'no')"
        awk '{{gsub(\">@\",\">\", $0); print $0}}' $(readlink -f {input}) > {output}
        """
                
# Rule to perform Bowtie alignment to pabA, extract aligned reads, and format to FASTA
rule bowtie_align:
    input:
        fasta=join(FASTA_DIR, "{sample}_{rep}_clean_fixed.fa") 
    output:
        # Marked as temporary so Snakemake removes them after merge_aligned_fasta
        fa=temporary(join(PABA_ALIGNED_DIR, "{sample}_{rep}_aligned.fa")) 
    params:
        outdir=PABA_ALIGNED_DIR # Use the new pabA specific output directory
    shell:
        """
        mkdir -p {params.outdir}
        
        # Perform Bowtie alignment to the pabA index
        # Pipe output to samtools view to filter for aligned reads (-F 4)
        # Pipe to awk to extract only the read ID (>$1) and sequence ($10)
        bowtie \\
            -f \\
            -v 0 \\
            -p 10 \\
            --sam \\
            {INDEX_PABA} \\
            {input.fasta} \\
            | samtools view -F 4 - \\
            | awk '{{print ">"$1"\\n"$10}}' \\
            > {output.fa}
        """

# Rule to merge aligned FASTA files from different replicates for each sample
rule merge_aligned_fasta:
    input:
        # Collect all individual aligned FASTA files for a given sample
        expand(join(PABA_ALIGNED_DIR, "{{sample}}_{rep}_aligned.fa"), rep=[1, 2, 3])
    output:
        # Output the single merged FASTA file for the sample
        join(PABA_ALIGNED_DIR, "{sample}_aligned_merged.fa")
    shell:
        """
        echo "Merging aligned FASTA files for {wildcards.sample}"
        # Concatenate all input FASTA files into one merged output file
        cat {input} > {output}
        """

# New rule to extract read length and first nucleotide from merged FASTA files
rule extract_read_info:
    input:
        merged_fa=join(PABA_ALIGNED_DIR, "{sample}_aligned_merged.fa")
    output:
        length_txt=join(PABA_ALIGNED_DIR, "{sample}_aligned_merged_length.txt"),
        first_nt_txt=join(PABA_ALIGNED_DIR, "{sample}_aligned_merged_1st_nt.txt")
    shell:
        """
        echo "Extracting read length and first nucleotide for {wildcards.sample}"
        # Extract read lengths using bioawk
        bioawk -c fastx '{{print length($seq)}}' {input.merged_fa} > {output.length_txt}
        # Extract first nucleotide using bioawk
        bioawk -c fastx '{{print substr($seq, 1, 1)}}' {input.merged_fa} > {output.first_nt_txt}
        """