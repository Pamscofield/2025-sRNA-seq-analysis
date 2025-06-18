# Project: IRT sRNA-seq Alignment and Analysis
# Author: 潘晓清 Xiaoqing Pan
# Date: 2025-06-18
# Description: This Snakemake workflow handles cleaning raw FASTA reads (from Novogene),
#              aligning them to both Aspergillus fumigatus and pabA references.
#              The AFU pipeline generates sorted BAMs, merged BAMs, and strand-specific BDG files.
#              The pabA pipeline generates aligned FASTA files, merged FASTA, and extracts
#              read length and first nucleotide information.
# Version: 8.10.8

from os.path import join, expanduser
import glob
from datetime import datetime
import yaml # For loading configuration from YAML file

# Load configuration from config.yaml manually
# Ensure 'config.yaml' is in the same directory as this Snakefile
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

# Define input directory for clean FASTA reads from config
FASTA_DIR = config["FASTA_DIR"]

# Output directory for Aspergillus fumigatus (AFU) alignment products from config
AFU_ALIGNED_DIR = config["AFU_ALIGNED_DIR"]

# Output directory for pabA alignment products from config
PABA_ALIGNED_DIR = config["PABA_ALIGNED_DIR"]

# Define the A. fumigatus index prefix from config
INDEX_AFU = config["INDEX_AFU"]

# Define the pabA index prefix from config
INDEX_PABA = config["INDEX_PABA"]

# Get list of samples from the clean FASTA files that actually exist
SAMPLES, = glob_wildcards(join(FASTA_DIR, '{sample}_1_clean.fa'))

# Define the rule all to specify all final outputs for both pipelines
rule all:
    input:
        # AFU Pipeline Outputs (BAM and BDG files)
        expand(join(AFU_ALIGNED_DIR, "{sample}_{rep}_sorted.bam"),
               sample=SAMPLES,
               rep=[1, 2, 3]),
        expand(join(AFU_ALIGNED_DIR, "{sample}_{rep}_sorted.bam.bai"),
               sample=SAMPLES,
               rep=[1, 2, 3]),
        expand(join(AFU_ALIGNED_DIR, "merge_bam", "{sample}_aligned_merged.bam"),
               sample=SAMPLES),
        expand(join(AFU_ALIGNED_DIR, "merge_bam", "{sample}_aligned_merged_minus.bdg"),
               sample=SAMPLES),
        expand(join(AFU_ALIGNED_DIR, "merge_bam", "{sample}_aligned_merged_plus.bdg"),
               sample=SAMPLES),

        # PABA Pipeline Outputs (FASTA and analysis files)
        expand(join(PABA_ALIGNED_DIR, "{sample}_aligned_merged.fa"),
               sample=SAMPLES),
        expand(join(PABA_ALIGNED_DIR, "{sample}_aligned_merged_length.txt"),
               sample=SAMPLES),
        expand(join(PABA_ALIGNED_DIR, "{sample}_aligned_merged_1st_nt.txt"),
               sample=SAMPLES)


# Rule to fix FASTA headers (common for both pipelines)
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
                
# --- AFU Alignment and BAM Processing Pipeline ---

# Rule to perform Bowtie alignment to A. fumigatus, convert to bam, and index
rule bowtie_align_afu:
    input:
        fasta=join(FASTA_DIR, "{sample}_{rep}_clean_fixed.fa") 
    output:
        bam=join(AFU_ALIGNED_DIR, "{sample}_{rep}_sorted.bam"),
        bai=join(AFU_ALIGNED_DIR, "{sample}_{rep}_sorted.bam.bai")
    params:
        outdir=AFU_ALIGNED_DIR # Use the AFU_ALIGNED_DIR for all outputs of this rule
    shell:
        """
        mkdir -p {params.outdir}
        
        # Define temporary file paths within the shell command for clarity and easier deletion
        TEMP_SAM={params.outdir}/{wildcards.sample}_{wildcards.rep}.sam
        TEMP_UNSORTED_BAM={params.outdir}/{wildcards.sample}_{wildcards.rep}.bam

        # 1. Perform Bowtie alignment to the A. fumigatus index
        # Output directly to a temporary SAM file
        bowtie \\
            --sam \\
            -v 0 \\
            -p 10 \\
            -f {INDEX_AFU} \\
            {input.fasta} > $TEMP_SAM &&

        # 2. Convert SAM to unsorted BAM
        samtools view -bS $TEMP_SAM > $TEMP_UNSORTED_BAM &&

        # 3. Sort the BAM file and output to the final sorted BAM file
        samtools sort $TEMP_UNSORTED_BAM -o {output.bam} &&

        # 4. Index the sorted BAM file
        samtools index {output.bam} &&

        # Remove temporary SAM and unsorted BAM files
        rm $TEMP_SAM
        rm $TEMP_UNSORTED_BAM
        """

# Rule to merge aligned BAM files from different replicates for each sample (AFU-specific)
rule merge_aligned_bam_afu:
    input:
        # Collect all individual aligned BAM files for a given sample
        expand(join(AFU_ALIGNED_DIR, "{{sample}}_{rep}_sorted.bam"), rep=[1, 2, 3])
    output:
        # Output the single merged BAM file for the sample in the new subfolder
        bam=join(AFU_ALIGNED_DIR, "merge_bam", "{sample}_aligned_merged.bam"),
        # IMPORTANT: Add the .bai index as an output here!
        bai=join(AFU_ALIGNED_DIR, "merge_bam", "{sample}_aligned_merged.bam.bai")
    params:
        outdir_merge_bam=join(AFU_ALIGNED_DIR, "merge_bam") # Specific output directory for merged BAMs
    shell:
        """
        echo "Merging aligned BAM files for {wildcards.sample}"
        mkdir -p {params.outdir_merge_bam}
        # Concatenate all input BAM files into one merged output file
        samtools merge {output.bam} {input} &&
        # Index the newly merged BAM file so it can be used by deepTools
        samtools index {output.bam}
        """

# Rule to convert merged BAM files to strand-specific BDG files (AFU-specific)
rule bam_to_bdg_afu:
    input:
        bam=join(AFU_ALIGNED_DIR, "merge_bam", "{sample}_aligned_merged.bam"),
        # Explicitly declare the index as an input, even if not directly used in command args
        bai=join(AFU_ALIGNED_DIR, "merge_bam", "{sample}_aligned_merged.bam.bai") 
    output:
        minus_bdg=join(AFU_ALIGNED_DIR, "merge_bam", "{sample}_aligned_merged_minus.bdg"),
        plus_bdg=join(AFU_ALIGNED_DIR, "merge_bam", "{sample}_aligned_merged_plus.bdg")
    params:
        outdir_merge_bam=join(AFU_ALIGNED_DIR, "merge_bam") # Ensure the output directory exists
    shell:
        """
        echo "Converting {wildcards.sample}_aligned_merged.bam to BDG files."
        mkdir -p {params.outdir_merge_bam} # Ensure output directory exists

        # bamCoverage for reverse strand (which corresponds to --filterRNAstrand forward in deepTools)
        # Note: deepTools' --filterRNAstrand forward option is for reads mapped to the reverse strand of the gene.
        # Check your library prep if your definition of forward/reverse differs.
        bamCoverage \\
            --filterRNAstrand forward \\
            --normalizeUsing CPM \\
            -b {input.bam} \\
            -o {output.minus_bdg} \\
            -bs 1 \\
            -of bedgraph \\
            -p 12 &&

        # bamCoverage for forward strand (which corresponds to --filterRNAstrand reverse in deepTools)
        # Note: deepTools' --filterRNAstrand reverse option is for reads mapped to the forward strand of the gene.
        bamCoverage \\
            --filterRNAstrand reverse \\
            --normalizeUsing CPM \\
            -b {input.bam} \\
            -o {output.plus_bdg} \\
            -bs 1 \\
            -of bedgraph \\
            -p 12
        """

# --- PABA Alignment and FASTA Processing Pipeline ---

# Rule to perform Bowtie alignment to pabA, extract aligned reads, and format to FASTA
rule bowtie_align_paba:
    input:
        fasta=join(FASTA_DIR, "{sample}_{rep}_clean_fixed.fa") 
    output:
        # Marked as temporary so Snakemake removes them after merge_aligned_fasta
        fa=temp(join(PABA_ALIGNED_DIR, "{sample}_{rep}_aligned.fa")) 
    params:
        outdir=PABA_ALIGNED_DIR # Use the PABA_ALIGNED_DIR for outputs of this rule
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

# Rule to merge aligned FASTA files from different replicates for each sample (PABA-specific)
rule merge_aligned_fasta_paba:
    input:
        # Collect all individual aligned FASTA files for a given sample from PABA_ALIGNED_DIR
        expand(join(PABA_ALIGNED_DIR, "{{sample}}_{rep}_aligned.fa"), rep=[1, 2, 3])
    output:
        # Output the single merged FASTA file for the sample to PABA_ALIGNED_DIR
        join(PABA_ALIGNED_DIR, "{sample}_aligned_merged.fa")
    shell:
        """
        echo "Merging aligned FASTA files for {wildcards.sample}"
        # Concatenate all input FASTA files into one merged output file
        cat {input} > {output}
        """

# Rule to extract read length and first nucleotide from merged FASTA files (PABA-specific)
rule extract_read_info_paba:
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
