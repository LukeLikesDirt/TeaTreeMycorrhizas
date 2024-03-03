#!/bin/bash

## Script: Quality filter and denoise Illumina paired-end reads.
## Credit: This script is adapted from https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline
## Author: Luke Florence.
## Date: 4th November 2023.
## Software: VSEARCH v2.22.1: https://github.com/torognes/vsearch

# Constants and subdirectories
readonly THREADS=8
readonly MAXEE=0.5
readonly MINLEN=80
readonly MAXN=0
readonly QMAX=41
readonly ITS_EXTRACTED_DIR="../../data/bioinformatics/04.ITS_extracted"
readonly QUAL_FILTERED_DIR="../../data/bioinformatics/05.Quality_filtered"
readonly DENOISED_DIR="../../data/bioinformatics/06.Denoised"

# Create subdirectories
mkdir -p "$QUAL_FILTERED_DIR" "$DENOISED_DIR"

# Quality filter function
quality_filter() {

    for file in "$ITS_EXTRACTED_DIR"/*.fastq.gz; do
    
        local filename=$(basename "$file")

        # Quality filter
        vsearch \
            --fastq_filter "${file}" \
            --fastq_maxee $MAXEE \
            --fastq_minlen $MINLEN \
            --fastq_maxns $MAXN \
            --fastq_qmax $QMAX \
            --fastqout "$QUAL_FILTERED_DIR/${filename%.fastq.gz}.fastq" \

    done
}

# Denoise function
denoise() {

    for file in "$QUAL_FILTERED_DIR"/*.fastq; do

        local filename=$(basename "$file")

        vsearch \
            --cluster_unoise "$file" \
            --sizeout \
            --minsize 1 \
            --fasta_width 0 \
            --relabel "$filename". \
            --centroids "$DENOISED_DIR/${filename%.fastq}.fasta" \
            --uc "$DENOISED_DIR/${filename%.fastq}.derep.uc"
    done
}

###############################################################################
### Main script ###############################################################
###############################################################################

# Activate conda environment
conda activate shell

# Run functions
quality_filter
denoise

# Deactivate conda environment
conda deactivate
