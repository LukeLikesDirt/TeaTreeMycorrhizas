#!/bin/bash

## Script: Remove chimeras using de novo and reference-based chimera detection in VSEARCH.
## Purpose: To remove chimeras prior to clustering and generating the OTU table.
## Credit: This script is adapted from https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline
## Author: Luke Florence.
## Date: 4th November 2023.
## Software: VSEARCH v2.22.1: https://github.com/torognes/vsearch

# CHANGE ME: Absolute path to the project directory
readonly PROJECT_PATH="/path/to/project/directory"

# Constants and subdirectories
readonly THREADS=8
readonly MAXEE=1
readonly MINLEN=80
readonly MAXN=0
readonly QMAX=41
readonly ITS_EXTRACTED_DIR="$PROJECT_PATH/data/bioinformatics/04.ITS_extracted"
readonly QUAL_FILTERED_DIR="$PROJECT_PATH/data/bioinformatics/05.Quality_filtered"
readonly DEREPLICATED_DIR="$PROJECT_PATH/data/bioinformatics/06.Dereplicated"

# Create subdirectories
mkdir -p "$QUAL_FILTERED_DIR" "$DEREPLICATED_DIR"

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
            --fastaout "$QUAL_FILTERED_DIR/${filename%.fastq.gz}.fasta" \

    done
}

# Dereplicate within samples function
dereplicate_within_samples() {

    for file in "$QUAL_FILTERED_DIR"/*.fasta; do

        local filename=$(basename "$file")

        vsearch \
            --derep_fulllength "$file" \
            --strand plus \
            --sizeout \
            --fasta_width 0 \
            --relabel "$filename". \
            --output "$DEREPLICATED_DIR/$filename" \
            --uc "$DEREPLICATED_DIR/${filename%.fasta}.derep.uc"
    done
}

### Main script ###############################################################

# Activate conda environment
conda activate shell

# Run functions
quality_filter
dereplicate_within_samples

# Deactivate conda environment
conda deactivate
