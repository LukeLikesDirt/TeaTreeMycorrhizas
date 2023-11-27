#!/bin/bash

## Script: Remove chimeras using de novo and reference-based chimera detection in VSEARCH.
## Credit: This script is adapted from https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline
## Author: Luke Florence.
## Date: 4th November 2023.
##
## Software:
## --------
## VSEARCH v2.22.1: https://github.com/torognes/vsearch
## perl v5.32.1: https://www.perl.org/get.html
##
## Script Overview:
## ---------------
##   (1) Dereplicate across samples.
##   (2) De novo chimera detection.
##   (3) Reference-based chimera detection.
##   (4) Extract all dereplicated non-chimeric sequences.
##
## Pre-requisites:
## ---------------
## By this point reads should be quality processed and dereplicated within
## samples.
##
## There should be one file named 'all.fasta' formatted for VSEARCH. That is,
## each unique sequence should occupy two lines, formatted as follows:
##  (1) <sample_name>.fasta.<unique_ID>;size=<number_of_reads>;
##  (2) <the actual sequence>
##
## There should also be a mapping file named 'map.pl'. The mapping file can be
## found in this repository or on the VSEARCH wiki linked at the top of the 
## page. Otherwise here is a link to wiki with an example VSEARCH pipeline that
## doesn't use a map.pl file: https://github.com/torognes/vsearch/wiki/Alternative-VSEARCH-pipeline

## Constants and subdirectories
readonly THREADS=8                                                                               ## Set the number of threads
readonly IDENTITY=0.97                                                                           ## Set the identity threshold for preclustering
readonly MAP_SCRIPT="../../code/bioinformatics/map.pl"                                   ## Map file for fasta reconstruction
readonly REFERENCE_SEQS="../../data/bioinformatics/07.Reference_dataset/ITS2/ITS2.fasta" ## Path to UNITE reference dataset                          
readonly DENOISED_DIR="../../data/bioinformatics/06.Denoised_uchime3"                    ## Path to denoised fasta file
readonly CHIMERA_FILTERED_DIR="../../data/bioinformatics/08.Chimeras_filtered_uchime3"   ## Path for chimera filtered fasta file

## Create subdirectory
mkdir -p "$CHIMERA_FILTERED_DIR"                                                  

## Function for dereplication across samples, de novo and referenced-based chimera filtering, and mapping reads back to dereplicated non-chimeric sequences
chimera_filter() {

    # Remove any downstream files if re-running the analysis

    local files_to_remove=("all.denovo.nonchimeras.fasta" "all.derep.fasta" "all.derep.uc" "all.fasta" "all.nonchimeras.derep.fasta" "all.nonchimeras.fasta" "all.preclustered.fasta" "all.preclustered.uc" "all.ref.nonchimeras.fasta")

    for file in "${files_to_remove[@]}"; do

        rm -f "$CHIMERA_FILTERED_DIR/$file"

    done

    # Merge samples into one file
    cat "$DENOISED_DIR"/*.fasta > "$CHIMERA_FILTERED_DIR/all.fasta"

    # Dereplicate across samples

    vsearch \
        --derep_fulllength "$CHIMERA_FILTERED_DIR/all.fasta" \
        --sizein --sizeout \
        --fasta_width 0 \
        --uc "$CHIMERA_FILTERED_DIR/all.derep.uc" \
        --output "$CHIMERA_FILTERED_DIR/all.derep.fasta"

    # Precluster at 97% before chimera detection

    vsearch \
        --cluster_size "$CHIMERA_FILTERED_DIR/all.derep.fasta" \
        --threads "$THREADS" \
        --id "$IDENTITY" \
        --strand plus \
        --sizein --sizeout \
        --fasta_width 0 \
        --uc "$CHIMERA_FILTERED_DIR/all.preclustered.uc" \
        --centroids "$CHIMERA_FILTERED_DIR/all.preclustered.fasta"

     # De novo chimera detection
    
    vsearch \
        --uchime3_denovo "$CHIMERA_FILTERED_DIR/all.preclustered.fasta" \
        --sizein --sizeout \
        --fasta_width 0 \
        --nonchimeras "$CHIMERA_FILTERED_DIR/all.denovo.nonchimeras.fasta"    
    
    # Reference-based chimera detection

    vsearch \
        --uchime_ref "$CHIMERA_FILTERED_DIR/all.denovo.nonchimeras.fasta" \
        --threads "$THREADS" \
        --db "$REFERENCE_SEQS" \
        --sizein --sizeout \
        --fasta_width 0 \
        --nonchimeras "$CHIMERA_FILTERED_DIR/all.ref.nonchimeras.fasta"

    # Extract all non-chimeric, dereplicated sequences across samples
    perl "$MAP_SCRIPT" "$CHIMERA_FILTERED_DIR/all.derep.fasta" "$CHIMERA_FILTERED_DIR/all.preclustered.uc" "$CHIMERA_FILTERED_DIR/all.ref.nonchimeras.fasta" > "$CHIMERA_FILTERED_DIR/all.nonchimeras.derep.fasta"

    # Extract all non-chimeric, dereplicated sequences within samples
    perl "$MAP_SCRIPT" "$CHIMERA_FILTERED_DIR/all.fasta" "$CHIMERA_FILTERED_DIR/all.derep.uc" "$CHIMERA_FILTERED_DIR/all.nonchimeras.derep.fasta" > "$CHIMERA_FILTERED_DIR/all.nonchimeras.fasta"

}

###############################################################################
### Main script ###############################################################
###############################################################################

## Activate the conda environment
conda activate shell

## Dereplicate across samples and remove chimeras
chimera_filter

## Deactivate the conda environment
conda deactivate
