#!/bin/bash

## Script: To cluster OTUs prior to taxonomic assignment.
## Credit: This script is adapted from https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline.
## Author: Luke Florence.
## Date: 5th November 2023.
## Software: VSEARCH v2.22.1 - https://github.com/torognes/vsearch

## Constants and file paths
readonly THREADS=8                                                                    ## Set the number of threads
readonly IDENTITY=0.97                                                                ## Set the identity threshold for clustering
readonly CHIMERA_FILTERED_DIR="../../data/bioinformatics/08.Chimeras_filtered_uchime3" ## Path to chimera filtered fasta file
readonly CLUSTERED_DIR="../../data/bioinformatics/09.Clustered_uchime3"               ## Path to clustered fasta file and OTU table
mkdir -p "$CLUSTERED_DIR"                                                             ## Make the clustered subdirectory


## Function for clustering OTUs and formatting the OTU table
cluster_OTUs() {

    ## Cluster at 98% and generate OTU table

    vsearch \
        --cluster_size "$CHIMERA_FILTERED_DIR/all.nonchimeras.fasta" \
        --threads "$THREADS" \
        --id "$IDENTITY" \
        --strand plus \
        --sizein --sizeout \
        --fasta_width 0 \
        --uc "$CHIMERA_FILTERED_DIR/all.clustered.uc" \
        --relabel_sha \
        --centroids "$CLUSTERED_DIR/OTUs.fasta" \
        --otutabout "$CLUSTERED_DIR/OTUs.txt"

    ## Format OTU table
    ## Rename the header in the file
    sed -i '1s/#OTU ID/OTU_ID/' "$CLUSTERED_DIR/OTUs.txt"
    ## Convert to .csv and save to the output directory
    sed -e 's/\s\+/,/g' "$CLUSTERED_DIR/OTUs.txt" > "$OUTPUT_DIR/OTUs.csv"

    printf '\nNumber of unique sequences and OTUs\n'
    printf '    Unique non-chimeric sequence: %s\n' "$(grep -c "^>" "$CHIMERA_FILTERED_DIR/all.nonchimeras.fasta")"
    printf '    Clustered OTUs: %s\n' "$(grep -c "^>" "$CLUSTERED_DIR/OTUs.fasta")"

}

# Function to perform sequence matching table for post-clustering curation with LULU
match_list() {

    vsearch \
        --usearch_global "$CLUSTERED_DIR/OTUs.fasta" \
        --db "$CLUSTERED_DIR/OTUs.fasta" \
        --self --id .84 \
        --iddef 1 --userout "$CLUSTERED_DIR/match_list.txt" \
        -userfields query+target+id \
        --maxaccepts 0 --query_cov 0.9 --maxhits 10
}


###############################################################################
## Main script ################################################################
###############################################################################

## Activate the conda environment
conda activate shell

## Cluster OTUs
cluster_OTUs

## Create a match list
match_list

## Deactivate the conda environment
conda deactivate