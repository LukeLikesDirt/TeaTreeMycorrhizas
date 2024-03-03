#!/bin/bash

## Script: Taxonomic assignment with BLASTn
## Author: Luke Florence
## Date: 5th November 2023
## Software: BLAST v2.14.1 - https://blast.ncbi.nlm.nih.gov/Blast.cgi

## Constants and subdirectories
readonly THREADS=8
readonly REFERENCE_SEQUENCES="../../data/bioinformatics/07.Reference_dataset/ITS2/ITS2"
readonly OTU_FASTA="../../data/bioinformatics/09.Clustered/OTUs.fasta"
readonly TAXA_DIR="../../data/bioinformatics/10.Taxonomy"
mkdir -p "$TAXA_DIR"

## Log function
log() {
    local timestamp
    timestamp=$(date)
    printf "\n%s %s\n\n" "$1" "$timestamp"
}

## BLAST function
my_blast() {
    log 'BLAST best hit starting at'

    # blast ten hits
    blastn \
        -task blastn \
        -outfmt "6 qseqid sseqid pident length slen mismatch gapopen qstart qend sstart send evalue bitscore" \
        -strand both \
        -query "$OTU_FASTA" \
        -db "$REFERENCE_SEQUENCES" \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -out "$TAXA_DIR/BLAST_best_10.txt" \
        -num_threads "$THREADS"
        
}

## Reformat taxa table function
reformat_taxa_tables() {
    sed -i '1s/^/OTU_ID;abundance\treference;kingdom;phylum;class;order;family;genus;species\tpident\tlength\tslen\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n/' "$TAXA_DIR/BLAST_best_hit.txt"
    sed 's/[[:space:]]\{1,\}/;/g' "$TAXA_DIR/BLAST_best_hit.txt" > "$TAXA_DIR/BLAST_best_hit.csv"

    sed -i '1s/^/OTU_ID;abundance\treference;kingdom;phylum;class;order;family;genus;species\tpident\tlength\tslen\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n/' "$TAXA_DIR/BLAST_best_10.txt"
    sed 's/[[:space:]]\{1,\}/;/g' "$TAXA_DIR/BLAST_best_10.txt" > "$TAXA_DIR/BLAST_best_10.csv"
}

###############################################################################
## Main script ################################################################
###############################################################################
log 'Starting at:'

## Activate the conda environment
conda activate shell

## BLAST
my_blast
reformat_taxa_tables

## Deactivate the conda environment
conda deactivate

log 'Finished at:'
