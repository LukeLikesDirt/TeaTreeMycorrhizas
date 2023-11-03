#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --partition=day

## Script: Taxonomic assignment with BLASTn
## Author: Luke Florence.
## Date: 5th November 2023.
## Software: BLAST v2.14.1 - https://blast.ncbi.nlm.nih.gov/Blast.cgi

# CHANGE ME: Absolute path to the project directory
readonly PROJECT_DIR="/path/to/project/directory"

# Constants and subdirectories
readonly THREADS=8
readonly REFERENCE_SEQUENCES="$PROJECT_DIR/data/06.Reference_dataset/ITS2/ITS2"
readonly OTU_FASTA="$PROJECT_DIR/data/08.Clustered/OTUs.fasta"
readonly TAXA_DIR="$PROJECT_DIR/data/09.Taxonomy"
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

    # Best hit
    blastn \
        -task blastn \
        -outfmt "6 qseqid sseqid pident length slen mismatch gapopen qstart qend sstart send evalue bitscore" \
        -query "$OTU_FASTA" \
        -db "$REFERENCE_SEQUENCES" \
        -strand both \
        -evalue 0.001 \
        -word_size 7 \
        -reward 1 \
        -penalty -1 \
        -gapopen 1 \
        -gapextend 2 \
        -max_target_seqs 1 \
        -max_hsps 1 \
        -out "$TAXA_DIR/BLAST_best_hit.txt" \
        -num_threads "$THREADS"

    log 'BLAST best ten hits starting at'

    # Best 10
    blastn \
        -task blastn \
        -outfmt "6 qseqid sseqid pident length slen mismatch gapopen qstart qend sstart send evalue bitscore" \
        -query "$OTU_FASTA" \
        -db "$REFERENCE_SEQUENCES" \
        -strand both \
        -evalue 0.001 \
        -word_size 7 \
        -reward 1 \
        -penalty -1 \
        -gapopen 1 \
        -gapextend 2 \
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

## Main script

log 'Starting at:'

## Activate the conda environment
conda activate shell

## BLAST
my_blast
reformat_taxa_tables

## Deactivate the conda environment
conda deactivate

log 'Finished at:'
