#!/bin/bash

## Script:  Prepare UNITE dataset for chimera detection and taxonomic assignment
## Author:  Luke Florence.
## Date:    4th October 2023.

## Software:
## --------
## ITSx v1.1.3: https://microbiology.se/software/itsx/
## HMMER v3.1b2: http://hmmer.org/

## Constants
readonly CPUS=8              ## Set the number of CPUS
readonly REGION="ITS1,ITS2"  ## ITSx: The target ITS region: SSU, ITS1, 5.8S, ITS2, LSU, all, none. Outputs only the actual ITS sequences (ITS1, ITS2) by default.
readonly TAXA="Fungi"        ## ITSx: The taxonomic group to extract: ALL, Fungi, Metazoa, Viridiplantae, Streptophyta, Rhodophyta, Protozoa, Chromista, Bacteria, Archaea, none. Extracts all taxa by default.

## File paths
readonly UNITE_DIR="../../data/bioinformatics/07.Reference_dataset"   ## Path to the UNITE fasta file
readonly ITS1_DIR="$UNITE_DIR/ITS1"                                           ## Path to the ITS1 fasta file
readonly ITS2_DIR="$UNITE_DIR/ITS2"                                           ## Path to the ITS2 fasta file

## Create the output subdirectories if they don't exist
mkdir -p "$UNITE_DIR" "$ITS1_DIR" "$ITS2_DIR"

# Validate paths to data subdirectories
if [ ! -d "$UNITE_DIR" ]; then
  echo "Error: UNITE_DIR path '$UNITE_DIR' does not exist"
  exit 1
fi
if [ ! -d "$ITS1_DIR" ]; then
  echo "Error: ITS1_DIR path '$ITS1_DIR' does not exist"
  exit 1
fi
if [ ! -d "$ITS2_DIR" ]; then
  echo "Error: ITS2_DIR path '$ITS2_DIR' does not exist"
  exit 1
fi

## UNITE dataset information: 
## CHANGE ME to the URL and file name for the desired UNITE version: https://unite.ut.ee/repository.php
## The URL to the latest release (25.07.2023) of UNITE general release in dynamic files for fungi with singletons set as RefS
UNITE_URL="https://files.plutof.ut.ee/public/orig/86/E8/86E80475EDB915AC7173E82787BC0B73463A7690C30A91D715CF9BA0D51059BD.tgz"
## UNITE dataset file names: 
##    - The stable version ends in '.fasta' at the development version ends in '_dev.fasta'.
##    - The 'prepare_UNITE_dataset' function uses the stable version by default.
##    - Amend "$UNITE_version" from '.fatsa' to '_dev.fasta' if you want to use the development version.
UNITE_DATASET="$UNITE_DIR/sh_general_release_dynamic_25.07.2023"
UNITE_VERSION=".fasta"
UNITE_REFORMATTED_DATASET="$UNITE_DIR/UNITE_reformatted.fasta"

## Log function
log() {
    local timestamp
    timestamp=$(date)
    printf "\n%s %s\n\n" "$1" "$timestamp"
}

## Function to prepare the UNITE dataset for ITS extraction
prepare_UNITE_dataset() {

    ## Retrieve UNITE
    log 'Downloading UNITE at:'
    wget -O "$UNITE_DIR/UNITE_compressed.tgz" "$UNITE_URL"

    ## Decompress UNITE
    log 'Decompressing UNITE at:'
    tar -zxvf "$UNITE_DIR/UNITE_compressed.tgz" -C "$UNITE_DIR/"

    ## Reformat to remove lowercase and spaces
    ## Also replace "|k__" with ";k__" to help with taxa table formating in the "07.Assign_taxonomy.sh" script
    log 'Reformatting UNITE at:'
    awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' \
       "$UNITE_DATASET$UNITE_VERSION" | tr -d ' ' | sed 's/|k__/;k__/g' > \
       "$UNITE_REFORMATTED_DATASET"
    # Remove the intermediate files
    rm "$UNITE_DIR/UNITE_compressed.tgz"
    rm ""$UNITE_DATASET"_dev.fasta"
    rm ""$UNITE_DATASET".fasta"
    # Rename the reformatted UNITE file
    mv "$UNITE_REFORMATTED_DATASET" "$UNITE_DATASET$UNITE_VERSION"

    ## Extract the ITS1 and ITS2 subregions from the UNITE dataset
    log 'Extracting the ITS region from UNITE at:'

    ITSx \
      -i "$UNITE_DATASET$UNITE_VERSION" \
      --complement T \
      --save_regions "$REGION" \
      --graphical F \
      --positions T \
      -E 1e-1 \
      -t "$TAXA" \
      --cpu "$CPUS" \
      --preserve T \
      -o "$UNITE_DIR/ITSx"

    ## Build UNITE ITS1 and ITS2 datasets for input into BLAST
    log 'Building blast databases for the ITS1 and ITS2 subregions at:'

    mv "$UNITE_DIR/ITSx.ITS1.fasta" "$ITS1_DIR/ITS1.fasta"
    mv "$UNITE_DIR/ITSx.ITS2.fasta" "$ITS2_DIR/ITS2.fasta"

    makeblastdb \
       -in "$ITS1_DIR/ITS1.fasta" \
       -out "$ITS1_DIR/ITS1" \
       -dbtype 'nucl' \
       -hash_index

    makeblastdb \
       -in "$ITS2_DIR/ITS2.fasta" \
       -out "$ITS2_DIR/ITS2" \
       -dbtype 'nucl' \
       -hash_index

}

###############################################################################
### Main script ###############################################################
###############################################################################

log 'Starting at:'

## Activate the conda environment
conda activate shell

## Execute the function
prepare_UNITE_dataset

## Deactivate the conda environment
conda deactivate

log 'Finished at:'
