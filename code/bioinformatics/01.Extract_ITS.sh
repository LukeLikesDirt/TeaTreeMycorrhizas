#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00:00
#SBATCH --partition=week

# Script: Trim primers, quality truncate, and extract ITS region from Illumina paired-end reads.
# Purpose: Prepare Illumina paired-end reads targeting the ITS2 region for denoising with DADA2.
# Author: Luke Florence
# Date: 28th October 2023

# Software:
# --------
# Cutadapt v4.4: https://cutadapt.readthedocs.io/en/stable/
# ITSxpress v2.0.0: https://github.com/USDA-ARS-GBRU/itsxpress
# Trimmomatic v0.36: http://www.usadellab.org/cms/?page=trimmomatic

# Script Overview:
# ---------------
# This script performs the following tasks:
#   (1) Trim the primer sequences with Cutadapt
#   (2) Quality truncate reads with Trimmomatic
#   (3) Extract the ITS region with ITSxpress

# CHANGE ME: Absolute path to the project directory
readonly PROJECT_PATH="/path/to/project/directory/"

# Constants
readonly NUM_THREADS=8                      # Number of threads
readonly NUM_CORES=8                        # Number of cores
readonly FILE_EXT=".fastq.gz"               # File extension for the raw data
readonly PRIMER_FWD="GTGAATCATCGAATCTTTGAA" # Cutadapt: Forward primer ITS86F
readonly PRIMER_REV="TCCTCCGCTTATTGATATGC"  # Cutadapt: Reverse primer ITS4
readonly OVERLAP_FWD=19                     # Cutadapt: Overlap for the forward primer
readonly OVERLAP_REV=18                     # Cutadapt: Overlap for the reverse primer
readonly ERROR_RATE=0.1                     # Cutadapt: Maximum error rate
readonly WINDOW=4                           # Trimmomatic: The sliding window size for averaging quality scores
readonly QUAL=10                            # Trimmomatic: The quality threshold for sliding window trimming
readonly ITS_REGION="ITS2"                  # ITSxpress: The target ITS region to be extracted. Choose from ITS1, ITS2, All.
readonly TAXA="All"                         # ITSxpress: The target taxa to be extracted. Choose from 'Alveolata', 'Bryophyta', 'Bacillariophyta', 'Amoebozoa', 'Euglenozoa', 'Fungi', 'Chlorophyta', 'Rhodophyta', 'Phaeophyceae', 'Marchantiophyta', 'Metazoa', 'Oomycota', 'Haptophyceae', 'Raphidophyceae', ' Rhizaria', 'Synurophyceae', 'Tracheophyta', 'Eustigmatophyceae', 'All'.
readonly MIN_LEN=80                         # Multiple functions: Minimum length of the reads

# Data subdirectories
readonly RAW_DATA="$PROJECT_PATH/data/01.Raw_data"
readonly PRIMERS_TRIMMED="$PROJECT_PATH/data/02.Primers_trimmed"
readonly QUALITY_TRIMMED="$PROJECT_PATH/data/03.Quality_truncated"
readonly ITS_EXTRACTED="$PROJECT_PATH/data/04.ITS_extracted"

# Create data subdirectories if they do not exist
mkdir -p "$RAW_DATA" "$PRIMERS_TRIMMED" "$QUALITY_TRIMMED" "$ITS_EXTRACTED"

# Validate paths to data subdirectories
if [ ! -d "$RAW_DATA" ]; then
  echo "Error: RAW_DATA path '$RAW_DATA' does not exist"
  exit 1
fi
if [ ! -d "$PRIMERS_TRIMMED" ]; then
  echo "Error: PRIMERS_TRIMMED path '$PRIMERS_TRIMMED' does not exist" 
  exit 1
fi
if [ ! -d "$QUALITY_TRIMMED" ]; then
  echo "Error: QUALITY_TRIMMED path '$QUALITY_TRIMMED' does not exist"
  exit 1
fi
if [ ! -d "$ITS_EXTRACTED" ]; then
  echo "Error: ITS_EXTRACTED path '$ITS_EXTRACTED' does not exist"
  exit 1
fi

# Log function
log() {
    local timestamp
    timestamp=$(date)
    printf "\n%s %s\n\n" "$1" "$timestamp"
}

# Function to construct reverse-complement sequences sourced from https://github.com/Mycology-Microbiology-Center/GSMc/blob/main/02.Extract_ITS.sh
RC () {
  echo "$1" | tr "[ATGCUatgcuNnYyRrSsWwKkMmBbDdHhVv]" "[TACGAtacgaNnRrYySsWwMmKkVvHhDdBb]" | rev
}
# Reverse complement of the forward and reverse primer sequences
readonly PRIMER_FWD_RC=$( RC "$PRIMER_FWD") # Cutadapt: Reverse complement of the forward primer
readonly PRIMER_REV_RC=$( RC "$PRIMER_REV") # Cutadapt: Reverse complement of the reverse primer

# Function to trim the primer sequences from the forward and reverse reads
trim_primers() {

  log "Trimming primer sequences from the forward reads"

  for fwd_file in "$RAW_DATA"/*R1$FILE_EXT; do

    local filename=$(basename "$fwd_file")
    
    cutadapt \
      -a "$PRIMER_FWD;min_overlap=$OVERLAP_FWD"..."$PRIMER_REV_RC;min_overlap=$OVERLAP_REV" \
      -e "$ERROR_RATE" \
      --minimum-length "$MIN_LEN" \
      --discard-untrimmed \
      --cores "$NUM_CORES" \
      -o "$PRIMERS_TRIMMED/$filename" "$fwd_file" \
      >> "$PRIMERS_TRIMMED/log.txt"
  done

  log "Trimming primer sequences from the reverse reads"  

  for rev_file in "$RAW_DATA"/*R2$FILE_EXT; do

    local filename=$(basename "$rev_file")

    cutadapt \
      -a "$PRIMER_REV;min_overlap=$OVERLAP_REV"..."$PRIMER_FWD_RC;min_overlap=$OVERLAP_FWD" \
      -e "$ERROR_RATE" \
      --minimum-length "$MIN_LEN" \
      --discard-untrimmed \
      --cores "$NUM_CORES" \
      -o "$PRIMERS_TRIMMED/$filename" "$rev_file" \
      >> "$PRIMERS_TRIMMED/log.txt"

  done

  log "Primer trimming complete"
}

# Function to trim the quality of the forward and reverse reads
quality_trim() {

  log "Quality trimming reads"

  for fwd_file in "$PRIMERS_TRIMMED"/*R1$FILE_EXT; do

    rev_file="${fwd_file/_R1/_R2}"

    local fwd_filename=$(basename "$fwd_file")
    local rev_filename=$(basename "$rev_file")

    trimmomatic PE -threads "$NUM_THREADS" -phred33 \
    "$fwd_file" "$rev_file" \
    "$QUALITY_TRIMMED/$fwd_filename" "$QUALITY_TRIMMED"/unpaired_"$fwd_filename" \
    "$QUALITY_TRIMMED/$rev_filename" "$QUALITY_TRIMMED"/unpaired_"$rev_filename" \
    SLIDINGWINDOW:4:"$QUAL" MINLEN:"$MIN_LEN"

done

}

# Function to extract the ITS region from the forward and reverse reads
extract_ITS(){
    log "Starting ITS extraction"

    for fwd_file in "$QUALITY_TRIMMED"/*_R1$FILE_EXT; do

    rev_file="${fwd_file/_R1/_R2}"

    local rename_fwd_file=$(cut -d- -f1 <<< "$(basename "$fwd_file")")
    local rename_rev_file=$(cut -d- -f1 <<< "$(basename "$rev_file")")

    log "Extracting the ITS from $rename_files"

    itsxpress \
      --fastq "$fwd_file" \
      --fastq2 "$rev_file" \
      --region "$ITS_REGION" \
      --taxa "$TAXA" \
      --log "$ITS_EXTRACTED/logfile.txt" \
      --outfile "$ITS_EXTRACTED/$rename_fwd_file"_R1"$FILE_EXT" \
      --outfile2 "$ITS_EXTRACTED/$rename_rev_file"_R2"$FILE_EXT" \
      --threads "$NUM_THREADS"
      
    done

    log "ITS extraction complete"
}

# Execute the main script

log "Starting at:"

# Activate Conda environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh
conda activate shell

# Run the functions
trim_primers
quality_trim
extract_ITS

# Deactivate Conda environment
conda deactivate

log "Finished at:"
