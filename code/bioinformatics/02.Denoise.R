#
# Creadit: This script is adapted from https://benjjneb.github.io/dada2/tutorial.html
#
# Script Purpose:
# ---------------
# This R script is designed for the quality filtering and denoising of Illumina
# paired-end amplicons targeting the ITS2 using DADA2. Subsequent chimera
# detection and removal will be done in VSEARCH. I therefore convert the DADA2
# 'rds' sequence table to a 'fasta' file formatted for VESEARCH. Reads are not
# quality-filtered in this script because I have already quality truncated and
# ITS2 extracted in Trimmomatic and ITSxpress, respectively. ITSxpress
# impliments quality filtering with maxEE = 2. If you want more stringent
# quality-filtering you add a filterAndTrim() step prior to denoising.
#
# Script Overview:
# ---------------
#   (1) Denoise and merge the R1 and R2 reads.
#   (2) Convert the merged sequence table to a fasta file for VSEARCH.
#
# Pre-processing:
# ---------------
# I'm processing reads targeting the ITS2 subregion. I have removed the primers
# with Cutadapt, quality truncated the reads with Trimmomatic and extracted the
# ITS2 using ITSxpress. I find preprocessing with Trimmomatic improves read
# recovery across that DADA2 pipeline. Preprocessing with ITSxpress also
# improves read recovery. But importantly, extracting the ITS is a necessary
# step to remove highly conserved non-informative flanking regions of the ITS,
# which will enhance clustering as these flanking regions likely contain random
# errors.

# Load required libraries:
library(dada2)
library(seqinr)
library(here)
library(tidyverse)

# Define the path to the data subdirectory containing the 16S fastq files
data_dir <- here("data/")
ITS_extracted_dir <- here("data/04.ITS_extracted")

# Create subdirectory for denoised files
dir.create(file.path(data_dir, "05.Denoised"))

# Assign the path and filenames for ITS extracted fastq.gz files.
fwd_file <- sort(list.files(ITS_extracted_dir, pattern = "R1.fastq.gz", full.names = TRUE))
rev_file <- sort(list.files(ITS_extracted_dir, pattern = "R2.fastq.gz", full.names = TRUE))

# Learn error rates for the forward and reverse reads
err_fwd <- learnErrors(fwd_file, multithread = TRUE)
err_rev <- learnErrors(rev_file, multithread = TRUE)
# Denoise the forward and reverse reads
dada_fwd <- dada(fwd_file, err = err_fwd, multithread = TRUE)
dada_rev <- dada(rev_file, err = err_rev, multithread = TRUE)
# Merge the forward and reverse reads
mergers <- mergePairs(dada_fwd, fwd_file, dada_rev, rev_file, verbose = TRUE)

# CONSTRUCT THE SEQUENCE TABLE:
seq_tab <- makeSequenceTable(mergers)
# Initialise an empty list to store the formatted data frames for each chunk
fasta_tab <- as.data.frame(seq_tab) %>%
  rownames_to_column('sample.names') %>%
  as_tibble() %>%
  pivot_longer(-sample.names, names_to = 'seq', values_to = 'size') %>%
  filter(size > 0) %>%
  ungroup() %>%
  mutate(seq.name = paste0(sample.names, '.fasta', '.',
                           row_number(), ';size=', size)) %>%
  select(seq.name, seq) %>%
  mutate(seq.name = str_replace(seq.name, "_R1", "")) %>%
  glimpse()

# Save the rds file for the merged sequence table
saveRDS(seq_tab, file.path(data_dir, '05.Denoised/all_seqtab.rds'))

# Write and save the fasta file for the merged sequence table
write.fasta(as.list(fasta_tab$seq), fasta_tab$seq.name,
            file.path(data_dir, '05.Denoised/all.fasta'),
            open = 'w', nbchar = 60, as.string = FALSE)
