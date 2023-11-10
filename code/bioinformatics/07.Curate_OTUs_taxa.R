
# (1) Filter fungi using an e-value threshold of e-50
# (2) Remove low abundance OTUs
# (3) Subset arbuscular mycorrhizal OTU and taxa tables
# (4) Subset ectomycorrhizal OTU and taxa tables

# Required packages
source("code/statistics/functions.R")
require(tidyverse)

# (1) Quality filter taxa at level fungi using e-value threshold of e-50

# Filter to fungi
taxa <- read.csv('data/bioinformatics/10.Taxonomy/BLAST_best_hit.csv',
                  header = TRUE, sep = ';', row.names = 'OTU_ID') %>%
  mutate_all(~gsub("^.*__", "", .)) %>%
  mutate(similarity = as.numeric(pident)) %>%
  mutate(coverage = (as.numeric(length) / as.numeric(slen)) * 100) %>%
  mutate(coverage = ifelse(coverage > 100, 100, coverage)) %>%
  mutate(evalue = as.numeric(evalue)) %>%
  filter(
    kingdom == 'Fungi' & similarity >= 85 & coverage >= 90 & evalue <= 1.00e-50
  ) %>%
  select(
    phylum, class, order, family, genus, species, evalue, similarity, coverage
  ) %>%
  glimpse()

# Filter OTUs to taxa table
otu <- read.table('data/bioinformatics/09.Clustered/OTUs.txt',
                  header = TRUE) %>%
  filter(OTU_ID %in% rownames(taxa)) %>%
  column_to_rownames(var = "OTU_ID") %>%
  glimpse()

# Number of quality filtered reads and OTUs annotated to fungi
sum(otu) # Number of reads after quality filtering = 8382739
nrow(otu) # Number of OTUs after quality filtering = 7204

# (2) Remove low abundance OTUs: filter low abundance OTUs that are likely
#     the result of errors, such as index switching, using a within sample
#     abundance threshold of 0.01%.

otu1 <- filter_low_abundance_otus(otu, threshold = 0.01)

# Number of reads and OTUs after low abundance filter
n_reads <- sum(otu1) # Number of reads after quality filtering = 8308200
n_OTUs <- nrow(otu1) # Number of OTUs before quality filtering = 1843

# Filter taxa to OTUs
taxa1 <- taxa %>%
  rownames_to_column(var = "OTU_ID") %>%
  filter(OTU_ID %in% rownames(otu1)) %>%
  column_to_rownames(var = "OTU_ID") %>%
  glimpse()
nrow(taxa1)

# Write filtered tables
taxa1 %>%
  rownames_to_column(var = 'OTU_ID') %>%
  write_csv('data/statistics/taxa_fungi.csv')
otu1 %>%
  rownames_to_column(var = 'OTU_ID') %>%
  write_csv('data/statistics/otu_fungi.csv')

# (3) Subset arbuscular mycorrhizal OTU and taxa tables to Glomeromycota
#     a similarity threshold of 85%: Defined by Tedersoo et al. (2022). Best
#     practices in metabarcoding of fungi. Molecular Ecology, 31(10), 2769-2795.
#     To conform with the most recent systematics of Glomeromycota, I will
#     re-annotated the genera Entrophospora (otherwise family incertae sedis,
#     Diversisporales) and Claroideoglomus (otherwise Claroideoglomeraceae,
#     Glomerales) to family Entrophosporaceae (Glomerales): Wijayawardene et al.
#     (2022). Outline of Fungi and fungus-like taxa – 2021. Mycosphere, 13(1),
#     53–453.

taxaAM <- taxa1 %>%
  filter(phylum == 'Glomeromycota') %>%
  # Update taxonomy
  mutate(genus = ifelse(genus == "Claroideoglomeraceae_gen_Incertae_sedis", "Entrophospora", genus),
         genus = ifelse(genus == "Claroideoglomus", "Entrophospora", genus),
         family = ifelse(genus == "Entrophospora", "Entrophosporaceae", family),
         order = ifelse(genus == "Entrophospora", "Glomerales", order)) %>%
  glimpse()
unique(taxaAM$genus)
range(taxa1$similarity)
range(taxa1$coverage)

otuAM <- otu %>%
  rownames_to_column(var = "OTU_ID") %>%
  filter(OTU_ID %in% rownames(taxaAM)) %>%
  column_to_rownames(var = "OTU_ID") %>%
  glimpse()

# AM quality filtered reads and OTUs
sum(otuAM) # Number of reads after quality filtering = 7979
nrow(otuAM) # Number of OTUs after quality filtering = 58
7979/n_reads * 100
58/n_OTUs * 100

# Export quality filtered AM OTUs and taxa tables
otuAM %>%
  rownames_to_column(var = 'OTU_ID') %>%
  write_csv('data/statistics/otuAM.csv')
taxaAM %>%
  rownames_to_column(var = 'OTU_ID') %>%
  write_csv('data/statistics/taxaAM.csv')

# (4) Subset ectomycorrhizal OTU and taxa tables based on FungalTraits of genera
#     Põlme et al. (2021). FungalTraits: a user-friendly traits database of
#     fungi and fungus-like stramenopiles. Fungal Diversity, 105, 1–16.

taxaEM <- taxa1 %>%
  rownames_to_column(var = "OTU_ID") %>%
  inner_join(read.csv("data/FungalTraits.csv", header = TRUE), by = 'genus') %>%              
  filter(primary_lifestyle == 'ectomycorrhizal') %>%
  filter(similarity >= 90) %>%
  column_to_rownames(var = "OTU_ID") %>%
  glimpse()

otuEM = otu %>%
  rownames_to_column(var = "OTU_ID") %>%
  filter(OTU_ID %in% rownames(taxaEM)) %>%
  column_to_rownames(var = "OTU_ID") %>%
  glimpse()

# EM quality filtered reads and OTUs
sum(otuEM) # Number of reads before quality filtering = 3017369
nrow(otuEM) # Number of OTUs before quality filtering = 143
3017369/n_reads * 100
143/n_OTUs * 100

# Export quality filtered EM OTU and taxa tables
otuEM %>%
  rownames_to_column(var = 'OTU_ID') %>%
  write_csv('data/statistics/otuEM.csv')
taxaEM %>%
  rownames_to_column(var = 'OTU_ID') %>%
  write_csv('data/statistics/taxaEM.csv')
