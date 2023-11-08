
# (1) Filter taxa to fungi and e-value threshold of e-50
# (2) Remove low abundance OTUs
# (3) Subset arbuscular mycorrhizal OTU and taxa tables
# (4) Subset ectomycorrhizal OTU and taxa tables

# Required packages
require(tidyverse)

# (1) Quality filter taxa at level fungi using e-value threshold of e-50 and
#     removing OTUs that are annotated to kingdoms other than Fungi

# Filter to fungi
taxa <- read.csv('data/bioinformatics/10.Taxonomy_original/BLAST_best_hit.csv',
                  header = TRUE, sep = ';', row.names = 'OTU_ID') %>%
  mutate_all(~gsub("^.*__", "", .)) %>%
  mutate(similarity = as.numeric(pident)) %>%
  mutate(coverage = (as.numeric(length) / as.numeric(slen)) * 100) %>%
  mutate(coverage = ifelse(coverage > 100, 100, coverage)) %>%
  filter(
    kingdom == 'Fungi' | similarity >= 85 | coverage >= 90 | evalue <= 1.00e-50
  ) %>%
  select(
    phylum, class, order, family, genus, species, evalue, similarity, coverage
  ) %>%
  glimpse()

taxa_denoised <- read.csv('data/bioinformatics/09.Taxonomy_uchime3/BLAST_best_hit.csv',
                 header = TRUE, sep = ';', row.names = 'OTU_ID') %>%
  mutate_all(~gsub("^.*__", "", .)) %>%
  mutate(similarity = as.numeric(pident)) %>%
  mutate(coverage = (as.numeric(length) / as.numeric(slen)) * 100) %>%
  mutate(coverage = ifelse(coverage > 100, 100, coverage)) %>%
  filter(
    kingdom == 'Fungi' | similarity >= 85 | coverage >= 90 | evalue <= 1.00e-50
  ) %>%
  select(
    phylum, class, order, family, genus, species, evalue, similarity, coverage
  ) %>%
  glimpse()

# Filter OTUs to taxa table
otu <- read.table('data/bioinformatics/09.Clustered_original/OTUs.txt',
                  header = TRUE) %>%
  filter(OTU_ID %in% rownames(taxa)) %>%
  column_to_rownames(var = "OTU_ID") %>%
  glimpse()

otu_denoised <- read.table('data/bioinformatics/08.Clustered_uchime3/OTUs.txt',
                  header = TRUE) %>%
  filter(OTU_ID %in% rownames(taxa_denoised)) %>%
  column_to_rownames(var = "OTU_ID") %>%
  glimpse()

# Number of quality filtered reads and OTUs annotated to fungi
sum(otu) # Number of reads after quality filtering = 8822879
nrow(otu) # Number of OTUs after quality filtering = 14968

sum(otu_denoised) # Number of reads after quality filtering = 8871022
nrow(otu_denoised) # Number of OTUs after quality filtering = 31098


# (2) Remove low abundance OTUs: filter low abundance OTUs that are likely
#     the result of errors, such as index switching, using a within sample
#     abundance threshold of 0.01%.

otu1 <- otu %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample') %>%
  as_tibble() %>%
  pivot_longer(-sample) %>%
  # Group by sample to account for differences in sequencing depth
  group_by(sample) %>%
  # Calculate relative abundance within samples
  mutate(rel_abund = 100 * value / sum(value)) %>%
  ungroup() %>%
  # Remove OTUs that have relative abundance < 0.01 % within individual samples
  mutate(count = replace(value, rel_abund < 0.01, 0)) %>%
  select(-c(value, rel_abund)) %>%
  pivot_wider(values_from = 'count') %>%
  column_to_rownames(var = 'sample') %>%
  t() %>%
  as.data.frame() %>%
  filter(rowSums(.) != 0) %>%
  glimpse()

otu_denoised1 <- otu_denoised %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample') %>%
  as_tibble() %>%
  pivot_longer(-sample) %>%
  # Group by sample to account for differences in sequencing depth
  group_by(sample) %>%
  # Calculate relative abundance within samples
  mutate(rel_abund = 100 * value / sum(value)) %>%
  ungroup() %>%
  # Remove OTUs that have relative abundance < 0.01 % within individual samples
  mutate(count = replace(value, rel_abund < 0.01, 0)) %>%
  select(-c(value, rel_abund)) %>%
  pivot_wider(values_from = 'count') %>%
  column_to_rownames(var = 'sample') %>%
  t() %>%
  as.data.frame() %>%
  filter(rowSums(.) != 0) %>%
  glimpse()

# Number of reads and OTUs after low abundance filter
sum(otu1) # Number of reads after quality filtering = 8716370
nrow(otu1) # Number of OTUs before quality filtering = 2283

sum(otu_denoised1) # Number of reads after quality filtering = 8751245
nrow(otu_denoised1) # Number of OTUs before quality filtering = 1972


taxa1 <- taxa %>%
  rownames_to_column(var = "OTU_ID") %>%
  filter(OTU_ID %in% rownames(otu1)) %>%
  column_to_rownames(var = "OTU_ID") %>%
  glimpse()
nrow(taxa1)

taxa_denoised1 <- taxa_denoised %>%
  rownames_to_column(var = "OTU_ID") %>%
  filter(OTU_ID %in% rownames(otu_denoised1)) %>%
  column_to_rownames(var = "OTU_ID") %>%
  glimpse()
nrow(taxa_denoised1)

# Write filtered tables to csv
taxa1 %>%
  rownames_to_column(var = 'OTU_ID') %>%
  write_csv('data/statistics/taxa_fungi.csv')
otu1 %>%
  rownames_to_column(var = 'OTU_ID') %>%
  write_csv('data/statistics/otu_fungi.csv')
taxa_denoised1 %>%
  rownames_to_column(var = 'OTU_ID') %>%
  write_csv('data/statistics/taxa_fungi_denoised.csv')
otu_denoised1 %>%
  rownames_to_column(var = 'OTU_ID') %>%
  write_csv('data/statistics/otu_fungi_denoised.csv')

# (3) Subset arbuscular mycorrhizal OTU and taxa tables to Glomeromycota
#     using e-value threshold of e-55 and to genus using a similarity threshold
#     of 85%: These thresholds are defined by Tedersoo et al. (2022). Best
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
range(taxaAM$evalue)
otuAM <- otu %>%
  rownames_to_column(var = "OTU_ID") %>%
  filter(OTU_ID %in% rownames(taxaAM)) %>%
  column_to_rownames(var = "OTU_ID") %>%
  glimpse()

# AM quality filtered reads and OTUs
sum(otuAM) # Number of reads after quality filtering = 7538
nrow(otuAM) # Number of OTUs after quality filtering = 65

taxaAM_denoised <- taxa_denoised1 %>%
  filter(phylum == 'Glomeromycota') %>%
  # Update taxonomy
  mutate(genus = ifelse(genus == "Claroideoglomeraceae_gen_Incertae_sedis", "Entrophospora", genus),
         genus = ifelse(genus == "Claroideoglomus", "Entrophospora", genus),
         family = ifelse(genus == "Entrophospora", "Entrophosporaceae", family),
         order = ifelse(genus == "Entrophospora", "Glomerales", order)) %>%
  glimpse()
unique(taxaAM_denoised$genus)

otuAM_denoised <- otu_denoised %>%
  rownames_to_column(var = "OTU_ID") %>%
  filter(OTU_ID %in% rownames(taxaAM)) %>%
  column_to_rownames(var = "OTU_ID") %>%
  glimpse()

# AM quality filtered reads and OTUs
sum(otuAM_denoised) # Number of reads after quality filtering = 7309
nrow(otuAM_denoised) # Number of OTUs after quality filtering = 46

# Export quality filtered AM OTUs and taxa tables
otuAM %>%
  rownames_to_column(var = 'OTU_ID') %>%
  write_csv('data/statistics/otuAM.csv')
taxaAM %>%
  rownames_to_column(var = 'OTU_ID') %>%
  write_csv('data/statistics/taxaAM.csv')
otuAM_denoised %>%
  rownames_to_column(var = 'OTU_ID') %>%
  write_csv('data/statistics/otuAM_denoised.csv')
taxaAM_denoised %>%
  rownames_to_column(var = 'OTU_ID') %>%
  write_csv('data/statistics/taxaAM_denoised.csv')


# (4) Subset ectomycorrhizal OTU and taxa tables based on FungalTraits of genera
#     Põlme et al. (2021). FungalTraits: a user-friendly traits database of
#     fungi and fungus-like stramenopiles. Fungal Diversity, 105, 1–16.

taxaEM <- inner_join(read.csv("data/statistics/taxa.csv",
                     header = TRUE),
                     read.csv("data/FungalTraits.csv",
                     header = TRUE),
                     by = 'genus') %>%
  column_to_rownames(var = "OTU_ID") %>%                 
  filter(primary_lifestyle == 'ectomycorrhizal') %>%
  filter(similarity >= 90 & coverage >= 90) %>%
  glimpse()

otuEM = read.csv('data/statistics/OTUs.csv', header = TRUE) %>%
  filter(OTU_ID %in% rownames(taxaEM)) %>%
  column_to_rownames(var = "OTU_ID") %>%
  glimpse()

# EM quality filtered reads and OTUs
sum(otuEM) # Number of reads before quality filtering = 2967469 - 3180169
nrow(otuEM) # Number of OTUs before quality filtering = 156 - 211

# VSEARCH maxEE 0.5
2967469/8592129*100
156/3325*100

# VSEARCH maxEE 2
3180169/9152687*100
211/3471*100

# Export quality filtered EM OTU and taxa tables
otuEM %>%
  rownames_to_column(var = 'OTU_ID') %>%
  write_csv('data/statistics/otuEM.csv')
taxaEM %>%
  rownames_to_column(var = 'OTU_ID') %>%
  write_csv('data/statistics/taxaEM.csv')
