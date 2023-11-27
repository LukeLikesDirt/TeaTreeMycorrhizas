
# (1) Filter fungi using taxon specific similarity thresholds
# (2) Remove low abundance OTUs
# (3) Subset arbuscular mycorrhizal OTU and taxa tables
# (4) Subset ectomycorrhizal OTU and taxa tables

# Required packages
source("code/statistics/functions.R")
require(tidyverse)

# (1) Filter taxa at according to similarity thresholds defined by Tedersoo et 
#     al. (2022). Best practices in metabarcoding of fungi. Molecular Ecology,
#     31(10), 2769-2795.

# Similarity filters
# Note: Most fungi can be filtered at level phylum, however, there are 
# exceptions for ectomycorrto high levels then rejoin the data.
similarity_EM <- read.csv('data/bioinformatics/similarity_ectomycorrhiza.csv')
similarity_order <- read.csv('data/bioinformatics/similarity_fungi.csv') %>%
  select(order, similarity) %>%
  filter(order != "")
similarity_class <- read.csv('data/bioinformatics/similarity_fungi.csv') %>%
  select(class, similarity) %>%
  filter(class != "")
similarity_phylum <- read.csv('data/bioinformatics/similarity_fungi.csv') %>%
  select(phylum, similarity) %>%
  group_by(phylum) %>%
  summarise(similarity = first(similarity)) %>%
  filter(phylum != "")

# Read in taxa and assign functional guilds
taxa <- read.csv('data/bioinformatics/10.Taxonomy/BLAST_best_hit.csv',
                 header = TRUE, sep = ';') %>%
  mutate_all(~gsub("^.*__", "", .)) %>%
  mutate(
    similarity = as.numeric(pident),
    coverage = (as.numeric(length) / as.numeric(slen)) * 100,
    coverage = ifelse(coverage > 100, 100, coverage),
    evalue = as.numeric(evalue),
    abundance = as.numeric(str_replace(abundance, "size=", ""))
  ) %>%
  select(
    OTU_ID, abundance, phylum, class, order, family, genus, species, evalue, similarity, coverage
  ) %>%
  left_join(read.csv("data/FungalTraits.csv", header = TRUE), by = "genus") %>%
  mutate(primary_lifestyle = case_when(
    # Most functional guild are assigned at the generic level but arbuscular 
    # mycorrhizal fungi are an obvious exception
    phylum == "Glomeromycota" ~ "arbuscular_mycorrhizal",
    TRUE ~ primary_lifestyle)
  ) %>%
  filter(coverage >= 90 & !is.na(primary_lifestyle) & primary_lifestyle != "") %>%
  glimpse()

# Filter by ectomycorrhizal thresholds
filtered_taxa_EM <- taxa %>%
  filter(ectomycorrhiza_lineage %in% similarity_EM$ectomycorrhiza_lineage) %>%
  left_join(similarity_EM, by = "ectomycorrhiza_lineage") %>%
  filter(similarity.x >= similarity.y) %>%
  select(-starts_with("similarity"))

# Filter by order level thresholds
filtered_taxa_order <- taxa %>%
  filter(
    !ectomycorrhiza_lineage %in% similarity_EM$ectomycorrhiza_lineage,
    order %in% similarity_order$order
  ) %>%
  left_join(similarity_order, by = "order") %>%
  filter(similarity.x >= similarity.y) %>%
  select(-starts_with("similarity"))

# Filter by class level thresholds
filtered_taxa_class <- taxa %>%
  filter(
    !ectomycorrhiza_lineage %in% similarity_EM$ectomycorrhiza_lineage,
    !order %in% similarity_order$order,
    class %in% similarity_class$class
  ) %>%
  left_join(similarity_class, by = "class") %>%
  filter(similarity.x >= similarity.y) %>%
  select(-starts_with("similarity"))

# Filter by phylum level thresholds
filtered_taxa_phylum <- taxa %>%
  filter(
    !ectomycorrhiza_lineage %in% similarity_EM$ectomycorrhiza_lineage,
    !order %in% similarity_order$order,
    !class %in% similarity_class$class,
    phylum %in% similarity_phylum$phylum
  ) %>%
  left_join(similarity_phylum, by = "phylum") %>%
  filter(similarity.x >= similarity.y) %>%
  select(-starts_with("similarity"))

# Join filtered taxa tables
taxa1 <- bind_rows(
  filtered_taxa_EM,
  filtered_taxa_order,
  filtered_taxa_class,
  filtered_taxa_phylum) %>%
  arrange(desc(abundance)) %>%
  left_join(taxa %>% select(OTU_ID, similarity), by = "OTU_ID") %>%
  column_to_rownames(var = "OTU_ID") %>%
  glimpse()

# Filter OTUs to taxa table
otu <- read.table('data/bioinformatics/09.Clustered/OTUs.txt',
                  header = TRUE) %>%
  filter(OTU_ID %in% rownames(taxa1)) %>%
  column_to_rownames(var = "OTU_ID") %>%
  glimpse()

# Number of quality filtered reads and OTUs annotated to fungi
sum(otu) # Number of reads after quality filtering = 6 683 265
nrow(otu) # Number of OTUs after quality filtering = 4 105

#### (2) Remove low abundance OTUs #############################################

# Here I filter at a low threshold (0.01% or 1 in 10 000) because AM fungi are
# not well represented in the dataset. I will apply a higher abundance threshold
# across AM and ECM after subsetting.

otu1 <- filter_low_abundance_otus(otu, threshold = 0.01)

# Number of reads and OTUs after within sample low abundance filter
sum(otu1) # Reads after quality and low abundance filtering = 6 642 663
nrow(otu1) # OTUs after quality and low abundance filtering = 1 233

# Filter taxa to OTUs
taxa2 <- taxa1 %>%
  rownames_to_column(var = "OTU_ID") %>%
  filter(OTU_ID %in% rownames(otu1)) %>%
  glimpse()
nrow(taxa2)

#### (3) Subset arbuscular mycorrhizal OTU and taxa tables #####################

# I assign AM fungal OTUs to Glomeromycota to using a similarity threshold of 
# 85% and coverage of 90%, as defined by Tedersoo et al. (2022). Best practices 
# in metabarcoding of fungi. Molecular Ecology, 31(10), 2769-2795.

taxaAM <- taxa2 %>%
  filter(phylum == 'Glomeromycota') %>%
  glimpse()
range(taxaAM$similarity)
range(taxaAM$coverage)

otuAM <- otu %>%
  rownames_to_column(var = "OTU_ID") %>%
  filter(OTU_ID %in% taxaAM$OTU_ID) %>%
  column_to_rownames(var = "OTU_ID") %>%
  filter_low_abundance_otus(., threshold = 0.1) %>%
  glimpse()

# AM quality filtered reads and OTUs
sum(otuAM) # Number of reads after quality filtering = 82 240
nrow(otuAM) # Number of OTUs after quality filtering = 69
nrow(taxaAM)

# Export filtered AM OTU table
otuAM %>%
  rownames_to_column(var = 'OTU_ID') %>%
  write_csv('data/statistics/otuAM.csv')

# Manually inspect best ten hits of AM fungi to ensure abundant OTUs do not have
# ambiguous assignments
read.csv('data/bioinformatics/10.Taxonomy/BLAST_best_10.csv',
         header = TRUE, sep = ';') %>%
  mutate_all(~gsub("^.*__", "", .)) %>%
  mutate(similarity = as.numeric(pident)) %>%
  mutate(coverage = (as.numeric(length) / as.numeric(slen)) * 100) %>%
  mutate(coverage = ifelse(coverage > 100, 100, coverage)) %>%
  mutate(evalue = as.numeric(evalue)) %>%
  filter(OTU_ID %in% taxaAM$OTU_ID) %>%
  mutate(abundance = as.numeric(gsub("size=", "", abundance))) %>%
  arrange(desc(abundance)) %>%
  select(
    OTU_ID, family, genus, evalue, similarity, coverage
  ) %>%
  as_tibble() %>%
  print(n = Inf)

# There are a few random best hits at level genus that I will amend here. 
# Annotations were stable at family and I should consider focusing on family.
taxaAM1 <- taxaAM %>%
  mutate(
        genus = ifelse(OTU_ID == "fc642738e2eee360753ec56354270cf5683bdd11",
                       "Rhizophagus", genus),
        genus = ifelse(OTU_ID == "92abe6bfbf4b7c5be5f2f87fc9ab227f19c2d982",
                       "Dominikia", genus),
        genus = ifelse(OTU_ID == "425c7ba80e82ba44ca5817191ecc5384c182cb62",
                       "Rhizophagus", genus),
        genus = ifelse(OTU_ID == "08122ec8192500939a6d7c9f63e059a65c7bbf60",
                       "Rhizophagus", genus),
        genus = ifelse(OTU_ID == "b64328e2eeb7f32582feac638f074f3eda035d82",
                       "Dominikia", genus),
        genus = ifelse(OTU_ID == "a50bb24649a17886724f2ea8cf997f58a8ed0364",
                       "Dominikia", genus),
        genus = ifelse(OTU_ID == "5923a04272abb45d4227c8cc3247afc564801586",
                       "Glomus", genus),
        genus = ifelse(OTU_ID == "fb8976d6bb5be40a9def5e886b86e1cd5fe773fe",
                       "Rhizophagus", genus),
        genus = ifelse(OTU_ID == "c1e2c7ffc6bb6ac7dddc1355ea163a52b822b5a2",
                       "Rhizophagus", genus),
        genus = ifelse(OTU_ID == "61a0d4a3944e325488fdfc5ab44e3b247e28541d",
                       "Glomus", genus),
        class = ifelse(genus == "Rhizophagus",
                       "Glomeromycetes", class),
        order = ifelse(genus == "Rhizophagus",
                       "Glomerales", order),
        family = ifelse(genus == "Rhizophagus",
                        "Glomeraceae", family),
        genus = ifelse(genus == "Claroideoglomeraceae_gen_Incertae_sedis",
                       "Claroideoglomus", genus)
      ) %>%
      write_csv('data/statistics/taxaAM.csv')

#### (4) Subset ectomycorrhizal OTU and taxa tables ############################

# OTUs are assigned to ECM fungi based on FungalTraits of genera (Põlme et 
# al. (2021). FungalTraits: a user-friendly traits database of  fungi and 
# fungus-like stramenopiles. Fungal Diversity, 105, 1–16)

taxaEM <- taxa2 %>%
  filter(primary_lifestyle == 'ectomycorrhizal') %>%
  glimpse()

otuEM = otu1 %>%
  rownames_to_column(var = "OTU_ID") %>%
  filter(OTU_ID %in% taxaEM$OTU_ID) %>%
  column_to_rownames(var = "OTU_ID") %>%
  filter_low_abundance_otus(., threshold = 0.1) %>%
  glimpse()
  
# EM quality filtered reads and OTUs
sum(otuEM) # Number of reads before quality filtering = 3 012 772
nrow(otuEM) # Number of OTUs before quality filtering = 112

# Export quality filtered EM OTU table
otuEM %>%
  rownames_to_column(var = 'OTU_ID') %>%
  write_csv('data/statistics/otuEM.csv')
# Manually inspect best ten hits of EM fungi to ensure abundant OTUs do not have
# ambiguous assignments
read.csv('data/bioinformatics/10.Taxonomy/BLAST_best_10.csv',
         header = TRUE, sep = ';') %>%
  mutate_all(~gsub("^.*__", "", .)) %>%
  mutate(similarity = as.numeric(pident)) %>%
  mutate(coverage = (as.numeric(length) / as.numeric(slen)) * 100) %>%
  mutate(coverage = ifelse(coverage > 100, 100, coverage)) %>%
  mutate(evalue = as.numeric(evalue)) %>%
  filter(OTU_ID %in% rownames(taxaEM)) %>%
  mutate(abundance = as.numeric(gsub("size=", "", abundance))) %>%
  select(
    OTU_ID, family, genus, evalue, similarity, coverage
  ) %>%
  as_tibble() %>%
  print(n = 999)
# Looks good
# Export filtered EM taxa table
taxaEM %>%
  write_csv('data/statistics/taxaEM.csv')

#### (5) Assess other guilds ##################################################

OTU_ID <- taxa2 %>%
  select(OTU_ID, primary_lifestyle) %>%
  inner_join(otu1 %>% rownames_to_column(var = "OTU_ID"), by = "OTU_ID") %>%
  pivot_longer(cols = -c(OTU_ID, primary_lifestyle), names_to = "sample") %>%
  select(OTU_ID)

otu_fungi <- taxa2 %>%
  select(OTU_ID, primary_lifestyle) %>%
  inner_join(otu1 %>% rownames_to_column(var = "OTU_ID"), by = "OTU_ID") %>%
  select(-OTU_ID) %>%
  pivot_longer(-primary_lifestyle, names_to = "sample") %>%
  group_by(primary_lifestyle) %>%
  mutate(rel_abund = 100 * value / sum(value)) %>%
  ungroup() %>%
  mutate(value = replace(value, rel_abund < 0.1, 0)) %>%
  bind_cols(OTU_ID) %>%
  select(OTU_ID, sample, value) %>%
  pivot_wider(id_cols = OTU_ID, names_from = "sample", values_from = "value") %>%
  column_to_rownames(var = "OTU_ID") %>%
  filter(rowSums(.) != 0)

taxa_fungi <- taxa2 %>%
  filter(OTU_ID %in% rownames(otu_fungi))

sum(otu_fungi) # 6 346 897
nrow(otu_fungi) # 401

# Write OTU and taxa tables
taxa_fungi %>%
  write_csv('data/statistics/taxa_fungi.csv')
otu_fungi %>%
  rownames_to_column(var = 'OTU_ID') %>%
  write_csv('data/statistics/otu_fungi.csv')

# Sapratrpohs

# Plant pathogens