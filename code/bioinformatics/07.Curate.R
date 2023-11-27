
# Script: Post-bioinformatics curation of OTUs and taxa

# (1) Filter fungi using taxon specific similarity thresholds
# (2) Remove low abundance OTUs
# (3) Subset arbuscular mycorrhizal OTU and taxa tables
# (4) Subset ectomycorrhizal OTU and taxa tables
# (5) Assess other dominant guilds

# Required packages
source("code/statistics/functions.R")
require(tidyverse)

#### (1) Filter fungi using taxon specific similarity thresholds ##############

# Filter taxa at according to similarity thresholds defined by Tedersoo et al.
# (2022). Best practices in metabarcoding of fungi. Molecular Ecology, 31(10),
#2769-2795.

# Similarity filters:
# Most fungi can be filtered at level phylum, however, there are exceptions for
# ectomycorrhizal lineages, as well as some Order and Class. I will filter from
# low to high levels and then rejoin the data.
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

# Number of quality filtered reads and OTUs annotated to guilds
sum(otu) # Number of reads after quality filtering = 6765543
nrow(otu) # Number of OTUs after quality filtering = 4105

#### (2) Remove low abundance OTUs #############################################

otu1 <- filter_low_abundance_otus(otu, threshold = 0.1)

# Number of reads and OTUs after within sample low abundance filter
sum(otu1) # Reads after low abundance filtering = 6589721
nrow(otu1) # OTUs after low abundance filtering = 469

# Filter taxa to OTUs
taxa2 <- taxa1 %>%
  rownames_to_column(var = "OTU_ID") %>%
  filter(OTU_ID %in% rownames(otu1)) %>%
  glimpse()
nrow(taxa2)

# Write OTU and taxa tables
taxa2 %>%
  write_csv('data/statistics/taxa_fungi.csv')
otu1 %>%
  rownames_to_column(var = 'OTU_ID') %>%
  write_csv('data/statistics/otu_fungi.csv')

#### (3) Subset arbuscular mycorrhizal OTU and taxa tables #####################

taxaAM <- taxa2 %>%
  filter(phylum == 'Glomeromycota') %>%
  glimpse()
range(taxaAM$similarity)
range(taxaAM$coverage)

otuAM <- otu %>%
  rownames_to_column(var = "OTU_ID") %>%
  filter(OTU_ID %in% taxaAM$OTU_ID) %>%
  column_to_rownames(var = "OTU_ID") %>%
  glimpse()

# AM quality filtered reads and OTUs
sum(otuAM) # Number of reads after quality filtering = 74600
nrow(otuAM) # Number of OTUs after quality filtering = 68
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

# There are a few randomish best hits at level genus that I will amend here. 
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

taxaEM <- taxa2 %>%
  filter(primary_lifestyle == 'ectomycorrhizal') %>%
  glimpse()

otuEM = otu %>%
  rownames_to_column(var = "OTU_ID") %>%
  filter(OTU_ID %in% taxaEM$OTU_ID) %>%
  column_to_rownames(var = "OTU_ID") %>%
  glimpse()
  
# EM quality filtered reads and OTUs
sum(otuEM) # Number of reads before quality filtering = 3012279
nrow(otuEM) # Number of OTUs before quality filtering = 80

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

# Sapratrpohs
taxaSAP <- taxa2 %>%
  filter(
    grepl('soil', primary_lifestyle) |
      grepl('litter', primary_lifestyle) |
      grepl('wood', primary_lifestyle)) %>%
  glimpse()
unique(taxaSAP$primary_lifestyle)

otuSAP = otu1 %>%
  rownames_to_column(var = "OTU_ID") %>%
  filter(OTU_ID %in% taxaSAP$OTU_ID) %>%
  column_to_rownames(var = "OTU_ID") %>%
  glimpse()

# EM quality filtered reads and OTUs
sum(otuSAP) # Number of reads before quality filtering = 2569111
nrow(otuSAP) # Number of OTUs before quality filtering = 169

# Plant pathogens
# Sapratrpohs
taxaPAT <- taxa2 %>%
  filter(grepl('pathogen', primary_lifestyle)) %>%
  glimpse()
unique(taxaPAT$primary_lifestyle)

otuPAT = otu1 %>%
  rownames_to_column(var = "OTU_ID") %>%
  filter(OTU_ID %in% taxaPAT$OTU_ID) %>%
  column_to_rownames(var = "OTU_ID") %>%
  glimpse()

# EM quality filtered reads and OTUs
sum(otuPAT) # Number of reads before quality filtering = 473722
nrow(otuPAT) # Number of OTUs before quality filtering = 55
