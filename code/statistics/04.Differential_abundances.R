
# Script: Deferential abundances of mycorrhizal fungal genera across tea tree
# ecotypes

# Required packages and functions
source("code/statistics/functions.R")
require(phyloseq)
require(ANCOMBC)
require(patchwork)
require(tidyverse)

# Metadata
data <- read.csv('data/statistics/metadata.csv', 
                header = TRUE, row.name = "sample",
                stringsAsFactors = TRUE) %>%
  glimpse()

# OTUs
otuAM <- read.csv('data/statistics/otuAM.csv',
                 header = T, row.names = 'OTU_ID')
otuEM <- read.csv('data/statistics/otuEM.csv',
                 header = T, row.names = 'OTU_ID')

# Taxa
taxaAM <- read.csv('data/statistics/taxaAM.csv',
                  header = TRUE, row.names = 'OTU_ID') %>%
  as.matrix()
taxaEM <- read.csv('data/statistics/taxaEM.csv',
                  header = T, row.names = 'OTU_ID') %>%
  as.matrix()

#### (1) ANCOMBC ########################################################################

##### (1a) AM Fungi #####

# Generate phyloseq object
psOTU_AM <- otu_table(otuAM, taxa_are_rows = TRUE)
psTAXA_AM <- tax_table(taxaAM)
psDATA_AM <- sample_data(data)
psAM <- phyloseq(psOTU_AM, psTAXA_AM, psDATA_AM)

# Run ANCOMBC
set.seed(1986)
outAM <- ancombc2(
  psAM,
  tax_level = 'genus',
  fix_formula = 'ecotype',
  rand_formula = NULL,
  p_adj_method = 'BH',
  pseudo = 0,
  pseudo_sens = TRUE,
  prv_cut = 0.05,
  lib_cut = 0,
  s0_perc = 0.05,
  group = NULL,
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = 0.05,
  n_cl = 1,
  verbose = FALSE,
  global = FALSE,
  pairwise = FALSE,
  dunnet = FALSE,
  trend = FALSE,
  iter_control = list(tol = 0.01, max_iter = 20, verbose = FALSE),
  em_control = list(tol = 1e-05, max_iter = 100),
  lme_control = lme4::lmerControl(),
  mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
  trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100)
)
outAM

# Summarise statistics
difabund_AM_df <- outAM$res %>%
  select(c(taxon, 
           'lfc' = lfc_ecotypeUpland,
           'se' = se_ecotypeUpland, 
           'W-stat' = W_ecotypeUpland, 
           'p' = p_ecotypeUpland, 
           'p_adj' = q_ecotypeUpland, 
           'diff' = diff_ecotypeUpland)) %>%
  arrange('genus') %>% as_tibble() %>% print(n = Inf)

# Export the statistics table
write.csv(difabund_AM_df, 'output/ANCOM_AM.csv', row.names = FALSE)
  
##### (1b) EM Fungi #####

# Generate phyloseq object: required to run ANCOM-BC
psOTU_EM <- otu_table(otuEM, taxa_are_rows = TRUE)
psTAXA_EM <- tax_table(taxaEM)
psMETA_EM <- sample_data(data)
psEM <- phyloseq(psOTU_EM, psTAXA_EM, psMETA_EM)

# Run ANCOMBC
set.seed(1986)
outEM = ancombc2(
  psEM,
  tax_level = 'genus',
  fix_formula = 'ecotype',
  rand_formula = NULL,
  p_adj_method = 'BH',
  pseudo = 0,
  pseudo_sens = TRUE,
  prv_cut = 0.05,
  lib_cut = 0,
  s0_perc = 0.05,
  group = NULL,
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = 0.05,
  n_cl = 1,
  verbose = FALSE,
  global = FALSE,
  pairwise = FALSE,
  dunnet = FALSE,
  trend = FALSE,
  iter_control = list(tol = 0.01, max_iter = 20, verbose = FALSE),
  em_control = list(tol = 1e-05, max_iter = 100),
  lme_control = lme4::lmerControl(),
  mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
  trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100)
)
outEM

# Summarise statistics
difabund_EM_df <- outEM$res %>%
  select(c(taxon, 
           'lfc' = lfc_ecotypeUpland, 
           'se' = se_ecotypeUpland, 
           'W-stat' = W_ecotypeUpland, 
           'p' = p_ecotypeUpland, 
           'p_adj' = q_ecotypeUpland, 
           'diff' = diff_ecotypeUpland)) %>%
  arrange('genus') %>% as_tibble() %>% print(n = Inf)

# Export the ANCOMBC output
write.csv(difabund_EM_df, 'output/ANCOM_EM.csv', row.names = FALSE)

##### (1c) EM Fungi: Exploration type #####

# Generate phyloseq object
taxaEM_explore <- read.csv('data/statistics/taxaEM.csv',
                          header = TRUE, row.names = 'OTU_ID') %>%
  select(genus = ectomycorrhiza_exploration_type) %>%
  as.matrix()
psTAXA_EM_explore <- tax_table(taxaEM_explore)
psEM_explore <- phyloseq(psOTU_EM, psTAXA_EM_explore, psMETA_EM)

# Run ANCOMBC
set.seed(1986)
outEM_explore <- ancombc2(
  psEM_explore,
  tax_level = 'genus',
  fix_formula = 'ecotype',
  rand_formula = NULL,
  p_adj_method = 'BH',
  pseudo = 0,
  pseudo_sens = TRUE,
  prv_cut = 0.05,
  lib_cut = 0,
  s0_perc = 0.05,
  group = NULL,
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = 0.05,
  n_cl = 1,
  verbose = FALSE,
  global = FALSE,
  pairwise = FALSE,
  dunnet = FALSE,
  trend = FALSE,
  iter_control = list(tol = 0.01, max_iter = 20, verbose = FALSE),
  em_control = list(tol = 1e-05, max_iter = 100),
  lme_control = lme4::lmerControl(),
  mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
  trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100)
)
outEM_explore

# Summarise statistics
difabund_EM_df_exploration_type <- outEM_explore$res %>%
  select(c('exploration_type' = taxon, 
           'lfc' = lfc_ecotypeUpland, 
           'se' = se_ecotypeUpland, 
           'W-stat' = W_ecotypeUpland, 
           'p' = p_ecotypeUpland, 
           'p_adj' = q_ecotypeUpland, 
           'diff' = diff_ecotypeUpland)) %>%
  arrange('exploration_type') %>% as_tibble() %>% print(n = Inf)

# Export ANCOM-BC output
write.csv(difabund_EM_df_exploration_type, 
          'output/ANCOM_exploration_type.csv', row.names = FALSE)

#### (2) Figures ###############################################################

##### (2a) Calculate relative abundances #####

# Plots will include all genera with relative abundance > 2 %

# Calculate relative abundances of AM genera
relabund_AM <- 
  (relabund(read.csv('data/statistics/otuAM.csv',
                    header = TRUE, row.names = 'OTU_ID')) * 100) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample') %>%
  as_tibble() %>%
  pivot_longer(-sample, names_to = 'OTU_ID', values_to = 'rel_abund') %>%
  inner_join(read.csv('data/statistics/taxaAM.csv', header = TRUE),
             by = 'OTU_ID') %>%
  select(sample, taxon = genus, rel_abund) %>%
  filter(!grepl("_gen_Incertae_sedis", taxon)) %>%
  group_by(taxon) %>%
  summarise(rel_abund = round(sum(rel_abund), digits = 1))

# Calculate relative abundances of EM genera
relabund_EM <- 
  (relabund(read.csv('data/statistics/otuEM.csv',
                     header = TRUE, row.names = 'OTU_ID')) * 100) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample') %>%
  as_tibble() %>%
  pivot_longer(-sample, names_to = 'OTU_ID', values_to = 'rel_abund') %>%
  inner_join(read.csv('data/statistics/taxaEM.csv', header = TRUE),
             by = 'OTU_ID') %>%
  select(sample, taxon = genus, rel_abund) %>%
  group_by(taxon)%>%
  summarise(rel_abund = round(sum(rel_abund), digits = 1))

##### (2b) Plot AM Fungi #####

daAMplot <- relabund_AM %>%
  inner_join(difabund_AM_df, by = 'taxon') %>%
  filter(rel_abund > 2) %>%
  filter(taxon != "Glomeraceae_gen_Incertae_sedis") %>%
  mutate(ecotype = case_when(
    lfc > 0 ~ 'Upland',
    lfc < 0 ~ 'Coastal')
  ) %>%
  mutate(label = case_when(
    p_adj <= 0.05 & p_adj > 0.01 ~ '*',
    p_adj <= 0.01 & p_adj > 0.001 ~ '**',
    p_adj <= 0.001 ~ '***',
    TRUE ~ ''
  )) %>%
  ggplot(aes(x = reorder(taxon, -lfc), lfc, 
             fill = ecotype,
             label = label)
  ) +
  geom_hline(
    yintercept = 0, linetype = "dashed",
    colour = 'black', linewidth = 0.5
  ) +
  geom_errorbar(
    aes(ymin = lfc - se, ymax = lfc + se),
    width = 0.35,
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(ymin = lfc - (1.96 * se), ymax = lfc + (1.96 * se)),
    width = 0,
    linewidth = 0.5,
    linetype = 'dotted'
  ) +
  geom_text(aes(y = lfc + se + 1.), size = 7) + # Adjust astrix position
  geom_point(
    size = 4,
    shape = 21) +
  scale_y_continuous(
    limits = c(-7, 4.5),
    breaks = c(-6, -4, -2, 0, 2, 4)
  ) +
  scale_fill_manual(
    limits = c('Upland', 'Coastal'),
    values = c( '#d95f02', '#1b9e77')) +
  xlab(NULL) +
  ylab(bquote(Log[2]~~fold~change)) +
  ggtitle('Arbuscular mycorrhizal') +
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = rel(1.2), colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1.2), colour = 'black'),
    text = element_text(size = 10),
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    ### Legend ###
    legend.position = c(1, 1), 
    legend.justification = c(1.05, 1.1),
    legend.title = element_blank(),
    legend.text = element_text(size = rel(0.8)),
    legend.key.size = unit(0.25, "cm"),
    legend.background = element_rect(fill = "NA",
                                     color = "black",
                                     linewidth = 0.3,
                                     linetype = 'dotted'),
    legend.margin = margin(t=-2, r=5, b=3, l=3),
    legend.key = element_blank()
  )
daAMplot

##### (2c) Plot EM Fungi #####

daEMplot <- relabund_EM %>%
  inner_join(difabund_EM_df, by = 'taxon') %>%
  filter(rel_abund > 2) %>%
  mutate(ecotype = case_when(
    lfc > 0 ~ 'Upland',
    lfc < 0 ~ 'Coastal')
  ) %>%
  mutate(label = case_when(
    p_adj <= 0.05 & p_adj > 0.01 ~ '*',
    p_adj <= 0.01 & p_adj > 0.001 ~ '**',
    p_adj <= 0.001 ~ '***',
    TRUE ~ ''
  )) %>%
  ggplot(aes(x = reorder(taxon, -lfc), lfc, 
             fill = ecotype,
             label = label)) +
  geom_hline(
    yintercept = 0, linetype = "dashed",
    colour = 'black', linewidth = 0.5
  ) +
  geom_errorbar(
    aes(ymin = lfc - se, ymax = lfc + se),
    width = 0.35,
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(ymin = lfc - (1.96 * se), ymax = lfc + (1.96 * se)),
    width = 0,
    linewidth = 0.5,
    linetype = 'dotted'
  ) +
  geom_text(aes(y = lfc + se + 1.2), size = 7) + # Adjust astrix position
  geom_point(
    size = 4,
    shape = 21) +
  scale_y_continuous(
    limits = c(-6.9, 4.5)
    ) +
  scale_fill_manual(
    limits = c('Upland', 'Coastal'),
    values = c( '#d95f02', '#1b9e77')) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle('Ectomycorrhizal') +
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1.2), colour = 'black'),
    text = element_text(size = 10),
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    legend.position = 'none'
   )
daEMplot

##### (2d) Join plots #####

DAplot = daAMplot + daEMplot
DAplot

ggsave('output/dif_abund.pdf', width = 5, height = 3.75)
ggsave('output/dif_abund.jpg', width = 5, height = 3.75)
ggsave('output/dif_abund.tiff', width = 5, height = 3.75)

