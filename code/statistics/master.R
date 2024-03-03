
# Script: Assess differences in mycorrhizal root colonisation across tea tree
# ecotypes

# Required packages and functions
require(lme4)
require(emmeans)
require(parameters)
require(performance)
require(DHARMa)
require(adespatial)
require(vegan)
require(phyloseq)
require(ANCOMBC)
require(cowplot)
require(writexl)
require(tidyverse)
source("code/statistics/functions.R")

### (1) Organise data ##########################################################

# Read in the metadata and colonisation measurements
data <- read.csv('data/statistics/metadata.csv', stringsAsFactors = TRUE,
                header = TRUE) %>%
  # Order levels
  mutate(ecotype = factor(ecotype, levels = c(
  "Upland", "Coastal"))
  ) %>%
  glimpse(.)

# Ensure that the 'group' is a factor
data$treatment <- factor(data$ecotype)

# Check the unique levels of the 'group'
unique_levels <- levels(data$ecotype)

# Set sum contrasts
contrasts(data$ecotype) <- contr.treatment(unique_levels)

# Grab the geo-coordinates for spatial autocorrelation analyses
coords = data %>%
  select(x = longitude, y = latitude)

# Calculate candidate weighing matrices for spatial autocorrelation analyses
candidates = listw.candidates(coords,
                              nb = c('del', 'rel', 'dnear'),
                              weights = c('binary', 'flin'))

# Read in the OTU tables
otu_AM <- read.csv('data/statistics/otu_AM.csv',
                   header = TRUE, row.names = 'OTU_ID')
otu_ECM <- read.csv('data/statistics/otu_ECM.csv',
                    header = TRUE, row.names = 'OTU_ID')

# Read in the taxa tables
taxa_AM <- read.csv('data/statistics/taxa_AM.csv',
                   header = TRUE, row.names = 'OTU_ID')
taxa_ECM <- read.csv('data/statistics/taxa_ECM.csv',
                    header = TRUE, row.names = 'OTU_ID')

# Generate phyloseq objects
ps_data <- sample_data(
  read.csv('data/statistics/metadata.csv', stringsAsFactors = TRUE,
           header = TRUE) %>%
    column_to_rownames(var = 'sample'))
ps_otu_AM <- otu_table(otu_AM, taxa_are_rows = TRUE)
ps_otu_ECM <- otu_table(otu_ECM, taxa_are_rows = TRUE)
ps_taxa_AM <- tax_table(
  taxa_AM %>%
    select('phylum', 'class', 'order', 'family', 'genus') %>%
    as.matrix()
  )
ps_taxa_ECM <- tax_table(
  taxa_ECM %>%
    select('phylum', 'class', 'order', 'family', 'genus') %>%
    as.matrix()
)
phyloseq_AM <- phyloseq(ps_otu_AM, ps_taxa_AM, ps_data)
phyloseq_ECM <- phyloseq(ps_otu_ECM, ps_taxa_ECM, ps_data)

### (2) Root colonisation ######################################################

#### (2a) AM colonisation ######################################################

# Fit models
binomial_AM_1 <- glmer(col_AM/weight ~ ecotype + (1 | site),
                       weights = weight,
                       family = 'binomial',
                       data = data)
binomial_AM_2 <- glmer(col_AM/weight ~ ecotype + (1 | site/sample),
                       weights = weight,
                       family = 'binomial',
                       data = data)

# Compare model performance: Model 2 is the best fit based on AIC and
# predictions based on RMSE
compare_performance(binomial_AM_1, binomial_AM_2)

# Model validation: No obvious issues
set.seed(1986)
res <- simulateResiduals(binomial_AM_2)
plot(res)
testSpatialAutocorrelation(
  res, x = data$longitude, y = data$latitude
)

# Test statistics summary
summary_col_AM <- summarise_results(
  binomial_AM_2, test_type = "z"
  ) %>%
  print(.)

# Estimated marginal means for plots
emmeans_col_AM <- emmeans(
  binomial_AM_2, ~ ecotype,
  type = 'response'
  ) %>%
  as.data.frame() %>%
  mutate(mean = prob * 100,
         upper_se = mean + (SE * 100),
         lower_se = mean - (SE * 100),
         upper_ci = asymp.UCL * 100,
         lower_ci = asymp.LCL * 100) %>%
  select(
    group = ecotype, mean, upper_se, lower_se, upper_ci, lower_ci
    ) %>%
  mutate_at(vars(-group), ~round(., 2)) %>%
  print(.)

# Jitter values for means plots
jitter_col_AM <- jitter_values(data$ecotype, binomial_AM_2) %>%
  mutate_at(vars(-group), ~ . * 100)

# Effect size coefficients for plots
effects_col_AM <- effect_size_coefficients(binomial_AM_2)

#### (2b) ECM colonisation #####################################################

# Fit models
binomial_ECM_1 <- glmer(col_ECM/weight ~ ecotype + (1 | site),
                        weights = weight,
                        family = 'binomial',
                        data = data)
binomial_ECM_2 <- glmer(col_ECM/weight ~ ecotype + (1 | site/sample),
                        weights = weight,
                        family = 'binomial',
                        data = data)

# Compare model performance: Model 2 is the best fit based on AIC predictions
# based on RMSE
compare_performance(binomial_ECM_1, binomial_ECM_2)

# Model validation: !!Refit with Moran eigen vectors!!
set.seed(1986)
res <- simulateResiduals(binomial_ECM_2)
plot(res)
testSpatialAutocorrelation(
  res, x = data$longitude, y = data$latitude
  )

# Identify significant Moran eigenvector maps
set.seed(1986)
w <- listw.select(res[["scaledResiduals"]], candidates, 
                 MEM.autocor = 'positive', method = 'MIR', nperm = 999)
w$best.id # Autocorrelation not detected

# Create data frame with signifcant MEMs
data_binomial_ECM <- cbind(data, w$best$MEM.select) %>%
  glimpse(.)

# Re-fit model accounting for spatial autocorrelation: The model failed to 
# converge so I'll try alternate optimisers
binomial_ECM_3 <- glmer(
  col_ECM/weight ~ ecotype + MEM6 + MEM4 + MEM2 + (1 | site/sample),
  weights = weight,
  family = 'binomial',
  data = data_binomial_ECM
  )

# Compare model performance: Model 3 is the best fit based on AIC predictions
# based on RMSE, which is promising
compare_performance(binomial_ECM_2, binomial_ECM_3)

# Model validation: No obvious issues
set.seed(1986)
res <- simulateResiduals(binomial_ECM_3)
plot(res)
testSpatialAutocorrelation(
  res, x = data$longitude, y = data$latitude
)

# Test statistics summary
summary_col_ECM <- summarise_results(
  binomial_ECM_3, test_type = "z"
  ) %>%
  print(.)

# Estimated marginal means for plots
emmeans_col_ECM <- emmeans(
  binomial_ECM_3, ~ ecotype,
  type = 'response'
  ) %>%
  as.data.frame() %>%
  mutate(mean = prob * 100,
         upper_se = mean + (SE * 100),
         lower_se = mean - (SE * 100),
         upper_ci = asymp.UCL * 100,
         lower_ci = asymp.LCL * 100) %>%
  select(
    group = ecotype, mean, upper_se, lower_se, upper_ci, lower_ci
  ) %>%
  mutate_at(vars(-group), ~round(., 2)) %>%
  print(.)

# Jitter values for means plots
jitter_col_ECM <- jitter_values(data$ecotype, binomial_ECM_3) %>%
  mutate_at(vars(-group), ~ . * 100)

# Effect size coefficients for plots
effects_col_ECM <- effect_size_coefficients(binomial_ECM_3)

### (3) Alpha diversity ########################################################

#### (3a) AM richness ##########################################################

# Compute OTU richness and Shannon diversity values
alpha_diversity_AM <- calculate_alpha_diversity(otu_AM) %>%
  left_join(., data, by = 'sample') %>%
  glimpse(.) %>%
  # Order levels
  mutate(ecotype = factor(ecotype, levels = c(
    "Upland", "Coastal"))
  ) %>%
  glimpse(.)

# Fit GLMMs
poisson_AM_1 <- glmer(richness ~ ecotype + (1|site),
                      family = poisson,
                      data = alpha_diversity_AM)
poisson_AM_2 <- glmer(richness ~ ecotype + (1|site/sample),
                      family = poisson,
                      data = alpha_diversity_AM)

# Compare model performance: Model 2 is the best fit based on AIC predictions
# based on RMSE
compare_performance(poisson_AM_1, poisson_AM_2)

# Model validation: No obvious issues
set.seed(1986)
res <- simulateResiduals(poisson_AM_2)
plot(res)
testSpatialAutocorrelation(
  res, x = alpha_diversity_AM$longitude, y = alpha_diversity_AM$latitude
  )

# Test statistics summary
summary_richness_AM <- summarise_results(
  poisson_AM_2, test_type = "z"
  ) %>%
  print(.)

# Estimated marginal means for plots: 
emmeans_richness_AM <- emmeans(
  poisson_AM_2, ~ ecotype,
  type = 'response'
  ) %>%
  as.data.frame() %>%
  mutate(mean = .[,2],
         upper_se = mean + SE,
         lower_se = mean - SE,
         upper_ci = asymp.UCL,
         lower_ci = asymp.LCL
         ) %>%
  select(
    group = ecotype, mean, upper_se, lower_se, upper_ci, lower_ci
  ) %>%
  mutate_at(vars(-group), ~round(., 2)) %>%
  print(.)

# Jitter values for means plots
jitter_richness_AM <- jitter_values(alpha_diversity_AM$ecotype, poisson_AM_2)

# Effect size coefficients for plots
effects_richness_AM <- effect_size_coefficients(poisson_AM_2) %>%
  print(.)

#### (3b) ECM richness #########################################################

# Compute OTU richness and Shannon diversity values
alpha_diversity_ECM <- calculate_alpha_diversity(otu_ECM) %>%
  left_join(., data, by = 'sample') %>%
  glimpse(.) %>%
  # Order levels
  mutate(ecotype = factor(ecotype, levels = c(
    "Upland", "Coastal"))
  ) %>%
  glimpse(.)

# Fit GLMMs
poisson_ECM_1 <- glmer(richness ~ ecotype + (1|site),
                      family = poisson,
                      data = alpha_diversity_ECM)
poisson_ECM_2 <- glmer(richness ~ ecotype + (1|site/sample),
                      family = poisson,
                      data = alpha_diversity_ECM)

# Compare model performance: Model 2 is the best fit based on AIC predictions
# based on RMSE
compare_performance(poisson_ECM_1, poisson_ECM_2)

# Model validation: No obvious issues
set.seed(1986)
res <- simulateResiduals(poisson_ECM_2)
plot(res)
testSpatialAutocorrelation(
  res, x = alpha_diversity_ECM$longitude, y = alpha_diversity_ECM$latitude
  )

# Test statistics summary
summary_richness_ECM <- summarise_results(
  poisson_ECM_2, test_type = "z"
  ) %>%
  print(.)

# Estimated marginal means for plots
emmeans_richness_ECM <- emmeans(
  poisson_ECM_2, ~ ecotype,
  type = 'response'
  ) %>%
  as.data.frame() %>%
  mutate(mean = .[,2],
         upper_se = mean + SE,
         lower_se = mean - SE,
         upper_ci = asymp.UCL,
         lower_ci = asymp.LCL
         ) %>%
  select(
    group = ecotype, mean, upper_se, lower_se, upper_ci, lower_ci
  ) %>%
  mutate_at(vars(-group), ~round(., 2)) %>%
  print(.)

# Jitter values for means plots
jitter_richness_ECM <- jitter_values(alpha_diversity_ECM$ecotype, poisson_ECM_2)

# Effect size coefficients for plots
effects_richness_ECM <- effect_size_coefficients(poisson_ECM_2) %>%
  print(.)

#### (3c) AM Shannon diversity ################################################

# See this GitHub issue for information on why I don't fit sample as a random
# effect: https://github.com/lme4/lme4/issues/767

# Fit GLMMs
gaussian_AM_1 <- lmer(shannon_diversity ~ ecotype + (1|site),
                      data = alpha_diversity_AM)

# Looks good
set.seed(1986)
res <- simulateResiduals(gaussian_AM_1)
plot(res)
testSpatialAutocorrelation(
  res, x = alpha_diversity_AM$longitude, y = alpha_diversity_AM$latitude
  )

# Summary of test statistics
summary_shannon_AM <- summarise_results(
  gaussian_AM_1, test_type = "t"
  ) %>%
  print(.)

#### (3d) ECM Shannon diversity ################################################

# Fit GLMMs
gaussian_ECM_1 <- lmer(shannon_diversity ~ ecotype + (1|site),
                      data = alpha_diversity_ECM)


# Borderline level of spatial autocorrelation
set.seed(1986)
res <- simulateResiduals(gaussian_ECM_1)
plot(res)
testSpatialAutocorrelation(
  res, x = alpha_diversity_ECM$longitude, y = alpha_diversity_ECM$latitude
  )

# Summary of test statistics
summary_shannon_ECM <- summarise_results(
  gaussian_ECM_1, test_type = "t"
  ) %>%
  print(.)

### (4) Beta diversity #########################################################

# Calculate relative abundances to account for differences in sequencing depth
# and reatin OTUs that occur in at least two singles

otu_relabund_AM <- prevelance_filter_relative_abundance(otu_AM, 2)
otu_relabund_ECM <- prevelance_filter_relative_abundance(otu_ECM, 2)

# AM fungi were not detected in two of the samples, and another three samples
# were lost during the prevalence filter. I need to make a metadata table to
# reflect this for downstream analyses.
data_AM <- data %>%
  filter(sample %in% rownames(otu_relabund_AM)) %>%
  print(.)

#### (4a) AM community composition #############################################

# Calculate Jaccard distances
dist_Jaccard_AM <- vegdist(otu_relabund_AM, method = 'jaccard', binary = TRUE)

# NMDS on Jaccard distances
set.seed(1986)
nmds_Jaccard_AM = metaMDS(dist_Jaccard_AM)
set.seed(1986)
nmds_Jaccard_AM = metaMDS(dist_Jaccard_AM, previous.best = nmds_Jaccard_AM)
nmds_Jaccard_AM
# Stress = 0.103

# Differences in community composition based on Jaccard distances
summary_adonis_Jaccard_AM <- adonis2(dist_Jaccard_AM ~ ecotype,
                             data = data_AM) %>%
  print(.)

# Calculate Bray Curtis distances
dist_Bray_AM <- vegdist(otu_relabund_AM, method = 'bray')

# NMDS on Bray Curtis distances
set.seed(1986)
nmds_Bray_AM = metaMDS(dist_Bray_AM, k = 3)
nmds_Bray_AM
# Stress = 0.059

# Differences in community composition based on Bray Curtis distances
summary_adonis_Bray_AM <- adonis2(dist_Bray_AM ~ ecotype,
                             data = data_AM) %>%
  print(.)

#### (4b) ECM community composition ############################################

# Calculate Jaccard distances
dist_Jaccard_ECM <- vegdist(otu_relabund_ECM, method = 'jaccard', binary = TRUE)

# NMDS on Jaccard distances
set.seed(1986)
nmds_Jaccard_ECM = metaMDS(dist_Jaccard_ECM, k = 3)
nmds_Jaccard_ECM
# Stress = 0.136

# Differences in community composition based on Jaccard distances
summary_adonis_Jaccard_ECM <- adonis2(dist_Jaccard_ECM ~ ecotype,
                             data = data) %>%
  print(.)

# Calculate Bray Curtis distances
dist_Bray_ECM <- vegdist(otu_relabund_ECM, method = 'bray')

# NMDS on Bray Curtis distances
set.seed(1986)
nmds_Bray_ECM = metaMDS(dist_Bray_ECM, k = 4)
nmds_Bray_ECM
# Stress = 0.096

# Differences in community composition based on Bray Curtis distances
summary_adonis_Bray_ECM <- adonis2(dist_Bray_ECM ~ ecotype,
                             data = data) %>%
  print(.)

### (5) Differential abundance ################################################

#### (5a) AM Fungi ####

# Run ANCOMBC
set.seed(1986)
ancom_AM <- ancombc2(
  phyloseq_AM,
  tax_level = 'genus',
  fix_formula = 'ecotype',
  rand_formula = NULL,
  p_adj_method = 'BH',
  pseudo = 0,
  pseudo_sens = TRUE,
  # This removes OTUs that are not present in at least 2 samples based on the
  # 40 samples in my dataset
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
  ) %>%
  print(.)

# Summarise statistics
summary_ancom_AM <- ancom_AM$res %>%
  select(c(taxon, 
           'lfc' = lfc_ecotypeUpland,
           'se' = se_ecotypeUpland, 
           'W-stat' = W_ecotypeUpland, 
           'p' = p_ecotypeUpland, 
           'p_adj' = q_ecotypeUpland, 
           'diff' = diff_ecotypeUpland)) %>%
  arrange('genus') %>%
  as_tibble() %>%
  print(n = Inf)

#### (5b) ECM Fungi ####

# Run ANCOMBC
set.seed(1986)
ancom_ECM = ancombc2(
  phyloseq_ECM,
  tax_level = 'genus',
  fix_formula = 'ecotype',
  rand_formula =  NULL,
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
  ) %>%
  print(.)
ancom_ECM

# Summarise statistics
summary_ancom_ECM <- ancom_ECM$res %>%
  select(c(taxon, 
           'lfc' = lfc_ecotypeUpland, 
           'se' = se_ecotypeUpland, 
           'W-stat' = W_ecotypeUpland, 
           'p' = p_ecotypeUpland, 
           'p_adj' = q_ecotypeUpland, 
           'diff' = diff_ecotypeUpland)) %>%
  arrange('genus') %>%
  as_tibble() %>%
  print(n = Inf)

### (6) Save test statistics ##################################################

# Create a list of data frames
summary_list <- list(
  colonisation_AM = summary_col_AM,
  colonisation_ECM = summary_col_ECM,
  richness_AM = summary_richness_AM,
  richness_ECM = summary_richness_ECM,
  shannon_div_AM = summary_shannon_AM,
  shannon_div_ECM = summary_shannon_ECM,
  Jaccard_dist_AM = summary_adonis_Jaccard_AM,
  Jaccard_dist_ECM = summary_adonis_Jaccard_ECM,
  BrayCurtis_dist_AM = summary_adonis_Bray_AM,
  BrayCurtis_dist_ECM = summary_adonis_Bray_ECM,
  ANCOMBC_AM = summary_ancom_AM,
  ANCOMBC_ECM = summary_ancom_ECM
)

# Write data frames to the Excel workbook
write_xlsx(summary_list, "output/test_statistics.xlsx")

### (7) Plots: Colonisation and diversity ####################################

#### (7a) Organise plot data #################################################

# Colonisation: Join means data
emmeans_col <- bind_rows(
  emmeans_col_AM %>%
    mutate(type = 'AM'),
  emmeans_col_ECM %>%
    mutate(type = 'ECM')
) %>%
  glimpse(.)

# Colonisation: Join jitter data
jitter_col <- bind_rows(
  jitter_col_AM %>%
    mutate(type = 'AM'),
  jitter_col_ECM %>%
    mutate(type = 'ECM')
) %>%
  glimpse(.)

# Colonisation: Join effect size data
# NOTE: I want Upland in the data frame for the legend but off the plot for
# clarity
effects_col <- bind_rows(
  effects_col_AM %>%
    mutate(type = 'AM'),
  effects_col_ECM %>%
    mutate(type = 'ECM')
  ) %>%
  mutate(
    Coefficient = case_when(
      ecotype == "Upland" ~ Coefficient * 100,
    TRUE ~ Coefficient
    ),
    CI_low = case_when(
      ecotype == "Upland" ~ CI_low * 100,
      TRUE ~ CI_low
    ),
    CI_high = case_when(
      ecotype == "Upland" ~ CI_high * 100,
      TRUE ~ CI_high
    ),
  ) %>%
  glimpse(.)

# Richness: Join means data
emmeans_richness <- bind_rows(
  emmeans_richness_AM %>%
    mutate(type = 'AM'),
  emmeans_richness_ECM %>%
    mutate(type = 'ECM')
  ) %>%
  glimpse(.)

# Richness: Join jitter data
jitter_richness <- bind_rows(
  jitter_richness_AM %>%
    mutate(type = 'AM'),
  jitter_richness_ECM %>%
    mutate(type = 'ECM')
) %>%
  glimpse(.)

# Richness: Join effect size data
effects_richness <- bind_rows(
  effects_richness_AM %>%
    mutate(type = 'AM'),
  effects_richness_ECM %>%
    mutate(type = 'ECM')
  ) %>%
  mutate(
    Coefficient = case_when(
      ecotype == "Upland" ~ Coefficient * 100,
      TRUE ~ Coefficient
      ),
    CI_low = case_when(
      ecotype == "Upland" ~ CI_low * 100,
      TRUE ~ CI_low
      ),
    CI_high = case_when(
      ecotype == "Upland" ~ CI_high * 100,
      TRUE ~ CI_high
      )
    ) %>%
  glimpse(.)

# Community composition: Extract nmds points and join to metadata
# AM
points_AM <- inner_join(
  data_AM,
  nmds_Bray_AM$points %>%
    as.data.frame(.) %>%
    rownames_to_column('sample'),
  by = 'sample'
) %>%
  glimpse(.)
# ECM
points_ECM <- inner_join(
  data,
  nmds_Bray_ECM$points %>%
    as.data.frame(.) %>%
    rownames_to_column('sample'),
  by = 'sample'
) %>%
  glimpse(.)

# Community composition: Calculate centroids
# AM
centroid_AM <- points_AM %>%
  group_by(ecotype) %>%
  summarise(
    NMDS1 = mean(MDS1),
    NMDS2 = mean(MDS2),
    .groups = 'drop'
  )
# ECM
centroid_ECM <- points_ECM %>%
  group_by(ecotype) %>%
  summarise(
    NMDS1 = mean(MDS1),
    NMDS2 = mean(MDS2),
    .groups = 'drop'
  )

#### (7b) Colonisation means #################################################

# Plot colonisation means
means_plot_col <- ggplot(
  emmeans_col,
  aes(x = type,
      ymin = lower_ci,
      lower = lower_se,
      middle = mean,
      upper = upper_se,
      ymax = upper_ci,
      fill = group)
  ) +
  geom_boxplot(
    aes(),
    stat = 'identity',
    position = position_dodge(width = 0.5),
    width = 0.45
  ) +
  geom_point(
    data = jitter_col,
    aes(x = type, y = value, fill = group), 
    position = position_jitterdodge(
      seed = 16, dodge.width = 0.5, jitter.width = 0.25
      ),
    alpha = 0.6,
    shape = 21,
    stroke = 0.2,
    size = 1.5
  ) +
  scale_fill_manual(
    values = c('#1b9e77', '#d95f02'),
    breaks = c('Coastal', 'Upland')
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.25),
    legend.position = 'none',
    axis.text.x = element_text(colour = 'black', size = rel(1)),
    axis.text.y = element_text(colour = 'black', size = rel(1)),
    axis.title = element_text(colour = 'black', size = rel(0.9)),
    axis.ticks = element_blank()
  ) +
  scale_y_continuous(
    limits = c(-0.5, max(jitter_col$mean) + 0.01),
    breaks = c(0, 20, 40, 60)
  ) +
  xlab(element_blank()) +
  ylab('Root colonisation (%)') +
  annotate('text', x = 1, y = 0,
           label = 'paste(\ italic(z),\' = 1.19 \')',
           parse = T, size = 3, hjust = 0.5, vjust = 1)  +
  annotate('text', x = 2, y = 0,
           label = 'paste(\ italic(z),\' = -9.16*** \')',
           parse = T, size = 3, hjust = 0.5, vjust = 1)
means_plot_col

#### (7c) Colonisation effect sizes ##########################################

effects_plot_col <- ggplot(
  effects_col,
  aes(
    x = type,
    y = Coefficient,
    ymin = CI_low,
    ymax = CI_high,
    group = ecotype,
    colour = ecotype)
  ) +
  geom_hline(
    yintercept = 0,
    color = '#d95f02',
    linetype = "dotted",
    linewidth = 1
  ) +
  geom_errorbar(
    color = "black",
    width = 0,
    linewidth = 0.7
  ) +
  geom_point(
    shape = 16,
    size = 4
  ) +
  scale_colour_manual(
    values = c('#d95f02', '#1b9e77'),
    breaks = c('Upland', 'Coastal')
  ) +
  scale_y_continuous(
    limits = c(-3, 2)
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.25),
    axis.text.x = element_text(colour = 'black', size = rel(1)),
    axis.text.y = element_text(colour = 'black', size = rel(1)),
    axis.title = element_text(colour = 'black', size = rel(0.9)),
    axis.ticks = element_blank(),
    # Modify legend
    legend.position = c(1, 1), 
    legend.justification = c(1.05, 1.1),
    legend.title = element_blank(),
    legend.text = element_text(size = rel(0.9)),
    legend.key.size = unit(0.25, "cm"),
    legend.background = element_rect(fill="NA",
                                     color="black",
                                     linewidth = 0.3,
                                     linetype = 'dotted'),
    legend.margin = margin(t=-2.5, r=5, b=2, l=3),
    legend.key = element_blank()
  ) +
  xlab(element_blank()) +
  ylab('Mean effect size (log-odds)')
effects_plot_col

# Join plots
colonisation_plots <- cowplot::plot_grid(
  means_plot_col,
  effects_plot_col,
  labels = c('', ''),
  ncol = 2
  )

#### (7d) Richness means #################################################

# Plot richness means
means_plot_richness <- ggplot(
  emmeans_richness,
  aes(x = type,
      ymin = lower_ci,
      lower = lower_se,
      middle = mean,
      upper = upper_se,
      ymax = upper_ci,
      fill = group)
) +
  geom_boxplot(
    aes(),
    stat = 'identity',
    position = position_dodge(width = 0.5),
    width = 0.45
  ) +
  geom_point(
    data = jitter_richness,
    aes(x = type, y = value, fill = group), 
    position = position_jitterdodge(
      seed = 16, dodge.width = 0.5, jitter.width = 0.25
    ),
    alpha = 0.6,
    shape = 21,
    stroke = 0.2,
    size = 1.5
  ) +
  scale_fill_manual(
    values = c('#1b9e77', '#d95f02'),
    breaks = c('Coastal', 'Upland')
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.25),
    legend.position = 'none',
    axis.text.x = element_text(colour = 'black', size = rel(1)),
    axis.text.y = element_text(colour = 'black', size = rel(1)),
    axis.title = element_text(colour = 'black', size = rel(0.9)),
    axis.ticks = element_blank()
  ) +
  scale_y_continuous(
    limits = c(0, max(jitter_richness$mean) + 0.01),
    breaks = c(0, 10, 20)
  ) +
  xlab(element_blank()) +
  ylab('OTU richness') +
  annotate('text', x = 1, y = 0.25,
           label = 'paste(\ italic(z),\' = 2.52** \')',
           parse = T, size = 3, hjust = 0.5, vjust = 1)  +
  annotate('text', x = 2, y = 0.25,
           label = 'paste(\ italic(z),\' = 3.21*** \')',
           parse = T, size = 3, hjust = 0.5, vjust = 1)
means_plot_richness

#### (7e) Richness effect sizes ##########################################

effects_plot_richness <- ggplot(
  effects_richness,
  aes(
    x = type,
    y = Coefficient,
    ymin = CI_low,
    ymax = CI_high,
    group = ecotype,
    colour = ecotype)
  ) +
  geom_hline(
    yintercept = 0,
    color = '#d95f02',
    linetype = "dotted",
    size = 1
  ) +
  geom_errorbar(
    color = "black",
    width = 0,
    linewidth = 0.7
  ) +
  geom_point(
    shape = 16,
    size = 4
  ) +
  scale_colour_manual(
    values = c('#d95f02', '#1b9e77'),
    breaks = c('Upland', 'Coastal')
  ) +
  scale_y_continuous(
    limits = c(-0.25, 1.75), 
    breaks = c(0, 1)
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.25),
    axis.text.x = element_text(colour = 'black', size = rel(1)),
    axis.text.y = element_text(colour = 'black', size = rel(1)),
    axis.title = element_text(colour = 'black', size = rel(0.9)),
    axis.ticks = element_blank(),
    # Modify legend
    legend.position = "none"
  ) +
  xlab(element_blank()) +
  ylab('Mean effect size (log-rates)')
effects_plot_richness

# Join plots
richness_plots <- cowplot::plot_grid(
  means_plot_richness,
  effects_plot_richness,
  labels = c('', ''),
  ncol = 2
)

#### (7f) Community composition ###############################################

# AM fungi: Bray-Curtis dissimilarity on NMDS #

# Plot the ordination
ordi_plot_AM <- ggplot(
  points_AM,
  aes(
    MDS1, MDS2, colour = ecotype, fill = ecotype)
  ) +
  stat_ellipse(
    geom = 'polygon',
    level = 0.95,
    alpha = 0.1, 
    show.legend = F,
    linewidth = 0.1
  ) +
  geom_point() +
  scale_colour_manual(breaks = c('Upland', 'Coastal'),
                      values = c('#d95f02', '#1b9e77', 'black'))  +
  scale_fill_manual(breaks = c('Upland', 'Coastal'),
                    values = c('#d95f02', '#1b9e77', 'black'))  +
  geom_point(
    centroid_AM, mapping = aes(NMDS1, NMDS2),
    fill = c('#d95f02', '#1b9e77'),
    colour = 'black', 
    shape = 21,
    size = 4,
    stroke = 0.5,
    show.legend = F
  ) +
  labs(
    x = 'NMDS1',
    y = 'NMDS2'
    ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.25),
    axis.text.x = element_text(colour = 'black', size = rel(1)),
    axis.text.y = element_text(colour = 'black', size = rel(1)),
    axis.ticks = element_blank(),
    axis.title = element_text(colour = 'black', size = rel(0.85)),
    # Modify legend
    legend.position = "none"
  ) +
  scale_y_continuous(
    limits = c(-2.5, 2.3),
    breaks = c(-2, 0, 2)
    ) +
  scale_x_continuous(
    limits = c(-2.6, 2.2),
    breaks = c(-2, 0, 2)
    ) +
  annotate('text', x = -2.6, y = 2.25,
           label = "AM",
           parse = T, size = 3.4, hjust = 0, vjust = 1)  +
  annotate('text', x = -2.6, y = -2.375,
           label = "paste(\ 'Stress = 0.059')",
           parse = T, size = 3, hjust = 0, vjust = 1)  +
  annotate('text', x = 2, y = -2.375,
           label = "paste(\ '' , \ italic(F) ['1,36'], \ ' = 3.70** \')",
           parse = T, size = 3, hjust = 0.9, vjust = 1)
ordi_plot_AM

# ECM fungi: Bray-Curtis dissimilarity on NMDS #

# Plot the ordination
ordi_plot_ECM <- ggplot(
  points_ECM,
  aes(
    MDS1, MDS2, colour = ecotype, fill = ecotype)
) +
  stat_ellipse(
    geom = 'polygon',
    level = 0.95,
    alpha = 0.1, 
    show.legend = F,
    linewidth = 0.1
  ) +
  geom_point() +
  scale_colour_manual(breaks = c('Upland', 'Coastal'),
                      values = c('#d95f02', '#1b9e77', 'black'))  +
  scale_fill_manual(breaks = c('Upland', 'Coastal'),
                    values = c('#d95f02', '#1b9e77', 'black'))  +
  geom_point(
    centroid_ECM, mapping = aes(NMDS1, NMDS2),
    fill = c('#d95f02', '#1b9e77'),
    colour = 'black', 
    shape = 21,
    size = 4,
    stroke = 0.5,
    show.legend = F
  ) +
  labs(
    x = 'NMDS1',
    y = 'NMDS2'
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.25),
    axis.text.x = element_text(colour = 'black', size = rel(1)),
    axis.text.y = element_text(colour = 'black', size = rel(1)),
    axis.title = element_text(colour = 'black', size = rel(0.85)),
    axis.ticks = element_blank(),
    # Modify legend
    legend.position = "none"
  ) +
  scale_y_continuous(
    limits = c(-1.8, 1.6),
    breaks = c(-1, 0, 1)
    ) +
  scale_x_continuous(
    limits = c(-1.6, 1.8),
    breaks = c(-1, 0, 1)
    ) +
  annotate('text', x = -1.6, y = 1.575,
            label = "ECM",
            parse = T, size = 3.4, hjust = 0, vjust = 1) +
  annotate('text', x = -1.55, y = -1.7,
           label = "paste(\ 'Stress = 0.096')",
           parse = T, size = 3, hjust = 0, vjust = 1) +
  annotate('text', x = 1.6, y = -1.7,
           label = "paste(\ '' , \ italic(F) ['1,38'], \ ' = 2.29** \')",
           parse = T, size = 3, hjust = 0.9, vjust = 1)
ordi_plot_ECM

# Join plots
ordination_plots <- cowplot::plot_grid(
  ordi_plot_AM,
  ordi_plot_ECM,
  labels = c('', ''),
  ncol = 2
)

#### (7g) Join plots ##########################################################

figure_5 <- cowplot::plot_grid(
  means_plot_col,
  effects_plot_col,
  means_plot_richness,
  effects_plot_richness,
  ordi_plot_AM,
  ordi_plot_ECM,
  labels = c('a', '', 'b', '', 'c', ''),
  ncol = 2
)

# Save plots
ggsave('output/figure_5.jpg', width = 5, height = 7.5)
ggsave('output/figure_5.tiff', width = 5, height = 7.5)
ggsave('output/figure_5.pdf', width = 5, height = 7.5)

### (8) Plot differential abundances ###############################################################

##### (2a) Organise data #####

# Plots will include all genera with relative abundance > 2 %

# Calculate AM genera with relative abundances > 2
rel_abund_AM <- otu_AM %>%
  rownames_to_column(var = "OTU_ID") %>%
  pivot_longer(-OTU_ID) %>%
  select(-name) %>%
  inner_join(.,
  taxa_AM %>%
    rownames_to_column(var = "OTU_ID") %>%
    select(OTU_ID, taxon = genus)
  ) %>%
  filter(value > 0) %>%
  mutate(
    total_abundance = sum(value)
  ) %>%
  group_by(taxon, total_abundance) %>%
  summarise(
    taxon_abundance = sum(value)
    ) %>%
  summarise(
    rel_abund = taxon_abundance / total_abundance * 100
  ) %>%
  filter(
    rel_abund >= 2,
    !grepl("gen_Incertae_sedis", taxon)
  )
  
# Calculate ECM genera with relative abundances > 2
rel_abund_ECM <- otu_ECM %>%
  rownames_to_column(var = "OTU_ID") %>%
  pivot_longer(-OTU_ID) %>%
  select(-name) %>%
  inner_join(
    .,
    taxa_ECM %>%
      rownames_to_column(var = "OTU_ID") %>%
      select(OTU_ID, taxon = genus)
  ) %>%
  filter(value > 0) %>%
  mutate(
    total_abundance = sum(value)
  ) %>%
  group_by(taxon, total_abundance) %>%
  summarise(
    taxon_abundance = sum(value)
  ) %>%
  summarise(
    rel_abund = taxon_abundance / total_abundance * 100
  ) %>%
  filter(
    rel_abund >= 2
  )

##### (2b) Plot AM Fungi #####

ancom_plot_AM <- rel_abund_AM %>%
  inner_join(summary_ancom_AM, by = 'taxon') %>%
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
    yintercept = 0, linetype = "dotted",
    colour = 'darkgrey', linewidth = 0.5
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
  geom_text(aes(y = lfc + se + 1.4), size = 7) + # Adjust astrix position
  geom_point(
    size = 4,
    shape = 21) +
  scale_y_continuous(
    limits = c(-7, 9),
    breaks = c(-6, -4, -2, 0, 2, 4, 6, 8)
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
    axis.text.x = element_text(
      angle = 45, hjust = 1, size = rel(1.2), colour = 'black'
      ),
    text = element_text(size = 10),
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    legend.position = 'none'
  )
ancom_plot_AM

##### (2c) Plot EM Fungi #####

ancom_plot_ECM <- rel_abund_ECM %>%
  inner_join(summary_ancom_ECM, by = 'taxon') %>%
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
    yintercept = 0, linetype = "dotted",
    colour = 'darkgrey', linewidth = 0.5
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
  geom_text(aes(y = lfc + se + 1.5), size = 7) + # Adjust astrix position
  geom_point(
    size = 4,
    shape = 21) +
  scale_y_continuous(
    limits = c(-7, 9)
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
ancom_plot_ECM

##### (2d) Join plots #####
require(patchwork)
ancom_plot = ancom_plot_AM + ancom_plot_ECM
ancom_plot

ggsave('output/figure_4.pdf', width = 5, height = 3.75)
ggsave('output/figure_4.jpg', width = 5, height = 3.75)
ggsave('output/figure_4.tiff', width = 5, height = 3.75)

