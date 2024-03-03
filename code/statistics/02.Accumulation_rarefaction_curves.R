
# Script: Rarefaction and species acculmilation curves

# Required packages and functions
source("code/statistics/functions.R")
require(patchwork)
require(scales)
require(vegan)
require(tidyverse)


# Metadata
data <- read.csv('data/statistics/metadata.csv',
                 header = T)
# Taxa
taxa <- read.csv('data/statistics/taxa_fungi.csv',
                 header = TRUE, row.names = 'OTU_ID')

# Fungal OTUs
otu <- read.csv('data/statistics/otu_fungi.csv', header = TRUE) %>%
  pivot_longer(-OTU_ID, names_to = 'sample')
# AM fungal OTUs
otu_AM = read.csv('data/statistics/otu_AM.csv', header = TRUE)  %>%
  pivot_longer(-OTU_ID, names_to = 'sample')
# EM fungal OTUs
otu_ECM = read.csv('data/statistics/otu_ECM.csv', header = TRUE) %>%
  pivot_longer(-OTU_ID, names_to = 'sample')
                 
#### (1) Calculate sequencing depth and OTU richness per sample ################

# Sequencing depth per sample
sample_depth <- otu %>%
  group_by(sample) %>%
  summarise(n_seqs = sum(value)) %>%
  arrange(n_seqs) %>%
  print(n = Inf)

# Sample depth range
sample_depth$n_seqs %>% range
# Minimum sampling depth = 51259
min_depth <- sample_depth$n_seqs[1]
# Minimum sampling depth = 380142
max_depth <- sample_depth$n_seqs[40]
# Median sampling depth
mid_depth <- median(sample_depth$n_seqs)

# Visualise sample depth and range
sample_depth %>%
  ggplot(aes(x = 1:nrow(.), y = n_seqs)) +
  geom_line() +
  geom_point()
sample_depth %>%
  ggplot(aes(x = 1:nrow(.), y = n_seqs, label = sample)) +
  geom_line() +
  geom_point() +
  geom_label()
sample_depth %>%
  ggplot(aes(x = 1, y = n_seqs)) +
  geom_jitter()

# OTU richness per sample
OTU_richness <- otu %>%
  filter(value != 0) %>%
  group_by(sample) %>%
  summarise(n_OTUs = n_distinct(OTU_ID)) %>%
  arrange(n_OTUs) %>%
  print(n = Inf)
# OTU richness range
OTU_richness$n_OTUs %>% range
# Minimum richness = 14
min_rich <- OTU_richness$n_OTUs[1]
# Max richness = 113
max_rich <- OTU_richness$n_OTUs[40]

#### (2) Species accumulation curves ###########################################

##### (2a) All fungi #####

# Subset coastal OTUs for all fungi
otu_coastal <- otu %>%
  # Subset data to coastal site
  inner_join(data, by = 'sample') %>%
  filter(ecotype == 'Coastal') %>%
  select(sample, OTU_ID, value) %>%
  # Remove OTUs that do not occur in this subset
  group_by(OTU_ID) %>%
  mutate(total = sum(value)) %>%
  filter(total != 0) %>%
  ungroup() %>%
  select(-total)%>%
  pivot_wider(names_from = "OTU_ID") %>%
  column_to_rownames(var = 'sample')

# Subset upland OTUs for all fungi
otu_upland <- otu %>%
  # Subset data to upland sites
  inner_join(data, by = 'sample') %>%
  filter(ecotype == 'Upland') %>%
  select(sample, OTU_ID, value) %>%
  # Remove OTUs that do not occur in this subset
  group_by(OTU_ID) %>%
  mutate(total = sum(value)) %>%
  filter(total != 0) %>%
  ungroup() %>%
  select(-total) %>%
  pivot_wider(names_from = 'OTU_ID') %>%
  column_to_rownames(var = 'sample')

# Coastal species accumulation curve
species_acc_coastal <- specaccum(otu_coastal, 'random')
plot(species_acc_coastal)
# Upland species accumulation curve
species_acc_upland <- specaccum(otu_upland, 'random')
plot(species_acc_upland)

# Check out this other option for estimating species accumulation
pool_coast <- poolaccum(otu_coastal)
plot(pool_coast)
pool_upland = poolaccum(otu_upland)
plot(pool_upland)

##### (2b) AM fungi #####

# Coastal OTUs for AM fungi
otu_AM_coastal <- otu_AM %>%
  # Append data required for sub-setting (here ecotype)
  inner_join(data, by = 'sample') %>%
  filter(ecotype == 'Coastal') %>%
  select(sample, OTU_ID, value) %>%
  # Remove OTUs that do not occur in this subset
  group_by(OTU_ID) %>%
  mutate(total = sum(value)) %>%
  filter(total != 0) %>%
  ungroup() %>%
  select(-total) %>%
  pivot_wider(names_from = 'OTU_ID') %>%
  column_to_rownames(var = 'sample')

# Upland OTUs for AM fungi
otu_AM_upland = otu_AM %>%
  # Append data required for sub-setting (here ecotype)
  inner_join(data, by = 'sample') %>%
  filter(ecotype == 'Upland') %>%
  select(sample, OTU_ID, value) %>%
  # Remove OTUs that do not occur in this subset
  group_by(OTU_ID) %>%
  mutate(total = sum(value)) %>%
  filter(total != 0) %>%
  ungroup() %>%
  select(-total) %>%
  pivot_wider(names_from = 'OTU_ID') %>%
  column_to_rownames(var = 'sample')

pool_AM_coast <- poolaccum(otu_AM_coastal)
plot(pool_AM_coast)
pool_AM_upland <- poolaccum(otu_AM_upland)
plot(pool_AM_upland)

##### (2c) EM fungi #####

# Coastal OTUs for EM fungi
otu_ECM_coastal <- otu_ECM %>%
  # Append data required for sub-setting (here ecotype)
  inner_join(data, by = 'sample') %>%
  filter(ecotype == 'Coastal') %>%
  select(sample, OTU_ID, value) %>%
  # Remove OTUs that do not occur in this subset
  group_by(OTU_ID) %>%
  mutate(total = sum(value)) %>%
  filter(total != 0) %>%
  ungroup() %>%
  select(-total) %>%
  pivot_wider(names_from = 'OTU_ID') %>%
  column_to_rownames(var = 'sample')

# Upland OTUs for AM fungi
otu_ECM_upland <- otu_ECM %>%
  # Append data required for sub-setting (here ecotype)
  inner_join(data, by = 'sample') %>%
  filter(ecotype == 'Upland') %>%
  select(sample, OTU_ID, value) %>%
  # Remove OTUs that do not occur in this subset
  group_by(OTU_ID) %>%
  mutate(total = sum(value)) %>%
  filter(total != 0) %>%
  ungroup() %>%
  select(-total) %>%
  pivot_wider(names_from= 'OTU_ID') %>%
  column_to_rownames(var = 'sample')

pool_ECM_coast <- poolaccum(otu_ECM_coastal)
plot(pool_ECM_coast)
pool_ECM_upland <- poolaccum(otu_ECM_upland)
plot(pool_ECM_upland)

#### (3) Plots accumulation curve #############################################

##### (3a) All fungi #####

# Max y-axis value for all fungi plots
ymax <- summary(pool_coast)$S %>%
  as_tibble() %>%
  select("97.5%") %>%
  max()

# Coastal plot for all fungi
plot_coast <- data.frame(summary(pool_coast)$S, check.names = FALSE) %>%
  select(N = N, S = S, lowCI = '2.5%', upCI = '97.5%', SD = Std.Dev) %>%
  ggplot(data = ., aes(x = N,
                       y = S)) +
  # Add confidence intervals
  geom_ribbon(aes(ymin = lowCI,
                  ymax = upCI),
              alpha = 0.5,
              colour = "gray70") +
  # Add observed richness line 
  geom_line() +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    axis.text = element_text(colour = 'black'),
    axis.title = element_text(size = rel(1)),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) + 
  scale_y_continuous(limits = c(0, ymax)) +
  labs(title = 'Coastal sites',
       x = NULL,
       y = 'Number of all fungal OTUs')

# Upland plot for all fungi
plot_upland <- data.frame(summary(pool_upland)$S, check.names = FALSE) %>%
  select(N = N, S = S, lowCI = '2.5%', upCI = '97.5%', SD = Std.Dev) %>%
  ggplot(data = ., aes(x = N,
                       y = S)) +
  # Add confidence intervals
  geom_ribbon(aes(ymin = lowCI,
                  ymax = upCI),
              alpha = 0.5,
              colour = "gray70") +
  # Add observed richness line 
  geom_line() +
  scale_y_continuous(limits = c(0, ymax),
                     labels = NULL) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    axis.text = element_text(colour = 'black'),
    axis.title = element_text(size = rel(1)),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(title = 'Upland sites',
       x = NULL,
       y = NULL)

##### (3b) AM fungi #####

# Max y-axis value for AM fungi plots
ymax <- summary(pool_AM_coast)$S %>%
  as_tibble() %>%
  select("97.5%") %>%
  max()

# Coastal plot for AM fungi
plot_AM_coast <- data.frame(summary(pool_AM_coast)$S, check.names = FALSE) %>%
  select(N = N, S = S, lowCI = '2.5%', upCI = '97.5%', SD = Std.Dev) %>%
  ggplot(data = ., aes(x = N,
                       y = S)) +
  # Add confidence intervals
  geom_ribbon(aes(ymin = lowCI,
                  ymax = upCI),
              alpha = 0.5,
              colour = "gray70") +
  # Add observed richness line 
  geom_line() +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    axis.text = element_text(colour = 'black'),
    axis.title = element_text(size = rel(1)),
    axis.ticks = element_blank()
  ) + 
  scale_y_continuous(limits = c(0, ymax)) +
  labs(x = NULL,
       y = 'Number of AM fungal OTUs')

# Upland plot for AM fungi
plot_AM_upland <- data.frame(summary(pool_AM_upland)$S, check.names = FALSE) %>%
  select(N = N, S = S, lowCI = '2.5%', upCI = '97.5%', SD = Std.Dev) %>%
  ggplot(data = ., aes(x = N,
                       y = S)) +
  # Add confidence intervals
  geom_ribbon(aes(ymin = lowCI,
                  ymax = upCI),
              alpha = 0.5,
              colour = "gray70") +
  # Add observed richness line 
  geom_line() +
  scale_y_continuous(limits = c(0, 80),
                     labels = NULL) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    axis.text = element_text(colour = 'black'),
    axis.title = element_text(size = rel(1)),
    axis.ticks = element_blank()
  ) +
  labs(x = NULL,
       y = NULL)

##### (3c) EM fungi #####

# Max y-axis value for EM fungi plots
ymax <- summary(pool_ECM_coast)$S %>%
  as_tibble() %>%
  select("97.5%") %>%
  max()

# Coastal plot for EM fungi
plot_ECM_coast <- data.frame(summary(pool_ECM_coast)$S, check.names = FALSE) %>%
  select(N = N, S = S, lowCI = '2.5%', upCI = '97.5%', SD = Std.Dev) %>%
  ggplot(data = ., aes(x = N,
                       y = S)) +
  # Add confidence intervals
  geom_ribbon(aes(ymin = lowCI,
                  ymax = upCI),
              alpha = 0.5,
              colour = "gray70") +
  # Add observed richness line 
  geom_line() +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    axis.text = element_text(colour = 'black'),
    axis.title = element_text(size = rel(1)),
    axis.ticks = element_blank()
  ) + 
  scale_y_continuous(limits = c(0, ymax)) +
  labs(x = 'Number of root samples',
       y = 'Number of EcM fungal OTUs')

# Upland plot for EM fungi
plot_ECM_upland <- data.frame(summary(pool_ECM_upland)$S, check.names = FALSE) %>%
  select(N = N, S = S, lowCI = '2.5%', upCI = '97.5%', SD = Std.Dev) %>%
  ggplot(data = ., aes(x = N,
                       y = S)) +
  # Add confidence intervals
  geom_ribbon(aes(ymin = lowCI,
                  ymax = upCI),
              alpha = 0.5,
              colour = "gray70") +
  # Add observed richness line 
  geom_line() +
  scale_y_continuous(limits = c(0, ymax),
                     labels = NULL) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    axis.text = element_text(colour = 'black'),
    axis.title = element_text(size = rel(1)),
    axis.ticks = element_blank()
  ) + 
  labs(x = 'Number of root samples',
       y = NULL)

##### (3d) Combine plots #####

plot_coast + plot_upland + plot_AM_coast + plot_AM_upland + 
  plot_ECM_coast + plot_ECM_upland + plot_layout(ncol = 2)

ggsave('output/curve_accumilation.pdf', width = 5, height = 7.5)
ggsave('output/curve_accumilation.jpg', width = 5, height = 7.5)
ggsave('output/curve_accumilation.tiff', width = 5, height = 7.5)

#### (4) Rarefaction curves ####################################################

otu_wide <- otu %>%
  pivot_wider(names_from = 'OTU_ID', values_from = 'value') %>%
  column_to_rownames(var = 'sample')

rare_curve <- rarecurve(otu_wide, step = 100)
map_dfr(rare_curve, bind_rows) %>%
  bind_cols(sample = rownames(otu_wide),.) %>%
  pivot_longer(-sample) %>%
  drop_na() %>%
  mutate(depth = as.numeric(str_replace(name, 'N', ''))) %>%
  select(-name) %>%
  ggplot(aes(x = depth, y = value, sample = sample)) +
  geom_vline(xintercept = min_depth, colour = 'grey') +
  # geom_vline(xintercept = low_depth, colour = 'grey', linetype = 'dotted') +
  geom_line() +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black'),
    axis.ticks = element_blank()
  ) + 
  scale_x_continuous(labels = comma_format(big.mark = ' '),
                     limits = c(0, max_depth),
                     breaks = c(0, 100000, 200000, 300000)) +
  scale_y_continuous(labels = comma_format(big.mark = ' ')) +
  xlab('Number of reads') +
  ylab('Number of OTUs')

ggsave('output/rarecurve.pdf', width = 4, height = 4)
ggsave('output/rarecurve.jpg', width = 4, height = 4)
ggsave('output/rarecurve.tiff', width = 4, height = 4)

# Clean the environment
rm(list = ls())
