
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
otuAM = read.csv('data/statistics/otuAM.csv', header = TRUE)  %>%
  pivot_longer(-OTU_ID, names_to = 'sample')
# EM fungal OTUs
otuEM = read.csv('data/statistics/otuEM.csv', header = TRUE) %>%
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
# Minimum sampling depth = 44 569
min_depth <- sample_depth$n_seqs[1]
# Minimum sampling depth = 375 649
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
# Minimum richness = 6
min_rich <- OTU_richness$n_OTUs[1]
# Max richness = 72
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
otuAM_coastal <- otuAM %>%
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
otuAM_upland = otuAM %>%
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

poolAM_coast <- poolaccum(otuAM_coastal)
plot(poolAM_coast)
poolAM_upland <- poolaccum(otuAM_upland)
plot(poolAM_upland)

##### (2c) EM fungi #####

# Coastal OTUs for EM fungi
otuEM_coastal <- otuEM %>%
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
otuEM_upland <- otuEM %>%
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

poolEM_coast <- poolaccum(otuEM_coastal)
plot(poolEM_coast)
poolEM_upland <- poolaccum(otuEM_upland)
plot(poolEM_upland)

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
  # geom_text(label = 'a', x = 3, y = 1600, colour = 'black') +
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
  # geom_text(label = 'b', x = 3, y = 1600, colour = 'black') +
  labs(title = 'Upland sites',
       x = NULL,
       y = NULL)

##### (3b) AM fungi #####

# Max y-axis value for AM fungi plots
ymax <- summary(poolAM_coast)$S %>%
  as_tibble() %>%
  select("97.5%") %>%
  max()

# Coastal plot for AM fungi
plotAM_coast <- data.frame(summary(poolAM_coast)$S, check.names = FALSE) %>%
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
  # geom_text(label = 'a', x = 3, y = 100, colour = 'black') +
  labs(x = NULL,
       y = 'Number of AM fungal OTUs')

# Upland plot for AM fungi
plotAM_upland <- data.frame(summary(poolAM_upland)$S, check.names = FALSE) %>%
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
  # geom_text(label = 'b', x = 3, y = 50, colour = 'black') +
  labs(x = NULL,
       y = NULL)

##### (3c) EM fungi #####

# Max y-axis value for EM fungi plots
ymax <- summary(poolEM_coast)$S %>%
  as_tibble() %>%
  select("97.5%") %>%
  max()

# Coastal plot for EM fungi
plotEM_coast <- data.frame(summary(poolEM_coast)$S, check.names = FALSE) %>%
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
  # geom_text(label = 'a', x = 3, y = 100, colour = 'black') +
  labs(x = 'Number of root samples',
       y = 'Number of EcM fungal OTUs')

# Upland plot for EM fungi
plotEM_upland <- data.frame(summary(poolEM_upland)$S, check.names = FALSE) %>%
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
  # geom_text(label = 'b', x = 3, y = 50, colour = 'black') +
  labs(x = 'Number of root samples',
       y = NULL)

##### (3d) Combine plots #####

plot_coast + plot_upland + plotAM_coast + plotAM_upland + 
  plotEM_coast + plotEM_upland + plot_layout(ncol = 2)

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
