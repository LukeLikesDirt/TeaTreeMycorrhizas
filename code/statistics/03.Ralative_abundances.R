
# Script: Calculate and plot relative abundances of mycorrhizal genera across
# tea tree ecotypes

# Required packages and functions
source("code/statistics/functions.R")
require(patchwork)
require(paletteer)
require(tidyverse)

# Metadata
data <- read.csv('data/statistics/metadata.csv', 
                 header = TRUE, stringsAsFactors = TRUE) %>%
  select(sample, ecotype) %>%
  glimpse()

# OTUs
otuAM = t(read.csv('data/statistics/otuAM.csv',
                   h = TRUE, row.names = 'OTU_ID')) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample') %>%
  as_tibble() %>%
  pivot_longer(-sample, names_to = 'OTU_ID', values_to = 'count')
# AM taxa
taxaAM = read.csv('data/statistics/taxaAM.csv', header = TRUE)

# EM OTUs
otuEM = t(read.csv('data/statistics/otuEM.csv',
                   h = TRUE, row.names = 'OTU_ID')) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample') %>%
  as_tibble() %>%
  pivot_longer(-sample, names_to = 'OTU_ID', values_to = 'count')
# EM taxa
taxaEM = read.csv('data/statistics/taxaEM.csv', header = TRUE)

#### AM: Relative abundance and richness ################################

# AM: Relative abundance of genera
relabund_AM <- inner_join(otuAM, taxaAM, by = 'OTU_ID') %>%
  select(genus, count) %>%
  mutate(rel_abund = relabund(count) * 100) %>%
  group_by(genus) %>%
  summarise(rel_abund = round(sum(rel_abund), digits = 1))
sum(relabund_AM$rel_abund)
# AM: Richness of genera
richness_AM <- inner_join(otuAM, taxaAM, by = 'OTU_ID') %>%
  filter(count != 0) %>%
  distinct(OTU_ID, .keep_all = TRUE) %>%
  select(OTU_ID, genus) %>%
  count(genus)
sum(richness_AM$n)

# AM: Relative abundance coastal
relabund_coastal_AM <-
  inner_join(otuAM, taxaAM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  select(ecotype, genus, count) %>%
  filter(ecotype == 'Coastal') %>%
  mutate(rel_abund = relabund(count) * 100) %>%
  group_by(genus) %>%
  summarise(rel_abund_coastal = round(sum(rel_abund), digits = 1)) %>%
  print(n = Inf)
sum(relabund_coastal_AM$rel_abund_coastal)
# AM: Richness coastal
richness_coastal_AM <- inner_join(otuAM, taxaAM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  filter(ecotype == 'Coastal') %>%
  filter(count != 0) %>%
  distinct(OTU_ID, .keep_all = TRUE) %>%
  select(OTU_ID, genus) %>%
  count(genus) %>%
  rename(n_coastal = n) %>%
  print(n = Inf)

# AM: Relative abundance upland
relabund_upland_AM <-
  inner_join(otuAM, taxaAM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  select(ecotype, genus, count) %>%
  filter(ecotype == 'Upland') %>%
  mutate(rel_abund = relabund(count) * 100) %>%
  group_by(genus) %>%
  summarise(rel_abund_upland = round(sum(rel_abund), digits = 1)) %>%
  print(n = Inf)
sum(relabund_upland_AM$rel_abund_upland)
# AM: Richness upland
richness_upland_AM <- inner_join(otuAM, taxaAM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  filter(ecotype == 'Upland') %>%
  filter(count != 0) %>%
  distinct(OTU_ID, .keep_all = TRUE) %>%
  select(OTU_ID, genus) %>%
  count(genus) %>%
  rename(n_upland = n) %>%
  print(n = Inf)

# Export relative abundance and richness table
left_join(relabund_AM, richness_AM, by = 'genus') %>%
  left_join(relabund_coastal_AM, by = 'genus') %>%
  left_join(richness_coastal_AM, by = 'genus') %>%
  left_join(relabund_upland_AM, by = 'genus') %>%
  left_join(richness_upland_AM, by = 'genus') %>%
  select(Genus = genus, 
         'Relative Abundance (All)' = rel_abund,
         'OTU Richness (All)' = n, 
         'Relative Abundance (Coastal)' = rel_abund_coastal,
         'OTU Richness (Coastal)' = n_coastal, 
         'Relative Abundance (Upland)' = rel_abund_upland,
         'OTU Richness (Upland)' = n_upland) %>%
  arrange('Relative Abundance') %>% 
  replace(is.na(.), 0) %>%
  write_csv('output/genera_AM.csv')

### AM relative abundance plot ###

# Combine genera Incertae sedis for stacked bar plots
incertae_sedis_AM <- relabund_AM %>%
  filter(str_ends(genus, "gen_Incertae_sedis")) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  mutate(genus = 'Genera incertae sedis') %>%
  select(genus, rel_abund) %>%
  glimpse()

# Define stacked bar plot order:
# I want genera in ascending order of relative abundance with 'other' and 
# 'Genera incertae sedis' at the end of the bars
stacked_orderAM <- relabund_AM %>%
  filter(!str_ends(genus, "gen_Incertae_sedis")) %>%
  arrange(rel_abund) %>%
  bind_rows(incertae_sedis_AM, .) %>%
  select(-rel_abund) %>%
  mutate(stacked_order = c('a', 'b', 'c', 'd', 'e', 'f', 'g')) %>%
  print(n = Inf)

# Genera incertae sedis by ecotype
isAM <-
  bind_rows(
    relabund_coastal_AM %>% mutate(ecotype = "Coastal",
                                   rel_abund = rel_abund_coastal),
    relabund_upland_AM %>% mutate(ecotype = "Upland",
                                  rel_abund = rel_abund_upland)
  ) %>%
  select(ecotype, genus, rel_abund) %>%
  filter(str_ends(genus, "gen_Incertae_sedis")) %>%
  mutate(genus = 'Genera incertae sedis') %>%
  glimpse()

# Define palette
palAM <- c((paletteer_d("ggthemes::Tableau_10", 6)), 'light grey') %>%
  rev()
# Plot
relabund_plot_AM <-
  bind_rows(
    relabund_coastal_AM %>% mutate(ecotype = "Coastal",
                                   rel_abund = rel_abund_coastal),
    relabund_upland_AM %>% mutate(ecotype = "Upland",
                                  rel_abund = rel_abund_upland)
  ) %>%
  select(ecotype, genus, rel_abund) %>%
  mutate(genus = ifelse(str_ends(genus, "gen_Incertae_sedis"),
                        'Genera incertae sedis', genus)) %>%
  inner_join(stacked_orderAM, by = 'genus') %>%
  ggplot(aes(ecotype, rel_abund/100, fill = stacked_order)) + 
  geom_bar(stat = 'identity') +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
        axis.text.x = element_text(hjust = 0.5, colour = 'black', size = rel(1.1)),
        axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black", size = rel(1.1)),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.title = element_blank(),
        legend.text = element_text(size = rel(0.8)),
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5, size = rel(0.9))
  ) +
  guides(fill = guide_legend(reverse = T, nrow = 2, byrow = F)) +
  coord_flip() +
  labs(title = 'Relative abundance', x = 'Arbuscular mycorrhizal', y= NULL) +
  scale_x_discrete(limits = c('Coastal', 'Upland')) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(
    values = palAM,
    limits = c('a', 'b', 'c', 'd', 'e', 'f', 'g'),
    labels = c('Glomeraceae [10]\nincertae sedis', 'Oehlia [2]',
               'Rhizoglomus [6]', 'Dominikia [8]', 'Glomus [21]',
               'Claroideoglomus [6]', 'Rhizophagus [15]')
    )
relabund_plot_AM

#### EM: Relative abundance and richness ################################

# EM: Relative abundance of genera
relabund_EM <- inner_join(otuEM, taxaEM, by = 'OTU_ID') %>%
  select(genus, count) %>%
  mutate(rel_abund = relabund(count) * 100) %>%
  group_by(genus) %>%
  summarise(rel_abund = round(sum(rel_abund), digits = 1))
sum(relabund_EM$rel_abund)
# EM: Richness of genera
richness_EM <- inner_join(otuEM, taxaEM, by = 'OTU_ID') %>%
  filter(count != 0) %>%
  distinct(OTU_ID, .keep_all = TRUE) %>%
  select(OTU_ID, genus) %>%
  count(genus)
sum(richness_EM$n)

# EM: Relative abundance coastal
relabund_coastal_EM <-
  inner_join(otuEM, taxaEM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  select(ecotype, genus, count) %>%
  filter(ecotype == 'Coastal') %>%
  mutate(rel_abund = relabund(count) * 100) %>%
  group_by(genus) %>%
  summarise(rel_abund_coastal = round(sum(rel_abund), digits = 1)) %>%
  print(n = Inf)
sum(relabund_coastal_EM$rel_abund_coastal)
# EM: Richness coastal
richness_coastal_EM <- inner_join(otuEM, taxaEM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  filter(ecotype == 'Coastal') %>%
  filter(count != 0) %>%
  distinct(OTU_ID, .keep_all = TRUE) %>%
  select(OTU_ID, genus) %>%
  count(genus) %>%
  rename(n_coastal = n) %>%
  print(n = Inf)

# EM: Relative abundance upland
relabund_upland_EM <-
  inner_join(otuEM, taxaEM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  select(ecotype, genus, count) %>%
  filter(ecotype == 'Upland') %>%
  mutate(rel_abund = relabund(count) * 100) %>%
  group_by(genus) %>%
  summarise(rel_abund_upland = round(sum(rel_abund), digits = 1)) %>%
  print(n = Inf)
sum(relabund_upland_EM$rel_abund_upland)
# EM: Richness upland
richness_upland_EM <- inner_join(otuEM, taxaEM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  filter(ecotype == 'Upland') %>%
  filter(count != 0) %>%
  distinct(OTU_ID, .keep_all = TRUE) %>%
  select(OTU_ID, genus) %>%
  count(genus) %>%
  rename(n_upland = n) %>%
  print(n = Inf)

# Export relative abundance and richness table
left_join(relabund_EM, richness_EM, by = 'genus') %>%
  left_join(relabund_coastal_EM, by = 'genus') %>%
  left_join(richness_coastal_EM, by = 'genus') %>%
  left_join(relabund_upland_EM, by = 'genus') %>%
  left_join(richness_upland_EM, by = 'genus') %>%
  select(Genus = genus, 
         'Relative Abundance (All)' = rel_abund,
         'OTU Richness (All)' = n, 
         'Relative Abundance (Coastal)' = rel_abund_coastal,
         'OTU Richness (Coastal)' = n_coastal, 
         'Relative Abundance (Upland)' = rel_abund_upland,
         'OTU Richness (Upland)' = n_upland) %>%
  arrange('Relative Abundance') %>% 
  replace(is.na(.), 0) %>%
  write_csv('output/genera_EM.csv')

### EM relative abundance plot ###

# Combine low abundance taxaon for stacked bar plots
other <- relabund_EM %>%
  filter(rel_abund < 2) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  mutate(genus = 'other')
other_filter <- relabund_EM %>%
  filter(rel_abund < 2)

# Define stacked bar plot order:
# I want genera in ascending order of relative abundance with 'other' and 
# 'Genera incertae sedis' at the end of the bars
stacked_orderEM <- relabund_EM %>%
  filter(rel_abund > 2) %>%
  arrange(rel_abund) %>%
  bind_rows(other, .) %>%
  select(-rel_abund) %>%
  mutate(stacked_order = c('a', 'b', 'c', 'd', 'e', 'f', 'g')) %>%
  print(n = Inf)

# Low relative abundance by ecotype
other_EM <-
  bind_rows(
    relabund_coastal_EM %>%
      filter(genus %in% other_filter$genus) %>% 
      summarise(rel_abund = sum(rel_abund_coastal)) %>% 
      mutate(genus = 'other') %>% mutate(ecotype = 'Coastal'),
    relabund_upland_EM %>% 
      filter(genus %in% other_filter$genus) %>%
      summarise(rel_abund = sum(rel_abund_upland)) %>%
      mutate(genus = 'other') %>% mutate(ecotype = 'Upland')
  ) %>%
  select(ecotype, genus, rel_abund)

# Define pallatte
palEM = c('#C03728', '#919C4C', '#FD8F24', '#F5C04A', '#6E9FAF', 
          '#E68C7C', 'light grey') %>%
  rev()

# Plot 
relabund_plot_EM <-
  bind_rows(
    relabund_coastal_EM %>% mutate(rel_abund = rel_abund_coastal) %>%
      mutate(ecotype = 'Coastal'),
    relabund_upland_EM %>% mutate(rel_abund = rel_abund_upland) %>%
      mutate(ecotype = 'Upland')
  ) %>%
  select(ecotype, genus, rel_abund) %>%
  filter(!genus %in% other_filter$genus) %>%
  bind_rows(other_EM) %>%
  inner_join(stacked_orderEM, by = 'genus') %>%
  ggplot(aes(ecotype, rel_abund/100, fill = stacked_order)) + 
  geom_bar(stat = 'identity') +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
        axis.text.x = element_text(hjust = 0.5, colour = 'black', size = rel(1.1)),
        axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black", size = rel(1.1)),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.title = element_blank(),
        legend.text = element_text(size = rel(0.8)),
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5, size = rel(0.9))
  ) +
  guides(fill = guide_legend(reverse = T, nrow = 2, byrow = F)) +
  coord_flip() +
  labs(x = 'Ectomycorrhizal', y= NULL) +
  scale_x_discrete(limits = c('Coastal', 'Upland')) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(
    values = palEM,
    limits = c('a', 'b', 'c', 'd', 'e', 'f', 'g'),
    labels = c('other [19]\n(< 2%)', 'Cenococcum [2]', 'Inocybe [6]',
               'Thelephora [6]', 'Sebacina [10]', 'Ruhlandiella [9]',
               'Tomentella [28]'))
relabund_plot_EM

#### Export figures ###########################################################

relabund_plot_AM + relabund_plot_EM + plot_layout(ncol = 1)
ggsave('output/rel_abund.pdf', width = 15, height = 15, unit = 'cm')
ggsave('output/rel_abund.jpg', width = 15, height = 15, unit = 'cm')
ggsave('output/rel_abund.tiff', width = 15, height = 15, unit = 'cm')

#### EM exploration type: Relative abundance and richness #####################

# EM: Relative abundance of exploration type
relabund_EM_explore <- inner_join(otuEM, taxaEM, by = 'OTU_ID') %>%
  select(exploration_type = ectomycorrhiza_exploration_type, count) %>%
  mutate(rel_abund = relabund(count) * 100) %>%
  group_by(exploration_type) %>%
  summarise(rel_abund = round(sum(rel_abund), digits = 1))
sum(relabund_EM_explore$rel_abund)
# EM: Richness of exploration type
richness_EM_explore <- inner_join(otuEM, taxaEM, by = 'OTU_ID') %>%
  filter(count != 0) %>%
  distinct(OTU_ID, .keep_all = TRUE) %>%
  select(OTU_ID, exploration_type = ectomycorrhiza_exploration_type) %>%
  count(exploration_type)
sum(richness_EM_explore$n)

# EM: Relative abundance coastal
relabund_coastal_EM_explore <-
  inner_join(otuEM, taxaEM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  select(ecotype, exploration_type = ectomycorrhiza_exploration_type, count) %>%
  filter(ecotype == 'Coastal') %>%
  mutate(rel_abund = relabund(count) * 100) %>%
  group_by(exploration_type) %>%
  summarise(rel_abund_coastal = round(sum(rel_abund), digits = 1)) %>%
  print(n = Inf)
sum(relabund_coastal_EM_explore$rel_abund_coastal)
# EM: Richness coastal
richness_coastal_EM_explore <- inner_join(otuEM, taxaEM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  filter(ecotype == 'Coastal') %>%
  filter(count != 0) %>%
  distinct(OTU_ID, .keep_all = TRUE) %>%
  select(OTU_ID, exploration_type = ectomycorrhiza_exploration_type) %>%
  count(exploration_type) %>%
  rename(n_coastal = n) %>%
  print(n = Inf)

# EM: Relative abundance upland
relabund_upland_EM_explore <-
  inner_join(otuEM, taxaEM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  select(ecotype, exploration_type = ectomycorrhiza_exploration_type, count) %>%
  filter(ecotype == 'Upland') %>%
  mutate(rel_abund = relabund(count) * 100) %>%
  group_by(exploration_type) %>%
  summarise(rel_abund_upland = round(sum(rel_abund), digits = 1)) %>%
  print(n = Inf)
sum(relabund_upland_EM_explore$rel_abund_upland)
# EM: Richness upland
richness_upland_EM_explore <- inner_join(otuEM, taxaEM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  filter(ecotype == 'Upland') %>%
  filter(count != 0) %>%
  distinct(OTU_ID, .keep_all = TRUE) %>%
  select(OTU_ID, exploration_type = ectomycorrhiza_exploration_type) %>%
  count(exploration_type) %>%
  rename(n_upland = n) %>%
  print(n = Inf)

# Export relative abundance and richness table
left_join(relabund_EM_explore, richness_EM_explore, by = 'exploration_type') %>%
  left_join(relabund_coastal_EM_explore, by = 'exploration_type') %>%
  left_join(richness_coastal_EM_explore, by = 'exploration_type') %>%
  left_join(relabund_upland_EM_explore, by = 'exploration_type') %>%
  left_join(richness_upland_EM_explore, by = 'exploration_type') %>%
  select(ectomycorrhiza_exploration_type = exploration_type, 
         'Relative Abundance (All)' = rel_abund,
         'OTU Richness (All)' = n, 
         'Relative Abundance (Coastal)' = rel_abund_coastal,
         'OTU Richness (Coastal)' = n_coastal, 
         'Relative Abundance (Upland)' = rel_abund_upland,
         'OTU Richness (Upland)' = n_upland) %>%
  arrange('Relative Abundance') %>% 
  replace(is.na(.), 0) %>%
  write_csv('output/exploration_type_EM.csv')

### EM relative abundance plot ###

other <- relabund_EM_explore %>%
  filter(rel_abund < 1.5) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  mutate(exploration_type = 'other')
other_filter <- relabund_EM_explore %>%
  filter(rel_abund < 1.5)

# Define stacked bar plot order:
# I want exploration type in ascending order of relative abundance with 'other' and 
# 'exploration type incertae sedis' at the end of the bars
stacked_orderEM <- relabund_EM_explore %>%
  filter(rel_abund > 1.5) %>%
  arrange(rel_abund) %>%
  bind_rows(other, .) %>%
  select(-rel_abund) %>%
  mutate(stacked_order = c('a', 'b', 'c', 'd', 'e')) %>%
  print(n = Inf)

# Low relative abundance taxa
other_EM <-
  bind_rows(
    relabund_coastal_EM_explore %>%
      filter(exploration_type %in% other_filter$exploration_type) %>% 
      summarise(rel_abund = sum(rel_abund_coastal)) %>% 
      mutate(exploration_type = 'other') %>% mutate(ecotype = 'Coastal'),
    relabund_upland_EM_explore %>% 
      filter(exploration_type %in% other_filter$exploration_type) %>%
      summarise(rel_abund = sum(rel_abund_upland)) %>%
      mutate(exploration_type = 'other') %>% mutate(ecotype = 'Upland')
  ) %>%
  select(ecotype, exploration_type, rel_abund)

# Define palette
palEM <- c('#0073C2', '#EFC000', '#CD534C', '#003C67', 'light grey') %>%
  rev()

# Plot 
relabund_plot_EM_explore <-
  bind_rows(
    relabund_coastal_EM_explore %>% mutate(rel_abund = rel_abund_coastal) %>%
      mutate(ecotype = 'Coastal'),
    relabund_upland_EM_explore %>% mutate(rel_abund = rel_abund_upland) %>%
      mutate(ecotype = 'Upland')
  ) %>%
  select(ecotype, exploration_type, rel_abund) %>%
  filter(!exploration_type %in% other_filter$exploration_type) %>%
  bind_rows(other) %>%
  inner_join(stacked_orderEM, by = 'exploration_type') %>%
  ggplot(aes(ecotype, rel_abund/100, fill = stacked_order)) + 
  geom_bar(stat = 'identity') +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
        axis.text.x = element_text(hjust = 0.5, colour = 'black', size = rel(1)),
        axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black", size = rel(1)),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.title = element_blank(),
        legend.text = element_text(size = rel(0.8)),
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5, size = rel(0.9))
  ) +
  guides(fill = guide_legend(reverse = T, nrow = 2, byrow = F)) +
  coord_flip() +
  labs(x = 'Ecotype', y= "Relative abundance") +
  scale_x_discrete(limits = c('Coastal', 'Upland')) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = palEM,
                    limits = c('a', 'b', 'c', 'd', 'e'),
                    labels = c('other (< 2%) [5]', 'Contact [3]',
                               'Short-distance, delicate [23]',
                               'Short-distance, coarse [11]',
                               'Medium-distance, smooth [38]'))
relabund_plot_EM_explore

ggsave('output/rel_abund_explore.pdf', width = 15, height = 10, unit = 'cm')
ggsave('output/rel_abund_explore.jpg', width = 15, height = 10, unit = 'cm')
ggsave('output/rel_abund_explore.tiff', width = 15, height = 10, unit = 'cm')
