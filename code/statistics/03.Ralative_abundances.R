
# Script: Calculate and plot relative abundances of mycorrhizal genera across
# tea tree ecotypes

# Required packages and functions
require(patchwork)
require(paletteer)
require(tidyverse)
source("code/statistics/functions.R")

# Metadata
data <- read.csv('data/statistics/metadata.csv', 
                 header = TRUE, stringsAsFactors = TRUE) %>%
  select(sample, ecotype) %>%
  glimpse()

# OTUs
otu_AM = t(read.csv('data/statistics/otu_AM.csv',
                   h = TRUE, row.names = 'OTU_ID')) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample') %>%
  as_tibble() %>%
  pivot_longer(-sample, names_to = 'OTU_ID', values_to = 'count')
# AM taxa
taxa_AM = read.csv('data/statistics/taxa_AM.csv', header = TRUE)

# EM OTUs
otu_ECM = t(read.csv('data/statistics/otu_ECM.csv',
                   h = TRUE, row.names = 'OTU_ID')) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample') %>%
  as_tibble() %>%
  pivot_longer(-sample, names_to = 'OTU_ID', values_to = 'count')
# EM taxa
taxa_ECM = read.csv('data/statistics/taxa_ECM.csv', header = TRUE)

#### AM: Relative abundance and richness ################################

# AM: Relative abundance of genera
relabund_AM <- inner_join(otu_AM, taxa_AM, by = 'OTU_ID') %>%
  select(genus, count) %>%
  mutate(rel_abund = relabund(count) * 100) %>%
  group_by(genus) %>%
  summarise(rel_abund = round(sum(rel_abund), digits = 1))
sum(relabund_AM$rel_abund)
# AM: Richness of genera
richness_AM <- inner_join(otu_AM, taxa_AM, by = 'OTU_ID') %>%
  filter(count != 0) %>%
  distinct(OTU_ID, .keep_all = TRUE) %>%
  select(OTU_ID, genus) %>%
  count(genus) %>%
  print(n = Inf)
sum(richness_AM$n)

# AM: Relative abundance coastal
relabund_coastal_AM <-
  inner_join(otu_AM, taxa_AM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  select(ecotype, genus, count) %>%
  filter(ecotype == 'Coastal') %>%
  mutate(rel_abund = relabund(count) * 100) %>%
  group_by(genus) %>%
  summarise(rel_abund_coastal = round(sum(rel_abund), digits = 1)) %>%
  print(n = Inf)
sum(relabund_coastal_AM$rel_abund_coastal)
# AM: Richness coastal
richness_coastal_AM <- inner_join(otu_AM, taxa_AM, by = 'OTU_ID') %>%
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
  inner_join(otu_AM, taxa_AM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  select(ecotype, genus, count) %>%
  filter(ecotype == 'Upland') %>%
  mutate(rel_abund = relabund(count) * 100) %>%
  group_by(genus) %>%
  summarise(rel_abund_upland = round(sum(rel_abund), digits = 1)) %>%
  print(n = Inf)
sum(relabund_upland_AM$rel_abund_upland)
# AM: Richness upland
richness_upland_AM <- inner_join(otu_AM, taxa_AM, by = 'OTU_ID') %>%
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

# Combine low abundance taxaon for stacked bar plots
other <- relabund_AM %>%
  filter(
    !str_ends(genus, "gen_Incertae_sedis"),
    rel_abund < 2
    ) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  mutate(genus = 'other')
other_filter_AM <- relabund_AM %>%
  filter(
    rel_abund < 2,
    !str_ends(genus, "gen_Incertae_sedis")
    )

# Low relative abundance by ecotype
other_AM <-
  bind_rows(
    relabund_coastal_AM %>%
      filter(genus %in% other_filter_AM$genus) %>% 
      summarise(rel_abund = sum(rel_abund_coastal)) %>% 
      mutate(genus = 'other') %>% mutate(ecotype = 'Coastal'),
    relabund_upland_AM %>% 
      filter(genus %in% other_filter_AM$genus) %>%
      summarise(rel_abund = sum(rel_abund_upland)) %>%
      mutate(genus = 'other') %>% mutate(ecotype = 'Upland')
  ) %>%
  select(ecotype, genus, rel_abund)

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
  filter(
    !str_ends(genus, "gen_Incertae_sedis"),
    !rel_abund < 2
    ) %>%
  arrange(rel_abund) %>%
  bind_rows(
    incertae_sedis_AM,
    other,
    .) %>%
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
palAM <- c((paletteer_d("ggthemes::Tableau_10", 5)), "light grey", "#636363") %>%
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
  filter(!genus %in% other_filter_AM$genus) %>%
  bind_rows(other_AM) %>%
  mutate(
    genus = ifelse(str_ends(genus, "Incertae_sedis"),
                   'Genera incertae sedis', genus)) %>%
  group_by(ecotype, genus) %>%
  summarise(
    rel_abund = sum(rel_abund)
  ) %>%
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
    labels = c('Glomeromycota [18]\nincertae sedis', 'other [8]\n(< 2%)',
               'Dominikia [5]', 'Glomus [17]', 'Claroideoglomus [6]',
               'Rhizophagus [15]', 'Entrophospora [3]')
    )

#### EM: Relative abundance and richness ################################

# EM: Relative abundance of genera
relabund_ECM <- inner_join(otu_ECM, taxa_ECM, by = 'OTU_ID') %>%
  select(genus, count) %>%
  mutate(rel_abund = relabund(count) * 100) %>%
  group_by(genus) %>%
  summarise(rel_abund = round(sum(rel_abund), digits = 1))
sum(relabund_ECM$rel_abund)
# EM: Richness of genera
richness_ECM <- inner_join(otu_ECM, taxa_ECM, by = 'OTU_ID') %>%
  filter(count != 0) %>%
  distinct(OTU_ID, .keep_all = TRUE) %>%
  select(OTU_ID, genus) %>%
  count(genus) %>%
  print(n = Inf)
sum(richness_ECM$n)

# EM: Relative abundance coastal
relabund_coastal_ECM <-
  inner_join(otu_ECM, taxa_ECM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  select(ecotype, genus, count) %>%
  filter(ecotype == 'Coastal') %>%
  mutate(rel_abund = relabund(count) * 100) %>%
  group_by(genus) %>%
  summarise(rel_abund_coastal = round(sum(rel_abund), digits = 1)) %>%
  print(n = Inf)
sum(relabund_coastal_ECM$rel_abund_coastal)
# EM: Richness coastal
richness_coastal_ECM <- inner_join(otu_ECM, taxa_ECM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  filter(ecotype == 'Coastal') %>%
  filter(count != 0) %>%
  distinct(OTU_ID, .keep_all = TRUE) %>%
  select(OTU_ID, genus) %>%
  count(genus) %>%
  rename(n_coastal = n) %>%
  print(n = Inf)

# EM: Relative abundance upland
relabund_upland_ECM <-
  inner_join(otu_ECM, taxa_ECM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  select(ecotype, genus, count) %>%
  filter(ecotype == 'Upland') %>%
  mutate(rel_abund = relabund(count) * 100) %>%
  group_by(genus) %>%
  summarise(rel_abund_upland = round(sum(rel_abund), digits = 1)) %>%
  print(n = Inf)
sum(relabund_upland_ECM$rel_abund_upland)
# EM: Richness upland
richness_upland_ECM <- inner_join(otu_ECM, taxa_ECM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  filter(ecotype == 'Upland') %>%
  filter(count != 0) %>%
  distinct(OTU_ID, .keep_all = TRUE) %>%
  select(OTU_ID, genus) %>%
  count(genus) %>%
  rename(n_upland = n) %>%
  print(n = Inf)

# Export relative abundance and richness table
left_join(relabund_ECM, richness_ECM, by = 'genus') %>%
  left_join(relabund_coastal_ECM, by = 'genus') %>%
  left_join(richness_coastal_ECM, by = 'genus') %>%
  left_join(relabund_upland_ECM, by = 'genus') %>%
  left_join(richness_upland_ECM, by = 'genus') %>%
  select(Genus = genus, 
         'Relative Abundance (All)' = rel_abund,
         'OTU Richness (All)' = n, 
         'Relative Abundance (Coastal)' = rel_abund_coastal,
         'OTU Richness (Coastal)' = n_coastal, 
         'Relative Abundance (Upland)' = rel_abund_upland,
         'OTU Richness (Upland)' = n_upland) %>%
  arrange('Relative Abundance') %>% 
  replace(is.na(.), 0) %>%
  write_csv('output/genera_ECM.csv')

### EM relative abundance plot ###

# Combine low abundance taxaon for stacked bar plots
other <- relabund_ECM %>%
  filter(rel_abund < 2) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  mutate(genus = 'other')
other_filter <- relabund_ECM %>%
  filter(rel_abund < 2)

# Define stacked bar plot order:
# I want genera in ascending order of relative abundance with 'other' and 
# 'Genera incertae sedis' at the end of the bars
stacked_orderEM <- relabund_ECM %>%
  filter(rel_abund > 2) %>%
  arrange(rel_abund) %>%
  bind_rows(other, .) %>%
  select(-rel_abund) %>%
  mutate(stacked_order = c('a', 'b', 'c', 'd', 'e', 'f', 'g')) %>%
  print(n = Inf)

# Low relative abundance by ecotype
other_ECM <-
  bind_rows(
    relabund_coastal_ECM %>%
      filter(genus %in% other_filter$genus) %>% 
      summarise(rel_abund = sum(rel_abund_coastal)) %>% 
      mutate(genus = 'other') %>% mutate(ecotype = 'Coastal'),
    relabund_upland_ECM %>% 
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
relabund_plot_ECM <-
  bind_rows(
    relabund_coastal_ECM %>% mutate(rel_abund = rel_abund_coastal) %>%
      mutate(ecotype = 'Coastal'),
    relabund_upland_ECM %>% mutate(rel_abund = rel_abund_upland) %>%
      mutate(ecotype = 'Upland')
  ) %>%
  select(ecotype, genus, rel_abund) %>%
  filter(!genus %in% other_filter$genus) %>%
  bind_rows(other_ECM) %>%
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
    labels = c('other [28]\n(< 2%)', 'Cenococcum [3]', 'Inocybe [8]',
               'Thelephora [6]', 'Sebacina [15]', 'Ruhlandiella [10]',
               'Tomentella [36]'))

#### Export figures ###########################################################

relabund_plot_AM + relabund_plot_ECM + plot_layout(ncol = 1)
ggsave('output/rel_abund.pdf', width = 15, height = 15, unit = 'cm')
ggsave('output/rel_abund.jpg', width = 15, height = 15, unit = 'cm')
ggsave('output/rel_abund.tiff', width = 15, height = 15, unit = 'cm')

#### EM exploration type: Relative abundance and richness #####################

# EM: Relative abundance of exploration type
relabund_ECM_explore <- inner_join(otu_ECM, taxa_ECM, by = 'OTU_ID') %>%
  select(exploration_type = ectomycorrhiza_exploration_type, count) %>%
  mutate(rel_abund = relabund(count) * 100) %>%
  group_by(exploration_type) %>%
  summarise(rel_abund = round(sum(rel_abund), digits = 1))
sum(relabund_ECM_explore$rel_abund)
# EM: Richness of exploration type
richness_ECM_explore <- inner_join(otu_ECM, taxa_ECM, by = 'OTU_ID') %>%
  filter(count != 0) %>%
  distinct(OTU_ID, .keep_all = TRUE) %>%
  select(OTU_ID, exploration_type = ectomycorrhiza_exploration_type) %>%
  count(exploration_type)
sum(richness_ECM_explore$n)

# EM: Relative abundance coastal
relabund_coastal_ECM_explore <-
  inner_join(otu_ECM, taxa_ECM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  select(ecotype, exploration_type = ectomycorrhiza_exploration_type, count) %>%
  filter(ecotype == 'Coastal') %>%
  mutate(rel_abund = relabund(count) * 100) %>%
  group_by(exploration_type) %>%
  summarise(rel_abund_coastal = round(sum(rel_abund), digits = 1)) %>%
  print(n = Inf)
sum(relabund_coastal_ECM_explore$rel_abund_coastal)
# EM: Richness coastal
richness_coastal_ECM_explore <- inner_join(otu_ECM, taxa_ECM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  filter(ecotype == 'Coastal') %>%
  filter(count != 0) %>%
  distinct(OTU_ID, .keep_all = TRUE) %>%
  select(OTU_ID, exploration_type = ectomycorrhiza_exploration_type) %>%
  count(exploration_type) %>%
  rename(n_coastal = n) %>%
  print(n = Inf)

# EM: Relative abundance upland
relabund_upland_ECM_explore <-
  inner_join(otu_ECM, taxa_ECM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  select(ecotype, exploration_type = ectomycorrhiza_exploration_type, count) %>%
  filter(ecotype == 'Upland') %>%
  mutate(rel_abund = relabund(count) * 100) %>%
  group_by(exploration_type) %>%
  summarise(rel_abund_upland = round(sum(rel_abund), digits = 1)) %>%
  print(n = Inf)
sum(relabund_upland_ECM_explore$rel_abund_upland)
# EM: Richness upland
richness_upland_ECM_explore <- inner_join(otu_ECM, taxa_ECM, by = 'OTU_ID') %>%
  inner_join(data, by = 'sample') %>%
  filter(ecotype == 'Upland') %>%
  filter(count != 0) %>%
  distinct(OTU_ID, .keep_all = TRUE) %>%
  select(OTU_ID, exploration_type = ectomycorrhiza_exploration_type) %>%
  count(exploration_type) %>%
  rename(n_upland = n) %>%
  print(n = Inf)

# Export relative abundance and richness table
left_join(relabund_ECM_explore, richness_ECM_explore, by = 'exploration_type') %>%
  left_join(relabund_coastal_ECM_explore, by = 'exploration_type') %>%
  left_join(richness_coastal_ECM_explore, by = 'exploration_type') %>%
  left_join(relabund_upland_ECM_explore, by = 'exploration_type') %>%
  left_join(richness_upland_ECM_explore, by = 'exploration_type') %>%
  select(ectomycorrhiza_exploration_type = exploration_type, 
         'Relative Abundance (All)' = rel_abund,
         'OTU Richness (All)' = n, 
         'Relative Abundance (Coastal)' = rel_abund_coastal,
         'OTU Richness (Coastal)' = n_coastal, 
         'Relative Abundance (Upland)' = rel_abund_upland,
         'OTU Richness (Upland)' = n_upland) %>%
  arrange('Relative Abundance') %>% 
  replace(is.na(.), 0) %>%
  write_csv('output/exploration_type_ECM.csv')

### EM relative abundance plot ###

other <- relabund_ECM_explore %>%
  filter(rel_abund < 2) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  mutate(exploration_type = 'other')
other_filter <- relabund_ECM_explore %>%
  filter(rel_abund < 2)

# Define stacked bar plot order:
# I want exploration type in ascending order of relative abundance with 'other' and 
# 'exploration type incertae sedis' at the end of the bars
stacked_orderEM <- relabund_ECM_explore %>%
  filter(rel_abund > 2) %>%
  arrange(rel_abund) %>%
  bind_rows(other, .) %>%
  select(-rel_abund) %>%
  mutate(stacked_order = c('a', 'b', 'c', 'd')) %>%
  print(n = Inf)

# Low relative abundance taxa
other_ECM <-
  bind_rows(
    relabund_coastal_ECM_explore %>%
      filter(exploration_type %in% other_filter$exploration_type) %>% 
      summarise(rel_abund = sum(rel_abund_coastal)) %>% 
      mutate(exploration_type = 'other') %>% mutate(ecotype = 'Coastal'),
    relabund_upland_ECM_explore %>% 
      filter(exploration_type %in% other_filter$exploration_type) %>%
      summarise(rel_abund = sum(rel_abund_upland)) %>%
      mutate(exploration_type = 'other') %>% mutate(ecotype = 'Upland')
  ) %>%
  select(ecotype, exploration_type, rel_abund)

# Define palette
palEM <- c('#0073C2', '#EFC000', '#CD534C', 'light grey') %>%
  rev()

# Plot 
relabund_plot_ECM_explore <-
  bind_rows(
    relabund_coastal_ECM_explore %>% mutate(rel_abund = rel_abund_coastal) %>%
      mutate(ecotype = 'Coastal'),
    relabund_upland_ECM_explore %>% mutate(rel_abund = rel_abund_upland) %>%
      mutate(ecotype = 'Upland')
  ) %>%
  select(ecotype, exploration_type, rel_abund) %>%
  filter(!exploration_type %in% other_filter$exploration_type) %>%
  bind_rows(other_ECM) %>%
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
                    limits = c('a', 'b', 'c', 'd'),
                    labels = c('other (< 2%) [14]',
                               'Short-distance, delicate [32]',
                               'Short-distance, coarse [13]',
                               'Medium-distance, smooth [47]'))
relabund_plot_ECM_explore

ggsave('output/rel_abund_explore.pdf', width = 15, height = 10, unit = 'cm')
ggsave('output/rel_abund_explore.jpg', width = 15, height = 10, unit = 'cm')
ggsave('output/rel_abund_explore.tiff', width = 15, height = 10, unit = 'cm')

# Clear all objects from the environment
rm(list = ls())
