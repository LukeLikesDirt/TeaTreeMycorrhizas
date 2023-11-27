
# Script: Assess mycorrhizal fungal beta diversity across ecotypes

# Required packages and functions
source("code/statistics/functions.R")
require(patchwork)
require(vegan)
require(tidyverse)

# Metadata
data = read.csv('data/statistics/metadata.csv', 
                header = TRUE, stringsAsFactors = TRUE) %>%
  select(ecotype, sample) %>%
  glimpse()

# OTUs
otuAM = t(read.csv('data/statistics/otuAM.csv',
                   header = TRUE, row.names = 'OTU_ID')) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample') %>%
  as_tibble() %>%
  pivot_longer(-sample, names_to = 'OTU_ID', values_to = 'value') %>%
  glimpse()
otuEM = t(read.csv('data/statistics/otuEM.csv',
                   header = TRUE, row.names = 'OTU_ID')) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample') %>%
  as_tibble() %>%
  pivot_longer(-sample, names_to = 'OTU_ID', values_to = 'value')
# Taxa
taxaAM = read.csv('data/statistics/taxaAM.csv', header = TRUE)
taxaEM = read.csv('data/statistics/taxaEM.csv', header = TRUE)


#### (1) Calculate relative abundances and remove singltons ####################

##### (1a) AM Fungi ######
relAM <- inner_join(data, otuAM, by = 'sample', multiple = 'all') %>%
  as_tibble() %>%
  select(c(sample, OTU_ID, value)) %>%
  group_by(sample, OTU_ID) %>%
  # Calculate prevalence and remove genera that occur once
  summarise(count = sum(value), .groups = 'drop') %>%
  # Remove zeros
  filter(count != 0) %>%
  # Calculate prevalence
  group_by(OTU_ID) %>%
  mutate(prevelance = sum(count/count)) %>%
  ungroup() %>%
  # Remove taxa that occur once
  filter(prevelance != 1) %>%
  select(-prevelance) %>%
  ungroup() %>%
  # Calculate relative abundances to account for differences in sequencing depth
  group_by(sample) %>%
  mutate(rel_abund = relabund(count)) %>%
  # Remove counts
  select(-count) %>%
  # Fromat OTU table of relative abundance at occur in at least two samples
  pivot_wider(names_from = OTU_ID, values_from = rel_abund) %>%
  replace(is.na(.), 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = 'sample')

##### (1b) EM Fungi #####
relEM <- inner_join(data, otuEM, by = 'sample', multiple = 'all') %>%
  as_tibble() %>%
  select(c(sample, OTU_ID, value)) %>%
  group_by(sample, OTU_ID) %>%
  ### Calculate prevalence and remove OTUs that occur once ###
  summarise(count = sum(value), .groups = 'drop') %>%
  # Remove zeros
  filter(count != 0) %>%
  # Calculate prevalence
  group_by(OTU_ID) %>%
  mutate(prevelance = sum(count/count)) %>%
  ungroup() %>%
  # Remove taxa that occur once
  filter(prevelance != 1) %>%
  select(-prevelance) %>%
  ungroup() %>%
  # Calculate relative abundances to account for differences in sequencing depth
  group_by(sample)%>%
  mutate(rel_abund = relabund(count)) %>%
  # Remove counts
  select(-count) %>%
  # Fromat OTU table of relative abundance at occur in at least two samples
  pivot_wider(names_from = OTU_ID, values_from = rel_abund) %>%
  replace(is.na(.), 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = 'sample')

#### (2) Ordinations and PERMANOVAs ############################################

##### (2a) AM Fungi #####

# Calculate Jaccard distances and NMDS points
distJ_AM = vegdist(relAM, method = 'jaccard', binary = T)
set.seed(1986)
nmdsJ_AM = metaMDS(distJ_AM, k = 3)
nmdsJ_AM
# Stress = 0.059

# Calculate Bray-Curtis distances and NMDS points
distB_AM = vegdist(relAM, method = 'bray')
set.seed(1986)
nmdsB_AM = metaMDS(distB_AM, k = 3)
set.seed(1986)
nmdsB_AM = metaMDS(distB_AM, k = 3, previous.best = nmdsB_AM)
nmdsB_AM
# Stress = 0.067

# NOTE: AM Fungi were not detected at all sites; generate data file for sites 
# with AMF
AMfilter = row.names(as.matrix(distJ_AM))
Sample_ID = data$sample
dataAM = data[(Sample_ID %in% AMfilter),]
glimpse(dataAM)

# NOTE: adonis2 does not adjust significant values when assessing random 
# effects. So, unlike the alpha diversity analysis, I have not included site
# as a random effect

# Differences between ecotypes using Jaccard distances
set.seed(1986)
adJ_AM = adonis2(distJ_AM ~ ecotype, data = dataAM)
adJ_AM
bdJ_AM = betadisper(distJ_AM, dataAM$ecotype)
set.seed(1986)
bdJp_AM = permutest(bdJ_AM, permutations = 999)
bdJp_AM

# Differences between ecotypes using Bray-Curtis distances
set.seed(1986)
adB_AM = adonis2(distB_AM ~ ecotype, data = dataAM)
adB_AM
bdB_AM = betadisper(distB_AM, dataAM$ecotype)
set.seed(1986)
bdBp_AM = permutest(bdB_AM, permutations = 999)
bdBp_AM

##### (2b) EM Fungi #####

# Calculate Jaccard distances and NMDS points
distJ_EM = vegdist(relEM, method = 'jaccard', binary = TRUE)
set.seed(1986)
nmdsJ_EM = metaMDS(distJ_EM, k = 4)
nmdsJ_EM
# Stress = 0.100

# Calculate Bray-Curtis distances and NMDS points
distB_EM = vegdist(relEM, method = 'bray')
set.seed(1986)
nmdsB_EM = metaMDS(distB_EM, k = 3, trymax = 40)
nmdsB_EM = metaMDS(distB_EM, k = 3, trymax = 40,
                   previous.best = nmdsB_EM)
nmdsB_EM
# Stress =  0.112

# Differences between ecotypes using Jaccard distances
set.seed(1986)
adJ_EM = adonis2(distJ_EM ~ ecotype, data = data)
adJ_EM
bdJ_EM = betadisper(distJ_EM, data$ecotype)
set.seed(1986)
bdJp_EM = permutest(bdJ_EM, permutations = 999)
bdJp_EM

# Differences between ecotypes using Bray-Curtis distances
set.seed(1986)
adB_EM = adonis2(distB_EM ~ ecotype, data = data)
adB_EM
bdB_EM = betadisper(distB_EM, data$ecotype)
set.seed(1986)
bdBp_EM = permutest(bdB_EM, permutations = 999)
bdBp_EM

#### (3) Plot ordination #######################################################

##### (3a) AM Fungi #####
points_AM = nmdsB_AM$points # extract nmds points
plot_AM = cbind(dataAM, points_AM) # bind to metadata
glimpse(plot_AM)

# Calculate centroids
centroid_AM = plot_AM %>%
  group_by(ecotype) %>%
  summarise(NMDS1 = mean(MDS1),
            NMDS2 = mean(MDS2), .groups = 'drop')

# Plot ordinatio
ordAMplot = ggplot(plot_AM, aes(MDS1, MDS2, colour = ecotype, fill = ecotype)) +
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
    fill = c('#1b9e77', '#d95f02'),
    colour = 'black', 
    shape = 21,
    size = 4,
    stroke = 0.5,
    show.legend = F
  ) +
  coord_fixed(
    xlim = c(-3.4, 3.3),
    ylim = c(-3.25, 3.45)
  ) +
  labs(x = 'NMDS1',
       y = 'NMDS2') +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    axis.text.x = element_text(colour = 'black', size = rel(1.2)),
    axis.text.y = element_text(colour = 'black', size = rel(1)),
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
  scale_y_continuous(breaks = c(-3, 0, 3)) +
  scale_x_continuous(breaks = c(-3, 0, 3)) +
  # ggtitle("Arbuscular mycorrhizal") +
  theme(plot.title = element_text(hjust = 0.5, size = rel(1))) +
  geom_text(label = 'e', x = -3.4, y = 3.45, colour = 'black') +
  annotate('text', x = -3.4, y = -3.25,
           label = "paste(\ 'Stress = 0.067')",
           parse = T, size = 3.2, hjust = 0, vjust = 0)  +
  annotate('text', x = 3.2, y = -3.25,
           label = "paste(\ '' , \ italic(F) ['39'], \ ' = 3.76*** \')",
           parse = T, size = 3.2, hjust = 0.9, vjust = 0.25)
ordAMplot

##### (3b) EM Fungi #####

points_EM = nmdsB_EM$points # extract nmds points
plot_EM = cbind(data, points_EM) # bind to metadata
glimpse(plot_EM)

centroid_EM = plot_EM %>%
  group_by(ecotype) %>%
  summarise(NMDS1 = mean(MDS1),
            NMDS2 = mean(MDS2), .groups = 'drop')

ordEMplot = ggplot(plot_EM, aes(MDS1, MDS2, colour = ecotype, fill = ecotype)) +
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
    centroid_EM, mapping = aes(NMDS1, NMDS2),
    fill = c('#1b9e77', '#d95f02'),
    colour = 'black', 
    shape = 21,
    size = 4,
    stroke = 0.5,
    show.legend = F
  ) +
  coord_fixed(
    xlim = c(-2.5, 2.5),
    ylim = c(-2.8, 2.2)
  ) +
  labs(x = 'NMDS1',
       y = NULL) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    legend.position = 'none',
    axis.text.x = element_text(colour = 'black', size = rel(1.2)),
    axis.text.y = element_text(colour = 'black', size = rel(1)),
    axis.ticks = element_blank()
  ) +
  scale_y_continuous(breaks = c(-2, 0, 2)) +
  scale_x_continuous(breaks = c(-2, 0, 2)) +
  # ggtitle("Ectomycorrhizal") +
  theme(plot.title = element_text(hjust = 0.5, size = rel(1))) +
  geom_text(label = 'f', x = -2.5, y = 2.15, colour = 'black') + 
  annotate('text', x = -2.5, y = -2.8,
           label = "paste(\ 'Stress = 0.112')",
           parse = T, size = 3.2, hjust = 0, vjust = 0) + 
  annotate('text', x = 2.55, y = -2.8,
           label = "paste(\ '' , \ italic(F) ['39'], \ ' = 2.35** \')",
           parse = T, size = 3.2, hjust = 1, vjust = 0.25) 
ordEMplot

##### (3c) Join plots #####

ordAMEM = ordAMplot + ordEMplot
ordAMEM

ggsave('output/ordination.pdf', width = 6, height = 3.25)
ggsave('output/ordination.jpg', width = 6, height = 3.25)
ggsave('output/ordination.tiff', width = 6, height = 3.25)
