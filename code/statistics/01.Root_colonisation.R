
# Script: Assess differences in mycorrhizal root colonisation across tea tree
# ecotypes

# Required packages
require(lme4)
require(emmeans)
require(performance)
require(report)
require(DHARMa)
require(patchwork)
require(ggplot2)
require(tidyverse)

# Metadata
data = read.csv('data/statistics/metadata.csv', stringsAsFactors = TRUE,
                header = TRUE)
glimpse(data)

#### (1) AM root colonisation ##################################################

# Fit models
amB.1 = glmer(colAM/weight ~ ecotype + (1 | site),
              weights = weight, family = 'binomial', data = data)
amB.2 = glmer(colAM/weight ~ ecotype + (1 | site/sample),
              weights = weight, family = 'binomial', data = data)

# Compare model performance
compare_performance(amB.1, amB.2) # Model amB.2 is the best fit

# Model validation
set.seed(1986)
resAM = simulateResiduals(amB.2)
plot(resAM) # Normality, homogeneity and dispersion look good
check_predictions(amB.2) # Predictions resemble observed data

# Estimated marginal means
am.emm = emmeans(amB.2, ~ ecotype, type = 'response') %>%
  as.data.frame() %>%
  mutate(mean = prob * 100,
         upper_se = mean + (SE * 100),
         lower_se = mean - (SE * 100),
         upper_ci = asymp.UCL * 100,
         lower_ci = asymp.LCL * 100) %>%
  select(group = ecotype, mean, upper_se, lower_se, upper_ci, lower_ci)
am.emm
summary(amB.2)
# z = 1.19, p = 0.235

# # Fitted values
AMfit = cbind(as.character(data$ecotype), fitted(amB.2) * 100) %>%
  as.data.frame() %>%
  rename(group = V1) %>%
  rename(value = V2)
AMfit$group = as.factor(AMfit$group)
AMfit$value = as.numeric(AMfit$value)
glimpse(AMfit)

am.jitter = AMfit %>%
  mutate(lower_ci = value) %>%
  mutate(lower_se = value) %>%
  mutate(mean = value) %>%
  mutate(upper_se = value) %>%
  mutate(upper_ci = value)
# Plot
colAMplot = ggplot(am.emm, aes(x = group, 
                               ymin = lower_ci,
                               lower = lower_se,
                               middle = mean,
                               upper = upper_se,
                               ymax = upper_ci,
                               fill = group)) +
  geom_boxplot(stat = 'identity') +
  # Here I use point instead of jitter so I can a set seed so points do not overlap error bars, which is distracting
  geom_point(
    am.jitter, 
    mapping = aes(x = group, y = value, fill = group), 
    position = position_jitter(seed = 16, width = 0.25),
    shape = 21,
    size = 2,
    alpha = 0.66
  ) +
  scale_fill_manual(
    values = c('#1b9e77', '#d95f02'),
    breaks = c('Coastal', 'Upland')
  )  +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    legend.position = 'none',
    axis.text.x = element_text(colour = 'black', size = rel(1.2)),
    axis.text.y = element_text(colour = 'black', size = rel(1)),
    axis.ticks = element_blank()
  ) +
  scale_y_continuous(
    limits = c(0, 73.75),
    breaks = c(0, 20, 40, 60)
  ) +
  xlab(element_blank()) + ylab('Root colonisation (%)') +
  ggtitle('Arbuscular mycorrhizal') +
  theme(plot.title = element_text(hjust = 0.5, size = rel(1))) +
  geom_text(label = 'a', x = 0.5, y = 73, colour = 'black') +
  annotate('text', x = 1.55, y = 73.5,
           label = 'paste(\ italic(z),\'-value = 1.19 \')',
           parse = T, size = 3.2, hjust = -0, vjust = 1)
colAMplot

#### (2) EM root colonisation ##################################################

# Fit models
emB.1 = glmer(colEM/weight ~ ecotype + (1 | site),
              weights = weight, family = 'binomial', data = data)
emB.2 = glmer(colEM/weight ~ ecotype + (1 | site/sample),
              weights = weight, family = 'binomial', data = data)

# Compare performance
compare_performance(emB.1, emB.2) # emB.2 looks like the best fit

# Validate model
set.seed(1986)
resEM = simulateResiduals(emB.2)
plot(resEM) # Normality, homogeneity and dispersion look good
check_predictions(emB.2) # Predictions resemble observed data

# Estimated marginal means
em.emm <- emmeans(emB.2, ~ ecotype, type = 'response') %>%
  as.data.frame() %>%
  mutate(mean = prob * 100,
         upper_se = mean + (SE * 100),
         lower_se = mean - (SE * 100),
         upper_ci = asymp.UCL * 100,
         lower_ci = asymp.LCL * 100) %>%
  select(group = ecotype, mean, upper_se, lower_se, upper_ci, lower_ci)
em.emm
summary(emB.2)
# z = 3.71, p = 0.000207

# Fitted values
EMfit = cbind(as.character(data$ecotype), fitted(emB.2) * 100) %>%
  as.data.frame() %>%
  rename(group = V1) %>%
  rename(value = V2)
EMfit$group = as.factor(EMfit$group)
EMfit$value = as.numeric(EMfit$value)
glimpse(EMfit)

# Work around to overlay fitted values as jitter point
em.jitter = EMfit %>%
  mutate(lower_ci = value) %>%
  mutate(lower_se = value) %>%
  mutate(mean = value) %>%
  mutate(upper_se = value) %>%
  mutate(upper_ci = value)

# Plot
colEMplot = ggplot(em.emm, aes(x = group, 
                               ymin = lower_ci,
                               lower = lower_se,
                               middle = mean,
                               upper = upper_se,
                               ymax = upper_ci,
                               fill = group)) +
  geom_boxplot(stat = 'identity') +
  scale_fill_manual(
    values = c('#1b9e77', '#d95f02'),
    breaks = c('Coastal', 'Upland')
  )  +
  # Here I use point instead of jitter so I can a set seed so points do not overlap error bars, which is distracting
  geom_point(
    em.jitter, 
    mapping = aes(x = group, y = value, fill = group), 
    position = position_jitter(seed = 26, width = 0.25),
    shape = 21,
    size = 2,
    alpha = 0.66
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    legend.position = 'none',
    axis.text.x = element_text(colour = 'black', size = rel(1.2)),
    axis.text.y = element_text(colour = 'black', size = rel(1)),
    axis.ticks = element_blank()
  ) +
  scale_y_continuous(
    limits = c(0, 73.75),
    breaks = c(0, 20, 40, 60)
  ) +
  xlab(element_blank()) + ylab(element_blank()) +
  ggtitle('Ectomycorrhizal') +
  theme(plot.title = element_text(hjust = 0.5, size = rel(1))) +
  geom_text(label = 'b', x = 0.5, y = 73, colour = 'black') +
  annotate('text', x = 1.4, y = 73.75,
           label = 'paste(\ italic(z),\'-value = 3.71*** \')',
           parse = T, size = 3.2, hjust = 0, vjust = 1)
colEMplot

#### (3) Join plot #############################################################

colAMEM = colAMplot + colEMplot
colAMEM 

ggsave('output/colonisation.pdf', width = 6, height = 3.25)
ggsave('output/colonisation.jpg', width = 6, height = 3.25)
ggsave('output/colonisation.tiff', width = 6, height = 3.25)
