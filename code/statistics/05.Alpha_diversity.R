
# Script: Assess mycorrhizal fungal alpha diversity across ecotypes

# Required packages and functions
source("code/statistics/functions.R")
require(glmmTMB)
require(lme4)
require(emmeans)
require(performance)
require(DHARMa)
require(paletteer)
require(patchwork)
require(tidyverse)

# Read in the data
data <- read.csv('data/statistics/metadata.csv',
                 stringsAsFactors = TRUE, header = TRUE)
otuAM <- read.csv('data/statistics/otuAM.csv',
                  header = TRUE, row.names = 'OTU_ID')
otuEM <- read.csv('data/statistics/otuEM.csv',
                  header = TRUE, row.names = 'OTU_ID')

#### (1) Assess AM diversity ###################################################

# Estimate alpha diversity
alpha_diversity_AM = t(otuAM) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample') %>%
  as_tibble() %>% 
  pivot_longer(-sample) %>%
  group_by(sample) %>%
  summarise(AM_richness = richness(value),
            AM_shannon = shannon_diversity(value)) %>%
  arrange(AM_richness) %>%
  print(n = 40)
# 17.5 % of samples have zeros richness so check for zero inflation

# Combine environmental data and diversity metrics
alpha_diversity_AM  = inner_join(data, alpha_diversity_AM , by = 'sample')

# Fit GLMERs for richness
amP.1 = glmer(AM_richness ~ ecotype + (1|site),
              family = poisson, data = alpha_diversity_AM)
amP.2 = glmer(AM_richness ~ ecotype + (1|site/sample),
              family = poisson, data = alpha_diversity_AM)

# Compare model performance
compare_performance(amP.1, amP.2)
# amP.2 has a better fit based on AIC and predictions based on RMSE

# Check assumptions
set.seed(1986)
res = simulateResiduals(amP.2)
plot(res)  # Normality, homogeneity and dispersion look good
testZeroInflation(res) # No significant issues with zero inflation
check_predictions(amP.2) # A little wide/over predicted at the peak but not bad

# Fit LMERs for Shannon diversity
amG.1 = glmmTMB(AM_shannon ~ ecotype + (1|site),
                family = gaussian(), data = alpha_diversity_AM)
amG.2 = glmmTMB(AM_shannon ~ ecotype + (1|site/sample),
                family = gaussian(), data = alpha_diversity_AM)

# Compare model performance
compare_performance(amG.1, amG.2) 
# Similar based on AIC but slightly better predilections based on 
# RMSE. I'll stick with the nested model for consistency.

# Check assumptions
set.seed(1986)
res = simulateResiduals(amG.2)
plot(res) # Normality, homogeneity and dispersion looks good
check_predictions(amG.2) # A bit squirrely

#### (2) Assess EM diversity ###################################################

# Estimate alpha diversity
alpha_diversity_EM = t(otuEM) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample') %>%
  as_tibble() %>% 
  pivot_longer(-sample) %>%
  group_by(sample) %>%
  summarise(EM_richness = richness(value),
            EM_shannon = shannon_diversity(value)) %>%
  arrange(EM_richness) %>%
  print(n = 40)

# Combine environmental data and diversity metrics
alpha_diversity_EM  = inner_join(data, alpha_diversity_EM , by = 'sample')

# Fit GLMERs for richness
emP.1 = glmer(EM_richness ~ ecotype + (1|site),
              family = poisson, data = alpha_diversity_EM)
emP.2 = glmer(EM_richness ~ ecotype + (1|site/sample),
              family = poisson, data = alpha_diversity_EM)

# Compare model performance
compare_performance(emP.1, emP.2)
# amP.2 has a better fit based on AIC and predictions based on RMSE

# Check assumptions
set.seed(1986)
res = simulateResiduals(emP.2)
plot(res) # Normality, homogeneity and dispersion look good
check_predictions(emP.2) # Pretty good predicitions

# Fit LMERs for Shannon diversity
emG.1 = glmmTMB(EM_shannon ~ ecotype + (1|site),
                family = gaussian(), data = alpha_diversity_EM)
emG.2 = glmmTMB(EM_shannon ~ ecotype + (1|site/sample),
                family = gaussian(), data = alpha_diversity_EM)

# Compare model performance
compare_performance(emG.1, emG.2)
# More or less the same. I'll stick with the nested model for constancy.

# Check assumptions
set.seed(1986)
res = simulateResiduals(emG.2)
plot(res) # Normality, homogeneity and dispersion looks good
check_predictions(emG.2) # Not too bad.

#### (3) Print results #########################################################

# AM richness: Estimated marginal means
amP.emm = emmeans(amP.2, ~ ecotype, type = "response") %>%
  as.data.frame() %>%
  mutate(mean = rate,
         upper_se = mean + SE,
         lower_se = mean - SE,
         upper_ci = asymp.UCL,
         lower_ci = asymp.LCL) %>%
  select(group = ecotype, mean, upper_se, lower_se, upper_ci, lower_ci)
amP.emm
# AM Shannon: Estimated marginal means
amG.emm = emmeans(amG.2, ~ ecotype, type = "response") %>%
  as.data.frame() %>%
  mutate(mean = emmean,
         upper_se = mean + SE,
         lower_se = mean - SE,
         upper_ci = upper.CL,
         lower_ci = lower.CL) %>%
  select(group = ecotype, mean, upper_se, lower_se, upper_ci, lower_ci)
amG.emm

summary(amP.2)
report::report(amP.2)
report::report_table(amP.2)
# z = 2.32; p = 0.020
summary(amG.2)
report::report(amG.2)
# z = 2.60; p = 0.009

# EM richness: Estimated marginal means
emP.emm = emmeans(emP.2, ~ ecotype, type = "response") %>%
  as.data.frame() %>%
  mutate(mean = rate,
         upper_se = mean + SE,
         lower_se = mean - SE,
         upper_ci = asymp.UCL,
         lower_ci = asymp.LCL) %>%
  select(group = ecotype, mean, upper_se, lower_se, upper_ci, lower_ci)
emP.emm
# EM Shannon: Estimated marginal means
emG.emm = emmeans(emG.2, ~ ecotype, type = "response") %>%
  as.data.frame() %>%
  mutate(mean = emmean,
         upper_se = mean + SE,
         lower_se = mean - SE,
         upper_ci = upper.CL,
         lower_ci = lower.CL) %>%
  select(group = ecotype, mean, upper_se, lower_se, upper_ci, lower_ci)
emG.emm

summary(emP.2)
report::report(emP.2)
report::report_table(emP.2)
# z = 3.48; p = 0.000505
summary(emG.2)
report::report(emG.2)
# z = 3.27; p = 0.00108

#### (4) Plot AM richness ######################################################

# Fitted values
AMfit = cbind(as.character(data$ecotype), fitted(amP.2)) %>%
  as.data.frame() %>%
  rename(group = V1) %>%
  rename(value = V2)
AMfit$group = as.factor(AMfit$group)
AMfit$value = as.numeric(AMfit$value)
glimpse(AMfit)

# Work around to overlay fitted oo GGPLOT
am.jitter = AMfit %>%
  mutate(lower_ci = value) %>%
  mutate(lower_se = value) %>%
  mutate(mean = value) %>%
  mutate(upper_se = value) %>%
  mutate(upper_ci = value)

# Plot
richAMplot <-
  ggplot(amP.emm, 
         aes(x = group, ymin = lower_ci, lower = lower_se, middle = mean,
             upper = upper_se, ymax = upper_ci, fill = group)) +
  geom_boxplot(stat = 'identity') +
  scale_fill_manual(
    values = c('#1b9e77', '#d95f02'),
    breaks = c('Coastal', 'Upland')
  )  +
  geom_point(
    am.jitter, 
    mapping = aes(x = group, y = value, fill = group), 
    position = position_jitter(seed = 16, width = 0.25),
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
    limits = c(0, 23),
    breaks = c(0, 10, 20)
  ) +
  xlab(element_blank()) + ylab('OTU richness') +
  # ggtitle('Arbuscular mycorrhizal') +
  theme(plot.title = element_text(hjust = 0.5, size = rel(1)),
        legend.position = 'none') +
  geom_text(label = 'c', x = 0.5, y = 21.75, colour = 'black') +
  annotate('text', x = 1.5, y = 23,
           label = 'paste(\ italic(z),\'-value = 2.32* \')',
           parse = T, size = 3.2, hjust = -0, vjust = 1)
richAMplot

#### (5) Plot ECM richness #####################################################

# Fitted values
EMfit = cbind(as.character((data$ecotype)), fitted(emP.2)) %>%
  as.data.frame() %>%
  rename(group = V1) %>%
  rename(value = V2)
EMfit$group = as.factor(EMfit$group)
EMfit$value = as.numeric(EMfit$value)
glimpse(EMfit)
# Work around to overlay fitted values as jitter point
EMjitter = EMfit %>%
  mutate(lower_ci = value) %>%
  mutate(lower_se = value) %>%
  mutate(mean = value) %>%
  mutate(upper_se = value) %>%
  mutate(upper_ci = value)
# Plot
richEMplot <-
  ggplot(emP.emm, 
         aes(x = group, ymin = lower_ci, lower = lower_se, middle = mean,
             upper = upper_se, ymax = upper_ci, fill = group)) +
  geom_boxplot(stat = 'identity') +
  scale_fill_manual(
    values = c('#1b9e77', '#d95f02'),
    breaks = c('Coastal', 'Upland')
  ) +
  geom_point(
    EMjitter, 
    mapping = aes(x = group, y = value, fill = group), 
    position = position_jitter(seed = 16, width = 0.25),
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
    limits = c(0, 23),
    breaks = c(0, 10, 20)
    ) +
  xlab(element_blank()) + ylab(NULL) +
  # ggtitle('Ectomycorrhizal') +
  theme(plot.title = element_text(hjust = 0.5, size = rel(1)),
        legend.position = 'none') +
  geom_text(label = 'd', x = 0.5, y = 21.75, colour = 'black') +
  annotate('text', x = 1.4, y = 23,
           label = 'paste(\ italic(z),\'-value = 3.48*** \')',
           parse = T, size = 3.2, hjust = -0, vjust = 1)
richEMplot

#### (6) Join plots ############################################################

richAMEM = richAMplot + richEMplot
richAMEM

ggsave('output/richness.pdf', width = 6, height = 3.25)
ggsave('output/richness.jpg', width = 6, height = 3.25)
ggsave('output/richness.tiff', width = 6, height = 3.25)
