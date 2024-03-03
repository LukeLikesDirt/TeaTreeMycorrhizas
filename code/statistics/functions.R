
### Summarise model coefficients ###############################################
summarise_results <- function(model, test_type, exponentiate = FALSE, standardize = NULL) {
  # Retrieve parameter summary data
  summary_data <- parameters(model, exponentiate = exponentiate, standardize = standardize) %>%
    as_tibble() %>%
    mutate(
      Parameter = case_when(
        Parameter == "(Intercept)" ~ "Upland ecotype (Intercept)",
        Parameter == "ecotypeCoastal" ~ "Coastal ecotype",
        TRUE ~ Parameter
      ),
      across(c(Coefficient, SE, CI, CI_low, CI_high, all_of(test_type)), ~round(., 2))
    )
  
  # Return formatted summary data
  summary_data
}


### Jitter values for ggplot ##################################################

jitter_values <- function(factor_variable, model) {
  # Create a tibble with group and fitted values
  result <- tibble(
    group = as.factor(factor_variable),
    value = as.numeric(fitted(model))
  ) 
  
  # Mutate columns to add jittered values
  result <- result %>%
    mutate(lower_ci = value,
           lower_se = value,
           mean = value,
           upper_se = value,
           upper_ci = value)
  
  return(result)
}

### Effect sizes for ggplot ###################################################

effect_size_coefficients <- function(model) {
  parameters(model, effects = "fixed") %>%
    as_tibble() %>%
    select(Coefficient, CI_low, CI_high) %>%
    # Here I slice the first two rows of the data frame to remove MEMs
    slice(1:2) %>%
    mutate(
      ecotype = c("Upland", "Coastal"),
      ecotype = factor(ecotype, levels = c("Upland", "Coastal"))
    )
}

#### GGPLOT THEME ##############################################################

MyTheme <- function() {

  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    axis.text = element_text(colour = 'black'),
    axis.title = element_text(size = rel(1)),
    axis.ticks = element_blank()
  )

}

#### Alpha diversity functions #################################################

calculate_alpha_diversity <- function(otu_table) {
  
  # Richness function
  calculate_richness <- function(x) {
    sum(x > 0)
  }
  
  # Shannon diversity function
  calculate_shannon_diversity <- function(x) {
    rabund = x[x > 0] / sum(x)
    -sum(rabund * log(rabund))
  }
  
  alpha_diversity <- t(otu_table) %>%
    as.data.frame() %>%
    rownames_to_column(var = 'sample') %>%
    as_tibble() %>% 
    pivot_longer(-sample) %>%
    group_by(sample) %>%
    summarise(
      richness = calculate_richness(value),
      shannon_diversity = calculate_shannon_diversity(value)
    ) %>%
    arrange(richness) %>%
    print(n = Inf)
  
  return(alpha_diversity)
}

### Relative abundance calculation and remove samples within a single sample ###

prevelance_filter_relative_abundance <- function(
    otu_table, prevelance_threshold
    ) {
  
  # Relative abundance calculation function
  calculate_relative_abundance <- function(x){
    
    x / sum(x)
    
  }
  
  t(otu_table) %>%
    as.data.frame() %>%
    rownames_to_column(var = 'sample') %>%
    as_tibble() %>%
    pivot_longer(-sample, names_to = 'OTU_ID', values_to = 'value') %>%
    # Remove zeros
    filter(value > 0) %>%
    # Count each occurrence of an OTU
    mutate(count = 1) %>%
    # Calculate prevalence
    group_by(OTU_ID) %>%
    mutate(prevalence = sum(count)) %>%
    ungroup() %>%
    # Remove taxa that occur once
    filter(prevalence >= prevelance_threshold) %>%
    select(-prevalence) %>%
    # Calculate relative abundances within samples to account for differences
    # in sequencing depth
    group_by(sample) %>%
    mutate(
      rel_abund = calculate_relative_abundance(count)
    ) %>%
    # Remove counts
    select(-c(count, value)) %>%
    # Format OTU table an OTU table of relative abundances of OTUs that occur in
    # at least two samples
    pivot_wider(names_from = OTU_ID, values_from = rel_abund) %>%
    replace(is.na(.), 0) %>%
    as.data.frame() %>%
    column_to_rownames(var = 'sample')
}

# Relative abundance
calculate_rel_abundance <- function(otu_table, taxa_table, taxon_rank) {
  library(dplyr)
  library(tidyr)
  library(rlang)  # Load rlang for unquoting
  
  if (is.character(otu_table) && is.character(taxa_table)) {
    # Paths are provided, read data from paths
    otu_tab <- read.csv(otu_table, header = TRUE, row.names = 'OTU_ID') %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column(var = 'sample') %>%
      pivot_longer(-sample, names_to = 'OTU_ID', values_to = 'value')
    
    taxa_tab <- read.csv(taxa_table, header = TRUE)
  } else if (is.data.frame(otu_table) && is.data.frame(taxa_table)) {
    # Data frames are provided directly
    otu_tab <- otu_table %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column(var = 'sample') %>%
      pivot_longer(-sample, names_to = 'OTU_ID', values_to = 'value')
    
    taxa_tab <- taxa_table
  } else {
    stop("Invalid input types. Provide either paths or data frames.")
  }
  
  relab <- inner_join(otu_tab, taxa_tab, by = 'OTU_ID') %>%
    group_by(sample) %>%
    mutate(rel_abund = value / sum(value)) %>%
    ungroup() %>%
    pivot_longer(cols = c('phylum', 'class', 'order', 'family', 'genus', 'species'),
                 names_to = 'taxon_rank',
                 values_to = 'taxon') %>%
    mutate(is_selected_rank = ifelse(taxon_rank == "genus" & !is.nan(rel_abund), TRUE, FALSE)) %>%
    filter(is_selected_rank) %>%
    group_by(sample, taxon) %>%
    summarise(rel_abund = sum(rel_abund), .groups = "drop") %>%
    group_by(taxon) %>%
    summarise(rel_abund = 100 * mean(rel_abund), .groups = "drop") %>%
    mutate(across(where(is.numeric), round, 1)) %>%
    print(n = Inf)
}

### Relative abundance #########################################################

relabund = function(x){
  
  x / sum(x)
  
}
