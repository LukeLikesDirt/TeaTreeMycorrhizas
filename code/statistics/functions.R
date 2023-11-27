
#### OTU curation #############################################################

# The OTU table should have OTU IDs as row names and sample IDs as column
# names. Example usage with a specified threshold of 0.01:
# filtered_OTUs <- filter_low_abundance_otus(otu_table, threshold = 0.01)

filter_low_abundance_otus <- function(OTUs, threshold = 0.01) {
  OTUs1 <- OTUs %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = 'sample') %>%
    as_tibble() %>%
    pivot_longer(-sample) %>%
    # Group by sample to account for differences in sequencing depth
    group_by(sample) %>%
    # Calculate relative abundance within samples
    mutate(rel_abund = 100 * value / sum(value)) %>%
    ungroup() %>%
    # Remove OTUs that have relative abundance < threshold within individual
    # samples
    mutate(count = replace(value, rel_abund < threshold, 0)) %>%
    select(-c(value, rel_abund)) %>%
    pivot_wider(values_from = 'count') %>%
    column_to_rownames(var = 'sample') %>%
    t() %>%
    as.data.frame() %>%
    filter(rowSums(.) != 0) %>%
    glimpse()
    
  return(OTUs1)
}

### VARIATION INFLATION FACTOR FUNCTION ########################################

# To cite the VIF function, use:
# Mixed effects models and extensions in ecology with R. (2009).
# Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM. Springer.

# VIF function
VIF <- function(x) {

  x <- as.data.frame(x)
  # VIF calculation
  form    <- formula(paste("fooy ~ ", paste(strsplit(names(x), " "),
                          collapse = " + ")))
  x       <- data.frame(fooy = 1 + rnorm(nrow(x)) ,x)
  lm_mod  <- lm(form, x)
  # End VIF calculation
  cat("\n\nVariance inflation factors\n\n")
  print(supVIF(lm_mod))

}

# VIF function dependency
supVIF <- function(mod) {

  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor) <- 0
    if (any(tmp_cor == 1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)

}

### STANDARDISE AND UNSTANDARDISE CONTINUOUS COVARIATES ###

# Function to standardise predictors
std <- function(x) {

  (x - mean(x)) / sd(x)

}

# Function to back-transforms standardised predictors
unstd <- function(x.std, x) {

  x.std * sd(x) + mean(x)

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

# Richness function
richness <- function(x) {
  
  sum(x > 0)
  
}

# Shannon diversity function
shannon_diversity <- function(x) {
  
  rabund = x[x > 0] / sum(x)
  -sum(rabund * log(rabund))
  
}

# Relative abundance
relabund = function(x){
  
  x / sum(x)
  
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
