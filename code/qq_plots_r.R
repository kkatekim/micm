library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)

# runs fisher's test on case/control values
create_fisher <- function(vector_variants_cases_controls, n_cases, n_controls) {
  results <- fisher.test(cbind(vector_variants_cases_controls, (c(n_cases, n_controls) - vector_variants_cases_controls)))
  return(c(pval=results$p.value, lower_CI=results$conf.int[1], upper_CI=results$conf.int[2], OR=results$estimate))
}

# matrix of potential pvalues calculated from n_cases and n_controls
create_fisher_lookup_table <- function(
  n_cases, n_controls, df, case_col, control_col)
{	
  maximum_possible_carriers <- max(rowSums(df %>% select(case_col, control_col)))
  # Predefine contingency table
  contingency_table <- matrix(0, nrow=2, ncol=2)
  Fisher_precompute <- matrix(0,
                              nrow=(min(n_cases, maximum_possible_carriers)+1),
                              ncol=(min(n_controls, maximum_possible_carriers)+1))
  for (i in 1:nrow(Fisher_precompute))
  {
    contingency_table[1,1] <- i - 1
    contingency_table[1,2] <- n_cases - contingency_table[1,1]
    for (j in 1:ncol(Fisher_precompute))
    {
      contingency_table[2,1] <- j - 1
      contingency_table[2,2] <- n_controls - contingency_table[2,1]
      Fisher_precompute[i,j] <- fisher.test(contingency_table)$p.value
    }
    cat(paste0(i,' of ', dim(Fisher_precompute)[1],' possibilities \n'))
  }
  return(Fisher_precompute)
}

# calculate expected p value based on permutations 
create_permutation_p_fisher <- function(
  n_permutations, df, n_cases, n_controls,
  Fisher_lookup, case_col, control_col)
{
  n_genes <- nrow(df)
  n <- n_cases + n_controls
  p_values <- rep(0, n_genes)
  # For each gene, find out how many mutations need to be allocated.
  mutations_for_allocation <- rowSums(df %>% select(case_col, control_col))
  apply_rhyper <- function(mutations_for_allocation, n_cases, n_controls, n_perm) {
    rhyper(n_perm, n_cases, n_controls, mutations_for_allocation)
  }
  for(i in 1:n_permutations) {
    n_cases_with_mut <- sapply(mutations_for_allocation, apply_rhyper, n_cases, n_controls, 1)
    n_controls_with_mut <- mutations_for_allocation - n_cases_with_mut
    p_values <- p_values + sort(Fisher_lookup[cbind(n_cases_with_mut+1, n_controls_with_mut+1)])
    if ((i %% 10000) == 0) {
      cat('permutation =',i,'\n')
    }
  }
  p_values <- p_values / n_permutations
  return(p_values)
}

# creates dt with exp and obs p values
create_fisher_and_perm <- function(dt, n_cases, n_controls, variant,
                                          case_col, control_col, n_permutations=1000,
                                          save_lookup=TRUE, lookup_filename='lookup.Rdata',
                                          read_lookup=FALSE, pathname='../data/proteinDomain/', 
                                          include_col=FALSE)
{	
  cat(paste0("\nn cases: ", n_cases, "\nn controls: ", n_controls, "\n"))
  print(head(dt %>% select(case_col, control_col)))
  #dt <- dt %>% filter(consequence_category==consequence)
  results <- t(apply(dt %>% select(case_col, control_col), 1, create_fisher, n_cases, n_controls))
  results_dt <- cbind(dt, results)
  # return(results_dt)
  if (read_lookup) {
    cat("Reading Fisher look-up...\n")
    if (file.exists(lookup_filename)) {
      load(lookup_filename)
    } else {
      cat("Specified Rdata file does not exist, creating Fisher look-up...\n")
      Fisher_lookup <- create_fisher_lookup_table(
        n_cases, n_controls, dt,
        case_col, control_col
      )
    }
  } else {
    cat("Creating Fisher look-up...\n")
    Fisher_lookup <- create_fisher_lookup_table(
      n_cases, n_controls, dt,
      case_col, control_col
    )
  }
  if(save_lookup) {
    cat("Saving Fisher look-up...\n")
    save(Fisher_lookup, file=lookup_filename)
  }
  cat("Creating permutation p-values...\n")
  p_perm <- create_permutation_p_fisher(
    n_permutations, dt,
    n_cases, n_controls, Fisher_lookup, case_col, control_col)
  p_obs <- -log10(results_dt %>% arrange(pval) %>% select(pval))$pval
  dt_plot <- data.table(
    p_perm = -log10(p_perm),
    pval = p_obs
    #labels = (results_dt %>%
    #arrange(pval) %>% select(gene_symbol))$gene_symbol
  )
  save(dt_plot, file=paste0(pathname, "dt_plot_", variant, n_permutations, ".Rdata"))
  return(list(dt_plot=dt_plot, results_dt=results_dt))
}

plot_qq <- function(dt_plot, mutation, n_permutations)
{
  p <- ggplot(dt_plot) +
    geom_point(aes(x=p_perm, y=pval, shape=".")) +
    geom_abline(intercept=0, slope=1, colour='blue', lwd=1.2, alpha=0.8) +
    theme_bw() +
    theme(
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=1),
          legend.position = "none") +
    labs(title=paste(mutation, "with", n_permutations, "permutations"), x='exp', y='obs') 
  ggsave(paste0("protein_domain_", mutation, "_", n_permutations, '.png'), plot=p, path='../data/proteinDomain/figures')
}

# n_case = 3864
# n_control = 7839
create_all_qq_plots <- function(filename, pathname, n_permutations)
{
  df <- read.table(filename, sep=',', header=TRUE)
  
  variants <- c("synonymous", "missense", "PTVs", "missense_mpc_2")
  
  for (variant in variants)
  {
    cleaned_df <- df %>% select(ends_with(variant)) %>% drop_na()
    case_col_name <- paste0("case_", variant)
    control_col_name <- paste0("control_", variant)
    lookup_filename <- paste0(pathname, "lookup_", variant, n_permutations, ".Rdata")
    results <- create_fisher_and_perm(df, n_cases=3864, n_controls=7839, variant=variant,
                                             case_col=case_col_name, control_col=control_col_name,
                                             n_permutations=n_permutations, save_lookup=FALSE,
                                             read_lookup=TRUE, lookup_filename=lookup_filename,
                                             pathname=pathname)
    
    dt_plot <- results[[1]]
    #save(dt_plot, file=paste0(pathname, "dt_plot_", variant, n_permutations, ".Rdata"))
    
    plot_qq(dt_plot, mutation = variant, n_permutations = n_permutations)
  }
}
create_all_qq_plots("../data/proteinDomain/all_protein_variants_fishers.csv", "../data/proteinDomain/", 50)