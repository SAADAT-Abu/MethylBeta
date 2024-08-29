#' Perform Beta Regression for DMP Analysis
#'
#' This function fits a beta regression model to each CpG site in the beta value matrix
#' and returns a dataframe with regression results, including log fold change and adjusted p-values.
#'
#' @param beta_matrix A matrix of beta values with CpG sites as rows and samples as columns.
#' @param sample_groups A vector indicating the phenotype group for each sample. First level should represent controls and second level should represent cases
#' @param chunk_size An integer specifying the number of CpG sites to process in each chunk. Default is 100.
#' @param num_cores An integer specifying the number of cores to use for parallel processing. Default is 2.
#' @return A dataframe containing regression results for each CpG site.
#' @examples
#' \dontrun{
#' results_df <- beta.fit(myNorm, myPD$Sample_Group, chunk_size = 100, num_cores = 10)
#' }
#' @export
beta.fit <- function(beta_matrix, sample_groups, chunk_size = 100, num_cores = 2) {

  # Clean CpG IDs to remove numeric prefixes followed by a dot
  rownames(beta_matrix) <- sub("^[0-9]+\\.", "", rownames(beta_matrix))

  # Get the number of CpG sites
  num_cpgs <- nrow(beta_matrix)

  # Function to process a single chunk
  process_chunk <- function(start_idx, end_idx) {
    chunk_results <- list()
    chunk <- beta_matrix[start_idx:end_idx, , drop = FALSE]

    for (cpg in rownames(chunk)) {
      beta <- chunk[cpg, ]
      data <- data.frame(beta = beta, phenotype = sample_groups)

      # Fit the beta regression model
      model <- betareg(beta ~ phenotype, data = data)
      summary_model <- summary(model)

      # Extract necessary statistics
      phenotype_coef <- coef(summary_model)$mean[2, 1]
      phenotype_se <- coef(summary_model)$mean[2, 2]
      phenotype_pval <- coef(summary_model)$mean[2, 4]

      phi <- coef(summary_model)$precision[1, 1]
      phi_se <- coef(summary_model)$precision[1, 2]
      phi_pval <- coef(summary_model)$precision[1, 4]

      # Calculating MeanDiff

      mean_diff <- mean(beta[sample_groups == levels(sample_groups)[2]]) - mean(beta[sample_groups == levels(sample_groups)[1]])
      
      # Store the results in a list
      chunk_results[[cpg]] <- c(
        Phi = phi,
        Phi_SE = phi_se,
        Phi_Pval = phi_pval,
        Phenotype_Coef = phenotype_coef,
        Phenotype_SE = phenotype_se,
        Phenotype_Pval = phenotype_pval,
        MeanDiff = mean_diff
      )
    }

    return(chunk_results)
  }

  # Create indices for chunk processing
  chunk_indices <- split(seq_len(num_cpgs), ceiling(seq_len(num_cpgs) / chunk_size))

  # Process chunks in parallel
  results <- mclapply(chunk_indices, function(indices) {
    start_idx <- indices[1]
    end_idx <- indices[length(indices)]
    process_chunk(start_idx, end_idx)
  }, mc.cores = num_cores)

  # Combine results from all chunks
  results_df <- do.call(rbind, unlist(results, recursive = FALSE)) %>% as.data.frame()

  # Apply multiple testing correction to phenotype p-values
  results_df$Phenotype_Pval_Adj <- p.adjust(results_df$Phenotype_Pval, method = "BH")

  # Calculate logFC
  results_df$logFC <- with(results_df, ifelse(Phenotype_Coef >= 0, log2(Phenotype_Coef + 1), -log2(abs(Phenotype_Coef) + 1)))
  rownames(results_df) <- sub("^[0-9]+\\.", "", rownames(results_df))
  return(results_df[,c(1:6,8,9,7)])
}

