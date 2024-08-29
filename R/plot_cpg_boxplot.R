#' Plot Methylation Levels for a Specific CpG Site
#'
#' This function creates a boxplot of methylation values for a specific CpG site,
#' comparing the values across two phenotype groups.
#'
#' @param beta_matrix A matrix of beta values with CpG sites as rows and samples as columns.
#' @param sample_groups A vector indicating the phenotype group for each sample.
#' @param cpg_id The name of the CpG site to plot.
#' @return A ggplot object representing the boxplot of methylation levels.
#' @examples
#' \dontrun{
#' plot_cpg_boxplot(myNorm, myPD$Sample_Group, "cg18478105")
#' }
#' @export
plot_cpg_boxplot <- function(beta_matrix, sample_groups, cpg_id) {
  # Check if the CpG ID exists in the matrix
  if (!(cpg_id %in% rownames(beta_matrix))) {
    stop("CpG ID not found in the beta matrix.")
  }

  # Extract methylation values for the specified CpG
  beta_values <- beta_matrix[cpg_id, ]

  # Create a dataframe for plotting
  plot_data <- data.frame(
    Methylation = beta_values,
    Phenotype = sample_groups
  )

  # Generate the boxplot
  ggplot(plot_data, aes(x = Phenotype, y = Methylation, fill = Phenotype)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(title = paste("Methylation Levels at", cpg_id),
         x = "Sample Group",
         y = "Methylation Value") +
    theme_minimal() +
    theme(legend.position = "none")
}
