#' Plot Heatmap of Significant Differentially Methylated Positions (DMPs)
#'
#' This function creates a heatmap of methylation values for significant DMPs,
#' allowing visualization of methylation patterns across samples.
#'
#' @param beta_matrix A matrix of beta values with CpG sites as rows and samples as columns.
#' @param dmps A dataframe of significant DMPs, typically filtered by logFC and adjusted p-value.
#' @return A heatmap displaying methylation levels of significant DMPs.
#' @examples
#' \dontrun{
#' filtered_dmps <- filter.DMPs(results_df, logFC_threshold = 1, pval_adj_threshold = 0.05)
#' plot_heatmap(myNorm, filtered_dmps)
#' }
#' @export
plot_heatmap <- function(beta_matrix, dmps) {
  # Subset the beta matrix for significant DMPs
  significant_beta_matrix <- beta_matrix[rownames(beta_matrix) %in% rownames(dmps), ]
  
  # Plot the heatmap
  pheatmap::pheatmap(significant_beta_matrix, cluster_rows = TRUE, cluster_cols = TRUE,
           show_rownames = TRUE, show_colnames = TRUE,
           main = "Heatmap of Significant DMPs")
}
