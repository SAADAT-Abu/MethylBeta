#' Filter Differentially Methylated Positions (DMPs)
#'
#' This function filters the results dataframe to identify significant DMPs based on
#' specified thresholds for log fold change (logFC) and adjusted p-value.
#'
#' @param results_df A dataframe containing regression results, including logFC and adjusted p-values.
#' @param filter_by Option to choose "logFC" or "MeanDiff" to filter DMPs. Default is "logFC"
#' @param logFC_threshold A numeric value specifying the minimum absolute log fold change for significance. Used when filter_by is set to "logFC". Default is 1.
#' @param meanDiff_threshold A numeric value specifying the minimum absolute mean difference for significance. Used when filter_by is set to "MeanDiff". Default is 0.1.
#' @param pval_adj_threshold A numeric value specifying the maximum adjusted p-value for significance. Default is 0.05.
#' @return A dataframe containing only the DMPs that meet the specified criteria.
#' @examples
#' \dontrun{
#' filtered_dmps <- filter.DMPs(results_df, filter_by = "logFC", logFC_threshold = 1, pval_adj_threshold = 0.05)
#' filtered_dmps <- filter.DMPs(results_df, filter_by = "MeanDiff", meanDiff_threshold = 0.1, pval_adj_threshold = 0.05)
#' }
#' @export

filter.DMPs <- function(results_df, filter_by = "logFC", logFC_threshold = 1, meanDiff_threshold = 0.1, pval_adj_threshold = 0.05) {
  # Ensure the input dataframe has the necessary columns
  required_columns <- c("logFC", "MeanDiff", "Phenotype_Pval_Adj")
  if (!all(required_columns %in% colnames(results_df))) {
    stop("The input dataframe must contain 'logFC', 'MeanDiff', and 'Phenotype_Pval_Adj' columns.")
  }

  # Filter the dataframe based on the chosen method
  if (filter_by == "logFC") {
    filtered_dmps <- results_df[
      abs(results_df$logFC) >= logFC_threshold & results_df$Phenotype_Pval_Adj <= pval_adj_threshold,
    ]
  } else if (filter_by == "MeanDiff") {
    filtered_dmps <- results_df[
      abs(results_df$MeanDiff) >= meanDiff_threshold & results_df$Phenotype_Pval_Adj <= pval_adj_threshold,
    ]
  } else {
    stop("Invalid filter_by option. Choose 'logFC' or 'MeanDiff'.")
  }

  # Return the filtered dataframe
  return(filtered_dmps)
}

# Example usage:
# filtered_dmps_logFC <- filter.DMPs(results_df_2, filter_by = "logFC", logFC_threshold = 1, pval_adj_threshold = 0.05)
# filtered_dmps_meanDiff <- filter.DMPs(results_df_2, filter_by = "MeanDiff", meanDiff_threshold = 0.1, pval_adj_threshold = 0.05)
