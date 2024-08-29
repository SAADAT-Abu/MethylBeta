#' Create a Volcano Plot of Differentially Methylated Positions (DMPs)
#'
#' This function generates a volcano plot to visualize the results of DMP analysis,
#' highlighting significant changes in methylation based on log fold change and adjusted p-value.
#'
#' @param results_df A dataframe containing regression results, including logFC and adjusted p-values.
#' @param logFC_threshold A numeric value specifying the minimum absolute log fold change for significance. Default is 1.
#' @param pval_adj_threshold A numeric value specifying the maximum adjusted p-value for significance. Default is 0.05.
#' @return A ggplot object representing the volcano plot.
#' @examples
#' \dontrun{
#' plot_volcano(results_df, logFC_threshold = 1, pval_adj_threshold = 0.05)
#' }
#' @export
plot_volcano <- function(results_df, logFC_threshold = 1, pval_adj_threshold = 0.05) {
  ggplot(results_df, aes(x = logFC, y = -log10(Phenotype_Pval_Adj))) +
    geom_point(aes(color = (abs(logFC) > logFC_threshold & Phenotype_Pval_Adj < pval_adj_threshold)), alpha = 0.6) +
    scale_color_manual(values = c("grey", "red")) +
    labs(title = "Volcano Plot of DMPs", x = "Log Fold Change", y = "-Log10 Adjusted p-value") +
    theme_minimal() +
    theme(legend.position = "none")
}
