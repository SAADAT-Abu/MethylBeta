% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_volcano.R
\name{plot_volcano}
\alias{plot_volcano}
\title{Create a Volcano Plot of Differentially Methylated Positions (DMPs)}
\usage{
plot_volcano(results_df, logFC_threshold = 1, pval_adj_threshold = 0.05)
}
\arguments{
\item{results_df}{A dataframe containing regression results, including logFC and adjusted p-values.}

\item{logFC_threshold}{A numeric value specifying the minimum absolute log fold change for significance. Default is 1.}

\item{pval_adj_threshold}{A numeric value specifying the maximum adjusted p-value for significance. Default is 0.05.}
}
\value{
A ggplot object representing the volcano plot.
}
\description{
This function generates a volcano plot to visualize the results of DMP analysis,
highlighting significant changes in methylation based on log fold change and adjusted p-value.
}
\examples{
\dontrun{
plot_volcano(results_df, logFC_threshold = 1, pval_adj_threshold = 0.05)
}
}
