# MethylBeta | DMP Analysis Package

![MethylBeta_readme](https://github.com/user-attachments/assets/da0b0748-c83e-46b5-9992-9ccbfa19751a)

## Overview

This R package provides tools for analyzing differentially methylated positions (DMPs) using beta regression. It includes functions for regression analysis, filtering significant DMPs, and visualizing results.

## Features

- **Beta Regression Analysis:** Perform beta regression on methylation data to identify DMPs.
- **Flexible Filtering:** Filter DMPs based on log fold change or mean difference.
- **Visualization Tools:** Create volcano plots, heatmaps and boxplots to visualize methylation data.

## Installation

To install the package, use the following commands in R:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install the package from GitHub
devtools::install_github("SAADAT-Abu/MethylBeta")
