# MethylBeta | DMP Analysis Package

![methyBeta](https://github.com/user-attachments/assets/5af6936b-de90-4b86-b98b-604c6e487f17)

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
