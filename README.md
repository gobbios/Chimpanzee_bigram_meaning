
This folder contains files to reproduce the analysis presented in the manuscript.
Since running the entire analysis is somewhat time consuming we modified the scripts such that by default we work with a small subset of data and fit the model with fewer iterations.
As such the scripts *illustrate* the analysis rather than allowing *replicating* it.

If you want to fit the complete model you need to modify two scripts (`00_data_prep.R`: deactivate subsetting of data set at approximately line 15; `01_model_fit.R`: fit model with large number of iterations at approximately line 10).


# Prerequisites

We assume that all the files listed below are present in your current working directory (including the folder `helpers`).

The only other R package required is the `cmdstanr` package (see https://mc-stan.org/cmdstanr/ for installation instructions).


# Files

## R scripts

- `00_data_prep.R`: data reading and wrangling

- `01_model_fit.R`: model fit

- `02_figures.R`: figures

- `03_tables.R`: tables

- helper functions for data wrangling, results processing and plotting 

  - `helpers/context_distribution_plot.R`
  
  - `helpers/distance_plot.R`
  
  - `helpers/pairwise_dist_from_stan.R`
  
  - `helpers/prep_standat.R`

## Data

- `calldata.csv`: raw data file

## Model

- `stanmodel.stan`: model code for the Stan model

