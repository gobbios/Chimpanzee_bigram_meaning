library(cmdstanr)

# require object 'standat', which is created in the script '00_data_prep.R'
# source("00_data_prep.R")

m <- cmdstanr::cmdstan_model("stanmodel.stan")
m$check_syntax(pedantic = TRUE)

# model with 4000 iterations
# s <- m$sample(data = standat, parallel_chains = 4, chains = 8, seed = 42, iter_sampling = 500, iter_warmup = 2000, refresh = 50, adapt_delta = 0.85)

# model with small number of iterations (188)
s <- m$sample(data = standat, parallel_chains = 4, chains = 4, seed = 47, iter_sampling = 47, iter_warmup = 43, refresh = 5, adapt_delta = 0.8)


# save model output to disk
# creates one csv file per chain
s$save_output_files(basename = "modres", timestamp = FALSE, random = FALSE)


