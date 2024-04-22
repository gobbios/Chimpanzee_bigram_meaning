# we need `xdata` and `standat` again, e.g. uncomment:
# source("00_data_prep.R")

# reread model results if not already in workspace
s <- cmdstanr::as_cmdstan_fit(list.files(pattern = "modres"))

# make helper functions available
source("helpers/distance_plot.R")
source("helpers/context_distribution_plot.R")
source("helpers/pairwise_dist_from_stan.R")

# recreate table 3
# note that the small data subset only contains 5 bigrams

bigrams <- as.character(sort(unique(xdata$call.name[xdata$combi.length==2])))
res <- data.frame(bigram = bigrams)


for (i in seq_len(nrow(res))) {
  xunits <- unlist(strsplit(as.character(res$bigram[i]), split = "_"))

  tempres1 <- distance_plot(stanmod = s, standat = standat,
                            targets_single = xunits, target_bigram = res$bigram[i],
                            return_pair_support = FALSE, return_prop_support = TRUE)
  tempres2 <- distance_plot(stanmod = s, standat = standat,
                            targets_single = xunits, target_bigram = res$bigram[i],
                            return_pair_support = TRUE)

  res$closest1[i] <- tempres1$call[nrow(tempres1)]
  res$dist1[i] <- tempres1$mean_dist[nrow(tempres1)]
  res$p.support.closest1[i] <- tempres1$closest_prop[nrow(tempres1)]
  res$second.closest[i] <- tempres1$call[nrow(tempres1) - 1]
  res$dist2[i] <- tempres1$mean_dist[nrow(tempres1) - 1]
  res$p.support.second[i] <- tempres1$closest_prop[nrow(tempres1) - 1]
  res$p.support.both.closest[i] <- tempres2
}

res


