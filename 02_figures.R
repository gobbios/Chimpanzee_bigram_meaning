# we need `xdata` and `standat` again, e.g. uncomment:
# source("00_data_prep.R")

# reread model results if not already in workspace
s <- cmdstanr::as_cmdstan_fit(list.files(pattern = "modres"))

# make helper functions available
source("helpers/distance_plot.R")
source("helpers/context_distribution_plot.R")
source("helpers/pairwise_dist_from_stan.R")


# exemplifying primary results plots with one bigram (HO_PN, figure 3 in main text)

# create a figure sized one A4 page wide and one 1/3 of an A4 page high (11.7 * 8.3 inches)
# this would correspond to approximately one line in figure 3 or 4

# note that the figure that results from the small data subset differs quite
#   dramatically from the one using the full data set in the manuscript
#   because it is based on fewer events, call types and callers
# (it just illustrates how the code works)

pdf(file = "fig3.pdf", width = 11.7, height = 2.8)
par(family = "serif", mfrow = c(1, 4))

context_distribution_overids(which_call = "HO", xdata = xdata,
                             stanmod = s, standat = standat, plot_model = FALSE)
context_distribution_overids(which_call = "PN", xdata = xdata,
                             stanmod = s, standat = standat, plot_model = FALSE)
context_distribution_overids(which_call = "HO_PN", xdata = xdata,
                             stanmod = s, standat = standat, plot_model = FALSE)

distance_plot(stanmod = s, standat = standat,
              target_bigram = "HO_PN",
              targets_single = c("HO", "PN"),
              xlim = c(0, 2))

dev.off()

