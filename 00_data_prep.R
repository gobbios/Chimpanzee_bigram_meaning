# A) read data from csv file ----
xdata <- read.csv("calldata.csv", stringsAsFactors = FALSE)
summary(xdata)
# 4322 observations
# length(unique(xdata$caller)) # 53 unique callers
# length(unique(xdata$call.name)) # 27 unique calls (single utterances or call combinations)
# length(unique(xdata$context)) # 99 unique contexts

head(xdata)

## select small subset of data for demonstration of how the code works
##   this speeds up model fitting and postprocessing substantially
## if you want to replicate the entire analysis: ignore the following line
xdata <- droplevels(xdata[xdata$demosubset == 1, ])




# B) reshape for modelling ------

## create response matrix (cols = contexts, rows = observation)
cn_nm <- sort(unique(unlist(strsplit(xdata$context, "_"))))
ymat <- matrix(ncol = length(cn_nm), nrow = nrow(xdata), 0)
colnames(ymat) <- cn_nm
rownames(xdata) <- NULL
for (i in seq_len(nrow(xdata))) {
  ymat[i, unlist(strsplit(xdata$context[i], "_"))] <- 1
}

head(ymat)

## add trials and response matrix to xdata
xdata$Ytrials <- rowSums(ymat)
xdata$Ymat <- ymat

# sanity check
all(rowSums(ymat) == xdata$context.length) # should be TRUE



# C) split ('dyadic') calls into their individual calls ----
call_by_dyad <- lapply(strsplit(xdata$call.name, "_"), function(x) {
  if (length(x) == 1) {
    return(rep(x[1], 2))
  } else {
    return(x)
  }
})
call_by_dyad <- t(do.call("cbind", call_by_dyad))
xdata$call1 <- call_by_dyad[, 1]
xdata$call2 <- call_by_dyad[, 2]

## factorize for further processing
xdata$caller <- as.factor(xdata$caller)
xdata$context <- as.factor(xdata$context)
xdata$call.name <- as.factor(xdata$call.name)


# D) prepare data for Stan model ----
source("helpers/prep_standat.R")
xdata <- droplevels(xdata)
standat <- prep_standat(xdata = xdata, col_context = "context", col_calltype = "call.name", col_id = "caller")



# E) sanity checks ----
all(names(standat$label_calltype)[standat$index_calltype] == xdata$call.name)
all(xdata$caller == levels(xdata$caller)[standat$index_subject])
all(xdata$Ymat == standat$Ymat)
all(colnames(xdata$Ymat) == colnames(standat$Ymat))

