#' convert data table with data to Stan-ready list
#'
#' @param xdata the data frame
#' @param col_context,col_calltype,col_id character, column names for context, call type, and caller
#'
#' @return a list

prep_standat <- function(xdata, col_context, col_calltype, col_id) {

  x <- droplevels(xdata)
  x[, col_context] <- as.factor(x[, col_context])
  x[, col_calltype] <- as.factor(x[, col_calltype])
  x[, col_id] <- as.factor(x[, col_id])
  x <- droplevels(x)


  ymat <- xdata$Ymat

  dat <- list(
    N = nrow(xdata),
    ncat = ncol(ymat),
    Ymat = ymat,
    Ytrials = rowSums(ymat),
    index_calltype = as.numeric(as.factor(x[, col_calltype])),
    n_calltypes = length(unique(x[, col_calltype])),

    index_subject = as.numeric(as.factor(x[, col_id])),
    n_subjects = length(unique(x[, col_id]))
  )

  names(dat$index_subject) <- as.character(as.factor(x[, col_id]))
  names(dat$index_calltype) <- as.character(as.factor(x[, col_calltype]))
  # add context labels for indexing
  lab4 <- seq_len(ncol(ymat))
  names(lab4) <- colnames(ymat)
  dat$label_context <- lab4
  # add calltype labels for indexing
  lab5 <- seq_len(nlevels(x[, col_calltype]))
  names(lab5) <- levels(x[, col_calltype])
  dat$label_calltype <- lab5

  dat$for_multinom <- 1

  # indices for combined distributions
  bigram_index <- grep("_", names(dat$label_calltype))
  combined <- matrix(ncol = 3, nrow = length(bigram_index), 0)
  colnames(combined) <- c("bigram_index", "call1", "call2")
  combined[, "bigram_index"] <- bigram_index
  aux <- unlist(lapply(strsplit(names(dat$label_calltype)[bigram_index], "_"), function(X)X[1]))
  combined[, "call1"] <- as.integer(dat$label_calltype[aux])
  aux <- unlist(lapply(strsplit(names(dat$label_calltype)[bigram_index], "_"), function(X)X[2]))
  combined[, "call2"] <- as.integer(dat$label_calltype[aux])
  dat$index_combined <- combined
  dat$n_combined <- nrow(combined)

  # helper objects (indices and labels)
  # pairwise calltypes (-> for distances)
  ct <- levels(x[, col_calltype])
  aux <- combn(ct, 2)
  dat$index_comps <- t(apply(aux, 2, function(x)sapply(x, function(y)which(ct == y))))
  dat$n_distcomps <- nrow(dat$index_comps)
  lab1 <- rep(0, dat$n_distcomps)
  names(lab1) <- aux[1, ]
  lab2 <- rep(0, dat$n_distcomps)
  names(lab2) <- aux[2, ]
  lab3 <- rep(0, dat$n_distcomps)
  names(lab3) <- paste(aux[1, ], aux[2, ], sep = "_@_")
  dat$lab_call_1 <- lab1
  dat$lab_call_2 <- lab2
  dat$lab_call_pair <- lab3


  dat$length <- dat$Ytrials # same as Ytrials

  dat
}

