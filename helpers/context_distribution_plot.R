#' plot context distribution of a specific call type or a call combination across contexts
#'
#' @param xdata data.frame with calling data (needs to include either
#'        \code{call.name} or \code{calltype} column). Also requires
#'        \code{caller}, \code{context}, \code{context.length} columns.
#' @param which_call character, what call?
#' @param ylims numeric of length 2, limits of y axis (default is
#'        \code{0, 1.15} to allow space for legend at the top)
#' @param add_mean logical, if \code{TRUE}: plot the mean (as grey cross)
#' @param stanmod,standat optional cmdstanr model environment and data list,
#'        default is \code{NULL} for both. If both are supplied, posterior
#'        medians and 90% CI are plotted alongside the observed data.
#' @param use_color logical, if \code{TRUE}: use colors for contexts
#' @param return_post_mat return the plotting data visibly
#' @param return_full_post return the full posterior for all contexts, this
#'        takes priority over \code{return_post_mat=TRUE}
#'
#' @return a plot and (invisibly) a matrix with the underlying data
#'         generated for the plot


context_distribution_overids <- function(xdata,
                                         which_call,
                                         ylims = c(0, 1.15),
                                         add_mean = FALSE,
                                         stanmod = NULL,
                                         standat = NULL,
                                         use_color = TRUE,
                                         cex.axis1 = 1,
                                         cex.axis2 = 1,
                                         cex.text = 1,
                                         cex.ylab = 1,
                                         cex.title = 1.3,
                                         pt.cex.data = 2,
                                         lwd.data = 2,
                                         pt.cex.mod = 0.8,
                                         lwd.mod = 1,
                                         return_post_mat = FALSE,
                                         return_full_post = FALSE,
                                         plot_model = TRUE
                                         ) {

  # return full posterior before starting plotting...
  if (return_full_post) {
    if (is.null(stanmod) || is.null(standat)) {
      stop("can't return posterior samples unless model and standata are ",
           "supplied via 'stanmod = <...>' and 'standat = <...>",
           call. = FALSE)
    }
    post_full <- stanmod$draws("pred_mat_rep_multiplied", format = "draws_matrix")
    wc_num <- which(names(standat$label_calltype) == which_call)
    post_full <- post_full[, grepl(paste0("[", wc_num, ","), colnames(post_full), fixed = TRUE)]
    post_full <- matrix(post_full, nrow = nrow(post_full))
    colnames(post_full) <- names(standat$label_context)
    return(post_full)
  }

  callname <- c("call.name", "calltype")[which(c("call.name", "calltype") %in% colnames(xdata))]

  # all contexts:
  all_contexts <- sort(unique(unlist(strsplit(as.character(xdata$context), "_"))))

  # define colors if so desired
  # defaults:
  context_col <- rep("black", length(all_contexts))
  mean_col <- rep("lightgrey", length(all_contexts))
  estimate_col <- rep("red", length(all_contexts))

  if (use_color) {
    # Cedric's color universe
    col_cats <- c("deeppink1","lightpink","lightblue", "blue", "blueviolet", "brown", "brown2",
                  "coral1", "darkgreen", "chartreuse1",  "gold1",
                  "gold4", "aquamarine")
    context_col <- colorRampPalette(col_cats)(length(all_contexts))
    estimate_col <- rep("black", length(all_contexts))
  }

  aux <- xdata[xdata[, callname] == which_call, ]
  tot <- nrow(aux)

  # account properly for events with with more than two events
  xx <- data.frame(do.call("rbind", apply(as.matrix(aux[, c("caller", "context", "context.length")]), 1, function(x) {
    temp1 <- unlist(lapply(x[2], strsplit, split = "_"))
    temp2 <- rep(x[1], length(temp1))
    temp3 <- rep(x[3], length(temp1))
    cbind(caller = temp2, context = temp1, context.length = temp3)
  }, simplify = FALSE)))
  xx$context.length <- as.numeric(as.character(xx$context.length))
  # xx$context <- factor(xx$context, levels = all_contexts)
  temp <- xx


  temp$context <- factor(temp$context, levels = all_contexts)
  xtab <- table(temp$caller, temp$context)
  prop_table <- table(temp$caller, temp$context)
  # getting proportions that sum to 1 per individual
  prop_table <- prop_table / rowSums(prop_table)
  # that's what the multinom model gives us

  # multiply (per individual!) with the ratio of events and number of recordings
  i=3
  for (i in seq_len(nrow(prop_table))) {
    n_recordings <- sum(aux$caller == rownames(prop_table)[i])
    n_events <- sum(xtab[i, ])
    prop_table[i, ] <- prop_table[i, ] * (n_events/n_recordings)
  }
  round(prop_table, 2)
  didoccur <- colSums(prop_table) > 0
  ninds <- nrow(prop_table)
  pdata <- apply(prop_table, 2, quantile, c(0.25, 0.5, 0.75), na.rm = TRUE)

  # sanity checks...
  if (any(prop_table > 1)) message("there were actually proportion values larger than one after multiplication")
  if (any(pdata > ylims[2])) warning("values outside of plotting range", call. = FALSE)

  plot(0, 0, type = "n", xlim = c(0.5, ncol(pdata) + 0.5), ylim = ylims, yaxs = "i",
       axes = FALSE, xlab = "", ylab = "")


  title(ylab = "Proportion", cex.lab = cex.ylab)
  title(main = which_call, cex.main = cex.title)

  # mtext('Proportion', side = 2, line = 3, cex = cex.ylab)
  # mtext(which_call, side = 3, line = 1, cex = cex.title, font = 2)

  axis(1, at = seq_len(ncol(pdata)), labels = FALSE, cex.axis = 1)

  text(seq_len(ncol(pdata)),
       par("usr")[3] - 1/20, labels = colnames(pdata), adj=1, srt = 45, xpd = NA,cex = cex.axis1)


  box()
  abline(h = 1)

  points(seq_len(ncol(pdata))[didoccur], pdata[2, didoccur], pch = 16, xpd = TRUE,
         cex = pt.cex.data, col = context_col[didoccur])
  segments(seq_len(ncol(pdata))[didoccur], pdata[1, didoccur],
           seq_len(ncol(pdata))[didoccur], pdata[3, didoccur],
           lwd = lwd.data, col = context_col[didoccur])

  pdata <- rbind(pdata, "mean" = colMeans(prop_table))
  pdata <- rbind(pdata, "est_post_median" = NA, "post_lower" = NA, "post_upper" = NA)
  if (add_mean) {
    points(seq_len(ncol(pdata))[didoccur], pdata["mean", ][didoccur],
           pch = 4, col = mean_col, cex = 0.5, lwd = 2)
  }

  axis(2, las = 1, cex.axis = cex.axis2)
  legend("topleft",
         legend = c(bquote(italic(N)[calls] == .(tot)),
                    bquote(italic(N)[individuals] == .(ninds))),
         bty = "n", ncol = 2, cex = cex.text)


  # modelled estimates if model and data are supplied
  if (!is.null(stanmod) && !is.null(standat)) {
    xx <- stanmod$summary("pred_mat_rep_multiplied", ~quantile(.x, probs = c(0.055, 0.5, 0.945)))
    x1 <- matrix(data.frame(xx)[, "X50."], ncol = standat$ncat, byrow = FALSE)
    rownames(x1) <- names(standat$label_calltype)
    colnames(x1) <- names(standat$label_context)
    pp <- x1[which_call, , drop = FALSE]
    pdata["est_post_median", ] <- pp[, colnames(pdata)]
    x1 <- matrix(data.frame(xx)[, "X5.5."], ncol = standat$ncat, byrow = FALSE)
    rownames(x1) <- names(standat$label_calltype)
    colnames(x1) <- names(standat$label_context)
    pp <- x1[which_call, , drop = FALSE]
    pdata["post_lower", ] <- pp[, colnames(pdata)]
    x1 <- matrix(data.frame(xx)[, "X94.5."], ncol = standat$ncat, byrow = FALSE)
    rownames(x1) <- names(standat$label_calltype)
    colnames(x1) <- names(standat$label_context)
    pp <- x1[which_call, , drop = FALSE]
    pdata["post_upper", ] <- pp[, colnames(pdata)]

    if (plot_model == TRUE) {
      points(seq_len(ncol(pdata)) + 0.1, pdata["est_post_median", ],
             pch = 8, col = estimate_col, cex = pt.cex.mod, lwd = lwd.mod)
      segments(seq_len(ncol(pdata)) + 0.1, pdata["post_lower", ],
               seq_len(ncol(pdata)) + 0.1, pdata["post_upper", ],
               col = estimate_col, lwd = 0.5)
    }


  }

  if (return_post_mat) {
    return(pdata)
  }

  invisible(pdata)
}
