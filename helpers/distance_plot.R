#' create distance plot of bigram to all individual/single calls
#'
#' @param stanmod,standat cmdstanr model environment and data list used to fit
#'        the model
#' @param targets_single character of length 2, indicating two calls to be
#'        visually highlighted (typically the two constituent calls
#'        from \code{target_bigram})
#' @param target_bigram character, the bigram from which distances are displayed
#' @param target_calls
#' @param emp_distances
#' @param xlims
#' @param add_legend
#' @param add_cor_hist
#' @param return_post_mat
#'
#' @return a plot

distance_plot <- function(stanmod,
                          standat,
                          targets_single,
                          target_bigram,
                          target_calls = NULL,
                          emp_distances = NULL,
                          xlims = c(0, 1.4),
                          add_legend = TRUE,
                          add_cor_hist = FALSE,
                          return_post_mat = FALSE,
                          return_pair_support = FALSE,
                          return_prop_support = FALSE,
                          cex.legend = 0.8,
                          cex.p.support = 0.7,
                          cex.text = 0.6,
                          cex.xlab = 1,
                          cex.title = 1
                          ) {
  if (add_cor_hist) xlims[1] <- -0.5

  if (is.null(target_calls)) {
    if (isTRUE(standat$for_multinom == 1)) {
      target_calls <- unique(c(names(standat$lab_call_1), names(standat$lab_call_2)))
      if (target_bigram %in% target_calls) {
        target_calls <- target_calls[-c(which(target_calls == target_bigram))]
      }

    } else {
      target_calls <- unique(c(names(standat$lab_call_1), names(standat$lab_call_2)))
      target_calls <- target_calls[nchar(target_calls) == 2]
    }

    target_calls <- target_calls[!grepl("_", target_calls)]
  }
  # distances between pairs (takes a second or two)...
  pairwisedists <- sapply(target_calls, FUN = pairwise_dist_from_stan,
                          stanfit = stanmod, standat = standat,
                          call2 = target_bigram,
                          print_pair_code = FALSE, simplify = FALSE)

  xx <- do.call("cbind", pairwisedists)
  xranks <- t(apply(xx, 1, rank))

  if (return_post_mat) {
    return(xranks)
  }

  #### create a matrix with the essential information to create a plot
  xtab <- data.frame(call = target_calls, from = target_bigram)
  # colnames(xtab) <- c('call', 'mean_dist', 'low95', 'high95', 'low89', 'high89', 'low50', 'high50')
  xtab$call <- target_calls
  xtab$mean_dist <- unlist(lapply(pairwisedists, mean))
  xtab[, c('low95', 'high95', 'low89', 'high89', 'low50', 'high50')] <- do.call(
    "rbind",
    lapply(pairwisedists, quantile, prob = c(0.025, 0.975, 0.055, 0.945, 0.25, 0.75))
  )
  # sort according to mean
  xtab <- xtab[order(xtab$mean_dist, decreasing = TRUE), ]
  # add formatting
  xtab$xcol <- "black"
  xtab$xcol[xtab$call %in% targets_single] <- "darkorange"
  xtab$font <- 1
  xtab$font[xtab$call %in% targets_single] <- 2
  xtab$vertcoord <- seq_len(nrow(xtab))


  plot(xtab$mean_dist, xtab$vertcoord, col = "black", cex = 2.5,
       axes = FALSE, ann = FALSE, pch = 18, xlim = xlims, xaxs = "i", ylim = c(0.5, nrow(xtab) + 1.5))

  mtext(text = xtab$call, side = 2, line = 0.5, at = xtab$vertcoord, las = 1, col = xtab$xcol, font = xtab$font, cex = 0.9)
  segments(x0 = xtab$low89, y0 = xtab$vertcoord, x1 = xtab$high89, y1 = xtab$vertcoord, col = 'darkblue', lwd=2)
  segments(x0 = xtab$low50, y0 = xtab$vertcoord, x1 = xtab$high50, y1 = xtab$vertcoord, col = 'cyan3', lwd=4)
  points(xtab$mean_dist, xtab$vertcoord, col = "black", cex = 2.5, pch =18)
  axis(1, cex.axis = cex.xlab)
  box()
  if (add_legend) {
    legend("topright", legend = c("", "posterior mean", "50% CI", "89% CI"),
           lty = c(NA, NA, 1, 1), lwd = c(NA, NA, 4, 2),
           col = c("white", "black", "cyan3", "darkblue"),
           pt.cex = 2, cex = cex.legend, bty = "n", title = "", pch = c(NA, 18, NA, NA))
  }
  if (!is.null(emp_distances)) {
    emp <- emp_distances[, target_bigram]
    points(emp[xtab$call], xtab$vertcoord, pch = 4, col = "red")
  }

  title(main = paste("Euclidian distances to", target_bigram),
        cex.main = cex.title)

  if (add_cor_hist) {
    pdata <- cor_plot_data(stanmod = stanmod,
                           standat = standat,
                           target_bigram = target_bigram,
                           target_calls = xtab$call)
    for (i in xtab$vertcoord) {
      aux <- pdata[[xtab$call[i]]]
      aux <- table(round(aux, 2))
      aux <- aux / max(aux)
      aux <- aux * 0.6
      aux <- cbind(xval = as.numeric(names(aux)), freq = aux, yval = i)
      aux <- cbind(aux, col = as.numeric(aux[, "xval"] > 0) + 1)
      aux[, "xval"] <- aux[, "xval"] / 4
      aux[, "xval"] <- aux[, "xval"] - 0.25
      # aux[, "xval"] <- aux[, "xval"] + 0.5
      segments(x0 = aux[, "xval"], y0 = aux[, "yval"],
               x1 = aux[, "xval"], y1 = aux[, "yval"] + aux[, "freq"],
               col = c("lightgrey", "red")[aux[, "col"]], xpd = TRUE)
    }

  }

  xtext <- paste0("prop samples where ", shQuote(targets_single[1]),
                  " and ", shQuote(targets_single[2]),
                  " have smallest distance: ", round(mean(rowSums(xranks[, targets_single]) <= 3), 3))

  legend("topleft",
         legend = xtext,
         bty = "n",
         cex = cex.text)

  abline(h = nrow(xtab) + 0.6)

  xtab$closest_prop <- colMeans(xranks == 1)[xtab$call]
  plot_props <- xtab$closest_prop
  plot_props[plot_props < 0.01] <- NA
  plot_props <- sprintf("%.3f", plot_props)
  plot_props[plot_props == "NA"] <- ""
  text(0.01, xtab$vertcoord, plot_props, cex = cex.p.support, adj = 0)


    xtab

  if (return_prop_support) return (xtab)
  if (return_pair_support) return(mean(rowSums(xranks[, targets_single]) <= 3))
}



cor_plot_data <- function(stanmod, standat, target_bigram, target_calls = NULL) {
  if (is.null(target_calls)) {
    target_calls <- unique(c(names(standat$lab_call_1), names(standat$lab_call_2)))
    target_calls <- target_calls[nchar(target_calls) == 2]
  }

  dmat <- stanmod$draws("pred_mat_rep_multiplied", format = "draws_matrix")
  # row index for target_bigram
  i <- standat$label_calltype[target_bigram]
  tmat <- dmat[, grepl(paste0("[", i, ","), colnames(dmat), fixed = TRUE)]
  # for correlations: https://stackoverflow.com/questions/27943070/row-wise-correlations-in-r
  cA <- tmat - rowMeans(tmat)
  sA <- sqrt(rowMeans(cA^2))

  out <- vector(mode = "list", length = length(target_calls))
  names(out) <- target_calls
  i = target_calls[1]
  for (i in target_calls) {
    index <- standat$label_calltype[i]
    tempmat <- dmat[, grepl(paste0("[", index, ","), colnames(dmat), fixed = TRUE)]
    cB <- tempmat - rowMeans(tempmat)
    sB <- sqrt(rowMeans(cB^2))
    out[[i]] <- as.numeric(rowMeans(cA * cB) / (sA * sB))
  }
  out
}
