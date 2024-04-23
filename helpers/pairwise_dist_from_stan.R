#' extract Euclidean distance between two call types from a Stan model fit
#'
#' @param stanfit,standat cmdstanr model environment and data list used to fit
#'        the model
#' @param call1,call2 character, call type names for which the distance is to
#'        be extracted
#' @param print_pair_code logical, print some information to the screen
#'
#' @return a vector with pair-wise distances for the designated pair of calls

pairwise_dist_from_stan <- function(stanfit,
                                    standat,
                                    call1,
                                    call2,
                                    print_pair_code = TRUE
                                    ) {
  sel <- c(which(names(standat$lab_call_1) == call1 & names(standat$lab_call_2) == call2),
           which(names(standat$lab_call_1) == call2 & names(standat$lab_call_2) == call1))
  if (print_pair_code) {
    cat(names(standat$lab_call_pair)[sel], "; which is pair number",
        shQuote(sel), "in standat$lab_call_pair", "\n")
  }
  as.numeric(stanfit$draws(variables = paste0("distvals[", sel, "]"),
                           format = "draws_matrix"))
}

