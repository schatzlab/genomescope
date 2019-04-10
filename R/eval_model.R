#' Evaluate distinct model forms, in order to resolve ambiguity of which peak is the homozygous peak
#'
#' @param kmer_hist_orig A data frame of the original histogram data (starting at 1 and with last position removed).
#' @param nls0,nls1 The nlsLM model objects to evaluate and compare.
#' @param p An integer corresponding to the ploidy.
#' @param round An integer corresponding to the iteration number (0, 1, 2, 3) for the fitting process.
#' @param foldername A character vector corresponding to the name of the output directory.
#' @param arguments A data frame of the user-specified inputs.
#' @return A list (nls, nlsscore) where nls is the nlsLM model object (with some additional components)
#' and nlsscore is the score (model RSSE) corresponding to the best fit (of the p forms).
#' @export
eval_model<-function(kmer_hist_orig, nls0, nls1, p, round, foldername, arguments) {
  nls0score = -1
  nls1score = -1

  ## Evaluate the score the nls0
  if (!is.null(nls0)) {
    nls0score = score_model(kmer_hist_orig, nls0, round+0.1, foldername)

    #if(VERBOSE) {cat(paste("nls0score$all:\t", nls0score$all[[1]], "\n"))}

    if (VERBOSE) {
      mdir = paste(foldername, "/round", round, ".1", sep="")
      dir.create(mdir, showWarnings=FALSE)
      report_results(kmer_prof_orig,kmer_prof_orig, k, p, (list(nls0, nls0score)) , mdir, arguments, TRUE)
    }
  }
  else {
    if (VERBOSE) {cat("nls0score failed to converge\n")}
  }

  ## Evaluate the score of nls1
  if (!is.null(nls1)) {
    nls1score = score_model(kmer_hist_orig, nls1, round+0.2, foldername)

    if(VERBOSE) {cat(paste("nls1score$all:\t", nls1score$all[[1]], "\n"))}

    if (VERBOSE) {
      mdir = paste(foldername, "/round", round, ".2", sep="")
      dir.create(mdir, showWarnings=FALSE)
      report_results(kmer_prof_orig, kmer_prof_orig, k, p, (list(nls1, nls1score)) , mdir, arguments, FALSE)
    }
  }
  else {
    if (VERBOSE) {cat("nls1score failed to converge\n")}
  }

  ## Return the better of the scores
  if (!is.null(nls0)) {
    if (!is.null(nls1)) {
      pdiff = abs(nls0score$all[[1]] - nls1score$all[[1]]) / max(nls0score$all[[1]], nls1score$all[[1]])

      if (pdiff < SCORE_CLOSE) {
        het0 = nls0$ahet
        het1 = nls1$ahet

        #if (het1 * SCORE_HET_FOLD_DIFFERENCE < het0) {
        if (het1 + 0.01 < het0) {
          if (VERBOSE) {cat("returning nls0, similar score, higher het\n")}
          return (list(nls1, nls1score))
        }
        #else if (het0 * SCORE_HET_FOLD_DIFFERENCE < het1) {
        else if (het0  + 0.01 < het1) {
          if (VERBOSE) {cat("returning nls1, similar score, higher het\n")}
          return (list(nls0, nls0score))
        }
      }

      if (nls0score$all[[1]] < nls1score$all[[1]]) {
        if (VERBOSE) {cat("returning nls0, better score\n")}
        return (list(nls0, nls0score))
      }
      else {
        if (VERBOSE) {cat("returning nls1, better score\n")}
        return (list(nls1, nls1score))
      }
    }
    else {
      if (VERBOSE) {cat("returning nls0, nls1 fail\n")}
      return (list(nls0, nls0score))
    }
  }

  if (VERBOSE) {cat("returning nls1 by default\n")}
  return (list(nls1, nls1score))
}
