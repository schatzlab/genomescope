#' Function to fit 2p peak model, with p forms
#'
#' @param kmer_hist_orig A data frame of the original histogram data (starting at 1 and with last position removed).
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @param y A numeric vector of the y-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @param k An integer corresponding to the kmer length.
#' @param p An integer corresponding to the ploidy.
#' @param estKmercov An integer corresponding to the estimated kmer coverage of the polyploid genome.
#' Set to -1 if not specified by user.
#' @param round An integer corresponding to the iteration number (0, 1, 2, 3) for the fitting process.
#' @param foldername A character vector corresponding to the name of the output directory.
#' @param arguments A data frame of the user-specified inputs.
#' @return A list (nls, nlsscore) where nls is the nlsLM model object (with some additional components)
#' and nlsscore is the score (model RSSE) corresponding to the best fit (of the p forms).
#' @export
estimate_Genome_peakp<-function(kmer_hist_orig, x, y, k, p, estKmercov, round, foldername, arguments) {
  p_to_num_topologies = c(1, 1, 1, 2, 5, 15)
  num_topologies = p_to_num_topologies[p]
  numofKmers = sum(as.numeric(x)*as.numeric(y))
  if (estKmercov==-1) {
    ## First we see what happens when we set the estimated kmer coverage to be the x-coordinate where the max peak occurs (typically the homozygous peak)
    estKmercov1  = x[which(y==max(y))][1]
  }
  else {
    ## We set the estimated kmer coverage to be the user specified value
    estKmercov1 = estKmercov
  }
  estLength1 = numofKmers/estKmercov1

  if (VERBOSE) {cat(paste("trying with kmercov: ", estKmercov1, "\n"))}

  nls0 = NULL
  for (top in 1:num_topologies) {
    nls1 = nls_peak(x, y, k, p, top, estKmercov1, estLength1, MAX_ITERATIONS)
    if (top < num_topologies || (estKmercov==-1 && p>=2)) { #if this is not the last evaluation
      nls0 = eval_model(kmer_hist_orig, nls0, nls1, p, round, foldername, arguments)[[1]]
    }

    if (VERBOSE) {print(summary(nls0))}

    if (estKmercov==-1 && p>=2) {
      ploidies = 2:p
      for (i in ploidies) {
        ## Next we see what happens when we set the estimated kmer coverage to be 1/i times the x-coordinate where the max peak occurs (2 <= i <= p)
        estKmercov2  = estKmercov1 / i
        estLength2 = numofKmers/estKmercov2

        if (VERBOSE) {cat(paste("trying with kmercov: ", estKmercov2, "\n"))}

        nls1 = nls_peak(x, y, k, p, top, estKmercov2, estLength2, MAX_ITERATIONS)

        if (VERBOSE) {print(summary(nls1))}

        if (i<p || top < num_topologies) { #if this is not the last evaluation
          nls0 = eval_model(kmer_hist_orig, nls0, nls1, p, round, foldername, arguments)[[1]]
        }
      }
      #nls0 = eval_model(kmer_hist_orig, nls0, nls1, p, round, foldername, arguments)[[1]]
    }
#    else {
#      return(eval_model(kmer_hist_orig, nls0, nls0, p, round, foldername, arguments))
#    }
#    
  }
  return(eval_model(kmer_hist_orig, nls0, nls1, p, round, foldername, arguments))
}
