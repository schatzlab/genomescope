#' Function to fit 2p peak model, with p forms
#'
#' @param kmer_hist_orig A data frame of the original histogram data (starting at 1 and with last position removed).
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @param y A numeric vector of the y-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @param k An integer corresponding to the kmer length.
#' @param p An integer corresponding to the ploidy.
#' @param topology An integer corresponding to the topology to use.
#' @param estKmercov An integer corresponding to the estimated kmer coverage of the polyploid genome.
#' Set to -1 if not specified by user.
#' @param round An integer corresponding to the iteration number (0, 1, 2, 3) for the fitting process.
#' @param foldername A character vector corresponding to the name of the output directory.
#' @param arguments A data frame of the user-specified inputs.
#' @return A list (nls, nlsscore) where nls is the nlsLM model object (with some additional components)
#' and nlsscore is the score (model RSSE) corresponding to the best fit (of the p forms).
#' @export
estimate_Genome_peakp<-function(kmer_hist_orig, x, y, k, p, topology, estKmercov, round, foldername, arguments) {
  if (topology==-1) {
    p_to_num_topologies = c(1, 1, 1, 2, 5, 16)
    num_topologies = p_to_num_topologies[p]
    topologies = 1:num_topologies
  }
  else {
    num_topologies = 1
    topologies = c(topology)
  }
  numofKmers = sum(as.numeric(x)*as.numeric(y))
  if (estKmercov==-1) {
    #In situations with low heterozygosity, the peak with highest amplitude typically corresponds to the homozygous peak (i.e. the p-th peak).
    #However, with increasing heterozygosity, the highest amplitude peak may be an earlier peak.
    #Thus, when setting the estimated kmer coverage, we will need to iterate through these possibilities.
    #num_peak_indices indicates how many possibilities we need to iterate through.
    num_peak_indices = p
    y_transform = as.numeric(x)**transform_exp*as.numeric(y)
    estKmercov1 = x[which(y_transform==max(y_transform))][1]
  }
  else {
    # When the user sets the estimated kmer coverage, we only need to iterate through one possibility
    num_peak_indices = 1
    ## We set the estimated kmer coverage to be the user specified value
    estKmercov1 = estKmercov
  }
  estLength1 = numofKmers/estKmercov1

  nls00 = NULL
  peak_indices = 1:num_peak_indices
  for (i in peak_indices) {
    nls0 = NULL
    top_count = 0
    ## We see what happens when we set the estimated kmer coverage to be 1/i times the x-coordinate where the max peak occurs (1 <= i <= p if the user doesn't set the estimated kmer coverage, and i=1 if they do)
    estKmercov2 = estKmercov1 / i
    estLength2 = numofKmers/estKmercov2

    if (VERBOSE) {cat(paste("trying with kmercov: ", estKmercov2, "\n"))}

    for (top in topologies) {
      if (VERBOSE) {cat(paste("trying with topology: ", top, "\n"))}
      top_count = top_count + 1
      nls1 = nls_peak(x, y, k, p, top, estKmercov2, estLength2, MAX_ITERATIONS)
      nls0 = eval_model(kmer_hist_orig, nls0, nls1, p, round, foldername, arguments)[[1]]
    }
    if (i < num_peak_indices) { #if this is not the last evaluation
      nls00 = eval_model(kmer_hist_orig, nls00, nls0, p, round, foldername, arguments)[[1]]
    }
  }

  return(eval_model(kmer_hist_orig, nls00, nls0, p, round, foldername, arguments))
}
