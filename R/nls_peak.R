#' Uses nlsLM to fit 2p peak model
#'
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @param y A numeric vector of the y-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @param k An integer corresponding to the kmer length.
#' @param p An integer corresponding to the ploidy.
#' @param estKmercov A numeric with the estimated average kmer coverage of the polyploid genome.
#' @param estLength A numeric with the estimated polyploid genome length.
#' @param max_iterations An integer corresponding to the maximum number iterations to use for nlsLM.
#' @return An nlsLM model object with some additional components.
#' @export
nls_peak<-function(x, y, k, p, estKmercov, estLength, max_iterations) {
  #Initiate variables
  model = NULL
  d_min = 0
  d_initial = 0.001
  d_max = 1
  r_min = 0
  r1_initial = 0.001
  r2_initial = 0.002
  r3_initial = 0.003
  r4_initial = 0.004
  r5_initial = 0.005
  r_initials = c(r1_initial, r2_initial, r3_initial, r4_initial, r5_initial)
  r_start = vector("list", p-1)
  if (p > 1) {
    names(r_start) = paste("r", 1:(p-1), sep="")
    for (i in 1:(p-1)) {
      r_start[[paste("r",i,sep="")]] = r_initials[i]
    }
  }
  r_max = 1
  kmercov_min = 0
  kmercov_initial = estKmercov
  kmercov_max = Inf
  bias_min = 0
  bias_initial = 0.5
  bias_max = Inf
  length_min = 0
  length_initial = estLength
  length_max = Inf

  #Determine what formula to use, based on p
  if (p==1) {
    r_text = ""
  } else {
    r_text = paste(paste(lapply(1:(p-1), function(x) paste("r", as.character(x), sep="")), collapse=", "), ", ")
  }
  formula = as.formula(paste("y ~ length*predict",p,"(",r_text, "k, d, kmercov, bias, x)",sep=""))

  if (VERBOSE) {cat("trying nlsLM algorithm (Levenberg-Marquardt)\n")}

  try(model <- nlsLM(formula = formula,
                     start   = c(list(d = d_initial), r_start, list(kmercov = kmercov_initial, bias = bias_initial, length = length_initial)),
                     lower   = c(c(d_min), rep(r_min, p-1), c(kmercov_min, bias_min, length_min)),
                     upper   = c(c(d_max), rep(r_max, p-1), c(kmercov_max, bias_max, length_max)),
                     control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

  if (!is.null(model))
  {
    model_sum    = summary(model)
    model$p      = p
    if (p==1) {
      model$hets = list(c(0, 0))
    } else {
      model$hets = lapply(1:(p-1), function(x) min_max1(model_sum$coefficients[paste('r', x, sep=""),]))
    }
    model$het = c(1-Reduce("*", 1-unlist(lapply(model$hets, '[[', 1))), 1-Reduce("*", 1-unlist(lapply(model$hets, '[[', 2))))
    model$homo = 1-model$het
    model$dups   = min_max(model_sum$coefficients['bias',])
    model$kcov   = min_max(model_sum$coefficients['kmercov',])
    model$mlen   = min_max(model_sum$coefficients['length',])
    model$md     = min_max1(model_sum$coefficients['d',])
    if (p==1) {
      model$ahets = list(c(0))
    } else {
      model$ahets = lapply(1:(p-1), function(x) model_sum$coefficients[paste('r', x, sep=""),][[1]])
    }
    model$ahet = 1-Reduce("*", 1-unlist(model$ahets))
    model$ahomo = 1-model$ahet
    model$adups = model_sum$coefficients['bias',][[1]]
    model$akcov = model_sum$coefficients['kmercov',][[1]]
    model$amlen = model_sum$coefficients['length',][[1]]
    model$amd   = model_sum$coefficients['d',][[1]]
  }

  print(model)

  return(model)
}
