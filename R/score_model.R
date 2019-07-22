#' Score nlsLM model by number and percent of residual errors after excluding sequencing errors
#'
#' @param kmer_hist_orig A data frame of the original histogram data (starting at 1 and with last position removed).
#' @param nls An nlsLM model object.
#' @param round An integer corresponding to the iteration number (0, 1, 2, 3) for the fitting process.
#' @param foldername A character vector corresponding to the name of the output directory.
#' @return A data frame where the variables "all", "full", and "unique" correspond to the model residual sum of square error (excluding sequencing error)
#' including x-values up to the end of kmer_hist_orig, up to (2p+1)*lambda, and up to (p+1)*lambda respectively.
#' Additionally, the variables "allscore", "fullscore", and "uniquescore" are the
#' corresponding percentage of kmers that are correctly modeled.
#' @export
score_model<-function(kmer_hist_orig, nls, round, foldername) {
  x = kmer_hist_orig[[1]]
  y = kmer_hist_orig[[2]]
  y_transform = as.numeric(x)**transform_exp*as.numeric(y)

  pred=predict(nls, newdata=data.frame(x))
  model_sum=summary(nls)
  p=nls$p
  kcovfloor = max(1, floor(min_max(model_sum$coefficients['kmercov',])[[1]]))

  ## Compute error rate, by counting kmers unexplained by model through first peak
  ## truncate errors as soon as it goes to zero, dont allow it to go back up
  error_xcutoff = kcovfloor
  error_xcutoff_ind = tail(which(x<=error_xcutoff),n=1)
  if (length(error_xcutoff_ind)==0) {error_xcutoff_ind=1}

  error_kmers = x[1:error_xcutoff_ind]**(-transform_exp)*(y_transform[1:error_xcutoff_ind] - pred[1:error_xcutoff_ind])

  first_zero = -1

  for (i in 1:error_xcutoff_ind) {
    if (first_zero == -1) {
      if (error_kmers[i] < 1.0) {
        first_zero = i
      }
    }
    else {
      error_kmers[i] = 0
    }
  }

  if (first_zero == -1) {
    first_zero = error_xcutoff_ind
  }

  #if (TRANSFORM) {
  #  y_fit = y_transform
  #} else {
  #  y_fit = y
  #}
  y_fit = y_transform

  ## The fit is residual sum of square error, excluding sequencing errors
  model_fit_all    = c(sum(as.numeric(y_fit[first_zero:length(y_fit)]                          - pred[first_zero:length(y_fit)])                           ** 2), first_zero, x[length(y_fit)])
  model_fit_full   = c(sum(as.numeric(y_fit[first_zero:(min(length(y_fit),(2*p+1)*kcovfloor))] - pred[first_zero:(min(length(y_fit), (2*p+1)*kcovfloor))]) ** 2), first_zero, (min(length(y_fit), (2*p+1)*kcovfloor)))
  model_fit_unique = c(sum(as.numeric(y_fit[first_zero:((p+1)*kcovfloor)]                  - pred[first_zero:((p+1)*kcovfloor)])                   ** 2), first_zero, ((p+1)*kcovfloor))

  ## The score is the percentage of kmers correctly modeled, excluding sequencing errors
  model_fit_allscore    = c(1-sum(abs(as.numeric(y_fit[first_zero:length(y_fit)]                           - pred[first_zero:length(y_fit)])))                           / sum(as.numeric(y_fit[first_zero:length(y_fit)])),                           first_zero, x[length(y_fit)])
  model_fit_fullscore   = c(1-sum(abs(as.numeric(y_fit[first_zero:(min(length(y_fit), (2*p+1)*kcovfloor))] - pred[first_zero:(min(length(y_fit), (2*p+1)*kcovfloor))]))) / sum(as.numeric(y_fit[first_zero:(min(length(y_fit), (2*p+1)*kcovfloor))])), first_zero, (min(length(y_fit), (2*p+1)*kcovfloor)))
  model_fit_uniquescore = c(1-sum(abs(as.numeric(y_fit[first_zero:((p+1)*kcovfloor)]                   - pred[first_zero:((p+1)*kcovfloor)])))                   / sum(as.numeric(y_fit[first_zero:((p+1)*kcovfloor)])),                   first_zero, ((p+1)*kcovfloor))

  fit = data.frame(all  = model_fit_all,      allscore  = model_fit_allscore,
                   full = model_fit_full,     fullscore = model_fit_fullscore,
                   unique = model_fit_unique, uniquescore = model_fit_uniquescore)

  return (fit)
}
