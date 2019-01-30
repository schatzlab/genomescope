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
  model = NULL

  if (VERBOSE) {cat("trying nls_peak standard algorithm\n")}

#try p=1 model
  if (p == 1) {
    try(model <- 
      nlsLM(y ~ length*predict1(k, d, kmercov, bias, x),
      start = list(d=0.001, kmercov=estKmercov, bias = 0.5, length = estLength),
      lower = c(0, 0, 0, 0),
      upper = c(1, Inf, Inf, Inf),
      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    if(class(model) == "try-error" || is.null(model)){

      if (VERBOSE) {cat("retrying nls_peak with port algorithm\n")}

      try(model <- 
      nlsLM(y ~ length*predict1(k, d, kmercov, bias, x),
      start = list(d=0, kmercov=estKmercov, bias = 0.5, length=estLength),
      lower = c(0, 0, 0, 0),
      upper = c(1, Inf, Inf, Inf),
      algorithm="port", control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)
    }

    if (!is.null(model))
    {
      model_sum   = summary(model)
      model$p     = 1
      model$het   = c(0,0)
      model$homo  = 1-model$het
      model$dups  = min_max(model_sum$coefficients['bias',])
      model$kcov  = min_max(model_sum$coefficients['kmercov',])
      model$mlen  = min_max(model_sum$coefficients['length',])
      model$md    = min_max1(model_sum$coefficients['d',])
      model$ahet  = 0
      model$ahomo = 1-model$ahet
      model$adups = model_sum$coefficients['bias',][[1]]
      model$akcov = model_sum$coefficients['kmercov',][[1]]
      model$amlen = model_sum$coefficients['length',][[1]]
      model$amd   = model_sum$coefficients['d',][[1]]
    }
  }

#try p=2 model
  if (p == 2) {
    try(model <- 
      nlsLM(y ~ length*predict2(r1, k, d, kmercov, bias, x),
      start = list(d=0.001, r1=0.001, kmercov=estKmercov, bias = 0.5, length = estLength/2),
      lower = c(0, 0, 0, 0, 0),
      upper = c(1, 1, Inf, Inf, Inf),
      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    if(class(model) == "try-error" || is.null(model)){

      if (VERBOSE) {cat("retrying nls_peak with port algorithm\n")}

      try(model <- 
      nlsLM(y ~ length*predict2(r1, k, d, kmercov, bias, x),
      start = list(d=0, r1=0, kmercov=estKmercov, bias = 0.5, length=estLength/2),
      lower = c(0, 0, 0, 0, 0),
      upper = c(1, 1, Inf, Inf, Inf),
      algorithm="port", control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)
    }

    if (!is.null(model))
    {
      model_sum   = summary(model)
      model$p     = 2
      #print('HIIIIIII')
      #print(model_sum$coefficients['r1',])
      model$het1  = min_max1(model_sum$coefficients['r1',])
      #print(model$het1)
      model$het   = model$het1
      model$homo  = 1-model$het
      model$dups  = min_max(model_sum$coefficients['bias',])
      model$kcov  = min_max(model_sum$coefficients['kmercov',])
      model$mlen  = min_max(model_sum$coefficients['length',])
      model$md    = min_max1(model_sum$coefficients['d',])
      model$ahet1 = model_sum$coefficients['r1',][[1]]
      model$ahet  = model$ahet1
      model$ahomo = 1-model$ahet
      model$adups = model_sum$coefficients['bias',][[1]]
      model$akcov = model_sum$coefficients['kmercov',][[1]]
      model$amlen = model_sum$coefficients['length',][[1]]
      model$amd   = model_sum$coefficients['d',][[1]]
    }
  }

#try p=3 model
  if (p==3) {
    try(model <- 
    nlsLM(y ~ length*predict3(r1, r2, k, d, kmercov, bias, x),
    start = list(d=0.001, r1=0.001, r2=0.001, kmercov=estKmercov, bias = 0.5, length = estLength/3),
    lower = c(0, 0, 0, 0, 0, 0),
    upper = c(1, 1, 1, Inf, Inf, Inf),
    control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    if(class(model) == "try-error" || is.null(model)){

      if (VERBOSE) {cat("retrying nls_peak with port algorithm\n")}

      try(model <- 
      nlsLM(y ~ length*predict3(r1, r2, k, d, kmercov, bias, x),
      start = list(d=0, r1=0, r2=0, kmercov=estKmercov, bias = 0.5, length=estLength/3),
      lower = c(0, 0, 0, 0, 0, 0),
      upper = c(1, 1, 1, Inf, Inf, Inf),
      algorithm="port", control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)
    }

    if (!is.null(model))
    {
      model_sum   = summary(model)
      model$p     = 3
      model$het1  = min_max1(model_sum$coefficients['r1',])
      model$het2  = min_max1(model_sum$coefficients['r2',])
      model$het   = model$het1 + model$het2
      model$homo  = 1-model$het
      model$dups  = min_max(model_sum$coefficients['bias',])
      model$kcov  = min_max(model_sum$coefficients['kmercov',])
      model$mlen  = min_max(model_sum$coefficients['length',])
      model$md    = min_max1(model_sum$coefficients['d',])
      model$ahet1 = model_sum$coefficients['r1',][[1]]
      model$ahet2 = model_sum$coefficients['r2',][[1]]
      model$ahet  = model$ahet1 + model$ahet2
      model$ahomo = 1-model$ahet
      model$adups = model_sum$coefficients['bias',][[1]]
      model$akcov = model_sum$coefficients['kmercov',][[1]]
      model$amlen = model_sum$coefficients['length',][[1]]
      model$amd   = model_sum$coefficients['d',][[1]]
      #if (is.null(model) || model_3$m$deviance() *10 < model$m$deviance() || (model_3$m$deviance() < model$m$deviance() && !(model_3$het > 2*model$het) && !(model_3$d > 2*model$d))){
      #  model <- model_3
      #  model$p = 3
      #}
    }
  }

#try p=4 model
  if (p==4) {
    try(model <- 
    nlsLM(y ~ length*predict4(r1, r2, r3, k, d, kmercov, bias, x),
    #enforce r1 <= r2 <= r3
    start = list(d=0.001, r1=0.001, r2 = 0.002, r3=0.003, kmercov=estKmercov, bias = 0.5, length = estLength/4),
    #enforce r1 <= r3 <= r2
    #start = list(d=0.001, r1=0.001, r2 = 0.003, r3=0.002, kmercov=estKmercov, bias = 0.5, length = estLength/4),
    lower = c(0, 0, 0, 0, 0, 0, 0),
    upper = c(1, 1, 1, 1, Inf, Inf, Inf),
    control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    if(class(model) == "try-error" || is.null(model)){

      if (VERBOSE) {cat("retrying nls_peak with port algorithm\n")}

      try(model <- 
      nlsLM(y ~ length*predict4(r1, r2, r3, k, d, kmercov, bias, x),
      #enforce r1 <= r2 <= r3
      start = list(d=0.001, r1=0.001, r2=0.002, r3=0.003, kmercov=estKmercov, bias = 0.5, length=estLength/4),
      #enforce r1 <= r3 <= r2
      #start = list(d=0.001, r1=0.001, r2=0.003, r3=0.002, kmercov=estKmercov, bias = 0.5, length=estLength/4),
      lower = c(0, 0, 0, 0, 0, 0, 0),
      upper = c(1, 1, 1, 1, Inf, Inf, Inf),
      algorithm="port", control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)
	  }

    if (!is.null(model))
    {
      model_sum   = summary(model)
      model$p     = 4
      model$het1  = min_max1(model_sum$coefficients['r1',])
      model$het2  = min_max1(model_sum$coefficients['r2',])
      #model$het2  = c(0,0)
      model$het3  = min_max1(model_sum$coefficients['r3',])
      #model$het4  = min_max1(model_sum$coefficients['r4',])
      #model$het4  = c(0,0)
      model$het   = model$het1 + model$het2 + model$het3 #+ model$het4
      model$homo  = 1-model$het
      model$dups  = min_max(model_sum$coefficients['bias',])
      model$kcov  = min_max(model_sum$coefficients['kmercov',])
      model$mlen  = min_max(model_sum$coefficients['length',])
      model$md    = min_max1(model_sum$coefficients['d',])
      model$ahet1 = model_sum$coefficients['r1',][[1]]
      model$ahet2 = model_sum$coefficients['r2',][[1]]
      #model$ahet2 = 0
      model$ahet3 = model_sum$coefficients['r3',][[1]]
      #model$ahet4 = model_sum$coefficients['r4',][[1]]
      #model$ahet4 = 0
      model$ahet  = model$ahet1 + model$ahet2 + model$ahet3 #+ model$ahet4
      model$ahomo = 1-model$ahet
      model$adups = model_sum$coefficients['bias',][[1]]
      model$akcov = model_sum$coefficients['kmercov',][[1]]
      model$amlen = model_sum$coefficients['length',][[1]]
      model$amd   = model_sum$coefficients['d',][[1]]

      #if (is.null(model) || model_4$m$deviance() *10 < model$m$deviance() || (model_4$m$deviance() < model$m$deviance() && !(model_4$het > 2*model$het) && !(model_4$d > 2*model$d))){
      #  model <- model_4
      #  model$p = 4
      #}
    }
  }

#try p=5 model
  if (p==5) {
    try(model <- 
    nlsLM(y ~ length*predict5(r1, r2, r3, r4, r5, r6, k, d, kmercov, bias, x),
    start = list(d=0.001, r1=0.001, r2=0.001, r3=0.001, r4=0.001, r5=0.001, r6=0.001, kmercov=estKmercov, bias = 0.5, length = estLength/5),
    lower = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    upper = c(1, 1, 1, 1, 1, 1, 1, Inf, Inf, Inf),
    control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    if(class(model) == "try-error" || is.null(model)){

      if (VERBOSE) {cat("retrying nls_peak with port algorithm\n")}

      try(model <- 
      nlsLM(y ~ length*predict5(r1, r2, r3, r4, r5, r6, k, d, kmercov, bias, x),
      start = list(d=0, r1=0, r2=0, r3=0, r4=0, r5=0, r6=0, kmercov=estKmercov, bias = 0.5, length=estLength/5),
      lower = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      upper = c(1, 1, 1, 1, 1, 1, 1, Inf, Inf, Inf),
      algorithm="port", control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)
	  }

    if (!is.null(model))
    {
      model_sum   = summary(model)
      model$p     = 5
      model$het1  = min_max1(model_sum$coefficients['r1',])
      model$het2  = min_max1(model_sum$coefficients['r2',])
      model$het3  = min_max1(model_sum$coefficients['r3',])
      model$het4  = min_max1(model_sum$coefficients['r4',])
      model$het5  = min_max1(model_sum$coefficients['r5',])
      model$het6  = min_max1(model_sum$coefficients['r6',])
      model$het   = model$het1 + model$het2 + model$het3 + model$het4 + model$het5 + model$het6
      model$homo  = 1-model$het
      model$dups  = min_max(model_sum$coefficients['bias',])
      model$kcov  = min_max(model_sum$coefficients['kmercov',])
      model$mlen  = min_max(model_sum$coefficients['length',])
      model$md    = min_max1(model_sum$coefficients['d',])
      model$ahet1 = model_sum$coefficients['r1',][[1]]
      model$ahet2 = model_sum$coefficients['r2',][[1]]
      model$ahet3 = model_sum$coefficients['r3',][[1]]
      model$ahet4 = model_sum$coefficients['r4',][[1]]
      model$ahet5 = model_sum$coefficients['r5',][[1]]
      model$ahet6 = model_sum$coefficients['r6',][[1]]
      model$ahet  = model$ahet1 + model$ahet2 + model$ahet3 + model$ahet4 + model$ahet5 + model$ahet6
      model$ahomo = 1-model$ahet
      model$adups = model_sum$coefficients['bias',][[1]]
      model$akcov = model_sum$coefficients['kmercov',][[1]]
      model$amlen = model_sum$coefficients['length',][[1]]
      model$amd   = model_sum$coefficients['d',][[1]]
      #if (is.null(model) || model_5$m$deviance() *10 < model$m$deviance() || (model_5$m$deviance() < model$m$deviance() && !(model_5$het > 2*model$het) && !(model_5$d > 2*model$d))){
      #  model <- model_5
      #  model$p = 5
      #}
    }
  }

#try p=6 model
  if (p==6) {
    try(model <- 
    nlsLM(y ~ length*predict6(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, k, d, kmercov, bias, x),
    start = list(d=0.001, r1=0.001, r2=0.001, r3=0.001, r4=0.001, r5=0.001, r6=0.001, r7=0.001, r8=0.001, r9=0.001, r10=0.001, kmercov=estKmercov, bias = 0.5, length = estLength/6),
    lower = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    upper = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, Inf, Inf, Inf),
    control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    if(class(model) == "try-error" || is.null(model)){

      if (VERBOSE) {cat("retrying nls_peak with port algorithm\n")}

      try(model <- 
      nlsLM(y ~ length*predict6(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, k, d, kmercov, bias, x),
      start = list(d=0, r1=0, r2=0, r3=0, r4=0, r5=0, r6=0, r7=0, r8=0, r9=0, r10=0, kmercov=estKmercov, bias = 0.5, length=estLength/6),
      lower = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      upper = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, Inf, Inf, Inf),
      algorithm="port", control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)
	  }

    if (!is.null(model))
    {
      model_sum    = summary(model)
      model$p      = 6
      model$het1   = min_max1(model_sum$coefficients['r1',])
      model$het2   = min_max1(model_sum$coefficients['r2',])
      model$het3   = min_max1(model_sum$coefficients['r3',])
      model$het4   = min_max1(model_sum$coefficients['r4',])
      model$het5   = min_max1(model_sum$coefficients['r5',])
      model$het6   = min_max1(model_sum$coefficients['r6',])
      model$het7   = min_max1(model_sum$coefficients['r7',])
      model$het8   = min_max1(model_sum$coefficients['r8',])
      model$het9   = min_max1(model_sum$coefficients['r9',])
      model$het10  = min_max1(model_sum$coefficients['r10',])
      model$het    = model$het1 + model$het2 + model$het3 + model$het4 + model$het5 + model$het6 + model$het7 + model$het8 + model$het9 + model$het10
      model$homo   = 1-model$het
      model$dups   = min_max(model_sum$coefficients['bias',])
      model$kcov   = min_max(model_sum$coefficients['kmercov',])
      model$mlen   = min_max(model_sum$coefficients['length',])
      model$md     = min_max1(model_sum$coefficients['d',])
      model$ahet1  = model_sum$coefficients['r1',][[1]]
      model$ahet2  = model_sum$coefficients['r2',][[1]]
      model$ahet3  = model_sum$coefficients['r3',][[1]]
      model$ahet4  = model_sum$coefficients['r4',][[1]]
      model$ahet5  = model_sum$coefficients['r5',][[1]]
      model$ahet6  = model_sum$coefficients['r6',][[1]]
      model$ahet7  = model_sum$coefficients['r7',][[1]]
      model$ahet8  = model_sum$coefficients['r8',][[1]]
      model$ahet9  = model_sum$coefficients['r9',][[1]]
      model$ahet10 = model_sum$coefficients['r10',][[1]]
      model$ahet   = model$ahet1 + model$ahet2 + model$ahet3 + model$ahet4 + model$ahet5 + model$ahet6 + model$ahet7 + model$ahet8 + model$ahet9 + model$ahet10
      model$ahomo = 1-model$ahet
      model$adups = model_sum$coefficients['bias',][[1]]
      model$akcov = model_sum$coefficients['kmercov',][[1]]
      model$amlen = model_sum$coefficients['length',][[1]]
      model$amd   = model_sum$coefficients['d',][[1]]
      #if (is.null(model) || model_6$m$deviance() *10 < model$m$deviance() || (model_6$m$deviance() < model$m$deviance() && !(model_6$het > 2*model$het) && !(model_6$d > 2*model$d))){
      #  model <- model_6
      #  model$p = 6
      #}
    }
  }

  print(model)

  return(model)
}
