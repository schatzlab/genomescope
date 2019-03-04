#' Uses nlsLM to fit 2p peak model
#'
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @param y A numeric vector of the y-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @param k An integer corresponding to the kmer length.
#' @param p An integer corresponding to the ploidy.
#' @param top An integer corresponding to the topology.
#' @param estKmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param estLength A numeric corresponding to the estimated polyploid genome length.
#' @param max_iterations An integer corresponding to the maximum number iterations to use for nlsLM.
#' @return An nlsLM model object with some additional components.
#' @export
nls_peak<-function(x, y, k, p, top, estKmercov, estLength, max_iterations) {
  #Initiate variables
  model = NULL
  d_min = 0
  d_initial = 0.001
  d_max = 1
  r_min = 0
  r_initial = 0.001
  if (top==0) {
    p_to_num_r = c(0, 1, 2, 4, 6, 10)
  } else {
    p_to_num_r = c(0, 1, 2, 3, 4, 5)
  }
  num_r = p_to_num_r[p]
  if (r_inits!=-1) {
    r_initials = unlist(lapply(strsplit(r_inits,","),as.numeric))
    if (length(r_initials)!=num_r) {
      stop("Incorrect number of initial rates supplied.")
    }
  } else {
    if (ALPHA_RATES) {
      r_initials = c(12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1) #This doesn't work yet.
    } else if (KMER_RATES) {
      r_initials = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
    } else {
      r_initials = rep(r_initial, num_r)
      #r_initials = c(0.010, 0.009, 0.008, 0.007, 0.006, 0.005, 0.004, 0.003, 0.002, 0.001)
    }
  }
  r_start = vector("list", num_r)
  if (p > 1) {
    names(r_start) = paste("r", 1:(num_r), sep="")
    for (i in 1:(num_r)) {
      r_start[[paste("r",i,sep="")]] = r_initials[i]
    }
  }
  r_max = 1
  kmercov_min = 0.5*estKmercov
  kmercov_initial = estKmercov
  kmercov_max = 2*estKmercov
  bias_min = 0
  bias_initial = 0.5
  bias_max = 5
  length_min = 0.5*estLength/p
  length_initial = estLength/p
  length_max = 2*estLength/p

  #Determine what formula to use, based on p
  if (p==1) {
    r_text = ""
  } else {
    r_text = paste(paste(lapply(1:(num_r), function(x) paste("r", as.character(x), sep="")), collapse=", "), ", ",sep="")
    r_param_text = paste(paste(lapply(1:(num_r), function(x) paste("r", as.character(x), "=params[", as.character(x+1), "]",sep="")), collapse="; "), "; ",sep="")
  }
  if (TRANSFORM) {
    #x = head(x,100)
    #y_transform = head(x,100)*head(y,100)
    y_transform = as.numeric(x)*as.numeric(y)
    #formula = as.formula(paste("y_transform ~ x*length*predict",p,"_",top,"(",r_text, "k, d, kmercov, bias, x)",sep=""))
    formula = as.formula(paste("function(params) {d = params[1]; ", r_param_text, "kmercov=params[", as.character(num_r+2), "]; bias=params[", as.character(num_r+3), "]; length = params[", as.character(num_r+4), "]; sum((y_transform - x*length*predict",p,"_",top,"(",r_text, "k, d, kmercov, bias, x))**2)}",sep=""))
  } else {
    #formula = as.formula(paste("y ~ length*predict",p,"_",top,"(",r_text, "k, d, kmercov, bias, x)",sep=""))
    formula = eval(parse(text = paste("function(params) {d = params[1]; ", r_param_text, "kmercov=params[", as.character(num_r+2), "]; bias=params[", as.character(num_r+3), "]; length = params[", as.character(num_r+4), "]; sum((y - length*predict",p,"_",top,"(",r_text, "k, d, kmercov, bias, x))**2)}",sep="")))
  }

  if (VERBOSE) {cat("trying nlsLM algorithm (Levenberg-Marquardt)\n")}

#  try(model <- nlsLM(formula = formula,
#                     start   = c(list(d = d_initial), r_start, list(kmercov = kmercov_initial, bias = bias_initial, length = length_initial)),
#                     lower   = c(c(d_min), rep(r_min, num_r), c(kmercov_min, bias_min, length_min)),
#                     upper   = c(c(d_max), rep(r_max, num_r), c(kmercov_max, bias_max, length_max)),
#                     control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

  model <- DEoptim(fn    = formula,
                   lower = c(c(d_min), rep(r_min, num_r), c(kmercov_min, bias_min, length_min)),
                   upper = c(c(d_max), rep(r_max, num_r), c(kmercov_max, bias_max, length_max)),
                   control = list(itermax = 500))

  print(model)
  stop()

  if (!is.null(model))
  {
    model_sum    = summary(model)
    model$p      = p
    model$top = top
    if (p==1) {
      model$hets = list(c(0, 0))
    } else {
      model$hets = lapply(1:(num_r), function(x) min_max1(model_sum$coefficients[paste('r', x, sep=""),]))
    }
    #model$het = c(1-Reduce("*", 1-unlist(lapply(model$hets, '[[', 1))), 1-Reduce("*", 1-unlist(lapply(model$hets, '[[', 2))))
    model$het = c(sum(sapply(model$hets, '[[', 1)), sum(sapply(model$hets, '[[', 2)))
    model$homo = 1-model$het
    model$dups   = min_max(model_sum$coefficients['bias',])
    model$kcov   = min_max(model_sum$coefficients['kmercov',])
    model$mlen   = min_max(model_sum$coefficients['length',])
    model$md     = min_max1(model_sum$coefficients['d',])
    if (p==1) {
      model$ahets = list(c(0))
    } else {
      model$ahets = lapply(1:(num_r), function(x) model_sum$coefficients[paste('r', x, sep=""),][[1]])
    }
    #model$ahet = 1-Reduce("*", 1-unlist(model$ahets))
    model$ahet = Reduce("+", model$ahets)
    model$ahomo = 1-model$ahet
    model$adups = model_sum$coefficients['bias',][[1]]
    model$akcov = model_sum$coefficients['kmercov',][[1]]
    model$amlen = model_sum$coefficients['length',][[1]]
    model$amd   = model_sum$coefficients['d',][[1]]
  }

  print(model)

  return(model)
}
