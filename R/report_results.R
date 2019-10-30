## Format numbers
###############################################################################

bp_format<-function(num) {paste(formatC(round(num),format="f",big.mark=",", digits=0), "bp",sep=" ")}

percentage_format<-function(num) {paste(signif(num,6)*100,"%",sep="")}

X_format<-function(num) {paste(signif(num,4),"X",sep="")}

#' Report results and make plots
#'
#' @param kmer_hist A data frame of the original histogram data (starting at 1 and going up to the max kmer coverage threshold).
#' @param kmer_hist_orig A data frame of the original histogram data (starting at 1 and with last position removed).
#' @param k An integer corresponding to the kmer length.
#' @param p An integer corresponding to the ploidy.
#' @param container A list (nls, nlsscore) where nls is the nlsLM model object (with some additional components)
#' and nlsscore is the score (model RSSE) corresponding to the best fit (of the p forms).
#' @param foldername A character vector corresponding to the name of the output directory.
#' @param arguments A data frame of the user-specified inputs.
#' @param IN_VERBOSE A boolean flag to designate whether report_results is being called in a VERBOSE block.
#' @export
report_results<-function(kmer_hist,kmer_hist_orig, k, p, container, foldername, arguments, IN_VERBOSE) {
  
  x=kmer_hist_orig[[1]]
  y_orig=kmer_hist_orig[[2]]
  y = as.numeric(x)**transform_exp*as.numeric(y_orig)
  kmer_hist_transform = kmer_hist_orig
  kmer_hist_transform$V2 = as.numeric(kmer_hist_transform$V1)**transform_exp * as.numeric(kmer_hist_transform$V2)
  model = container[[1]]

  #automatically zoom into the relevant regions of the plot, ignore first 15 positions
  xmax=length(x)
  start_orig=which(y_orig == min(y_orig[1:TYPICAL_ERROR]))
  start=which(y == min(y[1:TYPICAL_ERROR]))
  zoomx=x[start:(xmax-1)]
  zoomy_orig=y_orig[start_orig:(xmax-1)]
  zoomy=y[start:(xmax-1)]

  ## allow for a little space above max value past the noise
  y_limit_orig = max(zoomy_orig[start_orig:length(zoomy_orig)])*1.1
  y_limit = max(zoomy[start:length(zoomy)])*1.1

  x_limit_orig = which(y_orig == max(y_orig[start_orig:length(zoomx)])) * 3
  x_limit = which(y == max(y[start:length(zoomx)])) * 3

  if (min(zoomy_orig) > zoomy_orig[1]){
    x_limit_orig=max(which(zoomy_orig<zoomy_orig[1])[2],600)
  }
  if (min(zoomy) > zoomy[1]){
    x_limit=max(which(zoomy<zoomy[1])[2],600)
  }

  if (!is.null(model))
  {
    model_sum=summary(model)
    kcov = min_max(model_sum$coefficients['kmercov',])[1]
    x_limit_orig = max(kcov*(2*p+1.1), x_limit_orig)
    x_limit = max(kcov*(2*p+1.1), x_limit)
    if (model$top==0) {
      p_to_num_r = c(0, 1, 2, 4, 6, 10)
    } else {
      p_to_num_r = c(0, 1, 2, 3, 4, 5)
    }
  } else {
    if (topology==0) {
      p_to_num_r = c(0, 1, 2, 4, 6, 10)
    } else {
      p_to_num_r = c(0, 1, 2, 3, 4, 5)
    }
  }

  ## Uncomment this to enforce a specific number
  # x_limit=150

  ## Features to report
  het=c(-1,-1)
  homo=c(-1,-1)
  num_r = p_to_num_r[p]
  if (p > 1) {
    hets = lapply(1:(num_r), function(x) c(-1, -1))
    ahets = lapply(1:(num_r), function(x) -1)
  }
  amd = -1
  akcov = -1
  adups = -1
  amlen = -1
  atotal_len = -1
  top = -1
  total_len=c(-1,-1)
  repeat_len=c(-1,-1)
  unique_len=c(-1,-1)
  dups=c(-1,-1)
  error_rate=c(-1,-1)
  model_status="fail"

  model_fit_unique      = c(0,0,0)
  model_fit_full        = c(0,0,0)
  model_fit_all         = c(0,0,0)
  model_fit_allscore    = c(0,0,0)
  model_fit_fullscore   = c(0,0,0)
  model_fit_uniquescore = c(0,0,0)

  plot_size=2000
  font_size=1.2
  resolution=300

  ## Plot the distribution, and hopefully with the model fit
  ylabel_orig = "Frequency"
  if (transform_exp == 1) {
    ylabel_transform = "Coverage*Frequency"
  } else {
    ylabel_transform = paste("Coverage^", transform_exp, "*Frequency", sep="")
  }
  png(paste(foldername, "/", arguments$name_prefix, "linear_plot.png", sep=""),
  width=plot_size, height=plot_size, res=resolution)
  par(mar = c(5.1,4.1,6.1,2.1))
  plot(kmer_hist_orig, type="n", main="GenomeScope Profile\n\n\n",
  xlab="Coverage", ylab=ylabel_orig, ylim=c(0,y_limit_orig), xlim=c(0,x_limit_orig),
  cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
  #rect(0, 0, max(kmer_hist_orig[[1]])*1.1 , max(kmer_hist_orig[[2]])*1.1, col=COLOR_BGCOLOR)
  rect(0, 0, x_limit_orig*1.1 , y_limit_orig*1.1, col=COLOR_BGCOLOR)
  points(kmer_hist_orig, type="h", col=COLOR_HIST, lwd=2)
#  if(length(kmer_hist[,1])!=length(kmer_hist_orig[,1])){
#    abline(v=length(kmer_hist[,1]),col=COLOR_COVTHRES,lty="dashed", lwd=3)
#  }
  box(col="black")

  png(paste(foldername, "/", arguments$name_prefix, "transformed_linear_plot.png", sep=""),
  width=plot_size, height=plot_size, res=resolution)
  par(mar = c(5.1,4.1,6.1,2.1))
  plot(kmer_hist_transform, type="n", main="GenomeScope Profile\n\n\n",
  xlab="Coverage", ylab=ylabel_transform, ylim=c(0,y_limit), xlim=c(0,x_limit),
  cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
  #rect(0, 0, max(kmer_hist_orig[[1]])*1.1 , max(kmer_hist_orig[[2]])*1.1, col=COLOR_BGCOLOR)
  rect(0, 0, x_limit*1.1 , y_limit*1.1, col=COLOR_BGCOLOR)
  points(kmer_hist_transform, type="h", col=COLOR_HIST, lwd=2)
#  if(length(kmer_hist[,1])!=length(kmer_hist_orig[,1])){
#    abline(v=length(kmer_hist[,1]),col=COLOR_COVTHRES,lty="dashed", lwd=3)
#  }
  box(col="black")

  ## Make a second plot in log space over entire range
  png(paste(foldername, "/", arguments$name_prefix, "log_plot.png", sep=""),
  width=plot_size, height=plot_size, res=resolution)
  par(mar = c(5.1,4.1,6.1,2.1))
  plot(kmer_hist_orig, type="n", main="GenomeScope Profile\n\n\n",
  xlab="Coverage", ylab=ylabel_orig, log="xy",
  cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
  rect(1e-10, 1e-10, max(kmer_hist_orig[,1])*10 , max(kmer_hist_orig[,2])*10, col=COLOR_BGCOLOR)
  points(kmer_hist_orig, type="h", col=COLOR_HIST, lwd=2)
  if(length(kmer_hist[,1])!=length(kmer_hist_orig[,1])){
    abline(v=length(kmer_hist[,1]),col=COLOR_COVTHRES,lty="dashed", lwd=3)
  }
  box(col="black")

  png(paste(foldername, "/", arguments$name_prefix, "transformed_log_plot.png", sep=""),
  width=plot_size, height=plot_size, res=resolution)
  par(mar = c(5.1,4.1,6.1,2.1))
  plot(kmer_hist_transform, type="n", main="GenomeScope Profile\n\n\n",
  xlab="Coverage", ylab=ylabel_transform, log="xy",
  cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
  rect(1e-10, 1e-10, max(kmer_hist_transform[,1])*10 , max(kmer_hist_transform[,2])*10, col=COLOR_BGCOLOR)
  points(kmer_hist_transform, type="h", col=COLOR_HIST, lwd=2)
  if(length(kmer_hist[,1])!=length(kmer_hist_transform[,1])){
    abline(v=length(kmer_hist[,1]),col=COLOR_COVTHRES,lty="dashed", lwd=3)
  }
  box(col="black")


  if(!is.null(model))
  {
    x=kmer_hist[[1]]
    y=kmer_hist[[2]]
    y_transform = as.numeric(x)**transform_exp*as.numeric(y)

    ## The model converged!
    pred=predict(model, newdata=data.frame(x))

    ## Compute the genome characteristics
    model_sum=summary(model)
    #print(model_sum)

    ## save the model to a file
    capture.output(model_sum, file=paste(foldername,"/", arguments$name_prefix, "model.txt", sep=""))

    ## Identify key values
    top   = model$top
    hets  = model$hets
    het   = model$het
    ahets = model$ahets
    ahet  = model$ahet
    homo  = model$homo
    ahomo = model$ahomo

    dups = model$dups
    kcov = model$kcov
    mlen = model$mlen
    md   = model$md

    adups = model$adups
    akcov = model$akcov
    amlen = model$amlen
    amd   = model$amd

    ## Compute error rate, by counting kmers unexplained by model through first peak
    ## truncate errors as soon as it goes to zero, dont allow it to go back up
    error_xcutoff = max(1, floor(kcov[1]))
    error_xcutoff_ind = tail(which(x<=error_xcutoff),n=1)
    if (length(error_xcutoff_ind)==0) {error_xcutoff_ind=1}

    error_kmers = x[1:error_xcutoff_ind]**(-transform_exp)*(y_transform[1:error_xcutoff_ind] - pred[1:error_xcutoff_ind])

    first_zero = -1

    for (i in 1:error_xcutoff_ind)
    {
      if (first_zero == -1)
      {
        if (error_kmers[i] < 1.0)
        {
          first_zero = i
          if (VERBOSE) {cat(paste("Truncating errors at", i, "\n"))}
        }
      }
      else
      {
        error_kmers[i] = 0
      }
    }

    if (first_zero == -1)
    {
      first_zero = error_xcutoff_ind
      if (VERBOSE) {cat(paste("Truncating errors at", error_xcutoff_ind, "\n"))}
    }

    ## Rather than "0", set to be some very small number so log-log plot looks okay
    error_kmers = pmax(error_kmers, 1e-10)

    total_error_kmers = sum(as.numeric(error_kmers) * as.numeric(x[1:error_xcutoff_ind]))

    total_kmers = sum(as.numeric(x)*as.numeric(y))

    error_rate = 1-(1-(total_error_kmers/total_kmers))**(1/k)
    error_rate = c(error_rate, error_rate)

    total_len = (total_kmers-total_error_kmers)/(p*kcov)
    atotal_len = (total_kmers-total_error_kmers)/(p*akcov)

    ## find kmers that fit the p peak model (no repeats)
    if (p==1)
    {
      unique_hist = amlen*predict1_1_unique(k, amd, akcov, adups, x)
    }
    if (p==2)
    {
      unique_hist = amlen*predict2_1_unique(ahets[[1]], k, amd, akcov, adups, x)
    }
    if (p==3)
    {
      unique_hist = amlen*predict3_1_unique(ahets[[1]], ahets[[2]], k, amd, akcov, adups, x)
    }
    if (p==4)
    {
      if (top==0) {
        unique_hist = amlen*predict4_0_unique(ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], k, amd, akcov, adups, x)
      } else {
        unique_hist = eval(parse(text = paste("amlen*predict4_", top, "_unique(ahets[[1]], ahets[[2]], ahets[[3]], k, amd, akcov, adups, x)", sep="")))
      }
    }
    if (p==5)
    {
      if (top==0) {
        unique_hist = amlen*predict5_0_unique(ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], ahets[[6]], k, amd, akcov, adups, x)
      } else {
        unique_hist = eval(parse(text = paste("amlen*predict5_", top, "_unique(ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], k, amd, akcov, adups, x)", sep="")))
      }
    }
    if (p==6)
    {
      if (top==0) {
        unique_hist = amlen*predict6_0_unique(ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], ahets[[6]], ahets[[7]], ahets[[8]], ahets[[9]], ahets[[10]], k, amd, akcov, adups, x)
      } else {
        unique_hist = eval(parse(text = paste("amlen*predict6_", top, "_unique(ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], k, amd, akcov, adups, x)", sep="")))
      }
    }

    unique_hist_transform = x**transform_exp*unique_hist

    unique_kmers = sum(as.numeric(x)*as.numeric(unique_hist))
    repeat_kmers = max(0, total_kmers - unique_kmers - total_error_kmers)

    repeat_len=repeat_kmers/(p*kcov)
    if (repeat_kmers == 0) {
      unique_len = total_len
    } else {
      unique_len=unique_kmers/(p*kcov)
    }

    score = container[[2]]

    model_fit_allscore    = score$allscore
    model_fit_fullscore   = score$fullscore
    model_fit_uniquescore = score$uniquescore

    model_fit_all    = score$all
    model_fit_full   = score$full
    model_fit_unique = score$unique

    residual_transform = y_transform - pred
    residual = x**(-transform_exp)*residual_transform

    hetline_simple = paste0("heterozygosity: ", format(100*ahet, digits=3), "%")

    if (p==1) {
      hetline = paste0("a:", format(100*ahomo, digits=3), "%")
    }
    if (p==2) {
      hetline = paste0("aa:", format(100*ahomo,      digits=3), "% ",
                       "ab:", format(100*ahets[[1]], digits=3), "%")
    }
    if (p==3) {
      hetline = paste0("aaa:", format(100*ahomo,      digits=3), "% ",
                       "aab:", format(100*ahets[[1]], digits=3), "% ",
                       "abc:", format(100*ahets[[2]], digits=3), "%")
    }
    if (p==4) {
      if (top==0) {
        hetline = paste0("aaaa:", format(100*ahomo,      digits=3), "% ",
                         "aaab:", format(100*ahets[[1]], digits=3), "% ",
                         "aabb:", format(100*ahets[[2]], digits=3), "% ",
                         "aabc:", format(100*ahets[[3]], digits=3), "% ",
                         "abcd:", format(100*ahets[[4]], digits=3), "%")
      } else {
        hetline = paste0("aaaa:",                       format(100*ahomo,      digits=3), "% ",
                         switch(top, "aaab:", "aabb:"), format(100*ahets[[1]], digits=3), "% ",
                         "aabc:",                       format(100*ahets[[2]], digits=3), "% ",
                         "abcd:",                       format(100*ahets[[3]], digits=3), "%")
      }
    }
    if (p==5) {
      if (top==0) {
        hetline = paste0("aaaaa:", format(100*ahomo,      digits=3), "% ",
                         "aaaab:", format(100*ahets[[1]], digits=3), "% ",
                         "aaabb:", format(100*ahets[[2]], digits=3), "% ",
                         "aaabc:", format(100*ahets[[3]], digits=3), "% ",'\n',
                         "aabbc:", format(100*ahets[[4]], digits=3), "% ",
                         "aabcd:", format(100*ahets[[5]], digits=3), "% ",
                         "abcde:", format(100*ahets[[6]], digits=3), "%")
      } else {
        hetline = paste0("aaaaa:",                                                      format(100*ahomo,      digits=3), "% ",
                         switch(top, "aaaab:", "aaaab:", "aaabb:", "aaabb:", "aaabb:"), format(100*ahets[[1]], digits=3), "% ",
                         switch(top, "aaabc:", "aabbc:", "aaabc:", "aabcc:", "aabcc:"), format(100*ahets[[2]], digits=3), "% ",'\n',
                         switch(top, "aabcd:", "aabcd:", "aabcd:", "aabcd:", "abcdd:"), format(100*ahets[[3]], digits=3), "% ",
                         "abcde:",                                                      format(100*ahets[[4]], digits=3), "%")
      }
    }
    if (p==6) {
      if (top==0) {
        hetline = paste0("aaaaaa:", format(100*ahomo, digits=3), "% ",
                         "aaaaab:", format(100*ahets[[1]], digits=3), "% ",
                         "aaaabb:", format(100*ahets[[2]], digits=3), "% ",
                         "aaabbb:", format(100*ahets[[3]], digits=3), "% ",'\n',
                         "aaaabc:", format(100*ahets[[4]], digits=3), "% ",
                         "aaabbc:", format(100*ahets[[5]], digits=3), "% ",
                         "aabbcc:", format(100*ahets[[6]], digits=3), "% ",
                         "aaabcd:", format(100*ahets[[7]], digits=3), "% ",'\n',
                         "aabbcd:", format(100*ahets[[8]], digits=3), "% ",
                         "aabcde:", format(100*ahets[[9]], digits=3), "% ",
                         "abcdef:", format(100*ahets[[10]], digits=3), "%")
      } else {
        hetline = paste0("aaaaaa:", format(100*ahomo, digits=3), "% ",
                         switch(top, "aaaaab:", "aaaaab:", "aaaaab:", "aaaaab:", "aaaaab:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaabbb:", "aaabbb:", "aaabbb:"), format(100*ahets[[1]], digits=3), "% ",
                         switch(top, "aaaabc:", "aaaabc:", "aaabbc:", "aaabbc:", "aaabbc:", "aaaabc:", "aaaabc:", "aaabcc:", "aaabcc:", "aaabcc:", "aabbcc:", "aabbcc:", "aabbcc:", "aaabbc:", "aaabbc:", "aaabbc:"), format(100*ahets[[2]], digits=3), "% ",'\n',
                         switch(top, "aaabcd:", "aabbcd:", "aaabcd:", "aabccd:", "aabccd:", "aaabcd:", "aabbcd:", "aaabcd:", "aabcdd:", "aabcdd:", "aabbcd:", "aabcdd:", "aabcdd:", "aaabcd:", "aabccd:", "aabccd:"), format(100*ahets[[3]], digits=3), "% ",
                         switch(top, "aabcde:", "aabcde:", "aabcde:", "aabcde:", "abcdde:", "aabcde:", "aabcde:", "aabcde:", "aabcde:", "abcdee:", "aabcde:", "aabcde:", "abcdee:", "aabcde:", "aabcde:", "abcdde:"), format(100*ahets[[4]], digits=3), "% ",
                         "abcdef:", format(100*ahets[[5]], digits=3), "%")
      }
    }

    if (p >= 5) {
      hetline = hetline_simple
    }

    if (!IN_VERBOSE) {
      cat(paste0(hetline,"\n"))
    }

    dev.set(dev.next())

    ## Finish Linear Plot
    title(paste("\n\nlen:",  prettyNum(total_len[1], big.mark=","),
    "bp",
    " uniq:", format(100*(unique_len[1]/total_len[1]), digits=3),
    "% ", "\n",
    hetline, "\n",
    " kcov:", format(akcov, digits=3),
    " err:",   format(100*error_rate[1], digits=3),
    "% ",
    " dup:",  format(adups, digits=3),
    " ",
    " k:",   format(k, digits=3),
    " p:",   format(p, digits=3),
    sep=""),
    cex.main=.85)

    ## Mark the modes of the peaks
    abline(v=akcov * (1:(2*p)), col=COLOR_KMERPEAK, lty=2)

    ## Draw just the unique portion of the model
    if (!NO_UNIQUE_SEQUENCE) {
      lines(x, unique_hist, col=COLOR_pPEAK, lty=1, lwd=3)
    }
    lines(x, x**(-transform_exp)*pred, col=COLOR_2pPEAK, lwd=3)
    lines(x[1:error_xcutoff_ind], error_kmers, lwd=3, col=COLOR_ERRORS)

    if (VERBOSE) {
      lines(x, residual, col=COLOR_RESIDUAL, lwd=3)
    }

    ## Add legend
    if (NO_UNIQUE_SEQUENCE) {
      legend(.62 * x_limit_orig, 1.0 * y_limit_orig,
      legend=c("observed", "full model", "errors", "kmer-peaks"),
      lty=c("solid", "solid", "solid", "dashed"),
      lwd=c(3,3,3,2),
      col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
      bg="white")
    } else {
      legend(.62 * x_limit_orig, 1.0 * y_limit_orig,
      legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks"),
      lty=c("solid", "solid", "solid", "solid", "dashed"),
      lwd=c(3,3,3,3,2),
      col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
      bg="white")
    }

    dev.set(dev.next())

    ## Finish Linear Plot
    title(paste("\n\nlen:",  prettyNum(total_len[1], big.mark=","),
    "bp",
    " uniq:", format(100*(unique_len[1]/total_len[1]), digits=3),
    "% ", "\n",
    hetline, "\n",
    " kcov:", format(akcov, digits=3),
    " err:",   format(100*error_rate[1], digits=3),
    "% ",
    " dup:",  format(adups, digits=3),
    " ",
    " k:",   format(k, digits=3),
    " p:",   format(p, digits=3),
    sep=""),
    cex.main=.85)

    ## Mark the modes of the peaks
    abline(v=akcov * (1:(2*p)), col=COLOR_KMERPEAK, lty=2)

    ## Draw just the unique portion of the model
    if (!NO_UNIQUE_SEQUENCE) {
      lines(x, unique_hist_transform, col=COLOR_pPEAK, lty=1, lwd=3)
    }
    lines(x, pred, col=COLOR_2pPEAK, lwd=3)
    lines(x[1:error_xcutoff_ind], (x[1:error_xcutoff_ind]**transform_exp)*error_kmers, lwd=3, col=COLOR_ERRORS)

    if (VERBOSE) {
      lines(x, residual_transform, col=COLOR_RESIDUAL, lwd=3)
    }

    ## Add legend
    if (NO_UNIQUE_SEQUENCE) {
      legend(.62 * x_limit, 1.0 * y_limit,
      legend=c("observed", "full model", "errors", "kmer-peaks"),
      lty=c("solid", "solid", "solid", "dashed"),
      lwd=c(3,3,3,2),
      col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
      bg="white")
    } else {
      legend(.62 * x_limit, 1.0 * y_limit,
      legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks"),
      lty=c("solid", "solid", "solid", "solid", "dashed"),
      lwd=c(3,3,3,3,2),
      col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
      bg="white")
    }

    dev.set(dev.next())

    ## Finish Log plot
    title(paste("\n\nlen:",  prettyNum(total_len[1], big.mark=","),
    "bp",
    " uniq:", format(100*(unique_len[1]/total_len[1]), digits=3),
    "% ", "\n",
    hetline, "\n",
    " kcov:", format(akcov, digits=3),
    " err:",   format(100*error_rate[1], digits=3),
    "% ",
    " dup:",  format(adups, digits=3),
    " ",
    " k:",   format(k, digits=3),
    " p:",   format(p, digits=3),
    sep=""),
    cex.main=.85)

    ## Mark the modes of the peaks
    abline(v=akcov * (1:(2*p)), col=COLOR_KMERPEAK, lty=2)

    ## Draw just the unique portion of the model
    if (!NO_UNIQUE_SEQUENCE) {
      lines(x, unique_hist, col=COLOR_pPEAK, lty=1, lwd=3)
    }
    lines(x, x**(-transform_exp)*pred, col=COLOR_2pPEAK, lwd=3)
    lines(x[1:error_xcutoff_ind], error_kmers, lwd=3, col=COLOR_ERRORS)

    if (VERBOSE) {
      lines(x, residual, col=COLOR_RESIDUAL, lwd=3)
    }

    ## Add legend
    if(length(kmer_hist[,1])==length(kmer_hist_orig[,1]))
    {
      if (NO_UNIQUE_SEQUENCE) {
        legend(exp(.62 * log(max(x))), 1.0 * max(y),
        legend=c("observed", "full model", "errors", "kmer-peaks"),
        lty=c("solid", "solid", "solid", "dashed"),
        lwd=c(3,3,3,3),
        col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
        bg="white")
      } else {
        legend(exp(.62 * log(max(x))), 1.0 * max(y),
        legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks"),
        lty=c("solid", "solid", "solid", "solid", "dashed"),
        lwd=c(3,3,3,3,3),
        col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
        bg="white")
      }
    }
    else
    {
      if (NO_UNIQUE_SEQUENCE) {
        legend("topright",
        ##legend(exp(.62 * log(max(x))), 1.0 * max(y),
        legend=c("observed", "full model", "errors", "kmer-peaks","cov-threshold"),
        lty=c("solid", "solid", "solid", "dashed", "dashed"),
        lwd=c(3,3,3,2,3),
        col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_ERRORS, COLOR_KMERPEAK, COLOR_COVTHRES),
        bg="white")
      } else {
        legend("topright",
        ##legend(exp(.62 * log(max(x))), 1.0 * max(y),
        legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks","cov-threshold"),
        lty=c("solid", "solid", "solid", "solid", "dashed", "dashed"),
        lwd=c(3,3,3,3,2,3),
        col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_pPEAK, COLOR_ERRORS, COLOR_KMERPEAK, COLOR_COVTHRES),
        bg="white")
      }
    }

    dev.set(dev.next())

    ## Finish Log plot
    title(paste("\n\nlen:",  prettyNum(total_len[1], big.mark=","),
    "bp",
    " uniq:", format(100*(unique_len[1]/total_len[1]), digits=3),
    "% ", "\n",
    hetline, "\n",
    " kcov:", format(akcov, digits=3),
    " err:",   format(100*error_rate[1], digits=3),
    "% ",
    " dup:",  format(adups, digits=3),
    " ",
    " k:",   format(k, digits=3),
    " p:",   format(p, digits=3),
    sep=""),
    cex.main=.85)

    ## Mark the modes of the peaks
    abline(v=akcov * (1:(2*p)), col=COLOR_KMERPEAK, lty=2)

    ## Draw just the unique portion of the model
    if (!NO_UNIQUE_SEQUENCE) {
      lines(x, unique_hist_transform, col=COLOR_pPEAK, lty=1, lwd=3)
    }
    lines(x, pred, col=COLOR_2pPEAK, lwd=3)
    lines(x[1:error_xcutoff_ind], (x[1:error_xcutoff_ind]**transform_exp)*error_kmers, lwd=3, col=COLOR_ERRORS)

    if (VERBOSE) {
      lines(x, residual_transform, col=COLOR_RESIDUAL, lwd=3)
    }

    ## Add legend
    if(length(kmer_hist[,1])==length(kmer_hist_orig[,1]))
    {
      if (NO_UNIQUE_SEQUENCE) {
        legend(exp(.62 * log(max(x))), 1.0 * max(y),
        legend=c("observed", "full model", "errors", "kmer-peaks"),
        lty=c("solid", "solid", "solid", "dashed"),
        lwd=c(3,3,3,3),
        col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
        bg="white")
      } else {
        legend(exp(.62 * log(max(x))), 1.0 * max(y),
        legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks"),
        lty=c("solid", "solid", "solid", "solid", "dashed"),
        lwd=c(3,3,3,3,3),
        col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
        bg="white")
      }
    }
    else
    {
      if (NO_UNIQUE_SEQUENCE) {
        legend("topright",
        ##legend(exp(.62 * log(max(x))), 1.0 * max(y),
        legend=c("observed", "full model", "errors", "kmer-peaks","cov-threshold"),
        lty=c("solid", "solid", "solid", "dashed", "dashed"),
        lwd=c(3,3,3,2,3),
        col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_ERRORS, COLOR_KMERPEAK, COLOR_COVTHRES),
        bg="white")
      } else {
        legend("topright",
        ##legend(exp(.62 * log(max(x))), 1.0 * max(y),
        legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks","cov-threshold"),
        lty=c("solid", "solid", "solid", "solid", "dashed", "dashed"),
        lwd=c(3,3,3,3,2,3),
        col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_pPEAK, COLOR_ERRORS, COLOR_KMERPEAK, COLOR_COVTHRES),
        bg="white")
      }
    }

    model_status="done"

    if (!IN_VERBOSE) {
      cat(paste("Model converged het:", format(ahet, digits=3),
      " kcov:", format(akcov, digits=3),
      " err:", format(error_rate[1], digits=3),
      " model fit:", format(adups, digits=3),
      " len:", round(total_len[1]), "\n", sep=""))
    }
  }
  else
  {
    title("\nFailed to converge")
    dev.set(dev.next())
    title("\nFailed to converge")
    cat("Failed to converge.", file=paste(foldername,"/", arguments$name_prefix, "model.txt", sep=""))
    cat("Failed to converge.\n")
  }

  dev.off()
  dev.off()
  dev.off()
  dev.off()

  ## Write key values to summary file
  summaryFile <- paste(foldername,"/", arguments$name_prefix, "summary.txt",sep="")

  format_column_1 = "%-30s"
  format_column_2 = "%-18s"
  format_column_3 = "%-18s"

  cat(paste("GenomeScope version 2.0", sep=""), file=summaryFile, sep="\n")
  cat(paste("input file = ", arguments$input, sep=""), file=summaryFile, sep="\n", append=TRUE)
  cat(paste("output directory = ", arguments$output, sep=""), file=summaryFile, sep="\n", append=TRUE)
  cat(paste("p = ", p,sep=""), file=summaryFile, sep="\n", append=TRUE)
  cat(paste("k = ", k,sep=""), file=summaryFile, sep="\n", append=TRUE)
  if (arguments$name_prefix!="") {
    cat(paste("name prefix = ", substring(arguments$name_prefix,1,nchar(arguments$name_prefix)-1), sep=""), file=summaryFile, sep="\n", append=TRUE)
  }
  if (arguments$lambda!=-1) {
    cat(paste("initial kmercov estimate = ", arguments$lambda, sep=""), file=summaryFile, sep="\n", append=TRUE)
  }
  if (arguments$max_kmercov!=-1) {
    cat(paste("max_kmercov = ", arguments$max_kmercov, sep=""), file=summaryFile, sep="\n", append=TRUE)
  }
  if (VERBOSE) {
    cat(paste("VERBOSE set to TRUE", sep=""), file=summaryFile, sep="\n", append=TRUE)
  }
  if (NO_UNIQUE_SEQUENCE) {
    cat(paste("NO_UNIQUE_SEQUENCE set to TRUE", sep=""), file=summaryFile, sep="\n", append=TRUE)
  }
  if (topology!=0) {
    cat(paste("topology = ", topology, sep=""), file=summaryFile, sep="\n", append=TRUE)
  }
  if (d_init!=-1) {
    cat(paste("initial repetitiveness = ", d_init, sep=""), file=summaryFile, sep="\n", append=TRUE)
  }
  if (r_inits!=-1) {
    cat(paste("initial heterozygosities = ", r_inits, sep=""), file=summaryFile, sep="\n", append=TRUE)
  }
  if (transform_exp != 1) {
    cat(paste("TRANSFORM_EXP = ", transform_exp, sep=""), file=summaryFile, sep="\n", append=TRUE)
  }
  if (TESTING) {
    cat(paste("TESTING set to TRUE", sep=""), file=summaryFile, sep="\n", append=TRUE)
  }
  if (TRUE_PARAMS != -1) {
    cat(paste("TRUE_PARAMS = ", TRUE_PARAMS, sep=""), file=summaryFile, sep="\n", append=TRUE)
  }
  if (TRACE_FLAG) {
    cat(paste("TRACE_FLAG set to TRUE", sep=""), file=summaryFile, sep="\n", append=TRUE)
  }
  if (NUM_ROUNDS != 4) {
    cat(paste("NUM_ROUNDS = ", NUM_ROUNDS, sep=""), file=summaryFile, sep="\n", append=TRUE)
  }
  cat(paste("\n",sprintf(format_column_1,"property"),               sprintf(format_column_2,"min"),                              sprintf(format_column_3,"max"), sep=""),                                     file=summaryFile, sep="\n", append=TRUE)
  if (p==1)
  {
    cat(paste(sprintf(format_column_1,"Homozygous (a)"),            sprintf(format_column_2,percentage_format(homo[2])),         sprintf(format_column_3,percentage_format(homo[1])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
  }
  if (p==2)
  {
    cat(paste(sprintf(format_column_1,"Homozygous (aa)"),           sprintf(format_column_2,percentage_format(homo[2])),         sprintf(format_column_3,percentage_format(homo[1])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"Heterozygous (ab)"),         sprintf(format_column_2,percentage_format(hets[[1]][1])),         sprintf(format_column_3,percentage_format(hets[[1]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
  }
  if (p==3)
  {
    cat(paste(sprintf(format_column_1,"Homozygous (aaa)"),          sprintf(format_column_2,percentage_format(homo[2])),         sprintf(format_column_3,percentage_format(homo[1])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"Heterozygous (not aaa)"),    sprintf(format_column_2,percentage_format(het[1])),          sprintf(format_column_3,percentage_format(het[2])),  sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"aab"),                       sprintf(format_column_2,percentage_format(hets[[1]][1])),         sprintf(format_column_3,percentage_format(hets[[1]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"abc"),                       sprintf(format_column_2,percentage_format(hets[[2]][1])),         sprintf(format_column_3,percentage_format(hets[[2]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
  }
  if (p==4)
  {
    cat(paste(sprintf(format_column_1,"Homozygous (aaaa)"),                    sprintf(format_column_2,percentage_format(homo[2])),      sprintf(format_column_3,percentage_format(homo[1])),      sep=""), file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"Heterozygous (not aaaa)"),              sprintf(format_column_2,percentage_format(het[1])),       sprintf(format_column_3,percentage_format(het[2])),       sep=""), file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1, switch(top+1, "aaab", "aaab", "aabb")), sprintf(format_column_2,percentage_format(hets[[1]][1])), sprintf(format_column_3,percentage_format(hets[[1]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1, switch(top+1, "aabb", "aabc", "aabc")), sprintf(format_column_2,percentage_format(hets[[2]][1])), sprintf(format_column_3,percentage_format(hets[[2]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1, switch(top+1, "aabc", "abcd", "abcd")), sprintf(format_column_2,percentage_format(hets[[3]][1])), sprintf(format_column_3,percentage_format(hets[[3]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
    if (top == 0) {
      cat(paste(sprintf(format_column_1,"abcd"),                               sprintf(format_column_2,percentage_format(hets[[4]][1])), sprintf(format_column_3,percentage_format(hets[[4]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
    }
  }
  if (p==5)
  {
    cat(paste(sprintf(format_column_1,"Homozygous (aaaaa)"),                                                sprintf(format_column_2,percentage_format(ahomo)), sep=""),      file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"Heterozygous (not aaaaa)"),                                          sprintf(format_column_2,percentage_format(ahet)),  sep=""),      file=summaryFile, sep="\n", append=TRUE)
    #cat(paste(sprintf(format_column_1, switch(top+1,"aaaab", "aaaab", "aaaab", "aaabb", "aaabb", "aaabb")), sprintf(format_column_2,percentage_format(hets[[1]][1])), sprintf(format_column_3,percentage_format(hets[[1]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
    #cat(paste(sprintf(format_column_1, switch(top+1,"aaabb", "aaabc", "aabbc", "aaabc", "aabcc", "aabcc")), sprintf(format_column_2,percentage_format(hets[[2]][1])), sprintf(format_column_3,percentage_format(hets[[2]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
    #cat(paste(sprintf(format_column_1, switch(top+1,"aaabc", "aabcd", "aabcd", "aabcd", "aabcd", "abcdd")), sprintf(format_column_2,percentage_format(hets[[3]][1])), sprintf(format_column_3,percentage_format(hets[[3]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
    #cat(paste(sprintf(format_column_1, switch(top+1,"aabbc", "abcde", "abcde", "abcde", "abcde", "abcde")), sprintf(format_column_2,percentage_format(hets[[4]][1])), sprintf(format_column_3,percentage_format(hets[[4]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
    if (top == 0) {
      #cat(paste(sprintf(format_column_1,"aabcd"),                     sprintf(format_column_2,percentage_format(hets[[5]][1])),         sprintf(format_column_3,percentage_format(hets[[5]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
      #cat(paste(sprintf(format_column_1,"abcde"),                     sprintf(format_column_2,percentage_format(hets[[6]][1])),         sprintf(format_column_3,percentage_format(hets[[6]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    }
  }
  if (p==6)
  {
    cat(paste(sprintf(format_column_1,"Homozygous (aaaaaa)"),       sprintf(format_column_2,percentage_format(ahomo)), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"Heterozygous (not aaaaaa)"), sprintf(format_column_2,percentage_format(ahet)), sep=""),                 file=summaryFile, sep="\n", append=TRUE)
    #cat(paste(sprintf(format_column_1, switch(top+1, "aaaaab", "aaaaab:", "aaaaab:", "aaaaab:", "aaaaab:", "aaaaab:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaabbb:", "aaabbb:", "aaabbb:")),                    sprintf(format_column_2,percentage_format(hets[[1]][1])),         sprintf(format_column_3,percentage_format(hets[[1]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    #cat(paste(sprintf(format_column_1, switch(top+1, "aaaabb", "aaaabc:", "aaaabc:", "aaabbc:", "aaabbc:", "aaabbc:", "aaaabc:", "aaaabc:", "aaabcc:", "aaabcc:", "aaabcc:", "aabbcc:", "aabbcc:", "aabbcc:", "aaabbc:", "aaabbc:", "aaabbc:")),                    sprintf(format_column_2,percentage_format(hets[[2]][1])),         sprintf(format_column_3,percentage_format(hets[[2]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    #cat(paste(sprintf(format_column_1, switch(top+1, "aaabbb", "aaabcd:", "aabbcd:", "aaabcd:", "aabccd:", "aabccd:", "aaabcd:", "aabbcd:", "aaabcd:", "aabcdd:", "aabcdd:", "aabbcd:", "aabcdd:", "aabcdd:", "aaabcd:", "aabccd:", "aabccd:")),                    sprintf(format_column_2,percentage_format(hets[[3]][1])),         sprintf(format_column_3,percentage_format(hets[[3]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    #cat(paste(sprintf(format_column_1, switch(top+1, "aaaabc", "aabcde:", "aabcde:", "aabcde:", "aabcde:", "abcdde:", "aabcde:", "aabcde:", "aabcde:", "aabcde:", "abcdee:", "aabcde:", "aabcde:", "abcdee:", "aabcde:", "aabcde:", "abcdde:")),                    sprintf(format_column_2,percentage_format(hets[[4]][1])),         sprintf(format_column_3,percentage_format(hets[[4]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    #cat(paste(sprintf(format_column_1, switch(top+1, "aaabbc", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:")),                    sprintf(format_column_2,percentage_format(hets[[5]][1])),         sprintf(format_column_3,percentage_format(hets[[5]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    if (top==0) {
      #cat(paste(sprintf(format_column_1,"aabbcc"),                    sprintf(format_column_2,percentage_format(hets[[6]][1])),         sprintf(format_column_3,percentage_format(hets[[6]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
      #cat(paste(sprintf(format_column_1,"aaabcd"),                    sprintf(format_column_2,percentage_format(hets[[7]][1])),         sprintf(format_column_3,percentage_format(hets[[7]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
      #cat(paste(sprintf(format_column_1,"aabbcd"),                    sprintf(format_column_2,percentage_format(hets[[8]][1])),         sprintf(format_column_3,percentage_format(hets[[8]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
      #cat(paste(sprintf(format_column_1,"aabcde"),                    sprintf(format_column_2,percentage_format(hets[[9]][1])),         sprintf(format_column_3,percentage_format(hets[[9]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
      #cat(paste(sprintf(format_column_1,"abcdef"),                    sprintf(format_column_2,percentage_format(hets[[10]][1])),        sprintf(format_column_3,percentage_format(hets[[10]][2])), sep=""),               file=summaryFile, sep="\n", append=TRUE)
    }
  }
  cat(paste(sprintf(format_column_1,"Genome Haploid Length"), sprintf(format_column_2,bp_format(total_len[2])),                  sprintf(format_column_3,bp_format(total_len[1])), sep=""),                   file=summaryFile, sep="\n", append=TRUE)
  cat(paste(sprintf(format_column_1,"Genome Repeat Length"),  sprintf(format_column_2,bp_format(repeat_len[2])),                 sprintf(format_column_3,bp_format(repeat_len[1])), sep=""),                  file=summaryFile, sep="\n", append=TRUE)
  cat(paste(sprintf(format_column_1,"Genome Unique Length"),  sprintf(format_column_2,bp_format(unique_len[2])),                 sprintf(format_column_3,bp_format(unique_len[1])), sep=""),                  file=summaryFile, sep="\n", append=TRUE)
  cat(paste(sprintf(format_column_1,"Model Fit "),            sprintf(format_column_2,percentage_format(model_fit_allscore[1])), sprintf(format_column_3,percentage_format(model_fit_fullscore[1])), sep=""), file=summaryFile, sep="\n", append=TRUE)
  cat(paste(sprintf(format_column_1,"Read Error Rate"),       sprintf(format_column_2,percentage_format(error_rate[1])),         sprintf(format_column_3,percentage_format(error_rate[2])), sep=""),          file=summaryFile, sep="\n", append=TRUE)
  if (VERBOSE)
  {
    cat(paste("\nPercent Kmers Modeled (All Kmers) = ",  percentage_format(model_fit_allscore[1]),    " [", model_fit_allscore[2],    ", ", model_fit_allscore[3],    "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
    cat(paste("Percent Kmers Modeled (Full Model) = ",   percentage_format(model_fit_fullscore[1]),   " [", model_fit_fullscore[2],   ", ", model_fit_fullscore[3],   "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
    cat(paste("Percent Kmers Modeled (Unique Kmers) = ", percentage_format(model_fit_uniquescore[1]), " [", model_fit_uniquescore[2], ", ", model_fit_uniquescore[3], "]", sep=""), file=summaryFile, sep="\n", append=TRUE)

    cat(paste("\nModel RSSE (All Kmers) = ",  model_fit_all[1],    " [", model_fit_all[2],    ", ", model_fit_all[3],    "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
    cat(paste("Model RSSE (Full Model) = ",   model_fit_full[1],   " [", model_fit_full[2],   ", ", model_fit_full[3],   "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
    cat(paste("Model RSSE (Unique Model) = ", model_fit_unique[1], " [", model_fit_unique[2], ", ", model_fit_unique[3], "]", sep=""), file=summaryFile, sep="\n", append=TRUE)	
  }
  ## Finalize the progress
  progressFilename=paste(foldername, "/", arguments$name_prefix, "progress.txt",sep="")
  cat(model_status, file=progressFilename, sep="\n", append=TRUE)

  if (TESTING) {
    if (TRUE_PARAMS!=-1) {
      testingFile <- paste(foldername, "/", arguments$name_prefix, "SIMULATED_testing.tsv",sep="")
    } else {
      testingFile <- paste(foldername,"/SIMULATED_testing.tsv",sep="")
    }
    if (p==1) {
      cat(paste(amd, akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
    }
    if (p==2) {
      cat(paste(amd, ahets[[1]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
    }
    if (p==3) {
      if (TRUE_PARAMS!=-1) {
        true_params = unlist(lapply(strsplit(TRUE_PARAMS, ","), as.numeric))
        cat(paste(amd, ahets[[1]], ahets[[2]], akcov, adups, atotal_len, top, true_params[1], true_params[2], true_params[3], true_params[4], sep="\t"), file=testingFile, sep="\n", append=FALSE)
      } else {
        cat(paste(amd, ahets[[1]], ahets[[2]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
      }
    }
    if (p==4) {
      if (topology==0) {
        if (TRUE_PARAMS!=-1) {
          true_params = unlist(lapply(strsplit(TRUE_PARAMS, ","), as.numeric))
          cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], akcov, adups, atotal_len, top, true_params[1], true_params[2], true_params[3], true_params[4], true_params[5], true_params[6], sep="\t"), file=testingFile, sep="\n", append=FALSE)
        } else {
          cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
        }
      } else {
        if (TRUE_PARAMS!=-1) {
          true_params = unlist(lapply(strsplit(TRUE_PARAMS, ","), as.numeric))
          cat(paste(amd, switch(top, ahets[[1]], 0), switch(top, 0, ahets[[1]]), ahets[[2]], ahets[[3]], akcov, adups, atotal_len, top, true_params[1], switch(true_params[5], true_params[2], 0), switch(true_params[5], 0, true_params[2]), true_params[3], true_params[4], true_params[5], sep="\t"), file=testingFile, sep="\n", append=FALSE)
        } else {
          cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
        }
      }
    }
    if (p==5) {
      if (topology==0) {
        if (TRUE_PARAMS!=-1) {
          true_params = unlist(lapply(strsplit(TRUE_PARAMS, ","), as.numeric))
          cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], ahets[[6]], akcov, adups, atotal_len, top, true_params[1], true_params[2], true_params[3], true_params[4], true_params[5], true_params[6], true_params[7], true_params[8], sep="\t"), file=testingFile, sep="\n", append=FALSE)
        } else {
          cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], ahets[[6]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
        }
      } else {
        if (TRUE_PARAMS!=-1) {
          true_params = unlist(lapply(strsplit(TRUE_PARAMS, ","), as.numeric))
          cat(paste(amd, switch(top, ahets[[1]], ahets[[1]], 0, 0, 0), switch(top, 0, 0, ahets[[1]], ahets[[1]], ahets[[1]]), switch(top, ahets[[2]], 0, ahets[[2]], 0, 0), switch(top, 0, ahets[[2]], 0, 0, 0), switch(top, 0, 0, 0, ahets[[2]], ahets[[2]]), switch(top, ahets[[3]], ahets[[3]], ahets[[3]], ahets[[3]], 0), switch(top, 0, 0, 0, 0, ahets[[3]]), ahets[[4]], akcov, adups, atotal_len, top, true_params[1], switch(true_params[6], true_params[2], true_params[2], 0, 0, 0), switch(true_params[6], 0, 0, true_params[2], true_params[2], true_params[2]), switch(true_params[6], true_params[3], 0, true_params[3], 0, 0), switch(true_params[6], 0, true_params[3], 0, 0, 0), switch(true_params[6], 0, 0, 0, true_params[3], true_params[3]), switch(true_params[6], true_params[4], true_params[4], true_params[4], true_params[4], 0), switch(true_params[6], 0, 0, 0, 0, true_params[4]), true_params[5], true_params[6], sep="\t"), file=testingFile, sep="\n", append=FALSE)
        } else {
          cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
        }
      }
    }
    if (p==6) {
      if (topology==0) {
        if (TRUE_PARAMS!=-1) {
          true_params = unlist(lapply(strsplit(TRUE_PARAMS, ","), as.numeric))
          cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], ahets[[6]], ahets[[7]], ahets[[8]], ahets[[9]], ahets[[10]], akcov, adups, atotal_len, top, true_params[1], true_params[2], true_params[3], true_params[4], true_params[5], true_params[6], true_params[7], true_params[8], true_params[9], true_params[10], true_params[11], true_params[12], sep="\t"), file=testingFile, sep="\n", append=FALSE)
        } else {
          cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], ahets[[6]], ahets[[7]], ahets[[8]], ahets[[9]], ahets[[10]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
        }
      } else {
        if (TRUE_PARAMS!=-1) {
          true_params = unlist(lapply(strsplit(TRUE_PARAMS, ","), as.numeric))
          cat(paste(amd, ifelse(top %in% c(1,2,3,4,5), ahets[[1]], 0), ifelse(top %in% c(6,7,8,9,10,11,12,13), ahets[[1]], 0), ifelse(top %in% c(14,15,16), ahets[[1]], 0), ifelse(top %in% c(1,2,6,7), ahets[[2]], 0), ifelse(top %in% c(3,4,5,14,15,16), ahets[[2]], 0), ifelse(top %in% c(8,9,10), ahets[[2]], 0), ifelse(top %in% c(11,12,13), ahets[[2]], 0), ifelse(top %in% c(1,3,6,8,14), ahets[[3]], 0), ifelse(top %in% c(2,7,11), ahets[[3]], 0), ifelse(top %in% c(4,5,15,16), ahets[[3]], 0), ifelse(top %in% c(9,10,12,13), ahets[[3]], 0), ifelse(top %in% c(1,2,3,4,6,7,8,9,11,12,14,15), ahets[[4]], 0), ifelse(top %in% c(5,16), ahets[[4]], 0), ifelse(top %in% c(10,13), ahets[[4]], 0), ahets[[5]], akcov, adups, atotal_len, top, true_params[1], ifelse(true_params[7] %in% c(1,2,3,4,5), true_params[2], 0), ifelse(true_params[7] %in% c(6,7,8,9,10,11,12,13), true_params[2], 0), ifelse(true_params[7] %in% c(14,15,16), true_params[2], 0), ifelse(true_params[7] %in% c(1,2,6,7), true_params[3], 0), ifelse(true_params[7] %in% c(3,4,5,14,15,16), true_params[3], 0), ifelse(true_params[7] %in% c(8,9,10), true_params[3], 0), ifelse(true_params[7] %in% c(11,12,13), true_params[3], 0), ifelse(true_params[7] %in% c(1,3,6,8,14), true_params[4], 0), ifelse(true_params[7] %in% c(2,7,11), true_params[4], 0), ifelse(true_params[7] %in% c(4,5,15,16), true_params[4], 0), ifelse(true_params[7] %in% c(9,10,12,13), true_params[4], 0), ifelse(true_params[7] %in% c(1,2,3,4,6,7,8,9,11,12,14,15), true_params[5], 0), ifelse(true_params[7] %in% c(5,16), true_params[5], 0), ifelse(true_params[7] %in% c(10,13), true_params[5], 0), true_params[6], true_params[7], sep="\t"), file=testingFile, sep="\n", append=FALSE)
        } else {
          cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
        }
      }
    }
  }
}
