#!/bin/env Rscript

## GenomeScope: Fast Genome Analysis from Unassembled Short Reads
## This is the automated script for computing genome characteristics from a histogram file, k-mer size, and readlength

## Number of rounds before giving up
NUM_ROUNDS=4

## Coverage steps to trim off between rounds
START_SHIFT=5

## Typical cutoff for sequencing error
TYPICAL_ERROR = 15

## Max rounds on NLS
MAX_ITERATIONS=20

## Print out VERBOSEging messages (0/1)
VERBOSE = 0

## Suppress the warnings if the modeling goes crazy, those are in try/catch blocks anyways
options(warn=-1)

## Colors for plots
COLOR_BGCOLOR  = "light grey"
COLOR_HIST     = "#56B4E9"
COLOR_4PEAK    = "black"
COLOR_2PEAK    = "#F0E442"
COLOR_ERRORS   = "#D55E00"
COLOR_KMERPEAK = "black"


## Helper Function
###############################################################################

min_max <- function(table){
	return (c( abs(table[1]) - abs(table[2]) , abs(table[1])+abs(table[2])))
}



## Use nls to fit 4 peak model
###############################################################################

nls_4peak<-function(x, y, k, estKmercov, estLength, max_iterations){
	model4 = NULL

    if (VERBOSE) { cat("trying nls_4peak standard algorithm\n") }
	try(model4 <- nls(y ~ ((2.0 * (1-d) * (1-(1-r)^k))            * dnbinom(x, size = kmercov   / bias, mu = kmercov)     * length +
                          ((d*(1-(1-r)^k)^2) + (1-2*d)*((1-r)^k)) * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2) * length + 
                          (2*d*((1-r)^k)*(1-(1-r)^k))             * dnbinom(x, size = kmercov*3 / bias, mu = kmercov * 3) * length + 
                          (d*(1-r)^(2*k))                         * dnbinom(x, size = kmercov*4 / bias, mu = kmercov * 4) * length), 
                      start = list(d=0, r=0, kmercov=estKmercov, bias = 0.5, length=estLength),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

	if(class(model4) == "try-error"){
        if (VERBOSE) { cat("retrying nls_4peak with port algorithm\n") }
		try(model4 <- nls(y ~ ((2.0 * (1-d) * (1-(1-r)^k))            * dnbinom(x, size = kmercov   / bias, mu = kmercov)     * length + 
                              ((d*(1-(1-r)^k)^2) + (1-2*d)*((1-r)^k)) * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2) * length + 
                              (2*d*((1-r)^k)*(1-(1-r)^k))             * dnbinom(x, size = kmercov*3 / bias, mu = kmercov * 3) * length + 
                              (d*(1-r)^(2*k))                         * dnbinom(x, size = kmercov*4 / bias, mu = kmercov * 4) * length), 
                          start = list(d=0, r=0, kmercov=estKmercov, bias = 0.5, length=estLength),
                          algorithm = "port", control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)
	}

	return(model4)
}


## score model by number and percent of residual errors
###############################################################################

score_model<-function(kmer_hist_orig, nls, round, foldername){
  x = kmer_hist_orig[[1]]
  y = kmer_hist_orig[[2]]

  pred=predict(nls, newdata=data.frame(x))
  model_sum=summary(nls)
  kcov = min_max(model_sum$coefficients['kmercov',])
  
  ## Compute error rate, by counting kmers unexplained by model through first peak
  ## truncate errors as soon as it goes to zero, dont allow it to go back up
  error_xcutoff = floor(kcov[1])
  error_xcutoff_ind = which(x==error_xcutoff)

  error_kmers = y[1:error_xcutoff_ind] - pred[1:error_xcutoff_ind]

  first_zero = -1

  for (i in 1:error_xcutoff_ind)
  {
    if (first_zero == -1)
    {
      if (error_kmers[i] < 1.0)
      {
        first_zero = i
        if (VERBOSE) { cat(paste("Truncating errors at", i, "\n")) }
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
  }

  suby = y[first_zero:length(y)]
  subp = pred[first_zero:length(y)]
  residual = as.numeric(suby)-as.numeric(subp)

  residualslow = y
  for (i in 1:length(residualslow))
  {
    residualslow[[i]] = y[[i]] - pred[[i]]
  }
  
  str(head(y[first_zero:length(y)]))
  str(head(pred[first_zero:length(y)]))
  str(head(residual))
  cat(head(residualslow))
  cat(suby[[1]]-subp[[1]], suby[[2]]-subp[[2]], suby[[3]] - subp[[3]], suby[[4]]-subp[[4]])

  ## The fit is residual sum of square error, excluding sequencing errors
  model_fit_all         = c(sum(as.numeric(y[first_zero:length(y)] - pred[first_zero:length(y)]) ** 2), first_zero, x[length(y)])
  model_fit_full        = c(sum(as.numeric(y[first_zero:5*kcov[1]] - pred[first_zero:5*kcov[1]]) ** 2), first_zero, 5*kcov[1])
  model_fit_unique      = c(sum(as.numeric(y[first_zero:3*kcov[1]] - pred[first_zero:3*kcov[1]]) ** 2), first_zero, 3*kcov[1])

  ## The score is the percentage of unexplained kmers, excluding sequencing errors
  model_fit_allscore    = c(sum(abs(as.numeric(y[first_zero:length(y)] - pred[first_zero:length(y)]))) / sum(as.numeric(y[first_zero:length(y)])), first_zero, x[length(y)])
  model_fit_fullscore   = c(sum(abs(as.numeric(y[first_zero:5*kcov[1]] - pred[first_zero:5*kcov[1]]))) / sum(as.numeric(y[first_zero:5*kcov[1]])), first_zero, 5*kcov[1])
  model_fit_uniquescore = c(sum(abs(as.numeric(y[first_zero:3*kcov[1]] - pred[first_zero:3*kcov[1]]))) / sum(as.numeric(y[first_zero:3*kcov[1]])), first_zero, 3*kcov[1])

  fit = data.frame(all  = model_fit_all,      allscore  = model_fit_allscore,
                   full = model_fit_full,     fullscore = model_fit_fullscore, 
                   unique = model_fit_unique, uniquescore = model_fit_uniquescore)

  pdf(paste(foldername, "/score.", round, ".pdf", sep=""))
  plot(x, y, xlim=c(0,150), ylim=c(0,5e6), typ="n", col=COLOR_HIST)
  lines(x, y, typ="l",  xlim=c(0,150), ylim=c(0,5e6), lwd=2, col=COLOR_HIST)
  lines(x, pred, typ="l", lwd=2, col=COLOR_4PEAK)
  lines(x[first_zero:length(y)], residual, col="purple")
  lines(x, residualslow, col="green", lty=2)
  abline(v=c(first_zero, 3*kcov[1], 5*kcov[1]), col="red")
  abline(h=0, col="grey")
  abline(v=seq(0,150,5), col="grey")
  dev.off()

  return (fit)
}



## ## Pick n the two model forms, resolves ambiguity between which is the homozygous and which is the heterozygous peak
## ###############################################################################
## 
## eval_model<-function(kmer_hist_orig, nls1, nls2, round, foldername){
##     nls1score = -1
##     nls2score = -1
##     
##     ## Count the kmers in the observed data in the interesting range
##     ox = kmer_hist_orig[[1]]
##     oy = kmer_hist_orig[[2]]
## 
##     allkmers = sum(as.numeric(ox * oy))
##     #allkmers = sum(as.numeric(oy))
## 
##     if(VERBOSE){ cat(paste("allkmers:\t", allkmers, "\n"))}
## 
##     ## Evaluate the score the nls1
##     if (!is.null(nls1))
##     {
##       res1 <- predict(nls1, newdata=data.frame(ox))
##       if(VERBOSE) { cat(paste("nls1 kmers:\t", sum(as.numeric(ox*res1)), "\n")) }
##      # if(VERBOSE) { cat(paste("nls1 kmers:\t", sum(as.numeric(res1)), "\n")) }
## 
##       nls1score = sum(as.numeric(abs(ox*oy-ox*res1))) / allkmers
##       #nls1score = sum(as.numeric(abs(oy-res1))) / allkmers
##       if(VERBOSE){ cat(paste("nls1score:\t", nls1score, "\n"))}
## 
##       if (VERBOSE)
##       {
##         mdir = paste(foldername, "/round", round, ".1", sep="")
##         dir.create(mdir, showWarnings=FALSE)
##         report_results(kmer_prof_orig, k, (list(nls1, nls1score)) , mdir)
##       }
##     }
##     else
##     {
##       if (VERBOSE) { cat("nls1score failed to converge\n") }
##     }
## 
##     
##     ## Evaluate the score of nls2
##     if (!is.null(nls2))
##     {
##       res2 <- predict(nls2, newdata=data.frame(ox))
##       if(VERBOSE) { cat(paste("nls2 kmers:\t", sum(as.numeric(ox*res2)), "\n")) }
##       # if(VERBOSE) { cat(paste("nls2 kmers:\t", sum(as.numeric(res2)), "\n")) }
## 
##       nls2score = sum(as.numeric(abs(ox*oy-ox*res2))) / allkmers
##       #nls2score = sum(as.numeric(abs(oy-res2))) / allkmers
##       if(VERBOSE){ cat(paste("nls2score:\t", nls2score, "\n"))}
## 
##       if (VERBOSE)
##       {
##         mdir = paste(foldername, "/round", round, ".2", sep="")
##         dir.create(mdir, showWarnings=FALSE)
##         report_results(kmer_prof_orig, k, (list(nls2, nls2score)) , mdir)
##       }
##     }
##     else
##     {
##       if (VERBOSE) { cat("nls2score failed to converge\n") }
##     }
## 
## 
##     ## Return the better of the scores
##     if (!is.null(nls1))
##     {
##       if (!is.null(nls2))
##       {
##         if (nls1score < nls2score)
##         {
##           if (VERBOSE) { cat(paste("returning nls1, better score\n")) }
##           return (list(nls1, nls1score))
##         }
##         else
##         {
##           if (VERBOSE) { cat(paste("returning nls2, better score\n")) }
##           return (list(nls2, nls2score))
##         }
##       }
##       else
##       {
##         if (VERBOSE) { cat(paste("returning nls1, nls2 fail\n")) }
##         return (list(nls1, nls1score))
##       }
##     }
## 
##     if (VERBOSE) { cat(paste("returning nls2 by default\n")) }
##     return (list(nls2, nls2score))
## }

## Pick between the two model forms, resolves ambiguity between which is the homozygous and which is the heterozygous peak
###############################################################################

eval_model<-function(kmer_hist_orig, nls1, nls2, round, foldername){
    nls1score = -1
    nls2score = -1
    
    ## Evaluate the score the nls1
    if (!is.null(nls1))
    {
      nls1score = score_model(kmer_hist_orig, nls1, round+0.1, foldername)

      if(VERBOSE){ cat(paste("nls1score$all:\t", nls1score$all[[1]], "\n"))}

      if (VERBOSE)
      {
        mdir = paste(foldername, "/round", round, ".1", sep="")
        dir.create(mdir, showWarnings=FALSE)
        report_results(kmer_prof_orig, k, (list(nls1, nls1score)) , mdir)
      }
    }
    else
    {
      if (VERBOSE) { cat("nls1score failed to converge\n") }
    }

    
    ## Evaluate the score of nls2
    if (!is.null(nls2))
    {
      nls2score = score_model(kmer_hist_orig, nls2, round+0.2, foldername)

      if(VERBOSE){ cat(paste("nls2score$all:\t", nls2score$all[[1]], "\n"))}

      if (VERBOSE)
      {
        mdir = paste(foldername, "/round", round, ".2", sep="")
        dir.create(mdir, showWarnings=FALSE)
        report_results(kmer_prof_orig, k, (list(nls2, nls2score)) , mdir)
      }
    }
    else
    {
      if (VERBOSE) { cat("nls2score failed to converge\n") }
    }


    ## Return the better of the scores
    if (!is.null(nls1))
    {
      if (!is.null(nls2))
      {
        if (nls1score$all[[1]] < nls2score$all[[1]])
        {
          if (VERBOSE) { cat(paste("returning nls1, better score\n")) }
          return (list(nls1, nls1score))
        }
        else
        {
          if (VERBOSE) { cat(paste("returning nls2, better score\n")) }
          return (list(nls2, nls2score))
        }
      }
      else
      {
        if (VERBOSE) { cat(paste("returning nls1, nls2 fail\n")) }
        return (list(nls1, nls1score))
      }
    }

    if (VERBOSE) { cat(paste("returning nls2 by default\n")) }
    return (list(nls2, nls2score))
}


## Wrapper function to try fitting 4 peak model with 2 forms
###############################################################################

estimate_Genome_4peak2<-function(kmer_hist_orig, x, y, k, readlength, round, foldername){
	## First we see what happens when the max peak is the kmercoverage (typically the homozygous peak) for the plot
	numofReads   = sum(as.numeric(x*y))/(readlength-k+1) 
	estKmercov1  = x[which(y==max(y))][1]
	estCoverage1 = estKmercov1*readlength/(readlength-k)
	estLength1   = numofReads*readlength/estCoverage1

    if (VERBOSE) { cat(paste("trying with kmercov: ", estKmercov1, "\n")) }
	nls1    = nls_4peak(x, y, k, estKmercov1, estLength1, MAX_ITERATIONS)
    if (VERBOSE) { print(summary(nls1)) }

	## Second we half the max kmercoverage (typically the heterozygous peak)
	estKmercov2  = estKmercov1 / 2 ##2.5
	estCoverage2 = estKmercov2*readlength/(readlength-k)
	estLength2   = numofReads*readlength/estCoverage2

    if (VERBOSE) { cat(paste("trying with kmercov: ", estKmercov2, "\n")) }
	nls2 = nls_4peak(x, y, k, estKmercov2, estLength2, MAX_ITERATIONS)
    if (VERBOSE) { print(summary(nls2)) }

	return(eval_model(kmer_hist_orig, nls1, nls2, round, foldername))
}


## Format numbers
###############################################################################
bp_format<-function(num) {
  paste(formatC(round(num),format="d",big.mark=","), "bp",sep=" ")
}

percentage_format<-function(num) {
  paste(signif(num,6)*100,"%",sep="")
}
X_format<-function(num) {
  paste(signif(num,4),"X",sep="")  
}


## Report results and make plots
###############################################################################

report_results<-function(kmer_hist, k, container, foldername)
{
    x=kmer_hist[[1]]
    y=kmer_hist[[2]]
	#d=data.frame(x=x, y=y)

	#automatically zoom into the relevant regions of the plot, ignore first 15 positions
    xmax=length(x)
	start=which(y == min(y[1:TYPICAL_ERROR]))
	zoomx=x[start:xmax-1]
	zoomy=y[start:xmax-1]

    ## allow for a little space above max value past the noise
	y_limit = max(zoomy[start:length(zoomy)])*1.1
	
	x_limit = which(y == max(y[start:length(zoomx)])) * 3

	if (min(zoomy) > zoomy[1]){
		x_limit=max(which(zoomy<zoomy[1])[2],600)
	}

    if (!is.null(container[[1]]))
    {
       model_sum=summary(container[[1]])
       kcov = min_max(model_sum$coefficients['kmercov',])[1]
       x_limit = max(kcov*5.1, x_limit)
    }

    ## Uncomment this to enforce a specific number
    # x_limit=150

    ## Features to report
    het=c(-1,-1)
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
	png(paste(foldername, "/plot.png", sep=""),width=plot_size,height=plot_size, res=resolution)
	plot(kmer_hist, type="n", main="GenomeScope Profile\n", xlab="Coverage", ylab="Frequency", ylim=c(0,y_limit), xlim=c(0,x_limit),cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
    rect(-1e10, -1e10, max(kmer_hist[[1]])*1.1 , max(kmer_hist[[2]])*1.1, col=COLOR_BGCOLOR)
    points(kmer_hist, type="h", col=COLOR_HIST, lwd=2)
    box(col="black")

    ## Make a second plot in log space over entire range
	png(paste(foldername, "/plot.log.png", sep=""),width=plot_size,height=plot_size,res=resolution)
	plot(kmer_hist, type="n", main="GenomeScope Profile\n", xlab="Coverage", ylab="Frequency", log="xy",cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
    rect(1e-10, 1e-10, max(kmer_hist[[1]])*10 , max(kmer_hist[[2]])*10, col=COLOR_BGCOLOR)
	points(kmer_hist, type="h", col=COLOR_HIST, lwd=2)
    box(col="black")

	if(!is.null(container[[1]])) 
    { 
       ## The model converged!
       pred=predict(container[[1]], newdata=data.frame(x))

       ## Compute the genome characteristics
       model_sum=summary(container[[1]])
       
       ## save the model to a file
       capture.output(model_sum, file=paste(foldername,"/model.txt", sep=""))
       
       ## Identify key values
       het  = min_max(model_sum$coefficients['r',])
       dups = min_max(model_sum$coefficients['bias',])
       kcov = min_max(model_sum$coefficients['kmercov',])
       mlen = min_max(model_sum$coefficients['length',])
       md   = min_max(model_sum$coefficients['d',])
       
       ## Compute error rate, by counting kmers unexplained by model through first peak
       ## truncate errors as soon as it goes to zero, dont allow it to go back up
       error_xcutoff = floor(kcov[1])
       error_xcutoff_ind = which(x==error_xcutoff)

       error_kmers = y[1:error_xcutoff_ind] - pred[1:error_xcutoff_ind]

       first_zero = -1

       for (i in 1:error_xcutoff_ind)
       {
         if (first_zero == -1)
         {
           if (error_kmers[i] < 1.0)
           {
             first_zero = i
             if (VERBOSE) { cat(paste("Truncating errors at", i, "\n")) }
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
       }

       ## Rather than "0", set to be some very small number so log-log plot looks okay
       error_kmers = pmax(error_kmers, 1e-10)

       total_error_kmers = sum(as.numeric(error_kmers * x[1:error_xcutoff_ind]))
       total_kmers = sum(as.numeric(x*y))

       f1 <- function(x){
             i=seq(1,k) 
             h=(1-x)^(k-i)*x^i*choose(k,i)
             sum(h)*total_kmers-total_error_kmers
       }

       error_rate_root = try(uniroot(f1, c(0,1))$root)

       if (class(error_rate_root) == "try-error")
       {
         error_rate  = c(total_error_kmers/total_kmers/k, total_error_kmers/total_kmers/k)
       }
       else
       {
          error_rate  = c(error_rate_root, error_rate_root)
       }

       total_len = (total_kmers-total_error_kmers)/(2*kcov)
       
       ## find kmers that fit the 2 peak model (no repeats)
       unique_hist <- (2 * (1 - md[1]) * (1 - (1 - het[1])^k))                                * dnbinom(x, size = kcov[1]     / dups[1], mu = kcov[1])     * mlen[1] +
                      ((md[1] * (1 - (1 - het[1])^k)^2) + (1 - 2 * md[1]) * ((1 - het[1])^k)) * dnbinom(x, size = kcov[1] * 2 / dups[1], mu = kcov[1] * 2) * mlen[1]
       
       unique_kmers = sum(as.numeric(x*unique_hist))
       repeat_kmers = total_kmers - unique_kmers - total_error_kmers
       
       repeat_len=repeat_kmers/(2*kcov)
       unique_len=unique_kmers/(2*kcov)

       score = container[[2]]

       model_fit_allscore    = score$allscore
       model_fit_fullscore   = score$fullscore
       model_fit_uniquescore = score$uniquescore

       model_fit_all    = score$all
       model_fit_full   = score$full
       model_fit_unique = score$unique
       
       ## Finish Log plot
       title(paste("\nlen:",  prettyNum(total_len[1], big.mark=","), 
                   "bp uniq:", format(100*(unique_len[1]/total_len[1]), digits=3),
                   "% het:",  format(100*het[1], digits=3), 
                   "% kcov:", format(kcov[1], digits=3), 
                   " err:",   format(100*error_rate[1], digits=3), 
                   "% dup:",  format(dups[1], digits=3), sep=""), 
                   cex.main=.85)
       
       ## Mark the modes of the peaks
       abline(v=kcov[1] * c(1,2,3,4), col=COLOR_KMERPEAK, lty=2)
       
       ## Draw just the unique portion of the model
       lines(x, unique_hist, col=COLOR_2PEAK, lty=1, lwd=3)
       lines(x, pred, col=COLOR_4PEAK, lwd=3)
       lines(x[1:error_xcutoff_ind], error_kmers, lwd=3, col=COLOR_ERRORS)

       ## Add legend
       legend(exp(.65 * log(max(x))), 1.0 * max(y), 
              legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks"),
              lty=c("solid", "solid", "solid", "solid", "dashed"),
              lwd=c(3,3,3,3,3),
              col=c(COLOR_HIST, COLOR_4PEAK, COLOR_2PEAK, COLOR_ERRORS, COLOR_KMERPEAK),
              bg="white")

       dev.set(dev.next())

       ## Finish Linear Plot
       title(paste("\nlen:",  prettyNum(total_len[1], big.mark=","), 
                   "bp uniq:", format(100*(unique_len[1]/total_len[1]), digits=3),
                   "% het:",  format(100*het[1], digits=3), 
                   "% kcov:", format(kcov[1], digits=3), 
                   " err:",   format(100*error_rate[1], digits=3), 
                   "% dup:",  format(dups[1], digits=3), sep=""), 
                   cex.main=.85)

       ## Mark the modes of the peaks
       abline(v=kcov[1] * c(1,2,3,4), col=COLOR_KMERPEAK, lty=2)
       
       ## Draw just the unique portion of the model
       lines(x, unique_hist, col=COLOR_2PEAK, lty=1, lwd=3)
       lines(x, pred, col=COLOR_4PEAK, lwd=3)
       lines(x[1:error_xcutoff_ind], error_kmers, lwd=3, col=COLOR_ERRORS)

       ## Add legend
       legend(.65 * x_limit, 1.0 * y_limit, 
              legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks"),
              lty=c("solid", "solid", "solid", "solid", "dashed"),
              lwd=c(3,3,3,3,3),
              col=c(COLOR_HIST, COLOR_4PEAK, COLOR_2PEAK, COLOR_ERRORS, COLOR_KMERPEAK),
              bg="white")
       
       model_status="done"

       cat(paste("Model converged het:", format(het[1], digits=3), 
                 " kcov:", format(kcov[1], digits=3), 
                 " err:", format(error_rate[1], digits=3), 
                 " dup:", format(dups[1], digits=3), 
                 " len:", round(total_len[1]), "\n", sep="")) 
	}
    else
    {
      title("\nFailed to converge")
      dev.set(dev.next())
      title("\nFailed to converge")
      cat("Failed to converge")
    }

	dev.off()
	dev.off()

    ## Write key values to summary file
	summaryFile <- paste(foldername,"/summary.txt",sep="")
	
    format_column_1 = "%-30s"
    format_column_2 = "%-18s"
    format_column_3 = "%-18s"
    
    cat(paste("k = ", k,sep=""),                                                                                                                                                          file=summaryFile, sep="\n") 
    cat(paste("\n",sprintf(format_column_1,"property"), sprintf(format_column_2,"min"), sprintf(format_column_3,"max"), sep=""),                                                          file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"Heterozygosity"), sprintf(format_column_2,percentage_format(het[1])), sprintf(format_column_3,percentage_format(het[2])), sep=""),                 file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"Genome Haploid Length"), sprintf(format_column_2,bp_format(total_len[2])), sprintf(format_column_3,bp_format(total_len[1])), sep=""),              file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"Genome Repeat Length"), sprintf(format_column_2,bp_format(repeat_len[2])), sprintf(format_column_3,bp_format(repeat_len[1])), sep=""),             file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"Genome Unique Length"), sprintf(format_column_2,bp_format(unique_len[2])), sprintf(format_column_3,bp_format(unique_len[1])), sep=""),             file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"Read Duplication Level"), sprintf(format_column_2,X_format(dups[1])), sprintf(format_column_3,X_format(dups[2])), sep=""),                         file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"Read Error Rate"), sprintf(format_column_2,percentage_format(error_rate[1])), sprintf(format_column_3,percentage_format(error_rate[2])), sep=""),  file=summaryFile, sep="\n", append=TRUE)
    
    cat(paste("\nModel Score (All Kmers) = ",  model_fit_allscore[1],    " [", model_fit_allscore[2],    " ", model_fit_allscore[3],    "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
    cat(paste("Model Score (Full Model) = ",   model_fit_fullscore[1],   " [", model_fit_fullscore[2],   " ", model_fit_fullscore[3],   "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
    cat(paste("Model Score (Unique Kmers) = ", model_fit_uniquescore[1], " [", model_fit_uniquescore[2], " ", model_fit_uniquescore[3], "]", sep=""), file=summaryFile, sep="\n", append=TRUE)

    cat(paste("Model Error (All Kmers) = ",    model_fit_all[1],    " [", model_fit_all[2],    " ", model_fit_all[3],    "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
    cat(paste("Model Error (Full Model) = ",   model_fit_full[1],   " [", model_fit_full[2],   " ", model_fit_full[3],   "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
    cat(paste("Model Error (Unique Model) = ", model_fit_unique[1], " [", model_fit_unique[2], " ", model_fit_unique[3], "]", sep=""), file=summaryFile, sep="\n", append=TRUE)

    ## Finalize the progress
    progressFilename=paste(foldername,"/progress.txt",sep="")
	cat(model_status, file=progressFilename, sep="\n", append=TRUE)
}




## Main program starts here
###############################################################################

args<-commandArgs(TRUE)

if(length(args) < 4) {
	cat("USAGE: genomescope.R histogram_file k-mer_length read_length output_dir [verbose]\n")
} else{

    ## Load the arguments from the user
	histfile   <- args[[1]]
	k          <- as.numeric(args[[2]])
	readlength <- as.numeric(args[[3]])
	foldername <- args[[4]] 

    if ((length(args) == 5) && (as.numeric(args[[5]] == 1))) { VERBOSE = 1 }

    ## values for testing
    #histfile <- "~/build/genomescope/simulation/simulation_results/Arabidopsis_thaliana.TAIR10.26.dna_sm.toplevel.fa_het0.01_br1_rl100_cov100_err0.01_reads.fa21.hist"
    #k <- 21
    #readlength <- 100
    #foldername <- "~/build/genomescope/simulation/simulation_analysis/Arabidopsis_thaliana.TAIR10.26.dna_sm.toplevel.fa_het0.01_br1_rl100_cov100_err0.01_reads.fa21.hist"

    if (k > readlength) { stop("K cannot be greater than readlength") }

    cat(paste("GenomeScope analyzing ", histfile, " k=", k, " readlen=", readlength, " outdir=", foldername, "\n", sep=""))

	dir.create(foldername, showWarnings=FALSE)

	kmer_prof <- read.csv(file=histfile,sep=" ", header=FALSE) 
	kmer_prof <- kmer_prof[c(1:length(kmer_prof[,2])-1),] #get rid of the last position
    kmer_prof_orig <- kmer_prof

    ## Initialize the status
    progressFilename <- paste(foldername,"/progress.txt",sep="")
	cat("starting", file=progressFilename, sep="\n")

    ## try to find the local minimum between errors and the first (heterozygous) peak
    start <- which(kmer_prof[,2]==min(kmer_prof[1:TYPICAL_ERROR,2]))

    ## terminate after NUM_ROUND iterations, store best result so far in container
	round <- 0
	best_container <- list(NULL,0)

	while(round < NUM_ROUNDS) 
    {
        cat(paste("round", round, "trimming to", start, "trying 4peak model... "), file=progressFilename, sep="", append=TRUE)
        if (VERBOSE) { cat(paste("round", round, "trimming to", start, "trying 4peak model... \n")) }

        ## Reset the input trimming off low frequency error kmers
        max <- length(kmer_prof[,1])
        x <- kmer_prof[start:max,1]
        y <- kmer_prof[start:max,2]

		model_4peaks <- estimate_Genome_4peak2(kmer_prof_orig, x, y, k, readlength, round, foldername) 

        if (!is.null(model_4peaks[[1]])) { 
          cat(paste("converged. score: ", model_4peaks[[2]]$all[[1]]), file=progressFilename, sep="\n", append=TRUE)

          if (VERBOSE)
          {
            mdir = paste(foldername, "/round", round, sep="")
	        dir.create(mdir, showWarnings=FALSE)
	        report_results(kmer_prof_orig, k, model_4peaks, mdir)
          }
        } else {
          cat(paste("unconverged"), file=progressFilename, sep="\n", append=TRUE)
        }

		#check if this result is better than previous
        if (!is.null(model_4peaks[[1]]))
        {
          if ((is.null(best_container[[1]])) || (model_4peaks[[2]]$all[[1]] < best_container[[2]]$all[[1]])){
              best_container = model_4peaks
          }
        }

        ## Ignore a larger number of kmers as errors
        start <- start + START_SHIFT
		round <- round + 1
	}

    ## Report the results, note using the original full profile
	report_results(kmer_prof_orig, k, best_container, foldername)
}
