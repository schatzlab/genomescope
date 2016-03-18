#!/bin/env Rscript

## GenomeScope: Fast Genome Analysis from Unassembled Short Reads
## This is the automated script for reading in a histogram file given the k-mer size, readlength and an optional parameter for what you want the x-axis limit to be.

## Number of rounds before giving up
NUM_ROUNDS=5

## Coverage steps to trim off between rounds
START_SHIFT=5

## Typical cutoff for sequencing error
TYPICAL_ERROR = 15

## Amount of random noise to introduce
RAND_ERROR = 0.0001

## Max rounds on NLS
MAX_ITERATIONS=20

## Print out debugging messages (0/1)
DEBUG = 0

## Suppress the warnings if the modeling goes crazy, those are in try/catch blocks anyways
options(warn=-1)



## Helper Function
###############################################################################

min_max <- function(table){
	return (c( abs(table[1]) - abs(table[2]) , abs(table[1])+abs(table[2])))
}



## Use nls to fit 4 peak model
###############################################################################

nls_4peak<-function(x, y, k, estKmercov, estLength, max_iterations){
	model4 = NULL

    if (DEBUG) { cat("trying nls_4peak standard algorithm\n") }
	try(model4 <- nls(y ~ ((2.0 * (1-d) * (1-(1-r)^k))            * dnbinom(x, size = kmercov   / bias, mu = kmercov)     * length +
                          ((d*(1-(1-r)^k)^2) + (1-2*d)*((1-r)^k)) * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2) * length + 
                          (2*d*((1-r)^k)*(1-(1-r)^k))             * dnbinom(x, size = kmercov*3 / bias, mu = kmercov * 3) * length + 
                          (d*(1-r)^(2*k))                         * dnbinom(x, size = kmercov*4 / bias, mu = kmercov * 4) * length), 
                      start = list(d=0, r=0, kmercov=estKmercov, bias = 0.5, length=estLength),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

	if(class(model4) == "try-error"){
        if (DEBUG) { cat("retrying nls_4peak with port algorithm\n") }
		try(model4 <- nls(y ~ ((2.0 * (1-d) * (1-(1-r)^k))            * dnbinom(x, size = kmercov   / bias, mu = kmercov)     * length + 
                              ((d*(1-(1-r)^k)^2) + (1-2*d)*((1-r)^k)) * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2) * length + 
                              (2*d*((1-r)^k)*(1-(1-r)^k))             * dnbinom(x, size = kmercov*3 / bias, mu = kmercov * 3) * length + 
                              (d*(1-r)^(2*k))                         * dnbinom(x, size = kmercov*4 / bias, mu = kmercov * 4) * length), 
                          start = list(d=0, r=0, kmercov=estKmercov, bias = 0.5, length=estLength),
                          algorithm = "port", control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)
	}

	return(model4)
}



## Pick between the two model forms, resolves ambiguity between which is the homozygous and which is the heterozygous peak
###############################################################################

eval_model<-function(kmer_hist_orig, nls1, nls2){
    nls1score = -1
    nls2score = -1

    ox = kmer_hist_orig[[1]]
    oy = kmer_hist_orig[[2]]

    allkmers = sum(as.numeric(ox * oy))
    if(DEBUG){ cat(paste("allkmers: ", allkmers, "\n"))}

    if (!is.null(nls1))
    {
      res1 <- predict(nls1, newdata=data.frame(ox))
      if(DEBUG) { cat(paste("nls1 kmers: ", sum(as.numeric(ox*res1)), "\n")) }

      nls1score = sum(as.numeric(abs(ox*oy-ox*res1))) / allkmers
      if(DEBUG){ cat(paste("nls1score: ", nls1score, "\n"))}
    }
    else
    {
      if (DEBUG) { cat("nls1score failed to converge\n") }
    }

    if (!is.null(nls2))
    {
      res2 <- predict(nls2, newdata=data.frame(ox))
      if(DEBUG) { cat(paste("nls2 kmers: ", sum(as.numeric(ox*res2)), "\n")) }

      nls2score = sum(as.numeric(abs(ox*oy-ox*res2))) / allkmers
      if(DEBUG){ cat(paste("nls2score: ", nls2score, "\n"))}
    }
    else
    {
      if (DEBUG) { cat("nls2score failed to converge\n") }
    }

    if (!is.null(nls1))
    {
      if (!is.null(nls2))
      {
        if (nls1score < nls2score)
        {
          if (DEBUG) { cat(paste("returning nls1, better score\n")) }
          return (list(nls1, nls1score))
        }
      }
      else
      {
        if (DEBUG) { cat(paste("returning nls1, nls2 fail\n")) }
        return (list(nls1, nls1score))
      }
    }

    if (DEBUG) { cat(paste("returning nls2 by default\n")) }
    return (list(nls2, nls2score))
}


## Wrapper function to try fitting 4 peak model with 2 forms
###############################################################################

estimate_Genome_4peak2<-function(kmer_hist_orig, x, y, k, readlength, foldername){
	## First we see what happens when the max peak is the kmercoverage (typically the homozygous peak) for the plot
	numofReads   = sum(as.numeric(x*y))/(readlength-k+1) 
	estKmercov1  = x[which(y==max(y))][1]
	estCoverage1 = estKmercov1*readlength/(readlength-k)
	estLength1   = numofReads*readlength/estCoverage1

    if (DEBUG) { cat(paste("trying with kmercov: ", estKmercov1, "\n")) }
	nls1    = nls_4peak(x, y, k, estKmercov1, estLength1, MAX_ITERATIONS)
    if (DEBUG) { print(summary(nls1)) }

	## Second we half the max kmercoverage (typically the heterozygous peak)
	estKmercov2  = estKmercov1 / 2 ##2.5
	estCoverage2 = estKmercov2*readlength/(readlength-k)
	estLength2   = numofReads*readlength/estCoverage2

    if (DEBUG) { cat(paste("trying with kmercov: ", estKmercov2, "\n")) }
	nls2 = nls_4peak(x, y, k, estKmercov2, estLength2, MAX_ITERATIONS)
    if (DEBUG) { print(summary(nls2)) }

	return(eval_model(kmer_hist_orig, nls1, nls2))
}


## Wrapper function to try fitting 2 peak model with 2 forms
###############################################################################

estimate_Genome_2peak2<-function(x, y, k, readlength, foldername){
	## First we see what happens when the max peak is the kmercoverage (typically the homozygous peak) for the plot
	numofReads   = sum(as.numeric(x*y))/(readlength-k+1) 
	estKmercov1  = x[which(y==max(y))][1]
	estCoverage1 = estKmercov1*readlength/(readlength-k)
	estLength1   = numofReads*readlength/estCoverage1

	nls1=nls_2peak(x,y,k,estKmercov1,estLength1,MAX_ITERATIONS)

	## Second we double the max kmercoverage
	estKmercov2 = estKmercov1*2
	estCoverage2= estKmercov2*readlength/(readlength-k)
	estLength2 = numofReads*readlength/estCoverage2

	nls2=nls_2peak(x,y,k,estKmercov2,estLength2,MAX_ITERATIONS)

	control=nls_control(x,y,k,estKmercov1,estLength1,MAX_ITERATIONS)
	return(eval_model(nls1,nls2))
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

    ## Features to report
    het=c(-1,-1)
    total_len=c(-1,-1)
    repeat_len=c(-1,-1)
    unique_len=c(-1,-1)
    dups=c(-1,-1)
    error_rate=c(-1,-1)
    model_status="fail"

    ## Plot the distribution, and hopefully with the model fit
	pdf(paste(foldername, "/plot.pdf", sep=""))
	plot(kmer_hist, type="h", main="GenomeScope profile\n", xlab="Coverage", ylab="Frequency", ylim=c(0,y_limit), xlim=c(0,x_limit))

    ## Make a second plot in log space over entire range
	pdf(paste(foldername, "/plot.log.pdf", sep=""))
	plot(kmer_hist, type="h", main="GenomeScope profile\n", xlab="Coverage", ylab="Frequency", log="xy")

	if(!is.null(container[[1]])) 
    { 
       ## The model converged!
       res<-data.frame(x,pred=predict(container[[1]], newdata=data.frame(x)))
       lines(x,res$pred,col="Blue",lwd=3)

       dev.set(dev.next())
       lines(x,res$pred,col="Blue",lwd=3)

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
       
       ## Compute error rate 
       error_xcutoff = floor(kcov[1])
       error_xcutoff_ind = which(x==error_xcutoff)

       error_kmers = y[1:error_xcutoff_ind] - res$pred[1:error_xcutoff_ind]
       error_kmers = pmax(error_kmers, 0)

       lines(x[1:error_xcutoff_ind], error_kmers, lwd=3, col="orange")
       dev.set(dev.next())
       lines(x[1:error_xcutoff_ind], error_kmers, lwd=3, col="orange")

       error_kmers = sum(as.numeric(error_kmers * x[1:error_xcutoff_ind]))
       total_kmers = sum(as.numeric(x*y))

       f1 <- function(x){
             i=seq(1,k) 
             h=(1-x)^(k-i)*x^i*choose(k,i)
             sum(h)*total_kmers-error_kmers
       }

       error_rate_root = try( uniroot(f1, c(0,1))$root)

       if (class(error_rate_root) == "try-error")
       {
         error_rate  = c(error_kmers/total_kmers/k, error_kmers/total_kmers/k)
       }
       else
       {
          error_rate  = c(error_rate_root, error_rate_root)
       }

       total_len = (total_kmers-error_kmers)/(2*kcov)
       
       ## find kmers that fit the 2 peak model (no repeats)
       unique_hist <- (2 * (1 - md[1]) * (1 - (1 - het[1])^k))                                * dnbinom(x, size = kcov[1]     / dups[1], mu = kcov[1])     * mlen[1] +
                      ((md[1] * (1 - (1 - het[1])^k)^2) + (1 - 2 * md[1]) * ((1 - het[1])^k)) * dnbinom(x, size = kcov[1] * 2 / dups[1], mu = kcov[1] * 2) * mlen[1]
       
       unique_kmers = sum(as.numeric(x*unique_hist))
       repeat_kmers = total_kmers - unique_kmers - error_kmers
       
       repeat_len=repeat_kmers/(2*kcov)
       unique_len=unique_kmers/(2*kcov)
       
       ## Add key values to plot
       title(paste("\nhet:", format(het[1], digits=3), 
                   " kcov:", format(kcov[1], digits=3), 
                   " err:", format(error_rate[1], digits=3), 
                   " dup:", format(dups[1], digits=3), 
                   " len:", round(total_len[1]), sep=""), 
                   cex.main=.85)
       
       ## Mark the modes of the peaks
       abline(v=kcov[1] * c(1,2,3,4), col="green", lty=2)
       
       ## Draw just the unique portion of the model
       lines(x, unique_hist, col="red", lty=3, lwd=3)

       dev.set(dev.next())

       ## update the other plot
       title(paste("\nhet:", format(het[1], digits=3), 
                   " kcov:", format(kcov[1], digits=3), 
                   " err:", format(error_rate[1], digits=3), 
                   " dup:", format(dups[1], digits=3), 
                   " len:", round(total_len[1]), sep=""), 
                   cex.main=.85)
       
       ## Mark the modes of the peaks
       abline(v=kcov[1] * c(1,2,3,4), col="green", lty=2)
       
       ## Draw just the unique portion of the model
       lines(x, unique_hist, col="red", lty=3, lwd=3)
       
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

	cat(paste("property", "min", "max", sep="\t"),                                      file=summaryFile, sep="\n")
	cat(paste("k", k, k, sep="\t"),                                                     file=summaryFile, sep="\n", append=TRUE) 
	cat(paste("Heterozygosity", het[1], het[2], sep="\t"),                              file=summaryFile, sep="\n", append=TRUE)
	cat(paste("GenomeHaploidLen", round(total_len[1]), round(total_len[2]), sep="\t"),  file=summaryFile, sep="\n", append=TRUE)
	cat(paste("GenomeRepeatLen", round(repeat_len[1]), round(repeat_len[2]), sep="\t"), file=summaryFile, sep="\n", append=TRUE)
	cat(paste("GenomeUniqueLen", round(unique_len[1]), round(unique_len[2]), sep="\t"), file=summaryFile, sep="\n", append=TRUE)
	cat(paste("ReadDuplicationLevel", dups[1], dups[2], sep="\t"),                      file=summaryFile, sep="\n", append=TRUE)
	cat(paste("ReadErrorRate", error_rate[1], error_rate[2], sep="\t"),                 file=summaryFile, sep="\n", append=TRUE)
	cat(paste("ModelScore", container[2], container[2], sep="\t"),                      file=summaryFile, sep="\n", append=TRUE)

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

    if ((length(args) == 5) && (as.numeric(args[[5]] == 1))) { DEBUG = 1 }

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

    ## terminate after num iterations or if the score becomes good enough, store best result so far in container
	num <- 0
	best_container <- list(NULL,0)

	while(num < NUM_ROUNDS) 
    {
        cat(paste("round", num, "trimming to", start, "trying 4peak model... "), file=progressFilename, sep="", append=TRUE)
        if (DEBUG) { cat(paste("round", num, "trimming to", start, "trying 4peak model... \n")) }

        ## Reset the input trimming off low frequency error kmers
        max <- length(kmer_prof[,1])
        x <- kmer_prof[start:max,1]
        y <- kmer_prof[start:max,2]

		model_4peaks <- estimate_Genome_4peak2(kmer_prof_orig, x, y, k, readlength, foldername) 

        if (!is.null(model_4peaks[[1]])) { 
          cat(paste("converged. score: ", model_4peaks[[2]][[1]]), file=progressFilename, sep="\n", append=TRUE)
        } else {
          cat(paste("unconverged"), file=progressFilename, sep="\n", append=TRUE)
        }

		#check if this result is better than previous
        if (!is.null(model_4peaks[[1]]))
        {
          if ((is.null(best_container[[1]])) || (model_4peaks[[2]][[1]] < best_container[[2]][[1]])){
              best_container = model_4peaks
          }
        }

        ## If we didnt get a good score, try again with random noise added
        start <- start + START_SHIFT
		kmer_prof[,2] <- kmer_prof[,2] + rnorm(length(kmer_prof[,2]), sd = RAND_ERROR) 
		num <- num + 1
	}

    ## Report the results, note using the original full profile
	report_results(kmer_prof_orig, k, best_container, foldername)
}
