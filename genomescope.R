#!/bin/env Rscript

## GenomeScope: Fast Genome Analysis from Unassembled Short Reads
## This is the automated script for reading in a histogram file given the k-mer size, readlength and an optional parameter for what you want the x-axis limit to be.

## Threshold on score for stopping search
GOOD_SCORE_THRESHOLD=90

## Number of rounds before giving up
NUM_ROUNDS=5

## Coverage steps to trim off between rounds
START_SHIFT=10

## Amount of random noise to introduce
RAND_ERROR = 0.01

## Max rounds on NLS
MAX_ITERATIONS=50

## Print out debugging messages (0/1)
DEBUG = 1



## Helper Function
###############################################################################

min_max <- function(table){
	return (c( abs(table[1]) - abs(table[2]) , abs(table[1])+abs(table[2])))
}



## Use nls to fit 4 peak model
###############################################################################

nls_4peak<-function(x, y, k, estKmercov, estLength, max_iterations){
	model4 = NULL

    if (DEBUG) { print("trying standard algorithm") }
	try(model4 <- nls(y ~ ((2.0 * (1-d) * (1-(1-r)^k))*dnbinom(x,size=kmercov/(bias),mu=kmercov)*length +
                          ((d*(1-(1-r)^k)^2) + (1-2*d)*((1-r)^k))*dnbinom(x,size=kmercov*2/(bias),mu=kmercov*2)*length + 
                          (2*d*((1-r)^k)*(1-(1-r)^k))*dnbinom(x,size = kmercov*3/(bias), mu = kmercov*3)*length + 
                          (d*(1-r)^(2*k))*dnbinom(x, size = kmercov*4 / (bias), mu = kmercov*4)*length), 
                      start=list(d=0, r=0, kmercov=estKmercov, bias = 0.5, length=estLength),
                      control=list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

	if(class(model4) == "try-error"){
        if (DEBUG) { print("retrying with port algorithm") }
		try(model4 <- nls(y ~ ((2.0 * (1-d) * (1-(1-r)^k))*dnbinom(x,size=kmercov/(bias),mu=kmercov)*length + 
                              ((d*(1-(1-r)^k)^2) + (1-2*d)*((1-r)^k))*dnbinom(x,size=kmercov*2/(bias),mu=kmercov*2)*length + 
                              (2*d*((1-r)^k)*(1-(1-r)^k))*dnbinom(x,size = kmercov*3/(bias), mu = kmercov*3)*length + 
                              (d*(1-r)^(2*k))*dnbinom(x, size = kmercov*4 / (bias), mu = kmercov*4)*length), 
                          start=list(d=0,r=0,kmercov=estKmercov,bias = 0.5,length=estLength),
                          algorithm="port", control=list(minFactor=1e-12,maxiter=max_iterations)), silent = TRUE)
	}

	return(model4)
}


## Use nls to fit 2 peak model
###############################################################################

nls_2peak<-function(x, y, k, estKmercov, estLength, max_iterations){	
	model2 = NULL

    if (DEBUG) { print("trying standard algorithm") }
	try(model2 <- nls(y ~ 2*(1-exp(-r*k))*dnbinom(x,size=(kmercov/2)/(bias+1),mu=kmercov/2)*length + 
                          exp(-r*k)*dnbinom(x,size=kmercov/(bias+1),mu=kmercov)*length,
                      start=list(length=estLength, r=0.0, kmercov=estKmercov, bias=1),
                      control=list(minFactor=1e-12,maxiter=max_iterations)))

	if(class(model2) == "try-error"){
        if (DEBUG) { print("retrying with port algorithm") }
		try(model2 <- nls(y ~ 2*(1-exp(-r*k))*dnbinom(x,size=(kmercov/2)/(bias+1),mu=kmercov/2)*length + 
                              exp(-r*k)*dnbinom(x,size=kmercov/(bias+1),mu=kmercov)*length,
                          start=list(length=estLength,r=0.0,kmercov=estKmercov,bias=1),
                          algorithm="port", control=list(minFactor=1e-12,maxiter=max_iterations)))
	}

	return(model2)
}


## Use nls to fit 1 peak model as control
###############################################################################

nls_control<-function(x, y, k, estKmercov, estLength, max_iterations){
    control = NULL

    if (DEBUG) { print("trying standard algorithm") }
	control <- try(nls(y ~ (dnbinom(x,size=kmercov/(bias),mu=kmercov)*length), 
                            start=list(kmercov=estKmercov,bias = 0.5,length=estLength),
                            control=list(minFactor=1e-12,maxiter=max_iterations)), silent = TRUE)

	if(class(control) == "try-error"){
        if (DEBUG) { print("retrying with port algorithm") }
		control <- try(nls(y ~ (dnbinom(x,size=kmercov/(bias),mu=kmercov)*length),
                           start=list(d=0,r=0,kmercov=estKmercov,bias = 0.5,length=estLength),
                           algorithm="port",control=list(minFactor=1e-12,maxiter=max_iterations)), silent = TRUE)
	}

	return(control)
}


## Pick between the two model forms, resolves ambiguity between which is the homozygous and which is the heterozygous peak
###############################################################################

eval_model<-function(nls1, nls2, control){
 	if(!is.null(nls1[[1]]) ){##&& (is.null(nls2[[1]]) || nls1$m$deviance() > nls2$m$deviance())){
		return (list(nls1,(control$m$deviance()/nls1$m$deviance())))
	}else if (!is.null(nls2[[1]])){ 
		return (list(nls2,(control$m$deviance()/nls2$m$deviance())))
	}else{	#No Model converged
		return (list(NULL,0))
	}
}


## Wrapper function to try fitting 4 peak model with 2 forms
###############################################################################

estimate_Genome_4peak2<-function(x, y, k, readlength, foldername){
	## First we see what happens when the max peak is the kmercoverage (typically the homozygous peak) for the plot
	numofReads   = sum(as.numeric(x*y))/(readlength-k+1) 
	estKmercov1  = x[which(y==max(y))][1]
	estCoverage1 = estKmercov1*readlength/(readlength-k)
	estLength1   = numofReads*readlength/estCoverage1

    if (DEBUG) { print(paste("trying with kmercov: ", estKmercov1)) }
	nls1    = nls_4peak(x, y, k, estKmercov1, estLength1, MAX_ITERATIONS)
    if (DEBUG) { print(summary(nls1)) }

	## Second we half the max kmercoverage (typically the heterozygous peak)
	estKmercov2  = estKmercov1 / 2 ##2.5
	estCoverage2 = estKmercov2*readlength/(readlength-k)
	estLength2   = numofReads*readlength/estCoverage2

    if (DEBUG) { print(paste("trying with kmercov: ", estKmercov2)) }
	nls2 = nls_4peak(x, y, k, estKmercov2, estLength2, MAX_ITERATIONS)
    if (DEBUG) { print(summary(nls2)) }

    ## Now try a simple control model
    if (DEBUG) { print(paste("trying control with kmercov: ", estKmercov1)) }
	control = nls_control(x, y, k, estKmercov1, estLength1, MAX_ITERATIONS)
    if (DEBUG) { print(summary(control)) }

	return(eval_model(nls1, nls2, control))
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
	return(eval_model(nls1,nls2,control))
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
	start=which(y == min(y[1:15]))
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
       unique_hist <- (2 * (1 - md[1]) * (1 - (1 - het[1])^k)) * dnbinom(x, size = kcov[1]/dups[1], mu = kcov[1]) * mlen[1] +
                      ((md[1] * (1 - (1 - het[1])^k)^2) + (1 - 2 * md[1]) * ((1 - het[1])^k)) * dnbinom(x, size = kcov[1] * 2/(dups[1]), mu = kcov[1] * 2) * mlen[1]
       
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
	}
    else
    {
      title("\nFailed to converge")
      dev.set(dev.next())
      title("\nFailed to converge")
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
	cat("USAGE: genomescope.R histogram_file k-mer_length read_length output_dir\n")
} else{
    #	options(show.error.messages = FALSE)
    #	options(error=FALSE)
    #	options(showWarnCalls=FALSE)
    #	options(showNCalls=FALSE)
	
	histfile   <- args[[1]]
	k          <- as.numeric(args[[2]])
	readlength <- as.numeric(args[[3]])
	title      <- args[[4]] #We dont really need this.

    print(paste("running on", histfile))

    ## values for testing
    #histfile <- "~/build/genomescope/simulation/simulation_results/Arabidopsis_thaliana.TAIR10.26.dna_sm.toplevel.fa_het0.01_br1_rl100_cov100_err0.01_reads.fa21.hist"
    #k <- 21
    #readlength <- 100
    #title <- "~/build/genomescope/simulation/simulation_analysis/Arabidopsis_thaliana.TAIR10.26.dna_sm.toplevel.fa_het0.01_br1_rl100_cov100_err0.01_reads.fa21.hist"

	foldername = paste(title, "_results", sep="")
	dir.create(foldername, showWarnings=FALSE)

	kmer_prof <- read.csv(file=histfile,sep=" ", header=FALSE) 
	kmer_prof <- kmer_prof[c(1:length(kmer_prof[,2])-1),] #get rid of the last position
    kmer_prof_orig <- kmer_prof

    ## Initialize the status
    progressFilename <- paste(foldername,"/progress.txt",sep="")
	cat("starting", file=progressFilename, sep="\n")

    ## try to find the local minimum between errors and the first (heterozygous) peak
    start <- which(kmer_prof[,2]==min(kmer_prof[1:15,2]))

    ## terminate after num iterations or if the score becomes good enough, store best result so far in container
	num <- 0
	container <- c(-1,-1)

	while(num < NUM_ROUNDS) 
    {
        ## Reset the input trimming off low frequency error kmers
        max <- length(kmer_prof[,1])
        x <- kmer_prof[start:max,1]
        y <- kmer_prof[start:max,2]

        cat(paste("round", num, "trimming to", start, "trying 4peak model... "), file=progressFilename, sep="", append=TRUE)
        print(paste("round", num, "trimming to", start, "trying 4peak model... "))

		model_4peaks <- estimate_Genome_4peak2(x,y,k, readlength, foldername) 

        if (!is.null(container[[1]])) { 
          cat(paste("converged. score: ", container[[2]][[1]]), file=progressFilename, sep="\n", append=TRUE)
        } else {
          cat(paste("unconverged"), file=progressFilename, sep="\n", append=TRUE)
        }

        #cat(paste("trying 2peak model, round", num), file=progressFilename, sep="\n", append=TRUE)
		#model_2peaks=estimate_Genome_2peak2(d,x,y,k, readlength, foldername) 

		#check if this result is better than previous
		if ((num==0) || (model_4peaks[2][[1]]>container[2][[1]])){
			container = model_4peaks
		}

	    #if(model_2peaks[2][[1]]>container[2][[1]]){
	    #	container=model_2peaks
	    #	best_p=p
	    #}

		if (container[2][[1]] > GOOD_SCORE_THRESHOLD)
        { 
            ## Score is good enough, stop trying
			break
		}

        ## If we didnt get a good score, try again with random noise added
        start <- start + START_SHIFT
		kmer_prof[,2] <- kmer_prof[,2] + rnorm(length(kmer_prof[,2]), sd = RAND_ERROR) 
		num <- num+1
	}

    ## Report the results, note using the original full profile
	report_results(kmer_prof_orig, k, container, foldername)
}
