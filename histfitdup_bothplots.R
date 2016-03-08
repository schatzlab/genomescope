#!/bin/env Rscript

#Gregory Vurture + Fritz J Sedlazeck
#histfitdup.R
#This is the automated script for reading in a histogram file given the k-mer size, readlength and an optional parameter for what you want the x-axis limit to be.

min_max <- function(table){
	return (c( abs(table[1]) - abs(table[2]) , abs(table[1])+abs(table[2])))
}

nls_4peak<-function(d,x,y,k,estKmercov,estLength,max_iterations){
	nlspoisson=NULL
	try(nlspoisson<-nls(y~((2.0 * (1-d) * (1-(1-r)^k))*dnbinom(x,size=kmercov/(bias),mu=kmercov)*length+( (d*(1-(1-r)^k)^2) + (1-2*d)*((1-r)^k))*dnbinom(x,size=kmercov*2/(bias),mu=kmercov*2)*length+(2*d*((1-r)^k)*(1-(1-r)^k))*dnbinom(x,size = kmercov*3/(bias), mu = kmercov*3)*length+(d*(1-r)^(2*k))*dnbinom(x, size = kmercov*4 / (bias), mu = kmercov*4)*length),start=list(d=0,r=0,kmercov=estKmercov,bias = 0.5,length=estLength),control=list(minFactor=1e-12,maxiter=max_iterations)), silent = TRUE)
	if(class(nlspoisson) == "try-error"){
		try(nlspoisson<-nls(y~((2.0 * (1-d) * (1-(1-r)^k))*dnbinom(x,size=kmercov/(bias),mu=kmercov)*length+( (d*(1-(1-r)^k)^2) + (1-2*d)*((1-r)^k))*dnbinom(x,size=kmercov*2/(bias),mu=kmercov*2)*length+(2*d*((1-r)^k)*(1-(1-r)^k))*dnbinom(x,size = kmercov*3/(bias), mu = kmercov*3)*length+(d*(1-r)^(2*k))*dnbinom(x, size = kmercov*4 / (bias), mu = kmercov*4)*length),start=list(d=0,r=0,kmercov=estKmercov,bias = 0.5,length=estLength),algorithm="port",control=list(minFactor=1e-12,maxiter=max_iterations)), silent = TRUE)
	}
	return(nlspoisson)
}

nls_2peak<-function(d,x,y,k,estKmercov,estLength,max_iterations){	
	nlspoisson=NULL
	try({nlspoisson<-nls(y~2*(1-exp(-r*k))*dnbinom(x,size=(kmercov/2)/(bias+1),mu=kmercov/2)*length+exp(-r*k)*dnbinom(x,size=kmercov/(bias+1),mu=kmercov)*length,start=list(length=estLength,r=0.0,kmercov=estKmercov,bias=1),control=list(minFactor=1e-12,maxiter=max_iterations))})
	if(class(nlspoisson) == "try-error"){
		try({nlspoisson<-nls(y~2*(1-exp(-r*k))*dnbinom(x,size=(kmercov/2)/(bias+1),mu=kmercov/2)*length+exp(-r*k)*dnbinom(x,size=kmercov/(bias+1),mu=kmercov)*length,start=list(length=estLength,r=0.0,kmercov=estKmercov,bias=1),algorithm="port",control=list(minFactor=1e-12,maxiter=max_iterations))})
	}
	return(nlspoisson)
}

nls_control<-function(d,x,y,k,estKmercov,estLength,max_iterations){
	control<-try(nls(y~(dnbinom(x,size=kmercov/(bias),mu=kmercov)*length),start=list(kmercov=estKmercov,bias = 0.5,length=estLength),control=list(minFactor=1e-12,maxiter=max_iterations)), silent = TRUE)
	if(class(control) == "try-error"){
		control<-try(nls(y~(dnbinom(x,size=kmercov/(bias),mu=kmercov)*length),start=list(d=0,r=0,kmercov=estKmercov,bias = 0.5,length=estLength),algorithm="port",control=list(minFactor=1e-12,maxiter=max_iterations)), silent = TRUE)
	}
	return(control)
}

eval_model<-function(nls1,nls2,control){
 	if(!is.null(nls1[[1]]) ){##&& (is.null(nls2[[1]]) || nls1$m$deviance() > nls2$m$deviance())){
		return (list(nls1,(control$m$deviance()/nls1$m$deviance())))
	}else if (!is.null(nls2[[1]])){ 
		return (list(nls2,(control$m$deviance()/nls2$m$deviance())))
	}else{	#No Model converged
		return (list(NULL,0))
	}
}

estimate_Genome_4peak2<-function(d,x,y,k,readlength,foldername){
	max_iterations=100
	#We can calculate the numer of reads from the number of kmers, read-length, and kmersize
	#First we see what happens when the max peak is the kmercoverage for the plot
	numofReads = sum(as.numeric(x*y))/(readlength-k+1) 
	estKmercov1 = x[which(y==max(y))][1]
	estCoverage1= estKmercov1*readlength/(readlength-k)
	estLength1 = numofReads*readlength/estCoverage1
	
	nls1=nls_4peak(d,x,y,k,estKmercov1,estLength1,max_iterations)
	control=nls_control(d,x,y,k,estKmercov1,estLength1,max_iterations)
	#Second we half the max kmercoverage
	estKmercov2 = estKmercov1/2 ##2.5
	estCoverage2= estKmercov2*readlength/(readlength-k)
	estLength2 = numofReads*readlength/estCoverage2
	nls2=nls_4peak(d,x,y,k,estKmercov2,estLength2,max_iterations)
	return(eval_model(nls1,nls2,control))
}

estimate_Genome_2peak2<-function(d,x,y,k,readlength,foldername){
	max_iterations=100
	numofReads = sum(as.numeric(x*y))/(readlength-k+1) 
	estKmercov1 = x[which(y==max(y))][1]
	estCoverage1= estKmercov1*readlength/(readlength-k)
	estLength1 = numofReads*readlength/estCoverage1
	nls1=nls_2peak(d,x,y,k,estKmercov1,estLength1,max_iterations)
	control=nls_control(d,x,y,k,estKmercov1,estLength1,max_iterations)
	#Second we double the max kmercoverage
	estKmercov2 = estKmercov1*2
	estCoverage2= estKmercov2*readlength/(readlength-k)
	estLength2 = numofReads*readlength/estCoverage2
	nls2=nls_2peak(d,x,y,k,estKmercov2,estLength2,max_iterations)
	return(eval_model(nls1,nls2,control))
}

report_results<-function(p,container,foldername)
{
        x=p[[1]]
        y=p[[2]]
	d=data.frame(x=p[[1]],y=p[[2]])

	#automatically zoom into the relevant regions of the plot
	start=which(p[,2]==min(p[1:15,2]))
	zoomx=d$x[start:xmax-1]
	zoomy=d$y[start:xmax-1]
	
	y_limit=max(zoomy[start:length(zoomy)])*1.1
	
	x_limit= which(max(p[c(start:length(zoomx)),2])==p[,2]) * 3

	if(min(zoomy)>zoomy[1]){
		x_limit=max(which(zoomy<zoomy[1])[2],600)
	}

        ## Features to report
        het=c(-1,-1)
        size=c(-1,-1)
        dups=c(-1,-1)
        repeats=c(-1,-1)
        error_rate=c(-1,-1)
        model_status="fail"

        ## Plot the distribution, and hopefully with the model fit
	pdf(paste(foldername,"/plot.pdf",sep=""))
	plot(p,type="h",main="GenomeScope profile\n",xlab="Coverage",ylab="Frequency",ylim=c(0,y_limit),xlim=c(0,x_limit))

	if(!is.null(container[[1]])) { 
                ## The model converged!
	        res<-data.frame(x,pred=predict(container[[1]], newdata=data.frame(x)))
		lines(x,res$pred,col="Blue",lwd=3)
	
                ## Compute the genome characteristics
                model_sum=summary(container[[1]])

                het  = min_max(model_sum$coefficients['r',])
                dups = min_max(model_sum$coefficients['bias',])
                kcov = min_max(model_sum$coefficients['kmercov',])

                abline(v=kcov[1] * c(1,2,3,4), col="green", lty=2)

                error_xcutoff = floor(kcov[1])
                error_xcutoff_ind = which(x==error_xcutoff)

                error_kmers = sum(y[1:error_xcutoff_ind] - res$pred[1:error_xcutoff_ind])
                total_kmers = sum(as.numeric(x*y))

                error_rate  = c(error_kmers/total_kmers/k, error_kmers/total_kmers/k)
                size = (total_kmers-error_kmers)/(2*kcov)

                repeats=c(-1,-1)
                
                title(paste("\nhet:", format(het[1], digits=3), 
                            " kcov:", format(kcov[1], digits=3), 
                            " err:", format(error_rate[1], digits=3), 
                            " dup:", format(dups[1], digits=3), 
                            " len:", round(size[1]), sep=""))

                ## save the model to a file
                capture.output(model_sum, file=paste(foldername,"/model.txt", sep=""))

                model_status="done"
	}

	dev.off()

	summaryFile <- paste(foldername,"/summary.txt",sep="")
      
	cat(paste("property", "min", "max", sep="\t"),                                  file=summaryFile, sep="\n")
	cat(paste("k", k, k, sep="\t"),                                                 file=summaryFile, sep="\n", append=TRUE) 
	cat(paste("Heterozygosity(%)", het[1], het[2], sep="\t"),                       file=summaryFile, sep="\n", append=TRUE)
	cat(paste("Haploid Genome Size(bp)", round(size[1]), round(size[2]), sep="\t"), file=summaryFile, sep="\n", append=TRUE)
	cat(paste("Read Duplication Level(times)", dups[1], dups[2], sep="\t"),         file=summaryFile, sep="\n", append=TRUE)
	cat(paste("Repeats", repeats[1], repeats[2], sep="\t"),                         file=summaryFile, sep="\n", append=TRUE)
	cat(paste("Ratio", container[2], container[2], sep="\t"),                               file=summaryFile, sep="\n", append=TRUE)

        ## Finalize the progress
        fileConn<-file(paste(foldername,"/progress.txt",sep=""),open="w")
        writeLines(model_status, fileConn)
        close(fileConn)
}

args<-commandArgs(TRUE)

if(length(args) < 4)
{
	cat("USAGE: histfit.R histogram_file k-mer_length read_length title_of_dataset [xlim max value (optional)]\n")
} else{
#	options(show.error.messages = FALSE)
#	options(error=FALSE)
#	options(showWarnCalls=FALSE)
#	options(showNCalls=FALSE)
	
	histfile <- args[[1]]
	k<-as.numeric(args[[2]])
	readlength <-as.numeric(args[[3]])
	title <- args[[4]] #We dont really need this.

   #     histfile <- "~/build/genomescope/simulation/simulation_results/Arabidopsis_thaliana.TAIR10.26.dna_sm.toplevel.fa_het0.01_br1_rl100_cov100_err0.01_reads.fa21.hist"
   #     k <- 21
   #     readlength <- 100
   #     title <- "~/build/genomescope/simulation/simulation_analysis/Arabidopsis_thaliana.TAIR10.26.dna_sm.toplevel.fa_het0.01_br1_rl100_cov100_err0.01_reads.fa21.hist"

	foldername = paste(title,"_results",sep="")
	dir.create(foldername,showWarnings=FALSE)

	kmer_prof=read.csv(file=histfile,sep=" ", header=FALSE) 
	kmer_prof=kmer_prof[c(1:length(kmer_prof[,2])-1),] #get rid of the last position
	
	num=0
	container=c(-1,-1)
	best_p=kmer_prof
	p=kmer_prof

	while(num < 5) { ## terminate after num iterations or if the score becomes good enough

		d=data.frame(x=p[[1]],y=p[[2]])
		if(length(args)==5){
			xmax <- as.numeric(args[[5]])
		} else{
			xmax <- nrow(p)
		}
	
		start=which(p[,2]==min(p[1:15,2]))
		x=d$x[start:xmax-1]
		y=d$y[start:xmax-1]

		#estimates the model
		model_4peaks=estimate_Genome_4peak2(d,x,y,k, readlength, foldername) 
		model_2peaks=estimate_Genome_2peak2(d,x,y,k, readlength, foldername) 

		#check wich model works better:
		if(num==0 || model_4peaks[2][[1]]>container[2][[1]]){
			container=model_4peaks
			best_p=p
		}
		if(model_2peaks[2][[1]]>container[2][[1]]){
			container=model_2peaks
			best_p=p
		}
		if(container[2][[1]]>90){ #it does not make sense to continue as the result is good enough. 
			break
		}

                ## If we didnt get a good score, try again with random noise added
		p[,2]= p[,2] + rnorm(length(p[,2]), sd = 0.01) 
		num = num+1
	}

	report_results(kmer_prof,container,foldername)
}
