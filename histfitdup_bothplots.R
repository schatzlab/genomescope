

#Gregory Vurture + Fritz J Sedlazeck + Maria Nattestad
#!/bin/env Rscript
#histfitdup.R
#This is the automated script for reading in a histogram file given the k-mer size, readlength and an optional parameter for what you want the x-axis limit to be.

# install.packages("ggplot2", repos='http://cran.us.r-project.org')
# install.packages("graphics", repos='http://cran.us.r-project.org')
# install.packages("stats", repos='http://cran.us.r-project.org')
# install.packages("grDevices", repos='http://cran.us.r-project.org')

# source("http://bioconductor.org/biocLite.R")
# biocLite("rhdf5")


library(ggplot2)


min_max <- function(table){
	return (c( abs(table[1]) - abs(table[2]) , abs(table[1])+abs(table[2])))
}


estimate_Genome<-function(d,x,y,k,readlength,foldername){
max_iterations=100
tresh_control=1
	#We can calculate the numer of reads from the number of kmers, read-length, and kmersize
	#First we see what happens when the max peak is the kmercoverage for the plot
	numofReads = sum(as.numeric(x*y))/(readlength-k+1) 
	estKmercov1 = x[which(y==max(y))][1]
	estCoverage1= estKmercov1*readlength/(readlength-k)
	estLength1 = numofReads*readlength/estCoverage1
	nlspoisson1=NULL
	nlspoisson2=NULL
	nlspoisson3=NULL
	

	nlspoisson1<-try(nls(y~((2.0 * (1-d) * (1-(1-r)^k))*dnbinom(x,size=kmercov/(bias),mu=kmercov)*length+( (d*(1-(1-r)^k)^2) + (1-2*d)*((1-r)^k))*dnbinom(x,size=kmercov*2/(bias),mu=kmercov*2)*length+(2*d*((1-r)^k)*(1-(1-r)^k))*dnbinom(x,size = kmercov*3/(bias), mu = kmercov*3)*length+(d*(1-r)^(2*k))*dnbinom(x, size = kmercov*4 / (bias), mu = kmercov*4)*length),start=list(d=0,r=0,kmercov=estKmercov1,bias = 0.5,length=estLength1),control=list(minFactor=1e-12,maxiter=max_iterations)), silent = FALSE)
	if(class(nlspoisson1) == "try-error"){
		nlspoisson1<-try(nls(y~((2.0 * (1-d) * (1-(1-r)^k))*dnbinom(x,size=kmercov/(bias),mu=kmercov)*length+( (d*(1-(1-r)^k)^2) + (1-2*d)*((1-r)^k))*dnbinom(x,size=kmercov*2/(bias),mu=kmercov*2)*length+(2*d*((1-r)^k)*(1-(1-r)^k))*dnbinom(x,size = kmercov*3/(bias), mu = kmercov*3)*length+(d*(1-r)^(2*k))*dnbinom(x, size = kmercov*4 / (bias), mu = kmercov*4)*length),start=list(d=0,r=0,kmercov=estKmercov1,bias = 0.5,length=estLength1),algorithm="port",control=list(minFactor=1e-12,maxiter=max_iterations)), silent = FALSE)
	}
	
	control_1<-try(nls(y~(dnbinom(x,size=kmercov/(bias),mu=kmercov)*length),start=list(kmercov=estKmercov1,bias = 0.5,length=estLength1),control=list(minFactor=1e-12,maxiter=max_iterations)), silent = FALSE)
	if(class(control_1) == "try-error"){
		control_1<-try(nls(y~(dnbinom(x,size=kmercov/(bias),mu=kmercov)*length),start=list(d=0,r=0,kmercov=estKmercov1,bias = 0.5,length=estLength1),algorithm="port",control=list(minFactor=1e-12,maxiter=max_iterations)), silent = FALSE)
	}
	
	#Second we half the max kmercoverage
	estKmercov2 = estKmercov1/2.5
	estCoverage2= estKmercov2*readlength/(readlength-k)
	estLength2 = numofReads*readlength/estCoverage2
	nlspoisson2<-try(nls(y~((2.0 * (1-d) * (1-(1-r)^k))*dnbinom(x,size=kmercov/(bias),mu=kmercov)*length+( (d*(1-(1-r)^k)^2) + (1-2*d)*((1-r)^k))*dnbinom(x,size=kmercov*2/(bias),mu=kmercov*2)*length+(2*d*((1-r)^k)*(1-(1-r)^k))*dnbinom(x,size = kmercov*3/(bias), mu = kmercov*3)*length+(d*(1-r)^(2*k))*dnbinom(x, size = kmercov*4 / (bias), mu = kmercov*4)*length),start=list(d=0,r=0,kmercov=estKmercov2,bias = 0.5,length=estLength2),control=list(minFactor=1e-12,maxiter=max_iterations)), silent = FALSE)
	if(class(nlspoisson2) == "try-error"){
		nlspoisson2<-try(nls(y~((2.0 * (1-d) * (1-(1-r)^k))*dnbinom(x,size=kmercov/(bias),mu=kmercov)*length+( (d*(1-(1-r)^k)^2) + (1-2*d)*((1-r)^k))*dnbinom(x,size=kmercov*2/(bias),mu=kmercov*2)*length+(2*d*((1-r)^k)*(1-(1-r)^k))*dnbinom(x,size = kmercov*3/(bias), mu = kmercov*3)*length+(d*(1-r)^(2*k))*dnbinom(x, size = kmercov*4 / (bias), mu = kmercov*4)*length),start=list(d=0,r=0,kmercov=estKmercov2,bias = 0.5,length=estLength2),algorithm="port",control=list(minFactor=1e-12,maxiter=max_iterations)), silent = FALSE)
	}
	
    if(class(nlspoisson1) != "try-error"){ # && (control_1$m$deviance()/nlspoisson1$m$deviance())>tresh_control){ 								#Model 1 converged
		return (list(nlspoisson1,(control_1$m$deviance()/nlspoisson1$m$deviance())))
	}else if (class(nlspoisson2) != "try-error"){ # && (control_1$m$deviance()/nlspoisson2$m$deviance())>tresh_control) { 						#Model 2 converged
		return (list(nlspoisson2,(control_1$m$deviance()/nlspoisson2$m$deviance())))
	}else{
		fileConn<-file(paste(foldername,"/","progress.txt",sep=""),open="w")
		if(class(nlspoisson1) != "try-error" && (control_1$m$deviance()/nlspoisson1$m$deviance())<tresh_control){
			writeLines(paste("Error: Control",(control_1$m$deviance()/nlspoisson1$m$deviance()),sep='\t'), fileConn)
		}else if (class(nlspoisson2) != "try-error" && (control_1$m$deviance()/nlspoisson2$m$deviance())<tresh_control) { 						#Model 2 converged
			writeLines(paste("Error: Control",(control_1$m$deviance()/nlspoisson2$m$deviance()),sep='\t'), fileConn)
		}else {
			#No Model converged
			writeLines("Error: Model did not converge.", fileConn)
		}
		close(fileConn)
		return (lsit(NULL,control_1))
	}
}

gen_summary<-function(model){
	size=min_max(summary(model)$coefficients[5,])
	het=min_max(summary(model)$coefficients[2,])
	repeats=min_max(summary(model)$coefficients[1,])
	dups=min_max(summary(model)$coefficients[4,])
	fileConn_txt<-file(paste(foldername,"/","summary.txt",sep=""))
	writeLines(paste("property\tmin\tmax\n","k\t",k,"\nHeterozygosity(%)\t",het[1],"\t",het[2],"\nHaploid Genome Size(bp)\t",round(size[1]),"\t",round(size[2]),"\nRead Duplication Level(times)\t",dups[1],"\t",dups[2]," \nRepeats\t",repeats[1],"\t",repeats[2],"\n",sep=""), fileConn_txt)
	close(fileConn_txt)

	fileConn_html<-file(paste(foldername,"/","summary.html",sep=""))
	# writeLines(paste('<table style="font-family:sans-serif">','<tr class="active"><th>k</th><th>Heterozygosity</th><th>Haploid genome size (bp)</th><th>Read duplication level (X)</th><th>Repeats</th></tr>', '<tr><td>',k,'</td><td>',signif(het[1],4),'</td><td>',round(size[1]),'</td><td>',signif(dups[1]),'</td><td>',signif(repeats[1],4),'</td></tr>','<tr><td></td><td>',signif(het[2],4),'</td><td>',round(size[2]),'</td><td>',signif(dups[2],4),'</td><td>',signif(repeats[2],4),'</td></tr>' ,'</table>',sep=""),fileConn_html)
	writeLines(paste('<style> table, th, td {border: 5px solid #B9E2F0; } td { padding: 15px;} th {padding: 15px;} table {border-collapse: collapse; text-align: left;} </style><table style="font-family:sans-serif;" border:"black">','<tr><th>k</th><td>',k,'</td></tr><tr><th>Heterozygosity</th><td>',signif(het[1],4)*100,'% - ',signif(het[2],4)*100,'%</td></tr><tr><th>Haploid genome size (bp)</th><td>',formatC(round(size[1]),format="d",big.mark=","),' - ',formatC(round(size[2]),format="d",big.mark=","),'</td></tr><tr><th>Read duplication level</th><td>',1+signif(dups[1],4),'X - ',1+signif(dups[2],4),'X</td></tr><tr><th>Repeats</th><td>',signif(repeats[1],4)*100,'% - ',signif(repeats[2],4)*100,'%</td></tr>' ,'</table>',sep=""),fileConn_html)
	close(fileConn_html)
}



args<-commandArgs(TRUE)

if(length(args)==0)
{
	cat("USAGE: histfit.R histogram_file [k-mer length] [read length] [title of dataset] [xlim max value (optional)]\n")
} else{
	options(show.error.messages = FALSE)
	histfile <- args[[1]]
	k<-as.numeric(args[[2]])
	readlength <-as.numeric(args[[3]])
	title <- args[[4]] 

	############################################################
	# FOR TESTING:
# 	histfile <- "/Applications/XAMPP/htdocs/kmers/user_uploads/YV2HXhqNt2ve6ga2AtZT"
# 	k<-21
# 	readlength <- 100
# 	title <- "/Applications/XAMPP/htdocs/kmers/user_data/Maria_4"
	############################################################

	foldername = title
	dir.create(foldername,showWarnings=FALSE)
	p=read.csv(file=histfile,sep=" ") 
	p[,2]= p[,2]+rnorm(length(p[,2]), sd = 0.01) #edited FS
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
	container=estimate_Genome(d,x,y,k, readlength, foldername) 
	model = container[[1]]
	#Maria: thanks for the plot:
    data_for_plot <- data.frame(Coverage=d$x,Frequency=d$y)
    head(data_for_plot)
    # x = coverage
    # y = frequency


    
	theme_set(theme_gray(base_size = 50))
    
    bar_colour <- "#50BADE"
    prediction_line_colour <- "#000000"
	cols <- c("Frequency"=bar_colour,"Predicted_frequency"=prediction_line_colour)
    
	if(!is.null(model)){ #check if prediction worked
        # PLOT with prediction line:
        
        prediction <- data.frame(Coverage=x,Predicted_frequency=predict(model))

        data_for_plot <- merge(data_for_plot,prediction,all=TRUE)
        png(paste(foldername,"/","plot.png",sep=""),width=1800,height=1200)
       
        #automatically choose the ylim:
        y_lim=max(data_for_plot$Frequency[start:xmax-1],data_for_plot$Predicted_frequency[start:xmax-1])
        y_lim=y_lim*1.1
        x_lim=max(which(data_for_plot$Predicted_frequency>10))
        x_lim=x_lim*1.1
        #  print() is necessary when this is run from the command-line, otherwise it doesn't produce a plot
        print(ggplot(data_for_plot,aes(x=Coverage,y=Frequency)) + geom_bar(mapping=aes(fill="Frequency"),stat="identity",colour=bar_colour) + geom_line(aes(x=Coverage,y=Predicted_frequency,colour="Predicted_frequency"),size=1.5) + scale_fill_manual(name="Data",values=cols) + scale_colour_manual(name="Prediction",values=cols)+ coord_cartesian(ylim = c(0, y_lim),xlim=c(0,x_lim)))
        dev.off()

		#generates the summary report
		gen_summary(model)
	
		fileConn<-file(paste(foldername,"/","progress.txt",sep=""))
		writeLines("done", fileConn)
		close(fileConn)
	} else { # If prediction did not work, draw histogram without the prediction line
        
	    # PLOT without prediction line:
	    png(paste(foldername,"/","plot.png",sep=""),width=1800,height=1200)
	    #  print() is necessary when this is run from the command-line, otherwise it doesn't produce a plot
		y_lim=max(data_for_plot$Frequency[start:xmax-1])
   		y_lim=y_lim*1.1
   		x_lim=max(which(data_for_plot$Frequency[1:length(data_for_plot$Frequency)-1]>100))
   		x_lim=x_lim*1.1
        print(ggplot(data_for_plot,aes(x=Coverage,y=Frequency,fill="Frequency")) + geom_bar(stat="identity",colour=bar_colour) + scale_fill_manual(name="Data",values=cols)+coord_cartesian(ylim = c(0, y_lim),xlim=c(0,x_lim)))
	    dev.off()
	    fileConn<-file(paste(foldername,"/","progress.txt",sep=""))
	    writeLines("Error: Model did not converge", fileConn)
	    close(fileConn)
	}
    
	
	
}
	
	#################### OLD STUFF:#############################
#	if(class(nlspoisson1) != "try-error"&&class(nlspoisson2) != "try-error"){  #edited by FS
#		rse1=NULL
#		rse2=NULL
#		rse1=try(summary(nlspoisson1)$sigma)
#		rse2=try(summary(nlspoisson2)$sigma)
#		print("Both work")
#		res<-data.frame(x,pred=predict(nlspoisson2))
#		nlsD=coef(nlspoisson2)[1]
#		r=coef(nlspoisson2)[2]
#		nlsKmercov=coef(nlspoisson2)[3]
#		nlsBias=coef(nlspoisson2)[4]
#		nlsLength=coef(nlspoisson2)[5]
#		r = round(r,digits = 4)*100
#		nlsBias=round(nlsBias,digits=4)
#		nlsD=round(nlsD,digits=4)
#		nlsLength=format(nlsLength,scientific=TRUE)
#		print("Creating plot...")
#		p=p[1:nrow(p)-1,]
#		png(paste(foldername,"/","_NLS2_plot.png",sep=""))   #edited by FS 
#	##	png(paste(foldername,"/",title,"_NLS2_plot.png",sep=""))
#		#png("test_plot.png")
#		try({plot(p,type="h",main=paste(title," NLS 2 - K-mer Analysis (k=",k,")\n Heterozygosity = ",r,"% Haploid Genome Size = ",nlsLength,"\n Read Duplication Level = ",nlsBias," % of Repeats = ",nlsD,sep=""),xlab="Coverage",ylab="Frequency",ylim=c(0,max(y[start:length(y)])*1.5),xlim=c(0,xmax))
#		lines(x,res$pred,col="Blue")
#
#		legend("topright",c("Original Data","Chosen NLS"),col=c("black","blue"),lwd=1)
#		dev.off()})
#		res<-data.frame(x,pred=predict(nlspoisson1))
#		nlsLength=coef(nlspoisson1)[5]
#		r=coef(nlspoisson1)[2]
#		nlsKmercov=coef(nlspoisson1)[3]
#		nlsBias=coef(nlspoisson1)[4]
#		nlsD=coef(nlspoisson1)[1]
#
#	}
#	if(class(nlspoisson1) != "try-error"&&class(nlspoisson2) == "try-error"){  #edited by FS
#		print("Half the max didn't work, use max")
#		res<-data.frame(x,pred=predict(nlspoisson1))
#		nlsLength=coef(nlspoisson1)[5]
#		r=coef(nlspoisson1)[2]
#		nlsKmercov=coef(nlspoisson1)[3]
#		nlsBias=coef(nlspoisson1)[4]
#		nlsD=coef(nlspoisson1)[1]
#	}
#	if(class(nlspoisson1) == "try-error"&&class(nlspoisson2) != "try-error"){  #edited by FS
#		print("The max didn't work, use half the max")
#		res<-data.frame(x,pred=predict(nlspoisson2))
#		nlsLength=coef(nlspoisson2)[5]
#		r=coef(nlspoisson2)[2]
#		nlsKmercov=coef(nlspoisson2)[3]
#		nlsBias=coef(nlspoisson2)[4]
#		nlsD=coef(nlspoisson2)[1]
#	}
#	if(class(nlspoisson1) == "try-error"&&class(nlspoisson2) == "try-error"){ #edited by FS
#		 print("Error: Both models failed. Possible reason is we could not detect two curves")
#		 png(paste(foldername,"/","_plot.png",sep=""))   #edited by FS   
#		 try(plot(p,type="h",main="Error, model did not converge"))
#		dev.off()
#		stop("Model did not converge. Please check the plot file.")	
#	}
#	r = round(r,digits = 4)*100
#	nlsBias=round(nlsBias,digits=4)
#	nlsD=round(nlsD,digits=4)
#	nlsLength=format(nlsLength,scientific=TRUE)
#	print("Creating plot...")
#	p=p[1:nrow(p)-1,]
#	png(paste(foldername,"/","plot.png",sep=""))   #edited by FS	
##png(paste(foldername,"/",title,"_plot.png",sep=""))
#	#png("test_plot.png")
#	try({plot(p,type="h",main=paste(title," - K-mer Analysis (k=",k,")\n Heterozygosity = ",r,"% Haploid Genome Size = ",nlsLength,"\n Read Duplication Level = ",nlsBias," % of Repeats = ",nlsD,sep=""),xlab="Coverage",ylab="Frequency",ylim=c(0,max(y[start:length(y)])*1.5),xlim=c(0,xmax))
#	lines(x,res$pred,col="Blue")
#
#	legend("topright",c("Original Data","Chosen NLS"),col=c("black","blue"),lwd=1)
#	dev.off()})
#	print("Writing results to file...")	
#	##try({capture.output(summary(nlspoisson1),file=paste(foldername,"/","_results.txt",sep=""))})   #edited by FS
#        ##try({capture.output(summary(nlspoisson2),file=paste(foldername,"/","_results.txt",sep=""),append=TRUE)})   #edited by FS
#	##try({capture.output(anova(nlspoisson1,nlspoisson2),file=paste(foldername,"/","_results.txt",sep=""),append=TRUE)})  #edited by FS
#	
#
#	fileConn<-file(paste(foldername,"/","summary.txt",sep=""))
#	size=min_max(summary(nlspoisson2)$coefficients[5,])
#	het=min_max(summary(nlspoisson2)$coefficients[2,])
#	repeats=min_max(summary(nlspoisson2)$coefficients[1,])
#	dups=min_max(summary(nlspoisson2)$coefficients[4,])
#	##writeLines(paste("k\t",k,"\nHeterozygosity\t",r,"\nHaploid Genome Size\t",nlsLength,"\nRead Duplication Level\t",nlsBias," %\n of Repeats\t",nlsD,"\n\n",sep=""), fileConn)
#
#	writeLines(paste("property\tmin\tmax\n","k\t",k,"\nHeterozygosity(%)\t",het[1],"\t",het[2],"\nHaploid Genome Size(bp)\t",round(size[1]),"\t",round(size[2]),"\nRead Duplication Level(times)\t",dups[1],"\t",dups[2]," \nRepeats\t",repeats[1],"\t",repeats[2],"\n",sep=""), fileConn)
#
#	close(fileConn)
#	##try({capture.output(summary(nlspoisson1),file=paste(foldername,"/","_results.txt",sep=""),append=TRUE)})   #edited by FS
#       ##try({capture.output(summary(nlspoisson2),file=paste(foldername,"/","_results.txt",sep=""),append=TRUE)})   #edited by FS
#
#
#	try({capture.output(summary(nlspoisson1),file=paste(foldername,"/",title,"_results.txt",sep=""))})
#	try({capture.output(summary(nlspoisson2),file=paste(foldername,"/",title,"_results.txt",sep=""),append=TRUE)})
#
#}
#
#
