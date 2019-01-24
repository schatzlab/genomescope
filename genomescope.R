#!/usr/bin/env Rscript

## GenomeScope: Fast Genome Analysis from Unassembled Short Reads
## This is the automated script for computing genome characteristics
## from a kmer histogram file, k-mer size, and readlength

## Load libraries for non-linear least squares and argument parser
library('minpack.lm')
library('argparse')

## Load the genomescope library
library('genomescope')

## Number of rounds before giving up
NUM_ROUNDS=4

## Coverage steps to trim off between rounds
START_SHIFT=5

## Typical cutoff for sequencing error
TYPICAL_ERROR = 15

## Max rounds on NLS
MAX_ITERATIONS=200

## Overrule if two scores are within this percent (0.05 = 5%) but larger difference in het
SCORE_CLOSE = 0.20 #originally 0.20

## Overrule heterozygosity if there is a large difference in het rate
SCORE_HET_FOLD_DIFFERENCE = 10 #originally 10

## Suppress the warnings if the modeling goes crazy, those are in try/catch blocks anyways
options(warn=-1)

## Colors for plots
COLOR_BGCOLOR  = "light grey"
COLOR_HIST     = "#56B4E9"
COLOR_2pPEAK    = "black"
COLOR_pPEAK    = "#F0E442"
COLOR_ERRORS   = "#D55E00"
COLOR_KMERPEAK = "black"
COLOR_RESIDUAL = "purple"
COLOR_COVTHRES = "red"

library(minpack.lm)

## Given mean +/- stderr, report min and max value within 2 SE
###############################################################################

min_max <- function(table) {
  return (c(max(0,table[1] - 2*table[2]), table[1]+ 2*table[2]))
}

min_max1 <- function(table) {
  return (c(max(0,table[1] - 2*table[2]), min(1, table[1]+ 2*table[2])))
}

## Main program starts here
###############################################################################

parser <- ArgumentParser()
parser$add_argument("-v", "--version", action="store_true", default=FALSE, help="print the version and exit")
parser$add_argument("-i", "--input", help = "input histogram file")
parser$add_argument("-o", "--output", help = "output directory name")
parser$add_argument("-k", "--kmer_size", type = "integer", default = 21, help = "kmer size used to calculate kmer spectra [default 21]")
parser$add_argument("-p", "--ploidy", type = "integer", default = 2, help = "ploidy (1, 2, 3, 4, 5, or 6) for model to use [default 2]")
parser$add_argument("-l", "--lambda", "--kcov", "--kmercov", type = "integer", default=-1, help = "optional initial kmercov estimate for model to use")
parser$add_argument("-m", "--max_kmercov", type = "integer", default=-1, help = "optional maximum kmer coverage threshold (kmers with coverage greater than max_kmercov are ignored by the model)")
parser$add_argument("-n", "--name_prefix", default = "OUTPUT", help = "name prefix for output files [default OUTPUT]")
parser$add_argument("--verbose", action="store_true", default=FALSE, help = "optional flag to print messages during execution")

arguments <- parser$parse_args()
version_message <- "GenomeScope 2.0\n"

if (arguments$version) {
  cat(version_message)
  quit()
}

if (is.null(arguments$input) | is.null(arguments$output)) {
  cat("USAGE: genomescope.R -i input_histogram_file -k kmer_length -p ploidy -o output_dir\n")
  cat("OPTIONAL PARAMETERS: -l lambda -m max_kmercov -n 'name_prefix' --verbose\n")
  cat("HELP: genomescope.R --help\n")
} else {

  ## Load the arguments from the user
  histfile   <- arguments$input
  k          <- arguments$kmer_size
  p          <- arguments$ploidy
  foldername <- arguments$output
  estKmercov <- arguments$lambda
  max_kmercov <- arguments$max_kmercov
  VERBOSE <- arguments$verbose

  cat(paste("GenomeScope analyzing ", histfile, " k=", k, " p=", p, " outdir=", foldername, "\n", sep=""))

  dir.create(foldername, showWarnings=FALSE)

  kmer_prof <- read.csv(file=histfile,sep="", header=FALSE)

  minkmerx = 1;
  if (kmer_prof[1,1] == 0) {
    if (VERBOSE) {cat("Histogram starts with zero, reseting minkmerx\n")}
    minkmerx = 2;
  }

  kmer_prof <- kmer_prof[c(minkmerx:(length(kmer_prof[,2])-1)),] #get rid of the last position
  kmer_prof_orig <- kmer_prof

  ## Initialize the status
  progressFilename <- paste(foldername,"/", arguments$name_prefix, "_progress.txt",sep="")
  cat("starting", file=progressFilename, sep="\n")

  ## try to find the local minimum between errors and the first (heterozygous) peak
  start <- tail(which(kmer_prof[1:TYPICAL_ERROR,2]==min(kmer_prof[1:TYPICAL_ERROR,2])),n=1)

  maxCovIndex = -1

  ## Figure out which kmers to exclude, if any
  if(max_kmercov == -1) {
    maxCovIndex <- length(kmer_prof[,1])
    max_kmercov <- kmer_prof[maxCovIndex,1]
  }
  else {
    ## Figure out the index we should use for this coverage length
    x <- kmer_prof[,1]
    maxCovIndex <- length(x[x<=max_kmercov])
  }

  if (VERBOSE) {cat(paste("using max_kmercov:", max_kmercov, " with index:", maxCovIndex, "trying 2p peak model... \n"))}

  # terminate after NUM_ROUND iterations, store best result so far in container
  round <- 0
  best_container <- list(NULL,0)

  while(round < NUM_ROUNDS) {
    cat(paste("round", round, "trimming to", start, "trying 2p peak model... "), file=progressFilename, sep="", append=TRUE)
    if (VERBOSE) {cat(paste("round", round, "trimming to", start, "trying 2p peak model... \n"))}

    ## Reset the input trimming off low frequency error kmers
    kmer_prof=kmer_prof_orig[1:maxCovIndex,]
    x <- kmer_prof[start:maxCovIndex,1]
    y <- kmer_prof[start:maxCovIndex,2]

    model_peaks <- estimate_Genome_peakp(kmer_prof, x, y, k, p, estKmercov, round, foldername, arguments)

    if (!is.null(model_peaks[[1]])) {
      cat(paste("converged. score: ", model_peaks[[2]]$all[[1]]), file=progressFilename, sep="\n", append=TRUE)

      if (VERBOSE) {
        mdir = paste(foldername, "/round", round, sep="")
        dir.create(mdir, showWarnings=FALSE)
        report_results(kmer_prof,kmer_prof_orig, k, p, model_peaks, mdir, arguments)
      }
    }
    else {
     cat(paste("unconverged"), file=progressFilename, sep="\n", append=TRUE)
    }

    #check if this result is better than previous
    if (!is.null(model_peaks[[1]])) {
      if (is.null(best_container[[1]])) {
        if (VERBOSE) {cat("no previous best, updating best\n")}
        best_container = model_peaks
      }
      else {
        pdiff = abs(model_peaks[[2]]$all[[1]] - best_container[[2]]$all[[1]]) / max(model_peaks[[2]]$all[[1]], best_container[[2]]$all[[1]])

        if (pdiff < SCORE_CLOSE) {
          hetm = model_peaks[[1]]$ahet
          hetb = best_container[[1]]$ahet

          if (hetb * SCORE_HET_FOLD_DIFFERENCE < hetm) {
            if (VERBOSE) {cat("model has significantly higher heterozygosity but similar score, overruling\n")}
            best_container = model_peaks
          }
          else if (hetm * SCORE_HET_FOLD_DIFFERENCE < hetb) {
            if (VERBOSE) {cat("previous best has significantly higher heterozygosity and similar score, keeping\n")}
          }
          else if (model_peaks[[2]]$all[[1]] < best_container[[2]]$all[[1]]) {
            if (VERBOSE) {cat("score is marginally better but het rate is not extremely different, updating\n")}
            best_container = model_peaks
          }
        }
        else if (model_peaks[[2]]$all[[1]] < best_container[[2]]$all[[1]]) {
          if (VERBOSE) {cat("score is significantly better, updating\n")}
          best_container = model_peaks
        }
      }
    }

    ## Ignore a larger number of kmers as errors
    start <- start + START_SHIFT
    round <- round + 1
  }
  ## Report the results, note using the original full profile
  report_results(kmer_prof,kmer_prof_orig, k, p, best_container, foldername, arguments)
  if (!is.null(best_container[[1]])) {
    print('model score')
    print(best_container[[2]]$all[[1]])
    print(best_container[[1]]$m$deviance())
  }
