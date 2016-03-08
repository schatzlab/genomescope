#!/bin/sh

K=21
READLEN=100
INDIR=simulation_results
OUTDIR=simulation_analysis
MODELBIN=../histfitdup_bothplots.R

mkdir -p $OUTDIR

set -xv

for hist in Drosophila_melanogaster.BDGP6.dna_sm.main_chr.fa_het0.01_br1_rl100_cov100_err0.01_reads.fa21.hist
do
  echo "Processing $i"
  Rscript $MODELBIN $INDIR/$hist $K $READLEN $OUTDIR/$hist
done


