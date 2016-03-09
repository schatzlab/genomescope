#!/bin/sh

K=21
READLEN=100
INDIR=simulation_results
MODELBIN=../histfitdup_bothplots.R


if [ $(which parallel) ] 
then
  echo "Running results in parallel"
  /bin/ls $INDIR/*21.hist | parallel -t Rscript $MODELBIN {} $K $READLEN {}

else
  for hist in `/bin/ls $INDIR/*21.hist`
  do
    echo "Processing $hist"
    Rscript $MODELBIN $hist $K $READLEN $hist
  done
fi







