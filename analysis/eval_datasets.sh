#!/bin/sh

K=21
READLEN=100
MODELBIN=../genomescope.R


#for INDIR in simulation_results
#for INDIR in real_data
for INDIR in simulation_results_100x simulation_results_50x simulation_results_25x simulation_results_15x simulation_results_10x real_data ecoli_mix
do
  echo "Processing $INDIR"

  if [ $(which parallel) ] 
  then
    echo "Running results in parallel"
    #/bin/ls $INDIR/seabass*.hist | parallel -t Rscript $MODELBIN {} $K $READLEN {}
    #/bin/ls $INDIR/D*.hist | parallel -t Rscript $MODELBIN {} $K $READLEN {}
    /bin/ls $INDIR/*.hist | parallel -t Rscript $MODELBIN {} $K $READLEN {}_results 1000

  else
    for hist in `/bin/ls $INDIR/*21.hist`
    do
      echo "Processing $hist"
      Rscript $MODELBIN $hist $K $READLEN ${hist}_results
    done
  fi
done







