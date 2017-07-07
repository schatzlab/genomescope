#!/bin/sh

if [[ $# != 6 ]]
then
  echo "simhist.sh LEN DUP HET1 HET2 HETZ KMER\n";
  exit
fi

LEN=$1
DUP=$2
HET1=$3
HET2=$4
HETZ=$5
KMER=$6

dir=pp-$LEN-$DUP-$HET1-$HET2-$HETZ-$KMER
mkdir -p $dir
cd $dir

../simpolyploid.pl $LEN $DUP $HET1 $HET2 $HETZ $KMER > genome.fa

echo "\n"
echo "Observed Counts:"

jellyfish count -t 2 -m $KMER -s 100000000 -C genome.fa 
jellyfish histo mer_counts.jf


dwgsim -r 0 -F 0 -y 0 -C 10 -e 0 -E 0 genome.fa reads
jellyfish count -m 21 -s 1000000 -t 4 -C -o reads.jf reads.bwa.read*
jellyfish histo reads.jf > reads.histo
~/build/genomescope/genomescope.R reads.histo 21 100 reads.gs
