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

./simgenome.pl $LEN $DUP $HET1 $HET2 $HETZ $KMER > genome.fa


echo "\n"
echo "Observed Counts:"

jellyfish count -t 2 -m $KMER -s 100000000 -C genome.fa 
jellyfish histo mer_counts.jf
