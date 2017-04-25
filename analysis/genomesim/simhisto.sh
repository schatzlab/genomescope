#!/bin/sh

if [[ $# != 5 ]]
then
  echo "simhist.sh LEN DUP HET KMER PLOIDY\n";
  exit
fi

LEN=$1
DUP=$2
HET=$3
KMER=$4
PLOIDY=$5

./simgenome.pl $LEN $DUP $HET $KMER $PLOIDY > genome.fa


echo "\n"
echo "Observed Counts:"

jellyfish count -t 2 -m $KMER -s 100000000 -C genome.fa 
jellyfish histo mer_counts.jf
