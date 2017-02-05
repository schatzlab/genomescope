#!/bin/sh

if [[ $# != 4 ]]
then
  echo "simhist.sh LEN DUP HET KMER\n";
  exit
fi

LEN=$1
DUP=$2
HET=$3
KMER=$4

./simgenome.pl $LEN $DUP $HET $KMER > genome.fa


echo "\n"
echo "Observed Counts:"

jellyfish count -t 2 -m $KMER -s 100000000 -C genome.fa 
jellyfish histo mer_counts.jf
