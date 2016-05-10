#!/bin/sh

grep "Genome Haploid Length" */summary.txt | column -t | less
grep "Genome Unique Length" */summary.txt | column -t | less
grep "Genome Repeat Length" */summary.txt | column -t | less
grep Heterozygosity */summary.txt | rev | sort -t '_' -k7,7 -k8,8 | rev | column -t | less
grep Error */summary.txt | rev | sort -t '_' -k3,3 -k8,8 | rev | column -t | less
grep Duplication */summary.txt | rev | sort -t '_' -k6,6 -k8,8 | rev | column -t | less
grep 'Percent Kmers Modeled (All Kmers)' */summary.txt | column -t | less
