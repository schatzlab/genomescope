#!/bin/sh

grep GenomeHaploidLen */summary.txt | column -t | less
grep GenomeUniqueLen */summary.txt | column -t | less
grep GenomeRepeatLen */summary.txt | column -t | less
grep Heterozygosity */summary.txt | rev | sort -t '_' -k7,7 -k8,8 | rev | column -t | less
grep Error */summary.txt | rev | sort -t '_' -k3,3 -k8,8 | rev | column -t | less
grep ReadDuplicationLevel */summary.txt | rev | sort -t '_' -k6,6 -k8,8 | rev | column -t | less
