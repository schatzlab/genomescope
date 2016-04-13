/sonas-hs/schatz/hpc/home/fsedlaze/programs/samtools view $reads  |cut -f12 | sed 's/NM:i://g' | awk '{print $1/=100}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' > $reads'.seqerr'
