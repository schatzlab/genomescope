prefix=${reads/.fq/}
file=$prefix'.bwa-mem.mapped.sort.bam'
/home/schatz/fsedlaze/programs/picard-tools MarkDuplicates I=$file M=$file'.picard.dups' O=$file'.tmp'
