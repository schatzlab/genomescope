# GenomeScope: Fast genome analysis from unassembled short reads

Current developments in de novo assembly technologies have been focused on relatively simple genomes. Even the human genome, with a heterozygosity rate of only ~0.1% and 2n diploid structure, is significantly simpler than many other species, especially plants. However, genomics is rapidly advancing towards sequencing more complex species such as pineapple, sugarcane, or wheat that have much higher rates of heterozygosity (>1% for pineapple), much higher ploidy (8n for sugarcane), and much larger genomes (16Gbp for wheat).

One of the first goals when sequencing a new species is determining the overall characteristics of the genome structure, including the genome size, abundance of repetitive elements, and the rate of heterozygosity. These features are needed to study trends in genome evolution, and can inform the parameters that should be used for the individual assembly steps. They can also serve as an independent quality control during any analysis, such as quantifying the quality of an assembly, or measuring the expected number of heterozygous bases in the genome before mapping any variants.

We have developed an analytical model and open-source software package GenomeScope that can infer the global properties of a genome from unassembled sequenced data. GenomeScope uses the k-mer count distribution, e.g. from Jellyfish, and within seconds produces a report and several informative plots describing the genome properties. We validate the approach on simulated heterozygous genomes, as well as synthetic crosses of related strains of microbial and eukaryotic genomes with known reference genomes. GenomeScope was also applied to study the characteristics of several novel species, including pineapple, pear, the regenerative flatworm Macrostomum lignano, and the Asian sea bass.

## Getting Started

Before running GenomeScope, you must first compute the histogram of k-mer frequencies. We recommend the tool [jellyfish](https://academic.oup.com/bioinformatics/article/27/6/764/234905/A-fast-lock-free-approach-for-efficient-parallel)  that is available here: http://www.genome.umd.edu/jellyfish.html. After compiling jellyfish, you can run it like this:

    $ jellyfish count -C -m 21 -s 1000000000 -t 10 *.fastq -o reads.jf
  
Note you should adjust the memory (-s) and threads (-t) parameter according to your server. This example will use 10 threads and 1GB of RAM. The kmer length (-m) may need to be scaled if you have low coverage or a high error rate. We recommend using a kmer length of 21 (m=21) for most genomes, as this length is sufficiently long that most k-mers are not repetitive and is short enough that the analysis will be more robust to sequencing errors. Extremely large (haploid size >>10GB) and/or very repetitive genomes may benefit from larger kmer lengths to increase the number of unique k-mers. Accurate inferences requires a minimum amount of coverage, at least 25x coverage of the haploid genome or greater, otherwise the model fit will be poor or not converge. GenomeScope also requires relatively low error rate sequencing, such as Illumina sequencing, so that most k-mers do not have errors in them. For example, a 2% error corresponds to an error every 50bp on average, which is greater than the typical k-mer size used (k=21). Notably, raw single molecule sequencing reads from Oxford Nanopore or Pacific Biosciences, which currently average 5-15% error, are not supported as an error will occur on average every 6 to 20 bp. You should always use "canonical kmers" (-C) since the sequencing reads will come from both the forward and reverse strand of DNA.

Then export the kmer count histogram

    $ jellyfish histo -t 10 reads.jf > reads.histo

Again the thread count (-t) should be scaled according to your server. After you have the jellyfish histogram file, you can run GenomeScope within the online web tool, or at the command line.

### Running GenomeScope Online

Users may prefer to use the online version, which offers all of the same functionality within an easy to use web interface:
http://qb.cshl.edu/genomescope/

### Running GenomeScope on the Command Line

Command line users can run the modeling with the R script genomescope.R, making sure that Rscript is in your PATH (alternatively, edit the shebang line to point to the Rscript location)

    $ Rscript genomescope.R histogram_file k-mer_length read_length output_dir [kmer_max] [verbose]

The histogram_file (from jellyfish), k-mer_length, read_length, and output_dir are required parameters. The optional parameter kmer_max specifies the cutoff for excluding high frequency kmers from the analysis. We recommend setting kmer_max=1000, but it depends on the specific characteristics of your data. If the optional parameter verbose is set to 1, then GenomeScope will output more details on the model fitting as different parameters are used. The output plots and a text file of the inferred genome characteristics will be output to the specified output_dir directory.

For example, you can download the histogram from the Arabidopsis F1 described in the manuscript here:
https://raw.githubusercontent.com/schatzlab/genomescope/master/analysis/real_data/ara_F1_21.hist

Then run GenomeScope like this:
    
    $ /PATH/TO/Rscript /PATH/TO/genomescope.R ara_F1_21.hist 21 150 output
    
This should complete in less than 1 minute, and report:

    Model converged het:0.0104 kcov:22.2 err:0.0035 model fit:0.446 len:151975724

The plots and the full results will be in output directory, showing the estimated genome size to be 151.9Mbp and a 1.04% heterozygosity rate (the exact values may slightly differ due to the randomization within the modeling)


## Frequently Asked Questions (FAQ)

**> Q: Why didnt the model converge, or why does it give very different results than expected?**

A: The most common problem is you have too low of coverage for the model to confidently identify the peaks corresponding to the homozygous kmers. This can be recognized by a lack of any peaks in the kmer plots or that the model fit doesnt match the observed kmer profile very well. To correct this problem, first make sure that you have used the cannonically kmer counting mode (-C if you are using jellyfish). If this still fails, you can try slightly decreasing the kmer size used to 17 or 19. If all of these attempts still fail you will unfortunately need to generate additional sequencing data.

**> Q:  As mentioned in the Supplementary Notes and Figures 1.3.2 Genome Size Estimation, the haploid genome size is estimated by: "This estimate is revised by summing the total number of k-mers, except presumptive sequencing errors identified as in section 1.3.1, and dividing by the 2*λ, the estimated coverage for homozygous k-mers". If I understand it correctly, λ is the mean of a distribution, the estimated coverage for homozygous k-mers; and in the genomescope profile, kcov is the estimated coverage for heterozygous kmers. 
Could you explain how genomescope estimates haploid genome size, specifically why dividing by 2 times of the estimated coverage for homozygous k-mers?**

A: Thank you for writing. I think there is a bit of confusion over the variable names and how they relate to each other. The first thing to note is λ and kcov refer to the same value, just that we use λ in the written document and kcov in the code. The modeling tries to identify 4 peaks centered at λ, 2*λ, 3*λ, and 4*λ. These 4 peaks correspond to the mean coverage levels of the unique heterozygous, unique homozygous, repetitive heterozygous and repetitive homozygous sequences, respectively. So when it estimates the haploid genome size, it divides by the 2*λ, which is the average homozygous coverage, not "2 times the estimated coverage for homozygous k-mers" as you write.

The other potentially confusing aspect is what is meant by haploid genome size versus diploid genome size. We consider the haploid genome size to mean the span of one complete set of haploid chromosomes and the diploid genome size to be the span of both haploid copies (total DNA content in one diploid cell). In particular, in a human cell, the haploid genome size is about 3Gbp and the diploid genome size is about 6Gbp. If you sequence a total of 300Gbp for a human genome, that would be about 150Gbp (50x coverage) of the maternal haplotype and about 150Gbp (50x coverage) of the paternal haplotype. But since the heterozygosity rate in humans is so low, the main peak in the distribution would be centered around 100x. However, GenomeScope will still try to fit the 4 peaks, so should set the heterozygous kmer coverage λ equal to 50x, and thus the homozygous coverage to 2*λ = 100x. From this GenomeScope will compute the haploid genome size as the total amount of sequence data (300GB) divided by the homozygous coverage (100x) to report 3Gbp as expected. Kmers with higher coverage are naturally scaled as well: kmers that occur 200 or 300 times in the kmer profile (and thus are 2 or 3 copy repeats in the diploid genome, 4 or 6 times in the haploid sequences) are still scaled by 100x to contribute 2 or 3 copies to the estimate. Finally, note that if the two haplotypes have significantly different lengths, then the reported haploid genome size will be the average of the two.

**> Q: Can GenomeScope be used to estimate ploidy, or used with genomes with higher ploidy?**

A: No, GenomeScope is only appropriate for diploid genomes. In principle the model could be extend to higher levels of polyploidy by extending the model to consider additional peaks in the k-mer profile, but this is not currently supported. GenomeScope also does not support genomes that have uneven copy number of their chromosomes, such as aneuploid cancer genomes or even unequal numbers of sex chromosomes. In these scenarios the reported heterozygosity rate will represent the fraction of bases that are haploid (copy number 1) versus diploid (copy number 2) as well as any heterozygous positions in the other chromosomes.


## Resources

VCF files of the variants identified in the larger genomes are available here:
http://labshare.cshl.edu/shares/schatzlab/www-data/genomescope/vcf/
