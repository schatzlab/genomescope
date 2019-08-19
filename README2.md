# GenomeScope 2.0: Reference-free profiling of polyploid genomes
## T. Rhyker Ranallo-Benavidez, Kamil S. Jaron, and Michael C. Schatz

Many genomics analyses require first establishing a reference genome. However, *de novo* genome assembly is a complicated and computationally intensive process with many assumptions hidden to the user. A popular assessment prior to genome assembly is genome profiling, where the k-mer frequencies within the sequencing reads are analyzed to efficiently estimate major genome characteristics such as genome size, heterozygosity, and repetitiveness. However, current genome profiling tools are suited only for diploid genomes and use heuristic approaches.

We have developed GenomeScope 2.0, which applies classical insights from combinatorial theory to establish a detailed mathematical model of how k-mer frequencies will be distributed in heterozygous and polyploid genomes. GenomeScope 2.0 employs a polyploid-aware mixture model that, within seconds, accurately infers genome properties from unassembled sequencing data. GenomeScope 2.0 uses the k-mer count distribution, e.g. from KMC or Jellyfish, and produces a report and several informative plots describing the genome properties. We validate the approach on simulated polyploid data created using a generative model with parameters for genome size, heterozygosity, repetitiveness, ploidy, and sequencing coverage, and find GenomeScope 2.0 retains accuracy across a broad range of realistic and extreme parameter values. We also validate GenomeScope 2.0 by analyzing genuine sequence data from 11 diverse polyploid genomes with known genome characteristics. 

## Getting Started

Before running GenomeScope 2.0, you must first compute the histogram of k-mer frequencies. We recommend [KMC](https://academic.oup.com/bioinformatics/article/33/17/2759/3796399) that is available at http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=download, or [jellyfish](https://academic.oup.com/bioinformatics/article/27/6/764/234905/A-fast-lock-free-approach-for-efficient-parallel) that is available at http://www.genome.umd.edu/jellyfish.html.

After compiling KMC, you can count the k-mers with these commands:

    $ mkdir tmp
    $ ls *.fastq > FILES
    $ kmc -k21 -t10 -m64 -ci1 -cs10000 @FILES reads tmp/

Note you should adjust the memory (-m) and threads (-t) parameters according to your server. This example will use 10 threads and 64GB of RAM. The k-mer length (-k) may need to be scaled if you have low coverage or high error rate. The lower (-ci) and upper (-cs) bounds exclude k-mers with counts outside these boundaries. FILES is a file with a list of input files.

After compiling jellyfish, you can count the k-mers with this command:

    $ jellyfish count -C -m 21 -s 1000000000 -t 10 *.fastq -o reads.jf
  
Note you should adjust the memory (-s) and threads (-t) parameters according to your server. This example will use 10 threads and 1GB of RAM. The k-mer length (-m) may need to be scaled if you have low coverage or a high error rate. You should always use "canonical k-mers" (-C) since the sequencing reads will come from both the forward and reverse strand of DNA.

We recommend using a k-mer length of 21 for most genomes, as this length is sufficiently long that most k-mers are not repetitive and is short enough that the analysis will be more robust to sequencing errors. Extremely large (haploid size >>10GB) and/or very repetitive genomes may benefit from larger k-mer lengths to increase the number of unique k-mers. Accurate inferences requires a minimum amount of coverage, at least 15x coverage of the haploid genome or greater, otherwise the model fit will be poor or not converge. GenomeScope also requires relatively low error rate sequencing, such as Illumina sequencing, so that most k-mers do not have errors in them. For example, a 2% error corresponds to an error every 50bp on average, which is greater than the typical k-mer size used (k=21). Notably, raw single molecule sequencing reads from Oxford Nanopore or Pacific Biosciences, which currently average 5-15% error, are not supported as an error will occur on average every 6 to 20 bp.

After counting the k-mers, you then export the k-mer count histogram. With KMC:

    $ kmc_tools transform reads histogram reads.histo -cx10000

The upper bound (-cx) gives the cutoff for the histogram.
    
With jellyfish:

    $ jellyfish histo -t 10 reads.jf > reads.histo

Again the thread count (-t) should be scaled according to your server.

After you have the histogram file, you can run GenomeScope within the online web tool, or at the command line.

### Running GenomeScope Online

Users may prefer to use the online version, which offers all of the same functionality within an easy to use web interface:
http://genomescope.org/

### Running GenomeScope on the Command Line

Command line users can run the modeling with the R script genomescope.R, making sure that Rscript is in your PATH (alternatively, edit the shebang line to point to the Rscript location).

    $ genomescope.R -i histogram_file -o output_dir -k k-mer_length

The input histogram_file (from KMC or jellyfish), output_dir, and k-mer_length are required parameters. The optional parameter *-p ploidy* sets the ploidy of the model for GenomeScope to use. The optional parameter *-l lambda* sets the initial guess for the average k-mer coverage of the sequencing. The optional parameter *-n 'name_prefix'* sets the prefix for the output files. The optional parameter *-m max_kmercov* specifies the cutoff for excluding high frequence k-mers from the analysis. The output plots and a text file of the inferred genome characteristics will be output to the specified output_dir directory.

For example, you can download the histogram from the Arabidopsis F1 described in the manuscript here:
https://raw.githubusercontent.com/schatzlab/genomescope/master/analysis/real_data/ara_F1_21.hist

Then run GenomeScope like this:
    
    $ /PATH/TO/genomescope.R -i ara_F1_21.hist -o output -k 21
    
This should complete in less than 1 minute, and report:

    Model converged het:0.0104 kcov:22.2 err:0.0035 model fit:0.446 len:151975724

The plots and the full results will be in output directory, showing the estimated genome size to be 151.9Mbp and a 1.04% heterozygosity rate (the exact values may slightly differ due to the randomization within the modeling).

## Tutorial
A tutorial by Andrew Severin on running GenomeScope 1.0 is available here: http://gif.biotech.iastate.edu/genomescope

## Frequently Asked Questions (FAQ)

**> Q: Why didn't the model converge, or why does it give very different results than expected?**

A: The most common problem is that you have too low of coverage for the model to confidently identify the peaks corresponding to the homozygous k-mers. This can be recognized by a lack of any peaks in the k-mer plots or that the model fit doesn't match the observed k-mer profile very well. To correct this problem, first make sure that you have used the canonical k-mer counting mode (-C if you are using jellyfish). If this still fails, you can try slightly decreasing the k-mer size used to 19 or 17. If all of these attempts still fail you will unfortunately need to generate additional sequencing data.

**> Q:  As mentioned in the Supplementary Notes in Section 1.3.2 Genome Size Estimation in the GenomeScope 1.0 paper, the haploid genome size estimate "is revised by summing the total number of k-mers, except presumptive sequencing errors identified as in section 1.3.1, and dividing by 2λ, the estimated coverage for homozygous k-mers." If I understand it correctly, λ is the mean of a distribution, the estimated coverage for homozygous k-mers; and in the genomescope profile, kcov is the estimated coverage for heterozygous k-mers. Could you explain how genomescope estimates haploid genome size, specifically why dividing by 2 times of the estimated coverage for homozygous k-mers?**

A: Thank you for writing. I think there is a bit of confusion over the variable names and how they relate to each other. The first thing to note is λ and kcov refer to the same value, just that we use λ in the written document and kcov in the code. The modeling tries to identify 4 peaks centered at λ, 2λ, 3λ, and 4λ. These 4 peaks correspond to the mean coverage levels of the unique heterozygous, unique homozygous, repetitive heterozygous and repetitive homozygous sequences, respectively. So when it estimates the haploid genome size, it divides by 2λ, which is the average homozygous coverage, not "2 times the estimated coverage for homozygous k-mers" as you write.

The other potentially confusing aspect is what is meant by haploid genome size versus diploid genome size. We consider the haploid genome size to mean the span of one complete set of haploid chromosomes and the diploid genome size to be the span of both haploid copies (total DNA content in one diploid cell). In particular, in a human cell, the haploid genome size is about 3Gbp and the diploid genome size is about 6Gbp. If you sequence a total of 300Gbp for a human genome, that would be about 150Gbp (50x coverage) of the maternal haplotype and about 150Gbp (50x coverage) of the paternal haplotype. But since the heterozygosity rate in humans is so low, the main peak in the distribution would be centered around 100x. However, GenomeScope will still try to fit the 4 peaks, and will set the heterozygous k-mer coverage λ equal to 50x, and thus the homozygous coverage to 2λ = 100x. From this GenomeScope will compute the haploid genome size as the total amount of sequence data (300GB) divided by the homozygous coverage (100x) to report 3Gbp as expected. K-mers with higher coverage are naturally scaled as well: k-mers that occur 200 or 300 times in the k-mer profile (and thus are 2 or 3 copy repeats in the diploid genome, 4 or 6 times in the haploid sequences) are still scaled by 100x to contribute 2 or 3 copies to the estimate. Finally, note that if the two haplotypes have significantly different lengths, then the reported haploid genome size will be the average of the two.

**> Q: Can GenomeScope be used to estimate ploidy, or used with genomes with higher ploidy?**

A: Yes, GenomeScope 2.0 is appropriate for genomes up to a hexaploid (6) ploidy level. To estimate ploidy, you can use [Smudgeplot](https://github.com/KamilSJaron/smudgeplot), which is another tool we developed. GenomeScope currently does not support genomes that have uneven copy number of their chromosomes, however, such as aneuploid cancer genomes or even unequal numbers of sex chromosomes. In these scenarios the reported heterozygosity rate will represent the fraction of bases that are haploid (copy number 1) versus diploid (copy number 2) as well as any heterozygous positions in the other chromosomes.

**> Q: I tried out GenomeScope in order to obtain a genome size estimate for an organism that I am working on. Using the widely-used and cited k-mer method and calculation outlined in the XXX paper, I got an estimate of 869 - 919 Mbp and this is somewhat consistent with a c-value of 0.83pg from a different published study. However, with GenomeScope, I got an estimate of 650 Mbp instead. Would you or your team have any insights into why I am observing this discrepancy of more than 200Mbp in my estimations?**

While GenomeScope can automate most of the analysis, it does have an optional parameter (--max_kmercov) to filter out high frequency k-mers. This parameter is sometimes needed because we see in real samples that k-mers that occur 10,000 times or 100,000 times or greater may be artifacts such as phiX sequencing, or organelle sequencing (or other contamination). To account for these, the user will have to specify a cutoff value for GenomeScope to exclude any k-mers that occur above this threshold. On the other hand, for highly repetitive or huge genomes you may need to expand your k-mer histogram to include higher counter k-mers (jellyfish histograms are truncated at 10,000 by default). If you want to include these ultrahigh frequency k-mers, then you will have to regenerate the histogram from jellyfish or KMC with a higher upper bound (perhaps 100,000, 1 million, or even higher). This will likely increase the estimated genome size, but could also make the analysis more sensitive to any artifacts in the data. Unfortunately every project is a little bit different on how to best remove those artifacts. Please see the GenomeScope 1.0 paper for details (especially the section in the supplement on characterizing the high frequency k-mers in Arabidopsis).


## Resources

VCF files of the variants identified in the larger genomes are available here:
http://labshare.cshl.edu/shares/schatzlab/www-data/genomescope/vcf/

## References:
Vurture, GW, Sedlazeck, FJ, Nattestad, M, Underwood, CJ, Fang, H, Gurtowski, J, Schatz, MC (2017) *Bioinformatics* doi: https://doi.org/10.1093/bioinformatics/btx153
