## Changes to original GenomeScope
GenomeScope.R has been updated to include ploidy as a parameter (p = 2 to 6); it can now be run in the command line as:

    $ Rscript genomescope.R histogram_file k-mer_length ploidy read_length output_dir [kmer_max] [verbose]

In the analysis/genomesim folder, we have added project.py, mobius1, mobius2, mobius3, mobius4, mobius5, and mobius6. Project.py outputs the updated model which is compared against the true kmer profile obtained by jellyfish. The mobius files are used by project.py to output the updated model. See the New Model section of the final writeup for an explanation of the mobius functions.

Also in the analysis/genomesim folder, we have updated simgenome.pl, eval_error.pl, simhisto.sh, run_sweep.sh, and run_eval.sh to include ploidy as a parameter. To test our model against simulated data, run the following:

    $ ./run_sweep.sh 31
    $ ./run_eval.sh sweep.k31

The run_sweep.sh file creates the simulated kmer profiles as well as the profiles predict by the updated model of a given ploidy, duplication rate, and heterozgyosity. The ploidy ranges from 2 to 10, the duplication rate from 0 to 1 by 0.1, and the heterozygosity rate for 0.000, 0.001, 0.010, 0.050, 0.100, 0.200, and 0.300. This amounts to 385 total simulations. The run_eval.sh file will evaluate the 385 models against the true data, and evaluate the accuracies. These accuracies were used to create the results histogram file of the final writeup. For all 385 simulations, the accuracies were greater than 98.6%.

In the analysis/scripts folder, chromosomeMutator.py and parameteranalysis.py have been updated. The user should also download the Triticum aestivum genome into analysis/scripts from ftp://ftp.ensemblgenomes.org/pub/plants/release-34/fasta/triticum_aestivum/dna/Triticum_aestivum.TGACv1.dna.toplevel.fa.gz. Afterwards, the user should run:

    $ head -n100000 Triticum_aestivum.TGACv1.dna.toplevel.fa > T_aestivum100000.fa
    $ grep -v '>' T_aestivum100000.fa > T_aestivum100000_2.fa
    $ grep -v 'N' T_aestivum100000_2.fa > T_aestivum100000.fa
    $ python parameteranalysis.py 31 2 100 100 1 0.01 T_aestivum100000.fa 0.01
    $ python parameteranalysis.py 31 3 100 100 1 0.01 T_aestivum100000.fa 0.01
    $ python parameteranalysis.py 31 4 100 100 1 0.01 T_aestivum100000.fa 0.01
    $ python parameteranalysis.py 31 5 100 100 1 0.01 T_aestivum100000.fa 0.01
    $ jellyfish count -C -m 21 -s 1000000000 -t 10 -o reads2.jf <(zcat T_aestivum1000000.fa_p2_het0.01_br1_rl100_cov100_err0.01_reads.fa.gz)
    $ jellyfish count -C -m 21 -s 1000000000 -t 10 -o reads3.jf <(zcat T_aestivum1000000.fa_p3_het0.01_br1_rl100_cov100_err0.01_reads.fa.gz)
    $ jellyfish count -C -m 21 -s 1000000000 -t 10 -o reads4.jf <(zcat T_aestivum1000000.fa_p4_het0.01_br1_rl100_cov100_err0.01_reads.fa.gz)
    $ jellyfish count -C -m 21 -s 1000000000 -t 10 -o reads5.jf <(zcat T_aestivum1000000.fa_p5_het0.01_br1_rl100_cov100_err0.01_reads.fa.gz)
    $ jellyfish histo -t 10 reads2.jf > reads2.histo
    $ jellyfish histo -t 10 reads3.jf > reads3.histo
    $ jellyfish histo -t 10 reads4.jf > reads4.histo
    $ jellyfish histo -t 10 reads5.jf > reads5.histo
    $ Rscript genomescope.R reads2.histo 21 2 100 output2
    $ Rscript genomescope.R reads3.histo 21 3 100 output3
    $ Rscript genomescope.R reads4.histo 21 4 100 output4
    $ Rscript genomescope.R reads5.histo 21 5 100 output5

This will produce the results in the final report table. The code first creates a subset of the T. aestivum genome for testing, and removes the > symbol and any Ns. Then it runs parameteranalysis.py to produce the simulated reads from the T. aestivum genome. Then jellyfish count gets the kmer reads.jf files and jellyfish histo produces the reads.histo files. Finally, genomescope.R is run on the reads.histo files and the output is stored to different output files.

## Getting Started

Before running GenomeScope, you must first compute the histogram of k-mer frequencies. We recommend the tool [jellyfish](https://academic.oup.com/bioinformatics/article/27/6/764/234905/A-fast-lock-free-approach-for-efficient-parallel)  that is available here: http://www.genome.umd.edu/jellyfish.html. After compiling jellyfish, you can run it like this:

    $ jellyfish count -C -m 21 -s 1000000000 -t 10 *.fastq -o reads.jf
  
Note you should adjust the memory (-s) and threads (-t) parameter according to your server. This example will use 10 threads and 1GB of RAM. The kmer length (-m) may need to be scaled if you have low coverage or a high error rate. We recommend using a kmer length of 21 (m=21) for most genomes, as this length is sufficiently long that most k-mers are not repetitive and is short enough that the analysis will be more robust to sequencing errors. Extremely large (haploid size >>10GB) and/or very repetitive genomes may benefit from larger kmer lengths to increase the number of unique k-mers. Accurate inferences requires a minimum amount of coverage, at least 25x coverage of the haploid genome or greater, otherwise the model fit will be poor or not converge. GenomeScope also requires relatively low error rate sequencing, such as Illumina sequencing, so that most k-mers do not have errors in them. For example, a 2% error corresponds to an error every 50bp on average, which is greater than the typical k-mer size used (k=21). Notably, raw single molecule sequencing reads from Oxford Nanopore or Pacific Biosciences, which currently average 5-15% error, are not supported as an error will occur on average every 6 to 20 bp. You should always use "canonical kmers" (-C) since the sequencing reads will come from both the forward and reverse strand of DNA.

Then export the kmer count histogram

    $ jellyfish histo -t 10 reads.jf > reads.histo

Again the thread count (-t) should be scaled according to your server. After you have the jellyfish histogram file, you can run GenomeScope within the online web tool, or at the command line.

### Running GenomeScope on the Command Line

Command line users can run the modeling with the R script genomescope.R, making sure that Rscript is in your PATH (alternatively, edit the shebang line to point to the Rscript location)

    $ Rscript genomescope.R histogram_file k-mer_length ploidy read_length output_dir [kmer_max] [verbose]

The histogram_file (from jellyfish), k-mer_length, ploidy, read_length, and output_dir are required parameters. The optional parameter kmer_max specifies the cutoff for excluding high frequency kmers from the analysis. We recommend setting kmer_max=1000, but it depends on the specific characteristics of your data. If the optional parameter verbose is set to 1, then GenomeScope will output more details on the model fitting as different parameters are used. The output plots and a text file of the inferred genome characteristics will be output to the specified output_dir directory.
