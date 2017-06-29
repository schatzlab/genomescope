#!/usr/bin/perl -w
use strict;
use List::Util qw(any);

my $USAGE = "simgenome.pl GENOME_LEN PLOIDY DUP_RATE KMER_LEN HET_RATE > genome.fa\n";

die $USAGE if (scalar @ARGV < 5);

my $GENOME_LEN = shift @ARGV or die $USAGE;
my ($PLOIDY, $DUP_RATE, $KMER_LEN, @HET_RATE)  = @ARGV;
#my $DUP_RATE   = shift @ARGV;
#my @HET_RATE   = shift @ARGV;
#my $KMER_LEN   = shift @ARGV;
#my $PLOIDY     = shift @ARGV;

die "Dup rate cant be > 1\n" if ($DUP_RATE > 1.0);
die "Het rate cant be > 1\n" if ( any{$_ > 1.0} @HET_RATE);

print STDERR "Simulating G=$GENOME_LEN, d=$DUP_RATE, r=$HET_RATE, k=$KMER_LEN, p=$PLOIDY\n";


## Initialize DNA coding table: 0<->A, 1<->C, 2<->G, 3<->T
###############################################################################
my @DNA = ("A", "C", "G", "T");
my %DNAIDX;
for (my $i = 0; $i < scalar @DNA; $i++) { $DNAIDX{$DNA[$i]} = $i; }


## Simutate a random genome of length GENOME_LEN
###############################################################################

print STDERR " ... simulating progenitor\n";

my @genomes = ("")x$PLOIDY;

for (my $i = 0; $i < $GENOME_LEN; $i++) #for every base
{
  my $z = int(rand(4));
  my $b = $DNA[$z];
  $genomes[0] .= $b;
}


## Now duplicate the first $DUP_RATE percent of the genome
###############################################################################

my $dup = substr($genomes[0], 0, $DUP_RATE * $GENOME_LEN);
$genomes[0] .= $dup;
my $newlen = length($genomes[0]);

print STDERR " ... duplicated $DUP_RATE, newlen=$newlen\n";


## Add random heterozygosity to all haplotypes
###############################################################################

print STDERR " ... add heterozygosity r=$HET_RATE\n";

my @nummut = (0)x$PLOIDY; #number of mutations for each haplotype, used to calculate per-haplotype mutation rate from progenitor
for (my $j = 1; $j < $PLOIDY; $j++) #for every non-progenitor haplotype
{
  $genomes[$j] = $genomes[0]; #set the haplotype equal to the 'progenitor' haplotype
}

for (my $i = 0; $i < $newlen; $i++) #for every base
{
  for (my $j = 1; $j < $PLOIDY; $j++) # for every non-progenitor haplotype
  {
    my $r = rand();
    if ($r <= $HET_RATE) #if the base on this haplotype is a mutation
    {
      $nummut[$j]++; #increase number of mutations for this haplotype
      my $s = 1+int(rand(2)); #random integer from 1, 2, 3, corresponding to which base it mutates to
      my $b = substr($genomes[0], $i, 1); #the DNA base on the progenitor at this location
      my $n = $DNA[($DNAIDX{$b}+$s)%4]; #the mutated DNA base
      substr($genomes[$j], $i, 1) = $n; #set the base on this haplotype to the mutated DNA base
    }
  }
}


## Print the haplotypes
###############################################################################

for (my $j = 0; $j < $PLOIDY; $j++) #for every haplotype
{
  print ">$j\n";
  print $genomes[$j];
  print "\n";
}

my @mutrate = (0)x$PLOIDY;
for (my $j = 0; $j < $PLOIDY; $j++) #for every haplotype
{
  my $mutrate[$j] = sprintf("%0.02f", 100*($nummut[$j] / $newlen));
}

print STDERR "Simulated $newlen total bases, @nummut mutations (@mutrate%), G=$GENOME_LEN, d=$DUP_RATE, r=@HET_RATE, p=$PLOIDY\n";



## Print expected # kmers in each peak
###############################################################################

#print STDERR "\n";
#print STDERR "Model:\n";

#my $delta = int($GENOME_LEN * $DUP_RATE * (1.0-$HET_RATE)**(2*$KMER_LEN));
#my $gamma = int($GENOME_LEN * 2 * $DUP_RATE * (1.0-$HET_RATE)**($KMER_LEN) * (1-(1-$HET_RATE)**($KMER_LEN)));
#my $alpha = int($GENOME_LEN * (2 * (1-$DUP_RATE)*(1.0-(1.0-$HET_RATE)**($KMER_LEN)) + 2*$DUP_RATE*(1-(1-$HET_RATE)**($KMER_LEN))**(2)) + $gamma);
#my $beta  = int($GENOME_LEN * ((1-$DUP_RATE)*((1.0-$HET_RATE)**($KMER_LEN)) + $DUP_RATE * (1-(1-$HET_RATE)**($KMER_LEN))**(2)));

#print STDERR "1 $alpha\n";
#print STDERR "2 $beta\n";
#print STDERR "3 $gamma\n";
#print STDERR "4 $delta\n";
