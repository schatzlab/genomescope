#!/usr/bin/perl -w
use strict;

my $USAGE = "simgenome.pl GENOME_LEN DUP_RATE HET_RATE KMER_LEN > genome.fa\n";

die $USAGE if (scalar @ARGV != 4);

my $GENOME_LEN = shift @ARGV or die $USAGE;
my $DUP_RATE   = shift @ARGV;
my $HET_RATE   = shift @ARGV;
my $KMER_LEN   = shift @ARGV;

die "Dup rate cant be > 1\n" if ($DUP_RATE > 1.0);
die "Het rate cant be > 1\n" if ($HET_RATE > 1.0);

print STDERR "Simulating G=$GENOME_LEN, d=$DUP_RATE, r=$HET_RATE, k=$KMER_LEN\n";



## Initialize DNA coding table: 0<->A, 1<->C, 2<->G, 3<->T
###############################################################################
my @DNA = ("A", "C", "G", "T");
my %DNAIDX;
for (my $i = 0; $i < scalar @DNA; $i++) { $DNAIDX{$DNA[$i]} = $i; }


## Simutate a random genome of length GENOME_LEN
###############################################################################

print STDERR " ... simulating progenitor\n";

my $matgenome = "";

for (my $i = 0; $i < $GENOME_LEN; $i++)
{
  my $z = int(rand(4));
  my $b = $DNA[$z];
  $matgenome .= $b;
}


## Now duplicate the first $DUP_RATE percent of the genome
###############################################################################

my $dup = substr($matgenome, 0, $DUP_RATE * $GENOME_LEN);
$matgenome .= $dup;
my $newlen = length($matgenome);

print STDERR " ... duplicated $DUP_RATE, newlen=$newlen\n";

## Print the maternal haplotype
###############################################################################

print ">mat\n";
print $matgenome;
print "\n";


## Add random heterozygosity to the paternal haplotype
###############################################################################

print STDERR " ... add heterozygosity r=$HET_RATE\n";

my $nummut = 0;
my $patgenome = $matgenome;

for (my $i = 0; $i < $newlen; $i++)
{
  my $r = rand();
  if ($r <= $HET_RATE)
  {
    $nummut++;
    my $b = substr($patgenome, $i, 1);
    my $n = $DNA[($DNAIDX{$b}+1)%4];
    substr($patgenome, $i, 1) = $n;
  }
}


## Print the paternal haplotype
###############################################################################

print ">pat\n";
print $patgenome;
print "\n";

my $mutrate = sprintf("%0.02f", 100*($nummut / $newlen));
print STDERR "Simulated $newlen total bases, $nummut mutations ($mutrate%), G=$GENOME_LEN, d=$DUP_RATE, r=$HET_RATE\n";



## Print expected # kmers in each peak
###############################################################################

print STDERR "\n";
print STDERR "Model:\n";

my $delta = int($GENOME_LEN * $DUP_RATE * (1.0-$HET_RATE)**(2*$KMER_LEN));
my $gamma = int($GENOME_LEN * 2 * $DUP_RATE * (1.0-$HET_RATE)**($KMER_LEN) * (1-(1-$HET_RATE)**($KMER_LEN)));
my $alpha = int($GENOME_LEN * (2 * (1-$DUP_RATE)*(1.0-(1.0-$HET_RATE)**($KMER_LEN)) + 2*$DUP_RATE*(1-(1-$HET_RATE)**($KMER_LEN))**(2)) + $gamma);
my $beta  = int($GENOME_LEN * ((1-$DUP_RATE)*((1.0-$HET_RATE)**($KMER_LEN)) + $DUP_RATE * (1-(1-$HET_RATE)**($KMER_LEN))**(2)));

print STDERR "1 $alpha\n";
print STDERR "2 $beta\n";
print STDERR "3 $gamma\n";
print STDERR "4 $delta\n";
