#!/usr/bin/perl -w
use strict;

my $USAGE = "simgenome.pl GENOME_LEN DUP_RATE HET_RATE1 HET_RATE2 HET_RATE_ACROSS KMER_LEN > genome.fa\n";

die $USAGE if (scalar @ARGV != 6);

my $GENOME_LEN = shift @ARGV or die $USAGE;
my $DUP_RATE   = shift @ARGV;
my $HET_RATE1  = shift @ARGV;
my $HET_RATE2  = shift @ARGV;
my $HET_RATE_ACROSS = shift @ARGV;
my $KMER_LEN   = shift @ARGV;

die "Dup rate cant be > 1\n" if ($DUP_RATE > 1.0);
die "Hetrate1 cant be > 1\n" if ($HET_RATE1 > 1.0);
die "Hetrate2 cant be > 1\n" if ($HET_RATE2 > 1.0);
die "HetrateZ cant be > 1\n" if ($HET_RATE_ACROSS > 1.0);

print STDERR "Simulating G=$GENOME_LEN, d=$DUP_RATE, r1=$HET_RATE1, r2=$HET_RATE2, rZ=$HET_RATE_ACROSS k=$KMER_LEN\n";



## Initialize DNA coding table: 0<->A, 1<->C, 2<->G, 3<->T
###############################################################################
my @DNA = ("A", "C", "G", "T");
my %DNAIDX;
for (my $i = 0; $i < scalar @DNA; $i++) { $DNAIDX{$DNA[$i]} = $i; }


## Function for mutating a genome string
###############################################################################
sub mutate
{
  my $orig    = shift;
  my $hetrate = shift;

  my $new = $orig;
  my $l = length($new);
  
  for (my $i = 0; $i < $l; $i++)
  {
    my $r = rand();
    if ($r <= $hetrate)
    {
      my $b = substr($new, $i, 1);
      my $n = $DNA[($DNAIDX{$b}+1)%4];
      substr($new, $i, 1) = $n;
    }
  }

  return $new;
}


## Simutate a random genome of length GENOME_LEN
###############################################################################

print STDERR " ... simulating progenitor\n";

my $matgenome1 = "";

for (my $i = 0; $i < $GENOME_LEN; $i++)
{
  my $z = int(rand(4));
  my $b = $DNA[$z];
  $matgenome1 .= $b;
}


## Now duplicate the first $DUP_RATE percent of the genome
###############################################################################

my $dup = substr($matgenome1, 0, $DUP_RATE * $GENOME_LEN);
$matgenome1 .= $dup;
my $newlen = length($matgenome1);

print STDERR " ... duplicated $DUP_RATE, newlen=$newlen\n";



## Now create the heterozygous alleles
###############################################################################

print STDERR "generating paternal genome1 (r1=$HET_RATE1)\n";
my $patgenome1 = mutate($matgenome1, $HET_RATE1);

print STDERR "generating maternal genome2 (rZ=$HET_RATE_ACROSS)\n";
my $matgenome2 = mutate($matgenome1, $HET_RATE_ACROSS);

print STDERR "generating paternal genome2 (r2=$HET_RATE2)\n";
my $patgenome2 = mutate($matgenome2, $HET_RATE2);




## Print the genome sequences
###############################################################################

print ">mat1\n";
print $matgenome1;
print "\n";

print ">pat1\n";
print $patgenome1;
print "\n";

print ">mat2\n";
print $matgenome2;
print "\n";

print ">pat2\n";
print $patgenome2;
print "\n";






# ## Print expected # kmers in each peak
# ###############################################################################
# 
# print STDERR "\n";
# print STDERR "Model:\n";
# 
# my $delta = int($GENOME_LEN * $DUP_RATE * (1.0-$HET_RATE)**(2*$KMER_LEN));
# my $gamma = int($GENOME_LEN * 2 * $DUP_RATE * (1.0-$HET_RATE)**($KMER_LEN) * (1-(1-$HET_RATE)**($KMER_LEN)));
# my $alpha = int($GENOME_LEN * (2 * (1-$DUP_RATE)*(1.0-(1.0-$HET_RATE)**($KMER_LEN)) + 2*$DUP_RATE*(1-(1-$HET_RATE)**($KMER_LEN))**(2)) + $gamma);
# my $beta  = int($GENOME_LEN * ((1-$DUP_RATE)*((1.0-$HET_RATE)**($KMER_LEN)) + $DUP_RATE * (1-(1-$HET_RATE)**($KMER_LEN))**(2)));
# 
# print STDERR "1 $alpha\n";
# print STDERR "2 $beta\n";
# print STDERR "3 $gamma\n";
# print STDERR "4 $delta\n";
