#!/usr/bin/perl -w
use strict;

my $filename = shift @ARGV or die "eval_err.pl simhisto.out\n";

open FILE, $filename or die "Cant open $filename ($!)\n";

my $start = rindex($filename, "p")+1;
my $end = rindex($filename, "-");
my $length = $end - $start;
my $PlOIDY = int(substr($filename, $start, $length));
my @model = (0)x(2*$PLOIDY);
my @obs   = (0)x(2*$PLoIDY);

my $inobs = 0;

while(<FILE>)
{
  chomp;

  if (/^Model:/)
  {
    # The next 2p lines are the model output
    my @lines = ()x(2*$PLOIDY);
    for (my $i = 0; $i < 2*$PLOIDY; $i++)
    {
      $lines[$i] = <FILE>; chomp $lines[$i]; $model[$i] = (split/ /, $lines[$i])[1];
    }
  }
  elsif (/^Observed Counts:/)
  {
    # The next few lines are the jellyfish output
    $inobs = 1;
  }
  elsif ($inobs)
  {
    my ($f, $c) = split / /, $_;
    $obs[$f-1] = $c;
  }
}

print "model:\t".join("\t", @model)."\n";
print "obs:\t".join("\t", @obs)."\n";

my @diff = (0)x(2*$PLOIDY); # count of differences
my @diffp = (0)x(2*$PLOIDY); # percent of differences

my $sumkmer = 0;
my $sumdiff = 0;

for (my $i = 0; $i < 2*$PLOIDY; $i++)
{
  $diff[$i] = abs($model[$i] - $obs[$i]);
  $diffp[$i] = sprintf("%0.02f", ($obs[$i]) ? (100*$diff[$i] / $obs[$i]) : 0);
  $sumdiff += $diff[$i];
  $sumkmer += $obs[$i];
}

my $score = sprintf("%.03f", 100.0*(1-($sumdiff / $sumkmer)));

print "diff:\t".join("\t", @diff)."\n";
print "diffp:\t".join("\t", @diffp)."\n";
print "#$filename\tsumkmer\t$sumkmer\tsumdiff:\t$sumdiff\tscore:\t$score\n";
