#!/usr/bin/perl -w
use strict;

my $filename = shift @ARGV or die "eval_err.pl simhisto.out\n";

open FILE, $filename or die "Cant open $filename ($!)\n";

my @model = (0,0,0,0);
my @obs   = (0,0,0,0);

my $inobs = 0;

while(<FILE>)
{
  chomp;

  if (/^Model:/)
  {
    # The next 4 lines are the model output
    my $aline = <FILE>; chomp $aline; $model[0] = (split/ /, $aline)[1];
    my $bline = <FILE>; chomp $bline; $model[1] = (split/ /, $bline)[1];
    my $cline = <FILE>; chomp $cline; $model[2] = (split/ /, $cline)[1];
    my $dline = <FILE>; chomp $dline; $model[3] = (split/ /, $dline)[1];
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

print "model:\t$model[0]\t$model[1]\t$model[2]\t$model[3]\n";
print "obs:\t$obs[0]\t$obs[1]\t$obs[2]\t$obs[3]\n";

my @diff  = (0,0,0,0); # count of differences
my @diffp = (0,0,0,0); # percent of differences

my $sumkmer = 0;
my $sumdiff = 0;

for (my $i = 0; $i < 4; $i++)
{
  $diff[$i] = $model[$i] - $obs[$i];
  $diffp[$i] = sprintf("%0.02f", ($obs[$i]) ? (100*$diff[$i] / $obs[$i]) : 0);
  $sumdiff += $diff[$i];
  $sumkmer += $obs[$i];
}

my $score = sprintf("%.03f", 100.0*(1-($sumdiff / $sumkmer)));

print "diff:\t$diff[0]\t$diff[1]\t$diff[2]\t$diff[3]\n";
print "diffp:\t$diffp[0]\t$diffp[1]\t$diffp[2]\t$diffp[3]\n";
print "#$filename\tsumkmer\t$sumkmer\tsumdiff:\t$sumdiff\tscore:\t$score\n";
