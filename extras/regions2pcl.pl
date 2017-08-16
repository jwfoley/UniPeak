#!/usr/bin/perl

# Read region counts matrix and generate PCL file for microarray software

# usage: regions2pcl.pl < input_filename > output_filename
# Input files must be in region count matrix format, optionally with annotation column (first column after coordinates; no columns should come between kurtosis and the first experiment)
# Transformation of count values may be performed before this step as long as the format is correct

use strict;
use warnings;

my $kurtosis_col = 0;
my $annotated = 0;
my $cols = 0;
while (<STDIN>) {
  next if (/^\n/ || /^#/);
  chomp;
  my @row = split /\t/;
  if ($row[0] eq "") { # headers
  
    # set column number
    if ($cols == 0) {
      $cols = @row;
    } else {
      die "error: extra header row in line $.\n";
    }
    
    # find kurtosis
    for (my $i = 1; $i < @row; $i++) {
      if ($row[$i] eq "kurtosis") {
        if ($kurtosis_col == 0) {
          $kurtosis_col = $i;
        } else {
          die "error: kurtosis defined twice in line $.\n";
        }
      }
    }
    if ($kurtosis_col > 1) {$annotated = 1;}
  
    # make PCL header
    print join("\t", "region", ($annotated == 1 ? $row[1] : ""), "GWEIGHT", @row[($kurtosis_col + 1)..$#row]), "\n";
    print join("\t", "EWEIGHT", "", "");
    for (my $i = 0; $i < $#row - $kurtosis_col; $i++) {print "\t1";}
    print "\n";

  } else { # data row
    if ($cols == 0) {die "error: no header found before line $.\n";}
    if (@row != $cols) {die "error: wrong number of rows in line $.\n";}
    print join("\t", $row[0], ($annotated == 1 ? $row[1] : ""), 1, @row[($kurtosis_col + 1)..$#row]), "\n";
  }
}
