#!/usr/bin/perl

# make a contig sizes file from FASTA reference sequence
# usage:  fasta2contigs.pl < in.fasta > out.txt
#         cat in1.fasta in2.fasta in3.fasta | fasta2contigs.pl > out.txt

use strict;
use warnings;

my $contig = "contig";
my $size = 0;
while (<STDIN>) {
  chomp;
  if (/^>\s*(.+)$/) {
    # new contig
    if ($size) {
      print "$contig\t$size\n";
      $size = 0;  
    }
    $contig = $1;
  } else {
    die "error: bad character in line $.\n" if (/[^A-Za-z]/);
    $size += length;
  }
}
if ($size) {print "$contig\t$size\n";}

