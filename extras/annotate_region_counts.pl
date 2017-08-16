#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

# default column numbers (first column is 1)
use constant {
  CONTIG_COL => 3,
  STRAND_COL => 4,
  LEFT_COL => 5,
  RIGHT_COL => 6,
  NAME_COL => 13
};

# Print usage if no arguments
if (@ARGV == 0) {die '
Usage: annotate_region_counts.pl [arguments] annotation_file < infile > outfile

Optional arguments:
  --directional             these annotations are strand-specific
  --max-distance <integer>  maximum distance from feature, default 0
  --header <string>         name for new column
  
Column indices for annotation file (which column contains what):
  --name-col <positive integer>    feature name
  --contig-col <positive integer>  contig
  --strand-col <positive integer>  strand (+ or -), only used if directional
  --left-col <positive integer>    "left" end (start position on + strand)
  --right-col <positive integer>   "right" end (end position on + strand)
Note: the first column is number 1

';}

my ($directional, $max_dist, $annot_col_name, $name_col, $contig_col, $strand_col, $left_col, $right_col) = (0, 0, undef, NAME_COL, CONTIG_COL, STRAND_COL, LEFT_COL, RIGHT_COL);
GetOptions(
  "directional" => \$directional,
  "max-distance=i" => \$max_dist,
  "header=s" => \$annot_col_name,
  "name-col=i" => \$name_col,
  "contig-col=i" => \$contig_col,
  "strand-col=i" => \$strand_col,
  "left-col=i" => \$left_col,
  "right-col=i" => \$right_col
) or die "error: bad arguments\n";
my $annot_table_fname = shift @ARGV;
if (! $annot_table_fname) {die "error: must specify annotation filename\n";}

# read annotation file
my %annot;
open ANNOT, $annot_table_fname or die "error: could not read $annot_table_fname: $!\n";
$_ = <ANNOT>;
while (<ANNOT>) {
  next if (length == 1 || /^#/);
  chomp;
  my ($chr, $strand, $left, $right, $name) = (split /\t/)[$contig_col - 1, $strand_col - 1, $left_col - 1, $right_col - 1, $name_col - 1];
  if (! ($chr && $strand && $left && $right && $name)) {
    warn "warning: bad format in $annot_table_fname line $. - skipped\n";
    next;
  }
  if (! $directional) {$strand = '+';}
  
  if (! $annot{$chr}{$strand}{$name}) {
    $annot{$chr}{$strand}{$name}{left} = $left;
    $annot{$chr}{$strand}{$name}{right} = $right;
  } else {
    if ($left < $annot{$chr}{$strand}{$name}{left}) {
      $annot{$chr}{$strand}{$name}{left} = $left;
    }
    if ($right > $annot{$chr}{$strand}{$name}{right}) {
      $annot{$chr}{$strand}{$name}{right} = $right;
    }
  }
}
close ANNOT;

# header
print "# annotation_file=$annot_table_fname
# directional=$directional
# max_distance=$max_dist

";
while (<STDIN>) { 
  if (length == 1 || /^#/) {
    print $_;
  } else {
    print "\t$annot_col_name$_";
    last;
  }
}

# output
while (<STDIN>) {
  next if (length == 1 || /^#/);
  chomp;
  my @line = split /\t/;
  my $name = '';
  if ($line[0] =~ /(\w+):(\d+)-(\d+)/) {
    my $chr = $1;
    my ($strand, $left, $right);
    my @names;
    if ($2 < $3) {
      ($strand, $left, $right) = ('+', $2, $3);
    } elsif ($2 > $3) {
      ($strand, $left, $right) = ('-', $3, $2);
    } else {
      warn "warning: bad coordinates in line $. - removed\n";
      next;
    }
    if (! $directional) {$strand = '+';}
    foreach my $name (keys %{$annot{$chr}{$strand}}) {
      my $annot_start = $annot{$chr}{$strand}{$name}{left};
      my $annot_end = $annot{$chr}{$strand}{$name}{right};
      if (($left > $annot_start - $max_dist && $left <= $annot_end + $max_dist) || # left coord is inside region
        ($right > $annot_start - $max_dist && $right <= $annot_end + $max_dist) || # right coord is inside
        ($left < $annot_start - $max_dist && $right > $annot_end + $max_dist)) { # entire region is inside feature
        push @names, $name;
      }
    }
    if (@names) {
      $name = join ",", @names;
    }
  } else {
    warn "warning: bad format in line $. - skipped\n";
    next;
  }
  print join("\t", $line[0], $name, @line[1..$#line]), "\n";
}

