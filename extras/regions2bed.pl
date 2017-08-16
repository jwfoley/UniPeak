#!/usr/bin/perl

# Read region counts matrix and generate BED track for UCSC browser
# Displays only regions, not their counts

use strict;
use warnings;
use Getopt::Long;
use constant {
  FOR_COL => "0,0,255", # color for forward regions
  REV_COL => "255,0,0", # color for reverse regions
  NONDIR_COL => "191,0,191", # color for nondirectional regions
  TRACK_PRIORITY => 1,
  DEFAULT_NAME => "UniPeak regions"
};

# Print usage if no arguments
if (@ARGV == 0) {die '
Usage: regions2bed.pl [arguments]
Input files must be in region count matrix format

Optional arguments:
  --in <file>       input file (required)
  --out <file>      output file (required)
  --nondirectional  these annotations are not strand-specific (default false)
  --name <str>      track name
  --desc <str>      track description
  
';}

# Read arguments
my $infile;
my $outfile;
my $nondirectional;
my $name;
my $desc;
my $result = GetOptions (
  "in=s" => \$infile,
  "out=s" => \$outfile,
  "nondirectional" => \$nondirectional,
  "name=s" => \$name,
  "desc=s" => \$desc
);
die "error: could not parse command-line arguments\n" unless $result;
die "error: need to specify input file\n" unless $infile;
die "error: need to specify output file\n" unless $outfile;
if (! $name) {$name = DEFAULT_NAME;}
if ($name !~ /^".*"$/) {$name = '"' . $name . '"';}
if ($desc && $desc !~ /^".*"$/) {$desc = '"' . $desc . '"';}

# Generate output
open IN, $infile or die "error reading $infile: $!\n";
open OUT, ">$outfile" or die "error writing $outfile: $!\n";
print OUT "track name=$name ";
if ($desc) {print "description=$desc ";}
print OUT "visibility=\"pack\" priority=", TRACK_PRIORITY, " ", ($nondirectional ? "color=" . NONDIR_COL : "colorByStrand=\"" . FOR_COL . " " . REV_COL . "\""), "\n";
while (<IN>) {
  next if /^[#\t]/;
  chomp;
  if (/^(\w+?):(\d+?)-(\d+?)\t/) {
    my ($chr, $start, $end) = ($1, $2, $3);
    my @out = ($chr, undef, undef, "$chr:$start-$end");
    
    if ($nondirectional) {
      die "error: reverse region in nondirectional analysis\n" if $end < $start;
      @out[1, 2] = ($start - 1, $end);
    } else {
      $out[4] = 0;
      if ($end > $start) {
        @out[1, 2, 5] = ($start - 1, $end, '+');
      } elsif ($end < $start) {
        @out[1, 2, 5] = ($end - 1, $start, '-');
      } else {
        warn "warning: single-base region in line $.  - skipped\n";
        next;
      }
    }
    
    print OUT join("\t", @out), "\n";
  } else {
    warn "warning: bad format in line $.  - skipped\n";
  }
}
close IN;
close OUT;

