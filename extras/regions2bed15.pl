#!/usr/bin/perl

# Read region counts matrix and generate BED15 "microarray" track for UCSC browser

use strict;
use warnings;
use Getopt::Long;

use constant DEFAULT_SCALE => 1000; # highest value of color spectrum
use constant STEPS => 5; # number of colors in color spectrum
use constant TRACK_PRIORITY => 1;

# Print usage if no arguments
if (@ARGV == 0) {die '
Usage: regions2bed15.pl [arguments]
Input files must be in region count matrix format

Arguments:
  --input <str>   input file name (required)
  --output <str>  output file name (required)
  --name <str>    track name (optional)
  --desc <str>    track description (optional)
  --scale <num>   color scale maximum (optional, default 1000)
  
';}

# Read arguments
my $infile;
my $outfile;
my $name;
my $desc;
my $scale = DEFAULT_SCALE;
my $result = GetOptions (
  "input=s" => \$infile,
  "output=s" => \$outfile,
  "name=s" => \$name,
  "desc=s" => \$desc,
  "scale=f" => \$scale
);
die "error: could not parse command-line arguments\n" if (! $result);
die "error: need input filename\n" if (! $infile);
die "error: need output filename\n" if (! $outfile);
if (! $name) {
  $name = $outfile;
  if (my $last_slash = rindex($name, '/')) {$name = substr($name, $last_slash + 1, (length $name) - ($last_slash + 1));}
  my $first_dot = index($name, '.');
  if ($first_dot >= 0) {$name = substr($name, 0, $first_dot);}
}
my $step = $scale / STEPS;

# Write header
open IN, $infile or die "error: could not read $infile: $!\n";
open OUT, ">$outfile" or die "error: could not write $outfile: $!\n";
print OUT "# regions_file=$infile\n";
my @expts;
my $kurtosis_column = 0;
while (<IN>) {
  chomp;
  if (length $_ == 0) {
    next;
  } elsif (/^#/) {
    print OUT "$_\n";
  } elsif (/^\t/) {
    tr/",//d; # wait, what is this?
    @expts = split /\t/;
    for (my $i = 1; $i < @expts; $i++) {
      if ($expts[$i] eq "kurtosis") {
        if ($kurtosis_column) {die "error: kurtosis found twice\n";}
        $kurtosis_column = $i + 1;
      }
    }
    last;
  } else {
    die "error: did not find column headers\n";
  }
}
splice(@expts, 0, $kurtosis_column);   
my $num_expts = @expts;
print "found $num_expts experiments\n";
print OUT "track name=\"$name +\" ";
if ($desc) {print OUT "description=$desc ";}
print OUT "type=array expNames=\"", join(",", @expts), "\" expScale=$scale expStep=$step priority=TRACK_PRIORITY visibility=full\n";

# Read regions
my $reverse;
my $regions = 0;
while (<IN>) {
  chomp;
  next if (length $_ == 0);
  my @fields = split /\t/;
  if (@fields != @expts + $kurtosis_column) {die "error: wrong number of fields in line $.\n";}
  my $region = $fields[0];
  if ($region =~ /^(\w+?):(\d+?)-(\d+?)$/) {
    my ($chr, $start, $end) = ($1, $2, $3);
    if (! $reverse && $end < $start) {
      print OUT "track name=\"$name -\" ";
      if ($desc) {print OUT "description=$desc ";}
      print OUT "type=array expNames=\"", join(",", @expts), "\" expScale=$scale expStep=$step priority=TRACK_PRIORITY visibility=full\n";
      $reverse = 1;
    }
    splice(@fields, 0, $kurtosis_column - 1);
    my $kurtosis = shift @fields;
    my @out = ($chr, undef, undef, "k=" . $kurtosis, 0, undef, 0, 0, 0, 1, abs($start - $end) + 1, 0, $num_expts, join(",", 0..$#expts), join(",", @fields));
    if (! $reverse) {
      @out[1, 2, 5] = ($start - 1, $end, '+');
    } else {
      @out[1, 2, 5] = ($end - 1, $start, '-');
    }
    print OUT join("\t", @out), "\n";
    $regions++;
  } else {
    warn "warning: bad format in line $.  - skipped\n";
  }
}
close IN;
close OUT;
print "found $regions regions\n",
