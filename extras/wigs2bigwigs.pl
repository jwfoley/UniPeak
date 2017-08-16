#!/usr/bin/perl

# Convert UniPeak wiggle files into UCSC bigWig files
# Input wiggle files may be compressed with gzip or bzip2 (UniPeak options)
# Must put separate strands in separate files due to limitation of bigWig format
# Prints track headers to load in UCSC browser

use strict;
use warnings;
use Getopt::Long;
use constant {
  DEFAULT_SIZES_PATH => "somepath/somefile.chrom.sizes", # default path to chrom.sizes file
  DEFAULT_EXECUTABLE_PATH => "~/bin/wigToBigWig -clip", # default path to executable and any flags you prefer 
  DEFAULT_PREFIX => "http://www.stanford.edu/~yourusername/" # default prefix for data URLs (not used inside bigWig files)
};

my $gzip_regex = qr/.gz$/;
my $bzip2_regex = qr/.bz2$/;

# Print usage if no arguments
if (@ARGV == 0) {die "
Usage: $0 [arguments] file1.wig file2.wig file3.wig.gz ...

Optional arguments:
  --prefix <filename prefix>  prefix for URLs
  --sizes <filename>          path to chrom.sizes file
  --executable <filename>     path to executable
";}

# Read arguments
my ($prefix, $sizes_path, $executable_path) = (DEFAULT_PREFIX, DEFAULT_SIZES_PATH, DEFAULT_EXECUTABLE_PATH);
my $result = GetOptions (
  "prefix=s" => \$prefix,
  "sizes=s" => \$sizes_path,
  "executable=s" => \$executable_path,
);
if (! $result) {die "error: bad arguments\n";}
my @wig_files = @ARGV;
die "error: no input files\n" if (@wig_files == 0);

# Parse wig files
foreach my $infile (@wig_files) {
  my $out_root = $infile;
  $out_root =~ s/$gzip_regex$//;
  $out_root =~ s/$bzip2_regex$//;
  $out_root =~ s/.wig$//;
  open IN, ($infile =~ /$gzip_regex/ ? "gunzip < $infile |" : ($infile =~ /$bzip2_regex/ ? "bunzip2 < $infile |" : $infile)) or die "error reading $infile: $!\n";
  my ($out, $outfile);
  while (<IN>) {
    if (/^#/) { # comment line
      next;
      
    } elsif (/^track/) { # track header
      chomp (my $header = $_);
      if ($out) {close $out;}
      my $strand = "";
      if ($header =~ /name=".+([+-])"/) {
        $strand = $1;
      }
      $outfile = "$out_root$strand.bw";
      open $out, "| $executable_path stdin $sizes_path $outfile" or die "error: $!\n";
      
      $header =~ s/wiggle_0/bigWig/;
      print "$header bigDataUrl=$prefix$outfile\n";

    } elsif (length $_ > 0) { # profile line
      if ($out) {
        print $out $_;
      } else {
        die "error: no track defined at line $.\n";
      }
    }
  }
  if ($out) {close $out;}
  close IN;
}

