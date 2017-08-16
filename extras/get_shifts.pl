#!/usr/bin/perl

# automatically grab the recommended shift values from a list of strand_shift outputs and return them comma-delimited
# if multiple filenames are given together (comma-delimited), get the average of their shifts

use strict;
use warnings;

my @shifts;
foreach my $entry (@ARGV) {
  my @files = split /,/, $entry;
  my $sum;
  foreach my $file (@files) {
    open IN, &expand_tilde($file) or die "error reading $file: $!\n";
    my $success = 0;
    while (<IN>) {
      if (/^# best_shift=(\d+)/) {
        $sum += $1;
        $success = 1;
        last;
      }
    }
    close IN;
    die "error: didn't find best_shift line in $file\n" if (! $success);
  }
  push @shifts, int($sum / @files + 0.5);
}

print join(",", @shifts), "\n";

# http://docstore.mik.ua/orelly/perl/cookbook/ch07_04.htm
# Perl actually doesn't seem to have a better way
sub expand_tilde {
  my $filename = $_[0];
  $filename =~ s{ ^ ~ ( [^/]* ) }
              { $1
                    ? (getpwnam($1))[7]
                    : ( $ENV{HOME} || $ENV{LOGDIR}
                         || (getpwuid($>))[7]
                       )
  }ex;
  return $filename;
}
