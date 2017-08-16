#!/usr/bin/perl

use strict;
use warnings;
use Cwd;

my ($compress, $compress_string) = (1, ".gz"); # always gzip wiggle files; user can edit shell script if there is really some reason not to

# intro
print "
UniPeak script generator
This program automatically generates a Unix shell script to run UniPeak.
If a default value is given (in parentheses), leave blank to use it.
To save time, you can list your alignment files as command-line arguments.\n\n";

# check ARGV; can't use it if any of them are bad
foreach my $file (@ARGV) {die "error: couldn't find $file\n" if (! -e &expand_tilde($file));}

# get parameters
print "\n";
my $unipeak_dir;
while () {
  $unipeak_dir = &get_value("UniPeak home directory", &Cwd::abs_path(substr($0, 0, rindex($0, '/'))) . '/');
  if (-e &expand_tilde($unipeak_dir)) {last;}
  print "error: couldn't find directory $unipeak_dir\n";
}
if ((substr $unipeak_dir, -1, 1) ne '/') {$unipeak_dir .= '/';}

print "\n";
my $output_dir = &get_value("output directory", getcwd() . '/');
if ((substr $output_dir, -1, 1) ne '/') {$output_dir .= '/';}

print "\n";
my $contig_table;
while () {
  $contig_table = &get_value("contig table path", $unipeak_dir . "contig/hg19.txt");
  if (-e &expand_tilde($contig_table)) {last;}
  print "error: couldn't find file $contig_table\n";
}

print "\nNote: ChIP-seq and RNA-seq are generally not directional; 3SEQ is.\n";
my $directional = &get_bool("Are these experiments directional?", 0);

my $strand_shift = ! $directional;
if (! $directional) {$strand_shift = &get_bool("  Correct for strand shift?", $strand_shift);}

my $uses_background = &get_bool("\nAre any of these samples \"background\" (e.g. input, IgG)?", ! $directional);

#print "\nNote: you may skip this step if all alignment files are sorted\n"; # disabled till this is fixed
#my $make_wig = &get_bool("Convert alignments to wiggle format?", 1);
my $make_wig = 1;

#my $compress = 0;
#if ($make_wig) {
#  print "  "; # indent
#  $compress = &get_bool("Compress them?", $compress);
#}
#print "\n";
#my $compress = &get_bool("Compress wiggle files?", 1);
#my $compress_string = ($compress ? ".gz" : "");

print "\nNote: files may be in Eland, Corona Lite, BED, SAM, BAM, or UniPeak wiggle format;\nthey may be compressed with gzip or bzip2 if the file extension is correct\n";
my @align_files;
my @expt_names;
my @background;
if (@ARGV) {
  # read alignment files from command line
  foreach my $file (@ARGV) {
    if ($file =~ / /) {$file = &enquote($file);}
    print "$file\n";
    my $is_background = ($uses_background ? &get_bool("  Is this a \"background\" sample?", 0) : 0);
    push @align_files, $file;
    push @background, $is_background;
    print "  ";
    push @expt_names, &enquote(&get_value("optional name for this sample"));
  }
} else {
  # enter alignment files manually
  while () {
    my $file = &get_value("an alignment filename or leave blank if done");
    if ($file =~ / /) {$file = &enquote($file);}
        
    if ($file ne "") {
      if (-e &expand_tilde($file)) {
        my $is_background = ($uses_background ? &get_bool("  Is this a \"background\" sample?", 0) : 0);
        push @align_files, $file;
        push @background, $is_background;
        print "  ";
        push @expt_names, &enquote(&get_value("optional name for this sample"));
      } else {
        print "error: couldn't find file $file\n";
      }
    } else {
      if (@align_files == 0) {
        print "error: you must enter at least one alignment file\n";
      } else {
        my $has_non_background = 0;
        foreach my $state (@background) {
          if (! $state) {
            $has_non_background = 1;
            last;
          }
        }
        if (! $has_non_background) {
          print "error: you must enter at least one sample that is not \"background\"\n";
        } else {        
          last; # exit point
        }
      }
    }
  }
}
my @background_list;
for (my $i = 0; $i < @background; $i++) {
  if ($background[$i]) {push @background_list, $i + 1;}
}

print "\n";
my $use_defaults = &get_bool("Use default parameters?", 1);

my ($posterior_prob, $offset, $use_length, $mismatch_tolerance, $bandwidth, $region_threshold, $kurtosis_threshold, $hit_threshold, $mappable, $max_shift, $test_regions, $stringent_kurtosis_threshold);
my ($parse_args, $strand_shift_args, $regions_args) = ("", "", "");
if (! $use_defaults) {
  print "Custom parameters (see README.TXT for details, leave blank for default):\n";

  print "  ";
  my $posterior_prob = &get_value("alignment posterior probability threshold");
  if ($posterior_prob ne "") {
    die "bad posterior probability threshold\n" unless ($posterior_prob >= 0 && $posterior_prob <= 1);
    $parse_args .= " -p $posterior_prob";
  }
  
  print "  ";
  my $offset = &get_value("position offset, e.g. barcode length");
  if ($offset ne "") {$parse_args .= " -s $offset";}
  
  print "  ";
  my $use_length = &get_value("sequence length to use, if not entire sequence");
  if ($use_length ne "") {
    die "bad sequence length\n" if ($use_length < 0);
    $parse_args .= " -l $use_length";
  }
  
  print "  ";
  my $mismatch_tolerance = &get_value("mismatch tolerance");
  if ($mismatch_tolerance ne "") {
    die "bad mismatch tolerance\n" if ($mismatch_tolerance < 0);
    $parse_args .= " -i $mismatch_tolerance";
  }

  # copy parse args to later executables if not using convert_align
  if (! $make_wig) {
    $strand_shift_args = $parse_args;
    $regions_args = $parse_args;
  }
  
  print "  ";
  my $bandwidth = &get_value("KDE bandwidth");
  if ($bandwidth ne "") {
    die "bad bandwidth\n" unless ($bandwidth > 0);
    $regions_args .= " -b $bandwidth";
  }
  
  print "  ";
  my $region_threshold = &get_value("minimum fold-enrichment threshold for regions");
  if ($region_threshold ne "") {  
    die "bad region threshold\n" unless ($region_threshold > 0);
    $regions_args .= " -r $region_threshold";
  }
  
  print "  ";
  my $hit_threshold = &get_value("minimum average hits per experiment in a region");
  if ($hit_threshold ne "") {
    die "bad hit threshold\n" unless ($hit_threshold >= 0);
    $regions_args .= " -t $hit_threshold";
  }
  
  print "  ";
  my $mappable = &get_value("number of mappable positions, if not total of contig size");
  if ($mappable ne "") {
    die "bad mappable positions\n" unless ($mappable > 0);
    $regions_args .= " -m $mappable";
  }

  $strand_shift_args = $regions_args; # do this before kurtosis since that might be different
  print "  ";
  my $kurtosis_threshold = &get_value("maximum kurtosis threshold for output regions");
  if ($kurtosis_threshold ne "") {
    die "bad kurtosis threshold\n" unless ($kurtosis_threshold > 0);
    $regions_args .= " -k $kurtosis_threshold";
  }
  
  if ($strand_shift) {
    print "  ";
    my $stringent_kurtosis_threshold = &get_value("kurtosis threshold for strand shift calibration", $kurtosis_threshold);
    if (defined $stringent_kurtosis_threshold && $stringent_kurtosis_threshold ne "") {
      die "bad kurtosis threshold\n" unless ($stringent_kurtosis_threshold > 0);
      $strand_shift_args .= " -k $stringent_kurtosis_threshold";
    }

    print "  ";
    my $max_shift = &get_value("maximum strand shift to try");
    if ($max_shift ne "") {
      die "bad shift\n" unless ($max_shift > 0);
      $strand_shift_args .= " -x $max_shift";
    }
    
    print "  ";
    my $test_regions = &get_value("number of regions to use for recalibrating strand shift");
    if ($test_regions ne "") {
      die "bad region number\n" unless ($test_regions > 0);
      $strand_shift_args .= " -g $test_regions";
    }
  }  
}


print "\n";
my $output_script = $output_dir . "unipeak.sh";
get_value("name of script to generate", $output_script);

# report parameters
print "\n
PARAMETERS
UniPeak home directory: $unipeak_dir
Output directory: $output_dir
Alignment files:";
if ($uses_background) {
  for (my $i = 0; $i < @align_files; $i++) {
    print "\n  ", $background[$i] ? "Background" : "Sample", ": $align_files[$i]";
  }
} else {
  print "\n  ";
  print join("\n  ", @align_files);
}
print "
Contig table: $contig_table
" . ($make_wig ? "Make " . ($compress ? "compressed " : "") . "wiggle files\n" : "")
. ($directional ? "Directional" : "Non-directional" . ($strand_shift ? ", calibrate strand shift" : "")) . "\n";
if (! $use_defaults) {
  print "Additional arguments:\n";
  if ($make_wig && $parse_args ne "") {print "  convert_align$parse_args\n";}
  if ($strand_shift && $strand_shift_args ne "") {print "  strand_shift$strand_shift_args\n";}
  if ($regions_args ne "") {print "  regions$regions_args\n";}
}
print "\n";

# write output
if (! -e &expand_tilde($output_dir)) {mkdir &expand_tilde($output_dir) or die "error: couldn't create $output_dir: $!\n";}
$contig_table = &enquote($contig_table);
open OUT, ">" . &expand_tilde($output_script) or die "error writing $output_script  $!\n";
print OUT "#!/bin/bash\nset -e\n";

my @sorted_files;
if ($make_wig) {
  print OUT "\n# generate wiggle files\n";
  for (my $i = 0; $i < @align_files; $i++) {
    push @sorted_files, &enquote($output_dir . &get_root_name($align_files[$i]) . ".wig" . $compress_string);
    print OUT &enquote($unipeak_dir . "bin/convert_align") . "$parse_args -c $contig_table -o $sorted_files[-1]" . ($expt_names[$i] ne "" ? " -n $expt_names[$i] " : " ") . &enquote($align_files[$i]) . "\n";
  }
} else {
  @sorted_files = @align_files;
  foreach my $file (@sorted_files) {$file = &enquote($file);}
}

my @shift_files;
if ($strand_shift) {
  print OUT "\n# calibrate strand shift\n";
  foreach my $sorted_file (@sorted_files) {
    push @shift_files, &enquote($output_dir . &get_root_name($sorted_file) . "_shift.txt");
    print OUT &enquote($unipeak_dir . "bin/strand_shift") . "$strand_shift_args -c $contig_table -o $shift_files[-1] $sorted_file\n";
  }
}

print OUT "\n# run region caller\n";
print OUT &enquote($unipeak_dir . "bin/regions") . $regions_args . ($directional ? "" : " -D") . " -c $contig_table" . (@background_list ? " -e " . join(",", @background_list) : "") . ($directional ? "" : ($strand_shift ? " -s `" . &enquote($unipeak_dir . "extras/get_shifts.pl") . " " . join(" ", @shift_files) . "`": "")) . " -o " . &enquote($output_dir . "regions.txt") . " " . join(" ", @sorted_files) . "\n";

close OUT;
chmod 0755, &expand_tilde($output_script); # make executable
if (-e "/usr/bin/SetFile") {system("SetFile", "-t", "APPL", &expand_tilde($output_script));} # make clickable in OS X, I think?
print "
Successfully wrote $output_script
Execute it to run UniPeak!\n\n";


sub get_value {
  my ($phrase, $default) = @_;
  print "Enter $phrase" . (defined $default && $default ne "" ? " ($default)" : "") . ": ";
  chomp ($_ = <STDIN>);
  return ($_ eq "" && defined $default ? $default : $_);
}

sub get_bool {
  my ($phrase, $default) = @_;
  while () {
    print $phrase . (defined $default ? " (" . ($default ? "Y" : "N") . ")" : "") . ": ";
    chomp ($_ = <STDIN>);
    if ($_ eq "") {
      return $default;
    } else {
      if (/^n/i) {
        return 0;
      } elsif (/^y/i) {
        return 1;
      } else {
        print "please enter Y or N\n";
      }
    }
  }
}  

sub get_root_name {
  my $path = $_[0];
  my $last_slash = rindex($path, '/');
  my $first_dot = index($path, '.', $last_slash);
  if ($first_dot == -1) {$first_dot = length $path;}
  return substr($path, $last_slash + 1, $first_dot - $last_slash - 1);
}

sub enquote {
  my $string = $_[0];
  if ($string =~ / / && $string !~ /^".*"$/) {$string = "\"$string\"";}
  return $string;
}

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
