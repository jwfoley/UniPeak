#ifndef DEFAULTS_H
#define DEFAULTS_H

#include <cstdlib>
#include <string>
#include <cmath>

#include <boost/cstdint.hpp>

namespace unipeak {

  using namespace std;
  using namespace boost;

  const string VERSION = "1.0";

  // Data types
  typedef boost::uint16_t UShort;
  typedef boost::uint32_t ContigNo; // incomplete genomes have lots of contigs
  typedef boost::uint32_t Pos;
  typedef boost::uint32_t HitCount; // no memory savings if this is narrower than Pos
  typedef boost::uint64_t Counter;
  typedef double Score;

  // I/O parameters
  const string GZIP_SUFFIX = ".gz";
  const string BZIP2_SUFFIX = ".bz2";
  const string BAM_SUFFIX = ".bam"; 
  const string STDIN_STRING = "stdin"; // like UCSC
  const string STDOUT_STRING = "stdout"; // like UCSC

  // Alignment-parsing parameters
  const UShort DEFAULT_MISMATCH_TOLERANCE = 1;
  const double PROB_BASE_ERROR = 0.01;
  const double DEFAULT_PROB_THRESHOLD = 0.9;

  // UCSC output parameters
  const string FORWARD_COLOR = "0,0,255";
  const string REVERSE_COLOR = "255,0,0";
  const string NONDIRECTIONAL_COLOR = "191,0,191";
  const UShort ALIGN_PRIORITY = 3; // since there will be several alignment files per profile, it's useful to have them underneath rather than keep scrolling down
  const UShort PROFILE_PRIORITY = 2; // regions go on top

  // Peak-calling parameters
  const UShort DEFAULT_BANDWIDTH = 50;
  const Score DEFAULT_REGION_THRESHOLD = 25;
  const double DEFAULT_KURTOSIS_THRESHOLD = 50;
  const double DEFAULT_CORR_THRESHOLD = 0.3;
  const HitCount DEFAULT_HIT_THRESHOLD = 10; // minimum *average* hits per sample (threshold is this times number of experiments)

  // Strand-shift parameters
  const UShort DEFAULT_MIN_SHIFT = 25;
  const UShort DEFAULT_MAX_SHIFT = 150;
  const Counter DEFAULT_TEST_REGIONS = 1000;
  const UShort MAXIMA_SMOOTHENING_BANDWIDTH = 5;

}

#endif
