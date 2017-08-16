#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/format.hpp>
#include <tclap/CmdLine.h>

#include "defaults.hpp"
#include "filterstream.hpp"
#include "format.hpp"
#include "data.hpp"
#include "peakcall.hpp"

using namespace std;
using namespace boost;
using namespace TCLAP;
using namespace unipeak;


bool compareRegionCounts (const Region* const i, const Region* const j) {
  assert(i && j);
  return(i->sum() > j->sum());
}

int main (int argc, char* argv[]) {

  // get arguments
  string trackName;
  string contigFname;
  string outFname;
  UShort mismatchTolerance;
  UShort useLength;
  string offsetString;
  double probThreshold;
  UShort bandwidth;
  Pos mappable;
  double regionThreshold;
  double kurtosisThreshold;
  double corrThreshold;
  HitCount hitThreshold;
  UShort minShift;
  UShort maxShift;
  UShort nRegionsToTest;
  vector<string> alignFnames;
  try {
    CmdLine cmd("See README.TXT for more information", ' ', VERSION);
    ValueArg<Pos> mappableArg("m", "mappable", "number of mappable positions, default sum of contigs", false, 0, "positive integer", cmd);
    ValueArg<HitCount> hitThresholdArg("t", "hitThreshold", "minimum average hits per experiment in a region, default 10", false, DEFAULT_HIT_THRESHOLD, "non-negative integer", cmd);
    ValueArg<double> corrThresholdArg("u", "corrThreshold", "minimum correlation between strands' densities (non-directional only), default 0.3", false, DEFAULT_CORR_THRESHOLD, "number", cmd);
    ValueArg<double> kurtosisThresholdArg("k", "kurtosisThreshold", "maximum kurtosis threshold for regions, default 50", false, DEFAULT_KURTOSIS_THRESHOLD, "positive number", cmd);
    ValueArg<double> regionThresholdArg("r", "regionThreshold", "minimum fold enrichment threshold for regions, default 25", false, DEFAULT_REGION_THRESHOLD, "positive number", cmd);
    ValueArg<UShort> bandwidthArg("b", "bandwidth", "kernel bandwidth, default 50", false, DEFAULT_BANDWIDTH, "positive integer", cmd);
    ValueArg<UShort> nRegionsToTestArg("g", "testRegions", "number of regions to use for calibration, default 1000", false, DEFAULT_TEST_REGIONS, "non-negative integer", cmd);
    ValueArg<UShort> maxShiftArg("x", "maxShift", "maximum shift to try, default 150", false, DEFAULT_MAX_SHIFT, "positive integer", cmd);
    ValueArg<UShort> minShiftArg("n", "minShift", "minimum shift to try, default 25", false, DEFAULT_MIN_SHIFT, "positive integer", cmd);
    ValueArg<UShort> mismatchToleranceArg("i", "mismatches", "mismatch tolerance, default 2", false, DEFAULT_MISMATCH_TOLERANCE, "non-negative integer", cmd);
    ValueArg<UShort> useLengthArg("l", "length", "sequence length to use, default entire sequence", false, 0, "positive integer", cmd);
    ValueArg<string> offsetArg("s", "shift", "number of bases downstream to shift alignments (barcode length and/or strand shift), default 0", false, "", "comma-separated integers", cmd);
    ValueArg<double> probThresholdArg("p", "prob", "posterior probability threshold for alignments, default 0.9", false, DEFAULT_PROB_THRESHOLD, "probability", cmd);
    ValueArg<string> outFnameArg("o", "out", "output filename", false, "", "filename", cmd);
    ValueArg<string> contigFnameArg("c", "contigs", "contig table filename", true, "", "filename", cmd);
    UnlabeledMultiArg<string> alignFnamesArg("", ".wig files of alignments", true, "alignment filenames", cmd);
    cmd.parse(argc, argv);
    contigFname = contigFnameArg.getValue();
    outFname = outFnameArg.getValue();
    mismatchTolerance = mismatchToleranceArg.getValue();
    useLength = useLengthArg.getValue();
    offsetString = offsetArg.getValue();
    probThreshold = probThresholdArg.getValue();
    bandwidth = bandwidthArg.getValue();
    mappable = mappableArg.getValue();
    regionThreshold = regionThresholdArg.getValue();
    kurtosisThreshold = kurtosisThresholdArg.getValue();
    corrThreshold = corrThresholdArg.getValue();
    hitThreshold = hitThresholdArg.getValue();
    nRegionsToTest = nRegionsToTestArg.getValue();
    minShift = minShiftArg.getValue();
    maxShift = maxShiftArg.getValue();
    alignFnames = alignFnamesArg.getValue();
  } catch (ArgException& e) {
    cerr << "error: " << e.error() << " for arg " << e.argId() << "\n" << endl;
    exit(1);
  }

  // get offset(s)
  vector<short> offsets;
  if (! offsetString.empty()) {
    tokenizer<char_separator<char> > offsetTokenizer(offsetString, (char_separator<char>)(","));
    for (tokenizer<char_separator<char> >::const_iterator i = offsetTokenizer.begin(); i != offsetTokenizer.end(); ++i) {
      try {
        offsets.push_back(lexical_cast<short>(((string)(*i)).c_str()));
      } catch(bad_lexical_cast &) {
        cerr << "error: bad offset argument\n" << endl;
        exit(1);
      }
    }
    if (! (offsets.size() == alignFnames.size() || offsets.size() == 1)) {
      cerr << "error: wrong number of offset arguments\nmust have same number as alignment files or just one\n" << endl;
      exit(1);
    }
  }
  
  const ContigTable contigs = parseContigTable(contigFname);

  // create input streams
  MultiParse parse(alignFnames.size());
  MultiParse::iterator parseIter = parse.begin();
  vector<string>::const_iterator fnameIter = alignFnames.begin();
  vector<short>::const_iterator offsetIter = offsets.begin();
  Counter nTotalTags = 0;
  cerr << "reading alignment files..." << endl;
  while (parseIter != parse.end() && fnameIter != alignFnames.end()) {
    const short offset = (offsets.empty() ? 0 : *offsetIter);
    *parseIter = new NondirParseAlignStream(*fnameIter, &contigs, mismatchTolerance, useLength, offset, probThreshold); // polymorphism!
    const Counter tags = (*parseIter)->getExpectedTags();
    nTotalTags += tags;
    (*parseIter)->readAlign(); // get first count; may take a while if nondirectional
    cerr << "  " << (*parseIter)->exptName() << ": " << tags << " tags" << endl;
    ++parseIter;
    ++fnameIter;
    if (offsets.size() > 1) ++offsetIter;
  }

  // calculate background
  if (mappable == 0) mappable = contigs.getGenomeSize();
  cerr << nTotalTags << " usable tags at " << mappable << " mappable positions" << endl;
  double background = (double)nTotalTags / (double)mappable;
  cerr << "using background = " << background << " tags/position" << endl;

  // generate kernel
  const Kernel kernel(bandwidth, 1 / background);

  // read inputs
  cerr << "calling enriched regions for shift calibration... " << flush;
  vector<Region*> regions;
  ProfileBuffer buffer(&kernel, regionThreshold, kurtosisThreshold, -1, hitThreshold, true, vector<bool>(parse.size(), false), vector<double>(), &contigs, &regions);
  ContigNo contig = 0;
  while (contig < contigs.size()) { // each iteration is a new contig+strand but order isn't quite fixed
    
    const Pos posLimit = contigs.getSize(contig);
    Pos pos = 1;
    while (pos <= posLimit) { // not necessarily sequential; follows input to skip empty positions
      bool foundForwardHit = false;
      bool foundReverseHit = false;
      Pos nextPos = posLimit + 1;
      MultiParse::iterator parseIter = parse.begin();
      HitCountVec forwardHits(parse.size());
      HitCountVec reverseHits(parse.size());
      HitCountVec::iterator forwardHitIter = forwardHits.begin();
      HitCountVec::iterator reverseHitIter = reverseHits.begin();
      while (parseIter != parse.end() && (forwardHitIter != forwardHits.end() || reverseHitIter != reverseHits.end())) {     
        *forwardHitIter = 0;
        *reverseHitIter = 0;
        
        // count hits if at same position, and loop in case only one tag per line
        const Alignment& align = (*parseIter)->lastAlign();
        while (align.count != 0 && align.firstPos == pos && align.contig == contig) { 
          if (align.forward) {
            *forwardHitIter += align.count;
            foundForwardHit = true;
          } else {
            *reverseHitIter += align.count;
            foundReverseHit = true;
          }
          (*parseIter)->readAlign();
        }    
        
        // get next position, if on same contig     
        if (align.count != 0 && align.contig == contig && align.firstPos < nextPos) nextPos = align.firstPos;
        
        ++parseIter;
        ++forwardHitIter;
        ++reverseHitIter;
      }
       
      if (foundForwardHit) buffer.add(forwardHits, contig, pos, true);
      if (foundReverseHit) buffer.add(reverseHits, contig, pos, false);
      
      pos = nextPos; // if all inputs are done with this contig, it remains above posLimit and loop ends
    }

    buffer.flushContig();
    ++contig;
  }
  cerr << regions.size() << " found\n";
  
  // garbage collection
  for (MultiParse::iterator i = parse.begin(); i != parse.end(); ++i) if (*i) delete *i; // when is *i null? dunno, but this is always a safe way to delete pointers
  
  // partially sort regions by count (in decreasing order)
//  partial_sort(regions.begin(), (regions.size() > nRegionsToTest ? regions.begin() + nRegionsToTest : regions.end()), regions.end(), compareRegionCounts);
  sort(regions.begin(), regions.end(), compareRegionCounts);
  
  // get shift scores
  UShort nRegionsTested = 0;
  Counter nTagsInRegions = 0;
  vector<Counter> shiftFreq(maxShift + 1, 0);
//  for (vector<Region*>::const_iterator i = regions.begin(); nRegionsTested < nRegionsToTest && i != regions.end() && i != regions.begin() + nRegionsToTest; ++i) {
  for (vector<Region*>::const_iterator i = regions.begin(); nRegionsTested < nRegionsToTest && i != regions.end(); ++i) {
    if ((*i)->counts->size() > (unsigned)(2 * maxShift + 3) // large enough to try all shifts   
//          ((! control) || (double)((*i)->exptSums()[0] + (*i)->exptSums()[1]) / (double)((*i)->exptSums()[2] + (*i)->exptSums()[3]) >= controlThreshold)) { // if using control, exceeds threshold for relative enrichment
    ) {
      UShort bestShift = 0;
      double bestCorr = -1;
      for (UShort shift = 0; shift <= maxShift; ++shift) {
        const double corr = (*i)->strandCorr(shift);
        if (corr > bestCorr) { // if corr is somehow identical, this inequality means the lowest shift will be used
          bestShift = shift;
          bestCorr = corr;
        }
        
      }
      
      if (bestCorr >= corrThreshold) { // disqualify regions under threshold
        ++shiftFreq[bestShift];
        nTagsInRegions += (*i)->sum();
        ++nRegionsTested;
      }
    }  
//    delete *i; // garbage collection
    // TO DO: if reached end of sorted regions and still need more, sort more
  }
  
  // check number of regions used
  if (nRegionsTested == 0) {
    cerr << "error: no regions qualified with given settings" << endl << endl;
    exit(1);
  } else if (nRegionsTested < nRegionsToTest) {
    cerr << "warning: too few regions qualified with given settings" << endl;
  }
  cerr << "the top " << nRegionsTested << " qualified regions contained " << nTagsInRegions << " tags (" << format("%.1f") % (100 * (double)nTagsInRegions / (double)nTotalTags) << "%)" << endl;
  
    
  // smoothen distribution of maxima
  const Kernel modeKernel(MAXIMA_SMOOTHENING_BANDWIDTH);
  vector<double> modeDensity(shiftFreq.size(), 0);
   for (UShort i = 0; i < shiftFreq.size(); ++i) { // these loops are hacky
     for (int j = 0; j < (int)modeKernel.size(); ++j) {
      const int k = (int)i - (int)MAXIMA_SMOOTHENING_BANDWIDTH + j;
      if (k >= 0 && k < (int)modeDensity.size()) modeDensity[k] += (double)shiftFreq[i] * modeKernel[j];
    }
  } 
   
  // get best shift
  double bestDensity = 0;
  UShort bestShift = 0;
  for (UShort i = minShift; i < modeDensity.size() - MAXIMA_SMOOTHENING_BANDWIDTH; ++i) { // don't let them pile up at the ends
     if (modeDensity[i] > bestDensity) {
      bestShift = i;
      bestDensity = modeDensity[i];
    }
  }
  cerr << "estimated shift = " << bestShift << endl;
  
  // write output
  if (! outFname.empty()) {
    cerr << "writing output to " << outFname << "... " << flush;
    OutStream out(outFname);
    for (vector<string>::const_iterator i = alignFnames.begin(); i != alignFnames.end(); ++i) {
      out << "# align_file=" << *i << "\n";
    }
    if (! offsets.empty()) {
      if (offsets.size() == 1) {
        out << "# shift=" << offsets[0] << "\n";
      } else {
        out << "# shifts=";
        for (UShort i = 0; i < offsets.size() - 1; ++i) out << offsets[i] << ",";
        out << offsets[offsets.size() - 1] << "\n";
      }
    }  
    out << "# contig_table=" << contigFname << "\n" <<
      "# bandwidth=" << bandwidth << "\n" <<
      "# tags=" << nTotalTags << "\n" <<
      "# background=" << background << "\n" <<
      "# region_threshold=" << regionThreshold << "\n" <<
      "# kurtosis_threshold=" << kurtosisThreshold << "\n" <<
      "# corr_threshold=" << corrThreshold << "\n" <<
      "# hit_threshold=" << hitThreshold << "\n" <<
      "# regions_tested=" << nRegionsTested << "\n" <<
      "# tags_in_regions=" << nTagsInRegions << "\n" <<
      "# min_shift=" << minShift << "\n" <<
      "# best_shift=" << bestShift << "\n\n";
    // report mode frequences
    out << "shift\tregions\n";
    for (UShort i = 0; i < shiftFreq.size(); ++i) {
      if (shiftFreq[i] != 0) out << i << "\t" << shiftFreq[i] << "\n";
    }
//    // report smoothened mode distribution
//    out << "\n";
//    for (UShort i = 0; i < modeDensity.size(); ++i) {
//      out << i << "\t" << modeDensity[i] << "\n";
//    }
    out.close();
  }
  
  cerr << "done!\n" << endl;
  return 0;
}
