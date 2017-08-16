#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <deque>
#include <cassert>

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/format.hpp>

#include <tclap/CmdLine.h>

#include "defaults.hpp"
#include "data.hpp"
#include "filterstream.hpp"
#include "format.hpp"
#include "peakcall.hpp"


using namespace std;
using namespace boost;
using namespace TCLAP;
using namespace unipeak;

int main (int argc, char* argv[]) {

  // get arguments
  string regionsFname;
  string contigFname;
  string controlString;
  UShort mismatchTolerance;
  UShort useLength;
  string offsetString;
  double probThreshold;
  string assembly;
  UShort bandwidth;
  string profileFname;
  string trackName;
  Pos mappable;
  string coeffString;
  double regionThreshold;
  double kurtosisThreshold;
  double corrThreshold;
  double hitThreshold;
  bool outputPeaks;
  bool outputCorrs;
  vector<string> alignFnames;
  bool directional;
  bool quiet;
  try {
    CmdLine cmd("See README.TXT for more information", ' ', VERSION);
    SwitchArg quietArg("q", "quiet", "don't give contig-by-contig updates", cmd, false);
    SwitchArg directionalArg("D", "non-directional", "don't separate strands (data aren't from a directional protocol), default false", cmd, false);
    ValueArg<string> assemblyArg("a", "assembly", "reference assembly name, e.g. hg18", false, "", "name", cmd);
    ValueArg<string> trackNameArg("n", "name", "wiggle plot track name", false, "", "name", cmd);
    ValueArg<string> profileFnameArg("w", "wig", "output wiggle plot filename", false, "", "filename", cmd);
    ValueArg<Pos> mappableArg("m", "mappable", "number of mappable positions", false, 0, "positive integer", cmd);
    ValueArg<double> hitThresholdArg("t", "hitThreshold", "minimum average hits per experiment in a region, default 10", false, DEFAULT_HIT_THRESHOLD, "non-negative number", cmd);
    SwitchArg outputCorrsArg("y", "corr", "output strand correlations in addition to kurtosis boundaries (non-directional only), default false", cmd, false);
    ValueArg<double> corrThresholdArg("u", "corrThreshold", "minimum correlation between strands' densities (non-directional only), default 0.3", false, DEFAULT_CORR_THRESHOLD, "number", cmd);
    ValueArg<double> kurtosisThresholdArg("k", "kurtosisThreshold", "maximum kurtosis threshold for regions, default 50", false, DEFAULT_KURTOSIS_THRESHOLD, "non-negative number", cmd);
    ValueArg<double> regionThresholdArg("r", "regionThreshold", "minimum fold enrichment threshold for regions, default 25", false, DEFAULT_REGION_THRESHOLD, "positive number", cmd);
    ValueArg<string> coeffArg("z", "coeff", "coefficients to scale samples, or \'p\' to scale inversely to confident reads", false, "", "comma-separated numbers", cmd);
    ValueArg<UShort> bandwidthArg("b", "bandwidth", "kernel bandwidth, default 50", false, DEFAULT_BANDWIDTH, "positive integer", cmd);
    ValueArg<UShort> mismatchToleranceArg("i", "mismatches", "mismatch tolerance, default 2 (Eland extended only)", false, DEFAULT_MISMATCH_TOLERANCE, "non-negative integer", cmd);
    ValueArg<UShort> useLengthArg("l", "length", "sequence length to use, default entire sequence (Eland extended only)", false, 0, "positive integer", cmd);
    ValueArg<string> offsetArg("s", "shift", "number of bases downstream to shift alignments (barcode length and/or strand shift), default 0", false, "", "comma-separated integers", cmd);
    ValueArg<double> probThresholdArg("p", "prob", "posterior probability threshold for alignments", false, DEFAULT_PROB_THRESHOLD, "probability", cmd);
    ValueArg<string> controlArg("e", "exclude", "samples to exclude from peak-calling (e.g. negative controls)", false, "", "comma-separated positive integers", cmd);
    SwitchArg outputPeaksArg("f", "peaks", "output peak positions in addition to region boundaries, default false", cmd, false);
    ValueArg<string> regionsFnameArg("o", "out", "output regions filename", true, "", "filename", cmd);
    ValueArg<string> contigFnameArg("c", "contig", "contig table filename", true, "", "filename", cmd);
    UnlabeledMultiArg<string> alignFnamesArg("", ".wig files of alignments", true, "alignment filenames", cmd);
    cmd.parse(argc, argv);
    regionsFname = regionsFnameArg.getValue();
    contigFname = contigFnameArg.getValue();
    outputPeaks = outputPeaksArg.getValue();
    controlString = controlArg.getValue();
    mismatchTolerance = mismatchToleranceArg.getValue();
    useLength = useLengthArg.getValue();
    offsetString = offsetArg.getValue();
    probThreshold = probThresholdArg.getValue();
    assembly = assemblyArg.getValue();
    bandwidth = bandwidthArg.getValue();
    trackName = trackNameArg.getValue();
    profileFname = profileFnameArg.getValue();
    mappable = mappableArg.getValue();
    coeffString = coeffArg.getValue();
    regionThreshold = regionThresholdArg.getValue();
    kurtosisThreshold = kurtosisThresholdArg.getValue();
    outputCorrs = outputCorrsArg.getValue();
    corrThreshold = corrThresholdArg.getValue();
    hitThreshold = hitThresholdArg.getValue();
    alignFnames = alignFnamesArg.getValue();
    directional = ! directionalArg.getValue();
    quiet = quietArg.getValue();
  } catch (ArgException &e) {
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl << endl;
    exit(1);
  }
  if (! profileFname.empty() && trackName.empty()) trackName = getFnamePrefix(profileFname);
  if (directional) {
    if (corrThreshold != DEFAULT_CORR_THRESHOLD) cerr << "warning: correlation threshold is not used on strand-specific analysis" << endl;
    corrThreshold = -1;
    if (outputCorrs) cerr << "warning: strand correlations are not calculated for strand-specific analysis" << endl;
    outputCorrs = false;
  }
  
  // parse offset(s)
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
  
  // get control(s)
  vector<bool> control(alignFnames.size(), false);
  if (! controlString.empty()) {
    tokenizer<char_separator<char> > controlTokenizer(controlString, (char_separator<char>)(","));
    for (tokenizer<char_separator<char> >::const_iterator i = controlTokenizer.begin(); i != controlTokenizer.end(); ++i) {
      UShort index;
      try {
        index = lexical_cast<short>(((string)(*i)).c_str());
      } catch(bad_lexical_cast &) {
        cerr << "error: bad control index\n" << endl;
        exit(1);
      }
      if (index == 0) {
        cerr << "error: bad control index (first sample is 1)\n" << endl;
        exit(1);
      } else if (index > control.size()) {
        cerr << "error: bad control index (greater than number of samples)\n" << endl;
        exit(1);
      }      
      control[index - 1] = true;
    }
  }
  
  // scale hit threshold to number of non-control samples
  UShort nControl = 0;
  for (vector<bool>::const_iterator i = control.begin(); i != control.end(); ++i) if (*i) ++nControl; // this is necessary in case anyone enters the same index twice
  hitThreshold *= alignFnames.size() - nControl;
  
  // parse coefficients
  vector<double> coeffs;
  bool getProportionalCoeffs = false; // compute them later when we have read counts
  if (! coeffString.empty()) {
    if (coeffString == "p") {
      getProportionalCoeffs = true;
    } else {
    
      // read input
      tokenizer<char_separator<char> > coeffTokenizer(coeffString, (char_separator<char>)(","));
      for (tokenizer<char_separator<char> >::const_iterator i = coeffTokenizer.begin(); i != coeffTokenizer.end(); ++i) {
        try {
          coeffs.push_back(lexical_cast<double>(((string)(*i)).c_str()));
        } catch(bad_lexical_cast &) {
          cerr << "error: bad coeff argument\n" << endl;
          exit(1);
        }
      }
      if (coeffs.size() != alignFnames.size() - nControl) {
        cerr << "error: wrong number of coeff arguments\nmust have same number as non-control alignment files\n" << endl;
        exit(1);
      }
    }
  }
  
  const ContigTable contigs = parseContigTable(contigFname);
  
  // create input streams
  MultiParse parse(alignFnames.size());
  MultiParse::iterator parseIter;
  vector<string>::const_iterator fnameIter;
  vector<short>::const_iterator offsetIter;
  vector<bool>::const_iterator controlIter;
  Counter nNonControlTags = 0;
  Counter nTotalTags = 0;
  cerr << "reading alignment files..." << endl;
  for (parseIter = parse.begin(), fnameIter = alignFnames.begin(), offsetIter = offsets.begin(), controlIter = control.begin(); parseIter != parse.end() && fnameIter != alignFnames.end() && controlIter != control.end(); ++parseIter, ++fnameIter, ++controlIter) {
    const short offset = (offsets.empty() ? 0 : *offsetIter);
    if (directional) *parseIter = new ParseAlignStream(*fnameIter, &contigs, mismatchTolerance, useLength, offset, probThreshold); else *parseIter = new NondirParseAlignStream(*fnameIter, &contigs, mismatchTolerance, useLength, offset, probThreshold); // polymorphism!
    const Counter tags = (*parseIter)->getExpectedTags();
    if (! *controlIter) nNonControlTags += tags;
    nTotalTags += tags;
    (*parseIter)->readAlign(); // get first count; may take a while if nondirectional
    cerr << "  " << (*parseIter)->exptName() << ": " << tags << " tags" << endl;
    if (offsets.size() > 1) ++offsetIter;
  }
  
  // calculate background
  if (mappable == 0) mappable = contigs.getGenomeSize();
  cerr << nTotalTags << " usable tags";
  if (nNonControlTags != nTotalTags) cerr << ", " << nNonControlTags << " not from negative controls,";
  cerr << " at " << mappable << " mappable positions" << endl;
  double background = (double)nNonControlTags / (double)mappable;
  if (directional) background /= 2;
  cerr << "using background = " << background << " tags/position";
  if (directional) cerr << " on each strand";
  cerr << endl;
  
  // calculate scaling coefficients, scaled so that the modified total non-control read count is the same
  if (getProportionalCoeffs) {
    MultiParse::const_iterator i;
    vector<bool>::const_iterator j;
    assert(parse.size() == control.size());
    for (i = parse.begin(), j = control.begin(); i != parse.end() && j != control.end(); ++i, ++j) if (! *j) coeffs.push_back((double)nNonControlTags / ((double)((*i)->getExpectedTags()) * (double)(parse.size() - nControl)));
  } else if (! coeffs.empty()) {
    double scaledSum = 0;
    assert(coeffs.size() == parse.size() - nControl);
    vector<double>::const_iterator i;
    MultiParse::const_iterator j;
    vector<bool>::const_iterator k;
    for (i = coeffs.begin(), j = parse.begin(), k = control.begin(); i != coeffs.end() && j != parse.end() && k != control.end(); ++i, ++j, ++k) if (! *k) scaledSum += *i * (double)((*j)->getExpectedTags());
    for (vector<double>::iterator i = coeffs.begin(); i != coeffs.end(); ++i) *i *= nNonControlTags / scaledSum;
  }
  
  // generate kernel
  const Kernel kernel(bandwidth, 1 / background);
  
  // generate header
  stringstream commonHeader;
  for (vector<string>::const_iterator i = alignFnames.begin(); i != alignFnames.end(); ++i) {
    commonHeader << "# align_file=" << *i << "\n";
  }
  if (! offsets.empty()) {
    if (offsets.size() == 1) {
      commonHeader << "# shift=" << offsets[0] << "\n";
    } else {
      commonHeader << "# shifts=";
      for (UShort i = 0; i < offsets.size() - 1; ++i) commonHeader << offsets[i] << ",";
      commonHeader << offsets[offsets.size() - 1] << "\n";
    }
  }    
  commonHeader << "# contig_table=" << contigFname << "\n";
  commonHeader << "# bandwidth=" << bandwidth << "\n";

  if (parse.size() > 1) {
    commonHeader << "# tags=" << parse[0]->getExpectedTags();
    for (UShort i = 1; i < parse.size(); ++i) commonHeader << "," << parse[i]->getExpectedTags();
    commonHeader << "\n";
  }
  
  if (nControl > 0) {
    commonHeader << "# control=";
    bool first = true; // now this is hacky; just look up the join function
    for (UShort i = 0; i < control.size(); ++i) if (control[i]) {
      if (first) first = false; else commonHeader << ",";
      commonHeader << i + 1;
    }
    commonHeader << "\n";
  }
  
  if (! coeffs.empty()) {
    commonHeader << "# coeffs=" << coeffs[0];
    for (UShort i = 1; i < coeffs.size(); ++i) commonHeader << "," << coeffs[i];
    commonHeader << "\n";
  }
  
  commonHeader << "# total_tags=" << nNonControlTags << "\n";
  commonHeader << "# background=" << background << "\n";

  FormatOutStream* profileOut = 0; // ProfileBuffer knows what to do with null pointer
  if (! profileFname.empty()) {
    if (! assembly.empty()) {
      profileOut = new FormatOutStream(&contigs, profileFname, trackName, (directional ? DIRECTIONAL_WIG_FORMAT : NONDIRECTIONAL_WIG_FORMAT), PROFILE_PRIORITY, assembly);
    } else {
      profileOut = new FormatOutStream(&contigs, profileFname, trackName, (directional ? DIRECTIONAL_WIG_FORMAT : NONDIRECTIONAL_WIG_FORMAT), PROFILE_PRIORITY);
    }
  }
  if (profileOut) *profileOut << commonHeader.str();
  
  FormatOutStream regionsOut(&contigs, regionsFname, REGIONS_FORMAT);
  regionsOut.outputPeaks(outputPeaks);
  regionsOut.outputCorrs(outputCorrs);
  regionsOut << commonHeader.str();
  
  regionsOut << "# region_threshold=" << regionThreshold << "\n" <<
    "# kurtosis_threshold=" << kurtosisThreshold << "\n" <<
    "# corr_threshold=" << corrThreshold << "\n" <<
    "# hit_threshold=" << hitThreshold << "\n";
  
  // generate column names
//  regionsOut << "\tsize\tpositions sampled\tproportion\tSD"; // test
  if (outputPeaks) regionsOut << "\tpeak";
  if (outputCorrs) regionsOut << "\tcorrelation";
  regionsOut << "\tkurtosis"; // peak needs to be optional
  for (MultiParse::const_iterator i = parse.begin(); i != parse.end(); ++i) regionsOut << "\t" << (*i)->exptName();
  regionsOut << "\n";
  
  // read inputs
  cerr << "calling enriched regions..." << endl;
  vector<Region*> regions;
  ProfileBuffer forwardBuffer(&kernel, regionThreshold, kurtosisThreshold, corrThreshold, hitThreshold, true, control, coeffs, &contigs, &regions, profileOut);
  ProfileBuffer reverseBuffer(&kernel, regionThreshold, kurtosisThreshold, corrThreshold, hitThreshold, false, control, coeffs, &contigs, &regions, profileOut); // only used if directional, otherwise everything goes in forwardBuffer
  ContigNo contig = 0;
  bool forward = true; // only used if sorted by strand then contig
  while (contig < contigs.size()) { // each iteration is a new contig+strand but order isn't quite fixed
    if (! quiet) {
      cerr << "  " << contigs[contig];
      if (directional) {
        if (forward) {
          cerr << "+";
        } else {
          cerr << "-";
        }
      }
      cerr << "... " << flush;
    }
    
    const Pos posLimit = contigs.getSize(contig);
    Pos pos = 1;
    while (pos <= posLimit) { // not necessarily sequential; follows input to skip empty positions
      bool foundForwardHit = false;
      bool foundReverseHit = false;
      Pos nextPos = posLimit + 1;
      HitCountVec forwardHits(parse.size());
      HitCountVec reverseHits(parse.size());
      MultiParse::iterator parseIter;
      HitCountVec::iterator forwardHitIter;
      HitCountVec::iterator reverseHitIter;
      for (parseIter = parse.begin(), forwardHitIter = forwardHits.begin(), reverseHitIter = reverseHits.begin(); parseIter != parse.end() && (forwardHitIter != forwardHits.end() || reverseHitIter != reverseHits.end()); ++parseIter, ++forwardHitIter, ++reverseHitIter) {
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
      }
       
      if (foundForwardHit) forwardBuffer.add(forwardHits, contig, pos, true);
      if (foundReverseHit) {
        if (directional) reverseBuffer.add(reverseHits, contig, pos, false); else forwardBuffer.add(reverseHits, contig, pos, false);
      }
      
      // flush regions vector - not necessary to have them all in memory at once unless I think of a reason later
      if (regions.size() != 0) {
        for (vector<Region*>::const_iterator i = regions.begin(); i != regions.end(); ++i) {
          regionsOut.write(**i);
          delete *i;
        }
        regions.clear();
      }
      
      pos = nextPos; // if all inputs are done with this contig, it remains above posLimit and loop ends
    }

    const Counter forwardRegions = forwardBuffer.flushContig();
    const Counter reverseRegions = reverseBuffer.flushContig();
    
    // print summary
    if (! quiet) {
      if (forward) {
        cerr << forwardRegions << endl;
        if (directional && reverseRegions > 0) cerr << "  " << contigs[contig] << "-... " << reverseRegions << endl;
      } else {
        cerr << reverseRegions << endl;
      }
    }

    ++contig;    
    if (directional && contig == contigs.size() && forward) {
      contig = 0;
      forward = false;
    }
  }
  if (profileOut) delete profileOut;
  regionsOut.close();
  
  // summary
  cerr << forwardBuffer.nRegions() + reverseBuffer.nRegions() << " regions passed filters, " << forwardBuffer.nRegionRejects() + reverseBuffer.nRegionRejects() << " rejected" << endl;

  cerr << "tags in regions:" << endl;
  for (UShort i = 0; i < parse.size(); ++i) {
    ParseAlignStream& thisParse = *parse[i];
    if (thisParse.confidentTagCount() + thisParse.outOfBoundsTagCount() != thisParse.getExpectedTags()) cerr << "  warning: expected " << thisParse.getExpectedTags() << " tags in " << thisParse.exptName() << " but found " << thisParse.confidentTagCount() << "; results inaccurate" << endl;
    const Counter totalTagsInRegions = forwardBuffer.nTagsInRegions()[i] + reverseBuffer.nTagsInRegions()[i];
    cerr << "  " << thisParse.exptName() << ": " << totalTagsInRegions << " (" << format("%.1f") % (100 * (double)totalTagsInRegions / (double)thisParse.confidentTagCount()) << "%)" << endl;
    delete parse[i];
  }
 
  cerr << "\nAll done!\n" << endl;
  return 0;
}

