#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/format.hpp>

#include <tclap/CmdLine.h>

#include "defaults.hpp"
#include "data.hpp"
#include "filterstream.hpp"
#include "format.hpp"

using namespace std;
using namespace boost;
using namespace TCLAP;
using namespace unipeak;

int main (int argc, char* argv[]) {

  // get arguments
  string regionsInFname;
  string regionsOutFname;
  string contigFname;
  Pos extension;
  UShort mismatchTolerance;
  UShort useLength;
  string offsetString;
  double probThreshold;
  vector<string> alignFnames;
  bool directional;
  try {
    CmdLine cmd("See README.TXT for more information", ' ', VERSION);
    SwitchArg directionalArg("D", "non-directional", "don't separate strands (data aren't from a directional protocol), default false", cmd, false);
    ValueArg<UShort> mismatchToleranceArg("i", "mismatches", "mismatch tolerance, default 2 (Eland extended only)", false, DEFAULT_MISMATCH_TOLERANCE, "non-negative integer", cmd);
    ValueArg<UShort> useLengthArg("l", "length", "sequence length to use, default entire sequence (Eland extended only)", false, 0, "positive integer", cmd);
    ValueArg<string> offsetArg("s", "shift", "number of bases downstream to shift alignments (barcode length and/or strand shift), default 0", false, "", "comma-separated integers", cmd);
    ValueArg<double> probThresholdArg("p", "prob", "posterior probability threshold for alignments", false, DEFAULT_PROB_THRESHOLD, "probability", cmd);
    ValueArg<Pos> extensionArg("e", "extend", "extend region boundaries by a fixed amount", false, 0, "positive integer", cmd);
    ValueArg<string> regionsOutFnameArg("o", "out", "output regions filename", true, "", "filename", cmd);
    ValueArg<string> regionsInFnameArg("f", "in", "input regions filename", true, "", "filename", cmd);
    ValueArg<string> contigFnameArg("c", "contig", "contig table filename", true, "", "filename", cmd);
    UnlabeledMultiArg<string> alignFnamesArg("", ".wig files of alignments", true, "alignment filenames", cmd);
    cmd.parse(argc, argv);
    regionsInFname = regionsInFnameArg.getValue();
    regionsOutFname = regionsOutFnameArg.getValue();
    contigFname = contigFnameArg.getValue();
    extension = extensionArg.getValue();
    mismatchTolerance = mismatchToleranceArg.getValue();
    useLength = useLengthArg.getValue();
    offsetString = offsetArg.getValue();
    probThreshold = probThresholdArg.getValue();
    alignFnames = alignFnamesArg.getValue();
    directional = ! directionalArg.getValue();
  } catch (ArgException &e) {
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl << endl;
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
    if (directional) *parseIter = new ParseAlignStream(*fnameIter, &contigs, mismatchTolerance, useLength, offset, probThreshold); else *parseIter = new NondirParseAlignStream(*fnameIter, &contigs, mismatchTolerance, useLength, offset, probThreshold); // polymorphism!
    const Counter tags = (*parseIter)->getExpectedTags();
    nTotalTags += tags;
    (*parseIter)->readAlign(); // get first count; may take a while if nondirectional
    cerr << "  " << (*parseIter)->exptName() << ": " << tags << " tags" << endl;
    ++parseIter;
    ++fnameIter;
    if (offsets.size() > 1) ++offsetIter;
  }
  
  // generate header
  InStream regionsIn(regionsInFname);
  OutStream regionsOut(regionsOutFname);
  string line = regionsIn.readLine();
  while (line.empty() || line[0] == '#') { // copy previous header
    regionsOut << line << "\n";
    line = regionsIn.readLine();
  }
  if (line[0] != '\t') { // first line after header must be column names
    cerr << "error: bad format in " << regionsInFname << " line " << regionsIn.getLineNo() << "\n\n";
    exit(1);
  }
  if (extension != 0) regionsOut << "# region_extension=" << extension << "\n";
  for (vector<string>::const_iterator fnameIter = alignFnames.begin(); fnameIter != alignFnames.end(); ++fnameIter) regionsOut << "# extra_align_file=" << *fnameIter << "\n";
  if (! offsets.empty()) {
    if (offsets.size() == 1) {
      regionsOut << "# shift=" << offsets[0] << "\n";
    } else {
      regionsOut << "# shifts=";
      for (UShort i = 0; i < offsets.size() - 1; ++i) regionsOut << offsets[i] << ",";
      regionsOut << offsets[offsets.size() - 1] << "\n";
    }
  }
  
  // generate column names
  regionsOut << line;
  for (MultiParse::const_iterator parseIter = parse.begin(); parseIter != parse.end(); ++parseIter) regionsOut << "\t" << (*parseIter)->exptName();
  regionsOut << "\n";
  
  cerr << "processing regions... " << flush;
  vector<Counter> tagsInRegions(parse.size(), 0);
  Counter nRegions = 0;
  ContigNo lastContig = contigs.size();
  Pos lastRight = 0;
  while (regionsIn.good()) {
    const string line = regionsIn.readLine();
    if (line.empty()) continue;
    
    // parse region boundaries
    const size_t colon = line.find_first_of(":");
    const size_t dash = line.find_first_of("-");
    const size_t firstTab = line.find_first_of("\t");
    if (firstTab == string::npos || colon == string::npos || dash == string::npos || dash > firstTab || colon > dash) {
      cerr << "error: bad format in " << regionsInFname << " line " << regionsIn.getLineNo() << "\n\n";
      exit(1);
    }
    const ContigNo contig = contigs[line.substr(0, colon)];
    if (contig == contigs.size()) {
      cerr << "error: contig not in table in " << regionsInFname << " line " << regionsIn.getLineNo() << "\n\n";
      exit(1);
    }
    Pos start;
    Pos end;
    try {
      start = lexical_cast<Pos>(line.substr(colon + 1, dash - colon - 1));
      end = lexical_cast<Pos>(line.substr(dash + 1, firstTab - dash - 1));
    } catch(bad_lexical_cast &) {
      cerr << "error: bad format in " << regionsInFname << " line " << regionsIn.getLineNo() << "\n\n";
      exit(1);
    }
//    if (start > contigs.getSize(contig) || end > contigs.getSize(contig)) {
//      cerr << "error: region out of bounds in " << regionsInFname << " line " << regionsIn.getLineNo() << "\n\n";
//      exit(1);
//    }
    const bool forward = end >= start;
    if (! forward && ! directional) {
      cerr << "error: reverse regions in non-directional analysis\n\n";
      exit(1);
    }
    const Pos left = ((forward ? start : end) < extension ? 0 : (forward ? start : end) - extension); // don't let it go past the left end of the contig
    const Pos right = (forward ? end : start) + extension;
    if (contig == lastContig && left <= lastRight) {
      cerr << "error: " << line.substr(0, firstTab) << " overlaps previous region";
      if (extension != 0) cerr << "; decrease extension";
      cerr << "\n\n";
      exit(1);
    }
    regionsOut << line;
  
    // get counts from new align files, one at a time
    vector<Counter>::iterator tagsInRegionsIter = tagsInRegions.begin();
    for (MultiParse::iterator parseIter = parse.begin(); parseIter != parse.end(); ++parseIter) {
      const Alignment* align = &(*parseIter)->lastAlign(); // this is a pointer, not a reference, just in case the address changes (though it looks like it probably doesn't)
      while (align->count != 0 && (forward != align->forward || align->contig < contig || (align->contig == contig && align->firstPos < left))) align = &(*parseIter)->readAlign(); // skip all reads before region
      HitCount hits = 0;
      while (align->contig == contig && align->firstPos <= right) {
        hits += align->count;
        align = &(*parseIter)->readAlign();
      }
      
      regionsOut << "\t" << hits;
      *tagsInRegionsIter += hits;
      ++tagsInRegionsIter;
    }
    
    regionsOut << "\n";
    ++nRegions;
  }
  regionsIn.close();
  regionsOut.close();
  
  // summary
  cerr << nRegions << " in " << regionsInFname << endl;
  cerr << "tags in regions:" << endl;
  for (UShort i = 0; i < parse.size(); ++i) {
    ParseAlignStream& thisParse = *parse[i];
    cerr << "  " << thisParse.exptName() << ": " << tagsInRegions[i] << " (" << format("%.1f") % (100 * (double)tagsInRegions[i] / (double)thisParse.getExpectedTags()) << "%)" << endl;
    delete parse[i];
  } 
  
  cerr << "\nDone!\n" << endl;
  return(0);
}  
