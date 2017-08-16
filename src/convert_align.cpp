#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <ctime>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
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

int main (int argc, char* argv[]) {

  // get arguments
  string trackName;
  string outFname;
  string contigFname;
  UShort mismatchTolerance;
  UShort useLength;
  short offset;
  double probThreshold;
  string assembly;
  vector<string> alignFnames;
  bool directional;
  bool quiet;
  try {
    CmdLine cmd("See README.TXT for more information", ' ', VERSION);
    SwitchArg quietArg("q", "quiet", "don't give contig-by-contig updates", cmd, false);
    SwitchArg directionalArg("D", "non-directional", "combine strands (data aren't from a directional protocol), default false", cmd, false);
    ValueArg<UShort> mismatchToleranceArg("i", "mismatches", "mismatch tolerance, default 2", false, DEFAULT_MISMATCH_TOLERANCE, "non-negative integer", cmd);
    ValueArg<UShort> useLengthArg("l", "length", "sequence length to use, default entire sequence", false, 0, "positive integer", cmd);
    ValueArg<short> offsetArg("s", "shift", "number of bases downstream to shift alignments (barcode length + strand shift), default 0", false, 0, "integer", cmd);
    ValueArg<double> probThresholdArg("p", "prob", "posterior probability threshold for alignments, default 0.9", false, DEFAULT_PROB_THRESHOLD, "probability", cmd);
    ValueArg<string> assemblyArg("a", "assembly", "reference assembly name, e.g. hg18", false, "", "name", cmd);
    ValueArg<string> trackNameArg("n", "name", "track name", false, "", "name", cmd);
    ValueArg<string> contigFnameArg("c", "contigs", "contig table filename", true, "", "filename", cmd);
    ValueArg<string> outFnameArg("o", "output", "output file", true, "", "output filename", cmd);
    UnlabeledMultiArg<string> alignFnamesArg("input", "alignment files", true, "align filenames", cmd);
    cmd.parse(argc, argv);
    trackName = trackNameArg.getValue();
    outFname = outFnameArg.getValue();
    contigFname = contigFnameArg.getValue();
    mismatchTolerance = mismatchToleranceArg.getValue();
    useLength = useLengthArg.getValue();
    offset = offsetArg.getValue();
    probThreshold = probThresholdArg.getValue();
    assembly = assemblyArg.getValue();
    alignFnames = alignFnamesArg.getValue();
    directional = ! directionalArg.getValue();
    quiet = quietArg.getValue();
  } catch (ArgException& e) {
    cerr << "error: " << e.error() << " for arg " << e.argId() << "\n" << endl;
    exit(1);
  }
  
  const ContigTable contigs = parseContigTable(contigFname);

  // read alignments
  CountMap exptCounts(contigs.size());
  ParseAlignStream in(&contigs, mismatchTolerance, useLength, offset, probThreshold);
  for (vector<string>::const_iterator i = alignFnames.begin(); i != alignFnames.end(); ++i) {
    in.open(*i);
//    out.write(string("# original_file=" + *i + "\n")); // header
    cerr << "reading " << *i << "..." << endl;
    if (! quiet) cerr << "0 tags read";
    time_t lastUpdate = time(0);
    while (in.good()) {
//      Alignment align = in.readAlign();
//      if (! directional) align.forward = true; // firstPos will be less than lastPos but it shouldn't affect anything
//      exptCounts.add(align);
      exptCounts.add(in.readAlign());
      
      // status updates
      if (! quiet) {
        if (difftime(time(0), lastUpdate) > 1) { // every second
          cerr << "\r" << in.totalTagCount() << " tags read";
          lastUpdate = time(0);
        }
      }
    }
    in.close();
    if (! quiet) cerr << "\r";
    in.printSummary();
  }
  cerr << exptCounts.tagCount() << " usable tags" << endl;
  if (exptCounts.tagCount() == 0) {
    cerr << "error: nothing to do" << endl << endl;;
    exit(1);
  }
  
  // determine output mode
  if (trackName.empty()) trackName = getFnamePrefix(outFname);
  cerr << "writing " << outFname << " in ";
  UShort format;
  //  if (regex_search(outFname.c_str(), wig_suffix_regex)) { // restore this to allow BED again
    format = (directional ? DIRECTIONAL_WIG_FORMAT : NONDIRECTIONAL_WIG_FORMAT);
    cerr << "wiggle";
  //  } else {
  //    format = BED_FORMAT;
  //    cerr << "BED";
  //  }
  cerr << " format with " << (directional ? "separate strands" : "strands combined") << "... " << flush;

  // open output
  FormatOutStream out(&contigs);
  if (assembly.empty()) {
    out.open(outFname, trackName, format, ALIGN_PRIORITY);
  } else {
    out.open(outFname, trackName, format, ALIGN_PRIORITY, assembly);
  }
  
  // write header
  for (vector<string>::const_iterator i = alignFnames.begin(); i != alignFnames.end(); ++i) out.write(string("# original_file=" + *i + "\n"));  
  if (probThreshold != 0) out.write("# prob_threshold=" + lexical_cast<string>(probThreshold) + "\n");
  out.write("# tags=" + lexical_cast<string>(exptCounts.tagCount()) + "\n");

  // write output  
  if (directional) {
    for (CountMap::ConstIterator i = exptCounts.begin(); i != exptCounts.end(); ++i) out.write(*i);
  } else {
    for (CountMap::ConstNondirIterator i = exptCounts.nondirBegin(); i != exptCounts.nondirEnd(); ++i) out.write(*i);
  }  
  out.close();
  
  // test
  if (out.tagCount() != exptCounts.tagCount()) { // test
    cerr << "error: " << out.tagCount() << " " << exptCounts.tagCount() << endl;
    exit(1);
  }
  
  cerr << "done!\n" << endl;
//  cin.get(); // test: pause so I can check memory usage
  return 0;
}
