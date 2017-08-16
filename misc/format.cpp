#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>

#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/cstdint.hpp>

#include "bamtools/BamReader.h"

#include "defaults.hpp"
#include "data.hpp"
#include "filterstream.hpp"
#include "format.hpp"

namespace unipeak {

  using namespace std;
  using namespace boost;
  using namespace BamTools;

  ContigTable parseContigTable (const string& filename) {
    cerr << "reading " << filename << "... " << flush;
    ContigTable result;
    InStream fileStream(filename);
    while (fileStream.good()) {
      const string line = fileStream.readLine();
      if (line.length() && line[0] != '#') {
        cmatch fields;
        if (regex_search(line.c_str(), fields, CONTIG_TABLE_REGEX)) {
          const string contig = fields[1];
          const Pos size = lexical_cast<Pos>(((string)fields[2]).c_str());
          if (result[contig] != result.size()) {
            cerr << "error: " << contig << " defined twice in contig table\n" << endl;
            exit(1);
          } else {
            result.add(contig, size);
          }
        }
      }
    }
    fileStream.close();
    
    if (result.size() == 0) {
      cerr << "error: no contigs in table" << endl << endl;
      exit(1);
    } else {
      cerr << result.size() << " contigs" << endl;
    }
    
    return result;
  }      


  string getFnamePrefix (const string& name) {
    string result = name;
    const size_t lastSlash = result.rfind("/");
    if (lastSlash != string::npos) {result = result.substr(lastSlash + 1);}
    const size_t firstDot = result.find(".");
    if (firstDot != string::npos) {result = result.substr(0, firstDot);}
    return result;
  }

  BinomPosterior::BinomPosterior (const UShort& readLength) : coefficients_(new vector<double>(readLength + 1)) {
    const double c = PROB_BASE_ERROR / (1 - PROB_BASE_ERROR); 
    (*coefficients_)[0] = 1;
    for (UShort i = 1; i <= readLength; ++i) (*coefficients_)[i] = (*coefficients_)[i - 1] * c * (readLength - i + 1) / i;
  }

  BinomPosterior::~BinomPosterior () {
    delete coefficients_;
  }

  double BinomPosterior::prob (const UShort& mismatches, const vector<HitCount>& hits) const {
    assert(mismatches < coefficients_->size() && hits.size() <= coefficients_->size());
    const double numerator = (*coefficients_)[mismatches];
    double denominator = 0;
    for (UShort i = 0; i < hits.size(); ++i) denominator += hits[i] * (*coefficients_)[i];
    assert(denominator > 0);
    return numerator / denominator;
  }
  
  void ParseAlignStream::close () {
    if (format_ == BAM_FORMAT) {
      if (bam_) { 
        delete bam_;
        bam_ = 0;
      }
      if (bamAlign_) {
        delete bamAlign_;
        bamAlign_ = 0;
      }
    } else {
      InStream::close();
    }
  }

  void ParseAlignStream::open (const string& fnameArg) {
    close();
    if (ends_with(fnameArg, BAM_SUFFIX)) {
      bam_ = new BamReader;
      bam_->Open(fnameArg);
      bamAlign_ = new BamAlignment;
      bamDone_ = false;
      format_ = BAM_FORMAT;
      fname_ = fnameArg;
    } else {
      InStream::open(fnameArg);
      format_ = 0;
    }
    totalTagCounter_ = 0;
    rejectTagCounter_ = 0;
    outOfBoundsTagCounter_ = 0;
    confidentTagCounter_ = 0;
    expectedTags_ = 0;
    currentExptName_ = getFnamePrefix(fnameArg);
    align_->forward = true;
    align_->contig = 0;
    align_->firstPos = 0;
    align_->count = 0;
    align_->lastPos = 0;
    align_->seq.clear();
  }

  void ParseAlignStream::error (const string& message) const {
    cerr << "error: " << message << " in " << fname_ << " line " << lineNo_ << "\n" << endl; // BAM doesn't have lines...
    exit(1);
  }

  ParseAlignStream::ParseAlignStream () :
    bam_(0),
    bamAlign_(0),
    bamDone_(false),
    contigs_(0),
    mismatchTolerance_(DEFAULT_MISMATCH_TOLERANCE),
    useLength_(0),
    offset_(0),
    readLength_(0),
    probThreshold_(0),
    phredThreshold_(0),
    prob_(0),
    format_(0),
    totalTagCounter_(0),
    rejectTagCounter_(0),
    outOfBoundsTagCounter_(0),
    confidentTagCounter_(0),
    expectedTags_(0),
    align_(new Alignment)
  {}

  ParseAlignStream::ParseAlignStream (const ContigTable* const contigTableArg, const UShort& mismatchToleranceArg = DEFAULT_MISMATCH_TOLERANCE, const UShort& useLengthArg = 0, const short& offsetArg = 0, const double& probThresholdArg = 0) :
    bam_(0),
    bamAlign_(0),
    bamDone_(false),
    contigs_(contigTableArg),
    mismatchTolerance_(mismatchToleranceArg),
    useLength_(useLengthArg),
    offset_(offsetArg),
    readLength_(useLengthArg),
    probThreshold_(probThresholdArg),
    phredThreshold_(probThresholdArg == 0 ? 0 :
      probThresholdArg == 1 ? 255 :
      -10 * log10(1 - probThresholdArg)),
    prob_(useLengthArg == 0 ? 0 : new BinomPosterior(useLength_)),
    format_(0),
    totalTagCounter_(0),
    rejectTagCounter_(0),
    outOfBoundsTagCounter_(0),
    confidentTagCounter_(0),
    expectedTags_(0),
    align_(new Alignment)
  {}

  //ParseAlignStream::ParseAlignStream (const string& fnameArg, const ContigTable* const contigTableArg, const UShort& mismatchToleranceArg = DEFAULT_MISMATCH_TOLERANCE, const UShort& useLengthArg = 0, const short& offsetArg = 0, const double& probThresholdArg = 0) :
  ParseAlignStream::ParseAlignStream (const string& fnameArg, const ContigTable* const contigTableArg, const UShort& mismatchToleranceArg, const UShort& useLengthArg, const short& offsetArg, const double& probThresholdArg) :
    bam_(0),
    bamAlign_(0),
    bamDone_(false),
    contigs_(contigTableArg),
    mismatchTolerance_(mismatchToleranceArg),
    useLength_(useLengthArg),
    offset_(offsetArg),
    readLength_(useLengthArg),
    probThreshold_(probThresholdArg),
    phredThreshold_(probThresholdArg == 0 ? 0 :
      probThresholdArg == 1 ? 255 :
      -10 * log10(1 - probThresholdArg)),
    prob_(useLengthArg == 0 ? 0 : new BinomPosterior(useLength_)),
    format_(0),
    totalTagCounter_(0),
    rejectTagCounter_(0),
    outOfBoundsTagCounter_(0),
    confidentTagCounter_(0),
    expectedTags_(0),
    align_(new Alignment(true, contigTableArg->size(), 0, 0, "", 0))
  {
    open(fnameArg);
  }

  ParseAlignStream::~ParseAlignStream () {
    close();
    if (bam_) delete bam_;
    if (bamAlign_) delete bamAlign_;
    if (prob_) delete prob_;
    if (align_) delete align_; // added 20101012; did I miss this before, or does it break something?
  }

  void ParseAlignStream::printSummary () const {
    if (! currentExptName_.empty()) cerr << currentExptName_ << ": ";
    cerr << totalTagCounter_ << " tags";
    if (format_ == ELAND_MULTI_FORMAT || format_ == CORONA_FORMAT || format_ == SAM_FORMAT || format_ == BAM_FORMAT) { // must be one that allows probability filtering
      cerr << "\n  " << confidentTagCounter_;
      if (probThreshold_ != 0) cerr << " confidently mapped"; else cerr << " unique best";
      cerr <<
        " hits (" <<
        format("%.1f") % (100 * (double)confidentTagCounter_ / (double)totalTagCounter_) <<
        "%)";
      if (probThreshold_ != 0) cerr <<
        "\n  " <<
        rejectTagCounter_ <<
        " unique best hits rejected by filter (" <<
        format("%.1f") % (100 * (double)rejectTagCounter_ / (double)totalTagCounter_) <<
        "%)";  
    }
    if (outOfBoundsTagCounter_ > 0) cerr << "\n" << outOfBoundsTagCounter_ << " out of contig bounds";
    cerr << endl;
  }

  void ParseAlignStream::setContigTable (const ContigTable* const contigTableArg) {
    contigs_ = contigTableArg;
  }

  const Alignment& ParseAlignStream::lastAlign () const {
    return *align_;
  }

  const Alignment& ParseAlignStream::parseAlign (const string& line) {
    align_->count = 0;
    try {
      if ((line.empty() && format_ != BAM_FORMAT) || line[0] == '#') return *align_; // is # a universal comment?
      
      if (format_ == 0) {
        // determine format
        if (starts_with(line, UCSC_TRACK_PREFIX)) {
          cmatch matches;
          if (regex_search(line.c_str(), matches, NAME_REGEX1) || regex_search(line.c_str(), matches, NAME_REGEX2)) {
            currentExptName_ = matches[1];
          }
          if (regex_search(line.c_str(), WIG_REGEX)) {
            cmatch matches2;
            if (regex_search(currentExptName_.c_str(), matches2, DIRECTIONAL_WIG_NAME_REGEX)) {
              format_ = DIRECTIONAL_WIG_FORMAT;
              currentExptName_ = matches2[1];
              if (matches2[2] == "+") {
                align_->forward = true;
              } else {
                align_->forward = false;
              }
            } else {
              format_ = NONDIRECTIONAL_WIG_FORMAT;
              currentExptName_ = matches[1];
              align_->forward = true; // may not be necessary
            }
          } else {
            format_ = BED_FORMAT;
          }
          return *align_; // break with no return
        } else if (regex_search(line.c_str(), BED_REGEX)) {
          format_ = BED_FORMAT;
        } else if (regex_search(line.c_str(), ELAND_MULTI_REGEX)) {
          format_ = ELAND_MULTI_FORMAT;
        } else if (regex_search(line.c_str(), CORONA_REGEX)) {
          format_ = CORONA_FORMAT;
        } else if (regex_search(line.c_str(), SAM_HEADER_REGEX)) {
          format_ = SAM_FORMAT;
        } else {
          error();
        }
      }

      switch(format_) {
        case BED_FORMAT: {
          if (starts_with(line, UCSC_TRACK_PREFIX)) break; // ?
          vector<string> fields;
          split(fields, line, is_any_of("\t "));
          if (fields.size() < 6) {error();}
          const char parseStrand = fields[5][0];
          if (parseStrand == '+') {
            align_->forward = true;
          } else if (parseStrand == '-') {
            align_->forward = false;
          } else {
            error();
          }
          align_->contig = (*contigs_)[fields[0]];
          align_->seq = fields[3];
          if (align_->forward) {
            align_->firstPos = lexical_cast<Pos>(fields[1].c_str()) + 1;
            align_->lastPos = lexical_cast<Pos>(fields[2].c_str());
          } else {
            align_->firstPos = lexical_cast<Pos>(fields[2].c_str());
            align_->lastPos = lexical_cast<Pos>(fields[1].c_str()) + 1;
          }
          align_->count = 1;
          ++totalTagCounter_;
          break;
        }
        
        case ELAND_MULTI_FORMAT: { // Eland multi
          vector<string> fields;
          split(fields, line, is_any_of("\t"));
          if (fields.size() < 4) error();
          ++totalTagCounter_;
          if (fields[3] == "-") break;
          
          // check if unique best hit and get mismatches
          vector<HitCount> counts;
          UShort bestMismatches = 0;
          bool foundBestHit = false;
          bool unique = false;
          tokenizer<char_separator<char> > mismatchTokenizer(fields[2], (char_separator<char>)(":"));
          for (tokenizer<char_separator<char> >::const_iterator i = mismatchTokenizer.begin(); i != mismatchTokenizer.end(); ++i) {
            counts.push_back(lexical_cast<HitCount>(((string)(*i)).c_str()));
            if (! foundBestHit) {
              if (bestMismatches > mismatchTolerance_) break; // don't count if unique best fails mismatch tolerance
              if (counts.back() == 0) {
                ++bestMismatches;
              } else {
                foundBestHit = true;
                if (counts.back() == 1) {
                  unique = true;
                } else {
                  break; // don't bother to continue if we already know it's not unique
                }
              }
            }
          }
          if (! unique) break;
          
          // get sequence
          align_->seq = fields[1];
          string useSeq = align_->seq; // used only to find mismatches
          if (useLength_ != 0) { // no truncation set
            if (align_->seq.size() >= useLength_) {
              useSeq = align_->seq.substr(0, useLength_);
            } else {
              error("sequence shorter than requested length");
            }
          } else if (readLength_ != 0) { // read length previously established
            if (useSeq.size() != readLength_) error("different read length");
          } else { // need to establish read length
            readLength_ = align_->seq.size();
            if (probThreshold_ != 0) {
              assert(! prob_);
              prob_ = new BinomPosterior(readLength_);
            }
          } 
          const UShort nN = count(useSeq.begin(), useSeq.end(), 'N');
          
          // check probability
          if (prob_ && prob_->prob(bestMismatches, counts) < probThreshold_) {
            align_->count = 0;
            ++rejectTagCounter_;
    //        cerr << line << "\n"; for (vector<HitCount>::const_iterator i = counts.begin(); i != counts.end(); ++i) cerr << *i << " \n"; // test
            break;
          }
           
          // parse alignments
          tokenizer<char_separator<char> > hits(fields[3], (char_separator<char>)(","));
          for (tokenizer<char_separator<char> >::const_iterator i = hits.begin(); i != hits.end(); ++i) {
            string hitString = *i;
            
            // determine contig - if not listed, same as previous
            const size_t whereIsColon = hitString.find(":");
            if (whereIsColon != string::npos) {
              const string contigString = hitString.substr(0, whereIsColon);          
              const size_t whereIsSlash = contigString.find("/");
              align_->contig = (*contigs_)[getFnamePrefix(whereIsSlash == string::npos ? contigString : contigString.substr(whereIsSlash + 1))];
              hitString.erase(0, whereIsColon + 1); 
            }
            
            // determine position
            const size_t whereIsDir = hitString.find_first_of("FR");
            if (whereIsDir == string::npos) {error();}
            const Pos leftPos = lexical_cast<Pos>((hitString.substr(0, whereIsDir)).c_str());
            if (hitString[whereIsDir] == 'F') {
              align_->forward = true;
            } else if (hitString[whereIsDir] == 'R') {
              align_->forward = false;
            } else {
              assert(false);
            }
            hitString.erase(0, whereIsDir + 1);
            
            // determine mismatches
            UShort mismatches = 0;
            UShort lengthElapsed = 0;
            size_t whereIsMismatch = hitString.find_first_of("ACGTN");
            if (whereIsMismatch == string::npos) {
              const UShort value = lexical_cast<UShort>(hitString);
              if (value <= 2) mismatches = value; // if Eland is only reporting mismatch number, not position
            } else {
              while (whereIsMismatch != string::npos) {
                if (whereIsMismatch > 0) {
                  lengthElapsed += lexical_cast<UShort>((hitString.substr(0, whereIsMismatch)).c_str());
                }
                if (lengthElapsed >= useSeq.size()) {break;}
                if (hitString[whereIsMismatch] != 'N') {++mismatches;}
                ++lengthElapsed;
                hitString.erase(0, whereIsMismatch + 1);
                whereIsMismatch = hitString.find_first_of("ACGTN");
              }
              mismatches -= nN; // Eland doesn't count Ns as mismatches, but it reports them
            }

            // check if this is best hit
            if (mismatches == bestMismatches) {
              align_->firstPos = leftPos + (align_->forward ? 0 : readLength_ - 1);
              if (align_->contig < contigs_->size()) align_->count = 1;
              break;
            }
          }
          break;
        }
        
        case CORONA_FORMAT: {
          if (! starts_with(line, ">")) break;
          align_->seq = readLine();
          ++totalTagCounter_;
          
          const size_t alignPos = line.find(",");
          if (alignPos == string::npos) break;
          
          // get sequence
          if (useLength_ == 0) {
            if (readLength_ != 0) { // read length previously established
              if (align_->seq.size() - 1 != readLength_) {
                error("different read length");
              }
            } else { // need to establish read length
              readLength_ = align_->seq.size() - 1;
              if (probThreshold_ != 0) {
                assert(! prob_);
                prob_ = new BinomPosterior(readLength_);
              }
            }
          } else if (align_->seq.size() - 1 < useLength_) {
            error("sequence shorter than requested length");
          }
          
          UShort bestMismatches = mismatchTolerance_ + 1;
          vector<string> hits;
          const string hackString = line.substr(alignPos + 1); // hack
          split(hits, hackString, is_any_of(","));
          vector<HitCount> mismatchCounts(1, 0);
          for (vector<string>::const_iterator i = hits.begin(); i != hits.end(); ++i) {
            vector<string> fields;
            split(fields, *i, is_any_of("."));
            if (fields.size() != 3) {error();}
            const ContigNo contig = (*contigs_)[fields[0]];
            if (contig >= contigs_->size()) break;
            const UShort mismatches = lexical_cast<UShort>(fields[2].c_str());
            
            if (probThreshold_ != 0) {
              if (mismatchCounts.size() >= (UShort)(mismatches + 1)) { // the recast prevents a compiler warning
                ++mismatchCounts[mismatches];
              } else {
                while (mismatchCounts.size() < mismatches) mismatchCounts.push_back(0); // fill in zeroes between last count and this one, if necessary
                mismatchCounts.push_back(1);
              }
            }
            
            if (mismatches < bestMismatches) {
              align_->contig = contig;
              if (! starts_with(fields[1], "-")) {
                align_->forward = true;
                align_->firstPos = lexical_cast<Pos>(fields[1].c_str()) + 1; // really? + 1?
              } else {
                align_->forward = false;
                align_->firstPos = lexical_cast<Pos>(((fields[1]).substr(1)).c_str()) + 1; // again with the + 1?
              }
              bestMismatches = mismatches;
              align_->count = 1;
            } else if (mismatches == bestMismatches) {
              align_->count = 0;
              if (bestMismatches == 0) {break;} // non-unique
            }
          }
          
          if (prob_ && align_->count != 0 && prob_->prob(bestMismatches, mismatchCounts) < probThreshold_) {
            align_->count = 0;
            ++rejectTagCounter_;
    //        cerr << line << "\n"; for (vector<HitCount>::const_iterator i = mismatchCounts.begin(); i != mismatchCounts.end(); ++i) cerr << *i << " \n"; // test
          }
          break;
        }
          
        case DIRECTIONAL_WIG_FORMAT: {
          cmatch matches;
          if (isdigit(line[0])) { // position (this section is "optimized" since it's a bottleneck)
            if (align_->contig == contigs_->size()) break; // contig not defined, or not found in contig table (useful to skip these in case you use a reduced contig table)
            
            // get values
            const size_t delimiter = line.find_last_of("\t "); // why last? because the second number is usually shorter :-)
            if (delimiter == line.npos) error();
            align_->firstPos = lexical_cast<Pos>(line.substr(0, delimiter));
            align_->count = lexical_cast<HitCount>(line.substr(delimiter + (line[delimiter + 1] == '-' ? 2 : 1))); // strip minus sign

            align_->lastPos = align_->firstPos + (useLength_ == 0 ? 0 : (align_->forward ? useLength_ - 1 : -(useLength_ - 1))); // simulates read length if you tell it how long; o'wise zero length isn't the end of the world
            totalTagCounter_ += align_->count;          
          } else if (regex_search(line.c_str(), matches, WIG_CHR_REGEX)) {
            align_->contig = (*contigs_)[matches[1]];
            return *align_; // break with no return
          } else if (starts_with(line, UCSC_TRACK_PREFIX) && regex_search(line.c_str(), matches, NAME_REGEX1)) {
            align_->contig = contigs_->size();
            currentExptName_ = matches[1];
            cmatch matches2;
            // copied from above - make a protected method
            if (regex_search(currentExptName_.c_str(), matches2, DIRECTIONAL_WIG_NAME_REGEX)) {
              currentExptName_ = matches2[1];
              if (matches2[2] == "+") {
                align_->forward = true;
              } else {
                align_->forward = false;
              }
            } else {
              error("strand not defined");
            }
            return *align_; // break with no return
          } else {
            error();
          }
          break;
        }
        
        case NONDIRECTIONAL_WIG_FORMAT: {
          cmatch matches;
          if (isdigit(line[0])) { // position (this section is "optimized" since it's a bottleneck)
            if (align_->contig == contigs_->size()) break; // contig not defined, or not found in contig table (useful to skip these in case you use a reduced contig table)

            // get values
            const size_t delimiter = line.find_last_of("\t "); // why last? because the second number is usually shorter :-)
            if (delimiter == line.npos) error();
            align_->firstPos = lexical_cast<Pos>(line.substr(0, delimiter));
            align_->count = lexical_cast<HitCount>(line.substr(delimiter + (line[delimiter + 1] == '-' ? 2 : 1))); // strip minus sign

            align_->lastPos = align_->firstPos + (useLength_ == 0 ? 0 : useLength_ - 1); // simulates read length if you tell it how long; o'wise zero length isn't the end of the world
            totalTagCounter_ += align_->count;          
          } else if (regex_search(line.c_str(), matches, WIG_CHR_REGEX)) {
            align_->contig = (*contigs_)[matches[1]];
            return *align_; // break with no return
          } else if (starts_with(line, UCSC_TRACK_PREFIX) && regex_search(line.c_str(), matches, NAME_REGEX1)) {
            align_->contig = contigs_->size();
            currentExptName_ = matches[1];
            return *align_; // break with no return
          } else {
            error();
          }
          break;
        }        
        
        case SAM_FORMAT: {
          if (line[0] == '@') break;
          vector<string> fields;
          split(fields, line, is_any_of("\t"));
          if (fields.size() < 10) error();  
           
          // bitwise flag
          const boost::uint16_t flag = lexical_cast<boost::uint16_t>(fields[1]);
          if (flag & 0x0100) break; // non-primary
          ++totalTagCounter_;
          if (flag & (0x0004 + 0x0200)) break; // unmapped or failed QC
          
          // mismatches
          if (fields.size() > 11) {
            int32_t editDistance = 0;
            for (vector<string>::const_iterator i = fields.begin(); i != fields.end(); ++i) {
              cmatch matches;
              if (regex_search(i->c_str(), matches, SAM_EDITDISTANCE_REGEX)) {
                editDistance = lexical_cast<int32_t>(matches[1]);
                break;
              }
            }
            if (editDistance > mismatchTolerance_) break;
          }
          
          // MAPQ
          if (phredThreshold_ != 0 && phredThreshold_ != 255 && lexical_cast<boost::uint16_t>(fields[4]) < phredThreshold_) { // really int8 but that's a char
            ++rejectTagCounter_;
            break;
          }
          
          // contig
          align_->contig = (*contigs_)[fields[2]];
          if (align_->contig == contigs_->size()) break; // contig not found
          
          // rest of info
          align_->seq = fields[9];
          if (useLength_ != 0 && align_->seq.size() < useLength_) error("sequence shorter than requested length");          
          align_->forward = ! (flag & 0x0010);
          align_->firstPos = lexical_cast<Pos>(fields[3]) + (align_->forward ? 0 : align_->seq.size() - 1);
          align_->lastPos = align_->firstPos + (align_->forward ? 1 : -1) * ((useLength_ != 0 ? useLength_ : align_->seq.size()) - 1);
          align_->count = 1;
          break;
        }
        
        case BAM_FORMAT: {
          assert(bam_ && bam_->IsOpen() && bamAlign_);
        
          // puts next alignment in bamAlign_ or breaks on failure (end of input)
          if (! bam_->GetNextAlignment(*bamAlign_)) {
            bamDone_ = true;
            break;
          }

          // check if bad read
          if (! bamAlign_->IsPrimaryAlignment()) break;
          ++totalTagCounter_;
          if (bamAlign_->IsFailedQC() || ! bamAlign_->IsMapped()) break;
          
          // mismatches
          int32_t editDistance;
          if (bamAlign_->GetTag("NM", editDistance) && editDistance > mismatchTolerance_) break;
          
          // MAPQ
          if (phredThreshold_ != 0 && phredThreshold_ != 255 && bamAlign_->MapQuality < phredThreshold_) {
            ++rejectTagCounter_;
            break;
          }
          
          // contig
          align_->contig = (*contigs_)[bam_->GetReferenceData()[bamAlign_->RefID].RefName];
          if (align_->contig == contigs_->size()) break; // contig not found
          
          // rest of info
          align_->seq = bamAlign_->QueryBases;
          align_->forward = ! bamAlign_->IsReverseStrand();
          align_->firstPos = bamAlign_->Position + (align_->forward ? 1 : bamAlign_->Length); // BamAlignment.Position is 0-based even though SAM is 1-based
          align_->lastPos = align_->firstPos + (align_->forward ? 1 : -1) * ((useLength_ != 0 ? useLength_ : bamAlign_->Length) - 1);
          align_->count = 1;
          break;
        }

        default:
          assert(false);
          break;
      }

      if (align_->count != 0 && align_->contig != contigs_->size()) {    
        if (readLength_ != 0) align_->lastPos = align_->firstPos + (align_->forward ? readLength_ - 1 : -(readLength_ - 1));
        
        // apply offset
        if (offset_ != 0) {
          if ((align_->forward && ((int)align_->firstPos > -offset_)) || ((! align_->forward) && ((int)align_->lastPos > offset_))) { // avoid shifting off the left end (right end is already checked later)
            align_->firstPos += (align_->forward ? offset_ : -offset_);
            align_->lastPos += (align_->forward ? offset_ : -offset_);
          } else {
            outOfBoundsTagCounter_ += align_->count;
            align_->count = 0;
            return *align_;
          }
        }
            
        // make sure tag isn't out of bounds
        const Pos contigSize = contigs_->getSize(align_->contig);
        if (align_->firstPos == 0 || align_->firstPos > contigSize || align_->lastPos == 0 || align_->lastPos > contigSize) {
          outOfBoundsTagCounter_ += align_->count;
          align_->count = 0;
          return *align_;
        }
        
        confidentTagCounter_ += align_->count;
      }
    } catch(bad_lexical_cast &) {
      error();
    }
    return *align_;
  }    
  
  bool ParseAlignStream::good () const {
    if (format_ == BAM_FORMAT) return (bam_->IsOpen() && ! bamDone_); else return InStream::good();
  }
  
  string ParseAlignStream::readLine () {
    if (format_ == BAM_FORMAT) return ""; else return InStream::readLine();
  }

  const Alignment& ParseAlignStream::readAlign () {
    if (good()) {
      parseAlign(readLine());
    } else {
      align_->count = 0;
      align_->contig = contigs_->size();
    }
    
    // keep going till you hit a good line (not recursively, since that could make a stack millions of calls deep and apparently that is problematic)
    while ((align_->count == 0 || align_->contig == contigs_->size()) && good()) parseAlign(readLine());
    
    return *align_;
  }

  void ParseAlignStream::rewind () {
    if (fname_ == STDIN_STRING) error("can't rewind STDIN");
    if (format_ == BAM_FORMAT) {
      bam_->Rewind();
    } else {
      close();
      open(fname_); // couldn't easily find a better way for a boost::filtering_istream, and this is easy enough
    }
  }

  const string& ParseAlignStream::exptName () const {
    return currentExptName_;
  }

  const Counter& ParseAlignStream::totalTagCount () const {
    return totalTagCounter_;
  }

  const Counter& ParseAlignStream::rejectTagCount () const {
    return rejectTagCounter_;
  }

  const Counter& ParseAlignStream::outOfBoundsTagCount () const {
    return outOfBoundsTagCounter_;
  }

  const Counter& ParseAlignStream::confidentTagCount () const {
    return confidentTagCounter_;
  }

  const Counter& ParseAlignStream::getExpectedTags () {
    if (expectedTags_ == 0) {
      while (good()) {
        const string line = readLine();
        if ((! line.empty()) && line[0] == '#') {
//          if (offset_ == 0) { // only believe expected tags if no shift
            cmatch matches;
            
            // read tag count header
            if (regex_search(line.c_str(), matches, TAG_COUNT_REGEX)) {
              if (expectedTags_ == 0) expectedTags_ = lexical_cast<Counter>(matches[1]); else error("multiple tag count headers");
              if (expectedTags_ == 0) error("zero tag count"); // in case some joker puts that in the header
//            }
          }
        } else {
          parseAlign(line); // pass through to parser if not the line we're looking for; may also get track name
            
          if (expectedTags_ == 0) { // count tags the hard way
            ParseAlignStream tempParse(fname_, contigs_, mismatchTolerance_, useLength_, offset_, probThreshold_);
            while (tempParse.good()) tempParse.readAlign();
            expectedTags_ = tempParse.confidentTagCount();
          }
          
          break;
        }
      }
    }

    return expectedTags_;
  }    

  const UShort& ParseAlignStream::getFormat () const {
    return format_;
  }


  MultiParse::MultiParse (const UShort& size) : vector<ParseAlignStream*>(size, 0) {}

  bool MultiParse::good () const {
    for (const_iterator i = begin(); i != end(); ++i) {
      if (*i && (*i)->good()) {
        return true;
      }
    }
    return false;
  }


  StrandParseAlignStream::StrandParseAlignStream () : ParseAlignStream() {}

  StrandParseAlignStream::StrandParseAlignStream (const bool& forward) : forward_(forward) {}

  StrandParseAlignStream::StrandParseAlignStream (const bool& forward, const ContigTable* const contigTableArg, const UShort& mismatchToleranceArg = DEFAULT_MISMATCH_TOLERANCE, const UShort& useLengthArg = 0, const short& offsetArg = 0, const double& probThresholdArg = 0) : ParseAlignStream(contigTableArg, mismatchToleranceArg, useLengthArg, offsetArg, probThresholdArg), forward_(forward) {}

  StrandParseAlignStream::StrandParseAlignStream (const bool& forward, const string& fnameArg, const ContigTable* const contigTableArg, const UShort& mismatchToleranceArg, const UShort& useLengthArg, const short& offsetArg, const double& probThresholdArg) : ParseAlignStream(fnameArg, contigTableArg, mismatchToleranceArg, useLengthArg, offsetArg, probThresholdArg), forward_(forward) {}


  void StrandParseAlignStream::setDir (const bool& forward) {
    forward_ = forward;
  }

  const Alignment& StrandParseAlignStream::readAlign () {
    ParseAlignStream::readAlign();
    
    if (format_ == DIRECTIONAL_WIG_FORMAT) { // shortcuts for this format
      if (forward_ && ! align_->forward && align_->count > 0) { // avoids sending the forward-strand parser across the rest of the file if it's already done
        align_->count = 0;
        close();
      }
    }
      
    while (good() && align_->forward != forward_) ParseAlignStream::readAlign();
    
    return *align_;
  }


  void NondirParseAlignStream::open (const string& fnameArg) {
    forwardParse_->open(fnameArg);
    reverseParse_->open(fnameArg);
  }

  NondirParseAlignStream::NondirParseAlignStream () : forwardParse_(new StrandParseAlignStream(true)), reverseParse_(new StrandParseAlignStream(false)) {}

  NondirParseAlignStream::NondirParseAlignStream (const ContigTable* const contigTableArg, const UShort& mismatchToleranceArg = DEFAULT_MISMATCH_TOLERANCE, const UShort& useLengthArg = 0, const short& offsetArg = 0, const double& probThresholdArg = 0) : forwardParse_(new StrandParseAlignStream(true, contigTableArg, mismatchToleranceArg, useLengthArg, offsetArg, probThresholdArg)), reverseParse_(new StrandParseAlignStream(false, contigTableArg, mismatchToleranceArg, useLengthArg, offsetArg, probThresholdArg)) {}

  NondirParseAlignStream::NondirParseAlignStream (const string& fnameArg, const ContigTable* const contigTableArg, const UShort& mismatchToleranceArg, const UShort& useLengthArg, const short& offsetArg, const double& probThresholdArg) : forwardParse_(new StrandParseAlignStream(true, contigTableArg, mismatchToleranceArg, useLengthArg, offsetArg, probThresholdArg)), reverseParse_(new StrandParseAlignStream(false, contigTableArg, mismatchToleranceArg, useLengthArg, offsetArg, probThresholdArg)) {
    open(fnameArg);
  }

  NondirParseAlignStream::~NondirParseAlignStream () {
    delete forwardParse_;
    delete reverseParse_;
  }

  bool NondirParseAlignStream::whichLower () const { // returns true if forward parser is further in genome or at same position, false if reverse is further
    const Alignment& align1 = forwardParse_->lastAlign();
    const Alignment& align2 = reverseParse_->lastAlign();
    
    switch(2 * (align1.count == 0) + (align2.count == 0)) {
      case 0: { // both nonzero
        if (align1.contig == align2.contig) { // same contig
          return align1.firstPos <= align2.firstPos; // default is forward first
        } else {
          return align1.contig < align2.contig;
        }     
      }
      
      case 1: return true;
      case 2: return false;    
      case 3: return true;
      default: assert(false);
    }
  }

  bool NondirParseAlignStream::whichFurther () const { // returns true if forward parser is further in file, false otherwise
    return forwardParse_->getLineNo() >= reverseParse_->getLineNo();
  }

  void NondirParseAlignStream::printSummary () const {
    if (whichFurther()) {
      forwardParse_->printSummary();
    } else {
      reverseParse_->printSummary();
    }
  }

  void NondirParseAlignStream::setContigTable (const ContigTable* const contigTableArg) {
    forwardParse_->setContigTable(contigTableArg);
    reverseParse_->setContigTable(contigTableArg);
  }

  const Alignment& NondirParseAlignStream::lastAlign () const {
    return *align_;  
  }

  const Alignment& NondirParseAlignStream::readAlign () {
    // advance appropriate stream
    if (align_->forward) { // don't re-run the comparison, just check what the last one was
      forwardParse_->readAlign();
      if (forwardParse_->lastAlign().forward && forwardParse_->lastAlign().contig == align_->contig && forwardParse_->lastAlign().firstPos < align_->firstPos) {
        // test
        cerr << forwardParse_->lastAlign().firstPos << "\t" << align_->firstPos << "\n";
        forwardParse_->error("alignments out of order");
      }
    } else {
      reverseParse_->readAlign();
      if (! reverseParse_->lastAlign().forward && reverseParse_->lastAlign().contig == align_->contig && reverseParse_->lastAlign().firstPos < align_->firstPos) reverseParse_->error("alignments out of order");
    }
    
    // make sure both streams have been read
    if (forwardParse_->getLineNo() == 0) forwardParse_->readAlign();
    if (reverseParse_->getLineNo() == 0) reverseParse_->readAlign();
     
    // get new hit
    if (whichLower()) *align_ = forwardParse_->lastAlign(); else *align_ = reverseParse_->lastAlign();

    return *align_;
  }

  void NondirParseAlignStream::rewind () {
    forwardParse_->rewind();
    reverseParse_->rewind();
  }

  const string& NondirParseAlignStream::exptName () const {
    return (whichFurther() ? forwardParse_->exptName() : reverseParse_->exptName());
  }

  const Counter& NondirParseAlignStream::totalTagCount () const {
    return (whichFurther() ? forwardParse_->totalTagCount() : reverseParse_->totalTagCount());
  }

  const Counter& NondirParseAlignStream::rejectTagCount () const {
    return (whichFurther() ? forwardParse_->rejectTagCount() : reverseParse_->rejectTagCount());
  }

  const Counter& NondirParseAlignStream::outOfBoundsTagCount () const {
    return (whichFurther() ? forwardParse_->outOfBoundsTagCount() : reverseParse_->outOfBoundsTagCount());
  }

  const Counter& NondirParseAlignStream::confidentTagCount () const { 
    return (whichFurther() ? forwardParse_->confidentTagCount() : reverseParse_->confidentTagCount());
  }

  const Counter& NondirParseAlignStream::getExpectedTags () {
    if (expectedTags_ == 0) {
      if (whichFurther()) {
        expectedTags_ = forwardParse_->getExpectedTags();
      } else {
        expectedTags_ = reverseParse_->getExpectedTags();
      }
    }
    return expectedTags_;
  }

  const UShort& NondirParseAlignStream::getFormat () const {
    return (whichFurther() ? forwardParse_->getFormat() : reverseParse_->getFormat());
  }


  void FormatOutStream::setContigTable (const ContigTable* const contigTableArg) {
    contigs_ = contigTableArg;
  }

  void FormatOutStream::open (const string& fnameArg) {
    assert(format_ != 0);
    OutStream::open(fnameArg);
    if (format_ == BED_FORMAT) {
      trackHeader();
    }
  }

  void FormatOutStream::open (const string& fnameArg, const UShort& formatArg) {
    format_ = formatArg;
    open(fnameArg);
  }

  void FormatOutStream::open (const string& fnameArg, const string& expt_nameArg, const UShort& formatArg, const UShort& priorityArg) {
    exptName_ = expt_nameArg;
    format_ = formatArg;
    priority_ = priorityArg;
    open(fnameArg);
  }

  void FormatOutStream::open (const string& fnameArg, const string& expt_nameArg, const UShort& formatArg, const UShort& priorityArg, const string& assemblyArg) {
    assembly_ = assemblyArg;
    open(fnameArg, expt_nameArg, formatArg, priorityArg);
  }

  FormatOutStream::FormatOutStream () : contigs_(0), format_(0), priority_(0), tagCounter_(0), outputPeaks_(false), outputCorrs_(false) {}

  FormatOutStream::FormatOutStream (const ContigTable* const contigTableArg) : contigs_(contigTableArg), format_(0), priority_(0), tagCounter_(0), outputPeaks_(false), outputCorrs_(false) {}

  FormatOutStream::FormatOutStream (const ContigTable* const contigTableArg, const string& fnameArg, const UShort& formatArg) : contigs_(contigTableArg), format_(formatArg), priority_(0), tagCounter_(0), outputPeaks_(false), outputCorrs_(false) {
    open(fnameArg);
  }

  FormatOutStream::FormatOutStream (const ContigTable* const contigTableArg, const string& fnameArg, const string& expt_nameArg, const UShort& formatArg, const UShort& priorityArg) : contigs_(contigTableArg), format_(formatArg), priority_(priorityArg), exptName_(expt_nameArg), tagCounter_(0), outputPeaks_(false), outputCorrs_(false) {
    open(fnameArg);
  }

  FormatOutStream::FormatOutStream (const ContigTable* const contigTableArg, const string& fnameArg, const string& expt_nameArg, const UShort& formatArg, const UShort& priorityArg, const string& assemblyArg) : contigs_(contigTableArg), format_(formatArg), priority_(priorityArg), exptName_(expt_nameArg), assembly_(assemblyArg), tagCounter_(0), outputPeaks_(false), outputCorrs_(false) {
    open(fnameArg);
  }

  FormatOutStream::~FormatOutStream () {
    close();
  }
  
  void FormatOutStream::outputPeaks (const bool& arg) {
    outputPeaks_ = arg;
  }
  
  void FormatOutStream::outputCorrs (const bool& arg) {
    outputCorrs_ = arg;
  }

  const Counter& FormatOutStream::tagCount () const {
    return tagCounter_;
  }

  void FormatOutStream::write (const string& stringArg) {
    OutStream::write(stringArg);
  }

  void FormatOutStream::write (const Alignment& hit) {
    assert(good());
    if (hit.count != 0) {
      switch(format_) {
        case BED_FORMAT: {
          if (hit.forward) {
            *output_ << (*contigs_)[hit.contig] << "\t"
              << hit.firstPos - 1 << "\t"
              << hit.lastPos << "\t";
            if (hit.seq.empty()) {
              *output_ << "0";
            } else {
              *output_ << hit.seq;
            }
            *output_ << "\t"
              << "0\t"
              << "+\t"
              << "0\t"
              << "0\t"
              << FORWARD_COLOR;
          } else{
            *output_ << (*contigs_)[hit.contig] << "\t"
              << hit.lastPos - 1 << "\t"
              << hit.firstPos << "\t";
            if (hit.seq.empty()) {
              *output_ << "0";
            } else {
              *output_ << hit.seq;
            }
            *output_ << "\t"
              << "0\t"
              << "-\t"
              << "0\t"
              << "0\t"
              << REVERSE_COLOR;
          }
          *output_ << "\n";
          break;
        }

        case DIRECTIONAL_WIG_FORMAT: {
          if (contig_.empty() || forward_ != hit.forward) {
            forward_ = hit.forward;
            contig_.clear();
            trackHeader();
          }
          
          if (contig_.empty() || contig_ != (*contigs_)[hit.contig]) {
            contig_ = (*contigs_)[hit.contig];
            *output_ << "variableStep chrom=" << contig_ << "\n";
          }
          
          *output_ << hit.firstPos << " ";
          if (! forward_) {*output_ << "-";} // make reverse strand negative
          *output_ << hit.count << "\n";
          break;
        }
        
        case NONDIRECTIONAL_WIG_FORMAT: {
          if (contig_.empty()) {
            contig_.clear();
            trackHeader();
          }
          
          if (contig_.empty() || contig_ != (*contigs_)[hit.contig]) {
            contig_ = (*contigs_)[hit.contig];
            *output_ << "variableStep chrom=" << contig_ << "\n";
          }
          
          *output_ << hit.firstPos << " ";
          *output_ << hit.count << "\n";
          break;
        }        
        
        default: {
          assert(false);
          break;
        }
      }
      tagCounter_ += hit.count;
    }
  }

  void FormatOutStream::write (const PosCount& count) {
    assert(format_ != BED_FORMAT);
    write(Alignment(count.forward, count.contig, count.firstPos, 0, "", count.count));
  }

  void FormatOutStream::write (const PosScore& score) {
    if (score.score != 0) {
      switch(format_) {
        case BED_FORMAT: {
          cerr << "error: can't write profile as BED (" << fname_ << ")\n" << endl;
          exit(1);
          break;
        }
        
        case DIRECTIONAL_WIG_FORMAT: {
          if (contig_.empty() || forward_ != score.forward) {
            forward_ = score.forward;
            contig_.clear();
            trackHeader();
          }
          
          if (contig_.empty() || contig_ != (*contigs_)[score.contig]) {
            contig_ = (*contigs_)[score.contig];
            *output_ << "variableStep chrom=" << contig_ << "\n";
          }
          
          *output_ << score.pos << " ";
          if (! score.forward) {*output_ << "-";} // make reverse strand negative
          *output_ << score.score << "\n";
          break;
        }
        
        case NONDIRECTIONAL_WIG_FORMAT: {
          if (contig_.empty()) {
            contig_.clear();
            trackHeader();
          }
          
          if (contig_.empty() || contig_ != (*contigs_)[score.contig]) {
            contig_ = (*contigs_)[score.contig];
            *output_ << "variableStep chrom=" << contig_ << "\n";
          }
          
          *output_ << score.pos << " " << score.score << "\n";
          break;
        }
        
        default: {
          assert(false);
          break;
        }
      }
    }
  }

  void FormatOutStream::write (const Region& region) {
    assert(format_ == REGIONS_FORMAT && contigs_ && region.contig < contigs_->size());
    *output_ << (*contigs_)[region.contig] << ':';
    if (region.forward) {
      *output_ << region.left << '-' << region.left + region.counts->size() - 1;
    } else {
      *output_ << region.left + region.counts->size() - 1 << '-' << region.left;
    }
    
    if (outputPeaks_) *output_ << "\t" << region.peakPos;
    if (outputCorrs_) *output_ << "\t" << region.strandCorr();
    *output_ << "\t" << format("%.2f") % region.posKurtosis();
    
    const HitCountVec exptSums = region.exptSums();
    
    for (HitCountVec::const_iterator i = exptSums.begin(); i != exptSums.end(); ++i) {
      *output_ << "\t" << *i;
      tagCounter_ += *i;
    }
    
    *output_ << "\n";
  }    

  void FormatOutStream::trackHeader () {
    *output_ << "track name=\"" << exptName_;
    if (format_ == DIRECTIONAL_WIG_FORMAT) {
      if (forward_) {
        *output_ << " +";
      } else {
        *output_ << " -";
      }
    }
    *output_ << "\"";
    
    if (format_ == DIRECTIONAL_WIG_FORMAT)  {
      *output_ << " description=\"";
      if (forward_) {
        *output_ << exptName_;
      } else {
        *output_ << " ";
      }
      *output_ << "\"";
    }
    
    if (priority_ != 0) {
      *output_ << " priority=" << (int)priority_;
    }
        
    *output_ << " visibility=";
    switch(format_) {
      case BED_FORMAT: {
        *output_ << "squish itemRgb=on";
        break;
      }
      
      case DIRECTIONAL_WIG_FORMAT: {
        *output_ << "full type=wiggle_0 alwaysZero=on color=";
        if (forward_) {
          *output_ << FORWARD_COLOR;
        } else {
          *output_ << REVERSE_COLOR << " altColor=" << REVERSE_COLOR;
        }
        break;
      }
      
      case NONDIRECTIONAL_WIG_FORMAT: {
        *output_ << "full type=wiggle_0 alwaysZero=on color=" << NONDIRECTIONAL_COLOR;
        break;
      }
      
      default: {
        assert(false);
        break;
      }
    }
    
    if (! assembly_.empty()) *output_ << " db=" << assembly_;
    *output_ << "\n";
  }

}
