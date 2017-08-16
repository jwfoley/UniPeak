#ifndef FORMAT_H
#define FORMAT_H

#include <cstdlib>
#include <string>

#include <boost/regex.hpp>

#include "bamtools/BamReader.h"

#include "defaults.hpp"
#include "data.hpp"
#include "filterstream.hpp"

namespace unipeak {

  using namespace std;
  using namespace boost;
  using namespace BamTools;

  // note: regexes are used for identifying format, not for parsing, because that would be slow.
  const regex CONTIG_TABLE_REGEX("^(\\w+)\\W+(\\d+)"); // ^?
  const string UCSC_TRACK_PREFIX("track");
  const regex BED_REGEX("^[^\\s]+\\s\\d+\\s\\d+\\s[^\\s]*\\s[^\\s]*?\\s[+-]");
  const regex ELAND_MULTI_REGEX("^>.+?\\t[ACGTN\\.]+\\t(\\d+:\\d+:\\d+|RM|NM|QC)\\t.+");
  const regex CORONA_REGEX("^>\\d+_\\d+_\\d+_F3");
  const regex WIG_REGEX("type=wiggle_0");
  const regex WIG_CHR_REGEX("^variableStep chrom=(.+)$");
  const regex DIRECTIONAL_WIG_NAME_REGEX("(.+) ([+-])");
  const regex SAM_HEADER_REGEX("^@[A-Za-z][A-Za-z](\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$"); // from spec v1.4-r985
  const regex SAM_EDITDISTANCE_REGEX("^NM:i:(\\d+)$");
  const regex NAME_REGEX1("name=\"(.+?)\"");
  const regex NAME_REGEX2("name=(.+?) ");
  const regex TAG_COUNT_REGEX("# tags=(\\d+)");

  const UShort BED_FORMAT = 1;
  const UShort ELAND_MULTI_FORMAT = 2;
  const UShort CORONA_FORMAT = 3;
  const UShort SAM_FORMAT = 4;
  const UShort BAM_FORMAT = 5;
  const UShort DIRECTIONAL_WIG_FORMAT = 6;
  const UShort NONDIRECTIONAL_WIG_FORMAT = 7;
  const UShort REGIONS_FORMAT = 8;

  string getFnamePrefix (const string&);

  ContigTable parseContigTable (const string&); // filename

  class BinomPosterior { // calculate posterior probability of unique best hit being correct given list of other mappings (uses binomial error model)
    public:
      BinomPosterior (const UShort&); // read length
      ~BinomPosterior ();
      double prob (const UShort&, const vector<HitCount>&) const; // number of mismatches in candidate alignment, total number of alignments per number of mismatches
    
    protected:
      vector<double>* const coefficients_;
  };

  // could this store a ContigTable& instead of a ContigTable*?
  class ParseAlignStream : public InStream {
    public:
      ParseAlignStream ();
      ParseAlignStream (const ContigTable* const, const UShort&, const UShort&, const short&, const double&); // contig table (uses original, not copy), mismatchTolerance_, useLength_, offset_, probThreshold_
      ParseAlignStream (const string&, const ContigTable* const, const UShort& = DEFAULT_MISMATCH_TOLERANCE, const UShort& = 0, const short& = 0, const double& = 0);
      ~ParseAlignStream ();
      virtual void setContigTable (const ContigTable* const);
      virtual void open (const string&);
      virtual void close ();
      virtual bool good () const;
      virtual string readLine ();
      virtual const Alignment& parseAlign (const string&); // parse a given string
      virtual const Alignment& readAlign (); // read a line from the file and parse it - continue until a valid hit count is found, or end of file, i.e. lastAlign().count = 0 only before reading any alignments and after reading the entire file
      virtual const Alignment& lastAlign () const; // return the last alignment
      virtual void rewind ();
      virtual const string& exptName () const;
      virtual const Counter& totalTagCount () const;
      virtual const Counter& rejectTagCount () const;
      virtual const Counter& outOfBoundsTagCount () const;
      virtual const Counter& confidentTagCount () const;
      virtual const Counter& getExpectedTags (); // return expected tag count; find it if necessary
      virtual const UShort& getFormat () const;
      virtual void printSummary () const;
      virtual void error (const string& = "bad format") const;
      
    protected:
      BamReader* bam_;
      BamAlignment* bamAlign_;
      bool bamDone_;
      void parseBam ();
      const ContigTable* contigs_;
      const UShort mismatchTolerance_;
      const UShort useLength_;
      const short offset_;
      UShort readLength_;
      const double probThreshold_;
      const double phredThreshold_;
      const BinomPosterior* prob_; // null pointer if threshold is zero
      UShort format_;
      Counter totalTagCounter_;
      Counter rejectTagCounter_;
      Counter outOfBoundsTagCounter_;
      Counter confidentTagCounter_;
      Counter expectedTags_;
      string currentExptName_;
      Alignment* align_;
  };

  class MultiParse : public vector<ParseAlignStream*> {
    public:
      MultiParse (const UShort&);
      bool good () const;
  };

  class StrandParseAlignStream : public ParseAlignStream { // only returns alignments from one strand
    public:
      StrandParseAlignStream ();
      StrandParseAlignStream (const bool&);
      StrandParseAlignStream (const bool&, const ContigTable* const, const UShort&, const UShort&, const short&, const double&);
      StrandParseAlignStream (const bool&, const string&, const ContigTable* const, const UShort& = DEFAULT_MISMATCH_TOLERANCE, const UShort& = 0, const short& = 0, const double& = 0);
      void setDir (const bool&);
      const Alignment& readAlign ();
    protected:
      bool forward_;
  };

  class NondirParseAlignStream : public ParseAlignStream { // combines alignments from both strands (opens two filehandles)
    public:
      NondirParseAlignStream ();
      NondirParseAlignStream (const ContigTable* const, const UShort&, const UShort&, const short&, const double&);
      NondirParseAlignStream (const string&, const ContigTable* const, const UShort& = DEFAULT_MISMATCH_TOLERANCE, const UShort& = 0, const short& = 0, const double& = 0);
      ~NondirParseAlignStream ();
      void setContigTable (const ContigTable* const);   
      void open (const string&);
      const Alignment& readAlign ();
      const Alignment& lastAlign () const;
      void rewind ();
      const string& exptName () const;
      const Counter& totalTagCount () const;
      const Counter& rejectTagCount () const;
      const Counter& outOfBoundsTagCount () const;
      const Counter& confidentTagCount () const;
      const Counter& getExpectedTags (); // return expected tag count; find it if necessary (returns to same alignment)
      const UShort& getFormat () const;
      void printSummary () const;    
      
    protected:
      StrandParseAlignStream* const forwardParse_;
      StrandParseAlignStream* const reverseParse_;
      bool whichLower () const;
      bool whichFurther () const;
  };
      

  class FormatOutStream : public OutStream {
    public:
      FormatOutStream ();
      FormatOutStream (const ContigTable* const);
      FormatOutStream (const ContigTable* const, const string&, const UShort&); // contig table, fname, format
      FormatOutStream (const ContigTable* const, const string&, const string&, const UShort&, const UShort&); // contig table (uses original, not copy), fname, exptName, format, priority
      FormatOutStream (const ContigTable* const, const string&, const string&, const UShort&, const UShort&, const string&); // same, then assembly
      ~FormatOutStream ();
      void setContigTable (const ContigTable* const);
      void open (const string&); // filename
      void open (const string&, const UShort&); // filename, format
      void open (const string&, const string&, const UShort&, const UShort&); // filename, expt name, format, priority
      void open (const string&, const string&, const UShort&, const UShort&, const string&); // same, then assembly
      void write (const string&); // this doesn't inherit because it's overloaded
      void write (const Alignment&);
      void write (const PosCount&);
      void write (const PosScore&);
      void write (const Region&);
      void outputPeaks (const bool&); // write peaks along with regions?
      void outputCorrs (const bool&); // write correlations along with kurtosis?
      const Counter& tagCount () const;
      
    protected: 
      const ContigTable* contigs_;
      UShort format_;
      UShort priority_;
      string exptName_;
      string assembly_;
      bool forward_;
      string contig_; // or should this still be index. probably should.
      Counter tagCounter_;
      void trackHeader ();
      bool outputPeaks_;
      bool outputCorrs_;
  }; 

}      

#endif

