#ifndef DATA_H
#define DATA_H

#include <cstdlib>
#include <string>
#include <vector>
#include <map>

#include <boost/unordered_map.hpp>

#include "defaults.hpp"

namespace unipeak {

  using namespace std;

  // counts in more than one experiment
  typedef vector<Score> ScoreVec;
  typedef vector<HitCount> HitCountVec;
  
  double sum (const HitCountVec&);
  double mean (const ScoreVec::const_iterator&, const ScoreVec::const_iterator&); // mean of a vector between a begin and end iterator
  double sd (const ScoreVec::const_iterator&, const ScoreVec::const_iterator&); // standard deviation of a vector between a begin and end
  double sd (const ScoreVec::const_iterator&, const ScoreVec::const_iterator&, const double&); // slightly faster calculation if you already have the mean
  double corr (const ScoreVec::const_iterator&, const ScoreVec::const_iterator&, const ScoreVec::const_iterator&, const ScoreVec::const_iterator&); // Pearson correlation between two vectors (begin1, end1, begin2, end2); must be the same length, of course

  // compact version of Alignment for only the useful parts
  struct PosCount {
    PosCount ();
    PosCount (const bool&, const ContigNo&, const Pos&, const HitCount&);
    bool forward; // true = forward, false = reverse
    ContigNo contig;
    Pos firstPos;
    HitCount count;
  };

  struct Alignment : PosCount {
    Alignment ();
    Alignment (const bool&, const ContigNo&, const Pos&, const Pos&, const string&, const HitCount&);  // forward, contig, firstPos, lastPos, seq, count
    PosCount asPosCount ();
    Pos lastPos;
    string seq;
  };

  struct PosScore {
    PosScore ();
    PosScore (const bool&, const ContigNo&, const Pos&, const Score&);
    bool forward;
    ContigNo contig;
    Pos pos;
    Score score;
  };

  struct Region {
    bool forward;
    ContigNo contig;
    Pos left;
    Pos peakPos;
    Score peakScore;
    UShort nExpts;
    vector<HitCountVec*>* counts;
    pair<ScoreVec*, ScoreVec*> scores;
    Region ();
    Region (const bool&, const ContigNo&, const Pos&, const UShort&); // last arg is number of expts
    ~Region ();
    void addPos (HitCountVec* const, const Score&, const Score&); // hits, forward score, reverse score
    HitCount sum () const; // returns zero if region undefined, though other functions will die
    vector<HitCount> exptSums () const;
    double posMean () const;
  //  double posVariance () const;
    double posKurtosis () const; // non-excess sample kurtosis of positions of tags in a region; returns nan if applicable
    double strandCorr (const UShort& = 0) const; // correlation between strands' density profiles; optional strand shift; returns NaN if less than 3 positions being compared (NaN returns false in all comparisons)
  };

  class ContigTable {
    protected:
      typedef boost::unordered_map<string, ContigNo> LookupTable;

    public:
      ContigTable ();
      ~ContigTable ();
      void add (const string&, const Pos&); // name, size
      ContigNo size () const;
      bool empty () const;
      ContigNo operator[] (const string&) const; // get contig index by name; returns size() if not found
      string operator[] (const ContigNo&) const; // get contig name by index
      Pos getSize (const ContigNo&) const; // size of contig by index
      Pos getSize (const string&) const; // size of contig by name
      Pos getGenomeSize () const;
      vector<string>::const_iterator begin () const; // for iterating over names in order (maybe this class ought to have its own iterator, or at least a typedef so it isn't fragile)
      vector<string>::const_iterator end () const;
    
    protected:
      vector<string>* const names_;
      vector<Pos>* const sizes_;
      LookupTable* const indexLookup_;
      Pos genomeSize_;
  };

  class CountMap {
    protected:
      typedef map<Pos, HitCount> PosCountMap;
      typedef vector<PosCountMap*> ContigCountMap;
      typedef pair<ContigCountMap*, ContigCountMap*> StrandCountMap;

    public:
      CountMap (const ContigNo&); // number of contigs
      ~CountMap ();
      void add (const Alignment&);
      const Counter& tagCount () const;
      class ConstIterator { // constant "forward" iterator: increments position, then contig, then strand
        public:
          ConstIterator (); // sets to end
          ConstIterator (const CountMap&); // sets to begin
          ConstIterator& operator= (const ConstIterator&);
          bool operator== (const ConstIterator&) const;
          bool operator!= (const ConstIterator&) const;
          ConstIterator& operator++ ();
          const PosCount& operator* () const;
          const PosCount* operator-> () const;
        protected:
          bool forward_;
          const StrandCountMap* strandMap_;
          const ContigCountMap* contigMap_;
          ContigCountMap::const_iterator contigIter_;
          ContigNo contig_;
          const PosCountMap* posMap_;
          PosCountMap::const_iterator posIter_;
          PosCount result_; // this is probably unnecessary
          bool end_; // shortcut
          void getNextContig ();
      };
      ConstIterator begin () const; // not an iterator& because you might use it in an expression
      const ConstIterator& end () const; // a const iterator& so you can't use it in an expression, but doesn't need to be recomputed every time you check for != this->end()
      
      class ConstStrandIterator : public ConstIterator { // like ConstIterator but only iterates over one strand
        public:
          ConstStrandIterator (); // sets to end
          ConstStrandIterator (const CountMap&, const bool&); // beginning of selected strand
          ConstStrandIterator& operator= (const ConstStrandIterator&);
          bool operator== (const ConstStrandIterator&) const;
          bool operator!= (const ConstStrandIterator&) const;
      };
      ConstStrandIterator forwardBegin () const;
      ConstStrandIterator reverseBegin () const;
      ConstStrandIterator forwardEnd () const; // same as reverseBegin, just for ease of internal programming
      const ConstStrandIterator& reverseEnd () const;
      
      class ConstNondirIterator { // increments strand, then position, then contig (basically strandless)
        public:
          ConstNondirIterator (); // sets to end
          ConstNondirIterator (const CountMap&); // sets to begin
          ~ConstNondirIterator ();
          bool operator== (const ConstNondirIterator&) const;
          bool operator!= (const ConstNondirIterator&) const;
          ConstNondirIterator& operator++ ();
          PosCount operator* () const; // can't directly point anywhere if the two need to be combined
  //        const PosCount* operator-> () const; // this is problematic since operator* doesn't point anywhere
        protected:
          const CountMap* const parent_;
          CountMap::ConstStrandIterator* const forwardIter_;
          CountMap::ConstStrandIterator* const reverseIter_;
          UShort which () const;
          bool end_;
      };
      ConstNondirIterator nondirBegin () const;
      const ConstNondirIterator& nondirEnd () const;
      
    protected:
      StrandCountMap hits_;
      const ContigNo nContigs_;
      PosCountMap* findPosMap (const ContigNo&, const bool&) const;
      Counter tagCounter_; 
      const ConstIterator* const endIter_;
      const ConstStrandIterator* const strandEndIter_;
      const ConstNondirIterator* const nondirEndIter_;
  };

}

#endif
