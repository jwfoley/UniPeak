#ifndef PEAKCALL_H
#define PEAKCALL_H

#include <cstdlib>
#include <vector>
#include <deque>

#include "defaults.hpp"
#include "data.hpp"
#include "filterstream.hpp"
#include "format.hpp"
#include "kernel.hpp"

namespace unipeak {

  using namespace std;

  class ProfileBuffer {
    protected:
      struct BufferPos {
        HitCountVec* hits;
        Score forwardScore;
        Score reverseScore;
        BufferPos ();
  //      ~BufferPos ();
      };
      typedef deque<BufferPos> WindowDeq; // fixed to window size
      WindowDeq* const buffer_;
      Region* region_;
      vector<Region*>* const regionsOut_;
      const ScoreVec* const kernel_; // must have odd size; middle element is center
      const UShort bandwidth_;
      const Score regionThreshold_;
      const double kurtosisThreshold_;
      const double corrThreshold_;
      const double hitThreshold_;
      const bool forward_;
      ContigNo contig_;
      Pos bufferPos_; // *center* of buffer
      Pos lastPos_; // last *processed* position (not in buffer anymore)
      bool regionIsHit_; // true if region contains any hits at all (it doesn't always)
      Counter nRegions_;
      Counter nRegionRejects_;
      vector<Counter>* const nTagsInRegions_;
      Counter nContigRegions_;
      Counter nNonzeroPos_; // probably only useful to debug/optimize
      const vector<bool>* const control_;
      const vector<double>* const coeffs_;
      const UShort nExpt_;
      const ContigTable* const contigs_;
      FormatOutStream* const profileOut_;
  //    UShort sampledPositions(const HitCountVec& input);
      void processPosition (const Pos&, const Score&, const Score&, HitCountVec* const);
      void processRegion ();
    
    public:
      ProfileBuffer (
        const Kernel* const, // kernel
        const Score&, // region threshold
        const double&, // kurtosis threshold
        const double&, // strand correlation threshold
        const double&, // hit count threshold
        const bool&, // strand (fixed for buffer's lifetime - create a second buffer for other strand if necessary)
        const vector<bool>&, // indicates whether each sample is a negative control or not; also gives number of samples
        const vector<double>&, // scaling coefficients, empty if not scaling
        const ContigTable* const,
        vector<Region*>* const, // regions out
        FormatOutStream* const = 0 // profile out
      );
      ~ProfileBuffer ();
      const Counter& nRegions () const;
      const Counter& nRegionRejects() const;
      const vector<Counter>& nTagsInRegions () const;
      const Counter& nContigRegions () const;
      Counter flushContig (); // returns number of regions on contig
      void add (const HitCountVec&, const ContigNo&, const Pos&, const bool&);
  };

}

#endif
