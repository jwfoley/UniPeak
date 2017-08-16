#include <cstdlib>
#include <vector>
#include <deque>
#include <cassert>

#include "defaults.hpp"
#include "data.hpp"
#include "filterstream.hpp"
#include "format.hpp"
#include "peakcall.hpp"

namespace unipeak {

  using namespace std;

  double kernelFunction (const double& x) {
    return 3 * (1 - pow(x, 2)) / 4; // Epanechnikov
  }

  //UShort ProfileBuffer::sampledPositions (const HitCountVec& input) {
  //  UShort result = 0;
  //  for (HitCountVec::const_iterator i = input.begin(); i != input.end(); ++i) if (*i) ++result; // counts non-null pointers, on the assumption that only nonzero positions get pointers
  //  return result;
  //}

  ProfileBuffer::BufferPos::BufferPos () : hits(0), forwardScore(0), reverseScore(0) {}

  // removed because it gets deleted elsewhere but be careful!
  //ProfileBuffer::BufferPos::~BufferPos () {
  //  if (hits) delete hits;
  //}

  void ProfileBuffer::processRegion () {
    const HitCountVec exptSums = region_->exptSums();
    // check non-control sum
    HitCount nonControlSum = 0;
    assert(exptSums.size() == control_->size());
    HitCountVec::const_iterator exptSumIter;
    vector<bool>::const_iterator controlIter;
    for (exptSumIter = exptSums.begin(), controlIter = control_->begin(); exptSumIter != exptSums.end() && controlIter != control_->end(); ++exptSumIter, ++controlIter) if (! *controlIter) nonControlSum += *exptSumIter;
    if (nonControlSum >= hitThreshold_ && (kurtosisThreshold_ == 0 || (region_->counts->size() > 1 && region_->posKurtosis() <= kurtosisThreshold_)) && (corrThreshold_ <= -1 || region_->strandCorr() >= corrThreshold_)) {
      assert(nTagsInRegions_->size() == exptSums.size());
      vector<Counter>::iterator tagsIter;   
      for (exptSumIter = exptSums.begin(), tagsIter = nTagsInRegions_->begin(); exptSumIter != exptSums.end() && tagsIter != nTagsInRegions_->end(); ++exptSumIter, ++tagsIter) *tagsIter += *exptSumIter;
      regionsOut_->push_back(region_); // remember to delete these!
      ++nRegions_;
      ++nContigRegions_;
    } else {
      delete region_; // garbage collection
      ++nRegionRejects_;
    }
    region_ = new Region(forward_, contig_, 0, nExpt_);
  }

  void ProfileBuffer::processPosition (const Pos& pos, const Score& forwardScore, const Score& reverseScore, HitCountVec* const counts) {
    assert(pos > lastPos_);
    const Score score = forwardScore + reverseScore;
    if (pos == lastPos_ + 1) { // single step
      if (region_->left != 0) { // currently in region
        if (score >= regionThreshold_) {
          region_->addPos(counts, forwardScore, reverseScore);
        } else {
          processRegion();
          delete counts; // counts can be deleted here
        }
      } else { // not in region
        if (score >= regionThreshold_) {
          region_->left = pos;
          region_->addPos(counts, forwardScore, reverseScore);
        } else {
          delete counts; // counts can be deleted here
        }
      }
      ++lastPos_;
    
    } else { // leap (over zeroes)
      if (region_->left != 0) processRegion();
      if (score >= regionThreshold_) region_->addPos(counts, forwardScore, reverseScore); else delete counts; // counts can be deleted here
    }
    
    if (score != 0) {
      if (profileOut_) profileOut_->write(PosScore(forward_, contig_, pos, score));
      ++nNonzeroPos_;
    }
    lastPos_ = pos;
  }

  ProfileBuffer::ProfileBuffer (
    const Kernel* const kernelArg, 
    const Score& regionThresholdArg, 
    const double& kurtosisThresholdArg,
    const double& corrThresholdArg, 
    const double& hitThresholdArg, 
    const bool& forwardArg,
    const vector<bool>& controlArg,
    const vector<double>& coeffArg,
    const ContigTable* const contigTableArg,
    vector<Region*>* const regionsOutArg,
    FormatOutStream* const profileOutArg
  ) : 
  //  buffer_(new WindowDeq(kernelArg->size(), make_pair(0x0, 0))),
    buffer_(new WindowDeq(kernelArg->size(), BufferPos())),
    region_(new Region(true, 0, 0, controlArg.size())),
    regionsOut_(regionsOutArg),
    kernel_(kernelArg),
    bandwidth_((kernelArg->size() - 1) / 2),
    regionThreshold_(regionThresholdArg),
    kurtosisThreshold_(kurtosisThresholdArg),
    corrThreshold_(corrThresholdArg),
    hitThreshold_(hitThresholdArg), // don't include controls
    forward_(forwardArg),
    contig_(0),
    bufferPos_(0),
    lastPos_(0),
    nRegions_(0),
    nRegionRejects_(0),
    nTagsInRegions_(new vector<Counter>(controlArg.size(), 0)),
    nContigRegions_(0),
    nNonzeroPos_(0),
    control_(new vector<bool>(controlArg)),
    coeffs_(new vector<double>(coeffArg)),
    nExpt_(controlArg.size()),
    contigs_(contigTableArg),
    profileOut_(profileOutArg) 
  {
    assert(kernelArg->size() > 0 && kernelArg->size() % 2 == 1);
    UShort nControl = 0;
    for (vector<bool>::const_iterator i = control_->begin(); i != control_->end(); ++i) if (*i) ++nControl;
    assert(coeffs_->empty() || nControl + coeffs_->size() == nExpt_);
//    for (WindowDeq::iterator i = buffer_->begin(); i != buffer_->end(); ++i) { // no, this is done by BufferPos::BufferPos()
//      i->hits = 0;
//      i->forwardScore = 0;
//      i->reverseScore = 0;
//    }
  }

  ProfileBuffer::~ProfileBuffer () {
    if (buffer_) delete buffer_;
    if (region_) delete region_;
    if (nTagsInRegions_) delete nTagsInRegions_;
    if (control_) delete control_;
    if (coeffs_) delete coeffs_;
  }

  const Counter& ProfileBuffer::nRegions () const {
    return nRegions_;
  }

  const Counter& ProfileBuffer::nRegionRejects () const {
    return nRegionRejects_;
  }

  const vector<Counter>& ProfileBuffer::nTagsInRegions () const {
    return *nTagsInRegions_;
  }

  const Counter& ProfileBuffer::nContigRegions () const {
    return nContigRegions_;
  }

  void ProfileBuffer::add (const HitCountVec& counts, const ContigNo& contig, const Pos& pos, const bool& forward) {
    assert(counts.empty() || counts.size() == nExpt_);

    if (contig != contig_) {
      flushContig();
      contig_ = contig;
      region_->contig = contig;
    }
        
    // determine overlap of old buffer and new
    assert(pos >= bufferPos_);
    UShort nStaticPos = buffer_->size(); // number of positions that are done being updated
    if (pos <= bufferPos_ + 2 * bandwidth_) nStaticPos = pos - bufferPos_; // check this
    
    // process finished positions and shift buffers
    if (bufferPos_ != 0) { // don't bother with empty starting buffer
      for (UShort i = 0; i < nStaticPos; ++i) {
        if (bufferPos_ + i > bandwidth_) {
          processPosition(bufferPos_ + i - bandwidth_, buffer_->front().forwardScore, buffer_->front().reverseScore, buffer_->front().hits); // don't go off left end; should never get past right end (position would have to be right end + bandwidth)
        buffer_->pop_front(); // kept in this order, avoids resizing the deque, I think
        buffer_->push_back(BufferPos());
        }
      }
    }
    
    // get non-control count sum
    double countSum = 0;
    HitCountVec::const_iterator countIter;
    vector<bool>::const_iterator controlIter;
    vector<double>::const_iterator coeffIter;
    if (coeffs_->empty()) { // unscaled version
      for (countIter = counts.begin(), controlIter = control_->begin(); countIter != counts.end() && controlIter != control_->end(); ++countIter, ++controlIter) if (! *controlIter) countSum += *countIter;
    } else { // scaled version
      for (countIter = counts.begin(), controlIter = control_->begin(), coeffIter = coeffs_->begin(); countIter != counts.end() && controlIter != control_->end() && coeffIter != coeffs_->end(); ++countIter, ++controlIter) if (! *controlIter) {
        countSum += (double)(*countIter) * *coeffIter;
        ++coeffIter;
      }
    }
    
    for (countIter = counts.begin(), controlIter = control_->begin(), coeffIter = coeffs_->begin(); countIter != counts.end() && controlIter != control_->end() && coeffIter != coeffs_->end(); ++countIter, ++controlIter) if (! *controlIter) countSum += *countIter;
          
    // add density from new position
    if (countSum != 0) {
      assert(buffer_->size() == kernel_->size());
      Kernel::const_iterator kernelIter;
      WindowDeq::iterator bufferIter;
      for (kernelIter = kernel_->begin(), bufferIter = buffer_->begin(); kernelIter != kernel_->end() && bufferIter != buffer_->end(); ++kernelIter, ++bufferIter) {
        if (forward) bufferIter->forwardScore += *kernelIter * countSum; else bufferIter->reverseScore += *kernelIter * countSum;
      }
      if ((*buffer_)[bandwidth_].hits) { // this position already hit
        assert((*buffer_)[bandwidth_].hits->size() == counts.size());
        HitCountVec::iterator resultIter = (*buffer_)[bandwidth_].hits->begin();
        for (HitCountVec::const_iterator thisIter = counts.begin(); thisIter != counts.end(); ++thisIter) {
          *resultIter += *thisIter;
          ++resultIter;
        }
      } else {
        (*buffer_)[bandwidth_].hits = new HitCountVec(counts); // this is where it finally gets created (no need to allocate new vectors in the program that calls this); be sure it's later deleted in all possible outcomes
      }
    }
    bufferPos_ = pos;
  }

  Counter ProfileBuffer::flushContig () {
    add(HitCountVec(), contig_, bufferPos_ + buffer_->size(), true); // effectively fills the whole thing with zeroes
    bufferPos_ = 0;
    lastPos_ = 0;
    const Counter result = nContigRegions_;
    nContigRegions_ = 0;
    return result;
  }

}
