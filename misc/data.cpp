#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cassert>
#include <limits>

#include "defaults.hpp"
#include "data.hpp"

namespace unipeak {

  using namespace std;

  double sum (const HitCountVec& hits) {
    double result = 0;
    for (HitCountVec::const_iterator i = hits.begin(); i != hits.end(); ++i) result += *i;
    return result;
  }  

  double mean (const ScoreVec::const_iterator& begin, const ScoreVec::const_iterator& end) {
    double sum = 0;
    for (ScoreVec::const_iterator i = begin; i != end; ++i) sum += *i;
    return sum / (double)(end - begin);
  }

  double sd (const ScoreVec::const_iterator& begin, const ScoreVec::const_iterator& end, const double& x_bar) {
    double ssr = 0;
    for (ScoreVec::const_iterator i = begin; i != end; ++i) ssr += pow(*i - x_bar, 2);
    return sqrt(ssr / ((double)(end - begin) - 1));
  }
  
  double sd (const ScoreVec::const_iterator& begin, const ScoreVec::const_iterator& end) {
    return sd(begin, end, mean(begin, end));
  }

  double corr (const ScoreVec::const_iterator& begin1, const ScoreVec::const_iterator& end1, const ScoreVec::const_iterator& begin2, const ScoreVec::const_iterator& end2) {
    const size_t size = end1 - begin1;
    assert(end2 == begin2 + size); // in this order to avoid warnings about unsigned integers

    const double mean1 = mean(begin1, end1);
    const double mean2 = mean(begin2, end2);
    const double sd1 = sd(begin1, end1, mean1);
    const double sd2 = sd(begin2, end2, mean2);
    
    double ssr = 0;
    ScoreVec::const_iterator iter1 = begin1;
    ScoreVec::const_iterator iter2 = begin2;
    
    for (size_t i = 0; i < size; ++i) {
      ssr += (*iter1 - mean1) * (*iter2 - mean2);
      ++iter1;
      ++iter2;
    }
    
    return ssr / (((double)size - 1) * sd1 * sd2);
  }


  PosCount::PosCount () : firstPos(0), count(0) {}

  PosCount::PosCount (const bool& forwardArg, const ContigNo& contigArg, const Pos& firstPosArg, const HitCount& countArg) : forward(forwardArg), contig(contigArg), firstPos(firstPosArg), count(countArg) {}

  Alignment::Alignment () : PosCount(), lastPos(0) {}

  Alignment::Alignment (const bool& forwardArg, const ContigNo& contigArg, const Pos& firstPosArg, const Pos& lastPosArg, const string& seqArg, const HitCount& countArg) : PosCount(forwardArg, contigArg, firstPosArg, countArg), lastPos(lastPosArg), seq(seqArg) {}

  PosCount Alignment::asPosCount () {
    return PosCount(forward, contig, firstPos, count);
  }


  PosScore::PosScore () : pos(0), score(0) {}

  PosScore::PosScore (const bool& forwardArg, const ContigNo& contigArg, const Pos& posArg, const Score& scoreArg) : forward(forwardArg), contig(contigArg), pos(posArg), score(scoreArg) {}


  Region::Region() : forward(true), contig(0), left(0), peakPos(0), peakScore(0), nExpts(0), counts(new vector<HitCountVec*>), scores(make_pair(new ScoreVec, new ScoreVec)) {}

  Region::Region(const bool& forwardArg, const ContigNo& contigArg, const Pos& leftArg, const UShort& nExptsArg) : forward(forwardArg), contig(contigArg), left(leftArg), peakPos(0), peakScore(0), nExpts(nExptsArg), counts(new vector<HitCountVec*>), scores(make_pair(new ScoreVec, new ScoreVec)) {}

  Region::~Region() {
    if (counts) {
      for (vector<HitCountVec*>::iterator i = counts->begin(); i != counts->end(); ++i) if (*i) delete *i;
      delete counts; // counts can be deleted here
    }
    if (scores.first) delete scores.first;
    if (scores.second) delete scores.second;
  }

  void Region::addPos (HitCountVec* const hits, const Score& forwardScore, const Score& reverseScore) {
    assert(! hits || hits->size() == nExpts);
    counts->push_back(hits);
    scores.first->push_back(forwardScore);
    scores.second->push_back(reverseScore);
    const Score score = forwardScore + reverseScore;
    if (peakPos == 0 || score > peakScore) {
      peakPos = left + scores.first->size() - 1; // assumes scores.first and scores.second are always same length
      peakScore = score;
    }
  }    

  HitCount Region::sum () const {
    HitCount result = 0;
    for (vector<HitCountVec*>::const_iterator i = counts->begin(); i != counts->end(); ++i) {
      if (*i) {
        for (HitCountVec::const_iterator j = (*i)->begin(); j != (*i)->end(); ++j) {
          result += *j;
        }
      }
    }
    return result;
  }     

  vector<HitCount> Region::exptSums () const {
    HitCountVec result(nExpts, 0);
    for (vector<HitCountVec*>::const_iterator i = counts->begin(); i != counts->end(); ++i) {
      if (*i) {
        assert((*i)->size() == result.size());
        HitCountVec::const_iterator inputIter = (*i)->begin();
        HitCountVec::iterator outputIter = result.begin();
        while (inputIter != (*i)->end()) {
          *outputIter += *inputIter;
          ++inputIter;
          ++outputIter;
        }
      }
    }
    return result;
  }

  double Region::posMean () const {
    assert(counts && ! counts->empty());
    HitCount count = 0;
    HitCount sum = 0;
    UShort pos = 0;
    for (vector<HitCountVec*>::const_iterator i = counts->begin(); i != counts->end(); ++i) {
      if (*i) {
        HitCount posCount = 0;
        for (HitCountVec::const_iterator j = (*i)->begin(); j != (*i)->end(); ++j) posCount += *j;
        count += posCount;
        sum += posCount * pos;
      }
      ++pos;
    }  
    return (double)sum / (double)count;
  }

  //double posVariance (const MultiCount& values) { // sample variance of positions of tags in a region
  //  assert (! values.empty());
  //  const double x_bar = posMean(values);
  //  double ssd = 0;
  //  double pos = 0;
  //  double tags = 0;
  //  for (MultiCount::const_iterator i = values.begin(); i != values.end(); ++i) {
  //    ssd += *i * pow(pos - x_bar, 2);
  //    tags += *i;
  //    ++pos;
  //  }
  //    return ssd / (tags - 1);
  //}  

  double Region::posKurtosis () const {
    assert(counts && ! counts->empty());
    double x_bar = posMean();
    HitCount count = 0;
    double sum2 = 0;
    double sum4 = 0;
    UShort pos = 0;
    for (vector<HitCountVec*>::const_iterator i = counts->begin(); i != counts->end(); ++i) {
      if (*i) {
        HitCount posCount = 0;
        for (HitCountVec::const_iterator j = (*i)->begin(); j != (*i)->end(); ++j) posCount += *j;    
        count += posCount;
        sum2 += (double)(posCount) * pow((double)pos - x_bar, 2);
        sum4 += (double)(posCount) * pow((double)pos - x_bar, 4);
      }
      ++pos;
    }
    return ((double)count - 1) * sum4 / pow(sum2, 2);
  }
  
  double Region::strandCorr (const UShort& shift) const {
//    assert(scores.first->size() == scores.second->size() && scores.first->size() > (unsigned)(2 * shift + 3));
    if (scores.first->size() == scores.second->size() && scores.first->size() > (unsigned)(2 * shift + 3)) { // test
      return corr(scores.first->begin(), scores.first->begin() + scores.first->size() - 2 * shift, scores.second->begin() + 2 * shift, scores.second->end());
    } else {
      return numeric_limits<double>::quiet_NaN();
//      cerr << scores.first->size() << " " << scores.second->size() << endl << endl;
//      exit(1);
    }
  }


  ContigTable::ContigTable () : names_(new vector<string>), sizes_(new vector<Pos>), indexLookup_(new LookupTable), genomeSize_(0) {}

  ContigTable::~ContigTable () {
    if (names_) delete names_;
    if (sizes_) delete sizes_;
    if (indexLookup_) delete indexLookup_;
  }

  void ContigTable::add (const string& name, const Pos& size) {
    if (indexLookup_->find(name) == indexLookup_->end()) {
      names_->push_back(name);
      sizes_->push_back(size);
      (*indexLookup_)[name] = names_->size() - 1;
      genomeSize_ += size;
    } else {
      cerr << "error: " << name << " defined twice in contig table" << endl << endl;
      exit(1);
    }
  }

  ContigNo ContigTable::size () const {
    assert(names_->size() == sizes_->size() && names_->size() == indexLookup_->size());
    return names_->size();
  }

  bool ContigTable::empty () const {
    return size() == 0;
  }

  ContigNo ContigTable::operator[] (const string& contig) const {
    const LookupTable::const_iterator iter = indexLookup_->find(contig);
    if (iter == indexLookup_->end()) {
      return size();
    } else {
      return iter->second;
    }
  }

  string ContigTable::operator[] (const ContigNo& indexArg) const {
    assert(indexArg < names_->size());
    return (*names_)[indexArg];
  }

  Pos ContigTable::getSize (const ContigNo& indexArg) const {
    if (indexArg < sizes_->size()) {
      return (*sizes_)[indexArg];
    } else {
      return 0;
    }
  }

  Pos ContigTable::getSize (const string& contig) const {
    return getSize((*this)[contig]);
  }

  Pos ContigTable::getGenomeSize () const {
    return genomeSize_;
  }

  vector<string>::const_iterator ContigTable::begin () const {
    return names_->begin();
  }

  vector<string>::const_iterator ContigTable::end () const {
    return names_->end();
  }

  CountMap::CountMap (const ContigNo& contigArg) : hits_(new ContigCountMap(contigArg, 0), new ContigCountMap(contigArg, 0)), nContigs_(contigArg), tagCounter_(0), endIter_(new ConstIterator), strandEndIter_(new ConstStrandIterator), nondirEndIter_(new ConstNondirIterator) {
    // seems hacky but apparently the only way
    for (ContigCountMap::iterator i = hits_.first->begin(); i != hits_.first->end(); ++i) {
      *i = new PosCountMap;
    }
    for (ContigCountMap::iterator i = hits_.second->begin(); i != hits_.second->end(); ++i) {
      *i = new PosCountMap;
    }
  }
    
  CountMap::~CountMap () {
    if (hits_.first) {
      for (ContigCountMap::const_iterator contig = hits_.first->begin(); contig != hits_.first->end(); ++contig) delete *contig;
      delete hits_.first;
    }
    if (hits_.second) {
      for (ContigCountMap::const_iterator contig = hits_.second->begin(); contig != hits_.second->end(); ++contig) delete *contig;
      delete hits_.second;
    }
    if (endIter_) delete endIter_;
    if (strandEndIter_) delete strandEndIter_;
    if (nondirEndIter_) delete nondirEndIter_;
  }

  CountMap::PosCountMap* CountMap::findPosMap (const ContigNo& contig, const bool& forward) const { // returns null pointer if not found
    if (contig >= nContigs_) {
      return 0;
    } else {
      if (forward) {
        assert(hits_.first->size() > contig); // yes, really >
        return (*hits_.first)[contig];
      } else {
        assert(hits_.second->size() > contig); // yes, really >
        return (*hits_.second)[contig];
      }
    }
  }

  void CountMap::add (const Alignment& hit) {
    if (hit.count && hit.contig < nContigs_) {
      PosCountMap* const posMap_ = findPosMap(hit.contig, hit.forward);
      if (posMap_) {
        PosCountMap::iterator i = posMap_->find(hit.firstPos);
        if (i == posMap_->end()) {
          (*posMap_)[hit.firstPos] = hit.count;
        } else {
          i->second += hit.count;
        }
        tagCounter_ += hit.count;
      }
    }
  }

  const Counter& CountMap::tagCount () const {
    return tagCounter_;
  }


  void CountMap::ConstIterator::getNextContig () {
    assert(! end_ && strandMap_ && contigMap_ && posMap_ && contigIter_ != contigMap_->end());
    while (posIter_ == posMap_->end()) { // make sure we don't just advance to an empty contig
      ++contigIter_;
      if (contigIter_ == contigMap_->end()) {
        if (forward_) {
          forward_ = false;
          result_.forward = false;
          contigMap_ = strandMap_->second;
          contigIter_ = contigMap_->begin();
          contig_ = 0;
          result_.contig = 0;
        } else {
          end_ = true;
          result_.count = 0;
          return; // final exit point
        }
      } else {
        ++contig_;
        ++result_.contig;
      } 
      posMap_ = *contigIter_;
      posIter_ = posMap_->begin();    
      result_.contig = contig_;
    }
  }

  CountMap::ConstIterator::ConstIterator () : forward_(false), strandMap_(0), contigMap_(0), posMap_(0), end_(true) {}

  CountMap::ConstIterator::ConstIterator (const CountMap& parent) : forward_(true), strandMap_(&parent.hits_), contigMap_(strandMap_->first), contigIter_(contigMap_->begin()), contig_(0), posMap_(*contigIter_), posIter_(posMap_->begin()), end_(false) {
    if (posIter_ == posMap_->end()) getNextContig();
    result_.contig = contig_;
    result_.forward = forward_;
    result_.firstPos = posIter_->first;
    result_.count = posIter_->second;
  }

  CountMap::ConstIterator& CountMap::ConstIterator::operator= (const ConstIterator& original) {
    strandMap_ = original.strandMap_;
    contigMap_ = original.contigMap_;
    contigIter_ = original.contigIter_;
    contig_ = original.contig_;
    posMap_ = original.posMap_;
    posIter_ = original.posIter_;
    result_ = original.result_;
    end_ = original.end_;
    return *this;
  }

  bool CountMap::ConstIterator::operator== (const ConstIterator& query) const {
    if (query.end_) {
      return(end_);
    } else {
      return ((! end_) &&
        // danger: doesn't check query.result_
        strandMap_ == query.strandMap_ && 
        contigMap_ == query.contigMap_ &&
        posMap_ == query.posMap_ &&
        contigIter_ == query.contigIter_ &&
        posIter_ == query.posIter_ &&
        contig_ == query.contig_
      );
    }
  }

  bool CountMap::ConstIterator::operator!= (const ConstIterator& query) const {
    return ! operator==(query);
  //  if (query.end_) {
  //    return(! end_);
  //  } else {
  //    return (end_ ||
  //      // danger: doesn't check query.result_
  //      strandMap_ != query.strandMap_ ||
  //      contigMap_ != query.contigMap_ ||
  //      posMap_ != query.posMap_ ||
  //      contigIter_ != query.contigIter_ ||
  //      posIter_ != query.posIter_ ||
  //      contig_ != query.contig_
  //    );
  //  }
  }

  CountMap::ConstIterator& CountMap::ConstIterator::operator++ () {
    assert(! end_ && strandMap_ && contigMap_ && posMap_ && contigIter_ != contigMap_->end() && posIter_ != posMap_->end());
    ++posIter_;
    if (posIter_ == posMap_->end()) {
      getNextContig();
    }
    result_.firstPos = posIter_->first;
    result_.count = posIter_->second;
    return *this;
  }

  const PosCount& CountMap::ConstIterator::operator* () const {
    return result_;
  }

  const PosCount* CountMap::ConstIterator::operator-> () const {
    return &result_;
  }

  CountMap::ConstIterator CountMap::begin () const {
    return ConstIterator(*this);
  }

  const CountMap::ConstIterator& CountMap::end () const {
    return *endIter_;
  }


  CountMap::ConstStrandIterator::ConstStrandIterator () : CountMap::ConstIterator::ConstIterator () {}

  CountMap::ConstStrandIterator::ConstStrandIterator (const CountMap& parent, const bool& forward): ConstIterator(parent) {
    // no initialization list since protected members are in base class, so this is a little hacky for the reverse case
    if (! forward) {
      forward_ = forward;
      contigMap_ = strandMap_->second;
      contigIter_ = contigMap_->begin();
      contig_ = 0;
      posMap_ = *contigIter_;
      posIter_ = posMap_->begin();
      end_ = false;
      if (posIter_ == posMap_->end()) getNextContig();
      result_.contig = contig_; 
      result_.forward = forward_;
      result_.firstPos = posIter_->first;
      result_.count = posIter_->second;
    }
  }

  CountMap::ConstStrandIterator& CountMap::ConstStrandIterator::operator= (const ConstStrandIterator& original) {
    strandMap_ = original.strandMap_;
    contigMap_ = original.contigMap_;
    contigIter_ = original.contigIter_;
    contig_ = original.contig_;
    posMap_ = original.posMap_;
    posIter_ = original.posIter_;
    result_ = original.result_;
    end_ = original.end_;
    return *this;
  }

  bool CountMap::ConstStrandIterator::operator== (const ConstStrandIterator& query) const {
    if (query.end_) {
      return end_;
    } else {
      return ((! end_) &&
        // danger: doesn't check query.result_
        strandMap_ == query.strandMap_ && 
        contigMap_ == query.contigMap_ &&
        posMap_ == query.posMap_ &&
        contigIter_ == query.contigIter_ &&
        posIter_ == query.posIter_ &&
        contig_ == query.contig_
      );
    }
  }

  bool CountMap::ConstStrandIterator::operator!= (const ConstStrandIterator& query) const {
    return ! operator==(query);
  }

  CountMap::ConstStrandIterator CountMap::forwardBegin () const {
    return ConstStrandIterator(*this, true);
  }

  CountMap::ConstStrandIterator CountMap::reverseBegin () const {
    return ConstStrandIterator(*this, false);
  }

  CountMap::ConstStrandIterator CountMap::forwardEnd () const {
    return ConstStrandIterator(*this, false);
  }

  const CountMap::ConstStrandIterator& CountMap::reverseEnd () const {
    return *strandEndIter_;
  }


  UShort CountMap::ConstNondirIterator::which () const {
    assert(parent_);
    switch(2 * (*forwardIter_ == parent_->forwardEnd()) + (*reverseIter_ == parent_->reverseEnd())) {
      case 0: { // neither done
        if ((*forwardIter_)->contig == (*reverseIter_)->contig) { // same contig
          if ((*forwardIter_)->firstPos == (*reverseIter_)->firstPos) {
            return 2;
          } else if ((*forwardIter_)->firstPos < (*reverseIter_)->firstPos) {
            return 0;
          } else if ((*forwardIter_)->firstPos > (*reverseIter_)->firstPos) {
            return 1;
          } else assert(false);
        
        } else if ((*forwardIter_)->contig < (*reverseIter_)->contig) { 
          return 0;
        } else if ((*forwardIter_)->contig > (*reverseIter_)->contig) {
          return 1;
        } else assert(false);
      }
      
      case 1: return 0;
      case 2: return 1;
      case 3: return 2;
      default: assert(false);
    }
  }

  CountMap::ConstNondirIterator::ConstNondirIterator () : parent_(0), forwardIter_(0), reverseIter_(0), end_(true) {}

  CountMap::ConstNondirIterator::ConstNondirIterator (const CountMap& parentArg) : parent_(&parentArg), forwardIter_(new CountMap::ConstStrandIterator(parentArg, true)), reverseIter_(new CountMap::ConstStrandIterator(parentArg, false)), end_(*forwardIter_ == parentArg.forwardEnd() && *reverseIter_ == parentArg.reverseEnd()) {} // remember this could still return the end

  CountMap::ConstNondirIterator::~ConstNondirIterator () {
    if (forwardIter_) delete forwardIter_;
    if (reverseIter_) delete reverseIter_;
  }

  bool CountMap::ConstNondirIterator::operator== (const ConstNondirIterator& query) const {
    if (query.end_) {
      return end_;
    } else {
      return *forwardIter_ == *query.forwardIter_ && *reverseIter_ == *query.reverseIter_;
    }
  }

  bool CountMap::ConstNondirIterator::operator!= (const ConstNondirIterator& query) const {
    return ! operator==(query);
  }

  CountMap::ConstNondirIterator& CountMap::ConstNondirIterator::operator++ () {
    assert(! end_);
    switch(which()) {
      case 0: {
        ++(*forwardIter_);
        break;
      }
      
      case 1: {
        ++(*reverseIter_);
        break;
      }
      
      case 2: {
        ++(*forwardIter_);
        ++(*reverseIter_);
        break;
      }
      
      default: assert(false);
    }
    if (*forwardIter_ == parent_->forwardEnd() && *reverseIter_ == parent_->reverseEnd()) end_ = true;
    return *this;
  }

  PosCount CountMap::ConstNondirIterator::operator* () const {
    switch(which()) {
      case 0: return **forwardIter_;
      case 1: {
        PosCount result = **reverseIter_;
        result.forward = true;
        return result;
      }
      case 2: {
        PosCount result = **forwardIter_;
        result.count += (**reverseIter_).count;
        return result;
      }
      default: assert(false);
    }
  }

  CountMap::ConstNondirIterator CountMap::nondirBegin () const {
    return ConstNondirIterator(*this);
  }

  const CountMap::ConstNondirIterator& CountMap::nondirEnd () const {
    return *nondirEndIter_;
  }

}
