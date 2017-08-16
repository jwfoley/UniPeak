#include <cstdlib>
#include <vector>
#include <cassert>

#include "defaults.hpp"
#include "kernel.hpp"

namespace unipeak {

  using namespace std;

  double Kernel::f (const double& x) const {
    return 3 * (1 - pow(x, 2)) / 4; // Epanechnikov
  }

  Kernel::Kernel (const UShort& bandwidthArg, const double& sumArg) : vector<double>(2 * bandwidthArg + 1) {
    double sum = 0;
    
    // calculate values
    iterator kernelIter = begin();
    for (int i = -bandwidthArg; i <= bandwidthArg; ++i) {
      assert(kernelIter != end());
      *kernelIter = f((double)i / (double)bandwidthArg);
      sum += *kernelIter;
      ++kernelIter;
    }
    
    // scale to desired sum
    const double scale = sumArg / sum;
    kernelIter = begin();
    while (kernelIter != end()) {
      *kernelIter *= scale;
      ++kernelIter;
    }
  }

}
