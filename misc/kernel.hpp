#ifndef KERNEL_H
#define KERNEL_H

#include <cstdlib>
#include <vector>

namespace unipeak {

  using namespace std;

  class Kernel : public vector<double> {
    protected:
      double f (const double&) const;
    
    public:
      Kernel (const UShort&, const double& = 1); // bandwidth, sum
  };

}

#endif
