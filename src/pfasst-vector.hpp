/*
 * Host based encapsulated base sweeper.
 */

#ifndef _PFASST_VECTOR_HPP_
#define _PFASST_VECTOR_HPP_

#include <algorithm>
#include <vector>

#include "pfasst-encapsulated.hpp"

using namespace std;

namespace pfasst {
  namespace encap {

    template<typename scalar>
    struct VectorEncapsulation : public vector<scalar>, public Encapsulation {
      VectorEncapsulation(int size) : vector<scalar>(size) { }
      void setval(scalar v) {
	std::fill(this->begin(), this->end(), v);
      }
      void copy(const Encapsulation* X) {
	const auto* x = dynamic_cast<const VectorEncapsulation*>(X);
	std::copy(x->begin(), x->end(), this->begin());
      }
      void saxpy(double a, const Encapsulation *X) {
	const auto& x = *dynamic_cast<const VectorEncapsulation*>(X);
	auto&       y = *this;
	for (int i=0; i<y.size(); i++)
	  y[i] += a * x[i];
      }
    };

    template<typename T>
    class VectorFactory : public EncapsulationFactory {
      int size;
    public:
      int dofs() { return size; }
      VectorFactory(const int size) : size(size) { }
      Encapsulation* create(const EncapType) {
	return new VectorEncapsulation<T>(size);
      }
    };

  }
}

#endif
