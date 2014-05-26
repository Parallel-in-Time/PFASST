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
      virtual unsigned int nbytes() const {
	return sizeof(scalar) * this->size();
      }
      void setval(scalar v) {
	for (int i=0; i<this->size(); i++)
	  (*this)[i] = v;
      }
      void copy(const Encapsulation* X) {
	const auto* x = dynamic_cast<const VectorEncapsulation*>(X);
	for (int i=0; i<this->size(); i++)
	  (*this)[i] = (*x)[i];
      }
      void saxpy(double a, const Encapsulation *X) {
	const auto* x = dynamic_cast<const VectorEncapsulation*>(X);
	for (int i=0; i<this->size(); i++)
	  (*this)[i] += a * (*x)[i];
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
