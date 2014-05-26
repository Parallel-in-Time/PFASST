/*
 * Host based encapsulated base sweeper.
 */

#ifndef _PFASST_ENCAPSULATED_HPP_
#define _PFASST_ENCAPSULATED_HPP_

#include <vector>

#include "pfasst-interfaces.hpp"
#include "pfasst-quadrature.hpp"

using namespace std;

namespace pfasst {
  namespace encap {

    typedef enum EncapType { solution, function } EncapType;

    //
    // encapsulation
    //

    class Encapsulation {
    public:
      virtual ~Encapsulation() { }

      // required for interp/restrict helpers
      virtual void interpolate(const Encapsulation *) {
	throw NotImplementedYet("mlsdc/pfasst");
      }
      virtual void restrict(const Encapsulation *) {
	throw NotImplementedYet("mlsdc/pfasst");
      }

      // required for time-parallel communications
      virtual void send() {
	throw NotImplementedYet("pfasst");
      }
      virtual void recv() {
	throw NotImplementedYet("pfasst");
      }

      // required for host based encap helpers
      virtual void setval(double) {
	throw NotImplementedYet("encap");
      }
      virtual void copy(const Encapsulation *) {
	throw NotImplementedYet("encap");
      }
      virtual void saxpy(double a, const Encapsulation *) {
	throw NotImplementedYet("encap");
      }
      // virtual void mat_apply(Encapsulation dst[], double a, matrix m,
      // 			   const Encapsulation src[]) {
      //   throw NotImplementedYet("encap");
      // }
    };

    class EncapsulationFactory {
    public:
      virtual Encapsulation* create(const EncapType) = 0;
    };

    template<typename T>
    class EncapsulatedSweeperMixin : public ISweeper {
      vector<T>                        nodes;
      shared_ptr<EncapsulationFactory> factory;

    public:
      void set_nodes(vector<T> nodes) {
	this->nodes = nodes;
      }
      const vector<T> get_nodes() const {
	return nodes;
      }

      void set_factory(EncapsulationFactory* factory) {
	this->factory = shared_ptr<EncapsulationFactory>(factory);
      }
      EncapsulationFactory* get_factory() const {
	return factory.get();
      }

      virtual void set_q0(const Encapsulation* q0) {
	throw NotImplementedYet("sweeper");
      }
      virtual Encapsulation* get_qend() {
	throw NotImplementedYet("sweeper");
	return NULL;
      }

      virtual void advance() {
	throw NotImplementedYet("sweeper");
      }
    };

    template<class T>
    class PolyInterpMixin : public T {
      virtual void interpolate(const ISweeper*) { }
      virtual void restrict(const ISweeper*) { }
    };

  }
}

#endif
