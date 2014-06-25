/*
 * Host based encapsulated base sweeper.
 */

#ifndef _PFASST_ENCAPSULATED_HPP_
#define _PFASST_ENCAPSULATED_HPP_

#include <vector>

#include "../config.hpp"
#include "../interfaces.hpp"
#include "../quadrature.hpp"

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

      // required for time-parallel communications
      virtual void send() {
	throw NotImplementedYet("pfasst");
      }
      virtual void recv() {
	throw NotImplementedYet("pfasst");
      }

      // required for host based encap helpers
      virtual void zero() {
	throw NotImplementedYet("encap");
      }
      virtual void copy(const Encapsulation *) {
	throw NotImplementedYet("encap");
      }
      virtual void saxpy(time a, const Encapsulation *) {
	throw NotImplementedYet("encap");
      }
      virtual void mat_apply(vector<Encapsulation*> dst, time a, matrix<time> m,
			     vector<Encapsulation*> src, bool zero=true) {
        throw NotImplementedYet("encap");
      }
    };

    class EncapFactory {
    public:
      virtual Encapsulation* create(const EncapType) = 0;
    };

  }
}

#endif
