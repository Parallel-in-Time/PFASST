/*
 * Host based encapsulated base sweeper.
 */

#ifndef _PFASST_ENCAPSULATED_HPP_
#define _PFASST_ENCAPSULATED_HPP_

#include <vector>

#include "../interfaces.hpp"
#include "../quadrature.hpp"

using namespace std;

namespace pfasst
{
  namespace encap
  {

    typedef enum EncapType { solution, function } EncapType;

    //
    // encapsulation
    //
    template<typename scalar, typename time>
    class Encapsulation
    {
      public:
        virtual ~Encapsulation() { }

        // required for time-parallel communications
        virtual void send()
        {
          throw NotImplementedYet("pfasst");
        }
        virtual void recv()
        {
          throw NotImplementedYet("pfasst");
        }

        // required for host based encap helpers
        virtual void setval(scalar)
        {
          throw NotImplementedYet("encap");
        }
        virtual void copy(const Encapsulation<scalar, time>*)
        {
          throw NotImplementedYet("encap");
        }
        virtual void saxpy(time a, const Encapsulation<scalar, time>*)
        {
          throw NotImplementedYet("encap");
        }
        virtual void mat_apply(vector<Encapsulation<scalar, time>*> dst, time a, matrix<time> m,
                               vector<Encapsulation<scalar, time>*> src, bool zero = true)
        {
          throw NotImplementedYet("encap");
        }
    };

    template<typename scalar, typename time>
    class EncapFactory
    {
      public:
        virtual Encapsulation<scalar, time>* create(const EncapType) = 0;
    };

  }
}

#endif
