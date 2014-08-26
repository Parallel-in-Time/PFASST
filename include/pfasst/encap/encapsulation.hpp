/*
 * Host based encapsulated base sweeper.
 */

#ifndef _PFASST_ENCAPSULATED_HPP_
#define _PFASST_ENCAPSULATED_HPP_

#include <vector>
#include <memory>

#include "../globals.hpp"
#include "../interfaces.hpp"
#include "../quadrature.hpp"

using namespace std;

namespace pfasst
{
  namespace encap
  {

    typedef enum EncapType { solution, function } EncapType;

    /**
     * Data/solution encapsulation.
     *
     * An Encapsulation provides basic functionality of the user data used by PFASST such as
     * mathematical operation @em axpy \\(y=ax+y\\) and packing/unpacking for message passing.
     * @tparam time time precision
     *     defaults to pfasst::time_precision
     */
    template<typename time = time_precision>
    class Encapsulation
    {
      public:
        //! @{
        virtual ~Encapsulation()
        {}
        //! @}

        //! @{
        // required for time-parallel communications
        virtual void post(ICommunicator* comm, int tag)
        {
          UNUSED(comm); UNUSED(tag);
        }

        virtual void send(ICommunicator* comm, int tag, bool blocking)
        {
          UNUSED(comm); UNUSED(tag); UNUSED(blocking);
          throw NotImplementedYet("pfasst");
        }

        virtual void recv(ICommunicator* comm, int tag, bool blocking)
        {
          UNUSED(comm); UNUSED(tag); UNUSED(blocking);
          throw NotImplementedYet("pfasst");
        }

        virtual void broadcast(ICommunicator* comm)
        {
          UNUSED(comm);
          throw NotImplementedYet("pfasst");
        }
        //! @}

        //! @{
        // required for host based encap helpers
        virtual void zero()
        {
          throw NotImplementedYet("encap");
        }
        virtual void copy(shared_ptr<const Encapsulation<time>>)
        {
          throw NotImplementedYet("encap");
        }
        //! @}

        //! @{
        /**
         * provides basic mathematical operation \\(y+=ax\\).
         *
         * This is the main mathematical operation applied by PFASST on the data structures.
         * Here, \\(a\\) is a time point and \\(x\\) another data structure (usually of the
         * same type).
         */
        virtual void saxpy(time a, shared_ptr<const Encapsulation<time>> x)
        {
          UNUSED(a); UNUSED(x);
          throw NotImplementedYet("encap");
        }

        /**
         * defines matrix-vector multiplication for this data type.
         */
        virtual void mat_apply(vector<shared_ptr<Encapsulation<time>>> dst, time a, MatrixT<time> mat,
                               vector<shared_ptr<Encapsulation<time>>> src, bool zero = true)
        {
          size_t ndst = dst.size();
          size_t nsrc = src.size();

          if (zero) {
            for (auto elem : dst) { elem->zero(); }
          }

          for (size_t n = 0; n < ndst; n++) {
            for (size_t m = 0; m < nsrc; m++) {
              auto s = mat(n, m);
              if (s != 0.0) {
                dst[n]->saxpy(a*s, src[m]);
              }
            }
          }
        }
        //! @}
    };

    template<typename time = time_precision>
    class EncapFactory
    {
      public:
        virtual shared_ptr<Encapsulation<time>> create(const EncapType) = 0;
    };

  }  // ::pfasst::encap
} // ::pfasst

#endif
