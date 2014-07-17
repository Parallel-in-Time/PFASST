/*
 * Host based encapsulated base sweeper.
 */

#ifndef _PFASST_ENCAPSULATED_HPP_
#define _PFASST_ENCAPSULATED_HPP_

#include <vector>
#include <memory>

#include "../interfaces.hpp"
#include "../quadrature.hpp"

using namespace std;

namespace pfasst
{
  namespace encap
  {

    typedef enum EncapType { solution, function } EncapType;

    /**
     * basic encapsulation.
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
        virtual ~Encapsulation() { }
        //! @}

        //! @{
        // required for time-parallel communications
        virtual void send()
        {
          throw NotImplementedYet("pfasst");
        }
        virtual void recv()
        {
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
        virtual void saxpy(time /*a*/, shared_ptr<const Encapsulation<time>> /*x*/)
        {
          throw NotImplementedYet("encap");
        }
        /**
         * defines matrix-vector multiplication for this data type.
         */
        virtual void mat_apply(vector<shared_ptr<Encapsulation<time>>> /*dst*/, time /*a*/, matrix<time> /*m*/,
			       vector<shared_ptr<Encapsulation<time>>> /*src*/, bool zero = true)
        {
	  (void) zero;
          throw NotImplementedYet("encap");
        }
        //! @}
    };

    template<typename time = time_precision>
    class EncapFactory
    {
      public:
        virtual shared_ptr<Encapsulation<time>> create(const EncapType) = 0;
    };

  }
}

#endif
