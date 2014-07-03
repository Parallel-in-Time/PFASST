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

    /**
     * @brief basic encapsulation
     * @details An Encapsulation provides basic functionality of the user data used by PFASST such
     *     as mathematical operation @em axpy @f$y=ax+y@f$ and packing/unpacking for message 
     *     passing.
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
        virtual void copy(const Encapsulation<time>*)
        {
          throw NotImplementedYet("encap");
        }
        //! @}

        //! @{
        /**
         * @brief provides basic mathematical operation @f$y+=ax@f$
         * @details This is the main mathematical operation applied by PFASST on the data 
         *     structures.
         *     Here, @f$a@f$ is a time point and @f$x@f$ another data structure (usually of the 
         *     same type).
         */
        virtual void saxpy(time a, const Encapsulation<time>*)
        {
          throw NotImplementedYet("encap");
        }
        /**
         * @brief defines matrix-vector multiplication for this data type
         */
        virtual void mat_apply(vector<Encapsulation<time>*> dst, time a, matrix<time> m,
                               vector<Encapsulation<time>*> src, bool zero = true)
        {
          throw NotImplementedYet("encap");
        }
        //! @}
    };

    template<typename time = time_precision>
    class EncapFactory
    {
      public:
        virtual Encapsulation<time>* create(const EncapType) = 0;
    };

  }
}

#endif
