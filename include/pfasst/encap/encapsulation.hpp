#ifndef _PFASST_ENCAPSULATED_HPP_
#define _PFASST_ENCAPSULATED_HPP_

#include <memory>
#include <vector>
using namespace std;

#include <Eigen/Dense>

template<typename scalar>
using Matrix = Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

#include "pfasst/interfaces.hpp"


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
        virtual ~Encapsulation();
        //! @}

        //! @{
        // required for time-parallel communications
        virtual void post(ICommunicator* comm, int tag);
        virtual void send(ICommunicator* comm, int tag, bool blocking);
        virtual void recv(ICommunicator* comm, int tag, bool blocking);
        virtual void broadcast(ICommunicator* comm);
        //! @}

        //! @{
        // required for host based encap helpers
        virtual void zero();
        virtual void copy(shared_ptr<const Encapsulation<time>>);
        //! @}

        virtual time norm0() const;

        //! @{
        /**
         * provides basic mathematical operation \\(y+=ax\\).
         *
         * This is the main mathematical operation applied by PFASST on the data structures.
         * Here, \\(a\\) is a time point and \\(x\\) another data structure (usually of the
         * same type).
         */
        virtual void saxpy(time a, shared_ptr<const Encapsulation<time>> x);

        /**
         * defines matrix-vector multiplication for this data type.
         */
        virtual void mat_apply(vector<shared_ptr<Encapsulation<time>>> dst,
                               time a, Matrix<time> mat,
                               vector<shared_ptr<Encapsulation<time>>> src,
                               bool zero = true);
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

#include "pfasst/encap/encapsulation_impl.hpp"

#endif
