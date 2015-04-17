/**
 * @file pfasst/encap/encapsulation.hpp
 * @since v0.1.0
 */
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
  /**
   * Encapsulations (short _encaps_) are the central data type for all PFASST++ algorithms.
   * Encaps are both representing the unknown variable as well as the result of evaluating the
   * right hand side of the problem equation(s).
   */
  namespace encap
  {
    //! @todo document EncapType; especially its intention and expected usage
    typedef enum EncapType { solution, function } EncapType;

    /**
     * Data/solution encapsulation.
     *
     * An Encapsulation provides basic functionality of the user data used by PFASST such as
     * mathematical operation _axpy_ \\( y=ax+y \\) and packing/unpacking for message passing.
     *
     * @tparam time time precision; defaults to pfasst::time_precision
     */
    template<typename time = time_precision>
    class Encapsulation
    {
      public:
        //! @{
        virtual ~Encapsulation();
        //! @}

        //! @{
        // required for host based encap helpers
        /**
         * Zeroes out all values of this data structure.
         */
        virtual void zero();

        /**
         * Copies values from @p other into this data structure.
         *
         * @param[in] other other data structure to copy data from
         */
        virtual void copy(shared_ptr<const Encapsulation<time>> other);
        //! @}

        //! @{
        /**
         * Computes the \\( 0 \\)-norm of the data structure's values.
         *
         * @returns \\( 0 \\)-norm of this data structure
         */
        virtual time norm0() const;
        //! @}

        //! @{
        /**
         * Provides basic mathematical operation \\( y+=ax \\).
         *
         * This is the main mathematical operation applied by PFASST on the data structures.
         * Here, \\( a \\) is a time point and \\( x \\) another data structure (usually of the
         * same type) and \\( y \\) is this data structure.
         *
         * @param[in] a time point to multiply
         * @param[in] x another data structure to scale-add onto this one
         */
        virtual void saxpy(time a, shared_ptr<const Encapsulation<time>> x);

        /**
         * Defines matrix-vector multiplication for this data type.
         *
         * This implements the matrix-vector multiplication of the form
         * \\( \\vec{y}=a M \\vec{x} \\).
         * Here, \\( a \\) is a time point to scale the matrix-vector multiplication with, \\( M \\)
         * is a matrix and \\( x \\) another data structure usually of the same type as this one.
         *
         * If @p zero is `true`, then @p dst is zeroed out before applying the matrix-vector
         * multiplication.
         * Elsewise the result of the matrix-vector multiplication is added onto @p dst.
         *
         * @param[in,out] dst  data structure to store the result in
         * @param[in]     a
         * @param[in]     mat
         * @param[in]     src
         * @param[in]     zero
         */
        virtual void mat_apply(vector<shared_ptr<Encapsulation<time>>> dst,
                               time a, Matrix<time> mat,
                               vector<shared_ptr<Encapsulation<time>>> src,
                               bool zero = true);
        //! @}

        //! @{
        // required for time-parallel communications
        /**
         * @todo Write documentation for this member function. How is it different to post/recv?
         *
         * @param[in] comm     communicator managing the processes to post ??? to/from ???
         * @param[in] tag      tag to distinguish overlapping communication
         */
        virtual void post(ICommunicator* comm, int tag);

        /**
         * Send values stored in this data structure.
         *
         * @param[in] comm     communicator managing the processes to send to
         * @param[in] tag      tag to distinguish overlapping communication
         * @param[in] blocking whether to use blocking or non-blocking send
         */
        virtual void send(ICommunicator* comm, int tag, bool blocking);

        /**
         * Receive values to store in this data structure.
         *
         * @param[in] comm     communicator managing the processes to receive from
         * @param[in] tag      tag to distinguish overlapping communication
         * @param[in] blocking whether to use blocking or non-blocking receive
         */
        virtual void recv(ICommunicator* comm, int tag, bool blocking);

        /**
         * Broadcasting this data structure to all processes in @p comm.
         *
         * @param[in] comm communicator managing the processes to send this data structure to
         */
        virtual void broadcast(ICommunicator* comm);
        //! @}
    };


    /**
     * Abstract interface of factory for creating Encapsulation objects.
     *
     * This factory is intendet to be instantiated once to create multiple Encapsulation objects
     * of the same type and with the same parameters later on through calls to
     * EncapFactory::create().
     *
     * Implementations may add values and parameters to the constructor of the factory storing these
     * values as an instance member and accessing them within their implementations of create().
     * For an example see pfasst::encap::VectorFactory.
     *
     * @tparam time time precision; defaults to pfasst::time_precision
     */
    template<typename time = time_precision>
    class EncapFactory
    {
      public:
        /**
         * Actual method to create Encapsulation object of specific type.
         *
         * @param[in] type encapsulation type of the requested Encapsulation object
         */
        virtual shared_ptr<Encapsulation<time>> create(const EncapType type) = 0;
    };
  }  // ::pfasst::encap
} // ::pfasst

#include "pfasst/encap/encapsulation_impl.hpp"

#endif
