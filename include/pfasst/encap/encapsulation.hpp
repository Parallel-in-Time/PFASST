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
   * Encaps can represent the unknown variable or the result of evaluating the right hand side of
   * the problem equation(s).
   */
  namespace encap
  {
    typedef enum EncapType { solution, function } EncapType;

    /**
     * Data/solution encapsulation.
     *
     * An Encapsulation provides basic mathematical functionality for the user's data.
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
         * Provides basic mathematical operation \\( y = ax + y \\).
         *
         * This is the main mathematical operation applied by PFASST on the data structures.
         * Here, \\( a \\) is a constant and \\( x \\) another data structure (usually of the same
         * type) and \\( y \\) is this data structure.
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
        /**
         * Prepare to receive a solution (MPI_IRecv).
         *
         * @param[in] comm     communicator managing the processes to post ??? to/from ???
         * @param[in] tag      tag to distinguish overlapping communication
         */
        virtual void post(ICommunicator* comm, int tag);

        /**
         * Send solution (MPI_Send).
         *
         * @param[in] comm     communicator managing the processes to send to
         * @param[in] tag      tag to distinguish overlapping communication
         * @param[in] blocking whether to use blocking or non-blocking send
         */
        virtual void send(ICommunicator* comm, int tag, bool blocking);

        /**
         * Receive solution (MPI_Recv).
         *
         * @param[in] comm     communicator managing the processes to receive from
         * @param[in] tag      tag to distinguish overlapping communication
         * @param[in] blocking whether to use blocking or non-blocking receive
         */
        virtual void recv(ICommunicator* comm, int tag, bool blocking);

        /**
         * Broadcast this data structure to all processes in @p comm.
         *
         * @param[in] comm communicator managing the processes to send this data structure to
         */
        virtual void broadcast(ICommunicator* comm);
        //! @}
    };


    /**
     * Abstract interface of factory for creating Encapsulation objects.
     *
     * This factory is intended to be instantiated once to create multiple Encapsulation objects
     * of the same type and with the same parameters later on through calls to
     * EncapFactory::create().
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
