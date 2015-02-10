#ifndef _PFASST_MPI_VECTOR_HPP_
#define _PFASST_MPI_VECTOR_HPP_

#include <mpi.h>

#include "vector.hpp"
#include "../mpi_communicator.hpp"


namespace pfasst
{
  namespace encap
  {
    using namespace pfasst::mpi;

    /**
     * MPI enabled vector.
     */
    template<typename scalar, typename time = time_precision>
    class MPIVectorEncapsulation
      : public VectorEncapsulation<scalar, time>
    {
        //! @{
        MPI_Request recv_request = MPI_REQUEST_NULL;
        MPI_Request send_request = MPI_REQUEST_NULL;
        //! @}

        //! @{
        inline MPICommunicator& as_mpi(ICommunicator* comm)
        {
          auto mpi = dynamic_cast<MPICommunicator*>(comm);
          assert(mpi);
          return *mpi;
        }
        //! @}

      public:
        //! @{
        MPIVectorEncapsulation(const size_t size);
        virtual ~MPIVectorEncapsulation() = default;
        //! @}

        //! @{
        virtual void post(ICommunicator* comm, int tag) override;
        virtual void recv(ICommunicator* comm, int tag, bool blocking) override;
        virtual void send(ICommunicator* comm, int tag, bool blocking) override;
        virtual void broadcast(ICommunicator* comm) override;
        //! @}
    };


    template<typename scalar, typename time = time_precision>
    class MPIVectorFactory
      : public VectorFactory<scalar, time>
    {
      public:
        MPIVectorFactory(const size_t size);
        virtual shared_ptr<Encapsulation<time>> create(const EncapType) override;
    };

  }  // ::pfasst::encap
}  // ::pfasst

#include "mpi_vector_impl.hpp"

#endif
