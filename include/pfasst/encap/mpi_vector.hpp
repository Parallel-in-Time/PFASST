/*
 * MPI enabled vector.
 */

#ifndef _PFASST_MPI_VECTOR_HPP_
#define _PFASST_MPI_VECTOR_HPP_

#include <cassert>

#include <mpi.h>

#include "vector.hpp"
#include "../mpi_communicator.hpp"

using namespace std;

namespace pfasst
{
  namespace encap
  {

    using namespace pfasst::mpi;

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
        MPIVectorEncapsulation(int size)
          : VectorEncapsulation<scalar, time>(size)
        {}

        virtual ~MPIVectorEncapsulation()
        {}
        //! @}

        //! @{
        void post(ICommunicator* comm, int tag) override
        {
          auto& mpi = as_mpi(comm);
          if (mpi.size() == 1) { return; }
          if (mpi.rank() == 0) { return; }

          int err = MPI_Irecv(this->data(), sizeof(scalar) * this->size(), MPI_CHAR,
                              (mpi.rank() - 1) % mpi.size(), tag, mpi.comm, &recv_request);
          if (err != MPI_SUCCESS) {
            throw MPIError();
          }
        }

        virtual void recv(ICommunicator* comm, int tag, bool blocking) override
        {
          auto& mpi = as_mpi(comm);
          if (mpi.size() == 1) { return; }
          if (mpi.rank() == 0) { return; }

          int err;
          if (blocking) {
            MPI_Status stat;
            err = MPI_Recv(this->data(), sizeof(scalar) * this->size(), MPI_CHAR,
                           (mpi.rank() - 1) % mpi.size(), tag, mpi.comm, &stat);
          } else {
            MPI_Status stat;
            err = MPI_Wait(&recv_request, &stat);
          }

          if (err != MPI_SUCCESS) {
            throw MPIError();
          }
        }

        virtual void send(ICommunicator* comm, int tag, bool blocking) override
        {
          auto& mpi = as_mpi(comm);
          if (mpi.size() == 1) { return; }
          if (mpi.rank() == mpi.size() - 1) { return; }

          int err = MPI_SUCCESS;
          if (blocking) {
            err = MPI_Send(this->data(), sizeof(scalar) * this->size(), MPI_CHAR,
                           (mpi.rank() + 1) % mpi.size(), tag, mpi.comm);
          } else {
            MPI_Status stat;
            err = MPI_Wait(&send_request, &stat);
            if (err != MPI_SUCCESS) {
              throw MPIError();
            }

            err = MPI_Isend(this->data(), sizeof(scalar) * this->size(), MPI_CHAR,
                            (mpi.rank() + 1) % mpi.size(), tag, mpi.comm, &send_request);
          }

          if (err != MPI_SUCCESS) {
            throw MPIError();
          }
        }

        virtual void broadcast(ICommunicator* comm) override
        {
          auto& mpi = as_mpi(comm);
          int err = MPI_Bcast(this->data(), sizeof(scalar) * this->size(), MPI_CHAR,
                              comm->size()-1, mpi.comm);

          if (err != MPI_SUCCESS) {
            throw MPIError();
          }
        }
        //! @}
    };

    template<typename scalar, typename time = time_precision>
    class MPIVectorFactory
      : public EncapFactory<time>
    {
        int size;
      public:
        int dofs() { return size; }
        MPIVectorFactory(const int size) : size(size) { }
        shared_ptr<Encapsulation<time>> create(const EncapType)
        {
          return make_shared<MPIVectorEncapsulation<scalar, time>>(size);
        }
    };

  }  // ::pfasst::encap
}  // ::pfasst

#endif
