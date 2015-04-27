#include "pfasst/encap/mpi_vector.hpp"

#include <cassert>
using namespace std;


namespace pfasst
{
  namespace encap
  {
    template<typename scalar, typename time>
    MPIVectorEncapsulation<scalar, time>::MPIVectorEncapsulation(const size_t size)
      : VectorEncapsulation<scalar, time>(size)
    {}

    template<typename scalar, typename time>
    void MPIVectorEncapsulation<scalar, time>::post(ICommunicator* comm, int tag)
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

    template<typename scalar, typename time>
    void MPIVectorEncapsulation<scalar, time>::recv(ICommunicator* comm, int tag, bool blocking)
    {
      auto& mpi = as_mpi(comm);
      if (mpi.size() == 1) { return; }

      int err;
      if (blocking) {
        MPI_Status stat;
        err = MPI_Recv(this->data(), sizeof(scalar) * this->size(), MPI_CHAR,
                       mpi.rank() == 0 ? mpi.size() - 1 : mpi.rank() - 1, tag, mpi.comm, &stat);
      } else {
        MPI_Status stat;
        err = MPI_Wait(&recv_request, &stat);
      }

      if (err != MPI_SUCCESS) {
        throw MPIError();
      }
    }

    template<typename scalar, typename time>
    void MPIVectorEncapsulation<scalar, time>::send(ICommunicator* comm, int tag, bool blocking)
    {
      auto& mpi = as_mpi(comm);
      if (mpi.size() == 1) { return; }

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

    template<typename scalar, typename time>
    void MPIVectorEncapsulation<scalar, time>::broadcast(ICommunicator* comm)
    {
      auto& mpi = as_mpi(comm);
      int err = MPI_Bcast(this->data(), sizeof(scalar) * this->size(), MPI_CHAR,
                          comm->size()-1, mpi.comm);

      if (err != MPI_SUCCESS) {
        throw MPIError();
      }
    }


    template<typename scalar, typename time>
    MPIVectorFactory<scalar, time>::MPIVectorFactory(const size_t size)
      : pfasst::encap::VectorFactory<scalar, time>(size)
    {}

    template<typename scalar, typename time>
    shared_ptr<Encapsulation<time>> MPIVectorFactory<scalar, time>::create(const EncapType)
    {
      return make_shared<MPIVectorEncapsulation<scalar, time>>(this->dofs());
    }
  }  // ::pfasst::encap
}  // ::pfasst
