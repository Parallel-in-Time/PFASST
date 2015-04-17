#include "pfasst/mpi_communicator.hpp"

#include "pfasst/logging.hpp"


namespace pfasst
{
  namespace mpi
  {
    MPIError::MPIError(const string& msg)
      : runtime_error(msg)
    {}

    const char* MPIError::what() const throw()
    {
      return (string("mpi error: ") + string(runtime_error::what())).c_str();
    }


    MPICommunicator::MPICommunicator()
    {}

    MPICommunicator::MPICommunicator(MPI_Comm comm)
    {
      set_comm(comm);
    }

    void MPICommunicator::set_comm(MPI_Comm comm)
    {
      this->comm = comm;
      MPI_Comm_size(this->comm, &(this->_size));
      MPI_Comm_rank(this->comm, &(this->_rank));

      shared_ptr<MPIStatus> status = make_shared<MPIStatus>();
      this->status = status;
      this->status->set_comm(this);
    }

    int MPICommunicator::size()
    {
      return this->_size;
    }

    int MPICommunicator::rank()
    {
      return this->_rank;
    }


    void MPIStatus::set_comm(ICommunicator* comm)
    {
      this->comm = comm;
      this->converged.resize(comm->size());

      this->mpi = dynamic_cast<MPICommunicator*>(comm); assert(this->mpi);
    }

    void MPIStatus::clear()
    {
      std::fill(converged.begin(), converged.end(), false);
    }

    void MPIStatus::set_converged(bool converged)
    {
      LOG(DEBUG) << "mpi rank " << this->comm->rank() << " set converged to " << converged;
      this->converged.at(this->comm->rank()) = converged;
    }

    bool MPIStatus::get_converged(int rank)
    {
      return this->converged.at(rank);
    }

    void MPIStatus::post()
    {
      // noop: send/recv for status info is blocking
    }

    void MPIStatus::send()
    {
      // don't send forward if: single processor run, or we're the last processor
      if (mpi->size() == 1) { return; }
      if (mpi->rank() == mpi->size() - 1) { return; }

      int iconverged = converged.at(mpi->rank()) ? 1 : 0;

      LOG(DEBUG) << "mpi rank " << this->comm->rank() << " status send " << iconverged;

      int err = MPI_Send(&iconverged, sizeof(int), MPI_INT,
                         (mpi->rank() + 1) % mpi->size(), 1, mpi->comm);

      if (err != MPI_SUCCESS) {
        throw MPIError();
      }
    }

    void MPIStatus::recv()
    {
      // don't recv if: single processor run, or we're the first processor
      if (mpi->size() == 1) { return; }
      if (mpi->rank() == 0) { return; }

      if (get_converged(mpi->rank()-1)) {
        LOG(DEBUG) << "mpi rank " << this->comm->rank() << " skipping status recv";
        return;
      }

      MPI_Status stat;
      int iconverged;
      int err = MPI_Recv(&iconverged, sizeof(iconverged), MPI_INT,
                         (mpi->rank() - 1) % mpi->size(), 1, mpi->comm, &stat);

      if (err != MPI_SUCCESS) {
        throw MPIError();
      }

      converged.at(mpi->rank()-1) = iconverged == 1 ? true : false;

      LOG(DEBUG) << "mpi rank " << this->comm->rank() << " status recv " << iconverged;
    }
  }  // ::pfasst::mpi
}  // ::pfasst
