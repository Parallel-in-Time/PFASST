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
      int len = 0;
      char buff[MPI_MAX_OBJECT_NAME];
      MPI_Comm_get_name(this->comm, buff, &len);
      if (len == 0) {
        this->_name = string("world");
      } else {
        this->_name = string(buff, len);
      }

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

    string MPICommunicator::name()
    {
      return this->_name;
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
      CLOG(DEBUG, "Controller") << "set converged to " << boolalpha << converged;
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

      int dest_rank = (mpi->rank() + 1) % mpi->size();
      CLOG(DEBUG, "Controller") << "sending status " << iconverged
                                << " to " << dest_rank << " of communicator " << mpi->name();

      int err = MPI_Send(&iconverged, sizeof(int), MPI_INT,
                         dest_rank, 1, mpi->comm);

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
        CLOG(DEBUG, "Controller") << "skipping status recv";
        return;
      }

      MPI_Status stat;
      int iconverged;
      int src_rank = (mpi->rank() - 1) % mpi->size();
      int err = MPI_Recv(&iconverged, sizeof(iconverged), MPI_INT,
                         src_rank, 1, mpi->comm, &stat);

      if (err != MPI_SUCCESS) {
        throw MPIError();
      }

      converged.at(mpi->rank()-1) = iconverged == 1 ? true : false;

      CLOG(DEBUG, "Controller") << "received status " << iconverged
                                << " from rank " << src_rank << " of communicator " << mpi->name();
    }
  }  // ::pfasst::mpi
}  // ::pfasst


MAKE_LOGGABLE(MPI_Status, mpi_status, os)
{
  os << "MPI_Status(source=" << mpi_status.MPI_SOURCE << ", "
                << "tag=" << mpi_status.MPI_TAG << ", "
                << "error=" << mpi_status.MPI_ERROR << ")";
  return os;
}
