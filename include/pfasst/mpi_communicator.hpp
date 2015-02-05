/*
 * Interfaces for SDC/MLSDC/PFASST algorithms.
 */

#ifndef _PFASST_MPI_COMMUNICATOR_HPP_
#define _PFASST_MPI_COMMUNICATOR_HPP_

#include <exception>
#include <vector>

#include <mpi.h>

#include "interfaces.hpp"
#include "logging.hpp"

using namespace std;

namespace pfasst
{
  namespace mpi
  {

    class MPIError
      : public exception
    {
      public:
        const char* what() const throw()
        {
          return "mpi error";
        }
    };

    class MPIStatus;


    class MPICommunicator
      : public ICommunicator
    {
        //! @{
        int _rank;
        int _size;
        //! @}

      public:
        //! @{
        MPI_Comm comm;
        //! @}

        //! @{
        MPICommunicator()
        {}

        MPICommunicator(MPI_Comm comm)
        {
          set_comm(comm);
        }
        //! @}

        //! @{
        void set_comm(MPI_Comm comm)
        {
          this->comm = comm;
          MPI_Comm_size(this->comm, &(this->_size));
          MPI_Comm_rank(this->comm, &(this->_rank));

          shared_ptr<MPIStatus> status = make_shared<MPIStatus>();
          this->status = status;
          this->status->set_comm(this);
        }

        int size() { return this->_size; }
        int rank() { return this->_rank; }
        //! @}
    };


    class MPIStatus
      : public IStatus
    {
        vector<bool> converged;
        MPICommunicator* mpi;

      public:

        virtual void set_comm(ICommunicator* comm)
        {
          this->comm = comm;
          this->converged.resize(comm->size());

          this->mpi = dynamic_cast<MPICommunicator*>(comm); assert(this->mpi);
        }

        virtual void clear() override
        {
          std::fill(converged.begin(), converged.end(), false);
        }

        virtual void set_converged(bool converged) override
        {
          LOG(DEBUG) << "mpi rank " << this->comm->rank() << " set converged to " << converged;
          this->converged.at(this->comm->rank()) = converged;
        }

        virtual bool get_converged(int rank) override
        {
          return this->converged.at(rank);
        }

        virtual void post()
        {
          // noop: send/recv for status info is blocking
        }

        virtual void send()
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

        virtual void recv()
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
    };

  }  // ::pfasst::mpi
}  // ::pfasst

#endif
