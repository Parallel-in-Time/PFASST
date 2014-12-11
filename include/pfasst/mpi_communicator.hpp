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


    class MPICommunicator
      : public ICommunicator
    {
        //! @{
        int _rank;
        int _size;
        vector<int> status;
        //! @}

      public:
        //! @{
        MPI_Comm comm;
        MPI_Win  status_window;
        //! @}

        //! @{
        MPICommunicator()
        {}

        MPICommunicator(MPI_Comm comm)
        {
          set_comm(comm);
        }

        virtual ~MPICommunicator()
        {
          MPI_Win_fence(MPI_MODE_NOSUCCEED, this->status_window);
          MPI_Win_free(&this->status_window);
        }
        //! @}

        //! @{
        void set_comm(MPI_Comm comm)
        {
          this->comm = comm;
          MPI_Comm_size(this->comm, &(this->_size));
          MPI_Comm_rank(this->comm, &(this->_rank));

          status.resize(this->_size);
          auto err = MPI_Win_create(status.data(), this->_size * sizeof(int),
                                    sizeof(int), MPI_INFO_NULL, this->comm,
                                    &this->status_window);
          if (err != MPI_SUCCESS) {
            throw MPIError();
          }
          MPI_Win_fence(MPI_MODE_NOPRECEDE, this->status_window);
        }

        virtual void set_status(bool converged) override
        {
          LOG(DEBUG) << "mpi rank " << this->_rank << " set converged to " << converged;
          this->status[this->_rank] = converged ? 1 : 0;
          for (int dst = 0; dst < this->_size; dst++) {
            if (dst == this->_rank) { continue; }
            auto err = MPI_Put(this->status.data()+this->_rank, 1, MPI_INT,
                               dst, this->_rank, 1, MPI_INT, this->status_window);
            if (err != MPI_SUCCESS) {
              throw MPIError();
            }
          }
        }

        virtual void fence_status() override
        {
          auto err = MPI_Win_fence(0, this->status_window);
          if (err != MPI_SUCCESS) {
            throw MPIError();
          }
        }

        int size() { return this->_size; }
        int rank() { return this->_rank; }
        //! @
    };

  }  // ::pfasst::mpi
}  // ::pfasst

#endif
