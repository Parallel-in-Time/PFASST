/*
 * Interfaces for SDC/MLSDC/PFASST algorithms.
 */

#ifndef _PFASST_MPI_COMMUNICATOR_HPP_
#define _PFASST_MPI_COMMUNICATOR_HPP_

#include <exception>

#include <mpi.h>

#include "interfaces.hpp"

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

        virtual ~MPICommunicator()
        {}
        //! @}

        //! @{
        void set_comm(MPI_Comm comm)
        {
          this->comm = comm;
          MPI_Comm_size(this->comm, &(this->_size));
          MPI_Comm_rank(this->comm, &(this->_rank));
        }

        int size() { return this->_size; }
        int rank() { return this->_rank; }
        //! @
    };

  }  // ::pfasst::mpi
}  // ::pfasst

#endif
