/*
 * Interfaces for SDC/MLSDC/PFASST algorithms.
 */

#ifndef _PFASST_MPI_COMMUNICATOR_HPP_
#define _PFASST_MPI_COMMUNICATOR_HPP_

#include <exception>

#include <mpi.h>

using namespace std;

namespace pfasst {
  namespace mpi {

    class MPIError : public exception {
    public:
      const char *what() const throw() {
	return "mpi error";
      }
    };


    class MPICommunicator : public ICommunicator {
      int _rank, _size;
    public:
      MPI_Comm comm;

      MPICommunicator(MPI_Comm comm) { set_comm(comm); }
      MPICommunicator() { }

      void set_comm(MPI_Comm comm)
      {
	this->comm = comm;
	MPI_Comm_size(this->comm, &_size);
	MPI_Comm_rank(this->comm, &_rank);
      }

      int size() { return _size; }
      int rank() { return _rank; }

    };

  }

}

#endif
