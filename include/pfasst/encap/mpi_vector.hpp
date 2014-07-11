/*
 * MPI enabled vector.
 */

#ifndef _PFASST_MPI_VECTOR_HPP_
#define _PFASST_MPI_VECTOR_HPP_

#include <mpi.h>

#include "vector.hpp"
#include "../mpi_communicator.hpp"

using namespace std;

namespace pfasst {
  namespace encap {

    using namespace pfasst::mpi;

    template<typename scalar, typename time>
    class MPIVectorEncapsulation : public VectorEncapsulation<scalar,time> {

    public:

      MPIVectorEncapsulation(int size) : VectorEncapsulation<scalar,time>(size) { }

      void post(ICommunicator* _comm)
      {
	// auto* comm = dynamic_cast<MPICommunicator*>(_comm);
      }

      void recv(ICommunicator* _comm)
      {
	auto& mpi = *dynamic_cast<MPICommunicator*>(_comm);

      }

      void send(ICommunicator* _comm, int tag, bool blocking)
      {
	auto& mpi = *dynamic_cast<MPICommunicator*>(_comm);

	if (blocking) {
	  int err = MPI_Send(this->data(), sizeof(scalar)*this->size(), MPI_CHAR,
			     (mpi.rank()+1) % mpi.size(), tag, mpi.comm);
	  if (err != MPI_SUCCESS)
	    throw MPIError();
	} else {
	}
      }

    };

    template<typename scalar, typename time>
    class MPIVectorFactory : public EncapFactory<scalar,time> {
      int size;
    public:
      int dofs() { return size; }
      MPIVectorFactory(const int size) : size(size) { }
      Encapsulation<scalar,time>* create(const EncapType) {
	return new MPIVectorEncapsulation<scalar,time>(size);
      }
    };



  }
}

#endif
