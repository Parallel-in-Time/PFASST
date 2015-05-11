#ifndef _PFASST_MPI_COMMUNICATOR_HPP_
#define _PFASST_MPI_COMMUNICATOR_HPP_

#include <stdexcept>
#include <vector>
using namespace std;

#include <mpi.h>

#include "pfasst/interfaces.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  namespace mpi
  {
    class MPIError
      : public runtime_error
    {
      public:
        MPIError(const string& msg="");
        virtual const char* what() const throw();
        static MPIError from_code(const int err_code);
    };


    // forward declare for MPICommunicator
    class MPIStatus;


    class MPICommunicator
      : public ICommunicator
    {
        //! @{
        int _rank;
        int _size;
        string _name;
        //! @}

      public:
        //! @{
        MPI_Comm comm;
        //! @}

        //! @{
        MPICommunicator();
        MPICommunicator(MPI_Comm comm);
        //! @}

        //! @{
        virtual void set_comm(MPI_Comm comm);
        virtual int size();
        virtual int rank();
        virtual string name();
        //! @}
    };


    class MPIStatus
      : public IStatus
    {
      protected:
        vector<bool> converged;
        MPICommunicator* mpi;

      public:
        virtual void set_comm(ICommunicator* comm);
        virtual void clear() override;
        virtual void set_converged(bool converged) override;
        virtual bool get_converged(int rank) override;
        virtual void post();
        virtual void send();
        virtual void recv();
    };
  }  // ::pfasst::mpi
}  // ::pfasst


inline MAKE_LOGGABLE(MPI_Status, mpi_status, os);

#include "pfasst/mpi_communicator_impl.hpp"

#endif
