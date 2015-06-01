#ifndef _PFASST_MPI_COMMUNICATOR_HPP_
#define _PFASST_MPI_COMMUNICATOR_HPP_

#include <stdexcept>
#include <vector>
using namespace std;

#include <mpi.h>

#include "pfasst/interfaces.hpp"
#include "pfasst/logging.hpp"


/**
 * creates and initializes a new _empty_ `MPI_Status` object
 *
 * An _empty_ `MPI_Status` is defined as an `MPI_Status` object with `MPI_ERROR` as `MPI_SUCCESS`,
 * `MPI_SOURCE` as `MPI_ANY_SOURCE` and `MPI_TAG` as `MPI_ANY_TAG`
 * (cf. MPI Standard v3; section 3.7.3).
 *
 * rationale: some MPI implementations don't initialize the members of `MPI_Status` correctly
 *
 * @returns _empty_ `MPI_Status` object
 *
 * @since v0.5.0
 *
 * @ingroup Utilities
 */
inline static MPI_Status MPI_Status_factory()
{
  MPI_Status stat;
  stat.MPI_ERROR = MPI_SUCCESS;
  stat.MPI_SOURCE = MPI_ANY_SOURCE;
  stat.MPI_TAG = MPI_ANY_TAG;
  return stat;
}


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

    /**
     * checks MPI error code
     *
     * In case @p err_code is not `MPI_SUCCESS` this throws MPIError with the error code looked up
     * to a descriptive string as defined by the MPI implementation.
     *
     * @param[in] err_code MPI error code as returned from an MPI call
     *
     * @since v0.5.0
     *
     * @ingroup Utilities
     */
    inline static void check_mpi_error(const int err_code)
    {
      if (err_code != MPI_SUCCESS) {
        throw MPIError::from_code(err_code);
      }
    }


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
        virtual void post(int tag) override;
        virtual void send(int tag) override;
        virtual void recv(int tag) override;
    };
  }  // ::pfasst::mpi
}  // ::pfasst


inline MAKE_LOGGABLE(MPI_Status, mpi_status, os);

#include "pfasst/mpi_communicator_impl.hpp"

#endif
