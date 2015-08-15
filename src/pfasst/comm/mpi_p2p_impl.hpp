#include "pfasst/comm/mpi_p2p.hpp"

#include <string>
using namespace std;

#include "pfasst/logging.hpp"


namespace pfasst
{
  namespace comm
  {
    string error_from_code(const int err_code)
    {
      char err_str[MPI_MAX_ERROR_STRING];
      int err_len = 0;
      int err = MPI_Error_string(err_code, err_str, &err_len);
      check_mpi_error(err);
      return string(err_str, err_len) + " (code=" + to_string(err_code) + ")";
    }


    MPI_Status MPI_Status_factory()
    {
      MPI_Status stat;
      stat.MPI_ERROR = MPI_SUCCESS;
      stat.MPI_SOURCE = MPI_ANY_SOURCE;
      stat.MPI_TAG = MPI_ANY_TAG;
      return stat;
    }

    void check_mpi_error(const int err_code)
    {
      if (err_code != MPI_SUCCESS) {
        string err_msg = error_from_code(err_code);
        CLOG(ERROR, "COMM_P2P") << "MPI encountered an error: " << err_msg;
        throw runtime_error("MPI encountered an error: " + err_msg);
      }
    }


    MpiP2P::MpiP2P(MPI_Comm comm)
      :   _comm(comm)
        , _stati(0)
    {
      log::add_custom_logger("COMM_P2P");

      // get communicator's size and processors rank
      MPI_Comm_size(this->_comm, &(this->_size));
      MPI_Comm_rank(this->_comm, &(this->_rank));

      // get communicator's name (if available)
      int len = -1;
      char buff[MPI_MAX_OBJECT_NAME];
      int err = MPI_Comm_get_name(this->_comm, buff, &len);
      check_mpi_error(err);
      if (len > 0) {
        this->_name = string(buff, len);
      }
    }

    size_t MpiP2P::get_size() const
    {
      assert(this->_size > 0);
      return this->_size;
    }

    size_t MpiP2P::get_rank() const
    {
      assert(this->_rank >= 0);
      return this->_rank;
    }

    string MpiP2P::get_name() const
    {
      return this->_name;
    }

    bool MpiP2P::is_first() const
    {
      return (this->get_rank() == this->get_root());
    }

    bool MpiP2P::is_last() const
    {
      return false;
    }

    template<class DataT>
    void MpiP2P::send(const DataT* const data, const int count, const int dest_rank, const int tag)
    {
      CLOG(ERROR, "COMM_P2P") << "blocking send of generic data types not implemented."
                              << " type: " << typeid(data).name();
      throw runtime_error("blocking send for generic data type");
    }

    template<>
    void MpiP2P::send(const double* const data, const int count, const int dest_rank, const int tag)
    {
      CLOG(DEBUG, "COMM_P2P") << "sending " << count << " double values with tag=" << tag << " to " << dest_rank;
      int err = MPI_Send(data, count, MPI_DOUBLE, dest_rank, tag, this->_comm);
      check_mpi_error(err);
    }

    template<class DataT>
    void MpiP2P::isend(const DataT* const data, const int count, const int dest_rank, const int tag)
    {
      CLOG(ERROR, "COMM_P2P") << "non-blocking send of generic data types not implemented."
                              << " type: " << typeid(data).name();
      throw runtime_error("non-blocking send for generic data type");
    }

    template<>
    void MpiP2P::isend(const double* const data, const int count, const int dest_rank, const int tag)
    {
      auto request_index = make_pair(dest_rank, tag);
      
    }

    template<class DataT>
    void MpiP2P::recv(DataT* data, const int count, const int dest_rank, const int tag)
    {
      CLOG(ERROR, "COMM_P2P") << "blocking receive of generic data types not implemented."
                              << " type: " << typeid(data).name();
      throw runtime_error(" blocking recv for generic data type");
    }

    template<>
    void MpiP2P::recv(double* data, const int count, const int dest_rank, const int tag)
    {
      this->_stati.push_back(MPI_Status_factory());
      CLOG(DEBUG, "COMM_P2P") << "receiving " << count << " double values with tag=" << tag << " to " << dest_rank;
      int err = MPI_Recv(data, count, MPI_DOUBLE, dest_rank, tag, this->_comm, &(this->_stati.back()));
      check_mpi_error(err);
    }

    template<class DataT>
    void MpiP2P::irecv(DataT* data, const int count, const int src_rank, const int tag)
    {
      CLOG(ERROR, "COMM_P2P") << "non-blocking receive of generic data types not implemented."
                              << " type: " << typeid(data).name();
      throw runtime_error("non-blocking receive for generic data type");
    }

    template<>
    void MpiP2P::irecv(double* data, const int count, const int src_rank, const int tag)
    {
      auto request_index = make_pair(src_rank, tag);
    }

    template<class DataT>
    void MpiP2P::bcast(DataT* data, const int count, const int root_rank)
    {
      CLOG(ERROR, "COMM_P2P") << "braodcast of generic data types not implemented."
                              << " type: " << typeid(data).name();
      throw runtime_error("broadcast for generic data type");
    }

    template<>
    void MpiP2P::bcast(double* data, const int count, const int root_rank)
    {
      CLOG(DEBUG, "COMM_P2P") << "braodcasting " << count << " double values from root " << root_rank;
      int err = MPI_Bcast(data, count, MPI_DOUBLE, root_rank, this->_comm);
      check_mpi_error(err);
    }
  }
}
