#include "pfasst/comm/communicator.hpp"

#include <exception>
#include <stdexcept>
using namespace std;

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  namespace comm
  {
    size_t Communicator::get_size() const
    {
      return 0;
    }

    size_t Communicator::get_rank() const
    {
      return 0;
    }

    size_t Communicator::get_root() const
    {
      return 0;
    }

    bool Communicator::is_first() const
    {
      return false;
    }

    bool Communicator::is_last() const
    {
      return false;
    }
    
    void Communicator::abort(const int& err_code)
    {
      UNUSED(err_code);
      std::abort();
    }

    template<class DataT>
    void Communicator::send(const DataT* const data, const int count, const int dest_rank, const int tag)
    {
      CLOG(ERROR, "COMM") << "blocking send of generic data types not implemented."
                          << " type: " << typeid(data).name();
      throw runtime_error("send for generic data type");
    }

    template<class DataT>
    void Communicator::isend(const DataT* const data, const int count, const int dest_rank, const int tag)
    {
      CLOG(ERROR, "COMM") << "non-blocking send of generic data types not implemented."
                          << " type: " << typeid(data).name();
      throw runtime_error("isend for generic data type");
    }

    template<class DataT>
    void Communicator::recv(DataT* data, const int count, const int dest_rank, const int tag)
    {
      CLOG(ERROR, "COMM") << "blocking receive of generic data types not implemented."
                          << " type: " << typeid(data).name();
      throw runtime_error("recv for generic data type");
    }

    template<class DataT>
    void Communicator::irecv(DataT* data, const int count, const int src_rank, const int tag)
    {
      CLOG(ERROR, "COMM") << "non-blocking receive of generic data types not implemented."
                          << " type: " << typeid(data).name();
      throw runtime_error("irecv for generic data type");
    }

    template<class DataT>
    void Communicator::bcast(DataT* data, const int count, const int root_rank)
    {
      CLOG(ERROR, "COMM") << "braodcast of generic data types not implemented."
                          << " type: " << typeid(data).name();
      throw runtime_error("bcast for generic data type");
    }
  }
}
