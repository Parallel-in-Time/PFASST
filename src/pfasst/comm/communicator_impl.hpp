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


    void Communicator::send(const double* const data, const int count, const int dest_rank, const int tag)
    {
      UNUSED(data); UNUSED(count); UNUSED(dest_rank); UNUSED(tag);
      throw runtime_error("not implemented: send of double");
    }

    void Communicator::send_status(const StatusDetail<double>* const data, const int count, const int dest_rank, const int tag)
    {
      UNUSED(data); UNUSED(count); UNUSED(dest_rank); UNUSED(tag);
      throw runtime_error("not implemented: send of status details");
    }


    void Communicator::isend(const double* const data, const int count, const int dest_rank, const int tag)
    {
      UNUSED(data); UNUSED(count); UNUSED(dest_rank); UNUSED(tag);
      throw runtime_error("not implemented: isend of double");
    }

    void Communicator::isend_status(const StatusDetail<double>* const data, const int count, const int dest_rank, const int tag)
    {
      UNUSED(data); UNUSED(count); UNUSED(dest_rank); UNUSED(tag);
      throw runtime_error("not implemented: isend of status details");
    }


    void Communicator::recv(double* data, const int count, const int src_rank, const int tag)
    {
      UNUSED(data); UNUSED(count); UNUSED(src_rank); UNUSED(tag);
      throw runtime_error("not implemented: recv of double");
    }

    void Communicator::recv_status(StatusDetail<double>* data, const int count, const int src_rank, const int tag)
    {
      UNUSED(data); UNUSED(count); UNUSED(src_rank); UNUSED(tag);
      throw runtime_error("not implemented: recv of status details");
    }


    void Communicator::irecv(double* data, const int count, const int src_rank, const int tag)
    {
      UNUSED(data); UNUSED(count); UNUSED(src_rank); UNUSED(tag);
      throw runtime_error("not implemented: irecv of double");
    }

    void Communicator::irecv_status(StatusDetail<double>* data, const int count, const int src_rank, const int tag)
    {
      UNUSED(data); UNUSED(count); UNUSED(src_rank); UNUSED(tag);
      throw runtime_error("not implemented: irecv of status details");
    }


    void Communicator::bcast(double* data, const int count, const int root_rank)
    {
      UNUSED(data); UNUSED(count); UNUSED(root_rank);
      throw runtime_error("not implemented: bcast of double");
    }
  }  // ::pfasst::comm
}  // ::pfasst
