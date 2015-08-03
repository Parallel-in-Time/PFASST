#include "pfasst/comm/interface.hpp"

#include "pfasst/exceptions.hpp"


namespace pfasst
{
  namespace comm
  {
    size_t Communicator::get_size() const
    {
      return 3;
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

    void Communicator::send(const int* data, const int count, const int dest_rank, const int tag)
    {
      throw NotImplementedYet("send for generic data type");
    }
    void Communicator::send(const double* data, const int count, const int dest_rank, const int tag)
    {
      throw NotImplementedYet("send for generic data type");
    }

    void Communicator::isend(const int* data, const int count, const int dest_rank, const int tag)
    {
      throw NotImplementedYet("isend for generic data type");
    }
    void Communicator::isend(const double* data, const int count, const int dest_rank, const int tag)
    {
      throw NotImplementedYet("isend for generic data type");
    }

    void Communicator::recv(int* data, const int count, const int dest_rank, const int tag)
    {
      throw NotImplementedYet("recv for generic data type");
    }
    void Communicator::recv(double* data, const int count, const int dest_rank, const int tag)
    {
      throw NotImplementedYet("recv for generic data type");
    }

    void Communicator::irecv(int* data, const int count, const int src_rank, const int tag)
    {
      throw NotImplementedYet("irecv for generic data type");
    }
    void Communicator::irecv(double* data, const int count, const int src_rank, const int tag)
    {
      throw NotImplementedYet("irecv for generic data type");
    }

    void Communicator::bcast(int* data, const int count, const int root_rank)
    {
      throw NotImplementedYet("bcast for generic data type");
    }
    void Communicator::bcast(double* data, const int count, const int root_rank)
    {
      throw NotImplementedYet("bcast for generic data type");
    }
  }
}
