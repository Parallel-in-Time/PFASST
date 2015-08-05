#ifndef _PFASST__COMMUNICATOR_HPP_
#define _PFASST__COMMUNICATOR_HPP_

#include <memory>
#include <type_traits>
using namespace std;


namespace pfasst
{
  namespace comm
  {
    class Communicator
    {
      public:
        Communicator() = default;
        Communicator(const Communicator& other) = default;
        Communicator(Communicator&& other) = default;
        ~Communicator() = default;
        Communicator& operator=(const Communicator& other) = default;
        Communicator& operator=(Communicator&& other) = default;

        virtual size_t get_size() const;
        virtual size_t get_rank() const;
        virtual size_t get_root() const;

        virtual bool is_first() const;
        virtual bool is_last() const;

        virtual void send(const int* data, const int count,
                          const int dest_rank, const int tag);
        virtual void send(const double* data, const int count,
                          const int dest_rank, const int tag);

        virtual void isend(const int* data, const int count,
                           const int dest_rank, const int tag);
        virtual void isend(const double* data, const int count,
                           const int dest_rank, const int tag);

        virtual void recv(int* data, const int count,
                          const int dest_rank, const int tag);
        virtual void recv(double* data, const int count,
                          const int dest_rank, const int tag);

        virtual void irecv(int* data, const int count,
                           const int src_rank, const int tag);
        virtual void irecv(double* data, const int count,
                           const int src_rank, const int tag);

        virtual void bcast(int* data, const int count, const int root_rank);
        virtual void bcast(double* data, const int count, const int root_rank);
    };
  }  // ::pfasst::comm
}  // ::pfasst

#include "pfasst/comm/communicator_impl.hpp"

#endif  // _PFASST__COMMUNICATOR_HPP_
