#ifndef _PFASST__COMM__COMMUNICATOR_HPP_
#define _PFASST__COMM__COMMUNICATOR_HPP_

#include <memory>
#include <type_traits>
using namespace std;


namespace pfasst
{
  namespace comm
  {
    class Communicator
      : public enable_shared_from_this<Communicator>
    {
      public:
        Communicator() = default;
        Communicator(const Communicator& other) = default;
        Communicator(Communicator&& other) = default;
        virtual ~Communicator() = default;
        Communicator& operator=(const Communicator& other) = default;
        Communicator& operator=(Communicator&& other) = default;

        virtual size_t get_size() const;
        virtual size_t get_rank() const;
        virtual size_t get_root() const;

        virtual bool is_first() const;
        virtual bool is_last() const;
        
        virtual void abort(const int& err_code);

        template<class DataT>
        void send(const DataT* const data, const int count, const int dest_rank, const int tag);

        template<class DataT>
        void isend(const DataT* const data, const int count, const int dest_rank, const int tag);

        template<class DataT>
        void recv(DataT* data, const int count, const int dest_rank, const int tag);

        template<class DataT>
        void irecv(DataT* data, const int count, const int src_rank, const int tag);

        template<class DataT>
        void bcast(DataT* data, const int count, const int root_rank);
    };
  }  // ::pfasst::comm
}  // ::pfasst

#include "pfasst/comm/communicator_impl.hpp"

#endif  // _PFASST__COMM__COMMUNICATOR_HPP_
