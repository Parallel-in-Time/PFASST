#include "fixtures/test_helpers.hpp"

#include <pfasst/logging.hpp>
#include <pfasst/comm/communicator.hpp>


class CommMock
  : public pfasst::comm::Communicator
{
  public:
    MOCK_CONST_METHOD0(get_size, size_t());
    MOCK_CONST_METHOD0(get_rank, size_t());
    MOCK_CONST_METHOD0(get_root, size_t());

    MOCK_CONST_METHOD0(is_first, bool());
    MOCK_CONST_METHOD0(is_last, bool());

    template<typename T>
    void send(const T* const data, const int count, const int dest_rank, const int tag)
    {
      LOG(WARNING) << "TODO: mock blocking send";
    }

    template<typename T>
    void isend(const T* const data, const int count, const int dest_rank, const int tag)
    {
      LOG(WARNING) << "TODO: mock non-blocking send";
    }

    template<typename T>
    void recv(T* data, const int count, const int src_rank, const int tag)
    {
      LOG(WARNING) << "TODO: mock blocking receive";
    }

    template<typename T>
    void irecv(T* data, const int count, const int src_rank, const int tag)
    {
      LOG(WARNING) << "TODO: mock non-blocking receive";
    }

    template<typename T>
    void bcast(T* data, const int count, const int root_rank)
    {
      LOG(WARNING) << "TODO: mock broadcast";
    }
};
