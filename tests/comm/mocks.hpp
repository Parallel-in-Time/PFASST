#include "fixtures/test_helpers.hpp"

#include <pfasst/logging.hpp>
#include <pfasst/controller/status.hpp>
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

    template<typename DataT>
    void send(const DataT* const data, const int count, const int dest_rank, const int tag)
    {
      LOG(WARNING) << "TODO: mock blocking send";
    }
    template<typename DataT>
    void send_status(const pfasst::StatusDetail<DataT>* const data, const int count, const int dest_rank, const int tag)
    {
      LOG(WARNING) << "TODO: mock blocking send";
    }

    template<typename DataT>
    void isend(const DataT* const data, const int count, const int dest_rank, const int tag)
    {
      LOG(WARNING) << "TODO: mock non-blocking send";
    }
    template<typename DataT>
    void isend_status(const pfasst::StatusDetail<DataT>* const data, const int count, const int dest_rank, const int tag)
    {
      LOG(WARNING) << "TODO: mock non-blocking send";
    }

    template<typename DataT>
    void recv(DataT* data, const int count, const int dest_rank, const int tag)
    {
      LOG(WARNING) << "TODO: mock blocking receive";
    }
    template<typename DataT>
    void recv_status(pfasst::StatusDetail<DataT>* data, const int count, const int dest_rank, const int tag)
    {
      LOG(WARNING) << "TODO: mock blocking receive";
    }

    template<typename DataT>
    void irecv(DataT* data, const int count, const int src_rank, const int tag)
    {
      LOG(WARNING) << "TODO: mock non-blocking receive";
    }
    template<typename DataT>
    void irecv_status(pfasst::StatusDetail<DataT>* data, const int count, const int src_rank, const int tag)
    {
      LOG(WARNING) << "TODO: mock non-blocking receive";
    }

    template<typename DataT>
    void bcast(DataT* data, const int count, const int root_rank)
    {
      LOG(WARNING) << "TODO: mock broadcast";
    }
    template<typename DataT>
    void bcast_status(pfasst::StatusDetail<DataT>* data, const int count, const int root_rank)
    {
      LOG(WARNING) << "TODO: mock broadcast";
    }
};
