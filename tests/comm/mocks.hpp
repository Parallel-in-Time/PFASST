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

    MOCK_METHOD4(send, void(const double* const data, const int count, const int dest_rank, const int tag));
    MOCK_METHOD4(send_status, void(const pfasst::StatusDetail<double>* const data, const int count, const int dest_rank, const int tag));

    MOCK_METHOD4(isend, void(const double* const data, const int count, const int dest_rank, const int tag));
    MOCK_METHOD4(isend_status, void(const pfasst::StatusDetail<double>* const data, const int count, const int dest_rank, const int tag));

    MOCK_METHOD4(recv, void(double* data, const int count, const int dest_rank, const int tag));
    MOCK_METHOD4(recv_status, void(pfasst::StatusDetail<double>* data, const int count, const int dest_rank, const int tag));

    MOCK_METHOD4(irecv, void(double* data, const int count, const int src_rank, const int tag));
    MOCK_METHOD4(irecv_status, void(pfasst::StatusDetail<double>* data, const int count, const int src_rank, const int tag));

    MOCK_METHOD3(bcast, void(double* data, const int count, const int root_rank));
};
