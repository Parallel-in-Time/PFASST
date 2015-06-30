#include "fixtures/test_helpers.hpp"

#include <pfasst/comm/interface.hpp>


class CommMock
  : public pfasst::comm::Communicator
{
  public:
    MOCK_METHOD4(send, void(const int* data, const int count,
                            const int dest_rank, const int tag));
    MOCK_METHOD4(send, void(const double* data, const int count,
                            const int dest_rank, const int tag));

    MOCK_METHOD4(isend, void(const int* data, const int count,
                             const int dest_rank, const int tag));
    MOCK_METHOD4(isend, void(const double* data, const int count,
                             const int dest_rank, const int tag));

    MOCK_METHOD4(recv, void(int* data, const int count,
                            const int dest_rank, const int tag));
    MOCK_METHOD4(recv, void(double* data, const int count,
                            const int dest_rank, const int tag));

    MOCK_METHOD4(irecv, void(int* data, const int count,
                             const int src_rank, const int tag));
    MOCK_METHOD4(irecv, void(double* data, const int count,
                             const int src_rank, const int tag));

    MOCK_METHOD3(bcast, void(int* data, const int count, const int root_rank));
    MOCK_METHOD3(bcast, void(double* data, const int count, const int root_rank));
};
