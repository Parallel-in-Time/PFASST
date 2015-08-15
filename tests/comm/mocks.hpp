#include "fixtures/test_helpers.hpp"

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
};
