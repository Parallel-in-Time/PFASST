#include "fixtures/test_helpers.hpp"

#include <memory>
using namespace std;

#include "pfasst/sweeper/interface.hpp"


template<
  typename precision,
  class EncapT
>
class SweeperMock
  : public pfasst::Sweeper<precision, EncapT>
{
  public:
    typedef precision precision_type;
    typedef EncapT    encap_type;

  public:
//     MOCK_METHOD0_T(controller, shared_ptr<Controller<precision, EncapT>>&());
//     MOCK_CONST_METHOD0_T(get_controller, shared_ptr<Controller<precision, EncapT>>());
};
