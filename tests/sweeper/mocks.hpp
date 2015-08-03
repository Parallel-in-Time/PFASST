#include "fixtures/test_helpers.hpp"

#include <memory>
using namespace std;

#include "pfasst/sweeper/interface.hpp"


template<
  class SweeperTrait,
  typename Enabled = void
>
class SweeperMock
  : public pfasst::Sweeper<SweeperTrait, Enabled>
{
  public:
    typedef          SweeperTrait traits;

  public:
    MOCK_METHOD0_T(quadrature, shared_ptr<IQuadrature<typename SweeperTrait::time_type>>&());
    MOCK_CONST_METHOD0_T(get_quadrature, const shared_ptr<IQuadrature<typename SweeperTrait::time_type>>());

    MOCK_METHOD0_T(status, shared_ptr<pfasst::Status<typename SweeperTrait::time_type>>&());
    MOCK_CONST_METHOD0_T(get_status, const shared_ptr<pfasst::Status<typename SweeperTrait::time_type>>());

    MOCK_METHOD0_T(encap_factory, shared_ptr<typename SweeperTrait::encap_type::factory_type>&());
    MOCK_CONST_METHOD0_T(get_encap_factory, const shared_ptr<typename SweeperTrait::encap_type::factory_type>());

    MOCK_METHOD0_T(setup, void());

    MOCK_METHOD0_T(pre_predict, void());
    MOCK_METHOD0_T(predict, void());
    MOCK_METHOD0_T(post_predict, void());

    MOCK_METHOD0_T(pre_sweep, void());
    MOCK_METHOD0_T(sweep, void());
    MOCK_METHOD0_T(post_sweep, void());

    MOCK_METHOD0_T(advance, void());
    MOCK_METHOD0_T(spread, void());
    MOCK_METHOD0_T(save, void());

    MOCK_METHOD1_T(reevaluate, void(const bool initial_only));
    MOCK_METHOD1_T(integrate, vector<shared_ptr<typename SweeperTrait::encap_type>>(const typename SweeperTrait::time_type& dt));

    MOCK_METHOD0_T(converged, bool());
};
