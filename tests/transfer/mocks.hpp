#include "fixtures/test_helpers.hpp"

#include <memory>
using namespace std;

#include "pfasst/transfer/interface.hpp"


template<
  class TransferTraits,
  typename Enabled = void
>
class TransferMock
  : public pfasst::Transfer<TransferTraits, Enabled>
{
  public:
    typedef          TransferTraits              traits;

  public:
    MOCK_METHOD2_T(interpolate_initial, void(const shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse,
                                             shared_ptr<typename TransferTraits::fine_sweeper_type> fine));
    MOCK_METHOD3_T(interpolate, void(const shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse,
                                     shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                                     const bool initial));
    MOCK_METHOD2_T(interpolate_data, void(const shared_ptr<typename TransferTraits::coarse_encap_type> coarse,
                                          shared_ptr<typename TransferTraits::fine_encap_type> fine));

    MOCK_METHOD2_T(restrict_initial, void(const shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                                    shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse));
    MOCK_METHOD3_T(restrict, void(const shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                                  shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse,
                                  const bool initial));
    MOCK_METHOD2_T(restrict_data, void(const shared_ptr<typename TransferTraits::fine_encap_type> fine,
                                       shared_ptr<typename TransferTraits::coarse_encap_type> coarse));

    MOCK_METHOD3_T(fas, void(const typename TransferTraits::fine_time_type& dt,
                             const shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                             shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse));
};
