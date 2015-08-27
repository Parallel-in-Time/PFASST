#ifndef _PFASST__TRANSFER__INTERFACE_HPP_
#define _PFASST__TRANSFER__INTERFACE_HPP_

#include <memory>
#include <type_traits>
using namespace std;

#include "pfasst/transfer/traits.hpp"


namespace pfasst
{
  template<
    class TransferTraits,
    typename Enabled = void
  >
  class Transfer
    : public enable_shared_from_this<Transfer<TransferTraits, Enabled>>
  {
    public:
      typedef          TransferTraits              traits;

      typedef typename traits::coarse_encap_traits coarse_encap_traits;
      typedef typename traits::coarse_encap_type   coarse_encap_type;
      typedef typename traits::coarse_time_type    coarse_time_type;
      typedef typename traits::coarse_spacial_type coarse_spacial_type;

      typedef typename traits::fine_encap_traits   fine_encap_traits;
      typedef typename traits::fine_encap_type     fine_encap_type;
      typedef typename traits::fine_time_type      fine_time_type;
      typedef typename traits::fine_spacial_type   fine_spacial_type;

      static_assert(is_convertible<coarse_time_type, fine_time_type>::value,
                    "Coarse Time Type must be convertible to Fine Time Type");
      static_assert(is_convertible<fine_time_type, coarse_time_type>::value,
                    "Fine Time Type must be convertible to Coarse Time Type");

    public:
      Transfer() = default;
      Transfer(const Transfer<TransferTraits, Enabled>& other) = default;
      Transfer(Transfer<TransferTraits, Enabled>&& other) = default;
      virtual ~Transfer() = default;
      Transfer<TransferTraits, Enabled>& operator=(const Transfer<TransferTraits, Enabled>& other) = default;
      Transfer<TransferTraits, Enabled>& operator=(Transfer<TransferTraits, Enabled>&& other) = default;

      virtual void interpolate_initial(const shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse,
                                       shared_ptr<typename TransferTraits::fine_sweeper_type> fine);
      virtual void interpolate(const shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse,
                               shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                               const bool initial = false);
      virtual void interpolate_data(const shared_ptr<typename TransferTraits::coarse_encap_type> coarse,
                                    shared_ptr<typename TransferTraits::fine_encap_type> fine);

      virtual void restrict_initial(const shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                                    shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse);
      virtual void restrict(const shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                            shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse,
                            const bool initial = false);
      virtual void restrict_data(const shared_ptr<typename TransferTraits::fine_encap_type> fine,
                                 shared_ptr<typename TransferTraits::coarse_encap_type> coarse);

      virtual void fas(const typename TransferTraits::fine_time_type& dt,
                       const shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                       shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse);
  };
}  // ::pfasst

#include "pfasst/transfer/transfer_impl.hpp"

#endif  // _PFASST__TRANSFER__INTERFACE_HPP_
