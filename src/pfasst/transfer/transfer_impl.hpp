#include "pfasst/transfer/transfer.hpp"

#include <memory>
#include <stdexcept>
using namespace std;

#include "pfasst/globals.hpp"


namespace pfasst
{
  template<class TransferTraits, typename Enabled>
  void
  Transfer<TransferTraits, Enabled>::interpolate_initial(const shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse,
                                                         shared_ptr<typename TransferTraits::fine_sweeper_type> fine)
  {
    UNUSED(coarse); UNUSED(fine);
    throw runtime_error("interpolation of initial values for generic Sweeper");
  }

  template<class TransferTraits, typename Enabled>
  void
  Transfer<TransferTraits, Enabled>::interpolate(const shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse,
                                                 shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                                                 const bool initial)
  {
    UNUSED(coarse); UNUSED(fine); UNUSED(initial);
    throw runtime_error("interpolation for generic Sweeper");
  }

  template<class TransferTraits, typename Enabled>
  void
  Transfer<TransferTraits, Enabled>::interpolate_data(const shared_ptr<typename TransferTraits::coarse_encap_type> coarse,
                                                      shared_ptr<typename TransferTraits::fine_encap_type> fine)
  {
    UNUSED(coarse); UNUSED(fine);
    throw runtime_error("interpolation for generic Encapsulations");
  }

  template<class TransferTraits, typename Enabled>
  void
  Transfer<TransferTraits, Enabled>::restrict_initial(const shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                                                      shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse)
  {
    UNUSED(coarse); UNUSED(fine);
    throw runtime_error("restriction of initial value for generic Sweeper");
  }

  template<class TransferTraits, typename Enabled>
  void
  Transfer<TransferTraits, Enabled>::restrict(const shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                                              shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse,
                                              const bool initial)
  {
    UNUSED(coarse); UNUSED(fine); UNUSED(initial);
    throw runtime_error("restriction for generic Sweeper");
  }

  template<class TransferTraits, typename Enabled>
  void
  Transfer<TransferTraits, Enabled>::restrict_data(const shared_ptr<typename TransferTraits::fine_encap_type> fine,
                                                   shared_ptr<typename TransferTraits::coarse_encap_type> coarse)
  {
    UNUSED(coarse); UNUSED(fine);
    throw runtime_error("restriction for generic Encapsulations");
  }

  template<class TransferTraits, typename Enabled>
  void
  Transfer<TransferTraits, Enabled>::fas(const typename TransferTraits::fine_time_type& dt,
                                         const shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                                         shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse)
  {
    UNUSED(dt); UNUSED(coarse); UNUSED(fine);
    throw runtime_error("FAS correction for generic Sweeper");
  }
}  // ::pfasst
