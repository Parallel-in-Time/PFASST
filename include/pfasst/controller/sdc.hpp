#ifndef _PFASST_SDC_HPP_
#define _PFASST_SDC_HPP_

#include "pfasst/controller/interface.hpp"


namespace pfasst
{
  /**
   * Vanilla SDC controller.
   */
  template<typename time = time_precision>
  class SDC
    : public Controller<time>
  {
    public:
      virtual void run();
  };
}  // ::pfasst

#include "pfasst/controller/sdc_impl.hpp"

#endif
