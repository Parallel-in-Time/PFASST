/**
 * @file controller/sdc.hpp
 * @since v0.1.0
 */
#ifndef _PFASST__CONTROLLER__SDC_HPP_
#define _PFASST__CONTROLLER__SDC_HPP_

#include "pfasst/controller/interface.hpp"


namespace pfasst
{
  /**
   * Vanilla SDC controller.
   *
   * @tparam time time precision
   *
   * @see Controller on how to set up the controller and feed it with sweepers to do the actual
   *   integration.
   *
   * @ingroup Controllers
   */
  template<typename time = time_precision>
  class SDC
    : public Controller<time>
  {
    public:
      /**
       * Run vanilla SDC.
       */
      virtual void run();
  };
}  // ::pfasst

#include "pfasst/controller/sdc_impl.hpp"

#endif  // _PFASST__CONTROLLER__SDC_HPP_
