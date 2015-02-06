#ifndef _PFASST_HPP_
#define _PFASST_HPP_

#include "pfasst/config.hpp"
#include "pfasst/logging.hpp"
#include "pfasst/interfaces.hpp"
#include "pfasst/quadrature.hpp"
#include "pfasst/sdc.hpp"
#include "pfasst/encap/encap_sweeper.hpp"

namespace pfasst
{
  inline static void init(int argc, char** argv, std::function<void()> opts=nullptr, std::function<void()> logs=nullptr)
  {
    if (opts) {
      opts();
    }
    config::init();
    log::start_log(argc, argv);
    if (logs) {
      logs();
    }
    config::read_commandline(argc, argv);
  }
} // ::pfasst

#endif
