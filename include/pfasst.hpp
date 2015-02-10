#ifndef _PFASST_HPP_
#define _PFASST_HPP_

#include "pfasst/config.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  inline static void init(int argc, char** argv,
                          std::function<void()> opts = nullptr,
                          std::function<void()> logs = nullptr)
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
