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
  inline static void init(int argc, char** argv)
  {
    SDC<>::enable_config_options(0);
    encap::EncapSweeper<>::enable_config_options(0);
    Quadrature::enable_config_options(0);
    config::init_config();
    log::start_log(argc, argv);
    config::read_commandline(argc, argv);
  }
} // ::pfasst

#endif
