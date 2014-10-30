#ifndef _PFASST_HPP_
#define _PFASST_HPP_

#include "pfasst/config.hpp"
#include "pfasst/interfaces.hpp"
#include "pfasst/quadrature.hpp"
#include "pfasst/sdc.hpp"

namespace pfasst
{
  static void init(int argc, char** argv)
  {
    SDC<>::enable_config_options(0);
    quadrature::Quadrature::enable_config_options(0);
    config::init_config();
    config::read_commandline(argc, argv);
  }
} // ::pfasst

#endif
