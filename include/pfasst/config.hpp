
#ifndef _PFASST_CONFIG_HPP_
#define _PFASST_CONFIG_HPP_

namespace pfasst
{
#ifndef PFASST_TIME_PRECISION
  using time = double;
#else
  using time = PFASST_TIME_PRECISION;
#endif
}

#endif
