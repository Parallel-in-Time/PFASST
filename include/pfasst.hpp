#ifndef _PFASST_HPP_
#define _PFASST_HPP_

#include <iostream>
using namespace std;

#include "pfasst/config.hpp"
#include "pfasst/logging.hpp"

#ifdef WITH_MPI
  #include <mpi.h>
#endif


namespace pfasst
{
  inline static void init(int argc, char** argv,
                          std::function<void()> opts = nullptr,
                          std::function<void()> logs = nullptr,
                          bool with_mpi = false)
  {
    if (opts) {
      opts();
    }
    config::init();
    if (with_mpi) {
#ifdef WITH_MPI
      MPI_Init(&argc, &argv);
#else
      cerr << "PFASST::init() : 'with_mpi' flag used without enabling MPI" << endl;
#endif
    }
    config::read_commandline(argc, argv);
    log::start_log(argc, argv);
    if (logs) {
      logs();
    }
  }
} // ::pfasst

#endif
