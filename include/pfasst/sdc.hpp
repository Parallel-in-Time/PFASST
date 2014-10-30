/*
 * Vanilla SDC controller.
 */

#ifndef _PFASST_SDC_HPP_
#define _PFASST_SDC_HPP_

#include <iostream>
using namespace std;

#include "controller.hpp"
#include "config.hpp"

namespace pfasst
{

  template<typename time = time_precision>
  class SDC
    : public Controller<time>
  {
    private:
      static void init_config_options(po::options_description& opts)
      {
        opts.add_options()
          ("num_iter", po::value<size_t>(), "number of iterations")
          ("num_steps", po::value<size_t>(), "number of time steps")
          ("delte_step", po::value<time>(), "width of one time step")
          ;
      }

    public:
      static void enable_config_options(size_t index = -1)
      {
        pfasst::config::Options::get_instance()
          .register_init_function("SDC Sweeper",
                                  function<void(po::options_description&)>(init_config_options),
                                  index);
      }

    public:
      void run()
      {
        auto sweeper = this->get_level(0);

        for (; this->get_time() < this->get_end_time(); this->advance_time()) {
          bool initial = this->get_step() == 0;
          for (this->set_iteration(0);
               this->get_iteration() < this->get_max_iterations();
               this->advance_iteration()) {
            bool predict = this->get_iteration() == 0;
            if (predict) {
              sweeper->predict(initial);
            } else {
              sweeper->sweep();
            }
            if (sweeper->converged()) {
              break;
            }
          }
          sweeper->advance();
        }
      }
  };

}  // ::pfasst

#endif
