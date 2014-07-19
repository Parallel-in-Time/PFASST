/*
 * Vanilla SDC controller.
 */

#ifndef _PFASST_SDC_HPP_
#define _PFASST_SDC_HPP_

#include <iostream>
using namespace std;

#include "controller.hpp"

namespace pfasst
{

  template<typename time = time_precision>
  class SDC
    : public Controller<time>
  {
    public:

      void run()
      {
        auto sweeper = this->get_level(0);

        for (; this->get_time() < this->get_end_time(); this->advance_time()) {
          bool initial = this->get_step() == 0;
	  for (this->set_iteration(0); this->get_iteration() < this->get_max_iterations(); this->advance_iteration()) {
	    bool predict = this->get_iteration() == 0;
	    if (predict) {
	      sweeper->predict(initial);
	    } else {
	      sweeper->sweep();
	    }
          }
          sweeper->advance();
        }
      }

  };

}

#endif
