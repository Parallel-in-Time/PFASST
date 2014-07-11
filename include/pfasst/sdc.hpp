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
	auto& steps = this->steps;
	auto& iters = this->iterations;

        for (steps.reset(); steps.valid(); steps.next()) {
          bool initial = this->get_step() == 0;
	  for (iters.reset(); iters.valid(); iters.next()) {
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
