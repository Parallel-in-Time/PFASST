/*
 * Vanilla SDC controller.
 */

#ifndef _PFASST_SDC_HPP_
#define _PFASST_SDC_HPP_

#include <memory>
#include <string>

#include "controller.hpp"
#include "quadrature.hpp"

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

        for (size_t nstep = 0; nstep < this->nsteps; nstep++) {
          time t = nstep * this->dt;

          sweeper->predict(t, this->dt, nstep == 0);
          for (size_t niter = 1; niter < this->niters; niter++) {
            sweeper->sweep(t, this->dt);
          }

          sweeper->advance();
        }
      }
  };


  template<typename sweeperT, 
           typename dataFactoryT,
           typename time = time_precision>
  SDC<time> sdc_factory(
      size_t niters
    , size_t nsteps
    , time dt
    , size_t ndofs
    , size_t nnodes
    , string node_type
    , time iv
  )
  {
    SDC<time> sdc;

    auto nodes = compute_nodes(nnodes, node_type);
    auto factory = make_shared<dataFactoryT>(ndofs);
    auto sweeper = make_shared<sweeperT>(ndofs);

    sweeper->set_nodes(nodes);
    sweeper->set_factory(factory);

    sdc.add_level(sweeper);
    sdc.set_duration(dt, nsteps, niters);
    sdc.setup();

    auto q0 = sweeper->get_state(0);
    sweeper->exact(q0, iv);

    return sdc;
  }

}

#endif
