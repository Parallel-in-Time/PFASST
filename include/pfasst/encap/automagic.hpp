
#ifndef _PFASST_ENCAP_AUTOMAGIC_HPP_
#define _PFASST_ENCAP_AUTOMAGIC_HPP_

#include <tuple>
#include <cassert>

#include "../quadrature.hpp"
#include "encap_sweeper.hpp"

using namespace std;

namespace pfasst
{
  namespace encap
  {

    template<typename time = time_precision>
    using AutoBuildTuple = tuple<
                             shared_ptr<pfasst::encap::EncapSweeper<time>>,
                             shared_ptr<pfasst::ITransfer<time>>,
                             shared_ptr<pfasst::encap::EncapFactory<time>>
                           >;

    template<typename ControllerT, typename BuildT, typename time = time_precision>
    void auto_build(ControllerT& c, vector<pair<size_t, quadrature::QuadratureType>> nodes, BuildT build)
    {
      for (size_t l = 0; l < nodes.size(); l++) {
        auto quad = quadrature::quadrature_factory<time>(get<0>(nodes[l]), get<1>(nodes[l]));
        AutoBuildTuple<time> tpl = build(l);
        auto sweeper = get<0>(tpl);
        auto transfer = get<1>(tpl);
        auto factory = get<2>(tpl);
        sweeper->set_quadrature(quad);
        sweeper->set_factory(factory);
        c.add_level(sweeper, transfer, false);
      }
    }

    template<typename ControllerT, typename initialT, typename time = time_precision>
    void auto_setup(ControllerT& c, initialT initial)
    {
      c.setup();
      for (size_t l = 0; l < c.nlevels(); l++) {
        auto isweeper = c.get_level(l);
        auto sweeper = dynamic_pointer_cast<pfasst::encap::EncapSweeper<time>>(isweeper);
        assert(sweeper);
        auto q0 = sweeper->get_start_state();
        initial(sweeper, q0);
      }
    }

  }  // ::pfasst::encap
} // ::pfasst

#endif
