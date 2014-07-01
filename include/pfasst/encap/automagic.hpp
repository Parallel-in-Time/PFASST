
#ifndef _PFASST_ENCAP_AUTOMAGIC_HPP_
#define _PFASST_ENCAP_AUTOMAGIC_HPP_

#include <tuple>

#include "../quadrature.hpp"
#include "encap_sweeper.hpp"

using namespace std;

namespace pfasst
{
  namespace encap
  {

    template<typename ScalarT, typename timeT>
    using AutoBuildTupleT = tuple<pfasst::encap::EncapSweeper<ScalarT, timeT>*,
                                  pfasst::ITransfer*,
                                  pfasst::encap::EncapFactory<ScalarT, timeT>*>;

    template<typename ScalarT, typename timeT, typename ControllerT, typename BuildT>
    void auto_build(ControllerT& c, vector<pair<int, string>> nodes, BuildT build)
    {
      for (unsigned int l = 0; l < nodes.size(); l++) {
        auto nds = pfasst::compute_nodes<timeT>(get<0>(nodes[l]), get<1>(nodes[l]));
        AutoBuildTupleT<ScalarT, timeT> tpl = build(l);
        auto* sweeper = get<0>(tpl);
        auto* transfer = get<1>(tpl);
        auto* factory = get<2>(tpl);
        sweeper->set_nodes(nds);
        sweeper->set_factory(factory);
        c.add_level(sweeper, transfer, false);
      }
    }

    template<typename ScalarT, typename timeT, typename ControllerT, typename InitialT>
    void auto_setup(ControllerT& c, InitialT initial)
    {
      c.setup();

      for (int l = 0; l < c.nlevels(); l++) {
        auto* isweeper = c.get_level(l);
        auto* sweeper = dynamic_cast<pfasst::encap::EncapSweeper<ScalarT, timeT>*>(isweeper);
        auto* q0 = sweeper->get_state(0);
        initial(sweeper, q0);
      }
    }

  }
}

#endif
