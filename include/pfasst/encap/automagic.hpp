
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

    template<typename time = time_precision>
    using auto_build_tuple = tuple<pfasst::encap::EncapSweeper<time>*,
          pfasst::ITransfer<time>*,
          pfasst::encap::EncapFactory<time>*>;

    template<typename time = time_precision, typename controllerT, typename BuildT>
    void auto_build(controllerT& c, vector<pair<size_t, string>> nodes, BuildT build)
    {
      for (size_t l = 0; l < nodes.size(); l++) {
        auto nds = pfasst::compute_nodes<time>(get<0>(nodes[l]), get<1>(nodes[l]));
        auto_build_tuple<time> tpl = build(l);
        auto* sweeper = get<0>(tpl);
        auto* transfer = get<1>(tpl);
        auto* factory = get<2>(tpl);
        sweeper->set_nodes(nds);
        sweeper->set_factory(factory);
        c.add_level(sweeper, transfer, false);
      }
    }

    template<typename time = time_precision, typename controllerT, typename initialT>
    void auto_setup(controllerT& c, initialT initial)
    {
      c.setup();
      for (size_t l = 0; l < c.nlevels(); l++) {
        auto* isweeper = c.get_level(l);
        auto* sweeper = dynamic_cast<pfasst::encap::EncapSweeper<time>*>(isweeper);
        auto* q0 = sweeper->get_state(0);
        initial(sweeper, q0);
      }
    }

  }
}

#endif
