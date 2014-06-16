
#ifndef _PFASST_ENCAP_AUTOMAGIC_HPP_
#define _PFASST_ENCAP_AUTOMAGIC_HPP_

#include <tuple>

#include "../quadrature.hpp"
#include "encapsulation.hpp"

using namespace std;

namespace pfasst {
  namespace encap {

    template<typename scalar, typename time>
    using auto_build_tuple = tuple<pfasst::encap::EncapsulatedSweeperMixin<scalar,time>*,
				   pfasst::ITransfer*,
				   pfasst::encap::EncapsulationFactory<scalar,time>*>;

    template<typename scalar, typename time, typename controllerT, typename buildT>
    void auto_build(controllerT& c, vector<pair<int,string>> nodes, buildT build) {
      for (unsigned int l=0; l<nodes.size(); l++) {
	auto nds = pfasst::compute_nodes<time>(get<0>(nodes[l]), get<1>(nodes[l]));
	auto_build_tuple<scalar,time> tpl = build(l);
	auto* sweeper = get<0>(tpl);
	auto* transfer = get<1>(tpl);
	auto* factory = get<2>(tpl);
	sweeper->set_nodes(nds);
	sweeper->set_factory(factory);
	c.add_level(sweeper, transfer, false);
      }
    }

    template<typename scalar, typename time, typename controllerT, typename initialT>
    void auto_setup(controllerT& c, initialT initial) {
      c.setup();
      for (int l=0; l<c.nlevels(); l++) {
	//	auto* sweeper = c.get_level<pfasst::encap::EncapsulatedSweeperMixin<scalar,time>>(l);
	auto* isweeper = c.get_level(l);
	auto* sweeper = dynamic_cast<pfasst::encap::EncapsulatedSweeperMixin<scalar,time>*>(isweeper);
	auto* q0 = sweeper->get_state(0);
	initial(sweeper, q0);
	// sweeper->exact(q0, 0.0);
      }
    }

  }
}

#endif
