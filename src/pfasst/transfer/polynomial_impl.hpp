#include "pfasst/transfer/polynomial.hpp"

#include <stdexcept>
#include <memory>
#include <vector>
using namespace std;

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"
#include "pfasst/quadrature.hpp"
#include "pfasst/encap/encapsulation.hpp"


namespace pfasst
{
  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::interpolate_initial(const shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse,
                                                                   shared_ptr<typename TransferTraits::fine_sweeper_type> fine)
  {
    CVLOG(1, "TRANS") << "interpolate initial value only";

    auto coarse_factory = coarse->get_encap_factory();
    auto fine_factory = fine->get_encap_factory();

    // c_delta = restrict(u_0^F) - u_0^C
    auto coarse_delta = coarse_factory->create();
    // c_delta = restrict(u_0)
    this->restrict_data(fine->get_initial_state(), coarse_delta);
    // c_delta -= c_0
    coarse_delta->scaled_add(-1.0, coarse->get_initial_state());

    // f_delta = interpolate(c_delta)
    auto fine_delta = fine_factory->create();
    // f_delta = interpolate(c_delta)
    this->interpolate_data(coarse_delta, fine_delta);

    // u_0^F = u_0^F - f_delta
    fine->initial_state()->scaled_add(-1.0, fine_delta);

    fine->reevaluate(true);
  }

  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::interpolate(const shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse,
                                                           shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                                                           const bool initial)
  {
    CVLOG(1, "TRANS") << "interpolate";

    if (coarse->get_quadrature()->left_is_node() || fine->get_quadrature()->left_is_node()) {
      CLOG(ERROR, "TRANS") << "Time interpolation with left time point as a node is still not supported.";
      throw runtime_error("time interpolation with left time point as node");
    }
    
    CLOG_IF(fine->get_quadrature()->get_num_nodes() != coarse->get_quadrature()->get_num_nodes(),
            WARNING, "TRANS") << "interpolation between different number of nodes not tested!";

    if (initial) {
      this->interpolate_initial(coarse, fine);
    }

    this->setup_tmat(fine->get_quadrature(), coarse->get_quadrature());
    CVLOG(1, "TRANS") << "tmat: " << this->tmat;

    // +1 here for additional value in states
    const size_t num_coarse_nodes = coarse->get_quadrature()->get_num_nodes() + 1;

    auto coarse_factory = coarse->get_encap_factory();
    auto fine_factory = fine->get_encap_factory();

    // compute coarse level correction
    vector<shared_ptr<fine_encap_type>> fine_deltas(num_coarse_nodes);
    generate(fine_deltas.begin(), fine_deltas.end(),
             [fine_factory]() { return fine_factory->create(); });

    auto coarse_delta = coarse_factory->create();
    // u_m^F = u_m^F - interpolate(u_m^C - prev_u_m^C)
    for (size_t m = 1; m < num_coarse_nodes; ++m) {
      coarse_delta->data() = coarse->get_states()[m]->get_data();
      coarse_delta->scaled_add(-1.0, coarse->get_previous_states()[m]);
      CVLOG(1, "TRANS") << "  cd["<<m<<"]: " << to_string(coarse_delta);
      this->interpolate_data(coarse_delta, fine_deltas[m]);
      CVLOG(1, "TRANS") << "  fd["<<m<<"]: " << to_string(fine_deltas[m]);
    }

    CVLOG(1, "TRANS") << "fine states and deltas before interpolation:";
    for (size_t m = 1; m < num_coarse_nodes; ++m) {
      CVLOG(1, "TRANS") << "  f["<<m<<"]:  " << to_string(fine->get_states()[m]);
      CVLOG(1, "TRANS") << "  fd["<<m<<"]: " << to_string(fine_deltas[m]);
    }
    // add coarse level correction onto fine level's states
    encap::mat_apply(fine->states(), 1.0, this->tmat, fine_deltas, false);

    CVLOG(1, "TRANS") << "fine states after interpolation:";
    for (auto& n : fine->states()) {
      CVLOG(1, "TRANS") << "  " << to_string(n);
    }

    // update function evaluations on fine level
    fine->reevaluate();
  }

  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::interpolate_data(const shared_ptr<typename TransferTraits::coarse_encap_type> coarse,
                                                                shared_ptr<typename TransferTraits::fine_encap_type> fine)
  {
    UNUSED(coarse); UNUSED(fine);
    throw runtime_error("interpolation for generic Encapsulations");
  }

  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::restrict_initial(const shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                                                                shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse)
  {
    CVLOG(1, "TRANS") << "restrict initial value only";

    this->restrict_data(fine->get_initial_state(), coarse->initial_state());
  }

  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::restrict(const shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                                                        shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse,
                                                        const bool initial)
  {
    CVLOG(1, "TRANS") << "restrict";

    if (coarse->get_quadrature()->left_is_node() || fine->get_quadrature()->left_is_node()) {
      CLOG(ERROR, "TRANS") << "Time restriction with left time point as a node is still not supported.";
      throw runtime_error("time restriction with left time point as node");
    }


    if (initial) {
      this->restrict_initial(fine, coarse);
    }

    // +1 here for additional value in states
    const auto coarse_nodes = coarse->get_quadrature()->get_nodes();
    const auto fine_nodes = fine->get_quadrature()->get_nodes();
    const size_t num_coarse_nodes = coarse->get_quadrature()->get_num_nodes() + 1;

    // this commented out stuff is probably required for non-equal sets of time nodes
//     const int factor = ((int)num_fine_nodes - 1) / ((int)num_coarse_nodes - 1);

    for (size_t m = 1; m < num_coarse_nodes; ++m) {
//       if (coarse_nodes[m] != fine_nodes[m * factor]) {
//         CLOG(ERROR, "TRANS") << "coarse nodes are not nested within fine ones."
//                              << "coarse: " << coarse_nodes << " fine: " << fine_nodes;
//         throw NotImplementedYet("non-nested nodes");
//       }
      this->restrict_data(fine->get_states()[m], coarse->states()[m]);
    }

    coarse->reevaluate();
  }

  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::restrict_data(const shared_ptr<typename TransferTraits::fine_encap_type> fine,
                                                             shared_ptr<typename TransferTraits::coarse_encap_type> coarse)
  {
    UNUSED(coarse); UNUSED(fine);
    throw runtime_error("restriction for generic Encapsulations");
  }

  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::fas(const typename TransferTraits::fine_time_type& dt,
                                                   const shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                                                   shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse)
  {
    CVLOG(1, "TRANS") << "compute FAS correction";

    const auto coarse_nodes = coarse->get_quadrature()->get_nodes();
    const auto fine_nodes = fine->get_quadrature()->get_nodes();
    const size_t num_coarse_nodes = coarse->get_quadrature()->get_num_nodes();
    const size_t num_fine_nodes = fine->get_quadrature()->get_num_nodes();

    if (num_fine_nodes != num_coarse_nodes) {
      CLOG(ERROR, "TRANS") << "FAS Correction for different time scales not supported yet.";
      throw runtime_error("fas correction for different time scales");
    } else if (!equal(coarse_nodes.cbegin(), coarse_nodes.cend(), fine_nodes.cbegin())) {
      CLOG(ERROR, "TRANS") << "FAS Correction for different time nodes not supported yet.";
      throw runtime_error("fas correction for different sets of nodes");
    }
    assert(coarse_nodes.size() == fine_nodes.size());

    auto coarse_factory = coarse->get_encap_factory();
    auto fine_factory = fine->get_encap_factory();

    vector<shared_ptr<coarse_encap_type>> fas(num_coarse_nodes + 1);
    generate(fas.begin(), fas.end(), [coarse_factory]() { return coarse_factory->create(); });

    const auto coarse_integral = coarse->integrate(dt);
    const auto fine_integral = fine->integrate(dt);

    for (size_t m = 0; m < num_coarse_nodes + 1; ++m) {
      this->restrict_data(fine_integral[m], fas[m]);
      fas[m]->scaled_add(-1.0, coarse_integral[m]);
      coarse->tau()[m]->data() = fas[m]->get_data();
    }
  }


  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::setup_tmat(const shared_ptr<quadrature::IQuadrature<typename TransferTraits::fine_time_type>> fine_quad,
                                                          const shared_ptr<quadrature::IQuadrature<typename TransferTraits::coarse_time_type>> coarse_quad)
  {
    if (this->tmat.rows() == 0) {
      auto coarse_nodes = coarse_quad->get_nodes();
      auto fine_nodes = fine_quad->get_nodes();

      coarse_nodes.insert(coarse_nodes.begin(), 0.0);
      fine_nodes.insert(fine_nodes.begin(), 0.0);

      CLOG_IF(!equal(coarse_nodes.cbegin(), coarse_nodes.cend(), fine_nodes.cbegin()),
              WARNING, "TRANS") << "interpolation for different sets of nodes not tested.";

      this->tmat = quadrature::compute_interp<typename TransferTraits::fine_time_type>(coarse_nodes,
                                                                                       fine_nodes);
    }
  }
}  // ::pfasst
