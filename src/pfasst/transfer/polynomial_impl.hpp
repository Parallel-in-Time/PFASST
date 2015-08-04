#include "pfasst/transfer/polynomial.hpp"

#include <memory>
#include <vector>
using namespace std;

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"
#include "pfasst/quadrature.hpp"
#include "pfasst/encap/interface.hpp"


namespace pfasst
{
  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::interpolate_initial(const shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse,
                                                                   shared_ptr<typename TransferTraits::fine_sweeper_type> fine)
  {
    auto coarse_factory = coarse->get_encap_factory();
    auto fine_factory = fine->get_encap_factory();

    // c_delta
    auto coarse_delta = coarse_factory->create();
    // c_delta = restrict(f_0)
    this->restrict_data(fine->get_initial_state(), coarse_delta);
    // c_delta -= c_0
    coarse_delta->scaled_add(-1.0, coarse->get_initial_state());

    // f_delta
    auto fine_delta = fine_factory->create();
    // f_delta = interpolate(c_delta)
    this->interpolate_data(coarse_delta, fine_delta);
    // f_0 -= f_delta
    fine->initial_state()->scaled_add(-1.0, fine_delta);

    fine->reevaluate(true);
  }

  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::interpolate(const shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse,
                                                           shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                                                           const bool initial)
  {
    if (initial) {
      this->interpolate_initial(coarse, fine);
    }

    this->setup_tmat(fine->get_quadrature()->get_nodes(), coarse->get_quadrature()->get_nodes());

    const size_t num_fine_nodes = fine->get_quadrature()->get_num_nodes();
    const size_t num_coarse_nodes = coarse->get_quadrature()->get_num_nodes();

    auto coarse_factory = coarse->get_encap_factory();
    auto fine_factory = fine->get_encap_factory();

    vector<shared_ptr<fine_encap_type>> fine_delta(num_coarse_nodes);
    generate(fine_delta.begin(), fine_delta.end(),
             [fine_factory]() { return fine_factory->create(); });
    auto coarse_delta = coarse_factory->create();

    for (size_t m = 0; m < num_coarse_nodes; ++m) {
      coarse_delta = coarse->get_states()[m];
      coarse_delta->scaled_add(-1.0, coarse->get_previous_states()[m]);
      this->interpolate_data(coarse_delta, fine_delta[m]);
    }

    encap::mat_apply(fine->states(), 1.0, this->tmat, fine_delta, true);

    fine->reevaluate();
  }

  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::interpolate_data(const shared_ptr<typename TransferTraits::coarse_encap_type> coarse,
                                                                shared_ptr<typename TransferTraits::fine_encap_type> fine)
  {
    UNUSED(coarse); UNUSED(fine);
    throw NotImplementedYet("interpolation for generic Encapsulations");
  }

  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::restrict_initial(const shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                                                                shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse)
  {
    this->restrict_data(fine->get_initial_state(), coarse->initial_state());
  }

  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::restrict(const shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                                                        shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse,
                                                        const bool initial)
  {
    if (initial) {
      this->restrict_initial(fine, coarse);
    }

    const auto coarse_nodes = coarse->get_quadrature()->get_nodes();
    const auto fine_nodes = fine->get_quadrature()->get_nodes();
    const size_t num_coarse_nodes = coarse->get_quadrature()->get_num_nodes();
    const size_t num_fine_nodes = fine->get_quadrature()->get_num_nodes();

    const int factor = ((int)num_fine_nodes - 1) / ((int)num_coarse_nodes - 1);

    for (size_t m = 0; m < num_coarse_nodes; ++m) {
      if (coarse_nodes[m] != fine_nodes[m * factor]) {
        CLOG(ERROR, "TRANS") << "coarse nodes are not nested within fine ones."
                             << "coarse: " << coarse_nodes << " fine: " << fine_nodes;
        throw NotImplementedYet("non-nested nodes");
      }
      this->restrict_data(fine->get_states()[m * factor], coarse->states()[m]);
    }

    coarse->reevaluate();
  }

  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::restrict_data(const shared_ptr<typename TransferTraits::fine_encap_type> fine,
                                                             shared_ptr<typename TransferTraits::coarse_encap_type> coarse)
  {
    UNUSED(coarse); UNUSED(fine);
    throw NotImplementedYet("restriction for generic Encapsulations");
  }

  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::fas(const typename TransferTraits::fine_time_type& dt,
                                                   const shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                                                   shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse)
  {
    const auto coarse_nodes = coarse->get_quadrature()->get_nodes();
    const auto fine_nodes = fine->get_quadrature()->get_nodes();
    const size_t num_coarse_nodes = coarse->get_quadrature()->get_num_nodes();
    const size_t num_fine_nodes = fine->get_quadrature()->get_num_nodes();

    auto coarse_factory = coarse->get_encap_factory();
    auto fine_factory = fine->get_encap_factory();

    vector<shared_ptr<coarse_encap_type>> coarse_integral(num_coarse_nodes);
    vector<shared_ptr<fine_encap_type>> fine_integral(num_fine_nodes);
    vector<shared_ptr<coarse_encap_type>> restricted_fine_integral(num_coarse_nodes);

    generate(coarse_integral.begin(), coarse_integral.end(),
             [coarse_factory]() { return coarse_factory->create(); });
    generate(fine_integral.begin(), fine_integral.end(),
             [fine_factory]() { return fine_factory->create(); });
    generate(restricted_fine_integral.begin(), restricted_fine_integral.end(),
             [coarse_factory]() { return coarse_factory->create(); });

    coarse_integral = coarse->integrate(dt);
    fine_integral = fine->integrate(dt);

    const int factor = ((int)num_fine_nodes - 1) / ((int)num_coarse_nodes - 1);

    for (size_t m = 0; m < num_coarse_nodes; ++m) {
      if (coarse_nodes[m] != fine_nodes[m * factor]) {
        CLOG(ERROR, "TRANS") << "coarse nodes are not nested within fine ones."
                             << "coarse: " << coarse_nodes << " fine: " << fine_nodes;
        throw NotImplementedYet("non-nested nodes");
      }
      this->restrict_data(fine_integral[m * factor], restricted_fine_integral[m]);
    }

    vector<shared_ptr<coarse_encap_type>> restricted_and_coarse(2 * num_coarse_nodes);
    for (size_t m = 0; m < num_coarse_nodes; ++m) {
      restricted_and_coarse[m] = restricted_fine_integral[m];
      restricted_and_coarse[num_coarse_nodes + m] = coarse_integral[m];
    }

    // TODO: FAS w.t.r. Q not S
    this->setup_fmat(num_coarse_nodes);

    encap::mat_apply(coarse->tau(), 1.0, fmat, restricted_and_coarse, true);
  }


  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::setup_tmat(const vector<typename TransferTraits::fine_time_type>& fine_nodes,
                                                          const vector<typename TransferTraits::coarse_time_type>& coarse_nodes)
  {
    if (this->tmat.rows() == 0) {
      this->tmat = quadrature::compute_interp<typename TransferTraits::fine_time_type>(fine_nodes,
                                                                                       coarse_nodes);
    }
  }

  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::setup_fmat(const size_t& num_coarse)
  {
    if (this->fmat.rows() == 0) {
      // XXX: +1 ?!
      this->fmat.resize(num_coarse + 1, 2 * (num_coarse + 1));
      this->fmat.fill(0.0);

      // XXX: m=1...num_coarse ?!
      for (size_t m = 1; m < num_coarse; m++) {
        this->fmat(m, m) = 1.0;
        this->fmat(m, num_coarse + m) = -1.0;

        // subtract 0-to-(m-1) FAS so resulting FAS is (m-1)-to-m FAS,
        //  which will be required in the sweeper logic
//         for (size_t n = 0; n < m; n++) {
//           this->fmat(m, n) = -1.0;
//           this->fmat(m, num_coarse + n) = 1.0;
//         }
      }
    }
  }
}  // ::pfasst
