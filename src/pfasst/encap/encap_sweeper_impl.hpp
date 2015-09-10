#include "pfasst/encap/encap_sweeper.hpp"

#include <algorithm>
#include <cassert>
using namespace std;

#include "pfasst/globals.hpp"
#include "pfasst/config.hpp"


namespace pfasst
{
  namespace encap
  {
    template<typename time>
    EncapSweeper<time>::EncapSweeper()
      :   quadrature(nullptr)
        , abs_residual_tol(0.0)
        , rel_residual_tol(0.0)
    {}

    template<typename time>
    shared_ptr<Encapsulation<time>> EncapSweeper<time>::get_state(size_t m) const
    {
      return this->state[m];
    }

    template<typename time>
    shared_ptr<Encapsulation<time>> EncapSweeper<time>::get_tau(size_t m) const
    {
      return this->fas_corrections[m];
    }

    template<typename time>
    shared_ptr<Encapsulation<time>> EncapSweeper<time>::get_saved_state(size_t m) const
    {
      return this->saved_state[m];
    }

    /**
     * @internals
     * Sets EncapSweeper::abs_residual_tol and EncapSweeper::rel_residual_tol from command line or 
     * config file.
     * In case the parameters are not found, takes EncapSweeper::abs_residual_tol or 
     * EncapSweeper::rel_residual_tol respectively as default values.
     *
     * @todo Consider using EncapSweeper::set_residual_tolerances() instead of direct member access
     *   to allow for easy integration of future consistency checks on setting the residual
     *   tolerances.
     * @endinternals
     */
    template<typename time>
    void EncapSweeper<time>::set_options()
    {
      this->abs_residual_tol = time(config::get_value<double>("abs_res_tol", this->abs_residual_tol));
      this->rel_residual_tol = time(config::get_value<double>("rel_res_tol", this->rel_residual_tol));
    }

    /**
     * @internals
     * Initializes members holding state information by repeating calls to EncapFactory::create()
     * for EncapSweeper::start_state, EncapSweeper::end_state, EncapSweeper::state,
     * EncapSweeper::saved_state and - if @p coarse is `true` - also for 
     * EncapSweeper::fas_corrections.
     * @endinternals
     */
    template<typename time>
    void EncapSweeper<time>::setup(bool coarse)
    {
      auto const nodes = this->quadrature->get_nodes();
      auto const num_nodes = this->quadrature->get_num_nodes();

      this->start_state = this->get_factory()->create(pfasst::encap::solution);
      this->end_state = this->get_factory()->create(pfasst::encap::solution);

      for (size_t m = 0; m < num_nodes; m++) {
        this->state.push_back(this->get_factory()->create(pfasst::encap::solution));
        if (coarse) {
          this->saved_state.push_back(this->get_factory()->create(pfasst::encap::solution));
        }
      }

      if (coarse) {
        for (size_t m = 0; m < num_nodes; m++) {
          this->fas_corrections.push_back(this->get_factory()->create(pfasst::encap::solution));
        }
      }
    }

    /**
     * This implementation simply copies the initial value (see EncapSweeper::get_start_state())
     * to all time nodes.
     */
    template<typename time>
    void EncapSweeper<time>::spread()
    {
      for (size_t m = 1; m < this->quadrature->get_num_nodes(); m++) {
        this->state[m]->copy(this->state[0]);
      }
    }

    template<typename time>
    void EncapSweeper<time>::save(bool initial_only)
    {
      // XXX: if !left_is_node, this is a problem...
      if (initial_only) {
        this->saved_state[0]->copy(state[0]);
      } else {
        for (size_t m = 0; m < this->saved_state.size(); m++) {
          this->saved_state[m]->copy(state[m]);
        }
      }
    }

    template<typename time>
    void EncapSweeper<time>::set_quadrature(shared_ptr<IQuadrature<time>> quadrature)
    {
      this->quadrature = quadrature;
    }

    template<typename time>
    shared_ptr<const IQuadrature<time>> EncapSweeper<time>::get_quadrature() const
    {
      return this->quadrature;
    }

    template<typename time>
    shared_ptr<Encapsulation<time>> EncapSweeper<time>::get_start_state() const
    {
      return this->start_state;
    }

    template<typename time>
    const vector<time> EncapSweeper<time>::get_nodes() const
    {
      return this->quadrature->get_nodes();
    }

    template<typename time>
    void EncapSweeper<time>::set_factory(shared_ptr<EncapFactory<time>> factory)
    {
      this->factory = factory;
    }

    template<typename time>
    shared_ptr<EncapFactory<time>> EncapSweeper<time>::get_factory() const
    {
      return factory;
    }

    template<typename time>
    shared_ptr<Encapsulation<time>> EncapSweeper<time>::get_end_state()
    {
      return this->end_state;
    }

    /**
     * @throws NotImplementedYet This function is required by EncapSweeper
     */
    template<typename time>
    void EncapSweeper<time>::advance()
    {
      throw NotImplementedYet("sweeper");
    }

    /**
     * @throws NotImplementedYet This function is required by EncapSweeper
     */
    template<typename time>
    void EncapSweeper<time>::reevaluate(bool initial_only)
    {
      UNUSED(initial_only);
      throw NotImplementedYet("sweeper");
    }

    /**
     * @throws NotImplementedYet This function is required by EncapSweeper
     */
    template<typename time>
    void EncapSweeper<time>::integrate(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const
    {
      UNUSED(dt); UNUSED(dst);
      throw NotImplementedYet("sweeper");
    }

    template<typename time>
    void EncapSweeper<time>::set_residual_tolerances(time abs_residual_tol, time rel_residual_tol,
                                                     int order)
    {
      this->abs_residual_tol = abs_residual_tol;
      this->rel_residual_tol = rel_residual_tol;
      this->residual_norm_order = order;
    }

    /**
     * @throws NotImplementedYet This function is required by EncapSweeper for residual computation
     */
    template<typename time>
    void EncapSweeper<time>::residual(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const
    {
      UNUSED(dt); UNUSED(dst);
      throw NotImplementedYet("residual");
    }

    /**
     * @internals
     * The maximum of the absolute and relative residuals at all time nodes is checked against the
     * predefined tolerances.
     *
     * In case both, EncapSweeper::abs_residual_tol and EncapSweeper::rel_residual_tol, are `0.0`,
     * this is a no-op.
     *
     * Otherwise, and if residuals were not yet initialized, repeating calls to
     * EncapFactory::create() are done.
     * @endinternals
     */
    template<typename time>
    bool EncapSweeper<time>::converged()
    {
      if (this->abs_residual_tol > 0.0 || this->rel_residual_tol > 0.0) {
        if (this->residuals.size() == 0) {
          for (size_t m = 0; m < this->get_nodes().size(); m++) {
            this->residuals.push_back(this->get_factory()->create(pfasst::encap::solution));
          }
        }
        this->residual(this->get_controller()->get_step_size(), this->residuals);
        vector<time> anorms, rnorms;
        for (size_t m = 0; m < this->get_nodes().size(); m++) {
          anorms.push_back(this->residuals[m]->norm0());
          rnorms.push_back(anorms.back() / this->get_state(m)->norm0());
        }
        auto amax = *std::max_element(anorms.begin(), anorms.end());
        auto rmax = *std::max_element(rnorms.begin(), rnorms.end());
        if (amax < this->abs_residual_tol || rmax < this->rel_residual_tol) {
          return true;
        }
      }
      return false;
    }

    template<typename time>
    void EncapSweeper<time>::post(ICommunicator* comm, int tag)
    {
      this->start_state->post(comm, tag);
    }

    template<typename time>
    void EncapSweeper<time>::send(ICommunicator* comm, int tag, bool blocking)
    {
      this->end_state->send(comm, tag, blocking);
    }

    template<typename time>
    void EncapSweeper<time>::recv(ICommunicator* comm, int tag, bool blocking)
    {
      this->start_state->recv(comm, tag, blocking);
      if (this->quadrature->left_is_node()) {
        this->state[0]->copy(this->start_state);
      }
    }

    template<typename time>
    void EncapSweeper<time>::broadcast(ICommunicator* comm)
    {
      if (comm->rank() == comm->size() - 1) {
        this->start_state->copy(this->end_state);
      }
      this->start_state->broadcast(comm);
    }


    template<typename time>
    EncapSweeper<time>& as_encap_sweeper(shared_ptr<ISweeper<time>> x)
    {
      shared_ptr<EncapSweeper<time>> y = dynamic_pointer_cast<EncapSweeper<time>>(x);
      assert(y);
      return *y.get();
    }

    template<typename time>
    const EncapSweeper<time>& as_encap_sweeper(shared_ptr<const ISweeper<time>> x)
    {
      shared_ptr<const EncapSweeper<time>> y = dynamic_pointer_cast<const EncapSweeper<time>>(x);
      assert(y);
      return *y.get();
    }

  }  // ::pfasst::encap
}  // ::pfasst
