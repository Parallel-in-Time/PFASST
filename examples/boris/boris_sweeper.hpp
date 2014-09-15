#ifndef _EXAMPLES__BORIS__BORIS_SWEEPER__HPP_
#define _EXAMPLES__BORIS__BORIS_SWEEPER__HPP_

#include <cstdlib>
#include <vector>
#include <map>
#include <cassert>

#include <Eigen/Dense>

#include <pfasst/encap/encap_sweeper.hpp>

#include "particle.hpp"
#include "particle_3d.hpp"
#include "physics.hpp"

using namespace std;
using namespace pfasst;
using namespace pfasst::encap;

template<typename coeff>
using Matrix = Eigen::Matrix<coeff, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

typedef map<pair<size_t, size_t>, double> error_map;


template<
  typename scalar,
  typename time = time_precision
>
class BorisSweeper
  : public EncapSweeper<time>
{
  public:
    typedef Particle3DEncapsulation<scalar, time> encap_type;
    typedef typename encap_type::position_type position_type;
    typedef typename encap_type::velocity_type velocity_type;
    typedef typename encap_type::acceleration_type acceleration_type;
    typedef Particle3DFactory<scalar, time> factory_type;
    typedef ElectricField<scalar, time, Particle3DEncapsulation> e_field_type;
    typedef MagneticField<scalar, time, Particle3DEncapsulation> b_field_type;

  private:
    shared_ptr<e_field_type> e_field;
    shared_ptr<b_field_type> b_field;
    error_map errors;

  protected:
    vector<shared_ptr<encap_type>> particles;
    vector<shared_ptr<encap_type>> previous_particles;
    vector<shared_ptr<encap_type>> tau_corrections;

    vector<acceleration_type> e_field_evals;

    vector<velocity_type> s_integrals;
    vector<position_type> ss_integrals;

    //! delta_nodes[m] = nodes[m] - nodes[m-1]
    vector<time> delta_nodes;

    Matrix<time> s_mat;
    Matrix<time> ss_mat;
    Matrix<time> qx_mat;
    Matrix<time> qt_mat;

  public:
    //! @{
    BorisSweeper()
      :   e_field(nullptr)
        , b_field(nullptr)
        , errors()
    {}

    BorisSweeper(const BorisSweeper<scalar, time>& other) = delete;
    BorisSweeper(BorisSweeper<scalar, time>&& other) = delete;

    virtual ~BorisSweeper()
    {}
    //! @}

    //! @{
    virtual void set_nodes(vector<time> nodes) override
    {
      EncapSweeper<time>::set_nodes(nodes);

      // compute delta nodes
      size_t nnodes = this->nodes.size();
      this->delta_nodes = vector<time>(nnodes, time(0.0));
      for (size_t m = 1; m < nnodes; m++) {
        delta_nodes[m] = nodes[m] - nodes[m - 1];
      }
    }

    virtual void set_state(shared_ptr<const Encapsulation<time>> u0, size_t m) override
    {
      shared_ptr<const encap_type> u0_cast = dynamic_pointer_cast<const encap_type>(u0);
      assert(u0_cast);
      this->set_state(u0_cast, m);
    }

    virtual void set_state(shared_ptr<const encap_type> u0, size_t m)
    {
      this->particles[m]->operator=(*(u0.get()));
    }

    virtual shared_ptr<Encapsulation<time>> get_state(size_t m) const override
    {
      return this->particles[m];
    }

    virtual shared_ptr<Encapsulation<time>> get_tau(size_t m) const override
    {
      return this->tau_corrections[m];
    }

    virtual shared_ptr<Encapsulation<time>> get_saved_state(size_t m) const override
    {
      return this->previous_particles[m];
    }
    //! @}

    //! @{
    virtual void set_e_field(shared_ptr<e_field_type> e_field)
    {
      this->e_field = e_field;
    }

    virtual shared_ptr<e_field_type> get_e_field()
    {
      return this->e_field;
    }

    virtual void set_b_field(shared_ptr<b_field_type> b_field)
    {
      this->b_field = const_pointer_cast<b_field_type>(b_field);
    }

    virtual shared_ptr<b_field_type> get_b_field()
    {
      return this->b_field;
    }
    //! @}

    //! @{
    virtual void exact(shared_ptr<Encapsulation<time>> q, time t)
    {
      shared_ptr<encap_type> q_cast = dynamic_pointer_cast<encap_type>(q);
      assert(q_cast);
      this->exact(q_cast, t);
    }

    virtual void exact(shared_ptr<encap_type> q, time t)
    {
      UNUSED(q); UNUSED(t);
      // TODO: implement exact solution for Boris
    }

    virtual void echo_error(time t, bool predict = false) const
    {
      UNUSED(t); UNUSED(predict);
      // TODO: implement echo_error
    }

    virtual error_map get_errors() const
    {
      return this->errors;
    }
    //! @}

    //! @{
    virtual void setup(bool coarse = false) override
    {
      vector<time> nodes = this->get_nodes();
      auto is_proper = this->get_is_proper();
      assert(nodes.size() >= 1);
      size_t nnodes = nodes.size();

      Matrix<time> q_mat = compute_quadrature(nodes, nodes, is_proper, QuadratureMatrix::Q);
      assert(q_mat.rows() == q_mat.cols());
      Matrix<time> qq_mat = Matrix<time>(nnodes, nnodes);
      qq_mat.fill(time(0.0));
      this->s_mat = Matrix<time>(nnodes, nnodes);
      this->s_mat.fill(time(0.0));
      this->qt_mat = Matrix<time>(nnodes, nnodes);
      this->qt_mat.fill(time(0.0));
      this->qx_mat = Matrix<time>(nnodes, nnodes);
      this->qx_mat.fill(time(0.0));
      this->ss_mat = Matrix<time>(nnodes, nnodes);
      this->ss_mat.fill(time(0.0));

      // compute QQ, S, SS, Q_T, Q_x
      // building rules for Q_E and Q_I:
      //  Q_E is striclty lower diagonal matrix with delta nodes of column index
      //  Q_I is lower diagonal matrix with first row and column all zero and delta nodes of 
      //      column index minus one
      // TODO: should use Boost's matrix product function here sometime
      for (size_t i = 0; i < nnodes; i++) {
        for (size_t j = 0; i < nnodes; j++) {
          // Q_T = 0.5 * (Q_E + Q_I)
          this->s_mat(i, j) = q_mat(i, j) - ((i > 0) ? q_mat(i - 1, j) : time(0.0));

          this->qt_mat(i, j) = time(0.5) * (((j < i) ? delta_nodes[j + 1] : time(0.0)) +
                                            ((j > 0 && j <= i) ? delta_nodes[j + 1] : time(0.0)));

          for (size_t k = 0; k < nnodes; k++) {
            qq_mat(i, j) += q_mat(i, k) * q_mat(k, j);

            // Q_x = Q_E * Q_T + 0.5 * (Q_E ∘ Q_E)
            //  (only first term, i.e. matrix product)
            this->qx_mat(i, j) += ((k < i) ? delta_nodes[k + 1] : time(0.0)) * this->qt_mat(k, j);
          }
          // continuation of Q_x
          //  (i.e. Hadamard product of Q_E)
          this->qx_mat(i, j) += 0.5 * ((j < i) ? delta_nodes[j + 1] * delta_nodes[j + 1] : time(0.0));

          this->ss_mat(i, j) = qq_mat(i, j) - ((i > 0) ? qq_mat(i - 1, j) : time(0.0));
        }
      }

      UNUSED(coarse);
      // TODO: implement setup for BorisSweeper
    }

    virtual void advance() override
    {
      this->set_state(const_pointer_cast<const encap_type>(this->particles.back()), 0);
    }

    virtual void evaluate(size_t m) override
    {
      UNUSED(m);
      // TODO: implement evaluate for BorisSweeper
    }

    virtual void predict(bool initial) override
    {
      UNUSED(initial);
      this->spread();
    }

    virtual void sweep() override
    {
      const auto   nodes  = this->get_nodes();
      const size_t nnodes = nodes.size();
      assert(nnodes >= 1);
      velocity_type c_k_term;

      // compute integrals
      for (size_t m = 1; m < nnodes; m++) {
        this->s_integrals[m].zero();
        this->ss_integrals[m].zero();
        for (size_t l = 1; l < nnodes; l++) {
          this->s_integrals[l] += this->previous_particles[l]
                                      ->accel().convert(dt<time>(this->s_mat(m, l)));
          this->ss_integrals[l] += this->previous_particles[l]
                                       ->accel().convert(dtdt<time>(this->ss_mat(m, l)));
        }
      }

      for (size_t m = 0; m < nnodes - 1; m++) {
        ::dt<time> ds(this->delta_nodes[m+1]);

        //// Update Position (explicit)
        //
        // x_{m+1}^{k+1} = x_{m}^{k}
        this->particles[m+1]->pos() = position_type(this->particles[m]->pos());
        //               + delta_nodes_{m} * v_{0}
        this->particles[m+1]->pos() += this->particles[0]->vel().convert(ds);
        //               + \sum_{l=1}^{m} s_{m+1,l}^{x} (f_{l}^{k+1} - f_{l}^{k})
        for (size_t l = 1; l < m; l++) {
          this->particles[m+1]->pos() += this->particles[l]
                                           ->accel().convert(dtdt<time>(this->qx_mat(m+1, l)));
          this->particles[m+1]->pos() -= this->previous_particles[l]
                                           ->accel().convert(dtdt<time>(this->qx_mat(m+1,l)));
        }
        //               + ss_integral[m]
        for (size_t l = 1; l < nnodes; l++) {
          this->particles[m+1]->pos() += this->ss_integrals[l];
        }

        // evaluate electric field with new position
        this->e_field_evals[m+1] = this->e_field->evaluate(this->particles, m+1, nodes[m]);

        //// Update Velocity (semi-implicit)
        //
        c_k_term.zero();  // reset
        //                 - delta_nodes_{m} / 2 * f_{m+1}^{k}
        c_k_term -= 0.5 * this->previous_particles[m+1]->accel().convert(ds);
        //                 - delta_nodes_{m} / 2 * f_{m}^{k}
        c_k_term -= 0.5 * this->previous_particles[m]->accel().convert(ds);
        //                 + s_integral[m]
        for (size_t l = 1; l < nnodes; l++) {
          c_k_term += this->s_integrals[l];
        }

        // doing Boris' magic
        this->boris_solve(nodes[m], nodes[m+1], ds, m, c_k_term);
      }
    }

    virtual void save(bool initial_only=false) override
    {
      if (initial_only) {
        this->previous_particles[0] = make_shared<encap_type>(this->particles[0]);
      } else {
        for (size_t m = 0; m < this->previous_particles.size(); m++) {
          this->previous_particles[m] = make_shared<encap_type>(this->particles[m]);
        }
      }
    }

    virtual void spread() override
    {
      for (size_t m = 1; m < this->particles.size(); ++m) {
        this->set_state(const_pointer_cast<const encap_type>(this->particles[0]), m);
      }
    }
    //! @}

    //! @{
    virtual void post(ICommunicator* comm, int tag) override
    {
      UNUSED(comm); UNUSED(tag);
      // TODO: implement BorisSweeper::post
    };

    virtual void send(ICommunicator* comm, int tag, bool blocking) override
    {
      UNUSED(comm); UNUSED(tag); UNUSED(blocking);
      // TODO: implement BorisSweeper::send
    }

    virtual void recv(ICommunicator* comm, int tag, bool blocking) override
    {
      UNUSED(comm); UNUSED(tag); UNUSED(blocking);
      // TODO: implement BorisSweeper::recv
    }

    virtual void broadcast(ICommunicator* comm) override
    {
      UNUSED(comm);
      // TODO: implement BorisSweeper::broadcast
    }
    //! @}

    //! @{
    virtual void boris_solve(const time tm, const time t_next, const ::dt<time> dt, const size_t m,
                             const velocity_type& c_k_term)
    {
      UNUSED(tm);
      velocity_type c_k_term_half = scalar(0.5) * c_k_term;
      ::dt<time> beta(this->particles[m]->charge() * dt.v 
                      / (scalar(2.0) * this->particles[m]->mass()));
      acceleration_type e_mean = (this->e_field_evals[m] - this->e_field_evals[m+1]) * scalar(0.5);

      // v^{-} = v^{k}
      velocity_type v_minus = this->particles[m]->vel();

      // first Boris' drift
      //         + \beta * E_{mean} + c^{k} / 2
      v_minus += e_mean.convert(beta) + c_k_term_half;

      // Boris' kick
      velocity_type boris_t = this->b_field->evaluate(this->particles, m, t_next).convert(beta);
      velocity_type boris_t_sqr = dot(boris_t, boris_t);
      velocity_type boris_s = velocity_type((scalar(2.0) * boris_t.u) / (scalar(1.0) + boris_t_sqr.u),
                                            (scalar(2.0) * boris_t.v) / (scalar(1.0) + boris_t_sqr.v),
                                            (scalar(2.0) * boris_t.w) / (scalar(1.0) + boris_t_sqr.w));
      velocity_type v_prime = v_minus + dot(boris_t, v_minus);

      // final Boris' drift
      velocity_type v_plus = v_minus + dot(v_prime, boris_s);

      this->particles[m+1]->vel() = v_plus + e_mean.convert(beta) + c_k_term_half;
    }
    //! @}
};

#endif  // _EXAMPLES__BORIS__BORIS_SWEEPER__HPP_
