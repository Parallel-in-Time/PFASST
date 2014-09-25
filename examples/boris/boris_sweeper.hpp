#ifndef _EXAMPLES__BORIS__BORIS_SWEEPER__HPP_
#define _EXAMPLES__BORIS__BORIS_SWEEPER__HPP_

#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <vector>
#include <map>

#include <Eigen/Core>

#include <pfasst.hpp>
#include <pfasst/encap/encap_sweeper.hpp>

#include "particle.hpp"
#include "particle_3d.hpp"
#include "physics.hpp"

using namespace std;
using namespace pfasst;
using namespace pfasst::encap;

template<typename coeff>
using Vector3d = Eigen::Array<coeff, 3, 1>;

template<typename coeff>
using Matrix3d = Eigen::Matrix<coeff, 3, 3, Eigen::RowMajor>;

template<typename coeff>
using Matrix = Eigen::Matrix<coeff, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

typedef pair<size_t, size_t> error_index;

template<typename scalar>
struct ParticleError
{
  scalar x;
  scalar y;
  scalar z;
  scalar u;
  scalar v;
  scalar w;
};

template<typename scalar>
using error_pair = pair<ParticleError<scalar>, scalar>;

template<typename scalar>
using error_map = map<error_index, error_pair<scalar>>;


template<typename scalar>
static Vector3d<scalar> dot(const Vector3d<scalar>& first, const Vector3d<scalar>& second)
{
  Vector3d<scalar> vec;
  vec[0] = first[1] * second[2] - first[2] * second[1];
  vec[1] = first[2] * second[0] - first[0] * second[2];
  vec[2] = first[0] * second[1] - first[1] * second[0];
  return vec;
}

template<typename scalar>
static scalar norm(const Vector3d<scalar>& vec)
{
  return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}


template<
  typename scalar,
  typename time,
  typename EnergyOperatorT
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
    typedef EnergyOperatorT energy_operator_type;
    typedef typename energy_operator_type::e_field_type e_field_type;
    typedef typename energy_operator_type::b_field_type b_field_type;

  private:
    energy_operator_type energy_op;
    scalar epsilon;
    error_map<scalar> errors;

  protected:
    vector<shared_ptr<encap_type>> particles;
    vector<shared_ptr<encap_type>> previous_particles;
    vector<shared_ptr<encap_type>> tau_corrections;

    vector<acceleration_type> energy_evals;
    size_t f_evals;

    vector<velocity_type> s_integrals;
    vector<position_type> ss_integrals;

    //! delta_nodes[m] = nodes[m] - nodes[m-1]
    vector<time> delta_nodes;

    Matrix<time> s_mat;
    Matrix<time> ss_mat;
    Matrix<time> sx_mat;
    Matrix<time> qx_mat;
    Matrix<time> qt_mat;

  public:
    //! @{
    BorisSweeper(const energy_operator_type& energy_operator, scalar epsilon = -1.0)
      :   energy_op(energy_operator)
        , epsilon(epsilon)
        , errors()
        , f_evals(0)
    {}

    BorisSweeper(const BorisSweeper<scalar, time, EnergyOperatorT>& other) = delete;
    BorisSweeper(BorisSweeper<scalar, time, EnergyOperatorT>&& other) = delete;

    virtual ~BorisSweeper()
    {
      cout << "f evals: " << this->f_evals << endl;
    }
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
    virtual const e_field_type& get_e_field() const
    {
      return this->energy_op.get_e_field();
    }

    virtual const b_field_type& get_b_field() const
    {
      return this->energy_op.get_b_field();
    }

    virtual void set_energy_operator(const energy_operator_type& e_operator)
    {
      this->energy_op = e_operator;
    }

    virtual const energy_operator_type& get_energy_operator() const
    {
      return this->energy_op;
    }
    //! @}

    //! @{
    virtual void exact(shared_ptr<Encapsulation<time>> q, time t)
    {
      shared_ptr<encap_type> q_cast = dynamic_pointer_cast<encap_type>(q);
      assert(q_cast);
      this->exact(q_cast, t);
    }

    virtual void exact(shared_ptr<encap_type> q, const time t)
    {
      this->exact(*(q.get()), t);
    }

    virtual void exact(encap_type& q, const time t)
    {
      typedef complex<scalar> C;
      C i(0.0, 1.0);
      cout << OUT::bold << OUT::red << "EXACT (t=" << t << "):" << OUT::reset << endl;
      auto initial = this->particles[0];
      scalar x0 = initial->pos().x,
             y0 = initial->pos().y,
             z0 = initial->pos().z,
             u0 = initial->vel().u,
             v0 = initial->vel().v,
             w0 = initial->vel().w,
             omega_e = this->get_e_field().omega_e,
             omega_b = this->get_b_field().omega_b;
      C c_epsilon = C(this->epsilon, 0.0);

      C omega_tilde = sqrt(-2.0 * c_epsilon) * omega_e;
//       cout << "  omega_tilde: " << omega_tilde << endl;
      q.pos().z = (z0 * cos(omega_tilde * (scalar)(t)) + w0 / omega_tilde * sin(omega_tilde * (scalar)(t))).real();

      C sqrt_in_omega = sqrt(pow(omega_b, 2) + 4.0 * epsilon * pow(omega_e, 2));
//       cout << "  sqrt in Omega: " << sqrt_in_omega << endl;
      C omega_minus = 0.5 * (omega_b - sqrt_in_omega);
//       cout << "  Omega-: " << omega_minus << endl;
      C omega_plus = 0.5 * (omega_b + sqrt_in_omega);
//       cout << "  Omega+: " << omega_plus << endl;

      C r_minus = (omega_plus * x0 - v0) / (omega_plus - omega_minus);
//       cout << "  R-: " << r_minus << endl;
      C r_plus = x0 - r_minus;
//       cout << "  R+: " << r_plus << endl;

      C i_minus = (omega_plus * y0 + u0) / (omega_plus - omega_minus);
//       cout << "  I-: " << i_minus << endl;
      C i_plus = y0 - i_minus;
//       cout << "  I+: " << i_plus << endl;

      C x_y_move = (r_plus + i * i_plus) * exp(- i * omega_plus * (scalar)(t))
                                               + (r_minus + i * i_minus) * exp(- i * omega_minus * (scalar)(t));
//       cout << "  w: " << x_y_move << endl;

      q.pos().x = x_y_move.real();
      q.pos().y = x_y_move.imag();

      q.vel().w = (- z0 * omega_tilde * sin(omega_tilde * (scalar)(t)) + w0 * cos(omega_tilde * (scalar)(t))).real();
      C u_v_move = (- i * omega_plus * (r_plus + i * i_plus)) * exp(-i * omega_plus * (scalar)(t))
                                 - (i * omega_minus * (r_minus + i * i_minus)) * exp(-i * omega_minus * (scalar)(t));
//       cout << "  w/dt: " << u_v_move << endl;
      q.vel().u = u_v_move.real();
      q.vel().v = u_v_move.imag();
      cout << "    " << q << endl;;
    }

    virtual void echo_error(const time t, bool predict = false)
    {
      auto end = this->particles.back();
      shared_ptr<encap_type> ex = make_shared<encap_type>();
      this->exact(ex, t);
      scalar e_ex = this->energy_op.evaluate({ex}, t);
      scalar e_end = this->energy_op.evaluate({this->particles.back()}, t);

      scalar e_err = e_ex - e_end;
      ParticleError<scalar> p_err;
      p_err.x = ex->pos().x - end->pos().x;
      p_err.y = ex->pos().y - end->pos().y;
      p_err.z = ex->pos().z - end->pos().z;
      p_err.u = ex->vel().u - end->vel().u;
      p_err.v = ex->vel().v - end->vel().v;
      p_err.w = ex->vel().w - end->vel().w;
      error_pair<scalar> e_pair(p_err, e_err);

      size_t n = this->get_controller()->get_step();
      size_t k = this->get_controller()->get_iteration();
      error_index nk(n, k);

      cout << "err: " << n << " " << k << " " << scientific
           << e_err << fixed
           << " (e: " << e_end
           << "; pos: " << end->pos()
           << "; vel: " << end->pos()
           << scientific
           << "; xe: " << p_err.x
           << "; ye: " << p_err.y
           << "; ze: " << p_err.z
           << "; ue: " << p_err.u
           << "; ve: " << p_err.v
           << "; we: " << p_err.w
           << ", " << boolalpha << predict << fixed << ")" << endl;

      this->errors.insert(pair<error_index, error_pair<scalar>>(nk, e_pair));
    }

    virtual error_map<scalar> get_errors() const
    {
      return this->errors;
    }
    //! @}

    //! @{
    virtual void setup(bool coarse = false) override
    {
      cout << OUT::bold << OUT::red << "## SETUP" << OUT::reset << endl;
      auto nodes = this->get_nodes();
      auto is_proper = this->get_is_proper();
      assert(nodes.size() >= 1);
      const size_t nnodes = nodes.size();

      for (size_t m = 0; m < nnodes; ++m) {
        this->particles.push_back(dynamic_pointer_cast<encap_type>(this->get_factory()->create(pfasst::encap::solution)));
        this->previous_particles.push_back(dynamic_pointer_cast<encap_type>(this->get_factory()->create(pfasst::encap::solution)));
        this->s_integrals.push_back(velocity_type());
        this->ss_integrals.push_back(position_type());
        this->energy_evals.push_back(acceleration_type());
        if (coarse) {
          this->tau_corrections.push_back(dynamic_pointer_cast<encap_type>(this->get_factory()->create(pfasst::encap::solution)));
        }
      }

      cout << "Delta Nodes:  [";
      for(auto n : delta_nodes) { cout << "  " << n; }
      cout << "  ]" << endl;

      auto q_mat_trimmed = compute_quadrature(nodes, nodes, is_proper, QuadratureMatrix::Q);
      decltype(q_mat_trimmed) q_mat(nnodes, nnodes);
      q_mat.fill(time(0.0));
      q_mat.block(1, 0, nnodes - 1, nnodes) = q_mat_trimmed;
      cout << "Q:" /*<< cout.precision(16)*/ << endl << q_mat << endl;
      auto qq_mat = q_mat * q_mat;
      cout << "QQ:" << endl << qq_mat << endl;
      this->s_mat = Matrix<time>(nnodes, nnodes);
      this->s_mat.fill(time(0.0));
      this->qt_mat = Matrix<time>(nnodes, nnodes);
      this->qt_mat.fill(time(0.0));
      this->qx_mat = Matrix<time>(nnodes, nnodes);
      this->qx_mat.fill(time(0.0));
      this->ss_mat = Matrix<time>(nnodes, nnodes);
      this->ss_mat.fill(time(0.0));
      this->sx_mat = Matrix<time>(nnodes, nnodes);
      this->sx_mat.fill(time(0.0));

      // building rules for Q_E and Q_I:
      //  Q_E is striclty lower diagonal matrix with delta nodes of column index
      //  Q_I is lower diagonal matrix with first row and column all zero and delta nodes of 
      //      column index minus one
      Matrix<time> qe_mat = Matrix<time>(nnodes, nnodes);
      qe_mat.fill(time(0.0));
      Matrix<time> qi_mat = Matrix<time>(nnodes, nnodes);
      qi_mat.fill(time(0.0));
      for (size_t i = 0; i < nnodes; ++i) {
        for (size_t j = 0; j < nnodes; ++j) {
          if (j < i)           { qe_mat(i, j) = delta_nodes[j + 1]; }
          if (j > 0 && j <= i) { qi_mat(i, j) = delta_nodes[j]; }
        }
      }
      cout << "Q_E:" << endl << qe_mat << endl;
      cout << "Q_I:" << endl << qi_mat << endl;

      // Q_T = 0.5 * (Q_E + Q_I)
      qt_mat = 0.5 * (qe_mat + qi_mat);
      cout << "Q_T:" << endl << qt_mat << endl;

      // Q_x = Q_E * Q_T + 0.5 * (Q_E ∘ Q_E)
      //  (only first term, i.e. matrix product)
      qx_mat = qe_mat * qt_mat;
      cout << "Q_x1 = " << endl << qx_mat << endl;

      // compute QQ, S, SS, Q_T, Q_x
      for (size_t i = 1; i < nnodes; i++) {
        for (size_t j = 0; j < nnodes; j++) {
          // continuation of Q_x
          //  (i.e. Hadamard product of Q_E)
          this->qx_mat(i, j) += 0.5 * qe_mat(i, j) * qe_mat(i, j);
        }

        this->s_mat.row(i) = q_mat.row(i) - q_mat.row(i - 1);
        this->ss_mat.row(i) = qq_mat.row(i) - qq_mat.row(i - 1);
        this->sx_mat.row(i) = qx_mat.row(i) - qx_mat.row(i - 1);
      }
      cout << "Q_x:" << endl << qx_mat << endl;
      cout << "S:" << endl << s_mat << endl;
      cout << "SS:" << endl << ss_mat << endl;
      cout << "S_x:" << endl << sx_mat << endl;

      UNUSED(coarse);
      // TODO: implement setup for BorisSweeper
      cout << OUT::yellow << "setup DONE" << OUT::reset << endl;
    }

    virtual void advance() override
    {
      cout << OUT::bold << OUT::red << "ADVAAAAAANCE" << OUT::reset << endl;
      this->set_state(const_pointer_cast<const encap_type>(this->particles.back()), 0);
      this->energy_evals.front() = this->energy_evals.back();
      cout << OUT::yellow << "advaaaaaance DONE" << OUT::reset << endl;
    }

    virtual void evaluate(size_t m) override
    {
      cout << OUT::bold << OUT::red << "EVALUATE" << OUT::reset << endl;
//       Vector3d<scalar> vel; vel.fill(scalar(0.0));
//       Vector3d<scalar> B; B.fill(scalar(0.0));
//       Vector3d<scalar> b; b.fill(scalar(0.0));
      Eigen::Matrix<scalar, 3, 6> mat;
      mat << - this->epsilon * pow(this->get_e_field().omega_e, 2), 0, 0, 0, this->get_b_field().omega_b, 0,
             0, - this->epsilon * pow(this->get_e_field().omega_e, 2), 0, - this->get_b_field().omega_b, 0, 0,
             0, 0, 2 * this->epsilon * pow(this->get_e_field().omega_e, 2), 0, 0, 0;
      Eigen::Matrix<scalar, 6, 1> u_vec;
      u_vec.block(0, 0, 3, 1) = this->particles[m]->pos().as_matrix().transpose();
      u_vec.block(3, 0, 3, 1) = this->particles[m]->vel().as_matrix().transpose();

      Eigen::Matrix<scalar, 3, 1> f; f.fill(scalar(0.0));
      f = mat * u_vec;
//       cout << "  mat * u_vec = f(x,v):" << endl
//            << mat << endl << " * " << endl
//            << u_vec << endl << " = " << endl
//            << f << endl;
      this->particles[m]->accel() = acceleration_type(f);
      time t = this->get_controller()->get_time() + this->get_controller()->get_time_step() * this->delta_nodes[m];
//       auto e = this->get_e_field().evaluate(this->particles, m, t);
//       vel = this->particles[m]->vel().as_matrix().transpose().array();
//       B = this->get_b_field().get_field_vector().transpose().array() / this->particles[m]->alpha();
//       b = dot(vel, B);
//       cout << OUT::red << "  B-Field" << OUT::reset << endl
//            << "    b: " << b.transpose() << " = " << vel.transpose() << " . (" << this->get_b_field().get_field_vector() << " / " << this->particles[m]->alpha() << " [=" << B.transpose() << "] )" << endl
//            << OUT::yellow << "  b-field DONE" << OUT::reset << endl;
//       this->particles[m]->accel() = this->particles[m]->alpha() * (e + acceleration_type(b));
      this->f_evals++;
      cout << "f(x,v,t=" << t << ",m=" << m << "): " << this->particles[m]->accel().as_matrix() << endl;
      cout << OUT::yellow << "evaluate DONE" << OUT::reset << endl;
    }

    virtual void predict(bool initial) override
    {
      UNUSED(initial);
      cout << OUT::bold << OUT::red << "PREEEEEEDICT(" << boolalpha << initial << fixed << ")" << OUT::reset << endl;
      this->spread();
      for (size_t m = 0; m < this->particles.size(); ++m) {
        this->evaluate(m);
      }
      this->save();
      cout << OUT::yellow << "preeeeeedict DONE" << OUT::reset << endl;
    }

    virtual void sweep() override
    {
      cout << OUT::bold << OUT::red << "SWEEEEEEEEEEEEEEEEEEEEEEP:" << OUT::reset << endl;

      cout << "  previous:" << endl;
      for(auto p : this->previous_particles) { cout << "    " << *(p.get()) << endl; }
      cout << "  current:" << endl;
      for(auto p : this->particles) { cout << "    " << *(p.get()) << endl; }

      const auto   nodes  = this->get_nodes();
      const size_t nnodes = nodes.size();
      assert(nnodes >= 1);
      time t  = this->get_controller()->get_time();
      time dt = this->get_controller()->get_time_step();
      {
        scalar e = this->energy_op.evaluate({this->particles.front()}, t);
        cout << OUT::blue << "  !! initial total Energy: " << OUT::bold << e << OUT::reset << endl;
      }
      velocity_type c_k_term;

      // compute integrals
      cout << "  integrals:" << scientific << endl;
      for(size_t i = 0; i < nnodes; ++i) {
        this->s_integrals[i].zero();
        this->ss_integrals[i].zero();
      }
      // starting at m=1 as m=0 will only add zeros
      for (size_t m = 1; m < nnodes; m++) {
        cout << "   m=" << m << endl;
        for (size_t l = 0; l < nnodes; l++) {
          cout << "    l=" << l << endl;
          this->s_integrals[l] += this->previous_particles[l]
                                      ->accel().convert(::dt<time>(dt * this->s_mat(m, l)));
          cout << "     s["<<l<<"] += [" << this->previous_particles[l]->accel() << "] * " << dt << " * " << this->s_mat(m, l) << endl;
          this->ss_integrals[l] += this->previous_particles[l]
                                       ->accel().convert(::dtdt<time>(dt * dt * this->ss_mat(m, l)));
          cout << "     ss["<<l<<"] += [" << this->previous_particles[l]->accel() << "] * " << dt * dt << " * " << this->ss_mat(m, l) << endl;
        }
      }
      cout << "   ==>" << endl << "   s:  [";
      for(auto i : this->s_integrals) { cout << "  " << i; }
      cout << "  ]" << endl
           << "   ss: [";
      for(auto i : this->ss_integrals) { cout << "  " << i; }
      cout << "  ]" << fixed << endl;

      cout << "   p(0): " << *(this->particles[0]) << endl;
      this->evaluate(0);
//       this->particles.front()->accel() = this->get_e_field().evaluate(this->particles, 0, nodes[0]);
      cout << "   E(0): " << this->particles.front()->accel() << endl;

      for (size_t m = 0; m < nnodes - 1; m++) {
        ::dt<time> ds(dt * this->delta_nodes[m+1]);
        cout << "   m=" << m << " to m=" << (m+1) << "; ds=" << ds.v << endl;

        cout << "    old pos:   " << this->particles[m]->pos() << endl ;
        //// Update Position (explicit)
        //
        // x_{m+1}^{k+1} = x_{m}^{k}
        this->particles[m+1]->pos() = position_type(this->particles[m]->pos());
        //               + delta_nodes_{m+1} * v_{0}
        this->particles[m+1]->pos() += this->particles[0]->vel().convert(ds);
        cout << "             + " << this->particles[0]->vel().convert(ds) << endl;
        //               + \sum_{l=1}^{m} sx_{m+1,l}^{x} (f_{l}^{k+1} - f_{l}^{k})
        for (size_t l = 0; l <= m; l++) {
          this->particles[m+1]->pos() += this->particles[l]
                                             ->accel().convert(dtdt<time>(ds.v * ds.v * this->sx_mat(m+1, l)));
          cout << "             + ([" << this->particles[l]->accel() << "] * " << ds.v * ds.v << " * " << this->sx_mat(m+1, l) << " [=" << this->particles[l]
                                             ->accel().convert(dtdt<time>(ds.v * ds.v * this->sx_mat(m+1, l))) << "])" << endl;
//           cout << "            [=" << this->particles[m+1]->pos() << "]" << endl;
          this->particles[m+1]->pos() -= this->previous_particles[l]
                                             ->accel().convert(dtdt<time>(ds.v * ds.v * this->sx_mat(m+1, l)));
          cout << "             - ([" << this->previous_particles[l]->accel() << "] * " << ds.v * ds.v << " * " << this->sx_mat(m+1, l) << " [=" << this->previous_particles[l]
                                             ->accel().convert(dtdt<time>(ds.v * ds.v * this->sx_mat(m+1, l))) << "])" << endl;
        }
//         cout << "            [=" << this->particles[m+1]->pos() << "]" << endl;
        //               + ss_integral[m]
        for (size_t l = 0; l < nnodes; l++) {
          cout << "             + " << this->ss_integrals[l] << endl;
          this->particles[m+1]->pos() += this->ss_integrals[l];
        }
        cout << "    new pos: => " << this->particles[m+1]->pos() << endl;

        // evaluate electric field with new position
        this->particles[m+1]->accel() = this->get_e_field().evaluate(this->particles, m+1, t + nodes[m]);
        cout << "    E(m+1):  " << this->particles[m+1]->accel() << endl;

        cout << "    old vel: " << this->particles[m+1]->vel() << endl;
        //// Update Velocity (semi-implicit)
        //
        c_k_term.zero();  // reset
        cout << "    c_k_term:" << endl;
        //                 - delta_nodes_{m} / 2 * f_{m+1}^{k}
        c_k_term -= 0.5 * this->previous_particles[m+1]->accel().convert(::dt<time>(1.0));
        cout << "              - 0.5 * " << this->previous_particles[m+1]->accel() << endl;
        //                 - delta_nodes_{m} / 2 * f_{m}^{k}
        c_k_term -= 0.5 * this->previous_particles[m]->accel().convert(::dt<time>(1.0));
        cout << "              - 0.5 * " << this->previous_particles[m]->accel() << endl;
        //                 + s_integral[m]
        for (size_t l = 0; l < nnodes; l++) {
          c_k_term += this->s_integrals[l] / ds.v;
          cout << "              + " << this->s_integrals[l] << " / " << ds.v << endl;
        }
        cout << "              => " << c_k_term << endl;

        // doing Boris' magic
        this->boris_solve(nodes[m], nodes[m+1], ds, m, c_k_term);

        cout << "    new vel: " << this->particles[m+1]->vel() << endl;
      }
      this->save();
      this->echo_error(t + dt);
      cout << OUT::yellow << "sweeeeeeeeeeeeeeeeeeeeeep DONE" << OUT::reset << endl;
    }

    virtual void save(bool initial_only=false) override
    {
      cout << OUT::bold << OUT::red << "SAVE!!!" << OUT::reset << endl;
      if (initial_only) {
        this->previous_particles[0] = make_shared<encap_type>(*(this->particles[0].get()));
      } else {
        for (size_t m = 0; m < this->previous_particles.size(); m++) {
          cout << "  m=" << m;
          this->previous_particles[m] = make_shared<encap_type>(*(this->particles[m].get()));
        }
        cout << endl;
      }
      cout << OUT::yellow << "save!!! DONE" << OUT::reset << endl;
    }

    virtual void spread() override
    {
      cout << OUT::bold << OUT::red << "SPREAAAAAD" << OUT::reset << endl;
      for (size_t m = 1; m < this->particles.size(); ++m) {
        this->set_state(const_pointer_cast<const encap_type>(this->particles[0]), m);
      }
      cout << OUT::yellow << "spreaaaaad DONE" << OUT::reset << endl;
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
    virtual void boris_solve(const time tm, const time t_next, const ::dt<time> ds, const size_t m,
                             const velocity_type& c_k_term)
    {
      cout << OUT::bold << OUT::red << "BORIS(tm=" << tm << ", dt=" << ds.v << "):" << OUT::reset << endl;
      UNUSED(t_next);
      Vector3d<scalar> c_k_term_half; c_k_term_half.fill(scalar(0.0));
      c_k_term_half = c_k_term.as_matrix().transpose().array() / scalar(2.0);
      ::dt<time> beta(this->particles[m]->alpha() / time(2.0) * ds.v);
//       cout << "  beta:   " << beta.v << endl;
      Vector3d<scalar> e_mean; e_mean.fill(scalar(0.0));
      e_mean = (this->get_e_field().evaluate(this->particles, m, tm).as_matrix().transpose().array() + this->particles[m+1]->accel().as_matrix().transpose().array()) * scalar(0.5);
//       cout << "  E_mean: " << e_mean.transpose() << endl;

      // first Boris' drift
      //   v^{-} = v^{k}
      Vector3d<scalar> v_minus; v_minus.fill(scalar(0.0));
      v_minus = this->particles[m]->vel().as_matrix().transpose();
      //           + \beta * E_{mean} + c^{k} / 2
      v_minus += e_mean * beta.v + c_k_term_half;
//       cout << "  drift:  " << v_minus.transpose() << " = " << this->particles[m]->vel() << " + " << e_mean.transpose() * beta.v << " + " << c_k_term_half.transpose() << endl;

      // Boris' kick
      Vector3d<scalar> boris_t; boris_t.fill(scalar(0.0));
      boris_t = beta.v * this->get_b_field().get_field_vector().transpose();
//       cout << "       t: " << boris_t.transpose() << " = " << beta.v << " * " << this->get_b_field().get_field_vector() << endl;
      Vector3d<scalar> v_prime; v_prime.fill(scalar(0.0));
      v_prime = v_minus + dot(v_minus, boris_t);
//       cout << "  kick:   " << v_prime.transpose() << " = " << v_minus.transpose() << " + " << v_minus.transpose() << " . " << boris_t.transpose() << " [= " << dot(v_minus, boris_t).transpose() << "]" << endl;

      // final Boris' drift
      scalar boris_t_sqr = pow(norm(boris_t), 2);
//       cout << "   |t|^2: " << boris_t_sqr << endl;
      Vector3d<scalar> boris_s; boris_s.fill(scalar(0.0));
      boris_s = (scalar(2.0) * boris_t) / (boris_t_sqr + scalar(1.0));
//       cout << "       s: " << boris_s.transpose() << " = (2.0 * " << boris_t.transpose() << ") / (" << boris_t_sqr << " + 1.0)" << endl;
      Vector3d<scalar> v_plus; v_plus.fill(scalar(0.0));
      v_plus = v_minus + dot(v_prime, boris_s);
//       cout << "  drift:  " << v_plus.transpose() << " = " << v_minus.transpose() << " + " << v_prime.transpose() << " . " << boris_s.transpose() << " [= " << dot(boris_s, v_prime).transpose() << "]" << endl;

      cout << "   norm difference: " << scientific << abs(norm(v_minus) - norm(v_plus)) << fixed << endl;
      assert(abs(norm(v_minus) - norm(v_plus)) <= 10e-8);

      this->particles[m+1]->vel() = velocity_type(v_plus + e_mean * beta.v + c_k_term_half);
//       cout << " new_vel: " << this->particles[m+1]->vel() << " = " << v_plus.transpose() << " + " << beta.v << " * " << e_mean.transpose() << " + " << c_k_term_half.transpose() << endl;
      cout << OUT::yellow << "boris DONE" << OUT::reset << endl;
    }
    //! @}
};

#endif  // _EXAMPLES__BORIS__BORIS_SWEEPER__HPP_
