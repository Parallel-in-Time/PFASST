#include <cmath>
#include <memory>
#include <iostream>
#include <iomanip>

#include <mpi.h>

#define WITH_MPI

#include <pfasst.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/controller/sdc.hpp>
#include <pfasst/controller/pfasst.hpp>
#include <pfasst/mpi_communicator.hpp>
#include <pfasst/encap/vector.hpp>
#include <pfasst/encap/implicit_sweeper.hpp>
#include <pfasst/encap/imex_sweeper.hpp>
#include <pfasst/encap/poly_interp.hpp>

extern "C" {
  double exact(double t, double x);
  double forcing(double t, double x);
  void take_step(double * values,
                 int      size,
                 double   t,
                 double   xstart,
                 double   deltaX,
                 double   deltaT,
                 double * matrix,
                 double * temp);
  void matvec_tridiag(double *x, double *g, int N, double* matrix);
  void compute_stencil(double   deltaX, double   deltaT, double * matrix);

  void interpolate_1D(const double * cvalues, double * fvalues, int      csize, int      fsize);
  void coarsen_1D(double * cvalues, const double * fvalues, int      csize, int      fsize);

}

using namespace std;
using pfasst::encap::Encapsulation;

template<typename time = pfasst::time_precision>
class ImplicitHeatSweeper
  : public pfasst::encap::IMEXSweeper<time>
//  : public pfasst::encap::ImplicitSweeper<time>
{
    double xstart = 0.0;
    double xstop = M_PI;

  public:
    // ImplicitHeatSweeper(int nspace)
    // {
    //   this->nspace = nspace;
    // }

    void exact(pfasst::encap::VectorEncapsulation<double>& u, time t)
    {
      auto dx = (this->xstop - this->xstart) / (u.size() - 1.0);
      for (int i=0; i<u.size(); i++) {
        u[i] = sin(i * dx)*cos(t);
      }
    }

    void exact(shared_ptr<Encapsulation<time>> u_encap, time t)
    {
      auto& u = pfasst::encap::as_vector<double, time>(u_encap);
      this->exact(u, t);
    }

    void echo_error(time t)
    {
      auto& qend = pfasst::encap::as_vector<double, time>(this->get_end_state());
      pfasst::encap::VectorEncapsulation<double> qex(qend.size());

      this->exact(qex, t);

      double max = 0.0;
      for (size_t i = 0; i < qend.size(); i++) {
        double d = abs(qend[i] - qex[i]);
      if (d > max) { max = d; }
      }

      auto n = this->get_controller()->get_step();
      auto k = this->get_controller()->get_iteration();


      ML_CLOG(INFO, "User", "step: " << n << " iter: " << k << " err: " << std::scientific << max);
    }

    void echo_residual()
    {
      vector<shared_ptr<Encapsulation<time>>> residuals;

      for (size_t m = 0; m < this->get_nodes().size(); m++) {
        residuals.push_back(this->get_factory()->create(pfasst::encap::solution));
      }
      this->residual(this->get_controller()->get_step_size(), residuals);

      vector<time> rnorms;
      for (auto r: residuals) {
        rnorms.push_back(r->norm0());
      }
      auto rmax = *std::max_element(rnorms.begin(), rnorms.end());

      auto n = this->get_controller()->get_step();
      auto k = this->get_controller()->get_iteration();

      ML_CLOG(INFO, "User", "step: " << n << " iter: " << k << " res: " << rmax);
    }

    void post_sweep() override
    {
      time t  = this->get_controller()->get_time();
      time dt = this->get_controller()->get_step_size();
      this->echo_error(t+dt);
      //this->echo_residual();
    }

    void f_impl_eval(shared_ptr<Encapsulation<time>> f_impl_encap,
                     shared_ptr<Encapsulation<time>> u_encap,
                     time t) override
    {
      auto& u = pfasst::encap::as_vector<double, time>(u_encap);
      auto& f_impl = pfasst::encap::as_vector<double, time>(f_impl_encap);

      double dx = (xstop - xstart) / (u.size() - 1.0);
      double dx2 = dx*dx;
      // double matrix[3];

      // compute_stencil(dx, 1.0, matrix);
      // matrix[1] -= 1.0;
      // matvec_tridiag(u.data(), f_impl.data(), f_impl.size(), matrix);

      for (int i=1; i<f_impl.size()-1; i++) {
        f_impl[i] = (u[i-1] - 2*u[i] + u[i+1])/dx2 + forcing(t, i*dx);
      }
      f_impl[0] = 0.0; f_impl[f_impl.size()-1] = 0.0;

      // cout << " YAR!!! " << endl;
      // for (int i=0; i<f_impl.size(); i++) {
      //   cout << i << " " << u[i] << " " << f_impl[i] << endl;
      // }

    }

    /**
     * @copybrief pfasst::encap::IMEXSweeper::impl_solve()
     */
    void impl_solve(shared_ptr<Encapsulation<time>> f_impl_encap,
                    shared_ptr<Encapsulation<time>> u_encap,
                    time t, time dt,
                    shared_ptr<Encapsulation<time>> rhs_encap) override
    {
      UNUSED(t);
      auto& u = pfasst::encap::as_vector<double, time>(u_encap);
      auto& f_impl = pfasst::encap::as_vector<double, time>(f_impl_encap);
      auto& rhs = pfasst::encap::as_vector<double, time>(rhs_encap);

      for (int i=0; i<u.size(); i++) {
        u[i] = rhs[i];
      }

      double dx = (xstop - xstart) / (u.size() - 1.0);
      double matrix[3];
      double temp[u.size()];

      take_step(u.data(), u.size(), t+dt, 0.0, dx, dt, matrix, temp);

      for (int i=0; i<u.size(); i++) {
        f_impl[i] = (u[i] - rhs[i]) / dt;
      }

      // double dx2 = dx*dx;
      // for (int i=1; i<f_impl.size()-1; i++) {
      //   f_impl[i] = (u[i-1] - 2*u[i] + u[i+1])/dx2 + forcing(t, i*dx);
      // }
      // f_impl[0] = 0.0; f_impl[f_impl.size()-1] = 0.0;
      // pfasst::encap::VectorEncapsulation<double> qex(u.size());
      // this->exact(qex, t+1.5*dt);

      // for (int i=0; i<f_impl.size(); i++) {
      //   cout << dt << " " << i << " " << u[i] << " " << f_impl[i] << " " << rhs[i] << " " << qex[i] << " " << qex[i]/u[i] << endl;
      // }
    }

    void f_expl_eval(shared_ptr<Encapsulation<time>> f_expl_encap,
                     shared_ptr<Encapsulation<time>> u_encap,
                     time t) //override
    {
      UNUSED(u_encap); UNUSED(t);
      auto& f_expl = pfasst::encap::as_vector<double, time>(f_expl_encap);
      for (int i=0; i<f_expl.size(); i++) {
        f_expl[i] = 0.0;
      }
    }

};


template<typename time = pfasst::time_precision>
class BilinearTransfer1D
  : public pfasst::encap::PolyInterpMixin<time>
{
    using Encapsulation = pfasst::encap::Encapsulation<double>;

  public:
    void interpolate(shared_ptr<Encapsulation> dst,
                     shared_ptr<const Encapsulation> src) override
    {
      auto& fine = pfasst::encap::as_vector<double, time>(dst);
      auto& crse = pfasst::encap::as_vector<double, time>(src);
      interpolate_1D(crse.data(), fine.data(), crse.size(), fine.size());
    }

    void restrict(shared_ptr<Encapsulation> dst,
                  shared_ptr<const Encapsulation> src) override
    {
      auto& fine = pfasst::encap::as_vector<double, time>(src);
      auto& crse = pfasst::encap::as_vector<double, time>(dst);
      coarsen_1D(crse.data(), fine.data(), crse.size(), fine.size());
    }
};



int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  pfasst::init(argc, argv);

  pfasst::mpi::MPICommunicator comm(MPI_COMM_WORLD);
  pfasst::PFASST<> pf;

  auto const quad_type = pfasst::quadrature::QuadratureType::GaussLobatto;

  auto transfer = make_shared<BilinearTransfer1D<>>();

  auto const nlevels   = pfasst::config::get_value<int>("nlevels", 1);
  auto const nnodes    = pfasst::config::get_value<int>("nnodes", 3);
  auto const nspace    = pfasst::config::get_value<int>("nspace", 8193);
  auto const nsteps    = pfasst::config::get_value<int>("nsteps", 16);
  auto const niters    = pfasst::config::get_value<int>("niters", 4);
  auto const dt        = pfasst::config::get_value<double>("dt", 0.1);

  auto quad    = pfasst::quadrature::quadrature_factory(nnodes, quad_type);
  auto factory = make_shared<pfasst::encap::VectorFactory<double>>(nspace);
  auto sweeper = make_shared<ImplicitHeatSweeper<double>>();

  sweeper->set_quadrature(quad);
  sweeper->set_factory(factory);

  pf.set_comm(&comm);
  pf.add_level(sweeper, transfer);

  if (nlevels > 1) {
    auto quad    = pfasst::quadrature::quadrature_factory((nnodes-1)/2+1, quad_type);
    auto factory = make_shared<pfasst::encap::VectorFactory<double>>((nspace-1)/2+1);
    auto sweeper = make_shared<ImplicitHeatSweeper<double>>();
    sweeper->set_quadrature(quad);
    sweeper->set_factory(factory);
    pf.add_level(sweeper, transfer);
  }

  pf.set_duration(0.0, nsteps*dt, dt, niters);
  pf.setup();

  auto q0 = sweeper->get_start_state();
  sweeper->exact(q0, 0.0);

  pf.run();

  return 0;
}
