#include "pfasst/controller/status.hpp"

#include <sstream>
using namespace std;


namespace pfasst
{
  template<typename precision>
  Status<precision>::Status()
  {}

  template<typename precision>
  Status<precision>::Detail::Detail()
    :   step(0)
      , iteration(0)
      , time(0.0)
      , dt(0.0)
      , state(pfasst::State::UNKNOWN)
      , residual(0.0)
  {}

  template<typename precision>
  size_t Status<precision>::get_step() const
  {
    return this->_detail.step;
  }

  template<typename precision>
  size_t& Status<precision>::step()
  {
    return this->_detail.step;
  }

  template<typename precision>
  size_t Status<precision>::get_iteration() const
  {
    return this->_detail.iteration;
  }

  template<typename precision>
  size_t& Status<precision>::iteration()
  {
    return this->_detail.iteration;
  }

  template<typename precision>
  size_t Status<precision>::get_max_iterations() const
  {
    return this->_detail.max_iterations;
  }

  template<typename precision>
  size_t& Status<precision>::max_iterations()
  {
    return this->_detail.max_iterations;
  }

  template<typename precision>
  precision Status<precision>::get_time() const
  {
    return this->_detail.time;
  }

  template<typename precision>
  precision& Status<precision>::time()
  {
    return this->_detail.time;
  }

  template<typename precision>
  precision Status<precision>::get_dt() const
  {
    return this->_detail.dt;
  }

  template<typename precision>
  precision& Status<precision>::dt()
  {
    return this->_detail.dt;
  }

  template<typename precision>
  precision Status<precision>::get_t_end() const
  {
    return this->_detail.t_end;
  }

  template<typename precision>
  precision& Status<precision>::t_end()
  {
    return this->_detail.t_end;
  }

  template<typename precision>
  State Status<precision>::get_state() const
  {
    return this->_detail.state;
  }

  template<typename precision>
  State& Status<precision>::state()
  {
    return this->_detail.state;
  }

  template<typename precision>
  precision Status<precision>::get_residual() const
  {
    return this->_detail.residual;
  }

  template<typename precision>
  precision& Status<precision>::residual()
  {
    return this->_detail.residual;
  }

  template<typename precision>
  void Status<precision>::send(shared_ptr<comm::Communicator> comm,
                               const int dest_rank, const int tag, const bool blocking)
  {
    if (blocking) {
      comm->send(&(this->_detail.residual), 1, dest_rank, tag);
      comm->send((int*)&(this->_detail.state), 1, dest_rank, tag + 1); // TODO: tag computation
    } else {
      comm->isend(&(this->_detail.residual), 1, dest_rank, tag);
      comm->isend((int*)&(this->_detail.state), 1, dest_rank, tag + 1); // TODO: tag computation
    }
  }

  template<typename precision>
  void Status<precision>::recv(shared_ptr<comm::Communicator> comm,
                               const int src_rank, const int tag, const bool blocking)
  {
    if (blocking) {
      comm->recv((int*)&(this->_detail.state), 1, src_rank, tag + 1); // TODO: tag computation
      comm->recv(&(this->_detail.residual), 1, src_rank, tag);
    } else {
      comm->irecv((int*)&(this->_detail.state), 1, src_rank, tag + 1); // TODO: tag computation
      comm->irecv(&(this->_detail.residual), 1, src_rank, tag);
    }
  }

  template<typename precision>
  void Status<precision>::bcast(shared_ptr<comm::Communicator> comm,
                                const int root_rank)
  {
    comm->bcast(&(this->_detail.residual), 1, root_rank);
    comm->bcast((int*)&(this->_detail.state), 1, root_rank);
  }

  template<typename precision>
  void Status<precision>::log(el::base::type::ostream_t& os) const
  {
    os << "Status("
       << "t=" << this->get_time()
       << ", dt=" << this->get_dt()
       << ", t_end=" << this->get_t_end()
       << ", k_max=" << this->get_max_iterations()
       << ")";
  }
}  // ::pfasst
