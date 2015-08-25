#ifndef _PFASST__CONTROLLER__STATUS_HPP_
#define _PFASST__CONTROLLER__STATUS_HPP_

#include <limits>
#include <memory>
#include <type_traits>
using namespace std;

#ifdef WITH_MPI
  #include <mpi.h>
#endif

#include "pfasst/logging.hpp"


namespace pfasst
{
#ifdef WITH_MPI
  static MPI_Datatype status_data_type;
#endif

  enum State : int {
    // overall state
    CONVERGED        =  0,
    FAILED           =  1,

    // iterating states
    PREDICTING       = 10,
    ITERATING        = 11,

    // coarse level
    PRE_ITER_COARSE  = 20,
    ITER_COARSE      = 21,
    POST_ITER_COARSE = 22,

    // fine level
    PRE_ITER_FINE    = 30,
    ITER_FINE        = 31,
    POST_ITER_FINE   = 32,

    // last
    UNKNOWN          = numeric_limits<int>::max(),
  };


  template<
    typename precision
  >
  struct StatusDetail {
    State     state          = pfasst::State::UNKNOWN;

    size_t    step           = 0;
    size_t    num_steps      = 0;
    size_t    iteration      = 0;
    size_t    max_iterations = 0;

    precision time           = 0.0;
    precision dt             = 0.0;
    precision t_end          = 0.0;
    precision abs_res_norm   = 0.0;
    precision rel_res_norm   = 0.0;
  };


  template<
    typename precision
  >
  class Status
    : public enable_shared_from_this<Status<precision>>
      , public el::Loggable
  {
    static_assert(is_arithmetic<precision>::value,
                  "precision type must be arithmetic");

    public:
      typedef precision precision_type;

#ifdef WITH_MPI
      static inline void create_mpi_datatype();
      static inline void free_mpi_datatype();
#endif

      static_assert(is_standard_layout<StatusDetail<precision>>::value,
                    "Status::Detail needs to be have standard layout for MPI derived datatype");
      StatusDetail<precision> _detail;

    public:
      Status();
      Status(const Status<precision>& other) = default;
      Status(Status<precision>&& other) = default;
      ~Status() = default;
      Status<precision>& operator=(const Status<precision>& other) = default;
      Status<precision>& operator=(Status<precision>&& other) = default;

      virtual size_t& step();
      virtual size_t  get_step() const;

      virtual size_t& num_steps();
      virtual size_t  get_num_steps() const;

      virtual size_t& iteration();
      virtual size_t  get_iteration() const;

      virtual size_t& max_iterations();
      virtual size_t  get_max_iterations() const;

      virtual precision& time();
      virtual precision  get_time() const;

      virtual precision& dt();
      virtual precision  get_dt() const;

      virtual precision& t_end();
      virtual precision  get_t_end() const;

      virtual State& state();
      virtual State  get_state() const;

      virtual precision& abs_res_norm();
      virtual precision  get_abs_res_norm() const;

      virtual precision& rel_res_norm();
      virtual precision  get_rel_res_norm() const;

      template<class CommT>
      void send(shared_ptr<CommT> comm, const int dest_rank, const int tag, const bool blocking);

      template<class CommT>
      void recv(shared_ptr<CommT> comm, const int src_rank, const int tag, const bool blocking);

      template<class CommT>
      void bcast(shared_ptr<CommT> comm, const int root_rank);

      virtual vector<string> summary() const;
      virtual void log(el::base::type::ostream_t& os) const;
  };
}  // ::pfasst


#include "pfasst/controller/status_impl.hpp"

#endif  // _PFASST__CONTROLLER__STATUS_HPP_
