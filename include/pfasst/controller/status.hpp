#ifndef _PFASST__CONTROLLER__STATUS_HPP_
#define _PFASST__CONTROLLER__STATUS_HPP_

#include <limits>
#include <memory>
#include <type_traits>
using namespace std;

#include "pfasst/comm/interface.hpp"


namespace pfasst
{
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
  class Status
    : public enable_shared_from_this<Status<precision>>
  {
    static_assert(is_arithmetic<precision>::value,
                  "precision type must be arithmetic");

    public:
      typedef precision precision_type;

    protected:
      // TODO: create derived data type and API for creating it
      struct Detail {
        Detail();
        Detail(const Detail& other) = default;
        Detail(Detail&& other) = default;
        virtual ~Detail() = default;
        Detail& operator=(const Detail& other) = default;
        Detail& operator=(Detail&& other) = default;

        size_t    step;
        size_t    iteration;

        precision time;
        precision dt;

        precision t_end;
        size_t    max_iterations;

        State     state;
        precision residual;
      } _detail;

    public:
      Status();
      Status(const Status<precision>& other) = default;
      Status(Status<precision>&& other) = default;
      ~Status() = default;
      Status<precision>& operator=(const Status<precision>& other) = default;
      Status<precision>& operator=(Status<precision>&& other) = default;

      virtual size_t& step();
      virtual size_t  get_step() const;

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

      virtual precision& residual();
      virtual precision  get_residual() const;

      virtual void send(shared_ptr<comm::Communicator> comm,
                        const int dest_rank, const int tag, const bool blocking);

      virtual void recv(shared_ptr<comm::Communicator> comm,
                        const int src_rank, const int tag, const bool blocking);

      virtual void bcast(shared_ptr<comm::Communicator> comm, const int root_rank);
  };
}  // ::pfasst


template<typename precision>
string to_string(const pfasst::Status<precision>& status);


#include "pfasst/controller/status_impl.hpp"

#endif  // _PFASST__CONTROLLER__STATUS_HPP_
