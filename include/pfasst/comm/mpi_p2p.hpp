#ifndef _PFASST__COMM__MPI_P2P_HPP_
#define _PFASST__COMM__MPI_P2P_HPP_

#include <list>
#include <map>
#include <memory>
#include <type_traits>
#include <utility>
using namespace std;

#include <mpi.h>

#include "pfasst/comm/communicator.hpp"


namespace pfasst
{
  namespace comm
  {
    static string error_from_code(const int err_code);
    static MPI_Status MPI_Status_factory();
    static void check_mpi_error(const int err_code);


    class MpiP2P
      : public Communicator
    {
      protected:
        int      _size = -1;
        int      _rank = -1;
        string   _name = "";

        MPI_Comm _comm;

        list<MPI_Status>                 _stati;
        map<pair<int, int>, MPI_Request> _requests;

      public:
        explicit MpiP2P(MPI_Comm comm = MPI_COMM_WORLD);
        MpiP2P(const MpiP2P& other) = default;
        MpiP2P(MpiP2P&& other) = default;
        virtual ~MpiP2P() = default;  // TODO: might need to implement destructor to clean up pending MPI  stati
        MpiP2P& operator=(const MpiP2P& other) = default;
        MpiP2P& operator=(MpiP2P&& other) = default;

        virtual size_t get_size() const override;
        virtual size_t get_rank() const override;

        virtual string get_name() const;

        virtual bool is_first() const override;
        virtual bool is_last() const override;

        template<class DataT>
        void send(const DataT* const data, const int count, const int dest_rank, const int tag);

        template<class DataT>
        void isend(const DataT* const data, const int count, const int dest_rank, const int tag);

        template<class DataT>
        void recv(DataT* data, const int count, const int dest_rank, const int tag);

        template<class DataT>
        void irecv(DataT* data, const int count, const int src_rank, const int tag);

        template<class DataT>
        void bcast(DataT* data, const int count, const int root_rank);
    };
  }  // ::pfasst::comm
}  // ::pfasst

#include "pfasst/comm/mpi_p2p_impl.hpp"

#endif  // _PFASST__COMM__MPI_P2P_HPP_
