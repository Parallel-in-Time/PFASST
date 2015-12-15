#include <algorithm>
#include <cassert>

#include "pfasst/encap/vector.hpp"

namespace pfasst
{
  namespace encap
  {
    template<typename scalar, typename time>
    VectorEncapsulation<scalar, time>::VectorEncapsulation(const size_t size)
      : vector<scalar>(size)
    {
      zero();
    }

    template<typename scalar, typename time>
    VectorEncapsulation<scalar, time>::VectorEncapsulation(const VectorEncapsulation<scalar, time>& other)
      : vector<scalar>(other)
    {}

    template<typename scalar, typename time>
    VectorEncapsulation<scalar, time>::VectorEncapsulation(const Encapsulation<time>& other)
      : VectorEncapsulation(dynamic_cast<const VectorEncapsulation<scalar, time>>(other))
    {}

    template<typename scalar, typename time>
    VectorEncapsulation<scalar, time>::VectorEncapsulation(VectorEncapsulation<scalar, time>&& other)
      : vector<scalar>(other)
    {}

    template<typename scalar, typename time>
    VectorEncapsulation<scalar, time>::VectorEncapsulation(Encapsulation<time>&& other)
      : VectorEncapsulation(dynamic_cast<VectorEncapsulation<scalar, time>&&>(other))
    {}

    template<typename scalar, typename time>
    VectorEncapsulation<scalar, time>::~VectorEncapsulation()
    {
#ifdef WITH_MPI
      if (this->send_request != MPI_REQUEST_NULL) {
        MPI_Status stat = MPI_Status_factory();
        ML_CLOG(DEBUG, "Encap", "waiting for open send request");
        int err = MPI_Wait(&(this->send_request), &stat);
        check_mpi_error(err);
        ML_CLOG(DEBUG, "Encap", "waited for open send request");
      }
      assert(this->recv_request == MPI_REQUEST_NULL);
      assert(this->send_request == MPI_REQUEST_NULL);
#endif
    }

    template<typename scalar, typename time>
    void VectorEncapsulation<scalar, time>::zero()
    {
      this->assign(this->size(), scalar(0.0));
    }

    template<typename scalar, typename time>
    void VectorEncapsulation<scalar, time>::copy(shared_ptr<const Encapsulation<time>> x)
    {
      shared_ptr<const VectorEncapsulation<scalar, time>> x_cast = \
        dynamic_pointer_cast<const VectorEncapsulation<scalar, time>>(x);
      assert(x_cast);
      this->copy(x_cast);
    }

    template<typename scalar, typename time>
    void VectorEncapsulation<scalar, time>::copy(shared_ptr<const VectorEncapsulation<scalar, time>> x)
    {
      std::copy(x->cbegin(), x->cend(), this->begin());
    }

    template<typename scalar, typename time>
    void VectorEncapsulation<scalar, time>::saxpy(time a, shared_ptr<const Encapsulation<time>> x)
    {
      shared_ptr<const VectorEncapsulation<scalar, time>> x_cast = \
        dynamic_pointer_cast<const VectorEncapsulation<scalar, time>>(x);
      assert(x_cast);

      this->saxpy(a, x_cast);
    }

    template<typename scalar, typename time>
    void VectorEncapsulation<scalar, time>::saxpy(time a, shared_ptr<const VectorEncapsulation<scalar, time>> x)
    {
      assert(this->size() == x->size());
      for (size_t i = 0; i < this->size(); i++)
      { this->data()[i] += a * x->data()[i]; }
    }

    template<typename scalar, typename time>
    void
    VectorEncapsulation<scalar, time>::mat_apply(vector<shared_ptr<Encapsulation<time>>> dst,
                                                 time a, Matrix<time> mat,
                                                 vector<shared_ptr<Encapsulation<time>>> src,
                                                 bool zero)
    {
      size_t ndst = dst.size();
      size_t nsrc = src.size();

      vector<shared_ptr<VectorEncapsulation<scalar, time>>> dst_cast(ndst), src_cast(nsrc);
      for (size_t n = 0; n < ndst; n++) {
        dst_cast[n] = dynamic_pointer_cast<VectorEncapsulation<scalar, time>>(dst[n]);
        assert(dst_cast[n]);
      }
      for (size_t m = 0; m < nsrc; m++) {
        src_cast[m] = dynamic_pointer_cast<VectorEncapsulation<scalar, time>>(src[m]);
        assert(src_cast[m]);
      }

      dst_cast[0]->mat_apply(dst_cast, a, mat, src_cast, zero);
    }

    template<typename scalar, typename time>
    void
    VectorEncapsulation<scalar, time>::mat_apply(vector<shared_ptr<VectorEncapsulation<scalar, time>>> dst,
                                                 time a, Matrix<time> mat,
                                                 vector<shared_ptr<VectorEncapsulation<scalar, time>>> src,
                                                 bool zero)
    {
      size_t ndst = dst.size();
      size_t nsrc = src.size();

      if (zero) { for (auto elem : dst) { elem->zero(); } }

      size_t ndofs = dst[0]->size();
      for (size_t i = 0; i < ndofs; i++) {
        for (size_t n = 0; n < ndst; n++) {
          assert(dst[n]->size() == ndofs);
          for (size_t m = 0; m < nsrc; m++) {
            assert(src[m]->size() == ndofs);
            dst[n]->data()[i] += a * mat(n, m) * src[m]->data()[i];
          }
        }
      }
    }

    template<typename scalar, typename time>
    time VectorEncapsulation<scalar, time>::norm0() const
    {
      return std::abs(*std::max_element(this->cbegin(), this->cend(),
                                        [](scalar a, scalar b) {return std::abs(a) < std::abs(b); } ));
    }


    template<typename scalar, typename time>
    VectorFactory<scalar, time>::VectorFactory(const size_t size)
      : size(size)
    {}

    template<typename scalar, typename time>
    size_t VectorFactory<scalar, time>::dofs() const
    {
      return size;
    }

    template<typename scalar, typename time>
    shared_ptr<Encapsulation<time>> VectorFactory<scalar, time>::create(const EncapType)
    {
      return make_shared<VectorEncapsulation<scalar, time>>(this->dofs());
    }


    template<typename scalar, typename time>
    VectorEncapsulation<scalar,time>& as_vector(shared_ptr<Encapsulation<time>> x)
    {
      typedef VectorEncapsulation<scalar,time> VectorT;
      shared_ptr<VectorT> y = dynamic_pointer_cast<VectorT>(x);
      assert(y);
      return *y.get();
    }

    template<typename scalar, typename time>
    const VectorEncapsulation<scalar,time>& as_vector(shared_ptr<const Encapsulation<time>> x)
    {
      typedef VectorEncapsulation<scalar,time> VectorT;
      shared_ptr<const VectorT> y = dynamic_pointer_cast<const VectorT>(x);
      assert(y);
      return *y.get();
    }

#ifdef WITH_MPI
    template<typename scalar, typename time>
    void VectorEncapsulation<scalar, time>::post(ICommunicator* comm, int tag)
    {
      auto& mpi = as_mpi(comm);
      if (mpi.size() == 1) { return; }
      if (mpi.rank() == 0) { return; }

      if (this->recv_request != MPI_REQUEST_NULL) {
        throw MPIError("a previous receive request is still open");
      }

      int src = (mpi.rank() - 1) % mpi.size();
      ML_CLOG(DEBUG, "Encap", "non-blocking receiving from rank " << src << " with tag=" << tag);
      int err = MPI_Irecv(this->data(), sizeof(scalar) * this->size(), MPI_CHAR,
                          src, tag, mpi.comm, &this->recv_request);
      check_mpi_error(err);
      ML_CLOG(DEBUG, "Encap", "non-blocking received from rank " << src << " with tag=" << tag);
    }

    template<typename scalar, typename time>
    void VectorEncapsulation<scalar, time>::recv(ICommunicator* comm, int tag, bool blocking)
    {
      auto& mpi = as_mpi(comm);
      if (mpi.size() == 1) { return; }
      if (mpi.rank() == 0) { return; }

      MPI_Status stat = MPI_Status_factory();
      int err = MPI_SUCCESS;

      if (blocking) {
        int src = (mpi.rank() - 1) % mpi.size();
        ML_CLOG(DEBUG, "Encap", "blocking receive from rank " << src << " with tag=" << tag);
        err = MPI_Recv(this->data(), sizeof(scalar) * this->size(), MPI_CHAR,
                       src, tag, mpi.comm, &stat);
        check_mpi_error(err);
        ML_CLOG(DEBUG, "Encap", "received blocking from rank " << src << " with tag=" << tag << ": " << stat);
      } else {
        if (this->recv_request != MPI_REQUEST_NULL) {
          ML_CLOG(DEBUG, "Encap", "waiting on last receive request");
          err = MPI_Wait(&(this->recv_request), &stat);
          check_mpi_error(err);
          ML_CLOG(DEBUG, "Encap", "waited on last receive request: " << stat);
        }
      }
    }

    template<typename scalar, typename time>
    void VectorEncapsulation<scalar, time>::send(ICommunicator* comm, int tag, bool blocking)
    {
      auto& mpi = as_mpi(comm);
      if (mpi.size() == 1) { return; }
      if (mpi.rank() == mpi.size() - 1) { return; }

      MPI_Status stat = MPI_Status_factory();
      int err = MPI_SUCCESS;
      int dest = (mpi.rank() + 1) % mpi.size();

      if (blocking) {
        ML_CLOG(DEBUG, "Encap", "blocking send to rank " << dest << " with tag=" << tag);
        err = MPI_Send(this->data(), sizeof(scalar) * this->size(), MPI_CHAR, dest, tag, mpi.comm);
        check_mpi_error(err);
        ML_CLOG(DEBUG, "Encap", "sent blocking to rank " << dest << " with tag=" << tag);
      } else {
        // got never in here
        ML_CLOG(DEBUG, "Encap", "waiting on last send request to finish");
        err = MPI_Wait(&(this->send_request), &stat);
        check_mpi_error(err);
        ML_CLOG(DEBUG, "Encap", "waited on last send request: " << stat);
        ML_CLOG(DEBUG, "Encap", "non-blocking sending to rank " << dest << " with tag=" << tag);
        err = MPI_Isend(this->data(), sizeof(scalar) * this->size(), MPI_CHAR,
                        dest, tag, mpi.comm, &(this->send_request));
        check_mpi_error(err);
        ML_CLOG(DEBUG, "Encap", "sent non-blocking to rank " << dest << " with tag=" << tag);
      }
    }

    template<typename scalar, typename time>
    void VectorEncapsulation<scalar, time>::broadcast(ICommunicator* comm)
    {
      auto& mpi = as_mpi(comm);
      ML_CLOG(DEBUG, "Encap", "broadcasting");
      int err = MPI_Bcast(this->data(), sizeof(scalar) * this->size(), MPI_CHAR,
                          comm->size()-1, mpi.comm);
      check_mpi_error(err);
      ML_CLOG(DEBUG, "Encap", "broadcasted");
    }
#endif

  }  // ::pfasst::encap
}  // ::pfasst
