#ifndef _PFASST_VECTOR_HPP_
#define _PFASST_VECTOR_HPP_

#include <memory>
#include <vector>
using namespace std;

#ifdef WITH_MPI
#include "pfasst/mpi_communicator.hpp"
using namespace pfasst::mpi;
#endif

#include "pfasst/encap/encapsulation.hpp"


namespace pfasst
{
  namespace encap
  {
    /**
     * @tparam scalar
     *     precision and numerical type of the data values
     * @tparam time
     *     precision of the time points; defaults to pfasst::time_precision
     */
    template<typename scalar, typename time = time_precision>
    class VectorEncapsulation
      : public vector<scalar>,
        public Encapsulation<time>
    {
      public:

        //! @{
        VectorEncapsulation(const size_t size);

        /**
         * Copy constuctor.
         *
         * @note delegated to sdt::vector<scalar>
         */
        VectorEncapsulation(const VectorEncapsulation<scalar, time>& other);

        /**
         * @throws std::bad_cast
         *     if `other` can not be transformed into pfasst::encap::VectorEncapsulation via
         *     `dynamic_cast`
         */
        VectorEncapsulation(const Encapsulation<time>& other);

        /**
         * Move constructor.
         *
         * @note delegated to std::vector<scalar>
         */
        VectorEncapsulation(VectorEncapsulation<scalar, time>&& other);

        /**
         * @throws std::bad_cast
         *     if `other` can not be transformed into pfasst::encap::VectorEncapsulation via
         *     `dynamic_cast`
         */
        VectorEncapsulation(Encapsulation<time>&& other);

        virtual ~VectorEncapsulation();
        //! @}

        //! @{
        virtual void zero() override;
        virtual void copy(shared_ptr<const Encapsulation<time>> x) override;
        virtual void copy(shared_ptr<const VectorEncapsulation<scalar, time>> x);
        //! @}

        //! @{
        virtual void saxpy(time a, shared_ptr<const Encapsulation<time>> x) override;
        virtual void saxpy(time a, shared_ptr<const VectorEncapsulation<scalar, time>> x);

        /**
         * @note In case any of the elements of `dst` or `src` can not be transformed via
         *     `dynamic_cast` into pfasst::encap::VectorEncapsulation std::abort is called.
         */
        virtual void mat_apply(vector<shared_ptr<Encapsulation<time>>> dst,
                               time a, Matrix<time> mat,
                               vector<shared_ptr<Encapsulation<time>>> src,
                               bool zero = true) override;
        virtual void mat_apply(vector<shared_ptr<VectorEncapsulation<scalar, time>>> dst,
                               time a, Matrix<time> mat,
                               vector<shared_ptr<VectorEncapsulation<scalar, time>>> src,
                               bool zero = true);

        /**
         * Maximum norm of contained elements.
         *
         * This uses std::max with custom comparison function.
         */
        virtual time norm0() const override;
        //! @}

#ifdef WITH_MPI
        //! @{
        MPI_Request recv_request = MPI_REQUEST_NULL;
        MPI_Request send_request = MPI_REQUEST_NULL;
        //! @}

        //! @{
        inline MPICommunicator& as_mpi(ICommunicator* comm)
        {
          auto mpi = dynamic_cast<MPICommunicator*>(comm);
          assert(mpi);
          return *mpi;
        }
        //! @}

        //! @{
        virtual void post(ICommunicator* comm, int tag) override;
        virtual void recv(ICommunicator* comm, int tag, bool blocking) override;
        virtual void send(ICommunicator* comm, int tag, bool blocking) override;
        virtual void broadcast(ICommunicator* comm) override;
        //! @}
#endif

    };

    /**
     * @tparam scalar
     *     precision and numerical type of the data values
     * @tparam time
     *     precision of the time points; defaults to pfasst::time_precision
     */
    template<typename scalar, typename time = time_precision>
    class VectorFactory
      : public EncapFactory<time>
    {
      protected:
        size_t size;

      public:
        VectorFactory(const size_t size);
        virtual shared_ptr<Encapsulation<time>> create(const EncapType) override;
        size_t dofs() const;
    };

    template<typename scalar, typename time = time_precision>
    VectorEncapsulation<scalar,time>& as_vector(shared_ptr<Encapsulation<time>> x);

    template<typename scalar, typename time = time_precision>
    const VectorEncapsulation<scalar,time>& as_vector(shared_ptr<const Encapsulation<time>> x);
  }  // ::pfasst::encap
}  // ::pfasst

#include "pfasst/encap/vector_impl.hpp"

#endif
