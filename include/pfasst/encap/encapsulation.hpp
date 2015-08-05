#ifndef _PFASST__ENCAP__INTERFACE_HPP_
#define _PFASST__ENCAP__INTERFACE_HPP_

#include <memory>
#include <type_traits>
using namespace std;

#include <leathers/push>
#include <leathers/all>
#include <Eigen/Dense>
#include <leathers/pop>
template<typename precision>
using Matrix = Eigen::Matrix<precision, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"
#include "pfasst/encap/traits.hpp"
#include "pfasst/comm/communicator.hpp"


namespace pfasst
{
  namespace encap
  {
    template<
      class EncapsulationTrait,
      typename Enabled = void
    >
    class Encapsulation;

    template<
      class EncapsulationTrait,
      typename Enabled = void
    >
    class EncapsulationFactory;


    //! a*x + y
    template<
      class EncapsulationTrait
    >
    shared_ptr<Encapsulation<EncapsulationTrait>>
    axpy(const typename EncapsulationTrait::time_type& a,
         const shared_ptr<Encapsulation<EncapsulationTrait>> x,
         const shared_ptr<Encapsulation<EncapsulationTrait>> y);

    //! x = x + a*matrix*y
    template<
      class EncapsulationTrait
    >
    void
    mat_apply(vector<shared_ptr<Encapsulation<EncapsulationTrait>>>& x,
              const typename EncapsulationTrait::time_type& a,
              const Matrix<typename EncapsulationTrait::time_type>& matrix,
              const vector<shared_ptr<Encapsulation<EncapsulationTrait>>>& y,
              const bool zero_vec_x = true);

    //! a*matrix*x
    template<
      class EncapsulationTrait
    >
    vector<shared_ptr<Encapsulation<EncapsulationTrait>>>
    mat_mul_vec(const typename EncapsulationTrait::time_type& a,
                const Matrix<typename EncapsulationTrait::time_type>& matrix,
                const vector<shared_ptr<Encapsulation<EncapsulationTrait>>>& x);

    template<
      class EncapsulationTrait
    >
    typename EncapsulationTrait::spacial_type
    norm0(const shared_ptr<Encapsulation<EncapsulationTrait>> x);


    template<
      class EncapsulationTrait,
      typename Enabled
    >
    class Encapsulation
      :   public enable_shared_from_this<Encapsulation<EncapsulationTrait>>
        , el::Loggable
    {
      public:
        typedef          EncapsulationTrait            traits;
        typedef typename traits::time_type             time_type;
        typedef typename traits::spacial_type          spacial_type;
        typedef typename traits::data_type             data_type;
        typedef          EncapsulationFactory<traits>  factory_type;

      static_assert(is_arithmetic<time_type>::value,
                    "time precision must be an arithmetic type");
      static_assert(is_arithmetic<spacial_type>::value,
                    "spacial precision must be an arithmetic type");
      static_assert(is_constructible<data_type>::value,
                    "Data Type must be constructible");
      static_assert(is_default_constructible<data_type>::value,
                    "Data Type must be default constructible");
      static_assert(is_destructible<data_type>::value,
                    "Data Type must be destructible");
      static_assert(is_assignable<data_type, data_type>::value,
                    "Data Type must be assignable");

      STATIC_WARNING(is_move_constructible<data_type>::value,
                     "Data Type should be move constructible");
      STATIC_WARNING(is_copy_constructible<data_type>::value,
                     "Data Type should be copy constructible");
      STATIC_WARNING(is_move_assignable<data_type>::value,
                     "Data Type should be move assignable");
      STATIC_WARNING(is_copy_assignable<data_type>::value,
                     "Data Type should be copy assignable");

      protected:
        data_type _data;

      public:
        Encapsulation() = default;
        Encapsulation(const typename EncapsulationTrait::data_type& data);
        Encapsulation(const Encapsulation<EncapsulationTrait>& other) = default;
        Encapsulation(Encapsulation<EncapsulationTrait>&& other) = default;
        ~Encapsulation() = default;
        Encapsulation<EncapsulationTrait>& operator=(const typename EncapsulationTrait::data_type& data);
        Encapsulation<EncapsulationTrait>& operator=(const Encapsulation<EncapsulationTrait>& other) = default;
        Encapsulation<EncapsulationTrait>& operator=(Encapsulation<EncapsulationTrait>&& other) = default;

        virtual       typename EncapsulationTrait::data_type& data();
        virtual const typename EncapsulationTrait::data_type& get_data() const;

        virtual void zero();
        //! this += a*y
        virtual void scaled_add(const typename EncapsulationTrait::time_type& a,
                                const shared_ptr<Encapsulation<EncapsulationTrait>> y);

        virtual typename EncapsulationTrait::spacial_type norm0() const;

        virtual void send(shared_ptr<comm::Communicator> comm, const int dest_rank, const int tag,
                          const bool blocking);
        virtual void recv(shared_ptr<comm::Communicator> comm, const int src_rank, const int tag,
                          const bool blocking);
        virtual void bcast(shared_ptr<comm::Communicator> comm, const int root_rank);

        virtual void log(el::base::type::ostream_t& os) const override;
    };


    template<
      class EncapsulationTrait,
      class Enabled
    >
    class EncapsulationFactory
    {
      public:
        typedef          Encapsulation<EncapsulationTrait> encap_type;
        typedef typename EncapsulationTrait::data_type     data_type;

        EncapsulationFactory() = default;
        EncapsulationFactory(const EncapsulationFactory<EncapsulationTrait>& other) = default;
        EncapsulationFactory(EncapsulationFactory<EncapsulationTrait>&& other) = default;
        virtual ~EncapsulationFactory() = default;
        EncapsulationFactory<EncapsulationTrait>&
        operator=(const EncapsulationFactory<EncapsulationTrait>& other) = default;
        EncapsulationFactory<EncapsulationTrait>&
        operator=(EncapsulationFactory<EncapsulationTrait>&& other) = default;

        shared_ptr<Encapsulation<EncapsulationTrait>> create();
    };
  }  // ::encap
}  // ::pfasst

#include "pfasst/encap/encapsulation_impl.hpp"

#endif  // _PFASST__ENCAP__INTERFACE_HPP_
