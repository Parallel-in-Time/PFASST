#ifndef _PFASST__ENCAP__VECTOR_HPP_
#define _PFASST__ENCAP__VECTOR_HPP_

#include <memory>
#include <type_traits>
#include <vector>
using namespace std;

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"
#include "pfasst/encap/interface.hpp"
#include "pfasst/comm/interface.hpp"


namespace pfasst
{
  namespace encap
  {
    template<
      class EncapsulationTrait
    >
    class Encapsulation<EncapsulationTrait,
                        typename enable_if<
                                   is_same<
                                     vector<typename EncapsulationTrait::spacial_type>,
                                     typename EncapsulationTrait::data_type
                                   >::value>::type>
      : public enable_shared_from_this<Encapsulation<EncapsulationTrait>>
    {
      public:
        typedef typename EncapsulationTrait::time_type            time_type;
        typedef typename EncapsulationTrait::spacial_type         spacial_type;
        typedef typename EncapsulationTrait::data_type            data_type;
        typedef          EncapsulationFactory<EncapsulationTrait> factory_type;

      protected:
        data_type _data;

      public:
        Encapsulation() = default;
        Encapsulation(const size_t size);
        Encapsulation(const typename EncapsulationTrait::data_type& data);
        Encapsulation(const Encapsulation<EncapsulationTrait>& other) = default;
        Encapsulation(Encapsulation<EncapsulationTrait>&& other) = default;
        virtual ~Encapsulation() = default;
        Encapsulation<EncapsulationTrait>& operator=(const typename EncapsulationTrait::data_type& data);
        Encapsulation<EncapsulationTrait>& operator=(const Encapsulation<EncapsulationTrait>& other) = default;
        Encapsulation<EncapsulationTrait>& operator=(Encapsulation<EncapsulationTrait>&& other) = default;

        virtual       typename EncapsulationTrait::data_type& data();
        virtual const typename EncapsulationTrait::data_type& get_data() const;

        virtual void zero();
        //! this += a*y
        virtual void scale_add(const typename EncapsulationTrait::time_type& a,
                               const shared_ptr<Encapsulation<EncapsulationTrait>> y);

        virtual typename EncapsulationTrait::spacial_type norm0() const;

        virtual void send(shared_ptr<comm::Communicator> comm, const int dest_rank, const int tag,
                          const bool blocking);
        virtual void recv(shared_ptr<comm::Communicator> comm, const int src_rank, const int tag,
                          const bool blocking);
        virtual void bcast(shared_ptr<comm::Communicator> comm, const int root_rank);
    };

    template<
      typename time_precision,
      typename spacial_precision
    >
    using VectorEncapsulation = Encapsulation<vector_encap_traits<time_precision, spacial_precision>>;


    template<
      class EncapsulationTrait
    >
    class EncapsulationFactory<EncapsulationTrait,
                               typename enable_if<
                                          is_same<
                                            vector<typename EncapsulationTrait::spacial_type>,
                                            typename EncapsulationTrait::data_type
                                          >::value>::type>
      : public enable_shared_from_this<EncapsulationFactory<EncapsulationTrait>>
    {
      protected:
        size_t _size;

      public:
        explicit EncapsulationFactory(const size_t size=0);
        EncapsulationFactory(const EncapsulationFactory<EncapsulationTrait>& other) = default;
        EncapsulationFactory(EncapsulationFactory<EncapsulationTrait>&& other) = default;
        virtual ~EncapsulationFactory() = default;
        EncapsulationFactory<EncapsulationTrait>& operator=(const EncapsulationFactory<EncapsulationTrait>& other) = default;
        EncapsulationFactory<EncapsulationTrait>& operator=(EncapsulationFactory<EncapsulationTrait>&& other) = default;

        virtual shared_ptr<Encapsulation<EncapsulationTrait>> create();

        virtual void set_size(const size_t& size);
        virtual size_t size() const;
    };
  }  // ::pfasst::encap
}  // ::pfasst

#include "pfasst/encap/vector_impl.hpp"

#endif  // _PFASST__ENCAP__VECTOR_HPP_
