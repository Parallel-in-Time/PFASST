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
  /**
   * Dealing with user data.
   */
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


    /**
     * Computes \\( a * x + y \\) .
     *
     * @tparam    EncapsulationTrait type traits for Encapsulations @p x and @p y
     * @param[in] a                  scalar factor to scale @p x with
     * @param[in] x                  first Encapsulation to be scaled by @p a
     * @param[in] y                  second Encapsulation to be added to \\( a*x \\)
     * @returns result of \\( a * x + y \\)
     *
     * @see pfasst::encap_traits
     *  for instances of @p EncapsulationTrait
     */
    template<
      class EncapsulationTrait
    >
    shared_ptr<Encapsulation<EncapsulationTrait>>
    axpy(const typename EncapsulationTrait::time_type& a,
         const shared_ptr<Encapsulation<EncapsulationTrait>> x,
         const shared_ptr<Encapsulation<EncapsulationTrait>> y);

    /**
     * Computes scaled matrix-vector product \\( x += aMy \\) .
     *
     * The matrix-vector product of @p mat and @p y is scaled by @p a and added onto @p x.
     *
     * @tparam        EncapsulationTrait type traits for Encapsulations @p x and @p y
     * @param[in,out] x                  target Encapsulation
     * @param[in]     a                  scalar factor to scale the matrix-vector product of @p mat and @p y with
     * @param[in]     mat                matrix
     * @param[in]     y                  second Encapsulation used as vector in \\( aMy \\)
     * @param[in]     zero_vec_x         if `true` @p x is zeroed before adding the result of the scaled matrix-vector
     *                                   product
     *
     * @see pfasst::encap_traits
     *  for instances of @p EncapsulationTrait
     * @see pfasst::encap::mat_mul_vec()
     *  for a version not modifying any input data
     */
    template<
      class EncapsulationTrait
    >
    void
    mat_apply(vector<shared_ptr<Encapsulation<EncapsulationTrait>>>& x,
              const typename EncapsulationTrait::time_type& a,
              const Matrix<typename EncapsulationTrait::time_type>& matrix,
              const vector<shared_ptr<Encapsulation<EncapsulationTrait>>>& y,
              const bool zero_vec_x = true);

    /**
     * Computes scaled matrix-vector product \\( aMx \\).
     *
     * @tparam    EncapsulationTrait type traits for Encapsulations @p x
     * @param[in] x                  data used as vector
     * @param[in] a                  scalar factor to scale the matrix-vector product of @p mat and @p x with
     * @param[in] mat                matrix
     * @returns result of \\( aMx \\)
     *
     * @see pfasst::encap_traits
     *  for instances of @p EncapsulationTrait
     * @see pfasst::encap::mat_apply()
     *  for an in-place version
     */
    template<
      class EncapsulationTrait
    >
    vector<shared_ptr<Encapsulation<EncapsulationTrait>>>
    mat_mul_vec(const typename EncapsulationTrait::time_type& a,
                const Matrix<typename EncapsulationTrait::time_type>& matrix,
                const vector<shared_ptr<Encapsulation<EncapsulationTrait>>>& x);

    /**
     * Computes maximums norm of @p x.
     *
     * Computes the maximums or infinity norm of @p x, e.g. \\( \\|x\\|_{\\inf} = max(|x_0|, \\dots, |x_n|) \\).
     *
     * @tparam    EncapsulationTrait type traits for Encapsulations @p x
     * @param[in] x                  data vector to compute infinity/maximums norm of
     * @returns maximums norm
     *
     * @see Encapsulation::norm0()
     */
    template<
      class EncapsulationTrait
    >
    typename EncapsulationTrait::spacial_type
    norm0(const shared_ptr<Encapsulation<EncapsulationTrait>> x);


    /**
     * Encapsulations are the way _PFASST_ can handle arbitrary user data.
     *
     * @tparam EncapsulationTrait type trait describing encapsulated data
     * @tparam Enabled            utility type for template specializations
     *
     * @see pfasst::encap_traits
     *  for instances of @p EncapsulationTrait
     */
    template<
      class EncapsulationTrait,
      typename Enabled
    >
    class Encapsulation
      :   public enable_shared_from_this<Encapsulation<EncapsulationTrait, Enabled>>
        , el::Loggable
    {
      public:
        typedef          EncapsulationTrait            traits;
        typedef typename traits::time_type             time_type;
        typedef typename traits::spacial_type          spacial_type;
        typedef typename traits::data_type             data_type;
        typedef          EncapsulationFactory<traits>  factory_type;

      //! @internal
      //! time_type must be an arithmetic type
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

      //! @cond static_warnings
      STATIC_WARNING(is_move_constructible<data_type>::value,
                     "Data Type should be move constructible");
      STATIC_WARNING(is_copy_constructible<data_type>::value,
                     "Data Type should be copy constructible");
      STATIC_WARNING(is_move_assignable<data_type>::value,
                     "Data Type should be move assignable");
      STATIC_WARNING(is_copy_assignable<data_type>::value,
                     "Data Type should be copy assignable");
      //! @endcond
      //! @endinternal

      protected:
        //! @{
        //! actual storage of encapsulated data
        data_type _data;
        //! @}

      public:
        /*
         * Note on non-default ctors:
         * Explicitly requesting default ctors and copy- and move-assignment operators works fine
         * with clang but triggers compile errors 'cannot be defaulted' with GCC.
         * Thus I had to manually specify these.
         * Move ctor must be noexcept to comply with ISO C++11.
         */
        //! @{
        Encapsulation();
        /**
         * Encapsulating existing data by copying.
         *
         * A copy of @p data will be encapsulated in the newly constructed Encapsulation.
         *
         * @param[in] data
         */
        Encapsulation(const typename EncapsulationTrait::data_type& data);
        Encapsulation(const Encapsulation<EncapsulationTrait, Enabled>& other);
        Encapsulation(Encapsulation<EncapsulationTrait, Enabled>&& other) noexcept;
        virtual ~Encapsulation() = default;
        Encapsulation<EncapsulationTrait, Enabled>& operator=(const typename EncapsulationTrait::data_type& data);
        Encapsulation<EncapsulationTrait, Enabled>& operator=(const Encapsulation<EncapsulationTrait, Enabled>& other);
        Encapsulation<EncapsulationTrait, Enabled>& operator=(Encapsulation<EncapsulationTrait, Enabled>&& other);
        //! @}

        //! @name Accessor
        //! @{
        /**
         * Accessor for encapsulated data for modification.
         *
         * @returns encapsulated data for modification, e.g. as an lvalue
         */
        virtual       typename EncapsulationTrait::data_type& data();
        /**
         * Read-only accessor for encapsulated data.
         *
         * @returns encapsulated data
         */
        virtual const typename EncapsulationTrait::data_type& get_data() const;

        /**
         * Computes maximums norm of encapsulated data.
         *
         * @returns maximums norm of underlying data
         *
         * @note The implementation is strongly dependent on the encapsulated data type.
         */
        virtual typename EncapsulationTrait::spacial_type norm0() const;

        /**
         * Streams string representation of Encapsulation.
         *
         * @param[in,out] os stream used for output
         *
         * @see [Documentation of easylogging++](https://github.com/easylogging/easyloggingpp#logging-your-own-class)
         *  for details on where this comes from
         */
        virtual void log(el::base::type::ostream_t& os) const override;
        //! @}

        //! @name Modification
        //! @{
        /**
         * Zeros the underlying data without changing its reserved and occupied space.
         */
        virtual void zero();

        /**
         * Adds @p y scaled by @p a onto this Encapsulation.
         *
         * @param[in] a scaling factor
         * @param[in] y other data to be scaled and added onto this one
         */
        virtual void scaled_add(const typename EncapsulationTrait::time_type& a,
                                const shared_ptr<Encapsulation<EncapsulationTrait>> y);
        //! @}

        //! @name Communication
        //! @{
        /**
         * Sending encapsulated data over communicator.
         *
         * @param[in] comm      Communicator used for sending
         * @param[in] dest_rank target processor of the data
         * @param[in] tag       accociation of the data
         * @param[in] blocking  `true` for blocking sending, `false` for non-blocking communication
         */
        virtual void send(shared_ptr<comm::Communicator> comm, const int dest_rank, const int tag, const bool blocking);

        /**
         * Receiving encapsulated data over communicator.
         *
         * @param[in] comm      Communicator used for sending
         * @param[in] src_rank  source processor of the data
         * @param[in] tag       accociation of the data
         * @param[in] blocking  `true` for blocking receiving, `false` for non-blocking communication
         */
        virtual void recv(shared_ptr<comm::Communicator> comm, const int src_rank, const int tag, const bool blocking);

        /**
         * Sending encapsulated data over communicator.
         *
         * @param[in] comm      Communicator used for broadcasting
         * @param[in] root_rank source processor of the data
         */
        virtual void bcast(shared_ptr<comm::Communicator> comm, const int root_rank);
        //! @}
    };


    /**
     * Utility to create specific Encapsulations with predefined defaults.
     *
     * @tparam EncapsulationTrait type trait describing encapsulated data
     * @tparam Enabled            utility type for template specializations
     *
     * @note Specializations for certain encapsulated data types may define additional member functions for setting
     *  and accessing default values for instantiated Encapsulations.
     */
    template<
      class EncapsulationTrait,
      class Enabled
    >
    class EncapsulationFactory
    {
      public:
        typedef          Encapsulation<EncapsulationTrait> encap_type;
        typedef typename EncapsulationTrait::data_type     data_type;

        //! @{
        EncapsulationFactory();
        EncapsulationFactory(const EncapsulationFactory<EncapsulationTrait>& other);
        EncapsulationFactory(EncapsulationFactory<EncapsulationTrait>&& other);
        virtual ~EncapsulationFactory() = default;
        EncapsulationFactory<EncapsulationTrait>& operator=(const EncapsulationFactory<EncapsulationTrait>& other);
        EncapsulationFactory<EncapsulationTrait>& operator=(EncapsulationFactory<EncapsulationTrait>&& other);
        //! @}

        //! @{
        /**
         * Instantiating a new Encapsulation.
         *
         * @returns a newly instantiated Encapsulation with implementation dependent defaults.
         */
        shared_ptr<Encapsulation<EncapsulationTrait>> create() const;
        //! @}
    };
  }  // ::encap
}  // ::pfasst


#include "pfasst/encap/encapsulation_impl.hpp"

#endif  // _PFASST__ENCAP__INTERFACE_HPP_
