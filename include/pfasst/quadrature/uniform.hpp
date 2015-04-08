/**
 * @file pfasst/quadrature/uniform.hpp
 * @since v0.3.0
 */
#ifndef _PFASST__QUADRATURE__UNIFORM_HPP_
#define _PFASST__QUADRATURE__UNIFORM_HPP_

#include "pfasst/quadrature/interface.hpp"


namespace pfasst
{
  namespace quadrature
  {
    /**
     * Quadrature handler for uniform distributed nodes.
     *
     * @tparam scalar precision of quadrature (i.e. `double`)
     *
     * @since v0.3.0
     */
    template<typename precision = pfasst::time_precision>
    class Uniform
      : public IQuadrature<precision>
    {
      protected:
        //! @{
        static const bool LEFT_IS_NODE = true;
        static const bool RIGHT_IS_NODE = true;
        //! @}

      public:
        //! @{
        /**
         * @throws invalid_argument if less than two nodes are requested
         */
        explicit Uniform(const size_t num_nodes);
        Uniform() = default;
        virtual ~Uniform() = default;
        //! @}

        //! @{
        virtual bool left_is_node() const override;
        virtual bool right_is_node() const override;
        //! @}

      protected:
        //! @{
        virtual void compute_nodes() override;
        //! @}
    };
  }  // ::pfasst::quadrature
}  // ::pfasst

#include "pfasst/quadrature/uniform_impl.hpp"

#endif  // _PFASST__QUADRATURE__UNIFORM_HPP_
