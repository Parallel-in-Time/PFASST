#ifndef _PFASST__QUADRATURE__UNIFORM_HPP_
#define _PFASST__QUADRATURE__UNIFORM_HPP_

#include <cassert>
#include <vector>

#include "../interfaces.hpp"
#include "interface.hpp"

using namespace std;


namespace pfasst
{
  namespace quadrature
  {
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
        explicit Uniform(const size_t num_nodes)
          : IQuadrature<precision>(num_nodes)
        {
          if (this->num_nodes < 2) {
            throw invalid_argument("Uniform quadrature requires at least two quadrature nodes.");
          }
          this->compute_nodes();
          this->compute_weights();
          this->compute_delta_nodes();
        }

        Uniform() = default;

        virtual ~Uniform() = default;
        //! @}

        //! @{
        virtual bool left_is_node() const { return LEFT_IS_NODE; }

        virtual bool right_is_node() const { return RIGHT_IS_NODE; }
        //! @}

      protected:
        //! @{
        virtual void compute_nodes()
        {
          this->nodes = vector<precision>(this->num_nodes, precision(0.0));
          for (size_t j = 0; j < this->num_nodes; j++) {
            this->nodes[j] = precision(j) / (this->num_nodes - 1);
          }
        }
        //! @}
    };
  }  // ::pfasst::quadrature
}  // ::pfasst

#endif  // _PFASST__QUADRATURE__UNIFORM_HPP_
