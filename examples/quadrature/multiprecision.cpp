/*
 * Quick example of how to use Boost multiprecision reals to compute quadrature weights.
 */


#include <iostream>
#include <iomanip>

#include<pfasst/quadrature.hpp>

#include <boost/multiprecision/cpp_dec_float.hpp>
using boost::multiprecision::cpp_dec_float_50;

int main(int argc, char *argv[])
{
  pfasst::quadrature::GaussLobatto<cpp_dec_float_50> quad(3);

  auto Q = quad.get_q_mat();
  std::cout << std::setprecision(std::numeric_limits<cpp_dec_float_50>::digits10)
            << Q << std::endl;
}
