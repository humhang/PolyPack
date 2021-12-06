
#include <iostream>
#include <array>

#include "Polynomial.h"
#include "PolyInterpolation.h"

#include <type_traits>




int main(){


  constexpr std::array<int, 3> roots {0, 1, 2};

  constexpr Polynomial<3, int> poly(roots, 1);

  constexpr auto poly2 = poly.derivative();

  constexpr auto poly3 = poly2 * 2;

  std::cout << poly3.evaluate(3);

  constexpr std::array<double, 3> x{-1, 0, 1};
  constexpr std::array<double, 3> f{1, 0, 1};
  constexpr Polynomial<2, double> poly_interpolation = lagrange_interpolate_to_polynomial(x, f);
  std::cout << poly_interpolation.evaluate(4);


  constexpr std::array<double, 3> interpolation_coefficients = lagrange_interpolation_coefficients(x, 4.0);
  std::cout << interpolation_coefficients[0] + interpolation_coefficients[2];


  constexpr std::array<double, 3> interp_coeff = calculate_interpolation_coefficients_from_cell_to_face_for_uniformly_spaced_points<3, double>(0);

  constexpr std::array<double, 3> reconstruction_coeff = calculate_reconstruction_coefficients_from_cell_to_face_for_uniformly_spaced_points<3, double>(0);

  constexpr auto c = WENO_Uniform_Reconstruction_Coefficients<2, double>::c;

  constexpr auto d = WENO_Uniform_Reconstruction_Coefficients<2, double>::d;

  return 0;
}
