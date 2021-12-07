
#include <iostream>
#include <array>

#include "Polynomial.h"
#include "PolyInterpolation.h"

#include <type_traits>




int main(){


  constexpr std::array<double, 3> roots {0, 1, 2};

  constexpr Polynomial<3, double> poly(roots, 1);

  constexpr auto poly2 = poly.derivative();

  constexpr auto poly3 = poly2 * 2;

  constexpr auto poly4 = poly3.derivative<2>();
  constexpr auto poly5 = poly4.integrate();
  constexpr double integral = poly4.integrate(-0.5, 0.5);


  constexpr std::array<double, 3> x{-1, 0, 1};
  constexpr std::array<double, 3> f{1, 0, 1};
  constexpr Polynomial<2, double> poly_interpolation = lagrange_interpolate_to_polynomial(x, f);


  constexpr std::array<double, 3> interpolation_coefficients = lagrange_interpolation_coefficients(x, 4.0);


  constexpr std::array<double, 3> interp_coeff = calculate_interpolation_coefficients_from_cell_to_face_for_uniformly_spaced_points<3, double>(0);

  constexpr std::array<double, 1> reconstruction_coeff = calculate_reconstruction_coefficients_from_cell_to_face_for_uniformly_spaced_points<1, double>(0);


  constexpr auto c = WENO_Uniform_Reconstruction_Coefficients<3, double>::c;

  constexpr auto d = WENO_Uniform_Reconstruction_Coefficients<3, double>::d;

  constexpr auto s = WENO_Uniform_Reconstruction_Coefficients<3, double>::s;


  return 0;
}
