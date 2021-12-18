
#include <iostream>
#include <array>

#include "Polynomial.h"
#include "LagrangeInterpolation.h"
#include "UniformGridInterpolationAndReconstruction.h"
#include "WENO.h"

#include <type_traits>

int main() {

  {
    /*
     * basis polynomial operations
     */
    constexpr std::array<double, 3> roots{0, 1, 2};
    constexpr Polynomial<3, double> poly(roots, 1);
    constexpr auto poly2 = poly.derivative();
    constexpr auto poly3 = poly2 * 2;
    constexpr auto poly4 = poly3.derivative<2>();
    constexpr auto poly5 = poly4.integrate();
    constexpr double integral = poly5.integrate(-0.5, 0.5);
  }

  {
    /*
     * calculate lagrangian interpolation
     */
    constexpr std::array<double, 3> x{-1, 0, 1};
    constexpr std::array<double, 3> f{1, 0, 1};
    constexpr Polynomial<2, double> poly_interpolation = lagrange_interpolate_to_polynomial(x, f);
    constexpr std::array<double, 3> interpolation_coefficients = lagrange_interpolation_coefficients(x, 4.0);
  }

  {
    /*
     * reconstruction from cell to face with uniformly spaced grids
     */
    constexpr std::array<double, 4> coeffs =
        calculate_reconstruction_coefficients_from_cell_to_face_for_uniformly_spaced_points<4, double>();
    // coeffs = {-1/12, 7/12, 7/12, -1/12}
    // this means u(1/2) = -1/12 * bar{u}(-1) + 7/12 * bar{u}(0) + 7/12 * bar{u}(1) - 1/12 * bar{u}(2)
  }

  {
    /*
     * calculate the coefficients for calculating derivatives at faces with four-points FD stencils.
     */
    constexpr std::array<double, 4> coeffs =
        calculate_interpolation_coefficients_from_cell_to_face_for_uniformly_spaced_points<4, double, 1>();
    // coeffs = {1/24, -9/8, 9/8, -1/24};
    // this means u'(1/2) = (1/dx) * (1/24 * u(-1) - 9/8 * u(0) + 9/8 * u(1) - 1/24 * u(2));

  }

  {
    /*
     * calculate the coefficients for calculating derivatives at faces with four-points FV stencils.
     */
    constexpr std::array<double, 4> coeffs =
        calculate_reconstruction_coefficients_from_cell_to_face_for_uniformly_spaced_points<4, double, 1>();
    // coeffs = {1/12, -5/4, 5/4, -1/12};
    // this means u'(1/2) = (1/dx) * (1/12 * bar{u}(-1) - 5/4 * bar{u}(0) + 5/4 * bar{u}(1) - 1/12 * bar{u}(2));
  }

  {
    /*
     * calculate the coefficients for calculating 1-st order derivatives at cell with five-points centered FD stencils.
     */
    constexpr std::array<double, 7> coeffs =
        calculate_interpolation_coefficients_from_cell_to_cell_for_uniformly_spaced_points<7, double, 1>();
    // coeffs = {1/12, -2/3, 0, 2/3, -1/12};
    // this means u'(0) = (1/dx) * (1/12 * u(-2) - 2/3 * u(-1)  + 2/3 * u(1) - 1/12 * u(2));

  }


  {
    /*
     * calculate the coefficients for calculating 2-nd order derivatives at cell with five-points centered FD stencils.
     */
    constexpr std::array<double, 5> coeffs =
        calculate_interpolation_coefficients_from_cell_to_cell_for_uniformly_spaced_points<5, double, 2>();
    // coeffs = {-1/12, 4/3, -5/2, 4/3, -1/12};
    // this means u''(0) = (1/dx^2) * (-1/12 * u(-2) + 4/3 * u(-1) - 5/2 * u(0) + 4/3 * u(1) - 1/12 * u(2));

  }

  {
    /*
     * WENO-JS coefficients
     */
    constexpr auto c = WENO_Uniform_Reconstruction_Coefficients<3, double>::c;
    constexpr auto d = WENO_Uniform_Reconstruction_Coefficients<3, double>::d;
    constexpr auto s = WENO_Uniform_Reconstruction_Coefficients<3, double>::s;
  }



  return 0;
}
