//
// Created by Hang Yu on 12/17/21.
//

#ifndef POLYPACK__UNIFORMGRIDINTERPOLATIONANDRECONSTRUCTION_H
#define POLYPACK__UNIFORMGRIDINTERPOLATIONANDRECONSTRUCTION_H

#include "Helper.h"
#include "Polynomial.h"
#include "LagrangeInterpolation.h"



/*
 * evaluate the interpolation coefficients to face i+1/2 with Npts cell values.
 * If Npts is odd, the stencils are
 * i-(Npts-1)/2+bias, ..., i, ..., i+(Npts-1)/2+bias,
 * If Npts is even, the stencils are
 * i-Npts/2+bias+1, ..., i, ..., i+Npts/2+bias,
 * It is possible to output the coefficients for derivative by specifying template argument OrderOfDerivative
 */
template<size_t Npts, typename T, unsigned int OrderOfDerivative = 0>
static constexpr std::array<T, Npts>
calculate_interpolation_coefficients_from_cell_to_face_for_uniformly_spaced_points(const int bias = 0) {

  static_assert(Npts > 0);

  std::array<T, Npts> stencils{};

  if constexpr(is_odd(Npts)) {
    for(unsigned int k = 0; k < Npts; ++k) {
      stencils[k] =  - (static_cast<T>(Npts)-1)/2+static_cast<T>(bias) + static_cast<T>(k);
    }
  } else {
    for(unsigned int k = 0; k < Npts; ++k) {
      stencils[k] =  -static_cast<T>(Npts)/2+static_cast<T>(bias)+1 + static_cast<T>(k);
    }
  }

  return lagrange_interpolation_coefficients<Npts, T, OrderOfDerivative>(stencils, 0.5);

}


/*
 * evaluate the interpolation coefficients to cell i with Npts cell values.
 * If Npts is odd, the stencils are
 * i-(Npts-1)/2+bias, ..., i, ..., i+(Npts-1)/2+bias,
 * If Npts is even, the stencils are
 * i-Npts/2+bias+1, ..., i, ..., i+Npts/2+bias,
 * It is possible to output the coefficients for derivative by specifying template argument OrderOfDerivative
 */
template<size_t Npts, typename T, unsigned int OrderOfDerivative = 0>
static constexpr std::array<T, Npts>
calculate_interpolation_coefficients_from_cell_to_cell_for_uniformly_spaced_points(const int bias = 0) {

  static_assert(Npts > 0);

  std::array<T, Npts> stencils{};

  if constexpr(is_odd(Npts)) {
    for(unsigned int k = 0; k < Npts; ++k) {
      stencils[k] =  - (static_cast<T>(Npts)-1)/2+static_cast<T>(bias) + static_cast<T>(k);
    }
  } else {
    for(unsigned int k = 0; k < Npts; ++k) {
      stencils[k] =  -static_cast<T>(Npts)/2+static_cast<T>(bias)+1 + static_cast<T>(k);
    }
  }

  return lagrange_interpolation_coefficients<Npts, T, OrderOfDerivative>(stencils, 0.0);

}


/*
 * evaluate the reconstruction coefficients to face i+1/2 with Npts CELL AVERAGED values.
 * If Npts is odd, the stencils are
 * i-(Npts-1)/2+bias, ..., i, ..., i+(Npts-1)/2+bias,
 * If Npts is even, the stencils are
 * i-Npts/2+bias+1, ..., i, ..., i+Npts/2+bias,
 * It is possible to output the coefficients for derivative by specifying template argument OrderOfDerivative
 */
template<unsigned int Npts, typename T, unsigned int OrderOfDerivative = 0>
static constexpr std::array<T, Npts>
calculate_reconstruction_coefficients_from_cell_to_face_for_uniformly_spaced_points(const int bias = 0) {

  static_assert(Npts > 0);

  /*
   * do lagrange interpolation for the integral of function that is to be reconstructed
   * the interpolation points are the faces, note that we have Npts+1 faces
   */
  std::array<T, Npts+1> faces{};

  if constexpr(is_odd(Npts)) {
    for(unsigned int k = 0; k < Npts+1; ++k) {
      faces[k] =  - (static_cast<T>(Npts)-1)/2+static_cast<T>(bias) + static_cast<T>(k) - 0.5;
    }
  } else {
    for(unsigned int k = 0; k < Npts+1; ++k) {
      faces[k] =  - static_cast<T>(Npts)/2+static_cast<T>(bias)+1 + static_cast<T>(k) - 0.5;
    }
  }
  const std::array<Polynomial<Npts, T>, Npts+1> basis_integral = lagrange_interpolation_basis_polynomial(faces);
  std::array<Polynomial<Npts-1, T>, Npts> basis {};

  basis[Npts-1] = basis_integral[Npts].derivative();
  for (unsigned int pt = Npts - 1; pt > 0; --pt) {
    basis[pt-1] = basis_integral[pt].derivative() + basis[pt];
  }

  std::array<T, Npts> coefficients{};

  for(unsigned int pt = 0; pt < Npts; ++pt) {

    coefficients[pt] = basis[pt].template derivative<OrderOfDerivative>().evaluate(0.5);

  }

  return coefficients;

}


#endif //POLYPACK__UNIFORMGRIDINTERPOLATIONANDRECONSTRUCTION_H
