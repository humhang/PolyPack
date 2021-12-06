//
// Created by Hang Yu on 12/4/21.
//
#ifndef POLYPACK__POLYINTERPOLATION_H
#define POLYPACK__POLYINTERPOLATION_H

#include <type_traits>
#include <array>

#include "Polynomial.h"
#include "Helper.h"

/*
 * give N points and values, output the Lagrange interpolation polynomial of order N-1
 */
template<size_t Npts, typename T>
static constexpr Polynomial<Npts-1, T>
lagrange_interpolate_to_polynomial(const std::array<T, Npts>& pts,
                                   const std::array<T, Npts>& vals) {

  static_assert(Npts > 0);
  /*
   * Implementation is shitty but since this is supposed to be evaluated at compiling time it should be okay.
   */
  Polynomial<Npts-1, T> poly;

  for(unsigned int pt = 0; pt < Npts; ++pt) {

    const Polynomial<Npts-1, T> this_poly(deleteOne(pts, pt), 1);

    const T coeff = vals[pt] / this_poly.evaluate(pts[pt]);

    poly = poly + coeff * this_poly;

  }

  return poly;

}

/*
 * given N points and the target point, output the coefficients of Lagrange interpolation so that
 * val_at_target = sum val_at_pts[i] * coefficients[i]
 */
template<size_t Npts, typename T>
static constexpr std::array<T, Npts>
lagrange_interpolation_coefficients(const std::array<T, Npts>& pts,
                                    const T target) {

  /*
   * yet another shitty implementation
   */
  std::array<T, Npts> coefficients{};

  for(unsigned int pt = 0; pt < Npts; ++pt) {

    const Polynomial<Npts - 1, T> poly = lagrange_interpolate_to_polynomial(pts, getArray_OneAtIndex<Npts, T>(pt));
    coefficients[pt] = poly.evaluate(target);
  }

  return coefficients;

}

/*
 * evaluate the interpolation coefficients to face i+1/2 with Npts cell values.
 * If Npts is odd, the stencils are
 * i-(Npts-1)/2+bias, ..., i, ..., i+(Npts-1)/2+bias,
 * If Npts is even, the stencils are
 * i-Npts/2+bias, ..., i, ..., i+Npts/2+bias-1,
 *
 */
template<unsigned int Npts, typename T>
static constexpr std::array<T, Npts>
calculate_interpolation_coefficients_from_cell_to_face_for_uniformly_spaced_points(const int bias = 0) {

  std::array<T, Npts> stencils{};

  if constexpr(is_odd(Npts)) {
    for(unsigned int k = 0; k < Npts; ++k) {
      stencils[k] =  - (static_cast<T>(Npts)-1)/2+static_cast<T>(bias) + static_cast<T>(k);
    }
  } else {
    for(unsigned int k = 0; k < Npts; ++k) {
      stencils[k] =  -static_cast<T>(Npts)/2+static_cast<T>(bias) + static_cast<T>(k);
    }
  }

  return lagrange_interpolation_coefficients(stencils, 0.5);

}


/*
 * evaluate the reconstruction coefficients to face i+1/2 with Npts CELL AVERAGED values.
 * If Npts is odd, the stencils are
 * i-(Npts-1)/2+bias, ..., i, ..., i+(Npts-1)/2+bias,
 * If Npts is even, the stencils are
 * i-Npts/2+bias, ..., i, ..., i+Npts/2+bias-1,
 *
 */
template<unsigned int Npts, typename T>
static constexpr std::array<T, Npts>
calculate_reconstruction_coefficients_from_cell_to_face_for_uniformly_spaced_points(const int bias = 0) {

  /*
   * do lagrange interpolation for the integral of function that is to be reconstructed
   * the interpolation points are the faces, note that we have Npts+1 faces
   */
  std::array<T, Npts+1> stencils{};

  if constexpr(is_odd(Npts)) {
    for(unsigned int k = 0; k < Npts+1; ++k) {
      stencils[k] =  - (static_cast<T>(Npts)-1)/2+static_cast<T>(bias) + static_cast<T>(k) - 0.5;
    }
  } else {
    for(unsigned int k = 0; k < Npts+1; ++k) {
      stencils[k] =  - static_cast<T>(Npts)/2+static_cast<T>(bias) + static_cast<T>(k) - 0.5;
    }
  }

  std::array<T, Npts> coefficients{};

  for(unsigned int pt = 0; pt < Npts; ++pt) {

    const Polynomial<Npts, T> poly_integral = lagrange_interpolate_to_polynomial(stencils, getArray_OneAtAndAfterIndex<Npts+1, T>(pt+1));
    const Polynomial<Npts-1, T> poly = poly_integral.derivative();

    coefficients[pt] = poly.evaluate(0.5);

  }

  return coefficients;

}

/*
 * 2r-1 order WENO reconstruction coefficients: c
 */
template<unsigned int r, typename T>
static constexpr std::array<std::array<T, r>, r>
get_WENO_Uniform_Reconstruction_Coefficients_c() {

  std::array<std::array<T, r>, r> coefficients_array {};
  if constexpr(is_odd(r)) {
    for(int k = -(static_cast<int>(r)-1)/2; k <= (static_cast<int>(r)-1)/2; ++k) {
      coefficients_array[k+(static_cast<int>(r)-1)/2] = calculate_reconstruction_coefficients_from_cell_to_face_for_uniformly_spaced_points<r, T>(k);
    }
  } else {
    for(int k = -static_cast<int>(r)/2+1; k <= static_cast<int>(r)/2; ++k) {
      coefficients_array[k+static_cast<int>(r)/2-1] = calculate_reconstruction_coefficients_from_cell_to_face_for_uniformly_spaced_points<r, T>(k);
    }
  }
  return coefficients_array;

}

/*
 * 2r-1 order WENO reconstruction coefficients: d
 *
 * This is by solving lower triangle system
 *
 * [c_[0][0]                                       ]  [d[0]   ]      [ b[0]   ]
 * |c_[0][1]     c_[1][0]                          |  |d[1]   |   =  | b[1]   |
 * |...          ...                               |  |...    |      |  ...   |
 * [c_[0][r-1]   c_[1][r-2]        c_[r-1][0]      ]  [d[r-1] ]      [ b[r-1] ]
 *
 *
 * where b is the coefficients of polynomial reconstruction with all 2r-1 points
 */
template<unsigned int r, typename T>
static constexpr std::array<T, r>
get_WENO_Uniform_Reconstruction_Coefficients_d(const std::array<std::array<T, r>, r>& c) {

  constexpr std::array<T, 2*r-1> b = calculate_reconstruction_coefficients_from_cell_to_face_for_uniformly_spaced_points<2*r-1, T>(0);

  std::array<T, r> d{};

  for(unsigned int k = 0; k < r; ++k) {
    T res = 0;
    for(unsigned int m = 0; m < k; ++m) {
      res += d[m] * c[m][k-m];
    }
    d[k] = (b[k] - res) / c[k][0];
  }

  return d;

}


template<unsigned int r, typename T>
struct WENO_Uniform_Reconstruction_Coefficients {
  static constexpr std::array<std::array<T, r>, r> c = get_WENO_Uniform_Reconstruction_Coefficients_c<r, T>();
  static constexpr std::array<T, r> d = get_WENO_Uniform_Reconstruction_Coefficients_d<r, T>(c);

};




#endif //POLYPACK__POLYINTERPOLATION_H
