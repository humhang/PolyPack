//
// Created by Hang Yu on 12/17/21.
//

#ifndef POLYPACK__WENO_H
#define POLYPACK__WENO_H

#include <cmath>

#include "Helper.h"
#include "Polynomial.h"
#include "LagrangeInterpolation.h"
#include "UniformGridInterpolationAndReconstruction.h"

/*
 * evaluate the reconstruction coefficients to face i+1/2 with Npts CELL AVERAGED values {u_i}
 * If Npts is odd, the stencils are
 * i-(Npts-1)/2+bias, ..., i, ..., i+(Npts-1)/2+bias,
 * If Npts is even, the stencils are
 * i-Npts/2+bias+1, ..., i, ..., i+Npts/2+bias,
 *
 * the smoothness of the polynomial reconstructed is a quadratic function of the stencils values {u_i}, which are
 * smoothness_indicator = \sum_i \sum_j s_{ij} u_i u_j
 * for i, j in the stencils
 * This function output the Npts by Npts coefficients lists
 */
template<unsigned int Npts, typename T>
static constexpr std::array<std::array<T, Npts>, Npts>
calculate_reconstruction_smoothness_indicator_coefficients_from_cell_for_uniformly_spaced_points(const int bias = 0) {

  static_assert(Npts > 0);

  /*
   * do lagrange interpolation for the integral of function that is to be reconstructed
   * the interpolation points are the faces, note that we have Npts+1 faces
   */
  std::array<T, Npts+1> faces{};

  if constexpr(is_odd(Npts)) {
    for(unsigned int k = 0; k < Npts+1; ++k) {
      faces[k] =  - (static_cast<T>(Npts)-1)/2+static_cast<T>(bias) + static_cast<T>(k) - half<T>::value;
    }
  } else {
    for(unsigned int k = 0; k < Npts+1; ++k) {
      faces[k] =  - static_cast<T>(Npts)/2+static_cast<T>(bias)+1 + static_cast<T>(k) - half<T>::value;
    }
  }
  const std::array<Polynomial<Npts, T>, Npts+1> basis_integral = lagrange_interpolation_basis_polynomial(faces);
  std::array<Polynomial<Npts-1, T>, Npts> basis {};

  basis[Npts-1] = basis_integral[Npts].derivative();
  for (unsigned int pt = Npts - 1; pt > 0; --pt) {
    basis[pt-1] = basis_integral[pt].derivative() + basis[pt];
  }

  std::array<std::array<T, Npts>, Npts> smoothness_indicator_coefficients {};
  for(unsigned int m = 0; m < Npts; ++m) {
    for(unsigned int n = 0; n < Npts; ++n) {

      T coeff = 0;
      constexpr_for<1, Npts, 1>([&basis, m, n, &coeff](auto i){
        coeff += (basis[m].template derivative<i>() * basis[n].template derivative<i>()).integrate(-half<T>::value, half<T>::value);
      });

      smoothness_indicator_coefficients[m][n] = coeff;

    }
  }

  return smoothness_indicator_coefficients;

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
    for(int k = -static_cast<int>(r)/2; k <= static_cast<int>(r)/2-1; ++k) {
      coefficients_array[k+static_cast<int>(r)/2] = calculate_reconstruction_coefficients_from_cell_to_face_for_uniformly_spaced_points<r, T>(k);
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

/*
 * 2r-1 order WENO reconstruction coefficients: s
 */
template<unsigned int r, typename T>
static constexpr std::array<std::array<std::array<T, r>, r>, r>
get_WENO_Uniform_Reconstruction_Coefficients_s() {

  std::array<std::array<std::array<T, r>, r>, r> coefficients_array {};
  if constexpr(is_odd(r)) {
    for(int k = -(static_cast<int>(r)-1)/2; k <= (static_cast<int>(r)-1)/2; ++k) {
      coefficients_array[k+(static_cast<int>(r)-1)/2] = calculate_reconstruction_smoothness_indicator_coefficients_from_cell_for_uniformly_spaced_points<r, T>(k);
    }
  } else {
    for(int k = -static_cast<int>(r)/2; k <= static_cast<int>(r)/2-1; ++k) {
      coefficients_array[k+static_cast<int>(r)/2] = calculate_reconstruction_smoothness_indicator_coefficients_from_cell_for_uniformly_spaced_points<r, T>(k);
    }
  }
  return coefficients_array;

}


template<unsigned int r, typename T>
struct WENO_Uniform_Reconstruction_Coefficients {

  static constexpr std::array<std::array<T, r>, r> c = get_WENO_Uniform_Reconstruction_Coefficients_c<r, T>();
  static constexpr std::array<T, r> d = get_WENO_Uniform_Reconstruction_Coefficients_d<r, T>(c);
  static constexpr std::array<std::array<std::array<T, r>, r>, r> s = get_WENO_Uniform_Reconstruction_Coefficients_s<r, T>();

};

template<unsigned int r, typename T>
static constexpr T do_WENO(const T array[2*r-1]) {

  constexpr T epsilon = 1e-6;

  T q[r], s[r], w[r];

  /*
   * candidates
   */
  for(int i = 0; i < r; ++i) {
    q[i] = 0.0;
    for(int j = 0; j < r; ++j) {
      q[i] += WENO_Uniform_Reconstruction_Coefficients<r, T>::c[i][j] * array[j+i];
    }
  }

  /*
   * smoothness indicator
   */
  for(int i = 0; i < r; ++i) {
    s[i] = 0.0;
    for(int j = 0; j < r; ++j) {
      for(int k = 0; k < r; ++k) {
        s[i] += WENO_Uniform_Reconstruction_Coefficients<r, T>::s[i][j][k] * array[j+i] * array[k+i];
      }
    }
  }

  /*
   * nonlinear weight
   */
  T wt = 0.0;
  for(int i = 0; i < r; ++i) {
    w[i] = WENO_Uniform_Reconstruction_Coefficients<r, T>::d[i] / std::pow(epsilon + s[i], 2);
    wt += w[i];
  }

  /*
   * output
   */
  T out = 0.0;
  for(int i = 0; i < r; ++i) {
    out += q[i] * w[i] / wt;
  }

  return out;

}


#endif //POLYPACK__WENO_H
