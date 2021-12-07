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
 * give N points, output N lagrange interpolation basis polynomials phi_i(x) of degree N-1
 * where phi_i(x_j) = delta_ij
 */
template<size_t Npts, typename T>
static constexpr std::array<Polynomial<Npts-1, T> , Npts>
lagrange_interpolation_basis_polynomial(const std::array<T, Npts>& pts) {

  static_assert(Npts > 0);
  std::array<Polynomial<Npts-1, T> , Npts> basis{};

  for(unsigned int pt = 0; pt < Npts; ++pt) {

    const Polynomial<Npts-1, T> this_poly(deleteOne(pts, pt), 1);
    const T coeff = 1 / this_poly.evaluate(pts[pt]);
    basis[pt] = coeff * this_poly;

  }

  return basis;

}


/*
 * give N points and values, output the Lagrange interpolation polynomial of order N-1
 * i.e., p(x) = \sum val_i * phi_i(x)
 */
template<size_t Npts, typename T>
static constexpr Polynomial<Npts-1, T>
lagrange_interpolate_to_polynomial(const std::array<T, Npts>& pts,
                                   const std::array<T, Npts>& vals) {

  static_assert(Npts > 0);

  const std::array<Polynomial<Npts-1, T> , Npts> basis = lagrange_interpolation_basis_polynomial(pts);

  Polynomial<Npts-1, T> poly;

  for(unsigned int pt = 0; pt < Npts; ++pt) {

    poly = poly + vals[pt] * basis[pt];

  }

  return poly;

}

/*
 * given N points and the target point, output the coefficients of Lagrange interpolation so that
 * p(x_target) = \sum val_i * phi_i(x_target)
 * the coefficients are phi_i(x_target)
 */
template<size_t Npts, typename T>
static constexpr std::array<T, Npts>
lagrange_interpolation_coefficients(const std::array<T, Npts>& pts,
                                    const T target) {

  static_assert(Npts > 0);

  const std::array<Polynomial<Npts-1, T> , Npts> basis = lagrange_interpolation_basis_polynomial(pts);

  std::array<T, Npts> coefficients{};

  for(unsigned int pt = 0; pt < Npts; ++pt) {

    coefficients[pt] = basis[pt].evaluate(target);

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

  static_assert(Npts > 0);

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
      faces[k] =  - static_cast<T>(Npts)/2+static_cast<T>(bias) + static_cast<T>(k) - 0.5;
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

    coefficients[pt] = basis[pt].evaluate(0.5);

  }

  return coefficients;

}


/*
 * evaluate the reconstruction coefficients to face i+1/2 with Npts CELL AVERAGED values {u_i}
 * If Npts is odd, the stencils are
 * i-(Npts-1)/2+bias, ..., i, ..., i+(Npts-1)/2+bias,
 * If Npts is even, the stencils are
 * i-Npts/2+bias, ..., i, ..., i+Npts/2+bias-1,
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
      faces[k] =  - (static_cast<T>(Npts)-1)/2+static_cast<T>(bias) + static_cast<T>(k) - 0.5;
    }
  } else {
    for(unsigned int k = 0; k < Npts+1; ++k) {
      faces[k] =  - static_cast<T>(Npts)/2+static_cast<T>(bias) + static_cast<T>(k) - 0.5;
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
        coeff += (basis[m].template derivative<i>() * basis[n].template derivative<i>()).integrate(-1.0/2.0, 1.0/2.0);
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
    for(int k = -static_cast<int>(r)/2+1; k <= static_cast<int>(r)/2; ++k) {
      coefficients_array[k+static_cast<int>(r)/2-1] = calculate_reconstruction_smoothness_indicator_coefficients_from_cell_for_uniformly_spaced_points<r, T>(k);
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




#endif //POLYPACK__POLYINTERPOLATION_H
