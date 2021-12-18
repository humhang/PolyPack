//
// Created by Hang Yu on 12/4/21.
//
#ifndef POLYPACK__LAGRANGEINTERPOLATION_H
#define POLYPACK__LAGRANGEINTERPOLATION_H

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
 * It is possible to output derivative phi_i^(r)(x_target) instead if OrderOfDerivative is given in argument parameter.
 */
template<size_t Npts, typename T, unsigned int OrderOfDerivative = 0>
static constexpr std::array<T, Npts>
lagrange_interpolation_coefficients(const std::array<T, Npts>& pts,
                                    const T target) {

  static_assert(Npts > 0);

  const std::array<Polynomial<Npts-1, T> , Npts> basis = lagrange_interpolation_basis_polynomial(pts);

  std::array<T, Npts> coefficients{};

  for(unsigned int pt = 0; pt < Npts; ++pt) {

    coefficients[pt] = basis[pt].template derivative<OrderOfDerivative>().evaluate(target);

  }

  return coefficients;

}




#endif //POLYPACK__LAGRANGEINTERPOLATION_H
