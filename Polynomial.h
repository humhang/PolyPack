//
// Created by Hang on 12/5/21.
//

#ifndef POLYPACK__POLYNOMIAL_H
#define POLYPACK__POLYNOMIAL_H

#include <type_traits>
#include <array>
#include <algorithm>

#include "Helper.h"

/*
 * N-deg polynomial over field T
 */
template<unsigned int N, typename T>
class Polynomial {

  static_assert(is_field<T>::value);

public:

  /*
   * friends
   */
  template<unsigned int, typename>
  friend class Polynomial;

  /*
   * friends operators
   */
  template<unsigned int N_, typename T_>
  friend constexpr Polynomial<N_, T_> operator*(const T_ &scalar, const Polynomial<N_, T_> &poly);

  /*
   * default constructor
   */
  constexpr Polynomial() : coeffs({}) {};

  /*
   * constructor from coefficients
   */
  explicit constexpr Polynomial(const std::array<T, N + 1> &_coeffs) : coeffs(_coeffs) {};

  /*
   * trivial copy constructor
   */
  constexpr Polynomial(const Polynomial<N, T> &poly) = default;

  /*
   * trivial move constructor
   */
  constexpr Polynomial(Polynomial<N, T> &&poly) noexcept = default;

  /*
   * trivial copy assignment operator
   */
  Polynomial<N, T> &operator=(const Polynomial<N, T> &poly) = default;

  /*
   * trivial move assignment operator
   */
  Polynomial<N, T> &operator=(Polynomial<N, T> &&poly) noexcept = default;

  /*
   * constructor from roots and the leading coefficient (highest degree), for degree > 0
   */
  template<unsigned int N_ = N,
      typename std::enable_if<(N_ > 0), void *>::type = nullptr>
  constexpr Polynomial(const std::array<T, N> &roots, T leading_term) : coeffs({}) {

    /*
     * calculate the coefficients recursively
     */
    std::array<T, N - 1> roots_poly_one_degree_less({});
    for (unsigned int n = 0; n < N - 1; ++n) {
      roots_poly_one_degree_less[n] = roots[n];
    }
    Polynomial<N - 1, T> poly_one_order_less(roots_poly_one_degree_less, leading_term);

    Polynomial<1, T> poly1(std::array<T, 2>{-roots[N - 1], 1});

    coeffs = std::move((poly_one_order_less * poly1).coeffs);

  }

  /*
   * constructor from roots and the leading coefficient (the highest degree), for degree == 0
   */
  template<unsigned int N_ = N,
      typename std::enable_if<(N_ == 0), void *>::type = nullptr>
  constexpr Polynomial(const std::array<T, N> &roots, T leading_term) : coeffs({}) {

    coeffs[0] = leading_term;

  }

  /*
   * polynomial multiplication
   */
  template<unsigned int M>
  constexpr Polynomial<N + M, T>
  operator*(const Polynomial<M, T> &poly2) const {

    Polynomial<N + M, T> poly_out{};

    for (unsigned int m = 0; m <= N; ++m) {
      for (unsigned int n = 0; n <= M; ++n) {
        const unsigned int k = m + n;
        poly_out.coeffs[k] += coeffs[m] * poly2.coeffs[n];
      }
    }

    return poly_out;
  }

  /*
   * polynomial addition
   */
  template<unsigned int M>
  constexpr Polynomial<std::max(N, M), T>
  operator+(const Polynomial<M, T> &poly2) const {

    constexpr unsigned int degree_out = std::max(N, M);
    Polynomial<degree_out, T> poly_out;
    for (unsigned int m = 0; m <= degree_out; m++) {
      if (m <= N)
        poly_out.coeffs[m] += coeffs[m];
      if (m <= M)
        poly_out.coeffs[m] += poly2.coeffs[m];
    }

    return poly_out;

  }

  /*
   * polynomial minus
   */
  template<unsigned int M>
  constexpr Polynomial<std::max(N, M), T>
  operator-(const Polynomial<M, T> &poly2) const {

    constexpr unsigned int degree_out = std::max(N, M);
    Polynomial<degree_out, T> poly_out;
    for (unsigned int m = 0; m <= degree_out; m++) {
      if (m <= N)
        poly_out.coeffs[m] += coeffs[m];
      if (m <= M)
        poly_out.coeffs[m] -= poly2.coeffs[m];
    }

    return poly_out;
  }

  /*
   * polynomial multiplication with scalar
   */
  constexpr Polynomial<N, T>
  operator*(const T scalar) const {

    Polynomial<N, T> poly_out;
    for (unsigned int n = 0; n <= N; ++n) {
      poly_out.coeffs[n] = scalar * coeffs[n];
    }

    return poly_out;
  }

  /*
   * evaluate recursively with Horner's method, for polynomial deg > 0
   */
  template<typename... Dummy, unsigned int N_ = N>
  constexpr typename std::enable_if<(N_ > 0), T>::type
  evaluate(const T t) const {
    static_assert(sizeof...(Dummy) == 0, "Do not specify template arguments!");

    Polynomial<N_ - 1, T> poly_one_less(deleteOne(coeffs, 0));
    T out = poly_one_less.evaluate(t) * t + coeffs[0];
    return out;
  }

  /*
   * evaluate recursively with Horner's method, for polynomial deg == 0
   */
  template<typename... Dummy, unsigned int N_ = N>
  constexpr typename std::enable_if<(N_ == 0), T>::type
  evaluate(const T /*t*/) const {
    static_assert(sizeof...(Dummy) == 0, "Do not specify template arguments!");

    return coeffs[0];
  }

  /*
   * yet another evaluation wrapper
   */
  constexpr T
  operator()(const T t) const {
    return evaluate(t);
  }

  /*
   * taking r-th derivative for polynomial, if r > 1 and degree > 0
   */
  template<unsigned int r = 1, typename... Dummy, unsigned int N_ = N>
  constexpr typename std::enable_if<((N_ > 0) && (r > 1)), Polynomial<minusOrZero(N_, r), T> >::type
  derivative() const {

    static_assert(sizeof...(Dummy) == 0, "Do not specify template arguments other than the number of derivative!");

    const auto poly_derivative_r_minus_one_times = derivative<r-1>();

    return poly_derivative_r_minus_one_times.derivative();

  }

  /*
   * taking r-th derivative for polynomial, if r = 1 and degree > 0
   */
  template<unsigned int r = 1, typename... Dummy, unsigned int N_ = N>
  constexpr typename std::enable_if<((N_ > 0) && (r == 1)), Polynomial<minusOrZero(N_, r), T> >::type
  derivative() const {

    static_assert(sizeof...(Dummy) == 0, "Do not specify template arguments other than the number of derivative!");

    Polynomial<N_ - 1, T> poly_out;
    for (unsigned int n = 0; n <= N_ - 1; ++n) {
      poly_out.coeffs[n] = static_cast<T>(n + 1) * coeffs[n + 1];
    }
    return poly_out;
  }

  /*
   * taking r-th derivative for polynomial, if r >= 0 and degree = 0
   */
  template<unsigned int r = 1, typename... Dummy, unsigned int N_ = N>
  constexpr typename std::enable_if<((N_ == 0) && (r >= 1)), Polynomial<0, T> >::type
  derivative() const {

    static_assert(sizeof...(Dummy) == 0, "Do not specify template arguments other than the number of derivative!");

    Polynomial<0, T> poly_out;
    poly_out.coeffs[0] = 0;
    return poly_out;
  }

  /*
   * taking r-th derivative for polynomial, if r = 0
   */
  template<unsigned int r = 1, typename... Dummy, unsigned int N_ = N>
  constexpr typename std::enable_if<((r == 0)), Polynomial<N_, T> >::type
  derivative() const {

    static_assert(sizeof...(Dummy) == 0, "Do not specify template arguments other than the number of derivative!");

    return *this;
  }


  /*
   * indefinite integral of the polynomial
   * TODO: maybe we can make this recursive just like others;
   */
  constexpr Polynomial<N + 1, T>
  integrate(const T constant = 0) const {

    Polynomial<N + 1, T> poly_out;
    poly_out.coeffs[0] = constant;
    for (unsigned int n = 1; n <= N + 1; ++n) {
      poly_out.coeffs[n] = coeffs[n - 1] / static_cast<T>(n);
    }

    return poly_out;

  }

  /*
   * definite integral of the polynomial
   */
  constexpr T
  integrate(const T lower, const T upper) const {

    const Polynomial<N + 1, T> poly = integrate();

    return poly(upper) - poly(lower);

  }

private:

  /*
   * p(x) = sum_{k=0}^N coeffs[k] * x^k
   */
  std::array<T, N + 1> coeffs;

};

template<unsigned int N, typename T>
constexpr Polynomial<N, T> operator*(const T &scalar, const Polynomial<N, T> &poly) {
  Polynomial<N, T> poly_out;
  for (unsigned int n = 0; n <= N; ++n) {
    poly_out.coeffs[n] = scalar * poly.coeffs[n];
  }

  return poly_out;
}



/*
namespace {

using TestType = Polynomial<3, double>;
static_assert(std::is_literal_type<TestType>::value);
static_assert(std::is_constructible<TestType>::value);
static_assert(std::is_move_constructible<TestType>::value);
static_assert(std::is_copy_constructible<TestType>::value);
static_assert(std::is_trivially_move_constructible<TestType>::value);
static_assert(std::is_trivially_copy_constructible<TestType>::value);

}
*/

#endif //POLYPACK__POLYNOMIAL_H
