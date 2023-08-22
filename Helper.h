//
// Created by Hang on 12/4/21.
//

#ifndef POLYPACK__HELPER_H
#define POLYPACK__HELPER_H

#include <type_traits>
#include <array>

#if __cplusplus >= 202002L
#include <concepts>
#endif

#ifdef HAS_BOOST
#include "boost/rational.hpp"
#endif


static constexpr bool is_odd(unsigned int n) {
  return n % 2;
}

static constexpr unsigned int minusOrZero(unsigned int m, unsigned int n) {
  return m > n ? m - n : 0;
}

template<size_t N, typename T>
static constexpr std::array<T, N-1> deleteOne(const std::array<T, N>& array, const unsigned int index_to_delete){

  std::array<T, N-1> array_out{};
  unsigned int count = 0;
  for(unsigned int n = 0; n < N; ++n) {
    if(n != index_to_delete) {
      array_out[count] = array[n];
      count ++;
    }
  }

  return array_out;

}

template<size_t N, typename T>
static constexpr std::array<T, N> getArray_OneAtIndex(const unsigned int index) {

  std::array<T, N> array_out{};
  array_out[index] = 1;

  return array_out;

}

template<size_t N, typename T>
static constexpr std::array<T, N> getArray_OneAtAndAfterIndex(const unsigned int index) {

  std::array<T, N> array_out{};
  for(unsigned int k = index; k < N; ++k) {
    array_out[k] = 1;
  }

  return array_out;

}

template <auto begin, auto end, auto increment, class F>
constexpr void constexpr_for(F&& f) {

  if constexpr (begin < end) {
    f(std::integral_constant<decltype(begin), begin>());
    constexpr_for<begin + increment, end, increment>(f);
  }

}


template<typename T>
struct half {

};

template<>
struct half<float> {
  static constexpr float value = 0.5f;
};

template<>
struct half<double> {
  static constexpr double value = 0.5;
};

template<>
struct half<long double> {
  static constexpr long double value = 0.5L;
};

#ifdef HAS_BOOST
template<typename int_type>
struct half<boost::rational<int_type>> {
  static constexpr boost::rational<int_type> value = {1, 2};
};
#endif


#endif //POLYPACK__HELPER_H
