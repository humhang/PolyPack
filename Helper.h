//
// Created by Hang on 12/4/21.
//

#ifndef POLYPACK__HELPER_H
#define POLYPACK__HELPER_H

#include <type_traits>
#include <array>

static constexpr bool is_odd(unsigned int n) {
  return n % 2;
}

static constexpr unsigned int minusOrZero(unsigned int m, unsigned int n) {
  unsigned int k = 0;
  if (m > n) {
    k = m - n;
  }
  return k;
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

template <auto Start, auto End, auto Inc, class F>
constexpr void constexpr_for(F&& f) {

  if constexpr (Start < End) {
    f(std::integral_constant<decltype(Start), Start>());
    constexpr_for<Start + Inc, End, Inc>(f);
  }

}



#endif //POLYPACK__HELPER_H
