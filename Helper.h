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

template<typename>
std::false_type test_has_operatorPlus(...);

template<typename T>
std::true_type test_has_operatorPlus(std::remove_reference_t<decltype(std::declval<T>()+std::declval<T>())>*,
                                     std::enable_if_t<std::is_same<std::remove_cv_t<std::remove_reference_t<decltype(std::declval<T>()+std::declval<T>())>>, T>::value, void*>);

template<typename T>
struct has_operatorPlus {
  static constexpr bool value = decltype(test_has_operatorPlus<T>(nullptr, nullptr))::value;
};

template<typename>
std::false_type test_has_operatorMultiply(...);

template<typename T>
std::true_type test_has_operatorMultiply(std::remove_reference_t<decltype(std::declval<T>()*std::declval<T>())>*,
                                     std::enable_if_t<std::is_same<std::remove_cv_t<std::remove_reference_t<decltype(std::declval<T>()*std::declval<T>())>>, T>::value, void*>);

template<typename T>
struct has_operatorMultiply {
  static constexpr bool value = decltype(test_has_operatorMultiply<T>(nullptr, nullptr))::value;
};


#endif //POLYPACK__HELPER_H
