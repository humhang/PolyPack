//
// Created by hang on 8/21/23.
//

#ifndef POLYPACK__FIELD_H
#define POLYPACK__FIELD_H

#include <type_traits>

#if __cplusplus >= 202002L
#include <concepts>
#endif

#if !(__cplusplus >= 202002L)
// c++ 20 not available

template<typename>
std::false_type test_has_operatorPlus(...);

template<typename T>
std::true_type test_has_operatorPlus(std::remove_reference_t<decltype(std::declval<T>()+(std::declval<T>()))>*,
                                     std::enable_if_t<std::is_same_v<std::remove_cv_t<std::remove_reference_t<decltype(std::declval<T>()+(std::declval<T>()))>>, T>, void*>);

template<typename T>
struct has_operatorPlus {
  static constexpr bool value = decltype(test_has_operatorPlus<T>(nullptr, nullptr))::value;
};

template<typename>
std::false_type test_has_operatorMinus(...);

template<typename T>
std::true_type test_has_operatorMinus(std::remove_reference_t<decltype(std::declval<T>()-(std::declval<T>()))>*,
                                      std::enable_if_t<std::is_same_v<std::remove_cv_t<std::remove_reference_t<decltype(std::declval<T>()-(std::declval<T>()))>>, T>, void*>);

template<typename T>
struct has_operatorMinus {
  static constexpr bool value = decltype(test_has_operatorMinus<T>(nullptr, nullptr))::value;
};

template<typename>
std::false_type test_has_operatorMultiply(...);

template<typename T>
std::true_type test_has_operatorMultiply(std::remove_reference_t<decltype(std::declval<T>()*(std::declval<T>()))>*,
                                     std::enable_if_t<std::is_same_v<std::remove_cv_t<std::remove_reference_t<decltype(std::declval<T>()*(std::declval<T>()))>>, T>, void*>);

template<typename T>
struct has_operatorMultiply {
  static constexpr bool value = decltype(test_has_operatorMultiply<T>(nullptr, nullptr))::value;
};

template<typename>
std::false_type test_has_operatorDivide(...);

template<typename T>
std::true_type test_has_operatorDivide(std::remove_reference_t<decltype(std::declval<T>()/(std::declval<T>()))>*,
                                    std::enable_if_t<std::is_same_v<std::remove_cv_t<std::remove_reference_t<decltype(std::declval<T>()/(std::declval<T>()))>>, T>, void*>);

template<typename T>
struct has_operatorDivide {
  static constexpr bool value = decltype(test_has_operatorDivide<T>(nullptr, nullptr))::value;
};


template<typename T>
struct is_field {
  static constexpr bool value =
      has_operatorPlus<T>::value &&
      has_operatorMinus<T>::value &&
      has_operatorMultiply<T>::value &&
      has_operatorDivide<T>::value &&
      std::is_same_v<std::decay_t<T>, T> &&
      !std::is_integral_v<T>;
};

template<typename T>
inline constexpr bool field = is_field<T>::value;

#else
// c++ 20 available
template<typename T>
concept field =
    !std::integral<T> &&
    std::same_as<T, std::decay_t<T>> &&
    requires (const T& a,
              const T& b) {

  {a + b} -> std::same_as<T>;
  {a - b} -> std::same_as<T>;
  {a * b} -> std::same_as<T>;
  {a / b} -> std::same_as<T>;

};
#endif


#endif //POLYPACK__FIELD_H
