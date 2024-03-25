//
// Created by hang on 8/21/23.
//

#ifndef POLYPACK__FIELD_H
#define POLYPACK__FIELD_H

#include <type_traits>
#include <concepts>

template<typename T>
concept field =
    !std::integral<T> &&
    std::same_as<T, std::decay_t<T>> &&
    requires (T x) {

  {x + x} -> std::same_as<T>;
  {x - x} -> std::same_as<T>;
  {x * x} -> std::same_as<T>;
  {x / x} -> std::same_as<T>;

};


#endif //POLYPACK__FIELD_H
