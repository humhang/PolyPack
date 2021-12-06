cmake_minimum_required(VERSION 3.20)
project(PolyPack)

# C++17 is necessary because std::array::operator[] is constexpr only after C++17
# If we replace std::array with our own array type, then C++14 should be sufficient.
# I might also replace SFINAE with if constexpr so let's keep C++17.
set(CMAKE_CXX_STANDARD 17)

add_executable(PolyPack main.cpp)
