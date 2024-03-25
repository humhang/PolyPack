cmake_minimum_required(VERSION 3.20)
project(PolyPack CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
find_package(Boost)

if(Boost_FOUND)
    message(STATUS "boost found, boost::integral will be used!")
endif()


add_executable(PolyPack main.cpp)
if(Boost_FOUND)
    target_include_directories(PolyPack PRIVATE ${Boost_INCLUDE_DIRS})
    target_compile_definitions(PolyPack PRIVATE HAS_BOOST)
endif()
