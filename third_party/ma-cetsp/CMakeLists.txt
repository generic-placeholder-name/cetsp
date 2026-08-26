cmake_minimum_required(VERSION 3.20)
project(MA_CETSP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

file(GLOB_RECURSE MA_CETSP_SOURCES CONFIGURE_DEPENDS "src/*.cpp")
add_executable(MA-CETSP ${MA_CETSP_SOURCES})
target_include_directories(MA-CETSP PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/Genetic"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/Genetic/Crossover"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/LocalSearch"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/Utils")

set(GUROBI_ROOT "$ENV{GUROBI_HOME}" CACHE PATH "Gurobi installation root")
find_path(GUROBI_INCLUDE_DIR
    NAMES gurobi_c++.h
    HINTS "${GUROBI_ROOT}/include")
find_library(GUROBI_LIBRARY
    NAMES gurobi130 gurobi120 gurobi110 gurobi100 gurobi95
    HINTS "${GUROBI_ROOT}/lib")

if(MSVC)
    find_library(GUROBI_CXX_RELEASE
        NAMES gurobi_c++md2022 gurobi_c++md2019 gurobi_c++md2017
        HINTS "${GUROBI_ROOT}/lib")
    find_library(GUROBI_CXX_DEBUG
        NAMES gurobi_c++mdd2022 gurobi_c++mdd2019 gurobi_c++mdd2017
        HINTS "${GUROBI_ROOT}/lib")
else()
    find_library(GUROBI_CXX_RELEASE
        NAMES gurobi_c++
        HINTS "${GUROBI_ROOT}/lib")
    set(GUROBI_CXX_DEBUG "${GUROBI_CXX_RELEASE}")
endif()

if(NOT GUROBI_INCLUDE_DIR OR NOT GUROBI_LIBRARY OR
   NOT GUROBI_CXX_RELEASE OR NOT GUROBI_CXX_DEBUG)
    message(FATAL_ERROR
        "Could not find Gurobi. Set GUROBI_ROOT to the installation root.")
endif()

target_include_directories(MA-CETSP PRIVATE "${GUROBI_INCLUDE_DIR}")
target_link_libraries(MA-CETSP PRIVATE
    "${GUROBI_LIBRARY}"
    optimized "${GUROBI_CXX_RELEASE}"
    debug "${GUROBI_CXX_DEBUG}")

if(MSVC)
    target_compile_options(MA-CETSP PRIVATE /W4 /permissive-)
else()
    target_compile_options(MA-CETSP PRIVATE -Wall -Wextra -Wpedantic)
endif()
