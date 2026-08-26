cmake_minimum_required(VERSION 3.16)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

project(CETSP_project LANGUAGES CXX)

include(CTest)

# === user options ===
if(DEFINED USE_ABSEIL_HASH_SET)
  if(NOT DEFINED CETSP_SET_BACKEND)
    if(USE_ABSEIL_HASH_SET)
      set(CETSP_SET_BACKEND "ABSEIL" CACHE STRING
        "Set backend: AUTO, ABSEIL, UNORDERED, or ORDERED")
    else()
      set(CETSP_SET_BACKEND "ORDERED" CACHE STRING
        "Set backend: AUTO, ABSEIL, UNORDERED, or ORDERED")
    endif()
    message(DEPRECATION
      "USE_ABSEIL_HASH_SET is deprecated; mapped it to "
      "CETSP_SET_BACKEND=${CETSP_SET_BACKEND}.")
  else()
    message(DEPRECATION
      "USE_ABSEIL_HASH_SET is deprecated and ignored because "
      "CETSP_SET_BACKEND is already set.")
  endif()
  unset(USE_ABSEIL_HASH_SET CACHE)
endif()

set(CETSP_SET_BACKEND "AUTO" CACHE STRING
  "Set backend: AUTO, ABSEIL, UNORDERED, or ORDERED")
set_property(CACHE CETSP_SET_BACKEND PROPERTY STRINGS
  AUTO ABSEIL UNORDERED ORDERED)
string(TOUPPER "${CETSP_SET_BACKEND}" CETSP_REQUESTED_SET_BACKEND)
set(CETSP_VALID_SET_BACKENDS AUTO ABSEIL UNORDERED ORDERED)
if(NOT CETSP_REQUESTED_SET_BACKEND IN_LIST CETSP_VALID_SET_BACKENDS)
  message(FATAL_ERROR
    "Invalid CETSP_SET_BACKEND='${CETSP_SET_BACKEND}'. "
    "Expected one of: AUTO, ABSEIL, UNORDERED, ORDERED.")
endif()
option(DEBUG "Enable debug build and define DEBUG macro" OFF)
option(CETSP_BUILD_POSTPROCESSING
  "Build the optional external-solver tour postprocessors" OFF)
option(CETSP_BUILD_MA_CETSP_INTEGRATION
  "Build the MA-CETSP seed exporter" OFF)

# Default build type when single-config generators are used
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  if(DEBUG)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type" FORCE)
  else()
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
  endif()
endif()

# === language settings ===
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

function(cetsp_configure_target target)
  target_compile_features(${target} PRIVATE cxx_std_20)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /permissive-)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
      target_compile_options(${target} PRIVATE -march=native)
    endif()
  endif()
endfunction()

# Keep the algorithm independent of the command-line application so unit tests
# and future frontends can link it directly.
add_library(cetsp_core STATIC
  src/solver.cpp
  src/circle_geometry.cpp
  src/merge.cpp
  src/reconstruct.cpp
)
add_library(cetsp::core ALIAS cetsp_core)
target_include_directories(cetsp_core
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
)
target_compile_features(cetsp_core PUBLIC cxx_std_20)
cetsp_configure_target(cetsp_core)

add_executable(${PROJECT_NAME} app/main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE cetsp_core)
cetsp_configure_target(${PROJECT_NAME})

# === Boost (headers used: geometry / rtree) ===
# Boost.Geometry and Boost.Geometry.Index are header-only for typical usage.
# Modern CMake installations provide Boost::boost; fall back to include dirs if necessary.
find_package(Boost REQUIRED)
find_package(Threads REQUIRED)
target_link_libraries(cetsp_core PRIVATE Threads::Threads)
if(TARGET Boost::boost)
  target_link_libraries(cetsp_core PUBLIC Boost::boost)
else()
  target_include_directories(cetsp_core PUBLIC ${Boost_INCLUDE_DIRS})
endif()

if(CETSP_BUILD_POSTPROCESSING)
  add_library(cetsp_postprocess STATIC
    src/postprocess/process_support.cpp
    src/postprocess/lkh.cpp
    src/postprocess/socp.cpp
    src/postprocess/pipeline.cpp
  )
  add_library(cetsp::postprocess ALIAS cetsp_postprocess)
  target_include_directories(cetsp_postprocess
    PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    PRIVATE
      "${CMAKE_CURRENT_SOURCE_DIR}/src/postprocess"
  )
  target_link_libraries(cetsp_postprocess PUBLIC cetsp_core)
  cetsp_configure_target(cetsp_postprocess)
endif()

if(CETSP_BUILD_MA_CETSP_INTEGRATION)
  add_executable(
    CETSP_ma_cetsp_seed_export
    tools/ma_cetsp/seed_export.cpp
  )
  target_link_libraries(CETSP_ma_cetsp_seed_export PRIVATE cetsp_core)
  cetsp_configure_target(CETSP_ma_cetsp_seed_export)
  if(CETSP_BUILD_POSTPROCESSING)
    add_executable(
      CETSP_ma_cetsp_result_check
      tools/ma_cetsp/result_check.cpp
    )
    target_link_libraries(
      CETSP_ma_cetsp_result_check PRIVATE cetsp_postprocess)
    cetsp_configure_target(CETSP_ma_cetsp_result_check)
  endif()
endif()

# === internal set backend ===
set(CETSP_RESOLVED_SET_BACKEND "${CETSP_REQUESTED_SET_BACKEND}")
set(cetsp_abseil_targets)

function(cetsp_select_abseil_targets output)
  if(TARGET absl::flat_hash_set)
    set(${output} absl::flat_hash_set PARENT_SCOPE)
  elseif(TARGET absl::container AND TARGET absl::hash)
    set(${output} absl::container absl::hash PARENT_SCOPE)
  elseif(TARGET abseil::abseil)
    set(${output} abseil::abseil PARENT_SCOPE)
  else()
    set(${output} "" PARENT_SCOPE)
  endif()
endfunction()

if(CETSP_REQUESTED_SET_BACKEND STREQUAL "AUTO" OR
   CETSP_REQUESTED_SET_BACKEND STREQUAL "ABSEIL")
  cetsp_select_abseil_targets(cetsp_abseil_targets)
  if(NOT cetsp_abseil_targets)
    find_package(absl CONFIG QUIET)
    cetsp_select_abseil_targets(cetsp_abseil_targets)
  endif()

  if(cetsp_abseil_targets)
    set(CETSP_RESOLVED_SET_BACKEND "ABSEIL")
  elseif(CETSP_REQUESTED_SET_BACKEND STREQUAL "ABSEIL")
    message(FATAL_ERROR
      "CETSP_SET_BACKEND=ABSEIL was requested, but no supported Abseil "
      "target was found. Install Abseil or select AUTO, UNORDERED, or ORDERED.")
  else()
    set(CETSP_RESOLVED_SET_BACKEND "UNORDERED")
    message(STATUS
      "No supported Abseil target found; AUTO is using std::unordered_set")
  endif()
endif()

if(CETSP_RESOLVED_SET_BACKEND STREQUAL "ABSEIL")
  target_link_libraries(cetsp_core PRIVATE ${cetsp_abseil_targets})
  target_compile_definitions(cetsp_core PRIVATE CETSP_HASH_SET_ABSEIL=1)
elseif(CETSP_RESOLVED_SET_BACKEND STREQUAL "UNORDERED")
  target_compile_definitions(cetsp_core PRIVATE CETSP_HASH_SET_UNORDERED=1)
elseif(CETSP_RESOLVED_SET_BACKEND STREQUAL "ORDERED")
  target_compile_definitions(cetsp_core PRIVATE CETSP_HASH_SET_ORDERED=1)
endif()

if(DEBUG)
  target_compile_definitions(cetsp_core PRIVATE DEBUG=1)
else()
  target_compile_definitions(cetsp_core PRIVATE DEBUG=0)
endif()

if(BUILD_TESTING)
  # This target is deliberately a consumer test: it receives no access to src/.
  add_executable(CETSP_tests tests/CETSPTests.cpp)
  target_link_libraries(CETSP_tests PRIVATE cetsp_core)
  cetsp_configure_target(CETSP_tests)
  add_test(NAME CETSP_tests COMMAND CETSP_tests)

  # Focused white-box regressions may include private implementation headers.
  add_executable(CETSP_internal_tests tests/CETSPInternalTests.cpp)
  target_include_directories(CETSP_internal_tests PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
  )
  target_link_libraries(CETSP_internal_tests PRIVATE cetsp_core)
  cetsp_configure_target(CETSP_internal_tests)
  add_test(NAME CETSP_internal_tests COMMAND CETSP_internal_tests)

  add_executable(LKH_parser_tests tests/LKHParserTests.cpp)
  target_include_directories(LKH_parser_tests PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src/postprocess"
  )
  if(TARGET Boost::boost)
    target_link_libraries(LKH_parser_tests PRIVATE Boost::boost)
  else()
    target_include_directories(LKH_parser_tests PRIVATE ${Boost_INCLUDE_DIRS})
  endif()
  cetsp_configure_target(LKH_parser_tests)
  add_test(NAME LKH_parser_tests COMMAND LKH_parser_tests)

  # This is intentionally not a CTest test: it is stochastic and comparatively
  # expensive. Run it explicitly when checking solution quality.
  add_executable(CETSP_benchmark tests/CETSPBenchmark.cpp)
  target_link_libraries(CETSP_benchmark PRIVATE cetsp_core)
  cetsp_configure_target(CETSP_benchmark)

  if(CETSP_BUILD_POSTPROCESSING)
    add_executable(
      CETSP_postprocess_tests
      tests/CETSPPostprocessTests.cpp
    )
    target_link_libraries(CETSP_postprocess_tests PRIVATE cetsp_postprocess)
    cetsp_configure_target(CETSP_postprocess_tests)
    add_test(NAME CETSP_postprocess_tests COMMAND CETSP_postprocess_tests)

    add_executable(
      CETSP_postprocess_benchmark
      tests/CETSPPostprocessBenchmark.cpp
    )
    target_link_libraries(
      CETSP_postprocess_benchmark PRIVATE cetsp_postprocess)
    cetsp_configure_target(CETSP_postprocess_benchmark)
  endif()

  if(CETSP_RESOLVED_SET_BACKEND STREQUAL "ORDERED")
    add_test(
      NAME CETSP_seed_reproducibility
      COMMAND ${CMAKE_COMMAND}
        "-DBENCHMARK=$<TARGET_FILE:CETSP_benchmark>"
        "-DINSTANCE=${CMAKE_CURRENT_SOURCE_DIR}/data/mennell/instances/bubbles4.txt"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/CheckReproducibility.cmake"
    )
  endif()
endif()

# === helpful status output ===
message(STATUS "Project: ${PROJECT_NAME}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "CETSP_SET_BACKEND requested = ${CETSP_REQUESTED_SET_BACKEND}")
message(STATUS "CETSP_SET_BACKEND resolved = ${CETSP_RESOLVED_SET_BACKEND}")
message(STATUS "DEBUG = ${DEBUG}")
message(STATUS "CETSP_BUILD_POSTPROCESSING = ${CETSP_BUILD_POSTPROCESSING}")
message(STATUS
  "CETSP_BUILD_MA_CETSP_INTEGRATION = ${CETSP_BUILD_MA_CETSP_INTEGRATION}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
if(Boost_VERSION)
  message(STATUS "Found Boost (version ${Boost_VERSION})")
endif()
