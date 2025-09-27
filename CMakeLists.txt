cmake_minimum_required(VERSION 3.16)
project(CETSP_project LANGUAGES CXX)

# === user options ===
option(USE_ABSEIL_HASH_SET "Use Abseil flat_hash_set (otherwise fallback / remove usage yourself)" ON)
option(DEBUG "Enable debug build and define DEBUG macro" OFF)

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

# Basic warnings / flags
if(MSVC)
  add_compile_options(/W4 /permissive-)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    # requested: compile with -march=native
    add_compile_options(-march=native)
  endif()
endif()

# === sources & include dirs ===
# Collect all .cpp files in src/ (adjust if you want explicit listing)
file(GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cpp")
if(NOT SRC_FILES)
  message(FATAL_ERROR "No source files found in src/ — make sure your .cpp files are in src/")
endif()

add_executable(${PROJECT_NAME} ${SRC_FILES})
target_include_directories(${PROJECT_NAME} PRIVATE "${CMAKE_SOURCE_DIR}/src")

# === Boost (headers used: geometry / rtree) ===
# Boost.Geometry and Boost.Geometry.Index are header-only for typical usage.
# Modern CMake installations provide Boost::boost; fall back to include dirs if necessary.
find_package(Boost REQUIRED)
if(TARGET Boost::boost)
  target_link_libraries(${PROJECT_NAME} PRIVATE Boost::boost)
else()
  target_include_directories(${PROJECT_NAME} PRIVATE ${Boost_INCLUDE_DIRS})
endif()

# === Abseil (optional, required if USE_ABSEIL_HASH_SET=ON) ===
if(USE_ABSEIL_HASH_SET)
  # Prefer the config package if Abseil was installed via CMake / package manager
  find_package(absl CONFIG REQUIRED)
  # flat_hash_set is part of absl::container and also requires hash implementation.
  # Different Abseil packaging can expose different targets; handle the common cases.
  if(TARGET absl::container AND TARGET absl::hash)
    target_link_libraries(${PROJECT_NAME} PRIVATE absl::container absl::hash)
  elseif(TARGET absl::flat_hash_map)  # defensive / unlikely, some installs export many targets
    target_link_libraries(${PROJECT_NAME} PRIVATE absl::flat_hash_map)
  elseif(TARGET abseil::abseil) # some package managers provide a single meta-target
    target_link_libraries(${PROJECT_NAME} PRIVATE abseil::abseil)
  else()
    message(FATAL_ERROR "Abseil found but expected CMake targets (absl::container and absl::hash or abseil::abseil) are missing. Please check your Abseil installation.")
  endif()
endif()

# === compile-time definitions ===
if(USE_ABSEIL_HASH_SET)
  target_compile_definitions(${PROJECT_NAME} PRIVATE USE_ABSEIL_HASH_SET=1)
else()
  target_compile_definitions(${PROJECT_NAME} PRIVATE USE_ABSEIL_HASH_SET=0)
endif()

if(DEBUG)
  target_compile_definitions(${PROJECT_NAME} PRIVATE DEBUG=1)
else()
  target_compile_definitions(${PROJECT_NAME} PRIVATE DEBUG=0)
endif()

# === helpful status output ===
message(STATUS "Project: ${PROJECT_NAME}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "USE_ABSEIL_HASH_SET = ${USE_ABSEIL_HASH_SET}")
message(STATUS "DEBUG = ${DEBUG}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
if(Boost_VERSION)
  message(STATUS "Found Boost (version ${Boost_VERSION})")
endif()
