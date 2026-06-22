cmake_minimum_required(VERSION 3.16)
project(AeroAnalyzerPro CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Native tuning is opt-out: handy for a personal build, breaks portability.
option(AERO_USE_NATIVE "Enable -march=native / /arch:AVX2" ON)

# Compiler flags. DELIBERATELY NO -ffast-math: it implies -ffinite-math-only,
# which lets the compiler delete the NaN/Inf watchdog checks this code relies on,
# and it makes the trim/iterative tolerances non-deterministic. See README.
if(MSVC)
  add_compile_options(/O2 /fp:precise /EHsc /W3 /utf-8)
  if(AERO_USE_NATIVE)
    add_compile_options(/arch:AVX2)
  endif()
else()
  add_compile_options(-O3 -Wall -Wextra -fno-math-errno)
  if(AERO_USE_NATIVE)
    add_compile_options(-march=native)
  endif()
endif()

find_package(OpenMP QUIET)

# Optional Eigen (intended backend for the Milestone-3 panel solver). Not needed
# to build the current scaffold, which ships a self-contained LU in linalg.h.
find_package(Eigen3 3.3 QUIET NO_MODULE)

file(GLOB CORE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
add_library(aerocore STATIC ${CORE_SRC})
target_include_directories(aerocore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(OpenMP_CXX_FOUND)
  target_link_libraries(aerocore PUBLIC OpenMP::OpenMP_CXX)
  message(STATUS "OpenMP: enabled (GA evaluation parallelized)")
else()
  message(STATUS "OpenMP: not found (GA runs single-threaded)")
endif()
if(TARGET Eigen3::Eigen)
  target_link_libraries(aerocore PUBLIC Eigen3::Eigen)
  target_compile_definitions(aerocore PUBLIC AERO_HAVE_EIGEN=1)
  message(STATUS "Eigen3: found")
endif()

add_executable(aeroanalyzer app/main.cpp)
target_link_libraries(aeroanalyzer PRIVATE aerocore)

# Offline viscous-surrogate generator (drives XFOIL); not needed at runtime.
add_executable(build_surrogate tools/build_surrogate.cpp)
target_link_libraries(build_surrogate PRIVATE aerocore)

enable_testing()
file(GLOB TEST_SRC ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
add_executable(unit_tests ${TEST_SRC})
target_include_directories(unit_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
target_link_libraries(unit_tests PRIVATE aerocore)
add_test(NAME unit_tests COMMAND unit_tests)
