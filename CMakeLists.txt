cmake_minimum_required (VERSION 3.9)

# only adjust major and minor version - patch is taken from git commit number
project("libSmootherCpp" VERSION 0.1.0)

# disallow in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)


# enable warnings (always good)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall")

# Add -O0 to remove optimizations when using gcc
IF(CMAKE_COMPILER_IS_GNUCC)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
ENDIF(CMAKE_COMPILER_IS_GNUCC)

# Enable Link time optimization
include(CheckIPOSupported)
check_ipo_supported(RESULT supported OUTPUT error)
if( supported )
    message(STATUS "IPO / LTO enabled")
    set(INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(STATUS "IPO / LTO not supported: <${error}>")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")

# configure the version.h.in file
add_custom_target(
    cm_config_version ALL
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/src/conf_version.h.in.sh
        "${CMAKE_CURRENT_SOURCE_DIR}/inc/cm/version.h.in"
        "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
        "${CMAKE_BINARY_DIR}/generated/inc/cm/version.h"
    VERBATIM)

add_custom_target(
    cm_config_build_time ALL
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/src/conf_build_time.h.in.sh
        "${CMAKE_CURRENT_SOURCE_DIR}/inc/cm/build_time.h.in"
        "${CMAKE_BINARY_DIR}/generated/inc/cm/build_time.h"
    VERBATIM)

add_subdirectory(../libSps ${CMAKE_BINARY_DIR}/libSps EXCLUDE_FROM_ALL)
set(JSON_BuildTests OFF CACHE INTERNAL "")
find_package(pybind11 REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(pybind11_json REQUIRED)

pybind11_add_module( libSmootherCpp src/main.cpp )
target_include_directories( libSmootherCpp PUBLIC inc )
target_link_libraries( libSmootherCpp PRIVATE Sps nlohmann_json::nlohmann_json )
target_compile_definitions( libSmootherCpp PRIVATE WITH_PYTHON )
add_dependencies( libSmootherCpp cm_config_version cm_config_build_time )



