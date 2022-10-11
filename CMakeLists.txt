cmake_minimum_required (VERSION 3.9)

# only adjust major and minor version - patch is taken from git commit number
project("libContactMapping" VERSION 0.1.0)

# set building python library as default
option(WITH_PYTHON "Build python library." ON)

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
    COMMAND ${CMAKE_SOURCE_DIR}/src/conf_version.h.in.sh
        "${CMAKE_SOURCE_DIR}/inc/cm/version.h.in"
        "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
        "${CMAKE_BINARY_DIR}/generated/inc/cm/version.h"
    VERBATIM)

add_custom_target(
    cm_config_build_time ALL
    COMMAND ${CMAKE_SOURCE_DIR}/src/conf_build_time.h.in.sh
        "${CMAKE_SOURCE_DIR}/inc/cm/build_time.h.in"
        "${CMAKE_BINARY_DIR}/generated/inc/cm/build_time.h"
    VERBATIM)

add_subdirectory(../libSps ${CMAKE_BINARY_DIR}/libSps EXCLUDE_FROM_ALL)

# create python module if python has been found
if(WITH_PYTHON)
    pybind11_add_module( ContactMapping src/main.cpp )
    target_include_directories( ContactMapping PUBLIC inc )
    target_link_libraries( ContactMapping PUBLIC libSps )
    add_dependencies(ContactMapping cm_config_version cm_config_build_time)
endif()



