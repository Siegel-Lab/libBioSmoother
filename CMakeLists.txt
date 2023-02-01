cmake_minimum_required (VERSION 3.17)

set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")

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

find_package (Python QUIET COMPONENTS Interpreter Development)

find_package(libsps REQUIRED PATHS ${Python_SITEARCH})
find_package(pybind11 REQUIRED PATHS ${Python_SITEARCH})

pybind11_add_module( libsmoothercpp src/main.cpp )
target_include_directories( libsmoothercpp PUBLIC inc )
target_include_directories( libsmoothercpp PUBLIC contrib/inc )
target_include_directories( libsmoothercpp PUBLIC ${CMAKE_BINARY_DIR}/generated/inc/ )
target_include_directories( libsmoothercpp SYSTEM PRIVATE ${STXXL_INCLUDE_DIRS} )
target_link_libraries( libsmoothercpp PRIVATE libsps )
target_compile_definitions( libsmoothercpp PRIVATE WITH_PYTHON )
add_dependencies( libsmoothercpp cm_config_version cm_config_build_time )



