cmake_minimum_required (VERSION 3.17)

set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")

# only adjust major and minor version - patch is taken from git commit number
project("libBioSmootherCpp" VERSION ${LIB_BIO_SMOOTHER_VERSION})

# disallow in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    # enable warnings (always good)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall")
    # throw error for undefined symbols
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-undefined")

    # Add -O0 to remove optimizations when using gcc
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
ENDIF(CMAKE_CXX_COMPILER_ID MATCHES "GNU")

# Enable Link time optimization
include(CheckIPOSupported)
check_ipo_supported(RESULT supported OUTPUT error)
if( supported )
    message(STATUS "IPO / LTO enabled")
    set(INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(STATUS "IPO / LTO not supported: <${error}>")
endif()


find_package (Python REQUIRED COMPONENTS Interpreter Development.Module)

# configure the version.h.in file
add_custom_target(
    cm_config_version ALL
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/src/conf_version.h.in.py
        "${CMAKE_CURRENT_SOURCE_DIR}/inc/cm/version.h.in"
        "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}"
        "${CMAKE_BINARY_DIR}/generated/inc/cm/version.h"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    VERBATIM)

add_custom_target(
    cm_config_build_time ALL
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/src/conf_build_time.h.in.py
        "${CMAKE_CURRENT_SOURCE_DIR}/inc/cm/build_time.h.in"
        "${CMAKE_BINARY_DIR}/generated/inc/cm/build_time.h"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    VERBATIM)


find_package(libsps REQUIRED PATHS ${Python_SITEARCH})
find_package(pybind11 REQUIRED PATHS ${Python_SITEARCH})

if(WITH_STXXL)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")
endif(WITH_STXXL)

pybind11_add_module( libbiosmoothercpp src/main.cpp )
target_include_directories( libbiosmoothercpp PUBLIC inc )
target_include_directories( libbiosmoothercpp PUBLIC contrib/inc )
target_include_directories( libbiosmoothercpp PUBLIC ${CMAKE_BINARY_DIR}/generated/inc/ )
target_link_libraries( libbiosmoothercpp PRIVATE libsps )
if(WITH_STXXL)
    target_include_directories( libbiosmoothercpp SYSTEM PRIVATE ${STXXL_INCLUDE_DIRS} )
    target_link_libraries( libbiosmoothercpp PRIVATE stxxl )
endif(WITH_STXXL)
target_compile_definitions( libbiosmoothercpp PRIVATE WITH_PYTHON )
target_compile_definitions( libbiosmoothercpp PRIVATE CXX_COMPILER_ID="${CMAKE_CXX_COMPILER_ID}" )
add_dependencies( libbiosmoothercpp cm_config_version cm_config_build_time )

file(GLOB format_files
     "inc/cm/*.h"
     "src/*.cpp"
)
add_custom_target(format 
                  COMMAND "clang-format" "-i" "--style=file" ${format_files}
                  COMMAND "black" "libbiosmoother/"
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} VERBATIM)


# BUILD WHEEL
# https://wheel.readthedocs.io/en/stable/user_guide.html#building-wheels

# python -m pip install build
# python -m build --wheel
