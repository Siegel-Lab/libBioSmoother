#include "cm/build_time.h"
#include "cm/include.h"
#include "cm/version.h"
#include "sps/abstract_index.h"
#include "sps/build_time.h"
#include "sps/version.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE( libContactMapping, m )
{
    // prevent creation of stxxl log files
    if( getenv( (char*)"STXXLLOGFILE" ) == nullptr )
        putenv( (char*)"STXXLLOGFILE=/dev/null" );
    if( getenv( (char*)"STXXLERRLOGFILE" ) == nullptr )
        putenv( (char*)"STXXLERRLOGFILE=/dev/null" );

    m.attr( "CM_VERSION" ) = CM_VERSION;
    m.attr( "CM_BUILD_TIME" ) = CM_BUILD_TIME;
    m.attr( "SPS_VERSION" ) = SPS_VERSION;
    m.attr( "SPS_BUILD_TIME" ) = SPS_BUILD_TIME;

    pybind11::class_<cm::ContactMapping>( m, "Computation" ).def( "all", &cm::ContactMapping::computeAll );
}