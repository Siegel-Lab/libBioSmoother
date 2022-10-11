#include "sps/abstract_index.h"
#include "cm/version.h"
#include "cm/build_time.h"

#if WITH_PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#endif


PYBIND11_MODULE( libContactMapping, m )
{
    // prevent creation of stxxl log files
    if( getenv( (char*)"STXXLLOGFILE" ) == nullptr )
        putenv( (char*)"STXXLLOGFILE=/dev/null" );
    if( getenv( (char*)"STXXLERRLOGFILE" ) == nullptr )
        putenv( (char*)"STXXLERRLOGFILE=/dev/null" );

    m.attr( "VERSION" ) = VERSION;
    m.attr( "BUILD_TIME" ) = BUILD_TIME;

}