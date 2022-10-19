#include "cm/build_time.h"
#include "cm/include.h"
#include "cm/version.h"
#include "sps/abstract_index.h"
#include "sps/build_time.h"
#include "sps/version.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace cm
{

class PyPartialQuarry : public PartialQuarry
{
  public:
    /* Inherit the constructors */
    using PartialQuarry::PartialQuarry;

    using ret_t = std::vector<std::array<double, 2>>;
    /* Trampoline (need one for each virtual function) */
    ret_t normalizeBinominalTestTrampoline( std::vector<std::array<size_t, 2>>& vFlatValues,
                                            size_t uiNumInteractionsTotal, size_t uiNumBinsInRowTotal,
                                            double fPAccept ) override
    {
        PYBIND11_OVERRIDE( ret_t, /* Return type */
                           PartialQuarry, /* Parent class */
                           normalizeBinominalTestTrampoline, /* Name of function in C++ (must match Python name) */
                           vFlatValues, /* Argument(s) */
                           uiNumInteractionsTotal, /* Argument(s) */
                           uiNumBinsInRowTotal, /* Argument(s) */
                           fPAccept /* Argument(s) */
        );
    }

    std::vector<std::string> colorPalette( std::string sPaletteName, std::string sColorLow,
                                           std::string sColorHigh ) override
    {
        PYBIND11_OVERRIDE( std::vector<std::string>, /* Return type */
                           PartialQuarry, /* Parent class */
                           colorPalette, /* Name of function in C++ (must match Python name) */
                           sPaletteName, /* Argument(s) */
                           sColorLow, /* Argument(s) */
                           sColorHigh /* Argument(s) */
        );
    }
};

class ContectMappingPublicist : public PartialQuarry
{
  public:
    using PartialQuarry::colorPalette;
    using PartialQuarry::normalizeBinominalTestTrampoline;
};

} // namespace cm

PYBIND11_MODULE( libPartialQuarry, m )
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

    pybind11::class_<cm::PartialQuarry, cm::PyPartialQuarry>( m, "PartialQuarry" ) //
        .def( pybind11::init<std::string>( ) ) //
        .def( "get_colors", &cm::PartialQuarry::getColors ) //
        .def( "normalizeBinominalTestTrampoline", &cm::ContectMappingPublicist::normalizeBinominalTestTrampoline ) //
        .def( "colorPalette", &cm::ContectMappingPublicist::colorPalette ) //
        ;

    pybind11::class_<cm::SpsInterface<true>>( m, "CachedSpsInterface" ) //
        .def( pybind11::init<std::string>( ) ) //
        ;

    pybind11::class_<cm::SpsInterface<false>>( m, "DiskSpsInterface" ) //
        .def( pybind11::init<std::string>( ) ) //
        ;
}