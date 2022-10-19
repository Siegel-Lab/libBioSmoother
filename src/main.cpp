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

class PyContactMapping : public ContactMapping
{
  public:
    /* Inherit the constructors */
    using ContactMapping::ContactMapping;

    using ret_t = std::vector<std::array<double, 2>>;
    /* Trampoline (need one for each virtual function) */
    ret_t normalizeBinominalTestTrampoline( std::vector<std::array<size_t, 2>>& vFlatValues,
                                            size_t uiNumInteractionsTotal, size_t uiNumBinsInRowTotal,
                                            double fPAccept ) override
    {
        PYBIND11_OVERRIDE( ret_t, /* Return type */
                           ContactMapping, /* Parent class */
                           normalizeBinominalTestTrampoline, /* Name of function in C++ (must match Python name) */
                           vFlatValues, /* Argument(s) */
                           uiNumInteractionsTotal, /* Argument(s) */
                           uiNumBinsInRowTotal, /* Argument(s) */
                           fPAccept /* Argument(s) */
        );
    }

    std::vector<std::string> colorPalette( std::string sPaletteName ) override
    {
        PYBIND11_OVERRIDE( std::vector<std::string>, /* Return type */
                           ContactMapping, /* Parent class */
                           colorPalette, /* Name of function in C++ (must match Python name) */
                           sPaletteName /* Argument(s) */
        );
    }
};

class ContectMappingPublicist : public ContactMapping
{
  public:
    using ContactMapping::normalizeBinominalTestTrampoline;
    using ContactMapping::colorPalette;
};

} // namespace cm

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

    pybind11::class_<cm::ContactMapping, cm::PyContactMapping>( m, "ContactMapping" ) //
        .def( pybind11::init<std::string>( ) ) //
        .def( "get_colors", &cm::ContactMapping::getColors ) //
        .def( "normalizeBinominalTestTrampoline", &cm::ContectMappingPublicist::normalizeBinominalTestTrampoline ) //
        .def( "colorPalette", &cm::ContectMappingPublicist::colorPalette ) //
        ;
}