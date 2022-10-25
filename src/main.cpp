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

template <bool CACHE> void exportSpsInterface( pybind11::module& m )
{
    pybind11::class_<cm::SpsInterface<CACHE>>( m, CACHE ? "CachedSpsInterface" : "DiskSpsInterface" ) //
        .def( pybind11::init<std::string>( ) ) //
        .def( "loaded", &cm::SpsInterface<CACHE>::loaded )
        .def( "clear_points_and_desc", &cm::SpsInterface<CACHE>::clearPointsAndDesc )
        .def( "insert", &cm::SpsInterface<CACHE>::insert, //
              pybind11::arg( "d" ), pybind11::arg( "o" ), pybind11::arg( "start" ), pybind11::arg( "end" ),
              pybind11::arg( "desc" ) = "" )
        .def( "generate", &cm::SpsInterface<CACHE>::generate, //
              pybind11::arg( "d" ), pybind11::arg( "o" ), pybind11::arg( "fac" ) = -1,
              pybind11::arg( "verbosity" ) = 1 );
}

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

    pybind11::class_<cm::BinCoord>( m, "BinCoord" ) //
        .def( pybind11::init<>( ) ) //
        .def_readwrite( "chr_x", &cm::BinCoord::sChromosomeX ) //
        .def_readwrite( "chr_y", &cm::BinCoord::sChromosomeY ) //
        .def_readwrite( "screen_x", &cm::BinCoord::uiScreenX ) //
        .def_readwrite( "screen_y", &cm::BinCoord::uiScreenY ) //
        .def_readwrite( "index_x", &cm::BinCoord::uiIndexX ) //
        .def_readwrite( "index_y", &cm::BinCoord::uiIndexY ) //
        .def_readwrite( "w", &cm::BinCoord::uiW ) //
        .def_readwrite( "h", &cm::BinCoord::uiH ) //
        ;

    pybind11::class_<cm::AxisCoord>( m, "AxisCoord" ) //
        .def( pybind11::init<>( ) ) //
        .def_readwrite( "chr", &cm::AxisCoord::sChromosome ) //
        .def_readwrite( "screen_pos", &cm::AxisCoord::uiScreenPos ) //
        .def_readwrite( "index_pos", &cm::AxisCoord::uiIndexPos ) //
        .def_readwrite( "size", &cm::AxisCoord::uiSize ) //
        ;

    pybind11::class_<cm::PartialQuarry, cm::PyPartialQuarry>( m, "PartialQuarry" ) //
        .def( pybind11::init<std::string>( ) ) //
        .def( pybind11::init<>( ) ) //
        .def( "set_session", &cm::PartialQuarry::setSession ) //
        .def( "get_value", &cm::PartialQuarry::getValue<pybind11::object> ) //
        .def( "set_value", &cm::PartialQuarry::setValue<int> ) //
        .def( "set_value", &cm::PartialQuarry::setValue<std::string> ) //
        .def( "set_value", &cm::PartialQuarry::setValue<bool> ) //
        .def( "set_value", &cm::PartialQuarry::setValue<double> ) //
        .def( "set_value", &cm::PartialQuarry::setValue<json> ) //
        .def( "set_value", &cm::PartialQuarry::setValue<std::vector<std::string>> ) //
        .def( "has_undo", &cm::PartialQuarry::hasUndo ) //
        .def( "undo", &cm::PartialQuarry::undo ) //
        .def( "has_redo", &cm::PartialQuarry::hasRedo ) //
        .def( "redo", &cm::PartialQuarry::redo ) //
        .def( "get_colors", &cm::PartialQuarry::getColors ) //
        .def( "get_bin_coords", &cm::PartialQuarry::getBinCoords ) //
        .def( "get_axis_coords", &cm::PartialQuarry::getAxisCoords ) //
        .def( "get_annotation", &cm::PartialQuarry::getAnnotation ) //
        .def( "get_drawing_area", &cm::PartialQuarry::getDrawingArea ) //
        .def( "get_dot", &cm::PartialQuarry::getDOT ) //
        .def( "normalizeBinominalTestTrampoline", &cm::ContectMappingPublicist::normalizeBinominalTestTrampoline ) //
        .def( "colorPalette", &cm::ContectMappingPublicist::colorPalette ) //
        ;

    exportSpsInterface<true>( m );
    exportSpsInterface<false>( m );
}