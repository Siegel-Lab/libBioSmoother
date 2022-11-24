#include "cm/annotation_index.h"
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
                                            std::vector<std::array<size_t, 2>>& vNumInteractionsTotal,
                                            size_t uiNumBinsInRowTotal,
                                            double fPAccept ) override
    {
        pybind11::gil_scoped_acquire acquire;
        PYBIND11_OVERRIDE( ret_t, /* Return type */
                           PartialQuarry, /* Parent class */
                           normalizeBinominalTestTrampoline, /* Name of function in C++ (must match Python name) */
                           vFlatValues, /* Argument(s) */
                           vNumInteractionsTotal, /* Argument(s) */
                           uiNumBinsInRowTotal, /* Argument(s) */
                           fPAccept /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    std::vector<double> normalizeCoolerTrampoline( std::vector<size_t>& vFlatValues, size_t uiAxisSize ) override
    {
        pybind11::gil_scoped_acquire acquire;
        PYBIND11_OVERRIDE( std::vector<double>, /* Return type */
                           PartialQuarry, /* Parent class */
                           normalizeCoolerTrampoline, /* Name of function in C++ (must match Python name) */
                           vFlatValues, /* Argument(s) */
                           uiAxisSize /* Argument(s) */
        );
    }

    std::vector<std::string> colorPalette( std::string sPaletteName, std::string sColorLow,
                                           std::string sColorHigh ) override
    {
        pybind11::gil_scoped_acquire acquire;
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
    using PartialQuarry::normalizeCoolerTrampoline;
};


void testFunc( )
{
    pybind11::gil_scoped_release release;
    std::cout << "Sleeping..." << std::endl;
    std::this_thread::sleep_for( std::chrono::seconds( 3 ) );
    std::cout << "I'm awake! I'm awake!" << std::endl;
}

} // namespace cm


template <bool CACHE> void exportSpsInterface( pybind11::module& m )
{

    pybind11::class_<cm::SpsInterface<CACHE>>( m, CACHE ? "CachedSpsInterface" : "DiskSpsInterface" ) //
        .def( pybind11::init<std::string, bool>( ) ) //
        .def( "loaded", &cm::SpsInterface<CACHE>::loaded )
        .def( "clear_points_and_desc", &cm::SpsInterface<CACHE>::clearPointsAndDesc )
        .def( "insert", &cm::SpsInterface<CACHE>::insert, //
              pybind11::arg( "d" ), pybind11::arg( "o" ), pybind11::arg( "start" ), pybind11::arg( "end" ) )
        .def( "generate", &cm::SpsInterface<CACHE>::generate, //
              pybind11::arg( "d" ), pybind11::arg( "o" ), pybind11::arg( "fac" ) = -1,
              pybind11::arg( "verbosity" ) = 1 )
        .def_readwrite( "anno", &cm::SpsInterface<CACHE>::vAnno );
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

    m.def( "test", &cm::testFunc );


    pybind11::class_<cm::PartialQuarry, cm::PyPartialQuarry>( m, "PartialQuarry" ) //
        .def( pybind11::init<std::string>( ) ) //
        .def( pybind11::init<>( ) ) //
        .def( "set_session", &cm::PartialQuarry::setSession ) //
        .def( "get_session", &cm::PartialQuarry::getSession ) //
        .def( "get_value", &cm::PartialQuarry::getValue<pybind11::object> ) //

        // @note order is relevant here -> the functions are tried in order
        //      if e.g. json would be first everything would ose the json overload
        //      also int, double, string and bool can sometimes convert into one another
        .def( "set_value", &cm::PartialQuarry::setValue<bool> ) //
        .def( "set_value", &cm::PartialQuarry::setValue<int> ) //
        .def( "set_value", &cm::PartialQuarry::setValue<double> ) //
        .def( "set_value", &cm::PartialQuarry::setValue<std::vector<std::string>> ) //
        .def( "set_value", &cm::PartialQuarry::setValue<std::string> ) //
        .def( "set_value", &cm::PartialQuarry::setValue<json> ) //

        .def( "has_undo", &cm::PartialQuarry::hasUndo ) //
        .def( "undo", &cm::PartialQuarry::undo ) //
        .def( "has_redo", &cm::PartialQuarry::hasRedo ) //
        .def( "redo", &cm::PartialQuarry::redo ) //
        .def( "get_palette", &cm::PartialQuarry::getPalette ) //
        .def( "get_bin_coords", &cm::PartialQuarry::getBinCoords ) //
        .def( "get_axis_coords", &cm::PartialQuarry::getAxisCoords ) //
        .def( "get_annotation", &cm::PartialQuarry::getAnnotation ) //
        .def( "get_drawing_area", &cm::PartialQuarry::getDrawingArea ) //
        .def( "get_displayed_annos", &cm::PartialQuarry::getDisplayedAnnos ) //
        .def( "get_heatmap", &cm::PartialQuarry::getHeatmap ) //
        .def( "get_heatmap_export", &cm::PartialQuarry::getHeatmapExport ) //
        .def( "get_dot", &cm::PartialQuarry::getDOT ) //
        .def( "get_background_color", &cm::PartialQuarry::getBackgroundColor ) //
        .def( "get_ticks", &cm::PartialQuarry::getTicks ) //
        .def( "get_tick_list", &cm::PartialQuarry::getTickList ) //
        .def( "get_canvas_size", &cm::PartialQuarry::getCanvasSize ) //
        .def( "get_tracks", &cm::PartialQuarry::getTracks ) //
        .def( "get_ranked_slices", &cm::PartialQuarry::getRankedSlices ) //
        .def( "get_min_max_tracks", &cm::PartialQuarry::getMinMaxTracks ) //
        .def( "get_bin_size", &cm::PartialQuarry::getBinSize ) //
        .def( "get_annotation_list", &cm::PartialQuarry::getAnnotationList ) //
        .def( "get_track_export", &cm::PartialQuarry::getTrackExport ) //
        .def( "get_track_export_names", &cm::PartialQuarry::getTrackExportNames ) //
        .def( "get_decay", &cm::PartialQuarry::getDecayCDS ) //
        .def( "interpret_name", &cm::PartialQuarry::interpretName ) //

        .def( "print_sizes", &cm::PartialQuarry::printSizes ) //
        .def( "cancel", &cm::PartialQuarry::cancel ) //
        .def( "get_error", &cm::PartialQuarry::getError ) //
        .def( "update_cds", &cm::PartialQuarry::updateCDS ) //
        .def( "save_session", &cm::PartialQuarry::saveSession ) //

        .def( "normalizeBinominalTestTrampoline", &cm::ContectMappingPublicist::normalizeBinominalTestTrampoline ) //
        .def( "normalizeCoolerTrampoline", &cm::ContectMappingPublicist::normalizeCoolerTrampoline ) //
        .def( "colorPalette", &cm::ContectMappingPublicist::colorPalette ) //
        ;

    pybind11::class_<cm::AnnotationDescIndex<DiskVecGenerator>>( m, ( "AnnoIndex" ) ) //
        .def( "add_intervals", &cm::AnnotationDescIndex<DiskVecGenerator>::addIntervals, pybind11::arg( "intervals" ),
              pybind11::arg( "verbosity" ) = 0 ) //
        ;

    exportSpsInterface<true>( m );
    exportSpsInterface<false>( m );
}