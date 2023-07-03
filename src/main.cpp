#include "cm/annotation_index.h"
#include "cm/build_time.h"
#include "cm/include.h"
#include "cm/version.h"
#include "sps/abstract_index.h"
#include "sps/build_time.h"
#include "sps/version.h"

#include <pybind11/functional.h>
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
    ret_t normalizeBinominalTestTrampoline( const std::vector<std::array<size_t, 2>>& vFlatValues,
                                            const std::vector<std::array<size_t, 2>>& vNumInteractionsTotal,
                                            const std::vector<std::array<size_t, 2>>& vNumBinsInteractingWith,
                                            size_t uiSamples,
                                            size_t uiMaxNumInteractingWith,
                                            double fPAccept,
                                            bool bIsCol,
                                            size_t uiGridHeight ) override
    {
        pybind11::gil_scoped_acquire acquire;
        PYBIND11_OVERRIDE( ret_t, /* Return type */
                           PartialQuarry, /* Parent class */
                           normalizeBinominalTestTrampoline, /* Name of function in C++ (must match Python name) */
                           vFlatValues, /* Argument(s) */
                           vNumInteractionsTotal, /* Argument(s) */
                           vNumBinsInteractingWith, /* Argument(s) */
                           uiSamples, /* Argument(s) */
                           uiMaxNumInteractingWith, /* Argument(s) */
                           fPAccept, /* Argument(s) */
                           bIsCol, /* Argument(s) */
                           uiGridHeight );
    }

    /* Trampoline (need one for each virtual function) */
    std::vector<double> normalizeCoolerTrampoline( const std::vector<size_t>& vFlatValues, size_t uiAxisSize ) override
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

pybind11::dict test_py_dict( size_t uiA, size_t uiB )
{
    using namespace pybind11::literals;
    pybind11::gil_scoped_acquire acquire;

    pybind11::dict xRet;
    for( size_t uiI = 0; uiI < uiA; uiI++ )
    {
        pybind11::list vL;
        for( size_t uiJ = 0; uiJ < uiB; uiJ++ )
            vL.append( uiJ );
        xRet[ std::to_string( uiI ).c_str( ) ] = vL;
    }
    return xRet;
}

std::map<std::string, std::vector<size_t>> test_cpp_dict( size_t uiA, size_t uiB )
{
    std::map<std::string, std::vector<size_t>> xRet;

    for( size_t uiI = 0; uiI < uiA; uiI++ )
    {
        std::vector<size_t> vL;
        vL.reserve( uiB );
        for( size_t uiJ = 0; uiJ < uiB; uiJ++ )
            vL.push_back( uiJ );
        xRet[ std::to_string( uiI ) ] = vL;
    }
    return xRet;
}

} // namespace cm


template <bool CACHE> void exportSpsInterface( pybind11::module& m )
{
    pybind11::class_<cm::SpsInterface<CACHE>, std::shared_ptr<cm::SpsInterface<CACHE>>>( m, CACHE ? "CachedIndex"
                                                                                                  : "Index" ) //
        .def( pybind11::init<std::string, bool>( ) ) //
        .def( pybind11::init<std::string>( ) ) //
        .def( "loaded", &cm::SpsInterface<CACHE>::loaded )
        .def( "clear_points_and_desc", &cm::SpsInterface<CACHE>::clearPointsAndDesc )
        .def( "insert",
              static_cast<void ( cm::SpsInterface<CACHE>::* )( std::vector<uint64_t>, std::vector<uint64_t>, int )>(
                  &cm::SpsInterface<CACHE>::insert ), //
              pybind11::arg( "start" ), pybind11::arg( "end" ), pybind11::arg( "val" ) )
        .def( "insert_bias",
              static_cast<void ( cm::SpsInterface<CACHE>::* )( std::vector<uint64_t>, double )>(
                  &cm::SpsInterface<CACHE>::insertBias ), //
              pybind11::arg( "start" ), pybind11::arg( "val" ) )
        .def( "generate", &cm::SpsInterface<CACHE>::generate, //
              pybind11::arg( "fac" ) = -1.0, pybind11::arg( "verbosity" ) = 1 )
        .def( "generate_bias", &cm::SpsInterface<CACHE>::generateBias, //
              pybind11::arg( "fac" ) = -1.0, pybind11::arg( "verbosity" ) = 1 )
        .def_readwrite( "anno", &cm::SpsInterface<CACHE>::vAnno )
        .def( "get_value", &cm::HasSession::getSessionValue<pybind11::object> ) //

        // @note order is relevant here -> the functions are tried in order
        //      if e.g. json would be first everything would ose the json overload
        //      also int, double, string and bool can sometimes convert into one another
        .def( "set_value", &cm::HasSession::setSessionValue<bool> ) //
        .def( "set_value", &cm::HasSession::setSessionValue<int> ) //
        .def( "set_value", &cm::HasSession::setSessionValue<double> ) //
        .def( "set_value", &cm::HasSession::setSessionValue<std::vector<std::string>> ) //
        .def( "set_value", &cm::HasSession::setSessionValue<std::string> ) //
        .def( "set_value", &cm::HasSession::setSessionValue<json> ) //
        .def( "save_session", &cm::HasSession::saveSession ) //
        ;
}


PYBIND11_MODULE( libbiosmoothercpp, m )
{
#ifdef WITH_STXXL
    // prevent creation of stxxl log files
    if( getenv( (char*)"STXXLLOGFILE" ) == nullptr )
        putenv( (char*)"STXXLLOGFILE=/dev/null" );
    if( getenv( (char*)"STXXLERRLOGFILE" ) == nullptr )
        putenv( (char*)"STXXLERRLOGFILE=/dev/null" );
#endif

    m.attr( "LIB_BIO_SMOOTHER_CPP_VERSION" ) = LIB_BIO_SMOOTHER_CPP_VERSION;
    m.attr( "LIB_BIO_SMOOTHER_CPP_BUILD_TIME" ) = LIB_BIO_SMOOTHER_CPP_BUILD_TIME;
    m.attr( "SPS_VERSION" ) = SPS_VERSION;
    m.attr( "SPS_BUILD_TIME" ) = SPS_BUILD_TIME;

#ifdef NDEBUG
    m.attr( "DEBUG" ) = false;
#else
    m.attr( "DEBUG" ) = true;
#endif

#ifdef WITH_STXXL
    m.attr( "WITH_STXXL" ) = true;
#else
    m.attr( "WITH_STXXL" ) = false;
#endif

    m.attr( "COMPILER_ID" ) = CXX_COMPILER_ID;

    m.def( "test", &cm::testFunc );

    // observation: both functions are roughly the same speed -> not worth optimizing
    m.def( "test_py_dict", &cm::test_py_dict );
    m.def( "test_cpp_dict", &cm::test_cpp_dict );

    pybind11::class_<cm::AxisCoord>( m, "AxisCoord" ) //
        .def_readwrite( "chr_idx", &cm::AxisCoord::uiChromosome ) //
        .def_readwrite( "idx_pos", &cm::AxisCoord::uiIndexPos ) //
        .def_readwrite( "idx_size", &cm::AxisCoord::uiIndexSize ) //
        ;

    pybind11::class_<cm::PartialQuarry, cm::PyPartialQuarry>( m, "PartialQuarry" ) //
        .def( pybind11::init<std::string>( ) ) //
        .def( pybind11::init<std::shared_ptr<cm::SpsInterface<false>>>( ) ) //
        .def( pybind11::init<>( ) ) //
        .def_readwrite( "allow_ctrl_c_cancel", &cm::PartialQuarry::bAllowCtrlCCancel ) //
        .def_readwrite( "verbosity", &cm::PartialQuarry::uiVerbosity ) //
        .def_readwrite( "index", &cm::PartialQuarry::pIndices ) //
        .def( "set_session", &cm::PartialQuarry::setSession ) //
        .def( "get_session", &cm::PartialQuarry::getSession ) //
        .def( "has_value", &cm::PartialQuarry::hasValue ) //
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

        .def( "copy_from", &cm::PartialQuarry::copyFrom ) //
        .def( "has_undo", &cm::PartialQuarry::hasUndo ) //
        .def( "undo", &cm::PartialQuarry::undo ) //
        .def( "has_redo", &cm::PartialQuarry::hasRedo ) //
        .def( "redo", &cm::PartialQuarry::redo ) //
        .def( "get_palette", &cm::PartialQuarry::getPalette ) //
        .def( "get_bin_coords", &cm::PartialQuarry::getBinCoords ) //
        .def( "get_axis_coords", &cm::PartialQuarry::getAxisCoords ) //
        .def( "get_axis_size", &cm::PartialQuarry::getAxisSize ) //
        .def( "get_annotation", &cm::PartialQuarry::getAnnotation ) //
        .def( "get_drawing_area", &cm::PartialQuarry::getDrawingArea ) //
        .def( "get_displayed_annos", &cm::PartialQuarry::getDisplayedAnnos ) //
        .def( "get_heatmap", &cm::PartialQuarry::getHeatmap ) //
        .def( "get_heatmap_export", &cm::PartialQuarry::getHeatmapExport ) //
        .def( "get_scaled", &cm::PartialQuarry::getScaled ) //
        .def( "get_combined", &cm::PartialQuarry::getCombined ) //
        .def( "get_dot", &cm::PartialQuarry::getDOT ) //
        .def( "get_background_color", &cm::PartialQuarry::getBackgroundColor ) //
        .def( "get_ticks", &cm::PartialQuarry::getTicks ) //
        .def( "get_contig_ticks", &cm::PartialQuarry::getContigTicks ) //
        .def( "get_tick_list", &cm::PartialQuarry::getTickList ) //
        .def( "get_contig_start_list", &cm::PartialQuarry::getContigStartList ) //
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
        .def( "get_palette_ticks", &cm::PartialQuarry::getPaletteTicks ) //
        .def( "get_longest_common_suffix", &cm::PartialQuarry::getLongestCommonSuffix ) //

        .def( "get_file_prefix", &cm::PartialQuarry::getFilePrefix ) //
        .def( "print_num_reads_and_overlays", &cm::PartialQuarry::printNumReadsNumOverlays ) //
        .def( "print_sizes", &cm::PartialQuarry::printSizes ) //
        .def( "cancel", &cm::PartialQuarry::cancel ) //
        .def( "get_error", &cm::PartialQuarry::getError ) //
        .def( "update_cds", &cm::PartialQuarry::updateCDS ) //
        .def( "clear_cache", &cm::PartialQuarry::clearCache ) //
        .def( "update_all", &cm::PartialQuarry::updateAll ) //
        .def( "save_session", &cm::PartialQuarry::saveSession ) //
        .def( "get_runtimes", &cm::PartialQuarry::getRuntimes ) //

        .def( "normalizeBinominalTestTrampoline", &cm::ContectMappingPublicist::normalizeBinominalTestTrampoline ) //
        .def( "normalizeCoolerTrampoline", &cm::ContectMappingPublicist::normalizeCoolerTrampoline ) //
        .def( "colorPalette", &cm::ContectMappingPublicist::colorPalette ) //
        ;

    pybind11::class_<cm::AnnotationDescIndex<DiskVecGenerator>>( m, ( "AnnoIndex" ) ) //
        .def( "add_intervals", &cm::AnnotationDescIndex<DiskVecGenerator>::addIntervals, pybind11::arg( "intervals" ),
              pybind11::arg( "divided" ),
              pybind11::arg( "verbosity" ) = 0 ) //
        .def( "get_categories", &cm::AnnotationDescIndex<DiskVecGenerator>::getCategories, pybind11::arg( "pos" ),
              pybind11::arg( "dividend" ), pybind11::arg( "relevant" ), pybind11::arg( "interval_coords" ) = false,
              pybind11::arg( "interval_count" ) = false ) //
        ;

    exportSpsInterface<true>( m );
    exportSpsInterface<false>( m );
}