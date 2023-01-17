#include "cm/partial_quarry.h"
#include <cmath>

#pragma once

namespace cm
{

bool PartialQuarry::setActiveCoverage( )
{
    size_t uiSize = getValue<json>( { "coverage", "list" } ).size( );
    for( size_t uiI = 0; uiI < 2; uiI++ )
    {
        vActiveCoverage[ uiI ].clear( );
        vActiveCoverage[ uiI ].reserve( uiSize );
    }

    std::array<std::set<std::string>, 2> vGroups;

    for( auto& rV : getValue<json>( { "coverage", "in_column" } ) )
    {
        CANCEL_RETURN;
        vActiveCoverage[ 0 ].push_back( rV.get<std::string>( ) );
    }
    for( auto& rV : getValue<json>( { "coverage", "in_row" } ) )
    {
        CANCEL_RETURN;
        vActiveCoverage[ 1 ].push_back( rV.get<std::string>( ) );
    }
    END_RETURN;
}

std::tuple<size_t, int64_t, size_t, size_t> PartialQuarry::makeHeapTuple( size_t uiMapQMin,
                                                                          size_t uiMapQMax,
                                                                          bool bCol,
                                                                          bool bSymPart,
                                                                          const AxisCoord& xCoords,
                                                                          int64_t iDataSetId,
                                                                          size_t uiStart,
                                                                          size_t uiEnd )
{
    const size_t uiXMin = bCol != bSymPart ? xCoords.uiIndexPos : uiStart;
    const size_t uiXMax = bCol != bSymPart ? xCoords.uiIndexPos + xCoords.uiIndexSize : uiEnd;
    const size_t uiYMin = bCol != bSymPart ? uiStart : xCoords.uiIndexPos;
    const size_t uiYMax = bCol != bSymPart ? uiEnd : xCoords.uiIndexPos + xCoords.uiIndexSize;
    const size_t uiCount = xIndices.count( iDataSetId, { uiYMin, uiXMin, uiMapQMax, uiFromAnnoFilter },
                                           { uiYMax, uiXMax, uiMapQMin, uiToAnnoFilter }, xIntersect, 0 );

    return std::make_tuple( uiCount, iDataSetId, uiStart, uiEnd );
}

size_t PartialQuarry::getCoverageFromRepl( size_t uiMapQMin,
                                           size_t uiMapQMax,
                                           const AxisCoord& xCoords,
                                           const json& xRep,
                                           size_t uiCoverageGetMaxBinSize,
                                           const std::vector<ChromDesc>& vChroms,
                                           bool bCol,
                                           bool bSymPart )
{
    std::string sChromName = vChroms[ xCoords.uiChromosome ].sName;
    // score, datasetId, start, end
    std::vector<std::tuple<size_t, int64_t, size_t, size_t>> vHeap;
    for( const ChromDesc& rDesc : vChroms )
        vHeap.push_back( makeHeapTuple(
            uiMapQMin, uiMapQMax, bCol, bSymPart, xCoords,
            xRep[ "ids" ][ bCol != bSymPart ? sChromName : rDesc.sName ][ bCol != bSymPart ? rDesc.sName : sChromName ]
                .get<int64_t>( ),
            0, rDesc.uiLength ) );

    std::make_heap( vHeap.begin( ), vHeap.end( ) );

    while( vHeap.size( ) > 0 &&
           std::get<3>( vHeap.back( ) ) - std::get<2>( vHeap.back( ) ) > std::max( uiCoverageGetMaxBinSize, 1ul ) )
    {
        auto xFront = vHeap.back( );
        vHeap.pop_back( );

        size_t uiCenter = ( std::get<2>( xFront ) + std::get<3>( xFront ) ) / 2;
        vHeap.push_back( makeHeapTuple( uiMapQMin, uiMapQMax, bCol, bSymPart, xCoords, std::get<1>( xFront ), uiCenter,
                                        std::get<3>( xFront ) ) );
        std::push_heap( vHeap.begin( ), vHeap.end( ) );

        vHeap.push_back( makeHeapTuple( uiMapQMin, uiMapQMax, bCol, bSymPart, xCoords, std::get<1>( xFront ),
                                        std::get<2>( xFront ), uiCenter ) );
        std::push_heap( vHeap.begin( ), vHeap.end( ) );


        std::pop_heap( vHeap.begin( ), vHeap.end( ) );
    }
    if( vHeap.size( ) > 0 )
        return std::get<0>( vHeap.back( ) );

    return 0;
}

bool PartialQuarry::setCoverageValues( )
{

    size_t uiMinuend = getValue<size_t>( { "settings", "normalization", "min_interactions", "val" } );

    size_t uiMapQMin = getValue<size_t>( { "settings", "filters", "mapping_q", "val_min" } );
    size_t uiMapQMax = getValue<size_t>( { "settings", "filters", "mapping_q", "val_max" } );
    bool bIncomplAlignment = getValue<bool>( { "settings", "filters", "incomplete_alignments" } );

    if( !( uiMapQMin == 0 && bIncomplAlignment ) )
        ++uiMapQMin;
    ++uiMapQMax;

    uiMapQMin = 255 - uiMapQMin;
    uiMapQMax = 255 - uiMapQMax;

    for( size_t uiJ = 0; uiJ < 2; uiJ++ )
    {
        vvCoverageValues[ uiJ ].clear( );
        vvCoverageValues[ uiJ ].reserve( vActiveCoverage[ uiJ ].size( ) );

        for( auto& sRep : vActiveCoverage[ uiJ ] )
        {
            vvCoverageValues[ uiJ ].emplace_back( );
            vvCoverageValues[ uiJ ].back( ).reserve( vAxisCords[ uiJ ].size( ) );
            auto xRep = getValue<json>( { "coverage", "by_name", sRep } );

            for( AxisCoord& xCoords : vAxisCords[ uiJ ] )
            {
                size_t uiVal;
                if( xCoords.uiChromosome != std::numeric_limits<size_t>::max( ) )
                {
                    CANCEL_RETURN;

                    std::string sChromName = vActiveChromosomes[ uiJ ][ xCoords.uiChromosome ].sName;
                    int64_t iDataSetId = xRep[ "ids" ][ sChromName ].get<int64_t>( );
                    if( iDataSetId != -1 )
                        uiVal =
                            xIndices.count( iDataSetId,
                                            { xCoords.uiIndexPos, 0, uiMapQMax, uiFromAnnoFilter },
                                            { xCoords.uiIndexPos + xCoords.uiIndexSize, 1, uiMapQMin, uiToAnnoFilter },
                                            xIntersect,
                                            0 );
                }
                else
                    uiVal = 0;

                if( uiVal > uiMinuend )
                    uiVal -= uiMinuend;
                else
                    uiVal = 0;

                vvCoverageValues[ uiJ ].back( ).push_back( uiVal );
            }
        }
    }
    END_RETURN;
}

double normCoverage( size_t uiNorm, size_t uiVal, size_t uiTotal, size_t uiW, size_t uiH )
{
    uint64_t uiS = (uint64_t)uiW * (uint64_t)uiH;
    switch( uiNorm )
    {
        case 0: // dont
            return (double)uiVal;
        case 1: // rpm
            return uiTotal == 0 ? 0 : 1000000.0 * (double)uiVal / (double)uiTotal;
        case 2: // rpk
            return uiTotal == 0 ? 0 : 1000.0 * (double)uiVal / (double)uiTotal;
        case 3: // rpmb
            return uiS == 0 ? 0 : 1000000.0 * (double)uiVal / (double)uiS;
        case 4: // rpkb
            return uiS == 0 ? 0 : 1000.0 * (double)uiVal / (double)uiS;
        default:
            throw std::logic_error( "invalid value for normalize_by_coverage" );
    }
}


bool PartialQuarry::setTracks( )
{
    using namespace pybind11::literals;
    pybind11::gil_scoped_acquire acquire;

    size_t uiDividend = getValue<size_t>( { "dividend" } );
    for( size_t uiI = 0; uiI < 2; uiI++ )
    {
        vvMinMaxTracks[ uiI ][ 0 ] = std::numeric_limits<double>::max( );
        vvMinMaxTracks[ uiI ][ 1 ] = std::numeric_limits<double>::min( );

        for( size_t uiX = 0; uiX < vAxisCords[ uiI ].size( ); uiX++ )
        {
            CANCEL_RETURN;
            for( size_t uiId = 0; uiId < vvCoverageValues[ uiI ].size( ); uiI++ )
            {
                auto uiVal = (double)vvCoverageValues[ uiI ][ uiId ][ uiX ];
                vvMinMaxTracks[ uiI ][ 0 ] = std::min( vvMinMaxTracks[ uiI ][ 0 ], uiVal );
                vvMinMaxTracks[ uiI ][ 1 ] = std::max( vvMinMaxTracks[ uiI ][ 1 ], uiVal );
            }
        }
        pybind11::list vChrs;
        pybind11::list vScreenPoss;
        pybind11::list vIndexStarts;
        pybind11::list vIndexEnds;
        pybind11::list vValues;
        pybind11::list vColors;
        pybind11::list vNames;

        size_t uiCnt = 0;

        for( size_t uiId = 0; uiId < vvCoverageValues[ uiI ].size( ); uiI++ )
        {
            pybind11::list vScreenPos;
            pybind11::list vIndexStart;
            pybind11::list vIndexEnd;
            pybind11::list vValue;
            std::string sChr = "";


            for( size_t uiX = 0; uiX < vAxisCords[ uiI ].size( ); uiX++ )
            {
                CANCEL_RETURN;
                auto& xCoord = vAxisCords[ uiI ][ uiX ];
                std::string sChromName = vActiveChromosomes[ uiI ][ xCoord.uiChromosome ].sName;
                if( sChr != "" && sChr != sChromName )
                {
                    vChrs.append( substringChr( sChr ) );

                    vScreenPoss.append( vScreenPos );
                    vScreenPos = pybind11::list( );

                    vIndexStarts.append( vIndexStart );
                    vIndexStart = pybind11::list( );

                    vIndexEnds.append( vIndexEnd );
                    vIndexEnd = pybind11::list( );

                    vValues.append( vValue );
                    vValue = pybind11::list( );

                    vColors.append( vColorPaletteAnnotation[ uiCnt % vColorPaletteAnnotation.size( ) ] );

                    vNames.append( vActiveCoverage[ uiI ][ uiId ] );
                }

                if( uiX == 0 )
                {
                    // zero position at start
                    vIndexStart.append( readableBp( xCoord.uiIndexPos * uiDividend ) );
                    vIndexEnd.append( readableBp( xCoord.uiIndexPos * uiDividend ) );
                    vScreenPos.append( xCoord.uiScreenPos );
                    vValue.append( vvMinMaxTracks[ uiI ][ 0 ] );
                }

                sChr = sChromName;
                auto uiVal = vvCoverageValues[ uiI ][ uiId ][ uiX ];

                // front corner
                vIndexStart.append( readableBp( xCoord.uiIndexPos * uiDividend ) );
                vIndexEnd.append( readableBp( ( xCoord.uiIndexPos + xCoord.uiIndexSize ) * uiDividend ) );
                vScreenPos.append( xCoord.uiScreenPos );
                vValue.append( uiVal );

                // rear corner
                vIndexStart.append( readableBp( xCoord.uiIndexPos * uiDividend ) );
                vIndexEnd.append( readableBp( ( xCoord.uiIndexPos + xCoord.uiIndexSize ) * uiDividend ) );
                vScreenPos.append( xCoord.uiScreenPos + xCoord.uiScreenSize );
                vValue.append( uiVal );

                if( uiX + 1 == vAxisCords[ uiI ].size( ) )
                {
                    // zero position at end
                    vIndexStart.append( readableBp( ( xCoord.uiIndexPos + xCoord.uiIndexSize ) * uiDividend ) );
                    vIndexEnd.append( readableBp( ( xCoord.uiIndexPos + xCoord.uiIndexSize ) * uiDividend ) );
                    vScreenPos.append( xCoord.uiScreenPos + xCoord.uiScreenSize );
                    vValue.append( vvMinMaxTracks[ uiI ][ 0 ] );
                }
            }

            vChrs.append( substringChr( sChr ) );
            vScreenPoss.append( vScreenPos );
            vIndexStarts.append( vIndexStart );
            vIndexEnds.append( vIndexEnd );
            vValues.append( vValue );
            vColors.append( vColorPaletteAnnotation[ uiCnt % vColorPaletteAnnotation.size( ) ] );
            vNames.append( vActiveCoverage[ uiI ][ uiId ] );

            ++uiCnt;
        }

        xTracksCDS[ uiI ] = pybind11::dict( "chrs"_a = vChrs,
                                            "screen_pos"_a = vScreenPoss,

                                            "index_start"_a = vIndexStarts,
                                            "index_end"_a = vIndexEnds,

                                            "values"_a = vValues,

                                            "colors"_a = vColors,
                                            "names"_a = vNames );
    }
    END_RETURN;
}

bool PartialQuarry::setTrackExport( )
{
    size_t uiDividend = getValue<size_t>( { "dividend" } );
    for( size_t uiI = 0; uiI < 2; uiI++ )
    {
        vTrackExport[ uiI ].clear( );
        vTrackExport[ uiI ].reserve( vAxisCords[ uiI ].size( ) );
        vTrackExportNames[ uiI ].clear( );
        vTrackExportNames[ uiI ].reserve( vvCoverageValues[ uiI ].size( ) );

        for( size_t uiId = 0; uiId < vvCoverageValues[ uiI ].size( ); uiI++ )
            vTrackExportNames[ uiI ].push_back( vActiveCoverage[ uiI ][ uiId ] );

        for( size_t uiX = 0; uiX < vAxisCords[ uiI ].size( ); uiX++ )
        {
            CANCEL_RETURN;
            std::vector<double> vValues;
            for( size_t uiId = 0; uiId < vvCoverageValues[ uiI ].size( ); uiI++ )
                vValues.push_back( vvCoverageValues[ uiI ][ uiId ][ uiX ] );

            auto& xCoord = vAxisCords[ uiI ][ uiX ];
            std::string sChromName = vActiveChromosomes[ uiI ][ xCoord.uiChromosome ].sName;
            vTrackExport[ uiI ].emplace_back( sChromName,
                                              xCoord.uiIndexPos * uiDividend,
                                              ( xCoord.uiIndexPos + xCoord.uiIndexSize ) * uiDividend,
                                              vValues );
        }
    }
    END_RETURN;
}

#if 0
bool PartialQuarry::setRankedSlicesCDS( )
{
    std::array<std::vector<size_t>, 2> vSorted;

    for( size_t uiI = 0; uiI < 2; uiI++ )
    {
        vSorted[ uiI ].reserve( vvCombinedCoverageValues[ uiI ].size( ) );
        for( size_t uiX = 0; uiX < vvCombinedCoverageValues[ uiI ].size( ); uiX++ )
            vSorted[ uiI ].push_back( uiX );
        std::sort( vSorted[ uiI ].begin( ), vSorted[ uiI ].end( ), [ & ]( size_t uiA, size_t uiB ) {
            return vvCombinedCoverageValues[ uiI ][ uiA ] < vvCombinedCoverageValues[ uiI ][ uiB ];
        } );
    }

    using namespace pybind11::literals;
    pybind11::gil_scoped_acquire acquire;
    size_t uiDividend = getValue<size_t>( { "dividend" } );

    for( size_t uiI = 0; uiI < 2; uiI++ )
    {
        pybind11::list vChrs;
        pybind11::list vIndexStart;
        pybind11::list vIndexEnd;
        pybind11::list vXs;
        pybind11::list vYs;
        pybind11::list vColors;
        pybind11::list vScoreA;
        pybind11::list vScoreB;

        for( size_t uiX = 0; uiX < vvCombinedCoverageValues[ uiI ].size( ); uiX++ )
        {
            CANCEL_RETURN;
            auto& xCoord = vAxisCords[ uiI ][ vSorted[ uiI ][ uiX ] ];
            auto uiVal = vvCombinedCoverageValues[ uiI ][ vSorted[ uiI ][ uiX ] ];
            auto uiA = vvFlatCoverageValues[ uiI ][ vSorted[ uiI ][ uiX ] ][ 0 ];
            auto uiB = vvFlatCoverageValues[ uiI ][ vSorted[ uiI ][ uiX ] ][ 1 ];

            std::string sChromName = vActiveChromosomes[ uiI ][ xCoord.uiChromosome ].sName;
            vChrs.append( substringChr( sChromName ) );
            vIndexStart.append( readableBp( xCoord.uiIndexPos * uiDividend ) );
            vIndexEnd.append( readableBp( ( xCoord.uiIndexPos + xCoord.uiIndexSize ) * uiDividend ) );

            vXs.append( uiX );
            vYs.append( uiVal );

            vColors.append( xCoord.bFiltered ? "grey" : ( uiI == 0 ? "#0072B2" : "#D55E00" ) );

            vScoreA.append( uiA );
            vScoreB.append( uiB );
        }

        vRankedSliceCDS[ uiI ] = pybind11::dict( "chrs"_a = vChrs,

                                                 "index_start"_a = vIndexStart,
                                                 "index_end"_a = vIndexEnd,

                                                 "xs"_a = vXs,
                                                 "ys"_a = vYs,

                                                 "colors"_a = vColors,

                                                 "score_a"_a = vScoreA,
                                                 "score_b"_a = vScoreB );
    }

    END_RETURN;
}
#endif

const decltype( PartialQuarry::vTrackExport[ 0 ] ) PartialQuarry::getTrackExport( bool bXAxis )
{
    update( NodeNames::TrackExport );
    return vTrackExport[ bXAxis ? 0 : 1 ];
}

const std::vector<std::string> PartialQuarry::getTrackExportNames( bool bXAxis )
{
    update( NodeNames::TrackExport );
    return vTrackExportNames[ bXAxis ? 0 : 1 ];
}


const pybind11::dict PartialQuarry::getTracks( bool bXAxis )
{
    update( NodeNames::Tracks );
    return xTracksCDS[ bXAxis ? 0 : 1 ];
}

#if 0
const pybind11::dict PartialQuarry::getRankedSlices( bool bXAxis )
{
    update( NodeNames::RankedSlicesCDS );
    return vRankedSliceCDS[ bXAxis ? 0 : 1 ];
}
#endif

const std::array<double, 2> PartialQuarry::getMinMaxTracks( bool bXAxis )
{
    update( NodeNames::Tracks );
    return vvMinMaxTracks[ bXAxis ? 0 : 1 ];
}


void PartialQuarry::regCoverage( )
{
    registerNode( NodeNames::ActiveCoverage,
                  ComputeNode{ .sNodeName = "active_coverage",
                               .fFunc = &PartialQuarry::setActiveCoverage,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "coverage", "in_column" },
                                                     { "coverage", "in_row" },
                                                     { "coverage", "cov_column_a" },
                                                     { "coverage", "cov_column_b" },
                                                     { "coverage", "cov_row_a" },
                                                     { "coverage", "cov_row_b" },
                                                     { "replicates", "in_column" },
                                                     { "replicates", "in_row" },
                                                     { "replicates", "cov_column_a" },
                                                     { "replicates", "cov_column_b" },
                                                     { "replicates", "cov_row_a" },
                                                     { "replicates", "cov_row_b" },
                                                     { "settings", "normalization", "normalize_by" },
                                                     { "coverage", "list" },
                                                     { "replicates", "list" } },
                               .vSessionsIncomingInPrevious = {} } );

    registerNode( NodeNames::CoverageValues,
                  ComputeNode{ .sNodeName = "coverage_values",
                               .fFunc = &PartialQuarry::setCoverageValues,
                               .vIncomingFunctions = { NodeNames::ActiveCoverage, NodeNames::AxisCoords,
                                                       NodeNames::IntersectionType, NodeNames::Symmetry },
                               .vIncomingSession = { { "settings", "filters", "mapping_q", "val_min" },
                                                     { "settings", "filters", "mapping_q", "val_max" },
                                                     { "settings", "replicates", "coverage_get_max_col" },
                                                     { "settings", "replicates", "coverage_get_max_row" },
                                                     { "settings", "replicates", "coverage_get_max_bin_size", "val" },
                                                     { "settings", "normalization", "min_interactions", "val" },
                                                     { "coverage", "by_name" },
                                                     { "replicates", "by_name" },
                                                     { "settings", "filters", "incomplete_alignments" } },
                               .vSessionsIncomingInPrevious = { { "dividend" } } } );

    registerNode( NodeNames::Tracks,
                  ComputeNode{ .sNodeName = "coverage_tracks",
                               .fFunc = &PartialQuarry::setTracks,
                               .vIncomingFunctions = { NodeNames::LCS, NodeNames::CombinedCoverageValues,
                                                       NodeNames::AnnotationColors },
                               .vIncomingSession = { { "settings", "normalization", "display_ice_remainder" } },
                               .vSessionsIncomingInPrevious = { { "dividend" } } } );

    registerNode( NodeNames::TrackExport,
                  ComputeNode{ .sNodeName = "track_export",
                               .fFunc = &PartialQuarry::setTrackExport,
                               .vIncomingFunctions = { NodeNames::Tracks },
                               .vIncomingSession = { },
                               .vSessionsIncomingInPrevious = { { "dividend" } } } );
}


} // namespace cm