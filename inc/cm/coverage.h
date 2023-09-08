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

std::tuple<size_t, int64_t, size_t, size_t> PartialQuarry::makeHeapTuple( bool bCol, bool bSymPart, const size_t uiFrom,
                                                                          const size_t uiTo, int64_t iDataSetId,
                                                                          size_t uiStart, size_t uiEnd )
{
    const coordinate_t uiXMin = bCol != bSymPart ? uiFrom : uiStart;
    const coordinate_t uiXMax = bCol != bSymPart ? uiTo : uiEnd;
    const coordinate_t uiYMin = bCol != bSymPart ? uiStart : uiFrom;
    const coordinate_t uiYMax = bCol != bSymPart ? uiEnd : uiTo;

    const size_t uiCount = indexCount( iDataSetId, uiYMin, uiXMin, uiYMax, uiXMax, !bSymPart );

    return std::make_tuple( uiCount, iDataSetId, uiStart, uiEnd );
}

size_t PartialQuarry::getMaxCoverageFromRepl( const size_t uiChromId, const size_t uiFrom, const size_t uiTo,
                                              const size_t uiRepl, size_t uiCoverageGetMaxBinSize, bool bCol,
                                              bool bSymPart )
{
    // score, datasetId, start, end
    std::vector<std::tuple<size_t, int64_t, size_t, size_t>> vHeap;
    for( size_t uiI = 0; uiI < vActiveChromosomes[ bCol ? 1 : 0 ].size( ); uiI++ )
    {
        size_t uiActualContigId = vActiveChromosomes[ bCol ? 1 : 0 ][ uiI ].uiActualContigId;
        const size_t uiDatasetId = getDatasetIdfromReplAndChr( uiRepl, bCol != bSymPart ? uiChromId : uiActualContigId,
                                                               bCol != bSymPart ? uiActualContigId : uiChromId );
        if( uiDatasetId != std::numeric_limits<size_t>::max( ) )
            vHeap.push_back( makeHeapTuple( bCol, bSymPart, uiFrom, uiTo, uiDatasetId, 0,
                                            vActiveChromosomes[ bCol ? 1 : 0 ][ uiI ].uiLength ) );
    }

    std::make_heap( vHeap.begin( ), vHeap.end( ) );

    while( vHeap.size( ) > 0 && std::get<3>( vHeap.back( ) ) - std::get<2>( vHeap.back( ) ) >
                                    std::max( uiCoverageGetMaxBinSize, (size_t)1 ) )
    {
        auto xFront = vHeap.back( );
        vHeap.pop_back( );

        size_t uiCenter = ( std::get<2>( xFront ) + std::get<3>( xFront ) ) / 2;
        vHeap.push_back(
            makeHeapTuple( bCol, bSymPart, uiFrom, uiTo, std::get<1>( xFront ), uiCenter, std::get<3>( xFront ) ) );
        std::push_heap( vHeap.begin( ), vHeap.end( ) );

        vHeap.push_back(
            makeHeapTuple( bCol, bSymPart, uiFrom, uiTo, std::get<1>( xFront ), std::get<2>( xFront ), uiCenter ) );
        std::push_heap( vHeap.begin( ), vHeap.end( ) );


        std::pop_heap( vHeap.begin( ), vHeap.end( ) );
    }
    if( vHeap.size( ) > 0 )
        return std::get<0>( vHeap.back( ) );

    return 0;
}
size_t PartialQuarry::getMaxCoverageFromRepl( const AxisCoord& xCoords, const size_t uiRep,
                                              size_t uiCoverageGetMaxBinSize, bool bCol, bool bSymPart )
{

    return getMaxCoverageFromRepl( vActiveChromosomes[ bCol ? 0 : 1 ][ xCoords.uiChromosome ].uiActualContigId,
                                   xCoords.uiIndexPos, xCoords.uiIndexPos + xCoords.uiIndexSize, uiRep,
                                   uiCoverageGetMaxBinSize, bCol, bSymPart );
}

size_t PartialQuarry::getCoverageFromRepl( const size_t uiChromId, const size_t uiFrom, const size_t uiTo,
                                           const size_t uiRepl, bool bCol, bool bSymPart )
{
    size_t uiRet = 0;
    for( size_t uiI = 0; uiI < vActiveChromosomes[ bCol ? 1 : 0 ].size( ); uiI++ )
    {
        size_t uiActualContigId = vActiveChromosomes[ bCol ? 1 : 0 ][ uiI ].uiActualContigId;
        const size_t uiDatasetId = getDatasetIdfromReplAndChr( uiRepl, bCol != bSymPart ? uiChromId : uiActualContigId,
                                                               bCol != bSymPart ? uiActualContigId : uiChromId );

        if( uiDatasetId != std::numeric_limits<size_t>::max( ) )
        {
            const size_t uiXMin = bCol != bSymPart ? uiFrom : 0;
            const size_t uiXMax = bCol != bSymPart ? uiTo : vActiveChromosomes[ bCol ? 1 : 0 ][ uiI ].uiLength;
            const size_t uiYMin = bCol != bSymPart ? 0 : uiFrom;
            const size_t uiYMax = bCol != bSymPart ? vActiveChromosomes[ bCol ? 1 : 0 ][ uiI ].uiLength : uiTo;
            uiRet += indexCount( uiDatasetId, uiYMin, uiXMin, uiYMax, uiXMax, !bSymPart );
        }
    }


    return uiRet;
}
size_t PartialQuarry::getCoverageFromRepl( const AxisCoord& xCoords, const size_t uiRepl, bool bCol, bool bSymPart )
{
    return getCoverageFromRepl( vActiveChromosomes[ bCol ? 0 : 1 ][ xCoords.uiChromosome ].uiActualContigId,
                                xCoords.uiIndexPos, xCoords.uiIndexPos + xCoords.uiIndexSize, uiRepl, bCol, bSymPart );
}

bool PartialQuarry::setCoverageValues( )
{

    size_t uiMinuend = getValue<size_t>( { "settings", "normalization", "min_interactions", "val" } );


    for( size_t uiJ = 0; uiJ < 2; uiJ++ )
    {
        vvCoverageValues[ uiJ ].clear( );
        vvCoverageValues[ uiJ ].reserve( vActiveCoverage[ uiJ ].size( ) );

        for( auto& sRep : vActiveCoverage[ uiJ ] )
        {
            vvCoverageValues[ uiJ ].emplace_back( );
            vvCoverageValues[ uiJ ].back( ).reserve( vAxisCords[ uiJ ].size( ) );
            size_t uiFstDatasetId = getValue<json>( { "coverage", "by_name", sRep, "first_dataset_id" } );

            for( AxisCoord& xCoords : vAxisCords[ uiJ ] )
            {
                size_t uiVal;
                if( xCoords.uiChromosome != std::numeric_limits<size_t>::max( ) )
                {
                    CANCEL_RETURN;

                    int64_t iDataSetId =
                        uiFstDatasetId + vActiveChromosomes[ uiJ ][ xCoords.uiChromosome ].uiActualContigId;

                    uiVal = index1DCount( iDataSetId, xCoords.uiIndexPos, xCoords.uiIndexPos + xCoords.uiIndexSize,
                                          uiJ == 0 );
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
    const bool bCenterTracks = getValue<bool>( { "settings", "interface", "center_tracks_on_bins" } );
    const bool bZeroTracksAtEnd = getValue<bool>( { "settings", "interface", "zero_track_at_ends" } );
    const bool bConnectOverBorder = getValue<bool>( { "settings", "interface", "connect_tracks_over_contig_borders" } );
    const size_t uiMaxChar = getValue<size_t>( { "settings", "interface", "axis_label_max_char", "val" } );

    size_t uiDividend = getValue<size_t>( { "dividend" } );
    for( size_t uiI = 0; uiI < 2; uiI++ )
    {
        pybind11::list vChrs;
        pybind11::list vScreenPoss;
        pybind11::list vIndexStarts;
        pybind11::list vIndexEnds;
        pybind11::list vValues;
        pybind11::list vColors;
        pybind11::list vNames;

        size_t uiCnt = 0;

        for( size_t uiJ = 0; uiJ < vTrackExportNames[ uiI ].size( ); uiJ++ )
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
                if( sChr != "" && sChr != sChromName && !bConnectOverBorder )
                {
                    vChrs.append( substringChr( sChr, uiMaxChar ) );

                    vScreenPoss.append( vScreenPos );
                    vScreenPos = pybind11::list( );

                    vIndexStarts.append( vIndexStart );
                    vIndexStart = pybind11::list( );

                    vIndexEnds.append( vIndexEnd );
                    vIndexEnd = pybind11::list( );

                    vValues.append( vValue );
                    vValue = pybind11::list( );

                    vColors.append( vColorPaletteAnnotation[ uiCnt % vColorPaletteAnnotation.size( ) ] );

                    vNames.append( vTrackExportNames[ uiI ][ uiJ ] );
                }
                if( bZeroTracksAtEnd && uiX == 0 )
                {
                    // zero position at start
                    vIndexStart.append( readableBp( xCoord.uiIndexPos * uiDividend ) );
                    vIndexEnd.append( readableBp( xCoord.uiIndexPos * uiDividend ) );
                    vScreenPos.append( xCoord.uiScreenPos );
                    vValue.append( 0 );
                }

                sChr = sChromName;
                double uiVal = vTrackPrecursor[ uiI ][ uiX ].vValues[ uiJ ];

                if( bCenterTracks )
                {
                    // center
                    vIndexStart.append( readableBp( xCoord.uiIndexPos * uiDividend ) );
                    vIndexEnd.append( readableBp( ( xCoord.uiIndexPos + xCoord.uiIndexSize ) * uiDividend ) );
                    vScreenPos.append( xCoord.uiScreenPos + xCoord.uiScreenSize / 2 );
                    vValue.append( uiVal );
                }
                else
                {
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
                }

                if( bZeroTracksAtEnd && uiX + 1 == vAxisCords[ uiI ].size( ) )
                {
                    // zero position at end
                    vIndexStart.append( readableBp( ( xCoord.uiIndexPos + xCoord.uiIndexSize ) * uiDividend ) );
                    vIndexEnd.append( readableBp( ( xCoord.uiIndexPos + xCoord.uiIndexSize ) * uiDividend ) );
                    vScreenPos.append( xCoord.uiScreenPos + xCoord.uiScreenSize );
                    vValue.append( 0 );
                }
            }

            vChrs.append( substringChr( sChr, uiMaxChar ) );
            vScreenPoss.append( vScreenPos );
            vIndexStarts.append( vIndexStart );
            vIndexEnds.append( vIndexEnd );
            vValues.append( vValue );
            vColors.append( vColorPaletteAnnotation[ uiCnt % vColorPaletteAnnotation.size( ) ] );
            vNames.append( vTrackExportNames[ uiI ][ uiJ ] );

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

bool PartialQuarry::setTrackPrecursor( )
{

    const std::string sNorm = getValue<std::string>( { "settings", "normalization", "normalize_by" } );
    const bool bGridSeqNormDisp =
        sNorm == "grid-seq" && getValue<bool>( { "settings", "normalization", "grid_seq_display_background" } );
    const bool bRadiclNormDisp =
        sNorm == "radicl-seq" && getValue<bool>( { "settings", "normalization", "radicl_seq_display_coverage" } );
    bool bIsCol;
    if( bGridSeqNormDisp )
        bIsCol = getValue<bool>( { "settings", "normalization", "grid_seq_axis_is_column" } );
    else if( bRadiclNormDisp )
        bIsCol = getValue<bool>( { "settings", "normalization", "radicl_seq_axis_is_column" } );
    const bool bIceShowBias = getValue<bool>( { "settings", "normalization", "ice_show_bias" } );

    for( size_t uiI = 0; uiI < 2; uiI++ )
    {
        vvMinMaxTracks[ uiI ][ 0 ] = std::numeric_limits<double>::max( );
        vvMinMaxTracks[ uiI ][ 1 ] = std::numeric_limits<double>::min( );
        const bool bDoIce = sNorm == "ice" && bIceShowBias;

        for( size_t uiX = 0; uiX < vAxisCords[ uiI ].size( ); uiX++ )
        {
            CANCEL_RETURN;
            for( size_t uiId = 0; uiId < vvCoverageValues[ uiI ].size( ); uiId++ )
            {
                auto uiVal = (double)vvCoverageValues[ uiI ][ uiId ][ uiX ];
                vvMinMaxTracks[ uiI ][ 0 ] = std::min( vvMinMaxTracks[ uiI ][ 0 ], uiVal );
                vvMinMaxTracks[ uiI ][ 1 ] = std::max( vvMinMaxTracks[ uiI ][ 1 ], uiVal );
            }
            if( ( bRadiclNormDisp || bGridSeqNormDisp ) && ( bIsCol == ( uiI == 0 ) ) )
            {
                double uiVal;
                if( bGridSeqNormDisp )
                    uiVal = (double)vBackgroundGridSeq[ uiX ];
                else if( bRadiclNormDisp )
                    uiVal = std::min( (double)vRadiclSeqCoverage[ 0 ][ uiX ][ 0 ],
                                      (double)vRadiclSeqNumNonEmptyBins[ 0 ][ uiX ][ 0 ] );
                vvMinMaxTracks[ uiI ][ 0 ] = std::min( vvMinMaxTracks[ uiI ][ 0 ], uiVal );
                if( bRadiclNormDisp )
                    uiVal = std::max( (double)vRadiclSeqCoverage[ 0 ][ uiX ][ 0 ],
                                      (double)vRadiclSeqNumNonEmptyBins[ 0 ][ uiX ][ 0 ] );
                vvMinMaxTracks[ uiI ][ 1 ] = std::max( vvMinMaxTracks[ uiI ][ 1 ], uiVal );
            }

            if( vFlat4C[ 1 - uiI ].size( ) > 0 )
            {
                auto uiVal = (double)vFlat4C[ 1 - uiI ][ uiX ];
                vvMinMaxTracks[ uiI ][ 0 ] = std::min( vvMinMaxTracks[ uiI ][ 0 ], uiVal );
                vvMinMaxTracks[ uiI ][ 1 ] = std::max( vvMinMaxTracks[ uiI ][ 1 ], uiVal );
            }
        }

        for( size_t uiX = 0; uiX < vSampleAxisCoords[ uiI ].size( ); uiX++ )
            for( size_t uiY = 0; uiY < 2; uiY++ )
                if( bDoIce && vSampleAxisCoords[ uiI ][ uiY ].uiIdx != SAMPLE_COORD &&
                    vIceSliceBias[ 0 ][ uiI ][ uiY ].size( ) > 0 && vInGroup[ uiY ].size( ) > 0 )
                {
                    auto uiVal = vIceSliceBias[ 0 ][ uiI ][ uiY ][ uiX ];
                    vvMinMaxTracks[ uiI ][ 0 ] = std::min( vvMinMaxTracks[ uiI ][ 0 ], uiVal );
                    vvMinMaxTracks[ uiI ][ 1 ] = std::max( vvMinMaxTracks[ uiI ][ 1 ], uiVal );
                }


        vTrackExportNames[ uiI ].clear( );
        vTrackExportNames[ uiI ].reserve( 20 );

        for( size_t uiId = 0; uiId < vvCoverageValues[ uiI ].size( ); uiId++ )
            vTrackExportNames[ uiI ].push_back( vActiveCoverage[ uiI ][ uiId ] );
        if( bGridSeqNormDisp && ( bIsCol == ( uiI == 0 ) ) )
            vTrackExportNames[ uiI ].push_back( "background RNA associated elements" );
        if( bRadiclNormDisp && ( bIsCol == ( uiI == 0 ) ) )
        {
            vTrackExportNames[ uiI ].push_back( "coverage A" );
            vTrackExportNames[ uiI ].push_back( "coverage B" );
            vTrackExportNames[ uiI ].push_back( "non-empty bins A" );
            vTrackExportNames[ uiI ].push_back( "non-empty bins B" );
        }
        if( vFlat4C[ 1 - uiI ].size( ) > 0 )
            vTrackExportNames[ uiI ].push_back( "Virtual 4C" );
        for( size_t uiY = 0; uiY < 2; uiY++ )
            if( bIceShowBias && vIceSliceBias[ 0 ][ uiI ][ uiY ].size( ) > 0 && vInGroup[ uiY ].size( ) > 0 )
                vTrackExportNames[ uiI ].push_back( std::string( "ICE Bias " ) + ( uiY == 0 ? "A" : "B" ) );


        vTrackPrecursor[ uiI ].clear( );
        vTrackPrecursor[ uiI ].reserve( vAxisCords[ uiI ].size( ) );


        for( size_t uiX = 0; uiX < vAxisCords[ uiI ].size( ); uiX++ )
        {
            CANCEL_RETURN;
            vTrackPrecursor[ uiI ].push_back( Track{ /*.xCoord=*/vAxisCords[ uiI ][ uiX ], /*.vValues=*/{} } );

            for( size_t uiId = 0; uiId < vvCoverageValues[ uiI ].size( ); uiId++ )
                vTrackPrecursor[ uiI ].back( ).vValues.push_back( (double)vvCoverageValues[ uiI ][ uiId ][ uiX ] );

            if( bGridSeqNormDisp && ( bIsCol == ( uiI == 0 ) ) )
                vTrackPrecursor[ uiI ].back( ).vValues.push_back( (double)vBackgroundGridSeq[ uiX ] );
            if( bRadiclNormDisp && ( bIsCol == ( uiI == 0 ) ) )
            {
                vTrackPrecursor[ uiI ].back( ).vValues.push_back( (double)vRadiclSeqCoverage[ 0 ][ uiX ][ 0 ] );
                vTrackPrecursor[ uiI ].back( ).vValues.push_back( (double)vRadiclSeqCoverage[ 0 ][ uiX ][ 1 ] );
                vTrackPrecursor[ uiI ].back( ).vValues.push_back( (double)vRadiclSeqNumNonEmptyBins[ 0 ][ uiX ][ 0 ] );
                vTrackPrecursor[ uiI ].back( ).vValues.push_back( (double)vRadiclSeqNumNonEmptyBins[ 0 ][ uiX ][ 1 ] );
            }

            if( vFlat4C[ 1 - uiI ].size( ) > 0 )
                vTrackPrecursor[ uiI ].back( ).vValues.push_back( (double)vFlat4C[ 1 - uiI ][ uiX ] );

            for( size_t uiY = 0; uiY < 2; uiY++ )
                if( bIceShowBias && vIceSliceBias[ 0 ][ uiI ][ uiY ].size( ) > 0 && vInGroup[ uiY ].size( ) > 0 )
                    vTrackPrecursor[ uiI ].back( ).vValues.push_back( (double)vIceSliceBias[ 0 ][ uiI ][ uiY ][ uiX ] );
        }
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

        for( size_t uiX = 0; uiX < vAxisCords[ uiI ].size( ); uiX++ )
        {
            CANCEL_RETURN;
            auto& xCoord = vAxisCords[ uiI ][ uiX ];
            std::string sChromName = vActiveChromosomes[ uiI ][ xCoord.uiChromosome ].sName;
            vTrackExport[ uiI ].emplace_back( sChromName,
                                              xCoord.uiIndexPos * uiDividend,
                                              ( xCoord.uiIndexPos + xCoord.uiIndexSize ) * uiDividend,
                                              vTrackPrecursor[ uiI ][ uiX ].vValues );
        }
    }
    END_RETURN;
}

bool PartialQuarry::setRankedSlicesCDS( )
{
    const std::string sAnno = getValue<std::string>( { "settings", "normalization", "grid_seq_annotation" } );
    auto uiFistAnnoIdx =
        hasValue( { "annotation", "by_name", sAnno } ) ? getValue<size_t>( { "annotation", "by_name", sAnno } ) : 0;
    const uint32_t uiDividend = getValue<uint32_t>( { "dividend" } );
    const size_t uiMaxChar = getValue<size_t>( { "settings", "interface", "axis_label_max_char", "val" } );
    std::array<std::vector<size_t>, 2> vSorted;

    for( size_t uiI = 0; uiI < 2; uiI++ )
    {
        vSorted[ uiI ].reserve( vGridSeqAnnoCoverage.size( ) );
        for( size_t uiX = 0; uiX < vGridSeqAnnoCoverage.size( ); uiX++ )
            vSorted[ uiI ].push_back( uiX );
        std::sort( vSorted[ uiI ].begin( ), vSorted[ uiI ].end( ), [ & ]( size_t uiA, size_t uiB ) {
            return vGridSeqAnnoCoverage[ uiA ][ uiI ] < vGridSeqAnnoCoverage[ uiB ][ uiI ];
        } );
    }

    using namespace pybind11::literals;
    pybind11::gil_scoped_acquire acquire;

    for( size_t uiI = 0; uiI < 2; uiI++ )
    {
        pybind11::list vChrs;
        pybind11::list vAnnoDesc;
        pybind11::list vSampleId;
        pybind11::list vAnnoIdx;
        pybind11::list vIndexStart;
        pybind11::list vIndexEnd;
        pybind11::list vXs;
        pybind11::list vYs;
        pybind11::list vColors;
        pybind11::list vScoreA;
        pybind11::list vScoreB;

        for( size_t uiX = 0; uiX < vGridSeqAnnoCoverage.size( ); uiX++ )
        {
            CANCEL_RETURN;
            auto uiVal = vGridSeqAnnoCoverage[ vSorted[ uiI ][ uiX ] ][ uiI ];

            const AnnoCoord& rSample = vGridSeqSamples[ vSorted[ uiI ][ uiX ] ];
            const std::string& sChromName = this->vActiveChromosomes[ 0 ][ rSample.uiChromosome ].sName;
            const size_t uiChromId = this->vActiveChromosomes[ 0 ][ rSample.uiChromosome ].uiActualContigId;

            const auto rIntervalIt = pIndices->vAnno.get( uiFistAnnoIdx + uiChromId, rSample.uiAnnoId );

            vChrs.append( substringChr( sChromName, uiMaxChar ) );
            vAnnoDesc.append( pIndices->vAnno.desc( *rIntervalIt ) );
            vSampleId.append( vSorted[ uiI ][ uiX ] );
            vAnnoIdx.append( rSample.uiAnnoId );
            vIndexStart.append( readableBp( rSample.uiIndexPos * uiDividend ) );
            vIndexEnd.append( readableBp( ( rSample.uiIndexPos + rSample.uiIndexSize ) * uiDividend ) );

            vXs.append( uiX );
            vYs.append( uiVal );

            vColors.append( vGridSeqFiltered[ vSorted[ uiI ][ uiX ] ][ uiI ] ? ( uiI == 0 ? "#0072B2" : "#D55E00" )
                                                                             : "grey" );
        }

        vRankedSliceCDS[ uiI ] = pybind11::dict( "chrs"_a = vChrs,

                                                 "index_start"_a = vIndexStart,
                                                 "index_end"_a = vIndexEnd,

                                                 "anno_desc"_a = vAnnoDesc,
                                                 "sample_id"_a = vSampleId,
                                                 "anno_idx"_a = vAnnoIdx,

                                                 "xs"_a = vXs,
                                                 "ys"_a = vYs,

                                                 "colors"_a = vColors );
    }

    END_RETURN;
}

const decltype( PartialQuarry::vTrackExport[ 0 ] )
PartialQuarry::getTrackExport( bool bXAxis, const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::TrackExport, fPyPrint );
    return vTrackExport[ bXAxis ? 0 : 1 ];
}

const std::vector<std::string>
PartialQuarry::getTrackExportNames( bool bXAxis, const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::TrackExport, fPyPrint );
    return vTrackExportNames[ bXAxis ? 0 : 1 ];
}


const pybind11::dict PartialQuarry::getTracks( bool bXAxis, const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::Tracks, fPyPrint );
    return xTracksCDS[ bXAxis ? 0 : 1 ];
}

const pybind11::dict PartialQuarry::getRankedSlices( bool bXAxis,
                                                     const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::RankedSlicesCDS, fPyPrint );
    return vRankedSliceCDS[ bXAxis ? 0 : 1 ];
}

const std::array<double, 2> PartialQuarry::getMinMaxTracks( bool bXAxis,
                                                            const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::Tracks, fPyPrint );
    return vvMinMaxTracks[ bXAxis ? 0 : 1 ];
}


void PartialQuarry::regCoverage( )
{
    registerNode( NodeNames::ActiveCoverage,
                  ComputeNode{ /*.sNodeName =*/"active_coverage",
                               /*.fFunc =*/&PartialQuarry::setActiveCoverage,
                               /*.vIncomingFunctions =*/{ },
                               /*.vIncomingSession =*/
                               { { "coverage", "in_column" },
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
                               /*.vSessionsIncomingInPrevious =*/{ },
                               /*bHidden =*/false } );

    registerNode( NodeNames::CoverageValues,
                  ComputeNode{ /*.sNodeName =*/"coverage_values",
                               /*.fFunc =*/&PartialQuarry::setCoverageValues,
                               /*.vIncomingFunctions =*/
                               { NodeNames::ActiveCoverage, NodeNames::AxisCoords, NodeNames::IntersectionType,
                                 NodeNames::Symmetry, NodeNames::MappingQuality, NodeNames::Directionality },
                               /*.vIncomingSession =*/
                               { { "settings", "replicates", "coverage_get_max_col" },
                                 { "settings", "replicates", "coverage_get_max_row" },
                                 { "settings", "replicates", "coverage_get_max_bin_size", "val" },
                                 { "settings", "normalization", "min_interactions", "val" },
                                 { "coverage", "by_name" },
                                 { "replicates", "by_name" } },
                               /*.vSessionsIncomingInPrevious =*/
                               { { "settings", "filters", "incomplete_alignments" }, { "dividend" } },
                               /*bHidden =*/false } );

    registerNode( NodeNames::Tracks,
                  ComputeNode{ /*.sNodeName =*/"coverage_tracks",
                               /*.fFunc =*/&PartialQuarry::setTracks,
                               /*.vIncomingFunctions =*/
                               { NodeNames::TrackPrecursor },
                               /*.vIncomingSession =*/
                               { { "settings", "interface", "zero_track_at_ends" },
                                 { "settings", "interface", "center_tracks_on_bins" },
                                 { "settings", "interface", "connect_tracks_over_contig_borders" },
                                 { "settings", "interface", "axis_label_max_char", "val" } },
                               /*.vSessionsIncomingInPrevious =*/
                               { { "dividend" } },
                               /*bHidden =*/false } );
    registerNode(
        NodeNames::TrackPrecursor,
        ComputeNode{ /*.sNodeName =*/"coverage_track_precurser",
                     /*.fFunc =*/&PartialQuarry::setTrackPrecursor,
                     /*.vIncomingFunctions =*/
                     { NodeNames::LCS, NodeNames::AnnotationColors, NodeNames::CoverageValues, NodeNames::Flat4C },
                     /*.vIncomingSession =*/
                     { { "settings", "normalization", "display_ice_remainder" },
                       { "settings", "normalization", "grid_seq_display_background" },
                       { "settings", "normalization", "ice_show_bias" },
                       { "settings", "normalization", "radicl_seq_display_coverage" } },
                     /*.vSessionsIncomingInPrevious =*/
                     { { "dividend" },
                       { "settings", "normalization", "grid_seq_axis_is_column" },
                       { "settings", "normalization", "radicl_seq_axis_is_column" },
                       { "settings", "normalization", "normalize_by" } },
                     /*bHidden =*/false } );

    registerNode( NodeNames::TrackExport,
                  ComputeNode{ /*.sNodeName =*/"track_export",
                               /*.fFunc =*/&PartialQuarry::setTrackExport,
                               /*.vIncomingFunctions =*/{ NodeNames::Tracks },
                               /*.vIncomingSession =*/{ },
                               /*.vSessionsIncomingInPrevious =*/
                               { { "settings", "normalization", "radicl_seq_display_coverage" },
                                 { "settings", "normalization", "radicl_seq_axis_is_column" },
                                 { "settings", "normalization", "ice_show_bias" },
                                 { "settings", "normalization", "normalize_by" },
                                 { "settings", "normalization", "grid_seq_axis_is_column" },
                                 { "settings", "normalization", "grid_seq_display_background" },
                                 { "dividend" } },
                               /*bHidden =*/false } );

    registerNode(
        NodeNames::RankedSlicesCDS,
        ComputeNode{
            /*.sNodeName =*/"ranked_slices_cds",
            /*.fFunc =*/&PartialQuarry::setRankedSlicesCDS,
            /*.vIncomingFunctions =*/{ NodeNames::RnaAssociatedGenesFilter },
            /*.vIncomingSession =*/{ { "settings", "interface", "axis_label_max_char", "val" } },
            /*.vSessionsIncomingInPrevious =*/
            { { "annotation", "by_name" }, { "settings", "normalization", "grid_seq_annotation" }, { "dividend" } },
            /*bHidden =*/false } );
}


} // namespace cm