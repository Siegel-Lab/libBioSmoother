#include "cm/partial_quarry.h"
#include <cmath>

#pragma once

namespace cm
{

void PartialQuarry::setActiveCoverage( )
{
    size_t uiSize = this->xSession[ "coverage" ][ "list" ].size( ) + this->xSession[ "replicates" ][ "list" ].size( );
    for( size_t uiI = 0; uiI < 2; uiI++ )
    {
        vActiveCoverage[ uiI ].clear( );
        vActiveCoverage[ uiI ].reserve( uiSize );

        for( size_t uiJ = 0; uiJ < 3; uiJ++ )
        {
            vInGroupCoverage[ uiI ][ uiJ ].clear( );
            vInGroupCoverage[ uiI ][ uiJ ].reserve( uiSize );
        }
    }

    std::array<std::array<std::set<std::pair<std::string, bool>>, 3>, 2> vGroups;
    std::set<std::pair<std::string, bool>> xA{ };
    std::set<std::pair<std::string, bool>> xB{ };


    for( bool bB : { true, false } )
    {
        std::string sX = bB ? "coverage" : "replicates";

        for( std::string sY : { "in_column", "cov_column_a", "cov_column_b" } )
            for( auto& rV : this->xSession[ sX ][ sY ] )
            {
                std::pair<std::string, bool> xP( rV.get<std::string>( ), bB );
                xA.insert( xP );
                if( sY == "cov_column_a" )
                    vGroups[ 0 ][ 0 ].insert( xP );
                if( sY == "cov_column_b" )
                    vGroups[ 0 ][ 1 ].insert( xP );
                if( sY == "in_column" )
                    vGroups[ 0 ][ 2 ].insert( xP );
            }

        for( std::string sY : { "in_row", "cov_row_a", "cov_row_b" } )
            for( auto& rV : this->xSession[ sX ][ sY ] )
            {
                std::pair<std::string, bool> xP( rV.get<std::string>( ), bB );
                xB.insert( xP );

                if( sY == "cov_row_a" )
                    vGroups[ 1 ][ 0 ].insert( xP );
                if( sY == "cov_row_b" )
                    vGroups[ 1 ][ 1 ].insert( xP );
                if( sY == "in_row" )
                    vGroups[ 1 ][ 2 ].insert( xP );
            }
    }

    for( auto& rX : xA )
        vActiveCoverage[ 0 ].push_back( rX );
    for( auto& rX : xB )
        vActiveCoverage[ 1 ].push_back( rX );

    for( size_t uiI = 0; uiI < 2; uiI++ )
        for( size_t uiJ = 0; uiJ < vActiveCoverage[ uiI ].size( ); uiJ++ )
            for( size_t uiK = 0; uiK < 3; uiK++ )
                if( vGroups[ uiI ][ uiK ].count( vActiveCoverage[ uiI ][ uiJ ] ) > 0 )
                    vInGroupCoverage[ uiI ][ uiK ].push_back( uiJ );
}

void PartialQuarry::setCoverageValues( )
{
    for( size_t uiJ = 0; uiJ < 2; uiJ++ )
    {
        vvCoverageValues[ uiJ ].clear( );
        vvCoverageValues[ uiJ ].reserve( vActiveCoverage[ uiJ ].size( ) );
        for( auto& sRep : vActiveCoverage[ uiJ ] )
        {
            vvCoverageValues[ uiJ ].emplace_back( );
            vvCoverageValues[ uiJ ].back( ).reserve( vAxisCords[ uiJ ].size( ) );
            auto& xRep = this->xSession[ sRep.second ? "coverage" : "replicates" ][ "by_name" ][ sRep.first ];

            bool bHasMapQ = xRep[ "has_map_q" ];
            bool bHasMultiMap = xRep[ "has_multimapping" ];

            size_t uiMapQMin =
                255 - this->xSession[ "settings" ][ "filters" ][ "mapping_q" ][ "val_min" ].get<size_t>( );
            size_t uiMapQMax =
                255 - this->xSession[ "settings" ][ "filters" ][ "mapping_q" ][ "val_max" ].get<size_t>( );

            for( AxisCoord& xCoords : vAxisCords[ uiJ ] )
            {
                std::array<size_t, 2> vVals;
                if( xCoords.sChromosome != "" )
                    for( size_t uiI = 0; uiI < ( sRep.second ? 1 : 2 ); uiI++ )
                    {
                        int64_t iDataSetId;
                        if( sRep.second )
                            iDataSetId = xRep[ "ids" ][ xCoords.sChromosome ].get<int64_t>( );
                        else
                            iDataSetId = xRep[ "ids" ][ xCoords.sChromosome ][ uiJ ? "col" : "row" ].get<int64_t>( );

                        if( iDataSetId != -1 )
                        {
                            if( bHasMapQ && bHasMultiMap )
                                vVals[ uiI ] = xIndices.getIndex<2, 1>( )->count(
                                    iDataSetId,
                                    { xCoords.uiIndexPos, uiMapQMax },
                                    { xCoords.uiIndexPos + xCoords.uiSize, uiMapQMin },
                                    xIntersect,
                                    0 );
                            else if( !bHasMapQ && bHasMultiMap )
                                vVals[ uiI ] =
                                    xIndices.getIndex<1, 1>( )->count( iDataSetId,
                                                                       { xCoords.uiIndexPos },
                                                                       { xCoords.uiIndexPos + xCoords.uiSize },
                                                                       xIntersect,
                                                                       0 );
                            else if( bHasMapQ && !bHasMultiMap )
                                vVals[ uiI ] = xIndices.getIndex<2, 0>( )->count(
                                    iDataSetId,
                                    { xCoords.uiIndexPos, uiMapQMax },
                                    { xCoords.uiIndexPos + xCoords.uiSize, uiMapQMin },
                                    xIntersect,
                                    0 );
                            else // if(!bHasMapQ && !bHasMultiMap)
                                vVals[ uiI ] =
                                    xIndices.getIndex<1, 0>( )->count( iDataSetId,
                                                                       { xCoords.uiIndexPos },
                                                                       { xCoords.uiIndexPos + xCoords.uiSize },
                                                                       xIntersect,
                                                                       0 );
                        }
                    }
                else
                    vVals = { 0, 0 };

                if( sRep.second )
                    vvCoverageValues[ uiJ ].back( ).push_back( vVals[ 0 ] );
                else
                    vvCoverageValues[ uiJ ].back( ).push_back( symmetry( vVals[ 0 ], vVals[ 1 ] ) );
            }
        }
    }
}

void PartialQuarry::setFlatCoverageValues( )
{
    for( size_t uiJ = 0; uiJ < 2; uiJ++ )
    {
        vvFlatCoverageValues[ uiJ ].clear( );
        if( vvCoverageValues[ uiJ ].size( ) > 0 )
        {
            vvFlatCoverageValues[ uiJ ].reserve( vvCoverageValues[ uiJ ][ 0 ].size( ) );
            for( size_t uiI = 0; uiI < vvCoverageValues[ uiJ ][ 0 ].size( ); uiI++ )
            {
                std::array<size_t, 2> vVal;
                for( size_t uiK = 0; uiK < 2; uiK++ )
                {
                    std::vector<size_t> vCollected;
                    vCollected.reserve( vInGroupCoverage[ uiJ ][ uiK ].size( ) );
                    for( size_t uiX : vInGroupCoverage[ uiJ ][ uiK ] )
                        vCollected.push_back( vvCoverageValues[ uiJ ][ uiX ][ uiI ] );

                    vVal[ uiK ] = getFlatValue( vCollected );
                }
                vvFlatCoverageValues[ uiJ ].push_back( getMixedValue( (double)vVal[ 0 ], (double)vVal[ 1 ] ) );
            }
        }
    }
}

void PartialQuarry::setTracks( )
{
    using namespace pybind11::literals;
    size_t uiDividend = this->xSession[ "dividend" ].get<size_t>( );
    for( size_t uiI = 0; uiI < 2; uiI++ )
    {
        vvMinMaxTracks[ uiI ][ 0 ] = std::numeric_limits<int64_t>::max( );
        vvMinMaxTracks[ uiI ][ 1 ] = std::numeric_limits<int64_t>::min( );
        
        for( size_t uiId : vInGroupCoverage[ uiI ][ 2 ] )
            for( size_t uiX = 0; uiX < vAxisCords[ uiI ].size( ); uiX++ )
            {
                auto uiVal = vvCoverageValues[ uiI ][ uiId ][ uiX ];
                vvMinMaxTracks[ uiI ][ 0 ] = std::min( vvMinMaxTracks[ uiI ][ 0 ], (int64_t)uiVal );
                vvMinMaxTracks[ uiI ][ 1 ] = std::max( vvMinMaxTracks[ uiI ][ 1 ], (int64_t)uiVal );
            }
        pybind11::list vChrs;
        pybind11::list vScreenPoss;
        pybind11::list vIndexStarts;
        pybind11::list vIndexEnds;
        pybind11::list vValues;
        pybind11::list vColors;
        pybind11::list vNames;

        size_t uiCnt = 0;

        for( size_t uiId : vInGroupCoverage[ uiI ][ 2 ] )
        {
            pybind11::list vScreenPos;
            pybind11::list vIndexStart;
            pybind11::list vIndexEnd;
            pybind11::list vValue;
            std::string sChr = "";


            for( size_t uiX = 0; uiX < vAxisCords[ uiI ].size( ); uiX++ )
            {
                auto& xCoord = vAxisCords[ uiI ][ uiX ];
                if( sChr != "" && sChr != xCoord.sChromosome )
                {
                    vChrs.append( sChr ); // @todo remove longest common suffix

                    vScreenPoss.append( vScreenPos );
                    vScreenPos = pybind11::list( );

                    vIndexStarts.append( vIndexStart );
                    vIndexStart = pybind11::list( );

                    vIndexEnds.append( vIndexEnd );
                    vIndexEnd = pybind11::list( );

                    vValues.append( vValue );
                    vValue = pybind11::list( );

                    vColors.append( vColorPaletteAnnotation[ uiCnt % vColorPaletteAnnotation.size( ) ] );

                    vNames.append( vActiveCoverage[ uiI ][ uiId ].first );
                }

                if( uiX == 0 )
                {
                    // zero position at start
                    vIndexStart.append( xCoord.uiIndexPos * uiDividend );
                    vIndexEnd.append( xCoord.uiIndexPos * uiDividend );
                    vScreenPos.append( xCoord.uiScreenPos );
                    vValue.append( vvMinMaxTracks[ uiI ][ 0 ] );
                }

                sChr = xCoord.sChromosome;
                auto uiVal = vvCoverageValues[ uiI ][ uiId ][ uiX ];

                // front corner
                vIndexStart.append( xCoord.uiIndexPos * uiDividend );
                vIndexEnd.append( ( xCoord.uiIndexPos + xCoord.uiSize ) * uiDividend );
                vScreenPos.append( xCoord.uiScreenPos );
                vValue.append( uiVal );

                // rear corner
                vIndexStart.append( xCoord.uiIndexPos * uiDividend );
                vIndexEnd.append( ( xCoord.uiIndexPos + xCoord.uiSize ) * uiDividend );
                vScreenPos.append( xCoord.uiScreenPos + xCoord.uiSize );
                vValue.append( uiVal );

                if( uiX + 1 == vAxisCords[ uiI ].size( ) )
                {
                    // zero position at end
                    vIndexStart.append( ( xCoord.uiIndexPos + xCoord.uiSize ) * uiDividend );
                    vIndexEnd.append( ( xCoord.uiIndexPos + xCoord.uiSize ) * uiDividend );
                    vScreenPos.append( xCoord.uiScreenPos + xCoord.uiSize );
                    vValue.append( vvMinMaxTracks[ uiI ][ 0 ] );
                }
            }

            vChrs.append( sChr ); // @todo remove longest common suffix
            vScreenPoss.append( vScreenPos );
            vIndexStarts.append( vIndexStart );
            vIndexEnds.append( vIndexEnd );
            vValues.append( vValue );
            vColors.append( vColorPaletteAnnotation[ uiCnt % vColorPaletteAnnotation.size( ) ] );
            vNames.append( vActiveCoverage[ uiI ][ uiId ].first );

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
}


const pybind11::dict PartialQuarry::getTracks( bool bXAxis )
{
    update( NodeNames::Tracks );
    return xTracksCDS[ bXAxis ? 0 : 1 ];
}

const std::array<int64_t, 2> PartialQuarry::getMinMaxTracks( bool bXAxis )
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
                                                     { "replicates", "cov_row_b" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::CoverageValues,
                  ComputeNode{ .sNodeName = "coverage_values",
                               .fFunc = &PartialQuarry::setCoverageValues,
                               .vIncomingFunctions = { NodeNames::ActiveCoverage, NodeNames::AxisCoords,
                                                       NodeNames::IntersectionType, NodeNames::Symmetry },
                               .vIncomingSession = { { "settings", "filters", "mapping_q", "val_min" },
                                                     { "settings", "filters", "mapping_q", "val_max" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::FlatCoverageValues,
                  ComputeNode{ .sNodeName = "flat_coverage",
                               .fFunc = &PartialQuarry::setFlatCoverageValues,
                               .vIncomingFunctions = { NodeNames::CoverageValues, NodeNames::BetweenGroup },
                               .vIncomingSession = { },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::Tracks,
                  ComputeNode{ .sNodeName = "coverage_tracks",
                               .fFunc = &PartialQuarry::setTracks,
                               .vIncomingFunctions = { NodeNames::CoverageValues },
                               .vIncomingSession = { },
                               .uiLastUpdated = uiCurrTime } );
}


} // namespace cm