#include "cm/partial_quarry.h"
#include <cmath>

#pragma once

namespace cm
{

void PartialQuarry::setActiveCoverage( )
{
    size_t uiSize = this->xSession[ "coverage" ][ "list" ].size( ) + this->xSession[ "replicates" ][ "list" ].size( );
    vActiveCoverage[ 0 ].clear( );
    vActiveCoverage[ 1 ].clear( );
    vActiveCoverage[ 0 ].reserve( uiSize );
    vActiveCoverage[ 1 ].reserve( uiSize );

    std::set<std::pair<std::string, bool>> xA{ };
    std::set<std::pair<std::string, bool>> xB{ };


    for( bool bB : { true, false } )
    {
        std::string sX = bB ? "coverage" : "replicates";

        for( std::string sY : { "in_column", "cov_column_a", "cov_column_b" } )
            for( auto& rV : this->xSession[ sX ][ sY ] )
                xA.insert( std::make_pair( rV.get<std::string>( ), bB ) );

        for( std::string sY : { "in_row", "cov_row_a", "cov_row_b" } )
            for( auto& rV : this->xSession[ sX ][ sY ] )
                xB.insert( std::make_pair( rV.get<std::string>( ), bB ) );
    }

    std::copy( xA.begin( ), xA.end( ), vActiveCoverage[ 0 ].begin( ) );
    std::copy( xB.begin( ), xB.end( ), vActiveCoverage[ 1 ].begin( ) );
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
        std::array<std::vector<size_t>, 2> vInGroup;
        vInGroup[ 0 ].reserve( vvCoverageValues[ uiJ ].size( ) * 3 );
        vInGroup[ 1 ].reserve( vvCoverageValues[ uiJ ].size( ) * 3 );

        for( size_t uiI = 0; uiI < vvCoverageValues[ uiJ ].size( ); uiI++ )
        {
            for( size_t uiK = 0; uiK < 2; uiK++ )
                if( this->xSession[ vActiveCoverage[ uiJ ][ uiI ].second ? "coverage" : "replicates" ]
                                  [ uiJ == 0 ? ( uiK == 1 ? "cov_column_a" : "cov_column_b" )
                                             : ( uiK == 1 ? "cov_row_a" : "cov_row_b" ) ]

                                        // @todo @fixme contains does nothing use a set instead
                                      .contains( vActiveCoverage[ uiJ ][ uiI ].first ) )
                    vInGroup[ uiK ].push_back( uiI );
        }

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
                    vCollected.reserve( vInGroup[ uiK ].size( ) );
                    for( size_t uiX : vInGroup[ uiK ] )
                        vCollected.push_back( vvCoverageValues[ uiJ ][ uiX ][ uiI ] );

                    vVal[ uiK ] = getFlatValue( vCollected );
                }
                vvFlatCoverageValues[ uiJ ].push_back( getMixedValue( (double)vVal[ 0 ], (double)vVal[ 1 ] ) );
            }
        }
    }
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
                               .vIncomingSession = { { "coverage", "cov_column_a" },
                                                     { "coverage", "cov_column_b" },
                                                     { "coverage", "cov_row_a" },
                                                     { "coverage", "cov_row_b" },
                                                     { "replicates", "cov_column_a" },
                                                     { "replicates", "cov_column_b" },
                                                     { "replicates", "cov_row_a" },
                                                     { "replicates", "cov_row_b" } },
                               .uiLastUpdated = uiCurrTime } );
}


} // namespace cm