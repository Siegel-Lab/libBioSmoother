#include "cm/partial_quarry.h"
#include <cmath>

#pragma once

namespace cm
{

void PartialQuarry::setActiveCoverage( )
{
    size_t uiSize = this->xSession[ "coverage" ][ "list" ].size( ) + this->xSession[ "replicates" ][ "list" ].size( );
    vActiveCoverage[ 0 ].reserve( uiSize );
    vActiveCoverage[ 0 ].clear( );
    vActiveCoverage[ 1 ].reserve( uiSize );
    vActiveCoverage[ 1 ].clear( );

    for( bool bB : { true, false } )
        for( auto& xName : this->xSession[ bB ? "coverage" : "replicates" ][ "list" ] )
        {
            std::string sName = xName.get<std::string>( );
            if( this->xSession[ bB ? "coverage" : "replicates" ][ "coverage" ][ sName ][ "in_column" ].get<bool>( ) ||
                this->xSession[ bB ? "coverage" : "replicates" ][ "coverage" ][ sName ][ "cov_column_a" ]
                    .get<bool>( ) ||
                this->xSession[ bB ? "coverage" : "replicates" ][ "coverage" ][ sName ][ "cov_column_b" ].get<bool>( ) )
                vActiveCoverage[ 0 ].push_back( std::make_pair( sName, bB ) );

            if( this->xSession[ bB ? "coverage" : "replicates" ][ "coverage" ][ sName ][ "in_row" ].get<bool>( ) ||
                this->xSession[ bB ? "coverage" : "replicates" ][ "coverage" ][ sName ][ "cov_row_a" ].get<bool>( ) ||
                this->xSession[ bB ? "coverage" : "replicates" ][ "coverage" ][ sName ][ "cov_row_b" ].get<bool>( ) )
                vActiveCoverage[ 1 ].push_back( std::make_pair( sName, bB ) );
        }
}

void PartialQuarry::setCoverageValues( )
{
    for( size_t uiJ = 0; uiJ < 2; uiJ++ )
    {
        vvCoverageValues[ uiJ ].reserve( vActiveCoverage[ uiJ ].size( ) );
        vvCoverageValues[ uiJ ].clear( );
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
        vInGroup[ uiJ ].reserve( vvCoverageValues[ uiJ ].size( ) );

        for( size_t uiI = 0; uiI < vvCoverageValues[ uiJ ].size( ); uiI++ )
        {
            for( size_t uiK = 0; uiK < 2; uiK++ )
                if( this->xSession[ vActiveCoverage[ uiJ ][ uiI ].second ? "coverage" : "replicates" ][ "coverage" ]
                                  [ vActiveCoverage[ uiJ ][ uiI ].first ]
                                  [ uiJ == 0 ? ( uiK == 1 ? "cov_column_a" : "cov_column_b" )
                                             : ( uiK == 1 ? "cov_row_a" : "cov_row_b" ) ]
                                      .get<bool>( ) )
                    vInGroup[ uiK ].push_back( uiI );
        }

        vvFlatCoverageValues[ uiJ ].reserve( vvCoverageValues[ uiJ ][ 0 ].size( ) );
        vvFlatCoverageValues[ uiJ ].clear( );
        for( size_t uiI = 0; uiI < vvFlatCoverageValues[ uiJ ].size( ); uiI++ )
            for( size_t uiJ = 0; uiJ < 2; uiJ++ )
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


void PartialQuarry::regCoverage( )
{
    registerNode( NodeNames::ActiveCoverage,
                  ComputeNode{ .sNodeName = "active_coverage",
                               .fFunc = &PartialQuarry::setActiveCoverage,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "coverage", "coverage" }, { "replicates", "coverage" } },
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
                               .vIncomingSession = { { "coverage", "coverage" }, { "replicates", "coverage" } },
                               .uiLastUpdated = uiCurrTime } );
}


} // namespace cm