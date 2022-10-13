#include "cm/computation.h"
#include <cmath>

#pragma once

namespace cm
{

void ContactMapping::setActiveCoverage( )
{
    for( bool bB : { true, false } )
        for( auto& xName : this->xSession[ bB ? "coverage" : "replicates" ][ "list" ] )
        {
            std::string sName = xName.get<std::string>( );
            if( this->xSession[ bB ? "coverage" : "replicates" ][ "by_name" ][ sName ][ "in_column" ].get<bool>( ) ||
                this->xSession[ bB ? "coverage" : "replicates" ][ "by_name" ][ sName ][ "norm_column" ].get<bool>( ) )
                vActiveCoverage[ 0 ].push_back( std::make_pair( sName, bB ) );

            if( this->xSession[ bB ? "coverage" : "replicates" ][ "by_name" ][ sName ][ "in_row" ].get<bool>( ) ||
                this->xSession[ bB ? "coverage" : "replicates" ][ "by_name" ][ sName ][ "norm_row" ].get<bool>( ) )
                vActiveCoverage[ 1 ].push_back( std::make_pair( sName, bB ) );
        }
}

void ContactMapping::setCoverageValues( )
{
    for( size_t uiJ = 0; uiJ < 2; uiJ++ )
        for( auto& sRep : vActiveCoverage[ uiJ ] )
        {
            vvCoverageValues[ uiJ ].emplace_back( );
            vvCoverageValues[ uiJ ].back( ).reserve( vAxisCords[ uiJ ].size( ) );
            auto& xRep = this->xSession[ sRep.second ? "coverage" : "replicates" ][ "by_name" ][ sRep.first ];

            bool bHasMapQ = xRep[ "has_map_q" ];
            bool bHasMultiMap = xRep[ "has_multimapping" ];

            size_t uiMapQMin = 255 - this->xRenderSettings[ "filters" ][ "mapping_q" ][ "val_min" ].get<size_t>( );
            size_t uiMapQMax = 255 - this->xRenderSettings[ "filters" ][ "mapping_q" ][ "val_max" ].get<size_t>( );

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
                                vVals[ uiI ] = xIndices.getIndex<2, 1>( ).count(
                                    iDataSetId,
                                    { xCoords.uiIndexPos, uiMapQMax },
                                    { xCoords.uiIndexPos + xCoords.uiSize, uiMapQMin },
                                    xIntersect,
                                    0 );
                            else if( !bHasMapQ && bHasMultiMap )
                                vVals[ uiI ] =
                                    xIndices.getIndex<1, 1>( ).count( iDataSetId,
                                                                      { xCoords.uiIndexPos },
                                                                      { xCoords.uiIndexPos + xCoords.uiSize },
                                                                      xIntersect,
                                                                      0 );
                            else if( bHasMapQ && !bHasMultiMap )
                                vVals[ uiI ] = xIndices.getIndex<2, 0>( ).count(
                                    iDataSetId,
                                    { xCoords.uiIndexPos, uiMapQMax },
                                    { xCoords.uiIndexPos + xCoords.uiSize, uiMapQMin },
                                    xIntersect,
                                    0 );
                            else // if(!bHasMapQ && !bHasMultiMap)
                                vVals[ uiI ] =
                                    xIndices.getIndex<1, 0>( ).count( iDataSetId,
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

void ContactMapping::setFlatCoverageValues( )
{
    for( size_t uiJ = 0; uiJ < 2; uiJ++ )
    {
        std::vector<size_t> vInGroup;
        vInGroup.reserve( vvCoverageValues[ uiJ ].size( ) );

        for( size_t uiI = 0; uiI < vvCoverageValues[ uiJ ].size( ); uiI++ )
        {
            if( this->xSession[ vActiveCoverage[ uiJ ][ uiI ].second ? "coverage" : "replicates" ][ "by_name" ]
                              [ vActiveCoverage[ uiJ ][ uiI ].first ][ uiJ == 0 ? "norm_column" : "norm_row" ]
                                  .get<bool>( ) )
                vInGroup.push_back( uiI );
        }

        vvFlatCoverageValues[ uiJ ].reserve( vvCoverageValues[ uiJ ][ 0 ].size( ) );
        for( size_t uiI = 0; uiI < vvFlatCoverageValues[ uiJ ].size( ); uiI++ )
            for( size_t uiJ = 0; uiJ < 2; uiJ++ )
            {
                std::vector<size_t> vCollected;
                vCollected.reserve( vInGroup.size( ) );
                for( size_t uiX : vInGroup )
                    vCollected.push_back( vvCoverageValues[ uiJ ][ uiX ][ uiI ] );

                vvFlatCoverageValues[ uiJ ].push_back( getFlatValue( vCollected ) );
            }
    }
}


} // namespace cm