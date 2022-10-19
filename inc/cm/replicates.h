#include "cm/computation.h"
#include <cmath>

#pragma once

namespace cm
{

void ContactMapping::setActiveReplicates( )
{
    auto& rList = this->xSession[ "replicates" ][ "list" ];
    vActiveReplicates.resize( rList.size( ) );
    vActiveReplicates.clear( );

    for( auto& xRepl : rList )
    {
        std::string sRepl = xRepl.get<std::string>( );
        if( this->xSession[ "replicates" ][ "in_group" ][ sRepl ][ "in_group_a" ].get<bool>( ) ||
            this->xSession[ "replicates" ][ "in_group" ][ sRepl ][ "in_group_b" ].get<bool>( ) )
            vActiveReplicates.push_back( sRepl );
    }
}


void ContactMapping::setIntersectionType( )
{
    std::string sRenderSetting = this->xSession[ "settings" ][ "filters" ][ "ambiguous_mapping" ].get<std::string>( );

    if( sRenderSetting == "enclosed" )
        xIntersect = sps::IntersectionType::enclosed;
    else if( sRenderSetting == "encloses" )
        xIntersect = sps::IntersectionType::encloses;
    else if( sRenderSetting == "overlaps" )
        xIntersect = sps::IntersectionType::overlaps;
    else if( sRenderSetting == "first" )
        xIntersect = sps::IntersectionType::first;
    else if( sRenderSetting == "last" )
        xIntersect = sps::IntersectionType::last;
    else if( sRenderSetting == "points_only" )
        xIntersect = sps::IntersectionType::points_only;
    else
        throw std::logic_error( "unknown ambiguous_mapping value" );
}

size_t ContactMapping::symmetry( size_t uiA, size_t uiB )
{
    switch( uiSymmetry )
    {
        case 0:
            return uiA;
        case 1:
            return std::min( uiA, uiB );
        case 2:
            return (size_t)std::abs( (int64_t)uiA - (int64_t)uiB );
        case 3:
        case 4:
            return uiA + uiB;
        default:
            throw std::logic_error( "unknown symmetry setting" );
            break;
    }
}

void ContactMapping::setBinValues( )
{
    vvBinValues.reserve( vActiveReplicates.size( ) );
    vvBinValues.clear( );

    for( std::string& sRep : vActiveReplicates )
    {
        vvBinValues.emplace_back( );
        vvBinValues.back( ).reserve( vBinCoords.size( ) );
        auto& xRep = this->xSession[ "replicates" ][ "by_name" ][ sRep ];

        bool bHasMapQ = xRep[ "has_map_q" ];
        bool bHasMultiMap = xRep[ "has_multimapping" ];

        size_t uiMapQMin = 255 - this->xSession[ "settings" ][ "filters" ][ "mapping_q" ][ "val_min" ].get<size_t>( );
        size_t uiMapQMax = 255 - this->xSession[ "settings" ][ "filters" ][ "mapping_q" ][ "val_max" ].get<size_t>( );

        for( std::array<BinCoord, 2>& vCoords : vBinCoords )
        {
            std::array<size_t, 2> vVals;
            for( size_t uiI = 0; uiI < 2; uiI++ )
                if( vCoords[ uiI ].sChromosomeX != "" )
                {
                    int64_t iDataSetId =
                        xRep[ "ids" ][ vCoords[ uiI ].sChromosomeX ][ vCoords[ uiI ].sChromosomeY ].get<int64_t>( );

                    if( iDataSetId != -1 )
                    {
                        if( bHasMapQ && bHasMultiMap )
                            vVals[ uiI ] = xIndices.getIndex<3, 2>( ).count(
                                iDataSetId,
                                { vCoords[ uiI ].uiIndexY, vCoords[ uiI ].uiIndexX, uiMapQMax },
                                { vCoords[ uiI ].uiIndexY + vCoords[ uiI ].uiH,
                                  vCoords[ uiI ].uiIndexX + vCoords[ uiI ].uiW, uiMapQMin },
                                xIntersect,
                                0 );
                        else if( !bHasMapQ && bHasMultiMap )
                            vVals[ uiI ] =
                                xIndices.getIndex<2, 2>( ).count( iDataSetId,
                                                                  { vCoords[ uiI ].uiIndexY, vCoords[ uiI ].uiIndexX },
                                                                  { vCoords[ uiI ].uiIndexY + vCoords[ uiI ].uiH,
                                                                    vCoords[ uiI ].uiIndexX + vCoords[ uiI ].uiW },
                                                                  xIntersect,
                                                                  0 );
                        else if( bHasMapQ && !bHasMultiMap )
                            vVals[ uiI ] = xIndices.getIndex<3, 0>( ).count(
                                iDataSetId,
                                { vCoords[ uiI ].uiIndexY, vCoords[ uiI ].uiIndexX, uiMapQMax },
                                { vCoords[ uiI ].uiIndexY + vCoords[ uiI ].uiH,
                                  vCoords[ uiI ].uiIndexX + vCoords[ uiI ].uiW, uiMapQMin },
                                xIntersect,
                                0 );
                        else // if(!bHasMapQ && !bHasMultiMap)
                            vVals[ uiI ] =
                                xIndices.getIndex<2, 0>( ).count( iDataSetId,
                                                                  { vCoords[ uiI ].uiIndexY, vCoords[ uiI ].uiIndexX },
                                                                  { vCoords[ uiI ].uiIndexY + vCoords[ uiI ].uiH,
                                                                    vCoords[ uiI ].uiIndexX + vCoords[ uiI ].uiW },
                                                                  xIntersect,
                                                                  0 );
                    }
                    else
                        vVals[ uiI ] = 0;
                }

            vvBinValues.back( ).push_back( symmetry( vVals[ 0 ], vVals[ 1 ] ) );
        }
    }
}

void ContactMapping::setInGroup( )
{
    std::string sInGroupSetting = this->xSession[ "settings" ][ "replicates" ][ "in_group" ].get<std::string>( );
    if( sInGroupSetting == "min" )
        iInGroupSetting = 0;
    else if( sInGroupSetting == "sum" )
        iInGroupSetting = 1;
    else if( sInGroupSetting == "dif" )
        iInGroupSetting = 2;
    else if( sInGroupSetting == "max" )
        iInGroupSetting = 3;
    else
        throw std::logic_error( "invalid value for in_group" );
}

void ContactMapping::setBetweenGroup( )
{
    std::string sBetwGroupSetting = this->xSession[ "settings" ][ "replicates" ][ "between_group" ].get<std::string>( );
    if( sBetwGroupSetting == "1st" )
        iBetweenGroupSetting = 0;
    else if( sBetwGroupSetting == "2nd" )
        iBetweenGroupSetting = 1;
    else if( sBetwGroupSetting == "sub" )
        iBetweenGroupSetting = 2;
    else if( sBetwGroupSetting == "min" )
        iBetweenGroupSetting = 3;
    else if( sBetwGroupSetting == "max" )
        iBetweenGroupSetting = 4;
    else if( sBetwGroupSetting == "dif" )
        iBetweenGroupSetting = 5;
    else if( sBetwGroupSetting == "sum" )
        iBetweenGroupSetting = 6;
    else
        throw std::logic_error( "invalid value for between_group" );
}


size_t ContactMapping::getFlatValue( std::vector<size_t> vCollected )
{
    size_t uiVal = 0;
    if( iInGroupSetting == 0 && vCollected.size( ) > 0 )
        uiVal = std::numeric_limits<size_t>::max( );

    for( size_t uiC : vCollected )
        switch( iInGroupSetting )
        {
            case 0:
                uiVal = std::min( uiVal, uiC );
                break;
            case 1:
                uiVal += uiC;
                break;
            case 2:
                for( size_t uiC2 : vCollected )
                    uiVal += (size_t)std::abs( (int64_t)uiC - (int64_t)uiC2 );
                break;
            case 3:
                uiVal = std::max( uiVal, uiC );
                break;
            default:
                throw std::logic_error( "invalid value for in_group" );
        }
    return uiVal;
}

double ContactMapping::getMixedValue( double uiA, double uiB )
{
    switch( iBetweenGroupSetting )
    {
        case 0:
            return uiA;
        case 1:
            return uiB;
        case 2:
            return uiA - uiB;
        case 3:
            return std::min( uiA, uiB );
        case 4:
            return std::max( uiA, uiB );
        case 5:
            return std::abs( uiA - uiB );
        case 6:
            return uiA + uiB;
        default:
            throw std::logic_error( "invalid value for between_group" );
    }
}

void ContactMapping::setFlatValues( )
{
    vvFlatValues.clear( );

    if( vvBinValues.size( ) > 0 )
    {
        vvFlatValues.reserve( vvBinValues[ 0 ].size( ) );
        std::array<std::vector<size_t>, 2> vInGroup;
        vInGroup[ 0 ].reserve( vActiveReplicates.size( ) );
        vInGroup[ 1 ].reserve( vActiveReplicates.size( ) );
        for( size_t uiI = 0; uiI < vActiveReplicates.size( ); uiI++ )
        {
            if( this->xSession[ "replicates" ][ "in_group" ][ vActiveReplicates[ uiI ] ][ "in_group_a" ].get<bool>( ) )
                vInGroup[ 0 ].push_back( uiI );
            if( this->xSession[ "replicates" ][ "in_group" ][ vActiveReplicates[ uiI ] ][ "in_group_b" ].get<bool>( ) )
                vInGroup[ 1 ].push_back( uiI );
        }

        for( size_t uiI = 0; uiI < vvFlatValues.size( ); uiI++ )
        {
            vvFlatValues.push_back( { } );
            for( size_t uiJ = 0; uiJ < 2; uiJ++ )
            {
                std::vector<size_t> vCollected;
                vCollected.reserve( vInGroup[ uiJ ].size( ) );
                for( size_t uiX : vInGroup[ uiJ ] )
                    vCollected.push_back( vvBinValues[ uiX ][ uiI ] );

                vvFlatValues.back( )[ uiJ ] = getFlatValue( vCollected );
            }
        }
    }
}

void ContactMapping::regReplicates( )
{
    registerNode( NodeNames::ActiveReplicates,
                  ComputeNode{ .sNodeName = "active_replicates",
                               .fFunc = &ContactMapping::setActiveReplicates,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "replicates", "in_group" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::IntersectionType,
                  ComputeNode{ .sNodeName = "intersection_type",
                               .fFunc = &ContactMapping::setIntersectionType,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "settings", "filters", "ambiguous_mapping" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::BinValues,
                  ComputeNode{ .sNodeName = "bin_values",
                               .fFunc = &ContactMapping::setBinValues,
                               .vIncomingFunctions = { NodeNames::BinCoords, NodeNames::ActiveReplicates,
                                                       NodeNames::IntersectionType, NodeNames::Symmetry },
                               .vIncomingSession = { { "settings", "filters", "mapping_q", "val_min" },
                                                     { "settings", "filters", "mapping_q", "val_max" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::InGroup,
                  ComputeNode{ .sNodeName = "in_group_setting",
                               .fFunc = &ContactMapping::setInGroup,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "settings", "replicates", "in_group" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::BetweenGroup,
                  ComputeNode{ .sNodeName = "between_group_setting",
                               .fFunc = &ContactMapping::setBetweenGroup,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "settings", "replicates", "between_group" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::FlatValues,
                  ComputeNode{ .sNodeName = "flat_bins",
                               .fFunc = &ContactMapping::setFlatValues,
                               .vIncomingFunctions = { NodeNames::BinValues, NodeNames::InGroup },
                               .vIncomingSession = { { "replicates", "in_group" } },
                               .uiLastUpdated = uiCurrTime } );
}

} // namespace cm