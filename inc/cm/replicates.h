#include "cm/partial_quarry.h"
#include <cmath>

#pragma once

namespace cm
{

bool PartialQuarry::setActiveReplicates( )
{
    auto& rList = this->xSession[ "replicates" ][ "list" ];
    vActiveReplicates.clear( );
    vActiveReplicates.reserve( rList.size( ) );

    std::set<std::string> xHave{ };
    std::array<std::set<std::string>, 2> vGroups;
    for( size_t uiI = 0; uiI < 2; uiI++ )
        for( auto& xRepl : this->xSession[ "replicates" ][ uiI == 0 ? "in_group_a" : "in_group_b" ] )
        {
            CANCEL_RETURN;
            xHave.insert( xRepl.get<std::string>( ) );
            vGroups[ uiI ].insert( xRepl.get<std::string>( ) );
        }

    for( const std::string& rS : xHave )
        vActiveReplicates.push_back( rS );

    vInGroup[ 0 ].clear( );
    vInGroup[ 0 ].reserve( vActiveReplicates.size( ) );
    vInGroup[ 1 ].clear( );
    vInGroup[ 1 ].reserve( vActiveReplicates.size( ) );
    for( size_t uiI = 0; uiI < vActiveReplicates.size( ); uiI++ )
        for( size_t uiJ = 0; uiJ < 2; uiJ++ )
        {
            CANCEL_RETURN;
            if( vGroups[ uiJ ].count( vActiveReplicates[ uiI ] ) > 0 )
                vInGroup[ uiJ ].push_back( uiI );
        }
    END_RETURN;
}


bool PartialQuarry::setIntersectionType( )
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
    END_RETURN;
}

size_t PartialQuarry::symmetry( size_t uiA, size_t uiB )
{
    switch( uiSymmetry )
    {
        case 0: // all
            return uiA;
        case 1: // sym
            return std::min( uiA, uiB );
        case 2: // asym
            if( uiA > uiB )
                return uiA - uiB;
            else
                return 0;
        case 3: // toptobot
        case 4: // bottotop
            return uiA + uiB;
        default:
            throw std::logic_error( "unknown symmetry setting" );
            break;
    }
}

bool PartialQuarry::setBinValues( )
{
    vvBinValues.clear( );
    vvBinValues.reserve( vActiveReplicates.size( ) );
    vActiveReplicatesTotal.clear( );
    vActiveReplicatesTotal.reserve( vActiveReplicates.size( ) );

    size_t uiMinuend = this->xSession[ "settings" ][ "normalization" ][ "min_interactions" ][ "val" ].get<size_t>( );

    for( std::string& sRep : vActiveReplicates )
    {
        vvBinValues.emplace_back( );
        vvBinValues.back( ).reserve( vBinCoords.size( ) );

        bool bHasMapQ = this->xSession[ "replicates" ][ "by_name" ][ sRep ][ "has_map_q" ];
        bool bHasMultiMap = this->xSession[ "replicates" ][ "by_name" ][ sRep ][ "has_multimapping" ];

        size_t uiMapQMin = this->xSession[ "settings" ][ "filters" ][ "mapping_q" ][ "val_min" ].get<size_t>( );
        size_t uiMapQMax = this->xSession[ "settings" ][ "filters" ][ "mapping_q" ][ "val_max" ].get<size_t>( );

        bool bIncomplAlignment = this->xSession[ "settings" ][ "filters" ][ "incomplete_alignments" ].get<bool>( );

        if( !( uiMapQMin == 0 && bIncomplAlignment ) )
            ++uiMapQMin;
        ++uiMapQMax;

        uiMapQMin = 255 - uiMapQMin;
        uiMapQMax = 255 - uiMapQMax;

        for( std::array<BinCoord, 2>& vCoords : vBinCoords )
        {
            CANCEL_RETURN;
            std::array<size_t, 2> vVals;
            for( size_t uiI = 0; uiI < 2; uiI++ )
            {
                if( vCoords[ uiI ].sChromosomeX != "" )
                {
                    size_t iDataSetId = this->xSession[ "replicates" ][ "by_name" ][ sRep ][ "ids" ]
                                                      [ vCoords[ uiI ].sChromosomeX ][ vCoords[ uiI ].sChromosomeY ]
                                                          .get<size_t>( );

                    if( bHasMapQ && bHasMultiMap )
                        vVals[ uiI ] = xIndices.getIndex<3, 2>( )->count(
                            iDataSetId,
                            { vCoords[ uiI ].uiIndexY, vCoords[ uiI ].uiIndexX, uiMapQMax },
                            { vCoords[ uiI ].uiIndexY + vCoords[ uiI ].uiIndexH,
                              vCoords[ uiI ].uiIndexX + vCoords[ uiI ].uiIndexW, uiMapQMin },
                            xIntersect,
                            0 );
                    else if( !bHasMapQ && bHasMultiMap )
                        vVals[ uiI ] =
                            xIndices.getIndex<2, 2>( )->count( iDataSetId,
                                                               { vCoords[ uiI ].uiIndexY, vCoords[ uiI ].uiIndexX },
                                                               { vCoords[ uiI ].uiIndexY + vCoords[ uiI ].uiIndexH,
                                                                 vCoords[ uiI ].uiIndexX + vCoords[ uiI ].uiIndexW },
                                                               xIntersect,
                                                               0 );
                    else if( bHasMapQ && !bHasMultiMap )
                        vVals[ uiI ] = xIndices.getIndex<3, 0>( )->count(
                            iDataSetId,
                            { vCoords[ uiI ].uiIndexY, vCoords[ uiI ].uiIndexX, uiMapQMax },
                            { vCoords[ uiI ].uiIndexY + vCoords[ uiI ].uiIndexH,
                              vCoords[ uiI ].uiIndexX + vCoords[ uiI ].uiIndexW, uiMapQMin },
                            xIntersect,
                            0 );
                    else // if(!bHasMapQ && !bHasMultiMap)
                        vVals[ uiI ] =
                            xIndices.getIndex<2, 0>( )->count( iDataSetId,
                                                               { vCoords[ uiI ].uiIndexY, vCoords[ uiI ].uiIndexX },
                                                               { vCoords[ uiI ].uiIndexY + vCoords[ uiI ].uiIndexH,
                                                                 vCoords[ uiI ].uiIndexX + vCoords[ uiI ].uiIndexW },
                                                               xIntersect,
                                                               0 );
                }
                else
                    vVals[ uiI ] = 0;

                if( vVals[ uiI ] > uiMinuend )
                    vVals[ uiI ] -= uiMinuend;
                else
                    vVals[ uiI ] = 0;
            }

            vvBinValues.back( ).push_back( symmetry( vVals[ 0 ], vVals[ 1 ] ) );

            size_t uiTot = this->xSession[ "replicates" ][ "by_name" ][ sRep ][ "total_reads" ].get<size_t>( );
            vActiveReplicatesTotal.push_back( symmetry(uiTot, uiTot) );
        }
    }
    END_RETURN;
}

bool PartialQuarry::setDecayValues( )
{
    vvDecayValues.clear( );
    vvDecayValues.reserve( vActiveReplicates.size( ) );

    size_t uiMinuend = this->xSession[ "settings" ][ "normalization" ][ "min_interactions" ][ "val" ].get<size_t>( );

    for( std::string& sRep : vActiveReplicates )
    {
        vvDecayValues.emplace_back( );
        vvDecayValues.back( ).reserve( vDistDepDecCoords.size( ) );

        bool bHasMapQ = this->xSession[ "replicates" ][ "by_name" ][ sRep ][ "has_map_q" ];
        bool bHasMultiMap = this->xSession[ "replicates" ][ "by_name" ][ sRep ][ "has_multimapping" ];

        size_t uiMapQMin = this->xSession[ "settings" ][ "filters" ][ "mapping_q" ][ "val_min" ].get<size_t>( );
        size_t uiMapQMax = this->xSession[ "settings" ][ "filters" ][ "mapping_q" ][ "val_max" ].get<size_t>( );

        bool bIncomplAlignment = this->xSession[ "settings" ][ "filters" ][ "incomplete_alignments" ].get<bool>( );

        if( !( uiMapQMin == 0 && bIncomplAlignment ) )
            ++uiMapQMin;
        ++uiMapQMax;

        uiMapQMin = 255 - uiMapQMin;
        uiMapQMax = 255 - uiMapQMax;

        for( std::array<DecayCoord, 2>& vCoords : vDistDepDecCoords )
        {
            CANCEL_RETURN;
            std::array<size_t, 2> vVals;
            for( size_t uiI = 0; uiI < 2; uiI++ )
            {
                if( vCoords[ uiI ].sChromosome != "" )
                {
                    size_t iDataSetId = this->xSession[ "replicates" ][ "by_name" ][ sRep ][ "ids" ]
                                                      [ vCoords[ uiI ].sChromosome ][ "dist_dep_dec" ]
                                                          .get<size_t>( );

                    if( bHasMapQ && bHasMultiMap )
                        vVals[ uiI ] = xIndices.getIndex<2, 1>( )->count(
                            iDataSetId,
                            { vCoords[ uiI ].uiFrom, uiMapQMax },
                            { vCoords[ uiI ].uiTo, uiMapQMin },
                            xIntersect,
                            0 );
                    else if( !bHasMapQ && bHasMultiMap )
                        vVals[ uiI ] =
                            xIndices.getIndex<1, 1>( )->count( iDataSetId,
                                                               { vCoords[ uiI ].uiFrom },
                                                               { vCoords[ uiI ].uiTo },
                                                               xIntersect,
                                                               0 );
                    else if( bHasMapQ && !bHasMultiMap )
                        vVals[ uiI ] = xIndices.getIndex<2, 0>( )->count(
                            iDataSetId,
                            { vCoords[ uiI ].uiFrom, uiMapQMax },
                            { vCoords[ uiI ].uiTo, uiMapQMin },
                            xIntersect,
                            0 );
                    else // if(!bHasMapQ && !bHasMultiMap)
                        vVals[ uiI ] =
                            xIndices.getIndex<1, 0>( )->count( iDataSetId,
                                                               { vCoords[ uiI ].uiFrom },
                                                               { vCoords[ uiI ].uiTo },
                                                               xIntersect,
                                                               0 );
                }
                else
                    vVals[ uiI ] = 0;

                if( vVals[ uiI ] > uiMinuend )
                    vVals[ uiI ] -= uiMinuend;
                else
                    vVals[ uiI ] = 0;
            }

            vvDecayValues.back( ).push_back( symmetry( vVals[ 0 ], vVals[ 1 ] ) );
        }
    }
    END_RETURN;
}

bool PartialQuarry::setInGroup( )
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
    else if( sInGroupSetting == "mean" )
        iInGroupSetting = 4;
    else
        throw std::logic_error( "invalid value for in_group" );
    END_RETURN;
}

bool PartialQuarry::setBetweenGroup( )
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
    else if( sBetwGroupSetting == "div" )
        iBetweenGroupSetting = 7;
    else
        throw std::logic_error( "invalid value for between_group" );
    END_RETURN;
}


size_t PartialQuarry::getFlatValue( std::vector<size_t> vCollected )
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
            case 4:
                std::sort( vCollected.begin( ), vCollected.end( ) );
                uiVal = vCollected[ vCollected.size( ) / 2 ];
                break;
            default:
                throw std::logic_error( "invalid value for in_group" );
        }
    return uiVal;
}

double PartialQuarry::getMixedValue( double uiA, double uiB )
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
        case 7:
            return uiA / uiB;
        default:
            throw std::logic_error( "invalid value for between_group" );
    }
}

bool PartialQuarry::setFlatValues( )
{
    vvFlatValues.clear( );

    if( vvBinValues.size( ) > 0 )
    {
        vvFlatValues.reserve( vvBinValues[ 0 ].size( ) );

        for( size_t uiI = 0; uiI < vBinCoords.size( ); uiI++ )
        {
            vvFlatValues.push_back( { 0, 0 } );
            for( size_t uiJ = 0; uiJ < 2; uiJ++ )
            {
                std::vector<size_t> vCollected;
                vCollected.reserve( vInGroup[ uiJ ].size( ) );
                for( size_t uiX : vInGroup[ uiJ ] )
                {
                    CANCEL_RETURN;
                    vCollected.push_back( vvBinValues[ uiX ][ uiI ] );
                }

                vvFlatValues.back( )[ uiJ ] = getFlatValue( vCollected );
            }
        }

        for( size_t uiJ = 0; uiJ < 2; uiJ++ )
        {
            std::vector<size_t> vCollected;
            vCollected.reserve( vInGroup[ uiJ ].size( ) );
            for( size_t uiX : vInGroup[ uiJ ] )
            {
                CANCEL_RETURN;
                vCollected.push_back( vActiveReplicatesTotal[ uiX ] );
            }

            vvFlatTotal[ uiJ ] = getFlatValue( vCollected );
        }
    }
    END_RETURN;
}

bool PartialQuarry::setFlatDecay( )
{
    vvFlatDecay.clear( );

    if( vvDecayValues.size( ) > 0 )
    {
        vvFlatDecay.reserve( vvDecayValues[ 0 ].size( ) );

        for( size_t uiI = 0; uiI < vDistDepDecCoords.size( ); uiI++ )
        {
            vvFlatDecay.push_back( { 0, 0 } );
            for( size_t uiJ = 0; uiJ < 2; uiJ++ )
            {
                std::vector<size_t> vCollected;
                vCollected.reserve( vInGroup[ uiJ ].size( ) );
                for( size_t uiX : vInGroup[ uiJ ] )
                {
                    CANCEL_RETURN;
                    vCollected.push_back( vvDecayValues[ uiX ][ uiI ] );
                }

                vvFlatDecay.back( )[ uiJ ] = getFlatValue( vCollected );
            }
        }
    }
    END_RETURN;
}

bool PartialQuarry::setDecayCDS( )
{
    using namespace pybind11::literals;
    pybind11::gil_scoped_acquire acquire;

    pybind11::list vChrs;
    pybind11::list vXs;
    pybind11::list vYs;
    pybind11::list vColors;

    size_t uiCnt = 0;
    for( size_t uiJ = 0; uiJ < 2; uiJ++ )
    {
        std::string sChr = "";
        pybind11::list vX;
        pybind11::list vY;
        for( size_t uiI : vSortedDistDepDecCoords[uiJ] )
        {
            std::string sChrCurr = vDistDepDecCoords[uiJ][uiI].sChromosome;
            if(sChr.size() > 0 && sChrCurr != sChr)
            {
                vChrs.append( substringChr(sChr) + (uiJ == 0 ? " A" : " B") );

                vXs.append(vX);
                vX = pybind11::list();

                vYs.append(vY);
                vY = pybind11::list();

                vColors.append( vColorPaletteAnnotation[ uiCnt % vColorPaletteAnnotation.size( ) ]);
            }
            sChr = sChrCurr;

            vX.append(vDistDepDecCoords[uiJ][uiI].uiFrom);
            vY.append(vvFlatDecay[uiI][uiJ]);

            vX.append(vDistDepDecCoords[uiJ][uiI].uiTo);
            vY.append(vvFlatDecay[uiI][uiJ]);

            CANCEL_RETURN;
        }

        vChrs.append( substringChr(sChr) + (uiJ == 0 ? " A" : " B") );
        vXs.append(vX);
        vYs.append(vY);
        vColors.append( vColorPaletteAnnotation[ uiCnt % vColorPaletteAnnotation.size( ) ]);

        ++uiCnt;
    }


    xDistDepDecCDS = pybind11::dict( "color"_a = vColors,
                                  "chr"_a = vChrs,
                                  "xs"_a = vXs,
                                  "ys"_a = vYs );
    END_RETURN;
}


const pybind11::dict PartialQuarry::getDecayCDS( )
{
    update( NodeNames::DecayCDS );
    return xDistDepDecCDS;
}

void PartialQuarry::regReplicates( )
{
    registerNode(
        NodeNames::ActiveReplicates,
        ComputeNode{ .sNodeName = "active_replicates",
                     .fFunc = &PartialQuarry::setActiveReplicates,
                     .vIncomingFunctions = { },
                     .vIncomingSession = { { "replicates", "in_group_a" }, { "replicates", "in_group_b" } } } );

    registerNode( NodeNames::IntersectionType,
                  ComputeNode{ .sNodeName = "intersection_type",
                               .fFunc = &PartialQuarry::setIntersectionType,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "settings", "filters", "ambiguous_mapping" } } } );

    registerNode( NodeNames::BinValues,
                  ComputeNode{ .sNodeName = "bin_values",
                               .fFunc = &PartialQuarry::setBinValues,
                               .vIncomingFunctions = { NodeNames::BinCoords, NodeNames::ActiveReplicates },
                               .vIncomingSession = {} } );

    registerNode( NodeNames::DecayValues,
                  ComputeNode{ .sNodeName = "decay_values",
                               .fFunc = &PartialQuarry::setDecayValues,
                               .vIncomingFunctions = { NodeNames::DecayCoords, NodeNames::ActiveReplicates },
                               .vIncomingSession = {} } );

    registerNode( NodeNames::InGroup,
                  ComputeNode{ .sNodeName = "in_group_setting",
                               .fFunc = &PartialQuarry::setInGroup,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "settings", "replicates", "in_group" } } } );

    registerNode( NodeNames::BetweenGroup,
                  ComputeNode{ .sNodeName = "between_group_setting",
                               .fFunc = &PartialQuarry::setBetweenGroup,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "settings", "replicates", "between_group" } } } );

    registerNode( NodeNames::FlatValues,
                  ComputeNode{ .sNodeName = "flat_bins",
                               .fFunc = &PartialQuarry::setFlatValues,
                               .vIncomingFunctions = { NodeNames::BinValues },
                               .vIncomingSession = {} } );

    registerNode( NodeNames::FlatDecay,
                  ComputeNode{ .sNodeName = "flat_decay",
                               .fFunc = &PartialQuarry::setFlatDecay,
                               .vIncomingFunctions = { NodeNames::DecayValues },
                               .vIncomingSession = {} } );

    registerNode( NodeNames::DecayCDS,
                  ComputeNode{ .sNodeName = "decay_cds",
                               .fFunc = &PartialQuarry::setDecayCDS,
                               .vIncomingFunctions = { NodeNames::FlatDecay, NodeNames::AnnotationColors },
                               .vIncomingSession = {} } );
}

} // namespace cm