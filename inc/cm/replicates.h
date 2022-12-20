#include "cm/partial_quarry.h"
#include <cmath>

#pragma once

namespace cm
{

bool PartialQuarry::setActiveReplicates( )
{
    auto rList = getValue<json>( { "replicates", "list" } );
    vActiveReplicates.clear( );
    vActiveReplicates.reserve( rList.size( ) );

    std::set<std::string> xHave{ };
    std::array<std::set<std::string>, 2> vGroups;
    for( size_t uiI = 0; uiI < 2; uiI++ )
        for( auto& xRepl : getValue<json>( { "replicates", uiI == 0 ? "in_group_a" : "in_group_b" } ) )
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
    std::string sRenderSetting = getValue<std::string>( { "settings", "filters", "ambiguous_mapping" } );

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

template <typename v_t> v_t PartialQuarry::symmetry( v_t uiA, v_t uiB )
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

    size_t uiMinuend = getValue<size_t>( { "settings", "normalization", "min_interactions", "val" } );

    for( std::string& sRep : vActiveReplicates )
    {
        vvBinValues.emplace_back( );
        vvBinValues.back( ).reserve( vBinCoords.size( ) );

        bool bHasMapQ = getValue<bool>( { "replicates", "by_name", sRep, "has_map_q" } );
        bool bHasMultiMap = getValue<bool>( { "replicates", "by_name", sRep, "has_multimapping" } );

        size_t uiMapQMin = getValue<size_t>( { "settings", "filters", "mapping_q", "val_min" } );
        size_t uiMapQMax = getValue<size_t>( { "settings", "filters", "mapping_q", "val_max" } );

        bool bIncomplAlignment = getValue<bool>( { "settings", "filters", "incomplete_alignments" } );

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
                    size_t iDataSetId =
                        getValue<size_t>( { "replicates", "by_name", sRep, "ids", vCoords[ uiI ].sChromosomeX,
                                            vCoords[ uiI ].sChromosomeY } );

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

            size_t uiTot = getValue<size_t>( { "replicates", "by_name", sRep, "total_reads" } );
            vActiveReplicatesTotal.push_back( symmetry( uiTot, uiTot ) );
        }
    }
    END_RETURN;
}

bool PartialQuarry::setDecayValues( )
{
    vvDecayValues.clear( );
    vvDecayValues.reserve( vActiveReplicates.size( ) );

    size_t uiMinuend = getValue<size_t>( { "settings", "normalization", "min_interactions", "val" } );
    size_t uiSamplesMin = getValue<size_t>( { "settings", "normalization", "ddd_samples", "val_min" } );
    size_t uiSamplesMax = getValue<size_t>( { "settings", "normalization", "ddd_samples", "val_max" } );
    double fQuantExcl = getValue<double>( { "settings", "normalization", "ddd_quantile", "val" } );

    for( std::string& sRep : vActiveReplicates )
    {
        vvDecayValues.emplace_back( );
        vvDecayValues.back( ).reserve( vDistDepDecCoords.size( ) );

        bool bHasMapQ = getValue<bool>( { "replicates", "by_name", sRep, "has_map_q" } );
        bool bHasMultiMap = getValue<bool>( { "replicates", "by_name", sRep, "has_multimapping" } );

        size_t uiMapQMin = getValue<size_t>( { "settings", "filters", "mapping_q", "val_min" } );
        size_t uiMapQMax = getValue<size_t>( { "settings", "filters", "mapping_q", "val_max" } );

        bool bIncomplAlignment = getValue<bool>( { "settings", "filters", "incomplete_alignments" } );

        if( !( uiMapQMin == 0 && bIncomplAlignment ) )
            ++uiMapQMin;
        ++uiMapQMax;

        uiMapQMin = 255 - uiMapQMin;
        uiMapQMax = 255 - uiMapQMax;

        for( std::array<DecayCoord, 2>& vCoords : vDistDepDecCoords )
        {
            CANCEL_RETURN;
            std::array<double, 2> vVals;
            for( size_t uiI = 0; uiI < 2; uiI++ )
            {
                if( vCoords[ uiI ].sChromosomeX != "" && vCoords[ uiI ].sChromosomeY != "" )
                {
                    size_t iDataSetId =
                        getValue<size_t>( { "replicates", "by_name", sRep, "ids", vCoords[ uiI ].sChromosomeX,
                                            vCoords[ uiI ].sChromosomeY } );
                    int64_t iChrX = getValue<int64_t>( { "contigs", "lengths", vCoords[ uiI ].sChromosomeX } );
                    int64_t iChrY = getValue<int64_t>( { "contigs", "lengths", vCoords[ uiI ].sChromosomeY } );

                    /* The two contigs span a rectange of size cW x cH.
                     * The bin has a bottom left corner with a distance of c form the diagonal (positive = to right, neg
                     * = to left)
                     *
                     */

                    int64_t iCornerPos = ( vCoords[ uiI ].iFrom + vCoords[ uiI ].iTo ) / 2;
                    // corner pos within contig rectangle
                    if( iCornerPos >= -iChrY && iCornerPos <= iChrX )
                    {
                        int64_t iBot = std::abs( iCornerPos );
                        int64_t iChrTopRightCornerPos = iChrX - iChrY;
                        int64_t iTop = iChrX + iChrY - std::abs( iChrTopRightCornerPos - iCornerPos );

                        int64_t iMyH = vCoords[ uiI ].iTo - vCoords[ uiI ].iFrom;
                        int64_t iH = std::max( (int64_t)1, ( ( iTop - iMyH ) - iBot ) / (int64_t)uiSamplesMax );

                        if( ( iTop - iMyH ) - iBot >= (int64_t)(uiSamplesMin-1) * iMyH )
                        {
                            std::vector<size_t> vvVals;
                            vvVals.reserve( uiSamplesMax );
                            for( int64_t iMyBot = iBot; iMyBot <= iTop - iMyH; iMyBot += iH )
                            {
                                int64_t iMyTop = iMyBot + iMyH;

                                assert( iMyBot >= iCornerPos );
                                size_t uiYs = (size_t)( iMyBot + iCornerPos ) / 2;
                                size_t uiXs = (size_t)( iMyBot - iCornerPos ) / 2;
                                assert( iMyTop >= iCornerPos );
                                size_t uiYe = std::max( uiYs + 1, (size_t)( iMyTop + iCornerPos ) / 2 );
                                size_t uiXe = std::max( uiXs + 1, (size_t)( iMyTop - iCornerPos ) / 2 );
                                assert( uiYe <= (size_t)iChrY );
                                assert( uiXe <= (size_t)iChrX );

                                if( bHasMapQ && bHasMultiMap )
                                    vvVals.push_back( xIndices.getIndex<3, 2>( )->count( iDataSetId,
                                                                                         { uiXs, uiYs, uiMapQMax },
                                                                                         { uiXe, uiYe, uiMapQMin },
                                                                                         xIntersect,
                                                                                         0 ) );
                                else if( !bHasMapQ && bHasMultiMap )
                                    vvVals.push_back( xIndices.getIndex<2, 2>( )->count(
                                        iDataSetId, { uiXs, uiYs }, { uiXe, uiYe }, xIntersect, 0 ) );
                                else if( bHasMapQ && !bHasMultiMap )
                                    vvVals.push_back( xIndices.getIndex<3, 0>( )->count( iDataSetId,
                                                                                         { uiXs, uiYs, uiMapQMax },
                                                                                         { uiXe, uiYe, uiMapQMin },
                                                                                         xIntersect,
                                                                                         0 ) );
                                else // if(!bHasMapQ && !bHasMultiMap)
                                    vvVals.push_back( xIndices.getIndex<2, 0>( )->count(
                                        iDataSetId, { uiXs, uiYs }, { uiXe, uiYe }, xIntersect, 0 ) );

                                if( vvVals.back( ) > uiMinuend )
                                    vvVals.back( ) -= uiMinuend;
                                else
                                    vvVals.back( ) = 0;
                            }

                            std::sort( vvVals.begin( ), vvVals.end( ) );
                            if( fQuantExcl >= 0.5 )
                                vVals[ uiI ] = (double)vvVals[ vvVals.size( ) / 2 ];
                            else
                            {
                                vVals[ uiI ] = 0;
                                for( size_t uiJ = vvVals.size( ) * fQuantExcl;
                                     uiJ < vvVals.size( ) * ( 1 - fQuantExcl );
                                     uiJ++ )
                                    vVals[ uiI ] += (double)vvVals[ uiJ ];
                                vVals[ uiI ] /= (size_t)( (double)vvVals.size( ) * ( 1.0 - 2.0 * fQuantExcl ) );
                            }
                        }
                        else
                            vVals[ uiI ] = 0;
                    }
                    else
                        vVals[ uiI ] = 0;
                }
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
    std::string sInGroupSetting = getValue<std::string>( { "settings", "replicates", "in_group" } );
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
    std::string sBetwGroupSetting = getValue<std::string>( { "settings", "replicates", "between_group" } );
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


template <typename v_t> v_t PartialQuarry::getFlatValue( std::vector<v_t> vCollected )
{
    v_t uiVal = 0;
    if( iInGroupSetting == 0 && vCollected.size( ) > 0 )
        uiVal = std::numeric_limits<v_t>::max( );

    for( v_t uiC : vCollected )
        switch( iInGroupSetting )
        {
            case 0:
                uiVal = std::min( uiVal, uiC );
                break;
            case 1:
                uiVal += uiC;
                break;
            case 2:
                for( v_t uiC2 : vCollected )
                    uiVal += uiC >= uiC2 ? uiC - uiC2 : uiC2 - uiC;
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
                std::vector<double> vCollected;
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

    size_t uiDividend = getValue<size_t>( { "dividend" } );

    for( size_t uiJ = 0; uiJ < 2; uiJ++ )
    {
        std::map<std::string, size_t> xColors;
        for( size_t uiI = 0; uiI < vDistDepDecCoords.size( ); uiI++ )
        {
            std::string sChrCurr =
                vDistDepDecCoords[ uiI ][ uiJ ].sChromosomeX + " - " + vDistDepDecCoords[ uiI ][ uiJ ].sChromosomeY;
            if( sChrCurr.size( ) > 0 )
            {
                pybind11::list vX;
                pybind11::list vY;
                int64_t iF = vDistDepDecCoords[ uiI ][ uiJ ].iFrom * (int64_t)uiDividend;
                int64_t iT = vDistDepDecCoords[ uiI ][ uiJ ].iTo * (int64_t)uiDividend;
                double fVal = 1000.0 * 1000.0 * vvFlatDecay[ uiI ][ uiJ ] / (double)( ( iT - iF ) * ( iT - iF ) );
                vX.append( iF );
                vY.append( fVal );

                vX.append( iT );
                vY.append( fVal );

                vChrs.append( substringChr( vDistDepDecCoords[ uiI ][ uiJ ].sChromosomeX ) + " - " +
                              substringChr( vDistDepDecCoords[ uiI ][ uiJ ].sChromosomeY ) +
                              ( uiJ == 0 ? ", Group A" : ", Group B" ) );
                vXs.append( vX );
                vYs.append( vY );
                if( xColors.count( sChrCurr ) == 0 )
                    xColors[ sChrCurr ] = xColors.size( ) % vColorPaletteAnnotation.size( );
                vColors.append( vColorPaletteAnnotation[ xColors[ sChrCurr ] ] );
            }

            CANCEL_RETURN;
        }
    }


    xDistDepDecCDS = pybind11::dict( "color"_a = vColors, "chr"_a = vChrs, "xs"_a = vXs, "ys"_a = vYs );
    END_RETURN;
}


const pybind11::dict PartialQuarry::getDecayCDS( )
{
    update( NodeNames::DecayCDS );
    return xDistDepDecCDS;
}

void PartialQuarry::regReplicates( )
{
    registerNode( NodeNames::ActiveReplicates,
                  ComputeNode{ .sNodeName = "active_replicates",
                               .fFunc = &PartialQuarry::setActiveReplicates,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "replicates", "in_group_a" },
                                                     { "replicates", "in_group_b" },
                                                     { "replicates", "list" } },
                               .vSessionsIncomingInPrevious = {} } );

    registerNode( NodeNames::IntersectionType,
                  ComputeNode{ .sNodeName = "intersection_type",
                               .fFunc = &PartialQuarry::setIntersectionType,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "settings", "filters", "ambiguous_mapping" } },
                               .vSessionsIncomingInPrevious = {} } );

    registerNode(
        NodeNames::BinValues,
        ComputeNode{ .sNodeName = "bin_values",
                     .fFunc = &PartialQuarry::setBinValues,
                     .vIncomingFunctions = { NodeNames::BinCoords, NodeNames::ActiveReplicates },
                     .vIncomingSession = { },
                     .vSessionsIncomingInPrevious = { { "settings", "normalization", "min_interactions", "val" },
                                                      { "replicates", "by_name" },
                                                      { "settings", "filters", "mapping_q", "val_min" },
                                                      { "settings", "filters", "mapping_q", "val_max" },
                                                      { "settings", "filters", "incomplete_alignments" } } } );

    registerNode(
        NodeNames::DecayValues,
        ComputeNode{ .sNodeName = "decay_values",
                     .fFunc = &PartialQuarry::setDecayValues,
                     .vIncomingFunctions = { NodeNames::DecayCoords, NodeNames::ActiveReplicates },
                     .vIncomingSession = { { "settings", "normalization", "ddd_samples", "val_min" },
                                           { "settings", "normalization", "ddd_samples", "val_max" },
                                           { "settings", "normalization", "ddd_quantile", "val" } },
                     .vSessionsIncomingInPrevious = { { "settings", "normalization", "min_interactions", "val" },
                                                      { "replicates", "by_name" },
                                                      { "settings", "filters", "mapping_q", "val_min" },
                                                      { "settings", "filters", "mapping_q", "val_max" },
                                                      { "settings", "filters", "incomplete_alignments" },
                                                      { "contigs", "lengths" } } } );

    registerNode( NodeNames::InGroup,
                  ComputeNode{ .sNodeName = "in_group_setting",
                               .fFunc = &PartialQuarry::setInGroup,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "settings", "replicates", "in_group" } },
                               .vSessionsIncomingInPrevious = {} } );

    registerNode( NodeNames::BetweenGroup,
                  ComputeNode{ .sNodeName = "between_group_setting",
                               .fFunc = &PartialQuarry::setBetweenGroup,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "settings", "replicates", "between_group" } },
                               .vSessionsIncomingInPrevious = {} } );

    registerNode( NodeNames::FlatValues,
                  ComputeNode{ .sNodeName = "flat_bins",
                               .fFunc = &PartialQuarry::setFlatValues,
                               .vIncomingFunctions = { NodeNames::BinValues },
                               .vIncomingSession = { },
                               .vSessionsIncomingInPrevious = {} } );

    registerNode( NodeNames::FlatDecay,
                  ComputeNode{ .sNodeName = "flat_decay",
                               .fFunc = &PartialQuarry::setFlatDecay,
                               .vIncomingFunctions = { NodeNames::DecayValues },
                               .vIncomingSession = { },
                               .vSessionsIncomingInPrevious = {} } );

    registerNode( NodeNames::DecayCDS,
                  ComputeNode{ .sNodeName = "decay_cds",
                               .fFunc = &PartialQuarry::setDecayCDS,
                               .vIncomingFunctions = { NodeNames::FlatDecay, NodeNames::AnnotationColors },
                               .vIncomingSession = { },
                               .vSessionsIncomingInPrevious = { { "dividend" } } } );
}

} // namespace cm