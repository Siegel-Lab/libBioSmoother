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

bool PartialQuarry::setDatasetIdPerRepl( )
{
    vvDatasetIdsPerReplAndChr.clear( );
    vvDatasetIdsPerReplAndChr.reserve( vActiveReplicates.size( ) );
    vvDatasetFiltersPerIndex.clear( );
    vvDatasetFiltersPerIndex.reserve( vActiveReplicates.size( ) );

    for( const std::string& sRep : vActiveReplicates )
    {
        vvDatasetIdsPerReplAndChr.push_back(
            getValue<size_t>( { "replicates", "by_name", sRep, "first_dataset_id" } ) );
        vvDatasetFiltersPerIndex.push_back( DatasetFilters{
            /*bHasMapQ:*/ !getValue<bool>( { "replicates", "by_name", sRep, "no_map_q" } ),
            /*bHasMultimappers:*/ !getValue<bool>( { "replicates", "by_name", sRep, "no_multi_map" } ),
            /*bHasStrand:*/ !getValue<bool>( { "replicates", "by_name", sRep, "no_strand" } ),
            /*bHasCategory:*/ !getValue<bool>( { "replicates", "by_name", sRep, "no_category" } ),
        } );
    }
    uiContigListSize = getValue<json>( { "contigs", "ploidy_list" } ).size( );
    END_RETURN;
}


bool PartialQuarry::setIntersectionType( )
{
    std::string sRenderSetting = getValue<std::string>( { "settings", "filters", "ambiguous_mapping" } );

    bOnlyMMRs = false;
    if( sRenderSetting == "enclosed" )
        xIntersect = sps::IntersectionType::enclosed;
    else if( sRenderSetting == "enclosed_only" )
    {
        xIntersect = sps::IntersectionType::enclosed;
        bOnlyMMRs = true;
    }
    else if( sRenderSetting == "encloses" )
        xIntersect = sps::IntersectionType::encloses;
    else if( sRenderSetting == "overlaps" )
        xIntersect = sps::IntersectionType::overlaps;
    else if( sRenderSetting == "overlaps_only" )
    {
        xIntersect = sps::IntersectionType::overlaps;
        bOnlyMMRs = true;
    }
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
        case 3: // mirror
            return uiA + uiB;
        default:
            throw std::logic_error( "unknown symmetry setting" );
            break;
    }
}

bool PartialQuarry::setBinValues( )
{
    bool bOneFilledBinAtLeast = false;
    for( size_t uiY = 0; uiY < NUM_COORD_SYSTEMS; uiY++ )
    {
        vvBinValues[ uiY ].clear( );
        vvBinValues[ uiY ].reserve( vActiveReplicates.size( ) );
        vActiveReplicatesTotal.clear( );
        vActiveReplicatesTotal.reserve( vActiveReplicates.size( ) );

        size_t uiMinuend = getValue<size_t>( { "settings", "normalization", "min_interactions", "val" } );

        for( size_t uiRepl = 0; uiRepl < vActiveReplicates.size( ); uiRepl++ )
        {
            vvBinValues[ uiY ].emplace_back( );
            vvBinValues[ uiY ].back( ).reserve( vBinCoordsSampled[ uiY ].size( ) );

            for( const std::array<BinCoord, 2>& vCoords : vBinCoordsSampled[ uiY ] )
            {
                std::array<size_t, 2> vVals;
                CANCEL_RETURN;
                for( size_t uiI = 0; uiI < 2; uiI++ )
                {
                    if( vCoords[ uiI ].uiChromosomeX != std::numeric_limits<size_t>::max( ) )
                    {
                        size_t iDataSetId = getDatasetIdfromReplAndChr(
                            uiRepl,
                            vActiveChromosomes[ 0 ][ vCoords[ uiI ].uiChromosomeX ].uiActualContigId,
                            vActiveChromosomes[ 1 ][ vCoords[ uiI ].uiChromosomeY ].uiActualContigId );
                        vVals[ uiI ] = indexCount( uiRepl, iDataSetId, vCoords[ uiI ].uiIndexY, vCoords[ uiI ].uiIndexX,
                                                   vCoords[ uiI ].uiIndexY + vCoords[ uiI ].uiIndexH,
                                                   vCoords[ uiI ].uiIndexX + vCoords[ uiI ].uiIndexW, uiI == 0 );
                        if( vVals[ uiI ] > 0 )
                            bOneFilledBinAtLeast = true;
                    }
                    else
                        vVals[ uiI ] = { };

                    vVals[ uiI ] = vVals[ uiI ] > uiMinuend ? vVals[ uiI ] - uiMinuend : 0;
                }
                vvBinValues[ uiY ].back( ).push_back( symmetry( vVals[ 0 ], vVals[ 1 ] ) );
            }

            const size_t uiTot =
                getValue<size_t>( { "replicates", "by_name", vActiveReplicates[ uiRepl ], "total_reads" } );
            vActiveReplicatesTotal.push_back( symmetry( uiTot, uiTot ) );
        }
    }

    if( !bOneFilledBinAtLeast && vvBinValues[ 0 ].size( ) > 0 )
        setError( "All bins are empty." );

    END_RETURN;
}

bool PartialQuarry::setDecayValues( )
{
    std::vector<size_t> vvVals;
    for( size_t uiY = 0; uiY < NUM_COORD_SYSTEMS; uiY++ )
    {
        vvDecayValues[ uiY ].clear( );
        vvDecayValues[ uiY ].reserve( vActiveReplicates.size( ) );

        size_t uiMinuend = getValue<size_t>( { "settings", "normalization", "min_interactions", "val" } );
        size_t uiSamplesMin = getValue<size_t>( { "settings", "normalization", "ddd_samples", "val_min" } );
        size_t uiSamplesMax = getValue<size_t>( { "settings", "normalization", "ddd_samples", "val_max" } );
        bool bAllSamples = getValue<bool>( { "settings", "normalization", "ddd_all_samples" } );
        if( bAllSamples )
            uiSamplesMax = 1;
        double fQuantExcl =
            ( 1.0 - getValue<double>( { "settings", "normalization", "ddd_quantile", "val" } ) / 100.0 ) / 2.0;

        for( size_t uiRepl = 0; uiRepl < vActiveReplicates.size( ); uiRepl++ )
        {
            vvDecayValues[ uiY ].emplace_back( );
            vvDecayValues[ uiY ].back( ).reserve( vDistDepDecCoords[ uiY ].size( ) );

            for( std::array<DecayCoord, 2>& vCoords : vDistDepDecCoords[ uiY ] )
            {
                CANCEL_RETURN;
                std::array<double, 2> vVals;
                for( size_t uiI = 0; uiI < 2; uiI++ )
                {
                    if( vCoords[ uiI ].uiChromosomeX != std::numeric_limits<size_t>::max( ) &&
                        vCoords[ uiI ].uiChromosomeY != std::numeric_limits<size_t>::max( ) )
                    {
                        size_t iDataSetId = getDatasetIdfromReplAndChr(
                            uiRepl,
                            vActiveChromosomes[ 0 ][ vCoords[ uiI ].uiChromosomeX ].uiActualContigId,
                            vActiveChromosomes[ 1 ][ vCoords[ uiI ].uiChromosomeY ].uiActualContigId );

                        int64_t iChrX = vActiveChromosomes[ 0 ][ vCoords[ uiI ].uiChromosomeX ].uiLength;
                        int64_t iChrY = vActiveChromosomes[ 1 ][ vCoords[ uiI ].uiChromosomeY ].uiLength;

                        /* The two contigs span a rectange of size cW x cH.
                         * The bin has a bottom left corner with a distance of c form the diagonal (positive = to
                         * right, neg = to left)
                         *
                         */

                        int64_t iCornerPos = ( vCoords[ uiI ].iFrom + vCoords[ uiI ].iTo ) / 2;
                        // corner pos within contig rectangle
                        if( iCornerPos >= -iChrY && iCornerPos <= iChrX )
                        {
                            int64_t iBot = std::abs( iCornerPos );
                            int64_t iChrTopRightCornerPos = iChrX - iChrY;
                            int64_t iTop = iChrX + iChrY - std::abs( iChrTopRightCornerPos - iCornerPos );

                            int64_t iMyH = std::max( (int64_t)1, vCoords[ uiI ].iTo - vCoords[ uiI ].iFrom );
                            int64_t iH = std::max( (int64_t)1, ( ( iTop - iMyH ) - iBot ) / (int64_t)uiSamplesMax );
                            if( bAllSamples )
                                iH = 1;

                            if( iTop - iBot >= (int64_t)( uiSamplesMin - 1 ) * iMyH )
                            {
                                vvVals.clear( );
                                vvVals.reserve( std::min( (int64_t)uiSamplesMax, (int64_t)( iTop - iBot ) ) );
                                if( bAllSamples )
                                    vvVals.reserve( iTop - iBot );
                                for( int64_t iMyBot = iBot; iMyBot < iTop; iMyBot += iH )
                                {
                                    int64_t iMyTop = iMyBot + iMyH;

                                    assert( iMyBot >= iCornerPos );
                                    size_t uiYs = (size_t)( iMyBot + iCornerPos ) / 2;
                                    size_t uiXs = (size_t)( iMyBot - iCornerPos ) / 2;
                                    assert( iMyTop >= iCornerPos );
                                    size_t uiYe = std::max( uiYs + 1, (size_t)( iMyTop + iCornerPos ) / 2 );
                                    size_t uiXe = std::max( uiXs + 1, (size_t)( iMyTop - iCornerPos ) / 2 );
                                    // assert( uiYe <= (size_t)iChrX );
                                    // assert( uiXe <= (size_t)iChrY );
                                    vvVals.push_back(
                                        indexCount( uiRepl, iDataSetId, uiXs, uiYs, uiXe, uiYe, uiI == 0 ) );

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

                vvDecayValues[ uiY ].back( ).push_back( symmetry( vVals[ 0 ], vVals[ 1 ] ) );
            }
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
    else if( sBetwGroupSetting == "div+" )
        iBetweenGroupSetting = 8;
    else
        throw std::logic_error( "invalid value for between_group" );
    END_RETURN;
}


template <typename v_t> v_t PartialQuarry::getFlatValue( std::vector<v_t> vCollected )
{
    v_t uiVal = 0;
    if( iInGroupSetting == 0 && vCollected.size( ) > 0 )
        uiVal = std::numeric_limits<v_t>::max( );

    if( iInGroupSetting == 4 )
    {
        if( vCollected.size( ) == 0 )
            return { };
        std::sort( vCollected.begin( ), vCollected.end( ) );
        return vCollected[ vCollected.size( ) / 2 ];
    }

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
            return uiB > 0 ? uiA / uiB : 0;
        case 8:
            return ( uiA + 0.000001 ) / ( uiB + 0.000001 );
        default:
            throw std::logic_error( "invalid value for between_group" );
    }
}

bool PartialQuarry::setFlatValues( )
{
    for( size_t uiY = 0; uiY < NUM_COORD_SYSTEMS; uiY++ )
    {
        vvFlatValues[ uiY ].clear( );

        if( vvBinValues[ uiY ].size( ) > 0 )
        {
            // take size from fst replicate (sizes are equal anyways)
            vvFlatValues[ uiY ].reserve( vvBinValues[ uiY ][ 0 ].size( ) );

            for( size_t uiI = 0; uiI < vBinCoordsSampled[ uiY ].size( ); uiI++ )
            {
                vvFlatValues[ uiY ].push_back( { 0, 0 } );
                for( size_t uiJ = 0; uiJ < 2; uiJ++ )
                {
                    std::vector<size_t> vCollected;
                    vCollected.reserve( vInGroup[ uiJ ].size( ) );
                    for( size_t uiX : vInGroup[ uiJ ] )
                    {
                        CANCEL_RETURN;
                        vCollected.push_back( vvBinValues[ uiY ][ uiX ][ uiI ] );
                    }

                    vvFlatValues[ uiY ].back( )[ uiJ ] = getFlatValue( vCollected );
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

                vvFlatTotal[ uiY ][ uiJ ] = getFlatValue( vCollected );
            }
        }
    }
    END_RETURN;
}

double PartialQuarry::ploidyCorrect( size_t uiVal, const BinCoord& rxPos )
{
    if( this->vPloidyCounts[ rxPos.uiChromIdX + rxPos.uiChromIdY * uiFullContigListSize ] == 0 )
        return 0;
    return ( (double)uiVal ) /
           ( (double)this->vPloidyCounts[ rxPos.uiChromIdX + rxPos.uiChromIdY * uiFullContigListSize ] );
}

bool PartialQuarry::setPloidyValues( )
{
    for( size_t uiY = 0; uiY < NUM_COORD_SYSTEMS; uiY++ )
    {
        vvPloidyValues[ uiY ].clear( );

        if( vvFlatValues[ uiY ].size( ) > 0 )
        {
            // take size from fst replicate (sizes are equal anyways)
            vvPloidyValues[ uiY ].reserve( vvFlatValues[ uiY ][ 0 ].size( ) );

            for( size_t uiI = 0; uiI < vBinCoordsSampled[ uiY ].size( ); uiI++ )
            {
                vvPloidyValues[ uiY ].push_back( {
                    ploidyCorrect( vvFlatValues[ uiY ][ uiI ][ 0 ], vBinCoordsSampled[ uiY ][ uiI ][ 0 ] ),
                    ploidyCorrect( vvFlatValues[ uiY ][ uiI ][ 1 ], vBinCoordsSampled[ uiY ][ uiI ][ 1 ] ),
                } );
            }
        }
    }
    END_RETURN;
}

bool PartialQuarry::setFlatDecay( )
{
    for( size_t uiY = 0; uiY < NUM_COORD_SYSTEMS; uiY++ )
    {
        vvFlatDecay[ uiY ].clear( );

        if( vvDecayValues[ uiY ].size( ) > 0 )
        {
            vvFlatDecay[ uiY ].reserve( vvDecayValues[ uiY ][ 0 ].size( ) );

            for( size_t uiI = 0; uiI < vDistDepDecCoords[ uiY ].size( ); uiI++ )
            {
                vvFlatDecay[ uiY ].push_back( { 0, 0 } );
                for( size_t uiJ = 0; uiJ < 2; uiJ++ )
                {
                    std::vector<double> vCollected;
                    vCollected.reserve( vInGroup[ uiJ ].size( ) );
                    for( size_t uiX : vInGroup[ uiJ ] )
                    {
                        CANCEL_RETURN;
                        vCollected.push_back( vvDecayValues[ uiY ][ uiX ][ uiI ] );
                    }

                    vvFlatDecay[ uiY ].back( )[ uiJ ] = getFlatValue( vCollected );
                }
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

    const size_t uiMaxChar = getValue<size_t>( { "settings", "interface", "axis_label_max_char", "val" } );

    size_t uiDividend = getValue<size_t>( { "dividend" } );

    for( size_t uiJ = 0; uiJ < 2; uiJ++ )
    {
        std::map<std::pair<size_t, size_t>, size_t> xColors;
        for( size_t uiI = 0; uiI < vDistDepDecCoords[ 0 ].size( ); uiI++ )
        {
            std::pair<size_t, size_t> xuiChrCurr = std::make_pair( vDistDepDecCoords[ 0 ][ uiI ][ uiJ ].uiChromosomeX,
                                                                   vDistDepDecCoords[ 0 ][ uiI ][ uiJ ].uiChromosomeY );
            pybind11::list vX;
            pybind11::list vY;
            int64_t iF = vDistDepDecCoords[ 0 ][ uiI ][ uiJ ].iFrom * (int64_t)uiDividend;
            int64_t iT = vDistDepDecCoords[ 0 ][ uiI ][ uiJ ].iTo * (int64_t)uiDividend;
            double fVal = 1000.0 * 1000.0 * vvFlatDecay[ 0 ][ uiI ][ uiJ ] / (double)( ( iT - iF ) * ( iT - iF ) );
            vX.append( iF );
            vY.append( fVal );

            vX.append( iT );
            vY.append( fVal );

            std::string sChromNameX =
                vActiveChromosomes[ uiJ ][ vDistDepDecCoords[ 0 ][ uiI ][ uiJ ].uiChromosomeX ].sName;
            std::string sChromNameY =
                vActiveChromosomes[ uiJ ][ vDistDepDecCoords[ 0 ][ uiI ][ uiJ ].uiChromosomeY ].sName;
            vChrs.append( substringChr( sChromNameX, uiMaxChar ) + " - " + substringChr( sChromNameY, uiMaxChar ) +
                          ( uiJ == 0 ? ", Group A" : ", Group B" ) );
            vXs.append( vX );
            vYs.append( vY );
            if( xColors.count( xuiChrCurr ) == 0 )
                xColors[ xuiChrCurr ] = xColors.size( ) % vColorPaletteAnnotation.size( );
            vColors.append( vColorPaletteAnnotation[ xColors[ xuiChrCurr ] ] );

            CANCEL_RETURN;
        }
    }


    xDistDepDecCDS = pybind11::dict( "color"_a = vColors, "chr"_a = vChrs, "xs"_a = vXs, "ys"_a = vYs );
    END_RETURN;
}


const pybind11::dict PartialQuarry::getDecayCDS( const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::DecayCDS, fPyPrint );
    return xDistDepDecCDS;
}

void PartialQuarry::regReplicates( )
{
    registerNode(
        NodeNames::ActiveReplicates,
        ComputeNode{ /*.sNodeName =*/"active_replicates_setting",
                     /*.sNodeDesc =*/"Extract the active replicates from the session json.",
                     /*.fFunc=*/&PartialQuarry::setActiveReplicates,
                     /*.vIncomingFunctions =*/{ },
                     /*.vIncomingSession =*/
                     { { "replicates", "in_group_a" }, { "replicates", "in_group_b" }, { "replicates", "list" } },
                     /*.vSessionsIncomingInPrevious =*/{ },
                     /*bHidden =*/true } );

    registerNode( NodeNames::DatasetIdPerRepl,
                  ComputeNode{ /*.sNodeName =*/"dataset_id_per_repl",
                               /*.sNodeDesc =*/"Compute the dataset id for each replicate and active contig.",
                               /*.fFunc=*/&PartialQuarry::setDatasetIdPerRepl,
                               /*.vIncomingFunctions =*/{ NodeNames::ActiveReplicates, NodeNames::ActiveChrom },
                               /*.vIncomingSession =*/{ { "replicates", "by_name" } },
                               /*.vSessionsIncomingInPrevious =*/{ { "contigs", "ploidy_list" } },
                               /*bHidden =*/true } );

    registerNode( NodeNames::IntersectionType,
                  ComputeNode{ /*.sNodeName =*/"intersection_type_setting",
                               /*.sNodeDesc =*/"Extract the intersection type from the session json.",
                               /*.fFunc=*/&PartialQuarry::setIntersectionType,
                               /*.vIncomingFunctions =*/{ },
                               /*.vIncomingSession =*/{ { "settings", "filters", "ambiguous_mapping" } },
                               /*.vSessionsIncomingInPrevious =*/{ },
                               /*bHidden =*/true } );

    registerNode( NodeNames::BinValues,
                  ComputeNode{ /*.sNodeName =*/"bin_values",
                               /*.sNodeDesc =*/"Compute the values of all bins.",
                               /*.fFunc=*/&PartialQuarry::setBinValues,
                               /*.vIncomingFunctions =*/
                               { NodeNames::BinCoords, NodeNames::DatasetIdPerRepl, NodeNames::MappingQuality,
                                 NodeNames::Directionality },
                               /*.vIncomingSession =*/
                               { { "settings", "normalization", "min_interactions", "val" } },
                               /*.vSessionsIncomingInPrevious =*/
                               { { "replicates", "by_name" } },
                               /*bHidden =*/false } );


    registerNode( NodeNames::DecayValues,
                  ComputeNode{ /*.sNodeName =*/"decay_values",
                               /*.sNodeDesc =*/"Compute the coverage of the DDD normalization.",
                               /*.fFunc=*/&PartialQuarry::setDecayValues,
                               /*.vIncomingFunctions =*/
                               { NodeNames::DecayCoords, NodeNames::DatasetIdPerRepl, NodeNames::MappingQuality,
                                 NodeNames::Directionality },
                               /*.vIncomingSession =*/
                               { { "settings", "normalization", "ddd_samples", "val_min" },
                                 { "settings", "normalization", "ddd_samples", "val_max" },
                                 { "settings", "normalization", "ddd_all_samples" },
                                 { "settings", "normalization", "ddd_quantile", "val" },
                                 { "settings", "normalization", "min_interactions", "val" } },
                               /*.vSessionsIncomingInPrevious =*/{ },
                               /*bHidden =*/false } );

    registerNode( NodeNames::InGroup,
                  ComputeNode{ /*.sNodeName =*/"in_group_setting",
                               /*.sNodeDesc =*/"Extract the in_group setting from the session json.",
                               /*.fFunc=*/&PartialQuarry::setInGroup,
                               /*.vIncomingFunctions =*/{ },
                               /*.vIncomingSession =*/{ { "settings", "replicates", "in_group" } },
                               /*.vSessionsIncomingInPrevious =*/{ },
                               /*bHidden =*/true } );

    registerNode( NodeNames::BetweenGroup,
                  ComputeNode{ /*.sNodeName =*/"between_group_setting",
                               /*.sNodeDesc =*/"Extract the between_group setting from the session json.",
                               /*.fFunc=*/&PartialQuarry::setBetweenGroup,
                               /*.vIncomingFunctions =*/{ },
                               /*.vIncomingSession =*/{ { "settings", "replicates", "between_group" } },
                               /*.vSessionsIncomingInPrevious =*/{ },
                               /*bHidden =*/true } );

    registerNode( NodeNames::FlatValues,
                  ComputeNode{ /*.sNodeName =*/"flat_bins",
                               /*.sNodeDesc =*/"Combine the bin values of the datapools.",
                               /*.fFunc=*/&PartialQuarry::setFlatValues,
                               /*.vIncomingFunctions =*/{ NodeNames::BinValues, NodeNames::InGroup },
                               /*.vIncomingSession =*/{ },
                               /*.vSessionsIncomingInPrevious =*/{ },
                               /*bHidden =*/false } );

    registerNode( NodeNames::PloidyValues,
                  ComputeNode{ /*.sNodeName =*/"ploidy_corr",
                               /*.sNodeDesc =*/"Compute the ploidy corrected bin values.",
                               /*.fFunc=*/&PartialQuarry::setPloidyValues,
                               /*.vIncomingFunctions =*/{ NodeNames::FlatValues },
                               /*.vIncomingSession =*/{ },
                               /*.vSessionsIncomingInPrevious =*/{ },
                               /*bHidden =*/false } );

    registerNode( NodeNames::FlatDecay,
                  ComputeNode{ /*.sNodeName =*/"flat_decay",
                               /*.sNodeDesc =*/"Combine the DDD coverage of the datapools.",
                               /*.fFunc=*/&PartialQuarry::setFlatDecay,
                               /*.vIncomingFunctions =*/{ NodeNames::DecayValues, NodeNames::InGroup },
                               /*.vIncomingSession =*/{ },
                               /*.vSessionsIncomingInPrevious =*/{ },
                               /*bHidden =*/false } );

    registerNode( NodeNames::DecayCDS,
                  ComputeNode{ /*.sNodeName =*/"decay_cds",
                               /*.sNodeDesc =*/"Generate a python ColumnDataSource representation of the DDD coverage.",
                               /*.fFunc=*/&PartialQuarry::setDecayCDS,
                               /*.vIncomingFunctions =*/{ NodeNames::FlatDecay, NodeNames::AnnotationColors },
                               /*.vIncomingSession =*/{ { "settings", "interface", "axis_label_max_char", "val" } },
                               /*.vSessionsIncomingInPrevious =*/{ { "dividend" } },
                               /*bHidden =*/false } );
}

} // namespace cm